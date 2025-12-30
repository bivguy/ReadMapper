from abc import ABC, abstractmethod
from ..models.readMapper import ReadMapperInput, ReadMapperOutput
from typing import IO, List
import io
from ..index.solutionIndex import SolutionIndexBuilder, MetricAccumulator, Metrics
from ..mmm_parser.parser import Parser 
from ..mmm_parser.readParser import ReadParser
from ..models.read import Read
from ..models.sam import SAM, SAMInput
from ..index.build_index import ReferenceIndexBuilder
from multiprocessing import Pool
from ..parallelization.batch_reads import process_read_pair_batch, _init_worker


class AReadMapper(ABC):
    @abstractmethod
    def mapReads(inputData: ReadMapperInput) -> ReadMapperOutput:
        pass

class ReadMapper(AReadMapper):
    def __init__(self):
        pass

    def mapReads(self, inputData: ReadMapperInput) -> ReadMapperOutput:
        # constants
        k = inputData.kmerSize  # k-mer size
        w = inputData.windowSize  # window size

        # Assume the files need to be opened here 
        outputFile: IO = open(inputData.outputLocation, "w")

        # reference file handling
        referenceFile : IO = io.TextIOWrapper(inputData.referenceGenome, encoding='utf-8')

        referenceStringHeaderLine = referenceFile.readline().strip('\n') # Skip header
        referenceStringHeader = referenceStringHeaderLine.split()[0][1:]
        referenceString = "".join(line.strip() for line in referenceFile)
        referenceFile.close()
        inputData.referenceGenome.close()

        samWriter : SAM = SAM(  # create samOutput
            referenceName=referenceStringHeader, 
            referenceSize = len(referenceString), 
            outputFile=outputFile)


        readFrontFile : IO = io.TextIOWrapper(inputData.readsOne, encoding='utf-8')
        readBackFile : IO = io.TextIOWrapper(inputData.readsTwo, encoding='utf-8')

        # create reference solution map and accumulator for metrics if ground truth is provided
        solutionIndexBuilder : SolutionIndexBuilder = SolutionIndexBuilder()
        if inputData.groundTruth:
            readSolutionFile: IO = open(inputData.groundTruth, "r")
            solutionMap : dict = solutionIndexBuilder.getSolutionMap(readSolutionFile)
            accumulator : MetricAccumulator = MetricAccumulator(solutionMap)
        else:
            readSolutionFile = None
            accumulator = None

       
        # create parser
        parserFront : Parser = Parser(readFile=readFrontFile, referenceFile=None)
        parserBack : Parser = Parser(readFile=readBackFile, referenceFile=None)

        readParser : ReadParser = ReadParser(parserFront, parserBack)
        readPairs : List[List[Read]] = readParser.parseAllReadPairs()

        # Build minimizer index from reference
        builder : ReferenceIndexBuilder = ReferenceIndexBuilder(referenceString, k=k, w=w)
        referenceIndex : dict = builder.build_index()

        # Prepare read pairs with their indices for batch processing
        indexed_read_pairs = list(enumerate(readPairs))
        num_processes = 8 
        batch_size = max(1, len(indexed_read_pairs) // (num_processes * 3))
        
        # Split into batches
        batches = []
        for i in range(0, len(indexed_read_pairs), batch_size):
            batch = indexed_read_pairs[i:i + batch_size]
            batches.append((batch))
        
        totalReads = 0
        mappedReads = 0
        with Pool(
            processes=num_processes,
            initializer=_init_worker,
            initargs=(referenceIndex, referenceString),
            ) as pool:
            batch_results = pool.map(process_read_pair_batch, batches)
            # Flatten the results
            for batch_alignments in batch_results:
                if accumulator:
                    accumulator.update(batch_alignments=batch_alignments)
                totalReads += len(batch_alignments)
                for a in batch_alignments:
                        input : SAMInput
                        if not a.mapped:
                            input = SAMInput(
                                    QNAME = a.readId, 
                                    FLAG= a.flag,
                                    RNAME = referenceStringHeader,
                                    POS = -1,
                                    MAPQ = a.mapq,
                                    RNEXT= a.rnext,
                                    PNEXT=  a.pnext,
                                    CIGAR= "*",
                                    TLEN = -1,
                                    QUAL = a.qual,
                                )
                        else:
                            mappedReads += 1
                            input = SAMInput(
                                    QNAME = a.readId, 
                                    FLAG= a.flag,
                                    RNAME = referenceStringHeader,
                                    POS = a.ref_start,
                                    MAPQ = a.mapq,
                                    CIGAR= a.cigar,
                                    RNEXT= a.rnext, # Ref. name of the mate/next read
                                    PNEXT=  a.pnext,  # Position of the mate/next read
                                    TLEN = a.ref_end - a.ref_start,
                                    QUAL= a.qual,
                                )
                        samWriter.WriteReadToSam(input)
        
        output : ReadMapperOutput = ReadMapperOutput(
            samOutput=outputFile,
            numberOfMappedReads= mappedReads
        )

        outputFile.close()
        readFrontFile.close()
        readBackFile.close()

        return output
