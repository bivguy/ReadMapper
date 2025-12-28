"""
Main.py:
Handles pipeline and CLI
"""

import argparse
from multiprocessing import Pool, cpu_count
from parallelization.batch_reads import process_read_pair_batch, _init_worker
from seed.minimizer import Minimizer
from index.build_index import ReferenceIndexBuilder
from constants.constants import KMERSIZE, WINDOWSIZE
from mmm_parser.parser import Parser 
from mmm_parser.readParser import ReadParser
from models.read import FullRead, Read
from extend.extender import Extender, Alignment
from models.sam import SAM, SAMInput
from typing import IO, List
from index.solutionIndex import SolutionIndexBuilder, MetricAccumulator, Metrics
import time
# import tracemalloc
import sys
import psutil
import os

def get_memory_usage():
    """Get current memory usage in bytes using psutil"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss

def parse_args():
    parser = argparse.ArgumentParser(description='Genome Read Mapper')
    
    # Required arguments
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-1', '--reads1', required=True, help='First paired-end reads FASTQ file')
    parser.add_argument('-2', '--reads2', required=True, help='Second paired-end reads FASTQ file')
    
    # Optional arguments
    parser.add_argument('-o', '--output', default='io/outputs/SAMOutputFile.SAM', help='Output SAM file')
    parser.add_argument('--truth', help='Ground truth file for metrics')
    parser.add_argument('-k', '--kmer', type=int, default=KMERSIZE, help='K-mer size')
    parser.add_argument('-w', '--window', type=int, default=WINDOWSIZE, help='Window size')
    parser.add_argument('-m', '--memory', action='store_true', help='Track memory usage')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    startTime = time.perf_counter()
    startCpuTime = time.process_time()
    
    if args.memory:
        # tracemalloc.start()
        baseline_memory = get_memory_usage()
        peak_memory = baseline_memory

    # constants
    k = args.kmer  # k-mer size
    w = args.window  # window size

    # if args.memory:
    #     baseline_current_tm, baseline_peak_tm = tracemalloc.get_traced_memory()

    # open file... should need CLI handling
    readFrontFile : IO = open(args.reads1, "r")
    readBackFile : IO = open(args.reads2, "r")
    if args.truth:
        readSolutionFile: IO = open(args.truth, "r")
    else:
        readSolutionFile = None
    referenceFile : IO = open(args.reference, "r")
    outputFile: IO = open(args.output, "w")

    referenceStringHeaderLine = referenceFile.readline().strip('\n') # Skip header
    referenceStringHeader = referenceStringHeaderLine.split()[0][1:]
    referenceString = "".join(line.strip() for line in referenceFile)
    referenceFile.close()
    # create samOutput
    samWriter : SAM = SAM(referenceName=referenceStringHeader, referenceSize = len(referenceString), outputFile=outputFile)

    # create parser
    parserFront : Parser = Parser(readFile=readFrontFile, referenceFile=None)
    parserBack : Parser = Parser(readFile=readBackFile, referenceFile=None)

    readParser : ReadParser = ReadParser(parserFront, parserBack)

    # solution builder
    solutionIndexBuilder : SolutionIndexBuilder = SolutionIndexBuilder()

    # reads : List[FullRead] = readParser.parseAllReads()
    readPairs : List[List[Read]] = readParser.parseAllReadPairs()

    # referenceString : str = parserFront.getReferenceString()
    solutionMap : dict = solutionIndexBuilder.getSolutionMap(readSolutionFile)
    accumulator : MetricAccumulator = MetricAccumulator(solutionMap)
    # Build minimizer index from reference
    builder : ReferenceIndexBuilder = ReferenceIndexBuilder(referenceString, k=k, w=w)
    referenceIndex : dict = builder.build_index()

    # Track Indexing Memory Usage
    if args.memory:
        indexing_memory = get_memory_usage()
        indexing_memory_used = indexing_memory - baseline_memory
        peak_memory = max(peak_memory, indexing_memory)
        
        # first_current_tm, first_peak_tm = tracemalloc.get_traced_memory()
        # first_part_memory_tm = first_current_tm - baseline_current_tm
        
        print(f"\nMemory Usage After Indexing ")
        print(f"Indexing (actual system memory): {indexing_memory_used / 10**6:.2f} MB")
        # print(f"Indexing (tracemalloc Python only): {first_part_memory_tm / 10**6:.2f} MB")
        print(f"Total memory used: {(indexing_memory - baseline_memory) / 10**6:.2f} MB")

    # index time print("hello wrld grild")
    indexTime = time.perf_counter()
    indexCpuTime = time.process_time()

    indexElapsedTime = indexTime - startTime
    indexCpuElapsedTime = indexCpuTime - startCpuTime

    # Prepare read pairs with their indices for batch processing
    indexed_read_pairs = list(enumerate(readPairs))
    # Determine optimal batch size and number of processes
    num_processes = 8  # Use all available CPU cores
    # print(f"Total CPU cores: {num_processes}")
    batch_size = max(1, len(indexed_read_pairs) // (num_processes * 3))
    
    # Split into batches
    batches = []
    for i in range(0, len(indexed_read_pairs), batch_size):
        batch = indexed_read_pairs[i:i + batch_size]
        batches.append((batch))

    total_reads = 0
    with Pool(
        processes=num_processes,
        initializer=_init_worker,
        initargs=(referenceIndex, referenceString),
        ) as pool:
        batch_results = pool.map(process_read_pair_batch, batches)
        # Flatten the results
        for batch_alignments in batch_results:
            accumulator.update(batch_alignments=batch_alignments)
            total_reads += len(batch_alignments)
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

    # finalize the metrics 
    metrics : Metrics =  accumulator.compute_final_metrics(total_reads_processed=total_reads)
    # Mapping and Cumulative Usage
    if args.memory:
        mapping_memory = get_memory_usage()
        mapping_memory_used = mapping_memory - indexing_memory
        peak_memory = max(peak_memory, mapping_memory)
        # latter_current_tm, latter_peak_tm = tracemalloc.get_traced_memory()
        # latter_part_only_tm = latter_current_tm - first_current_tm
        
        print(f"\nMemory Usage After Mapping")
        print(f"Mapping phase (actual system memory): {mapping_memory_used / 10**6:.2f} MB")
        # print(f"Mapping phase (tracemalloc Python only): {latter_part_only_tm / 10**6:.2f} MB")
        print(f"\nTotal Memory Usage")
        print(f"Total cumulative (actual): {(mapping_memory - baseline_memory) / 10**6:.2f} MB")
        # print(f"Total cumulative (tracemalloc): {latter_current_tm / 10**6:.2f} MB")
        print(f"Peak memory (actual): {(peak_memory - baseline_memory) / 10**6:.2f} MB")
        # print(f"Peak memory (tracemalloc): {latter_peak_tm / 10**6:.2f} MB")
    outputFile.close()
    endTime = time.perf_counter()
    endCpuTime = time.process_time()

    elapsedTime = endTime - startTime
    elapsedCpuTime = endCpuTime - startCpuTime
    
    mappingTime = endTime - indexTime
    mappingCpuTime = endCpuTime - indexCpuTime

    # compute the metrics now 
    if readSolutionFile:
        # solutionMap : dict = solutionIndexBuilder.getSolutionMap(readSolutionFile)
        # metrics : Metrics = solutionIndexBuilder.computeMetrics(alignments, solutionMap)
        print(metrics)

    print(f"Finished processing {total_reads} reads.")
    # Print the execution time
    print(f"The 'indexing time' part took {indexElapsedTime:.4f} seconds to execute.")
    print(f"The 'mapping' part took {mappingTime:.4f} seconds to execute.")
    print(f"The 'main' part took {elapsedTime:.4f} seconds to execute.")
    print()

    print(f"The 'indexing cpu time' part took {indexCpuElapsedTime:.4f} seconds to execute.")
    print(f"The 'mapping cpu' part took {mappingCpuTime:.4f} seconds to execute.")
    print(f"The 'main cpu' part took {elapsedCpuTime:.4f} seconds to execute.")

    reads_per_minute = (total_reads * 60) / elapsedTime
    print(f"\nPerformance: {reads_per_minute:.0f} reads per minute")


if __name__ == "__main__":
    main()