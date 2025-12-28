from seed.minimizer import Minimizer
from index.build_index import ReferenceIndexBuilder
from constants.constants import KMERSIZE, WINDOWSIZE
from parser.parser import Parser
from extend.extender import Extender
from extend.chainer import Chainer
from models.sam import SAM 
from typing import IO, List


class ReadStub:
    def __init__(self, s: str, identifier: str):
        self._s = s
        self.identifier = identifier
    def getFullRead(self) -> str:
        return self._s
    def getIdentifier(self) -> str:
        return self.identifier


def _read_fastq_sequences(path: str) -> List[ReadStub]:
    seqs: List[ReadStub] = []
    with open(path, "r") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            p = fh.readline()
            q = fh.readline()
            if not s or not p or not q:
                break
            seqs.append(ReadStub(s, h))
    return seqs


def main():
    # constants
    k = KMERSIZE  # k-mer size
    w = WINDOWSIZE  # window size

    # open files
    readFile1: IO = open("io/inputs/reads/short_reads_1_1000_subset.fastq", "r")
    readFile2: IO = open("io/inputs/reads/short_reads_2_1000_subset.fastq", "r")
    referenceFile: IO = open("io/inputs/ref/short_reads_ref_genome.fasta", "r")
    # outputFile : IO = open("io/outputs/out.sam", "w")

    # create parser (reference only)
    parserFront: Parser = Parser(readFile=None, referenceFile=referenceFile)
    referenceString: str = parserFront.getReferenceString()

    # load reads (single-end)
    reads1: List[ReadStub] = _read_fastq_sequences("io/inputs/reads/short_reads_1_1000_subset.fastq")
    reads2: List[ReadStub] = _read_fastq_sequences("io/inputs/reads/short_reads_2_1000_subset.fastq")

    # Build minimizer index from reference
    builder: ReferenceIndexBuilder = ReferenceIndexBuilder(referenceString, k=k, w=w)
    referenceIndex: dict = builder.build_index()

    # seed the reads
    extractor : Minimizer = Minimizer(k=k, w=w, reference_index=referenceIndex)
    chainer = Chainer()
    extender = Extender()
    alignments = []
    anchor_lst = []
    chain_lst = []

    # iterate through each read 
    for index in range(len(reads1)):

        # Extract minimizers from read        
        minimizers1 = extractor.extract(reads1[index].getFullRead(), seq_id= index*2)

        # Filter and lookup candidates
        anchors1 = extractor.filter_and_lookup(minimizers1, referenceIndex)
        anchor_lst.append(anchors1)

        chain1 = chainer.chain(anchors1)
        chain_lst.append(chain1)

        # Handles chaining, local alignment, selects and outputs best alignment
        alignment1 = extender.extend(reads1[index].getIdentifier(), reads1[index].getFullRead(), referenceString, anchors1)
        alignments.append(alignment1)


        # Extract minimizers from read
        minimizers2 = extractor.extract(reads2[index].getFullRead(), seq_id= (index*2) + 1)

        # Filter and lookup candidates
        anchors2 = extractor.filter_and_lookup(minimizers2, referenceIndex)
        anchor_lst.append(anchors2)

        chain2 = chainer.chain(anchors2)
        chain_lst.append(chain2)

        # Handles chaining and local alignment and selects and outputs best alignment
        alignment2 = extender.extend(reads2[index].getIdentifier(), reads2[index].getFullRead(), referenceString, anchors2)
        alignments.append(alignment2)


    counter = 0
    # samWriter : SAM = SAM(outputFile=outputFile)
    for i in range(len(alignments)):
        if alignments[i] is not None:
            check = len(alignments[i].cigar)
            if check <= 8:
                counter += 1

        print(alignments[i])
        # print(anchor_lst[i])
        # print(chain_lst[i])
        # print()

    print(counter)


if __name__ == "__main__":
    main()
