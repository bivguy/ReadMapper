from extend.extender import Extender, Alignment
from seed.minimizer import Minimizer
from models.read import Read
from typing import Tuple, List
from constants.constants import KMERSIZE, WINDOWSIZE

_REFERENCE_INDEX : dict
_REFERENCE_STRING : str

def _init_worker(referenceIndex : dict, referenceString : str):
    global _REFERENCE_INDEX,_REFERENCE_STRING 

    _REFERENCE_INDEX= referenceIndex
    _REFERENCE_STRING = referenceString

def compute_sam_flag(is_read1: bool, current: Alignment, mate: Alignment) -> int:
    """
    Calculates the SAM flag based on TA specifications.
    
    Bit meanings:
    1   = Read paired (Always True for this project)
    2   = Proper pair (Both reads mapped)
    4   = Read unmapped
    8   = Mate unmapped
    16  = Read reverse strand
    32  = Mate reverse strand
    64  = First in pair
    128 = Second in pair
    """
    flag = 0

    # 1. read paired (Always 1 for paired-end reads)
    flag += 1

    # 2. proper pair (Was the pair mapped?)
    if current.mapped and mate.mapped:
        flag += 2

    # 4. read unmapped
    if not current.mapped:
        flag += 4
    
    # 16. read reverse strand (Only if mapped)
    elif not current.strand_plus: 
        flag += 16

    # 8. mate unmapped
    if not mate.mapped:
        flag += 8
    # 32. mate reverse strand (Only if mapped)
    elif not mate.strand_plus:
        flag += 32

    # 64 vs 128. First or Second in pair
    if is_read1:
        flag += 64
    else:
        flag += 128

    return flag

def process_read_pair_batch(args):
    """Process a batch of read pairs in parallel"""
    global _REFERENCE_INDEX,_REFERENCE_STRING
    batch = args
    k = KMERSIZE
    w = WINDOWSIZE
    
    # Create instances for this process
    extractor = Minimizer(k=k, w=w, reference_index=_REFERENCE_INDEX)
    extender = Extender()
    alignments = []
    
    for i, readPair in batch:
        fReadSeq = readPair[0].getSequence()
        bReadSeq = readPair[1].getSequence()
        frontReadMinimizers = extractor.extract(
            fReadSeq, 
            seq_id=i
        )  
        backReadMinimizers = extractor.extract(
            bReadSeq, 
            seq_id=(i+1)*2
        )  

        frontReadAnchors = extractor.filter_and_lookup(frontReadMinimizers, _REFERENCE_INDEX)
        backReadAnchors = extractor.filter_and_lookup(backReadMinimizers, _REFERENCE_INDEX)

        frontReadAlignment = extender.extend(readPair[0].getIdentifier(), fReadSeq, _REFERENCE_STRING, frontReadAnchors)
        backReadAlignment = extender.extend(readPair[1].getIdentifier(), bReadSeq, _REFERENCE_STRING, backReadAnchors)

        frontReadAlignment.flag = compute_sam_flag(
                is_read1=True, 
                current=frontReadAlignment, 
                mate=backReadAlignment
            )

        backReadAlignment.flag = compute_sam_flag(
            is_read1=False, 
            current=backReadAlignment, 
            mate=frontReadAlignment
        )

        if not frontReadAlignment.mapped:
            frontReadAlignment.mapq = 0
        if not backReadAlignment.mapped:
            backReadAlignment.mapq = 0
        
        frontReadAlignment.qual = readPair[0].getQualityScore()
        backReadAlignment.qual = readPair[1].getQualityScore()
        
        alignments.extend([frontReadAlignment, backReadAlignment])
    
    return alignments
