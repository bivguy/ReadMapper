# cython: language_level=3
# cython: boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False

from typing import List, Tuple  # ok to import; only used for hints
cimport cython

# Import your existing Python classes/modules
from extend.extender import Extender, Alignment
from seed.minimizer import Minimizer
from models.read import Read
from constants.constants import KMERSIZE, WINDOWSIZE

# --- per-worker module globals (set once in initializer) ---
cdef object _REFERENCE_INDEX  # dict-like; keep as generic object
cdef unicode _REFERENCE_STRING
cdef object _MINIMIZER
cdef object _EXTENDER

@cython.profile(False)
cpdef void _init_worker(object referenceIndex, unicode referenceString):
    """
    Called once per worker process via multiprocessing.Pool(initializer=...)
    Builds and caches heavy, read-only objects in module globals.
    """
    global _REFERENCE_INDEX, _REFERENCE_STRING, _MINIMIZER, _EXTENDER
    _REFERENCE_INDEX  = referenceIndex
    _REFERENCE_STRING = referenceString
    _MINIMIZER = Minimizer(k=KMERSIZE, w=WINDOWSIZE, reference_index=_REFERENCE_INDEX)
    _EXTENDER  = Extender()
    return

cdef inline unsigned char _comp_base(unsigned char b) nogil:
    if b == 65: return 84   #   A->T
    elif b == 84: return 65 #   T->A
    elif b == 67: return 71 #   C->G
    elif b == 71: return 67 #   G->C
    else: return 78         #   N->N

cdef str rc(str seq):
    """
    Reverse complement, uppercase.
    """
    cdef bytes up = seq.upper().encode("ascii")
    cdef Py_ssize_t L = len(up)
    cdef bytearray out = bytearray(L)
    cdef Py_ssize_t i
    for i in range(L):
        out[i] = _comp_base(up[L - 1 - i])
    return out.decode("ascii")


cpdef int compute_sam_flag(bint is_read1, object current, object mate):
    """
    Calculates the SAM flag using C types for performance.
    """
    cdef int flag = 0
    
    # Pre-fetch attributes to C variables to avoid repeated Python lookups
    # (This assumes current and mate have these attributes accessible)
    cdef bint cur_mapped = current.mapped
    cdef bint mate_mapped = mate.mapped
    
    # 1. Read paired (Always 1 for paired-end reads)
    flag += 1

    # 2. Proper pair (Was the pair mapped?)
    if cur_mapped and mate_mapped:
        flag += 2

    # 4. Read unmapped
    if not cur_mapped:
        flag += 4
    # 16. Read reverse strand (Only if mapped)
    # Note: We only check strand if it is mapped
    elif not current.strand_plus: 
        flag += 16

    # 8. Mate unmapped
    if not mate_mapped:
        flag += 8
    # 32. Mate reverse strand (Only if mapped)
    elif not mate.strand_plus:
        flag += 32

    # 64 vs 128. First or Second in pair
    if is_read1:
        flag += 64
    else:
        flag += 128

    return flag

@cython.profile(False)
cpdef list process_read_pair_batch(object batch):
    """
    batch: list of (i, readPair) where readPair is (Read, Read).
    Uses per-worker globals set by _init_worker.
    """
    cdef Py_ssize_t n = len(batch)
    cdef list alignments = [None] * (2 * n)  # preallocate
    cdef Py_ssize_t t, idx = 0

    # bind globals to locals for faster attribute resolution
    cdef object extractor = _MINIMIZER
    cdef object extender  = _EXTENDER
    cdef object refIndex  = _REFERENCE_INDEX
    cdef unicode refStr   = _REFERENCE_STRING

    cdef int i
    cdef object readPair, fRead, bRead
    cdef unicode fReadSeq, bReadSeq
    cdef object frontReadMinimizers, backReadMinimizers
    cdef object frontReadAnchors,  backReadAnchors
    cdef object frontReadAlignment, backReadAlignment

    for t in range(n):
        i, readPair = batch[t]

        # unpack reads
        fRead = readPair[0]
        bRead = readPair[1]

        # extract sequences
        fReadSeq = fRead.getSequence()
        bReadSeq = bRead.getSequence()

        # seed --> lookup
        frontReadMinimizers = extractor.extract(fReadSeq, seq_id=i)
        backReadMinimizers  = extractor.extract(bReadSeq, seq_id=(i + 1) * 2)

        frontReadAnchors = extractor.filter_and_lookup(frontReadMinimizers, refIndex)
        backReadAnchors  = extractor.filter_and_lookup(backReadMinimizers,  refIndex)

        # extend
        frontReadAlignment = extender.extend(fRead.getIdentifier(), fReadSeq, refStr, frontReadAnchors)
        backReadAlignment  = extender.extend(bRead.getIdentifier(), bReadSeq, refStr, backReadAnchors)

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

        # Setting RNEXT and PNEXT fields
        # Update FRONT Read (Look at Back Read)
        if backReadAlignment.mapped:
            frontReadAlignment.rnext = "="                  # Mate is on the same ref
            frontReadAlignment.pnext = backReadAlignment.ref_start
        else:
            frontReadAlignment.rnext = "*"                  # Mate is unmapped
            frontReadAlignment.pnext = 0

        # Update BACK Read (Look at Front Read)
        if frontReadAlignment.mapped:
            backReadAlignment.rnext = "="                   # Mate is on the same ref
            backReadAlignment.pnext = frontReadAlignment.ref_start
        else:
            backReadAlignment.rnext = "*"                   # Mate is unmapped
            backReadAlignment.pnext = 0
        
        # set the SEQ fields
        frontReadAlignment.seq = fReadSeq
        backReadAlignment.seq = bReadSeq

        frontReadAlignment.qual = fRead.getQualityScore()
        backReadAlignment.qual = bRead.getQualityScore()

        # write into preallocated list
        alignments[idx] = frontReadAlignment; idx += 1
        alignments[idx] = backReadAlignment;  idx += 1

    return alignments
