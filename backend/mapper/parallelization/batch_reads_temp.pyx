# file: extend/worker.pyx
# cython: language_level=3
# cython: boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False

from typing import List, Tuple
cimport cython

# Import your existing classes (Python classes are fine to import this way)
from extend.extender import Extender, Alignment
from seed.minimizer import Minimizer
from models.read import Read

@cython.profile(False)
cpdef list process_read_pair_batch(object args):
    cdef object batch
    cdef int k, w
    cdef object referenceIndex
    cdef str referenceString
    batch, k, w, referenceIndex, referenceString = args

    cdef object extractor = Minimizer(k=k, w=w, reference_index=referenceIndex)
    cdef object extender = Extender()

    cdef Py_ssize_t n = len(batch)
    cdef list alignments = [None] * (2 * n)   # preallocate
    cdef Py_ssize_t t, idx = 0

    cdef int i
    cdef object readPair, fRead, bRead
    cdef str fReadSeq, bReadSeq
    cdef object frontReadMinimizers, backReadMinimizers
    cdef object frontReadAnchors, backReadAnchors
    cdef object frontReadAlignment, backReadAlignment

    for t in range(n):
        i, readPair = batch[t]
        fRead = readPair[0]
        bRead = readPair[1]

        fReadSeq = fRead.getSequence()
        bReadSeq = bRead.getSequence()

        frontReadMinimizers = extractor.extract(fReadSeq, seq_id=i)
        backReadMinimizers  = extractor.extract(bReadSeq, seq_id=(i + 1) * 2)

        frontReadAnchors = extractor.filter_and_lookup(frontReadMinimizers, referenceIndex)
        backReadAnchors  = extractor.filter_and_lookup(backReadMinimizers,  referenceIndex)

        frontReadAlignment = extender.extend(fRead.getIdentifier(), fReadSeq, referenceString, frontReadAnchors)
        backReadAlignment  = extender.extend(bRead.getIdentifier(), bReadSeq, referenceString, backReadAnchors)

        alignments[idx] = frontReadAlignment; idx += 1
        alignments[idx] = backReadAlignment;  idx += 1

    return alignments

