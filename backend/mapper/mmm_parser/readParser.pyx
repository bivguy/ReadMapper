# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, cdivision=True
# distutils: language = c

from typing import List
# from mmm_parser.parser cimport Parser        # cimport for faster attribute access if compiled
from models.read import FullRead, Read    # regular import, still Python-level

cdef class ReadParser:
    cdef object parserFront
    cdef object parserBack

    def __cinit__(self, object parserFront, object parserBack):
        self.parserFront = parserFront
        self.parserBack = parserBack

    cpdef object parseFullRead(self):
        """Return a FullRead combining one read from each parser (None if either done)."""
        cdef object frontRead = self.parserFront.parseNextRead()
        cdef object backRead = self.parserBack.parseNextRead()

        if frontRead is None or backRead is None:
            return None

        return FullRead(frontRead, backRead)

    cpdef list parseReadPair(self):
        """Return a [frontRead, backRead] pair, or [] at EOF."""
        cdef object frontRead = self.parserFront.parseNextRead()
        cdef object backRead = self.parserBack.parseNextRead()

        if frontRead is None or backRead is None:
            return []
        return [frontRead, backRead]

    cpdef list parseAllReads(self):
        """Return list of FullRead objects until EOF."""
        cdef list reads = []
        cdef object fullRead
        while True:
            fullRead = self.parseFullRead()
            if not fullRead:
                break
            reads.append(fullRead)
        return reads

    cpdef list parseAllReadPairs(self):
        """Return list of [frontRead, backRead] pairs until EOF."""
        cdef list readPairs = []
        cdef list readPair
        while True:
            readPair = self.parseReadPair()
            if not readPair:
                break
            readPairs.append(readPair)
        return readPairs