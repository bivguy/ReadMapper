# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, nonecheck=False, cdivision=True
# distutils: language = c

cimport cython

from constants.constants import KMERSIZE
from models.read import Read  # treated as a Python class

cdef class Parser:
    cdef object readFile
    cdef object referenceFile
    cdef public int kmerSize

    def __cinit__(self, readFile, referenceFile=None):
        self.readFile = readFile
        self.referenceFile = referenceFile
        self.kmerSize = KMERSIZE

    @cython.cfunc
    @cython.inline
    cdef str _rstrip_crlf(self, str s):
        cdef Py_ssize_t n = len(s)
        if n == 0:
            return s
        # strip '\n'
        if s[n-1] == '\n':
            n -= 1
        # strip '\r'
        if n > 0 and s[n-1] == '\r':
            n -= 1
        if n != len(s):
            return s[:n]
        return s

    cpdef str getNextReferenceKmer(self):
        """Read the next k-mer (length kmerSize) from referenceFile, skipping newlines."""
        cdef list res = []
        cdef int i = 0
        cdef str ch
        cdef object ref = self.referenceFile
        if ref is None:
            return ""
        while i < self.kmerSize:
            ch = ref.read(1)
            if not ch:
                break
            if ch == '\n' or ch == '\r':
                continue
            res.append(ch)
            i += 1
        return ''.join(res)

    cpdef object parseNextRead(self):
        """
        Parse the next FASTQ record (4 lines) and return a models.read.Read.
        Returns None at EOF.
        """
        cdef str identifier = ""
        cdef str sequence = ""
        cdef str qualityScore = ""
        cdef bint isFront = True
        cdef int i
        cdef object f = self.readFile
        cdef str line

        for i in range(4):
            line = f.readline()
            if not line:
                return None

            line = self._rstrip_crlf(line)

            if i == 0:
                # header line: "@<id>[/1|/2]" or with trailing 1/2
                # remove first char '@'
                if line:
                    identifier = line[1:]
                    # Determine front/back from the last char if it's '2'
                    isFront = False if identifier and identifier[-1] == '2' else True
                else:
                    identifier = ""
                    isFront = True
            elif i == 1:
                sequence = line
            elif i == 2:
                # '+' line ignored
                pass
            else:
                qualityScore = line

        return Read(identifier=identifier, sequence=sequence, qualityScore=qualityScore, isFront=isFront)
