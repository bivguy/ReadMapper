# cython: language_level=3
# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

# Optional: only needed if you want to accept numpy arrays directly in Python calls
import numpy as np
from libc.stdint cimport uint64_t, int32_t

cdef inline uint64_t _encode_nt_char(str ch) nogil:
    """
    Map nucleotide char -> {A:0, T:1, G:2, C:3}, default 0 for anything else.
    Accepts upper/lowercase.
    """
    # ch is a 1-length Python str; direct equality checks are fine
    if ch == 'A' or ch == 'a':
        return <uint64_t>0
    elif ch == 'T' or ch == 't':
        return <uint64_t>1
    elif ch == 'G' or ch == 'g':
        return <uint64_t>2
    elif ch == 'C' or ch == 'c':
        return <uint64_t>3
    else:
        return <uint64_t>0

cdef inline uint64_t _encode_nt_codepoint(Py_UCS4 cp) nogil:
    """
    Faster branch on codepoints when iterating a Python str with Py_UCS4.
    """
    # 'A' = 65, 'a' = 97, 'T' = 84, 't' = 116, 'G' = 71, 'g' = 103, 'C' = 67, 'c' = 99
    if cp == 65 or cp == 97:
        return <uint64_t>0
    elif cp == 84 or cp == 116:
        return <uint64_t>1
    elif cp == 71 or cp == 103:
        return <uint64_t>2
    elif cp == 67 or cp == 99:
        return <uint64_t>3
    else:
        return <uint64_t>0

cdef class Hash:
    cdef:
        uint64_t k
        uint64_t base
        uint64_t mod
        uint64_t power

    def __cinit__(self, int k):
        self.k = <uint64_t>k
        self.base = <uint64_t>4
        self.mod = <uint64_t>0xFFFFFFFFFFFFFFFF  # 2**64 - 1
        # base**(k-1) (done in Python int then cast to uint64_t; fits since we mask anyway)
        self.power = (<uint64_t>(4 ** max(k - 1, 0))) & self.mod

    #
    # Public API mirrors your Python class
    #

    cpdef uint64_t hash_sequence(self, str s):
        """
        Compute rolling-base hash over a nucleotide string.
        Unknown characters map to 0 (same as dict.get(c, 0)).
        """
        cdef:
            uint64_t h = 0
            uint64_t base = self.base
            uint64_t mod = self.mod
            Py_UCS4 cp
            Py_ssize_t i, n = len(s)

        # Iterate by codepoint to avoid per-iteration Python object allocations
        for i in range(n):
            cp = <Py_UCS4> s[i]
            h = (h * base + _encode_nt_codepoint(cp)) & mod
        return h

    cpdef uint64_t hash_binary(self, object binary_s):
        """
        Same as your numpy version, but accepts ANY 1-D integer-like buffer:
          - NumPy arrays (int8/int16/int32/int64/uint8/...)
          - Python array('i'), list of ints (will copy), etc.

        For best performance, pass a contiguous NumPy array (e.g., dtype=int32/int64).
        """
        cdef uint64_t h = 0
        cdef uint64_t base = self.base
        cdef uint64_t mod = self.mod

        # Create a typed memoryview (will zero-copy for contiguous NumPy arrays)
        cdef long[:] view = np.asarray(binary_s, dtype=np.int64)  # ensures numeric and contiguous
        cdef Py_ssize_t i, n = view.shape[0]
        for i in range(n):
            h = (h * base + <uint64_t>view[i]) & mod
        return h

    cpdef uint64_t update(self, uint64_t prev_hash, str out_char, str in_char):
        """
        Rolling update for a window shift: remove out_char, add in_char.
        """
        cdef uint64_t out_val = _encode_nt_char(out_char)
        cdef uint64_t in_val  = _encode_nt_char(in_char)
        cdef uint64_t h

        # Remove outgoing char: prev - out_val * base^(k-1)
        h = (prev_hash - (out_val * self.power)) & self.mod

        # Shift/add incoming
        h = (h * self.base + in_val) & self.mod
        return h

    cpdef uint64_t update_binary(self, uint64_t prev_hash, int out_val, int in_val):
        """
        Rolling update when you already have integer-encoded symbols.
        """
        cdef uint64_t h
        h = (prev_hash - (<uint64_t>out_val * self.power)) & self.mod
        h = (h * self.base + <uint64_t>in_val) & self.mod
        return h
