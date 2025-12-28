# cython: language_level=3
# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

from cpython.list cimport PyList_GET_SIZE
# from libc.stdlib cimport  bint
from typing import List, Tuple


Anchor = Tuple[int, int, bool]  # (r, q, same_strand)a

cdef class Chainer:
    def __cinit__(self):
        pass

    # TODO: try doing cdef for this function at some point
    def chain(self, anchors: List[Anchor]) -> List[Anchor]:
        """
        Exact diagonal grouping (no tolerance). Returns the largest group.
        Input anchors are (r, q, same_strand) where:
           r: int (reference_position)
           q: int (query_position)
           same_strand: bool
        Returns: list[tuple[int, int, bool]]
        """
        if anchors is None:
            return []
        if not anchors:
            return []

        # Ensure deterministic order: sort by (q, r)
        # Keep Python tuple sorting (fast in CPython, fine to call from Cython)
        cdef list anchors_list = list(anchors)
        anchors_list.sort(key=lambda a: (a[1], a[0]))

        cdef dict buckets = {}     # key: (same_strand: bool, diag_key: int) -> list[anchors]
        cdef object a
        cdef int r, q
        cdef bint same
        cdef tuple k
        cdef list L

        # Group by invariant:
        #   same_strand == True  -> key = q - r
        #   same_strand == False -> key = q + r
        for a in anchors_list:
            r = <int>a[0]
            q = <int>a[1]
            same = <bint>a[2]
            if same:
                k = (True, q - r)
            else:
                k = (False, q + r)

            L = <list>buckets.get(k)
            if L is None:
                L = []
                buckets[k] = L
            # Store as a plain Python tuple to match your original Anchor type
            L.append((r, q, bool(same)))

        # Pick largest exact-diagonal bucket
        cdef list best = []
        for L in buckets.values():
            if PyList_GET_SIZE(L) > PyList_GET_SIZE(best):
                best = L

        return best
