from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any
from .chainer import Chainer
from array import array
import numpy as np
import time


# TODO: 
# 1. have a reject high edit rates
# 2. Check the band size and filling
# 3. check the bucket size for window

# Public datatypes
@dataclass
class Alignment:
    readId: str             # the id of the read
    ref_start: int          # absolute reference start 
    ref_end: int            # absolute reference end
    strand_plus: bool       # True = '+', False = '-'
    cigar: str              # SAM CIGAR over the read
    mapped: bool            # if this represens a read that was able to be mapped to the reference 
    # New fields needed for SAM Output
    flag: int = 0           # The calculated bitwise flag
    mapq: int = 60          # Default high quality (60) if mapped, 0 if unmapped
    rnext: str = "*"        # Ref. name of the mate/next read
    pnext: int = 0          # Position of the mate/next read
    qual: str = "*"         # ASCII of Phred-scaled base quality +33


# Helper
def rc(seq: str) -> str:
    comp = {'A': 'T', 
            'T': 'A', 
            'G': 'C', 
            'C': 'G'}
    return ''.join(comp.get(c) for c in seq[::-1].upper())

def compress_cigar(ops: List[str]) -> str:
    if not ops:
        return "0M"
    out, run_c, run_o = [], 1, ops[0]
    for o in ops[1:]:
        if o == run_o:
            run_c += 1
        else:
            out.append(f"{run_c}{run_o}")
            run_o, run_c = o, 1
    out.append(f"{run_c}{run_o}")
    return ''.join(out)

def construct_extension_window(chain: List[Tuple], read_len: int, pad_bp: int) -> Tuple[int, int, bool]:
    # Project every anchor to a [start, end) span on the reference
    r0, q0, same = chain[0]

    if same:
        c = r0 - q0
        ref_lo = c - pad_bp
        ref_hi = c + read_len + pad_bp
        strand_plus = True
    else:
        c = r0 + q0
        ref_lo = (c - read_len) - pad_bp
        ref_hi = c + pad_bp
        strand_plus = False

    return ref_lo, ref_hi, strand_plus

# Extender implementation
class Extender:
    """
    1) Chains anchors (ref_pos, read_pos, same_strand) to pick a best colinear seed.
    2) Creates a small reference window around the chain.
    3) Semi-global (global-on-read) banded edit-distance alignment within that window.
       - costs: match=0, mismatch=1, gap=1
       - returns: Alignment with CIGAR and absolute coords.
    """
    def __init__(self, max_edit_rate: float = 0.40):
        self.chainer = Chainer()
        self.max_edit_rate = max_edit_rate

    # Public API
    def extend(self,
               readId: str,
               read: str,
               reference: str,
               anchors: List[Tuple[int, int, bool]]) -> Optional[Alignment]:
        """
        :param read:      read sequence (A/C/G/T/N)
        :param reference: whole reference string
        :param anchors:   list of (ref_pos, read_pos, same_strand)
        :return: Alignment or None if no good chain
        """
        startTime = time.perf_counter()
        invalidAlignment : Alignment = Alignment(
            readId=readId,
            ref_start=-1,
            ref_end=-1,
            strand_plus=False,
            cigar="",
            mapped=False
        )
        chain = self.chainer.chain(anchors)
        if not chain:
            return invalidAlignment

        read_len = len(read)
        ref_len = len(reference)
        min_pad = 10

        # 1) Window around chain (then clamp to reference)
        ref_lo, ref_hi, strand_plus = construct_extension_window(chain, read_len, min_pad)
        ref_lo = max(0, ref_lo)
        ref_hi = min(ref_len, ref_hi)
        if ref_hi <= ref_lo:
            return invalidAlignment

        # # 2) Orient read
        # q_seq = read if strand_plus else rc(read)
        t_seq = reference[ref_lo:ref_hi]
        diag = min_pad

        if strand_plus:
            q_seq = read
        else:
            q_seq = rc(read)

        # 3) Call DP
        res = self._banded_semiglobal(q_seq, t_seq, diag_est_local=diag, band=15)
        if not res:
            return invalidAlignment
        
        score = res["score"]                     # total edit distance over the whole read
        m = max(1, len(q_seq))                   # safety: avoid div-by-zero for empty reads
        edit_rate = score / m
        if edit_rate > self.max_edit_rate:
            return invalidAlignment

        # 4) Map back to absolute coords and build CIGAR in forward-read order
        ref_start = ref_lo + res["t_start"]
        ref_end   = ref_lo + res["t_end"]

        ops = res["ops"]
        if not strand_plus:
            ops = ops[::-1]   # express operations over the original (forward) read
        cigar = compress_cigar(ops)

        endTime = time.perf_counter()
        elapsedTime = endTime - startTime

        return Alignment(
            readId=readId,
            ref_start=ref_start,
            ref_end=ref_end,
            strand_plus=strand_plus,
            cigar=cigar,
            mapped=True
        )

    # Banded semi-global alignment
    def _banded_semiglobal(self, q: str, t: str, diag_est_local: Optional[int], band: Optional[int]):
        """
        Semi-global: global on read, local on reference.
        Edit distance costs: match 0, mismatch 1, gap 1.
        Banded around j ≈ i + diag_est_local.
        Returns dict with score, ops (list of 'M','I','D'), t_start/t_end, and t_steps for MD.
        """
        m, n = len(q), len(t)
        if m == 0:
            return {"score": 0, "ops": [], "t_start": 0, "t_end": 0, "t_steps": []}
        INF = 10**9

        # Row 0: free start anywhere on reference
        dp_prev = [0] * (n + 1)
        parents: list[tuple[int,int,list[int]]] = []

        edge_slack = 0

        for i in range(1, m + 1):
            dp_cur = [INF] * (n + 1)

            if diag_est_local is None or band is None:
                # Unbanded row
                j_lo, j_hi = 0, n
            else:
                center = max(0, min(n, i + diag_est_local))
                j_lo = max(0, center - (band + edge_slack))
                j_hi = min(n, center + (band + edge_slack))
                if j_lo > j_hi:
                    j_lo = j_hi = center

            row_dir = [0] * (j_hi - j_lo + 1)

            for j in range(j_lo, j_hi + 1):
                ins = dp_prev[j] + 1
                left = dp_cur[j - 1] + 1 if j > 0 else INF
                if j > 0:
                    qi = q[i - 1].upper()
                    tj = t[j - 1].upper()
                    sub = 0 if (qi == tj and qi != 'N') else 1
                    diag = dp_prev[j - 1] + sub
                else:
                    diag = INF

                if diag <= left and diag <= ins:
                    dp_cur[j] = diag
                    row_dir[j - j_lo] = 2  # 'M'
                elif left <= ins:
                    dp_cur[j] = left
                    row_dir[j - j_lo] = 1  # 'D'
                else:
                    dp_cur[j] = ins
                    row_dir[j - j_lo] = 0  # 'I'

            parents.append((j_lo, j_hi, row_dir))
            dp_prev = dp_cur

        # end anywhere on reference
        j_end = min(range(n + 1), key=lambda j: dp_prev[j])
        score = dp_prev[j_end]

        # backtrack
        i, j = m, j_end
        ops_rev, t_steps_rev = [], []
        while i > 0:
            j_lo, j_hi, row_dir = parents[i - 1]
            if j < j_lo or j > j_hi:
                # if we somehow step outside the stored band, fall back by recomputing
                # should be rare if band is reasonable, can figure out later
                return None
            else:
                d = row_dir[j - j_lo]

            if d == 0:        # 'I' from (i-1, j)
                ops_rev.append('I')
                t_steps_rev.append(j)
                i, j = i - 1, j
            elif d == 1:      # 'D' from (i, j-1)
                ops_rev.append('D')
                t_steps_rev.append(j)
                i, j = i, j - 1
            elif d == 2:      # 'M' from (i-1, j-1)
                ops_rev.append('M')
                t_steps_rev.append(j)
                i, j = i - 1, j - 1
            else:
                return None

        ops = ops_rev[::-1]
        t_steps = t_steps_rev[::-1]

        t_used = [s - 1 for (o, s) in zip(ops, t_steps) if o in ('M', 'D')]
        if t_used:
            t_start = max(0, min(t_used))
            t_end   = max(t_used) + 1
        else:
            t_start = t_end = j_end

        return {"score": score, "ops": ops, "t_start": t_start, "t_end": t_end, "t_steps": t_steps}