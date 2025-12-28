# index/fm_index.py
from collections import defaultdict
from typing import Dict, List, Tuple

class FMIndex:
    """
    Plain FM-index (SA + BWT + C + Occ checkpoints).
    Suitable for moderate references; for very large refs consider a C-extension.
    """
    def __init__(self, s: str, step: int = 128):
        assert s.endswith("$"), "Reference must end with sentinel '$'"
        self.s = s
        self.n = len(s)
        self.step = step
        self.sa = self._suffix_array(s)
        self.bwt = self._bwt_from_sa(s, self.sa)
        self.alphabet = sorted(set(self.bwt))
        self.C = self._build_C(self.bwt)
        self.occ_chk, self.step = self._build_occ(self.bwt, self.alphabet, step)

    @staticmethod
    def _suffix_array(s: str) -> List[int]:
        n = len(s)
        k = 1
        sa = list(range(n))
        rank = [ord(c) for c in s]
        tmp = [0] * n
        while True:
            sa.sort(key=lambda i: (rank[i], rank[i + k] if i + k < n else -1))
            tmp[sa[0]] = 0
            for i in range(1, n):
                a, b = sa[i - 1], sa[i]
                tmp[b] = tmp[a] + (
                    rank[b] != rank[a] or
                    (rank[b + k] if b + k < n else -1) != (rank[a + k] if a + k < n else -1)
                )
            rank, tmp = tmp, rank
            if rank[sa[-1]] == n - 1:
                break
            k <<= 1
        return sa

    @staticmethod
    def _bwt_from_sa(s: str, sa: List[int]) -> str:
        n = len(s)
        out = []
        for p in sa:
            out.append(s[p - 1] if p != 0 else s[-1])
        return "".join(out)

    @staticmethod
    def _build_C(bwt: str) -> Dict[str, int]:
        counts = defaultdict(int)
        for ch in bwt:
            counts[ch] += 1
        total = 0
        C = {}
        for ch in sorted(counts):
            C[ch] = total
            total += counts[ch]
        return C

    @staticmethod
    def _build_occ(bwt: str, alphabet: List[str], step: int):
        n = len(bwt)
        chk = {ch: [0] * ((n + step - 1) // step + 1) for ch in alphabet}
        run = {ch: 0 for ch in alphabet}
        for i, ch in enumerate(bwt):
            if i % step == 0:
                bi = i // step
                for a in alphabet:
                    chk[a][bi] = run[a]
            run[ch] += 1
        bi = n // step
        for a in alphabet:
            chk[a][bi] = run[a]
        return chk, step

    def _occ(self, ch: str, i: int) -> int:
        if i <= 0:
            return 0
        block = i // self.step
        base = self.occ_chk.get(ch, [0])[block]
        start = block * self.step
        rem = 0
        bs = self.bwt
        for k in range(start, i):
            if bs[k] == ch:
                rem += 1
        return base + rem

    def search(self, pat: str) -> Tuple[int, int]:
        if not pat:
            return (0, self.n - 1)
        l, r = 0, self.n - 1
        for ch in reversed(pat):
            if ch not in self.C:
                return (1, 0)
            l = self.C[ch] + self._occ(ch, l)
            r = self.C[ch] + self._occ(ch, r + 1) - 1
            if l > r:
                return (1, 0)
        return (l, r)

    def locate(self, l: int, r: int) -> List[int]:
        if l > r:
            return []
        return [self.sa[i] for i in range(l, r + 1)]