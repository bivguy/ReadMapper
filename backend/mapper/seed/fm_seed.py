from typing import List, Tuple

class FMSeedExtractor:
    """
    Simple k-mer seeding on top of FMIndex.
    seed_read(seq) -> List[(ref_pos, read_offset)]
    """
    def __init__(self, fm_index, k: int, step: int = None):
        self.fm = fm_index
        self.k = k
        self.step = step if step is not None else k

    def seed_read(self, seq: str) -> List[Tuple[int, int]]:
        hits = []
        n = len(seq)
        i = 0
        while i + self.k <= n:
            kmer = seq[i:i+self.k]
            l, r = self.fm.search(kmer)
            if l <= r:
                for pos in self.fm.locate(l, r):
                    hits.append((pos, i))
            i += self.step
        return hits