from seed.minimizer import Minimizer


class ReferenceIndexBuilder:
    """
    Builds minimizer index from reference genome for fast read mapping
    """
    
    def __init__(self, ref_seq: str, k: int = 15, w: int = 10):
        """
        Initialize index builder
        
        Args:
        k: k-mer size (default 15)
        w: minimizer window size (default 10)
        """
        self.k = k
        self.w = w
        self.ref_seq = ref_seq
        self.extractor = Minimizer(k, w)

    def build_index(self):
        # Build index
        ref_minimizers = self.extractor.extract(self.ref_seq)
        ref_index = {}
        for hash_val, pos, _, is_rev in ref_minimizers:
            if hash_val not in ref_index:
                ref_index[hash_val] = []
            ref_index[hash_val].append((pos, is_rev))

        return ref_index
        