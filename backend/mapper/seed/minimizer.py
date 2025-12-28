from collections import deque
from typing import Dict, List, Tuple

import numpy as np

from hashing.hash import Hash


class Minimizer:
    def __init__(self, k: int, w: int, reference_index: Dict[int, List[Tuple[int, bool]]] = None):
        self.k = k
        self.w = w
        self.hash = Hash(k)       
        self.reference_index = reference_index if reference_index else {}


    def extract(
            self, 
            seq: str, 
            # seq_id: str = "temp",
            seq_id: int = 0
            ) -> List[Tuple[int, int, int, bool]]:
        """
        Extract minimizers from sequence
        Returns: List of (hash, position, seq_id, is_reverse) tuples
        
        Args:
        seq: DNA sequence (handles both uppercase and lowercase)
        seq_id: ID to track which read this is from 
        """
        # Convert to uppercase to handle mixed case input
        seq = seq.upper()
        
        if len(seq) < self.k + self.w - 1:
            return []
    

        minimizers = []
        window = deque()
        
        # Initialize first window
        current_hash = self.hash.hash_sequence(seq[:self.k])
        
        for i in range(self.w):
            if i + self.k > len(seq):
                break
            
            if i > 0:
                current_hash = self.hash.update(current_hash, seq[i-1], seq[i+self.k-1])
            
            # Use min of forward and reverse complement
            kmer = seq[i:i+self.k]
            rev_comp = self._reverse_complement(kmer)
            rev_hash = self.hash.hash_sequence(rev_comp)
            
            if rev_hash < current_hash:
                window.append((rev_hash, i, True))  # True = reverse strand
            else:
                window.append((current_hash, i, False))  # False = forward strand
        
        # Get first minimizer
        # Window contains exactly w k-mers from positions 0 to w-1
        if window:
            min_hash, min_pos, is_rev = min(window) # Find the k-mer with the min hash value in this first window
            minimizers.append((min_hash, min_pos, seq_id, is_rev)) # Add this as our first minimizer to the results
            # Remember this position to avoid duplicates later
            # since same k-mer might be minimum in multiple consecutive windows
            last_min_pos = min_pos
        
        # Reset current has to window slide postion
        current_hash = self.hash.hash_sequence(seq[self.w-1:self.w-1+self.k])
        
        # Slide window
        for i in range(self.w, len(seq) - self.k + 1):
            # Update hash from position i-1 to position i
            # This removes character at seq[i-1] and adds character at seq[i+k-1]
            current_hash = self.hash.update(current_hash, seq[i-1], seq[i+self.k-1])
            
            kmer = seq[i:i+self.k] # Extract the actual k-mer string
            rev_comp = self._reverse_complement(kmer) # Get reverse complement
            rev_hash = self.hash.hash_sequence(rev_comp) # Hash the reverse
            
            # Choose the smaller hash (canonical form)
            # This ensures we match k-mers regardless of strand
            canonical_hash = min(current_hash, rev_hash)
            is_reverse = (rev_hash < current_hash)
            
            # Update slifeing window
            window.popleft()
            window.append((canonical_hash, i, is_reverse))
            
            # Find which k-mer in the current window has minimum hash
            min_hash, min_pos, is_rev = min(window)
            
            # Check if new postion
            if min_pos != last_min_pos:
                minimizers.append((min_hash, min_pos, seq_id, is_rev))
                last_min_pos = min_pos
        
        return minimizers
    
    def filter_and_lookup(self, kmers: List[Tuple[int, int, int, bool]],
                         reference_index: Dict[int, List[Tuple[int, bool]]] = None) -> List[Tuple[int, int, bool]]:
        """
        Filter minimizers against reference and lookup positions
        
        Args:
        kmers: List of hash, position, seq_id, is_reverse tuples from read
        reference_index: Reference index mapping hash -> Dict[int, List[Tuple[int, bool]]]
        Returns:
        List of (ref_pos, read_pos, same_strand) tuples for chaining
        """
        index = reference_index if reference_index is not None else self.reference_index
        
        if not index:
            return []
        
        candidates = []
        
        # Process each read minimizer
        for hash_val, read_pos, _, read_is_rev in kmers:
            # Filter: only process if hash exists in reference
            if hash_val in index:
                # Lookup get all reference positions for this hash
                for ref_pos, ref_is_rev in index[hash_val]:
                    # CHekc strand consistency
                    same_strand = (read_is_rev == ref_is_rev)
                    candidates.append((
                        ref_pos,      
                        read_pos,     
                        same_strand   
                    ))
        
        return candidates
    
    def _reverse_complement(self, seq: str) -> str:
        """
        Get reverse complement of DNA sequence
        Handles both uppercase and lowercase input
        
        Args:
        seq: DNA sequence string
            
        Returns:
        Reverse complement string
        """
        complement = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'A',
            'a': 'T', 't': 'A', 'g': 'C', 'c': 'G', 'n': 'A'
        }
        return ''.join(complement[c] for c in seq[::-1])