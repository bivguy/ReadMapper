import numpy as np


class Hash:
    def __init__(self, k: int):
        self.k = k
        self.base = 4
        self.mod = 2**64 - 1
        self.power = self.base**(k-1)

        self.encode = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

    def hash_sequence(self, s: str):
        h = 0
        for c in s:
            h = (h * self.base + self.encode.get(c, 0)) & self.mod
        return h

    def hash_binary(self, binary_s: np.ndarray):
        """
        If we wish to perform bit operations later on for speed up through usage of binary seqs, 
        however does not need to be implemented rn.
        """
        h = 0
        for val in binary_s:
            h = (h * self.base + int(val)) & self.mod
        return h

    def update(self, prev_hash: int, out_char: str, in_char: str) -> int:
        """
        Rolling hash update - removes leftmost character and adds rightmost character
        
        Args:
            prev_hash: The previous hash value
            out_char: Character leaving the window (leftmost)
            in_char: Character entering the window (rightmost)
            
        Returns:
            Updated hash value
        """
        out_val = self.encode.get(out_char, 0)
        in_val = self.encode.get(in_char, 0)
        
        # Remove contribution of outgoing character
        h = (prev_hash - out_val * self.power) & self.mod
        
        # Shift left and add incoming character
        h = (h * self.base + in_val) & self.mod
        
        return h
    
    def update_binary(self, prev_hash: int, out_val: int, in_val: int) -> int:
        """
        Rolling hash update - removes leftmost character and adds rightmost character
        
        Args:
            prev_hash: The previous hash value
            out_val: val leaving the window (leftmost)
            in_val: val entering the window (rightmost)
            
        Returns:
            Updated hash value
        """       
        # Remove contribution of outgoing character
        h = (prev_hash - out_val * self.power) & self.mod
        
        # Shift left and add incoming character
        h = (h * self.base + in_val) & self.mod
        
        return h