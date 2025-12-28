from typing import List, Tuple, Dict
from collections import defaultdict

Anchor = Tuple[int, int, bool]  # (r, q, same_strand)a

class Chainer:
    def __init__(self):
        pass

    def chain(self, anchors: List[Anchor]) -> List[Anchor]:
        """
        Exact diagonal grouping (no tolerance). Returns the largest group.
        Input anchors are (reference_position=r, query_position=q, same_strand).
        """
        if not anchors:
            return []

        # Ensure deterministic order: sort by (q, r)
        anchors = sorted(anchors, key=lambda a: (a[1], a[0]))

        buckets: Dict[Tuple[bool, int], List[Anchor]] = defaultdict(list)
        for r, q, same in anchors:
            key = (q - r) if same else (q + r)  # invariants: + => q-r,  - => q+r
            buckets[(same, key)].append((r, q, same))

        # Pick largest exact-diagonal bucket
        best: List[Anchor] = []
        for group in buckets.values():
            if len(group) > len(best):
                best = group
        return best
