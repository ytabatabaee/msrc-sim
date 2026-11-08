from __future__ import annotations

from .lineages import ARGLineage
from .segments import Segment, normalize


def union_lineage(a: ARGLineage, b: ARGLineage, parent_node: int, time: float, new_id: int) -> ARGLineage:
    """Merge all material, including disjoint material, into one ancestor."""
    bounds = sorted({x for s in a.segments + b.segments for x in (s.left, s.right)})
    segments = []
    for left, right in zip(bounds, bounds[1:]):
        if any(s.left <= left and right <= s.right for s in a.segments + b.segments):
            segments.append(Segment(left, right, parent_node))
    return ARGLineage(new_id, a.population, normalize(segments, max(s.right for s in segments)), time, a.background)
