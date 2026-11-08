from __future__ import annotations

from .lineages import ARGLineage
from .segments import Segment, interval_union, normalize


def union_lineage(a: ARGLineage, b: ARGLineage, parent_node: int | None, time: float, new_id: int) -> ARGLineage:
    """Merge carriers without changing already resolved/disjoint local roots.

    A common-ancestor event changes the local genealogy only where both
    lineages carry material. Disjoint material is carried forward with its
    existing node, avoiding artificial unary ancestry above a completed root.
    """
    bounds = sorted({x for s in a.segments + b.segments for x in (s.left, s.right)})
    segments = []
    for left, right in zip(bounds, bounds[1:]):
        in_a = [s for s in a.segments if s.left <= left and right <= s.right]
        in_b = [s for s in b.segments if s.left <= left and right <= s.right]
        if in_a and in_b and parent_node is not None:
            segments.append(Segment(left, right, parent_node))
        elif in_a:
            segments.append(Segment(left, right, in_a[0].node))
        elif in_b:
            segments.append(Segment(left, right, in_b[0].node))
    active = interval_union(tuple(a.active) + tuple(b.active))
    return ARGLineage(new_id, a.population, normalize(segments, max(s.right for s in segments)), time, a.background, active)
