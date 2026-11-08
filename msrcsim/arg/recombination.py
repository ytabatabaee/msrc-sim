from __future__ import annotations

from .lineages import ARGLineage
from .recombination_map import PiecewiseConstantRecombinationMap
from .segments import Segment, interval_union, split


def recombine(lineage: ARGLineage, position: float, left_id: int, right_id: int, time: float, rec_map: PiecewiseConstantRecombinationMap, background_model):
    left, right = split(lineage.segments, position)
    if not left or not right:
        raise ValueError("recombination position must split ancestral material")
    state_left, state_right = background_model.recombination_states(lineage.population, lineage.background, left, right, time)
    # The parental nodes become the current ancestry after the crossover; the
    # old node is connected to them by edges in the caller.
    def remap(parts, parent_id):
        out = []
        for s in parts:
            cuts = sorted({s.left, s.right} | {x for a, b in lineage.active for x in (max(a, s.left), min(b, s.right)) if max(a, s.left) < min(b, s.right)})
            for a, b in zip(cuts, cuts[1:]):
                active = any(x <= a and b <= y for x, y in lineage.active)
                out.append(Segment(a, b, parent_id if active else s.node))
        return tuple(out)
    left = remap(left, left_id)
    right = remap(right, right_id)
    active_left = []
    active_right = []
    for a, b in lineage.active:
        if a < position:
            active_left.append((a, min(b, position)))
        if b > position:
            active_right.append((max(a, position), b))
    return (ARGLineage(left_id, lineage.population, left, time, state_left, interval_union(active_left)),
            ARGLineage(right_id, lineage.population, right, time, state_right, interval_union(active_right)))
