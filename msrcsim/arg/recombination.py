from __future__ import annotations

from .lineages import ARGLineage
from .recombination_map import PiecewiseConstantRecombinationMap
from .segments import Segment, split


def recombine(lineage: ARGLineage, position: float, left_id: int, right_id: int, time: float, rec_map: PiecewiseConstantRecombinationMap, background_model):
    left, right = split(lineage.segments, position)
    if not left or not right:
        raise ValueError("recombination position must split ancestral material")
    state_left, state_right = background_model.recombination_states(lineage.population, lineage.background, left, right, time)
    # The parental nodes become the current ancestry after the crossover; the
    # old node is connected to them by edges in the caller.
    left = tuple(Segment(s.left, s.right, left_id) for s in left)
    right = tuple(Segment(s.left, s.right, right_id) for s in right)
    return (ARGLineage(left_id, lineage.population, left, time, state_left),
            ARGLineage(right_id, lineage.population, right, time, state_right))
