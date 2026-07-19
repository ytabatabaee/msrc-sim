from dataclasses import dataclass
from typing import FrozenSet

@dataclass(frozen=True)
class Lineage:
    lineage_id: int
    descendants: FrozenSet[str]
    arrangement: int

@dataclass(frozen=True)
class CoalescenceEvent:
    time: float
    left_id: int
    right_id: int
    topology_index: int
