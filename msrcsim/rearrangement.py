from dataclasses import dataclass
@dataclass(frozen=True)
class Rearrangement:
    rearrangement_id: str
    kind: str
    origin_branch: str
    origin_time_from_branch_start: int
    initial_copy_count: int=1
    selection_coefficient: float=0.0
