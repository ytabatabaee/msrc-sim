from dataclasses import dataclass

@dataclass(frozen=True)
class Rearrangement:
    rearrangement_id: str = "inv1"
    kind: str = "inversion"
    origin_branch: str = "root"
    initial_copies: int = 1
    selection_coefficient: float = 0.0
    recombination_rate: float = 0.0
    suppression_factor: float = 1.0

    @property
    def effective_recombination_rate(self) -> float:
        return self.recombination_rate*self.suppression_factor
