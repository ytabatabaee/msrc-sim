from __future__ import annotations
from dataclasses import dataclass

@dataclass
class Lineage:
    lineage_id: int
    descendants: frozenset[str]
    population_id: str
    arrangement: int
    node_age: float = 0.0
    newick: str | None = None

    def label(self) -> str:
        return next(iter(self.descendants)) if len(self.descendants)==1 else "".join(sorted(self.descendants))
