from __future__ import annotations

from typing import Protocol, Any


class BackgroundModel(Protocol):
    """Extension point for MSRC backgrounds.

    A future implementation may return several state-specific hazards or split a
    lineage into segment-associated states. The neutral implementation has one
    state and no independent switching events.
    """

    name: str

    def initial_state(self, population: str, segment: Any, taxon: str | None = None) -> Any: ...

    def coalescence_rate(self, population: str, state_a: Any, state_b: Any, effective_size: float, time: float) -> float: ...

    def recombination_states(self, population: str, state: Any, left: Any, right: Any, time: float,
                             position: float | None = None, parent_background: Any = None) -> tuple[Any, Any]: ...
