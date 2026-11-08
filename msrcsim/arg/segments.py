from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class Segment:
    """A half-open ancestral interval and the node currently ancestral to it."""

    left: float
    right: float
    node: int

    def __post_init__(self):
        if self.left < 0 or self.right <= self.left:
            raise ValueError("segments must be non-empty and have right > left")


def normalize(segments: list[Segment] | tuple[Segment, ...], length: float) -> tuple[Segment, ...]:
    ordered = sorted(segments, key=lambda s: (s.left, s.right, s.node))
    out: list[Segment] = []
    for s in ordered:
        if s.right > length + 1e-12:
            raise ValueError("ancestral segment exceeds chromosome length")
        if out and s.left < out[-1].right - 1e-12:
            raise ValueError("segments overlap within a lineage")
        if out and abs(out[-1].right - s.left) <= 1e-12 and out[-1].node == s.node:
            out[-1] = Segment(out[-1].left, s.right, s.node)
        else:
            out.append(s)
    return tuple(out)


def split(segments: tuple[Segment, ...], x: float) -> tuple[tuple[Segment, ...], tuple[Segment, ...]]:
    left: list[Segment] = []
    right: list[Segment] = []
    for s in segments:
        if s.right <= x:
            left.append(s)
        elif s.left >= x:
            right.append(s)
        else:
            left.append(Segment(s.left, x, s.node))
            right.append(Segment(x, s.right, s.node))
    return tuple(left), tuple(right)


def total_length(segments: tuple[Segment, ...]) -> float:
    return sum(s.right - s.left for s in segments)
