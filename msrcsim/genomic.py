from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class GenomicInterval:
    chrom: str
    start: float
    end: float
    interval_id: str

    def __post_init__(self) -> None:
        if self.start < 0.0 or self.end <= self.start:
            raise ValueError("genomic intervals require 0 <= start < end")

    def contains(self, position: float) -> bool:
        return self.start <= position < self.end


@dataclass(frozen=True)
class GenomicWindow:
    window_id: int
    chrom: str
    start: float
    end: float

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


def ordered_windows(chrom: str, length: float, *, window_size: float | None = None, count: int | None = None) -> list[GenomicWindow]:
    if length <= 0.0:
        raise ValueError("chromosome length must be positive")
    if count is None:
        if window_size is None or window_size <= 0.0:
            raise ValueError("provide a positive window_size or count")
        count = int((length + window_size - 1.0) // window_size)
    if count <= 0:
        raise ValueError("window count must be positive")
    if window_size is None:
        window_size = length / count
    out: list[GenomicWindow] = []
    for i in range(count):
        start = i * window_size
        end = min(length, (i + 1) * window_size)
        if start < length and end > start:
            out.append(GenomicWindow(i, chrom, float(start), float(end)))
    return out


def interval_at(position: float, intervals: list[GenomicInterval]) -> GenomicInterval | None:
    for interval in intervals:
        if interval.contains(position):
            return interval
    return None


def validate_contiguous_block_ids(rows: list[Mapping[str, Any]]) -> None:
    ids = [int(row["block_id"]) for row in rows]
    if not ids:
        return
    if ids != sorted(ids):
        raise ValueError("block_id values must be nondecreasing along ordered coordinates")
    unique = sorted(set(ids))
    if unique != list(range(unique[-1] + 1)):
        raise ValueError("block_id values must be contiguous from zero")


def sample_metadata_rows(config: Mapping[str, Any], taxa: tuple[str, ...], sampled_arrangements: Mapping[str, int]) -> list[dict[str, Any]]:
    sampling = config.get("sampling", {}) or {}
    rows = sampling.get("samples")
    if rows:
        out = []
        for i, row in enumerate(rows):
            item = dict(row)
            species = str(item.get("species", item.get("population", item.get("taxon", ""))))
            item.setdefault("sample", f"sample_{i}")
            item.setdefault("species", species)
            if species in sampled_arrangements:
                item.setdefault("arrangement", int(sampled_arrangements[species]))
            out.append(item)
        return out
    return [
        {"sample": str(taxon), "species": str(taxon), "population": str(taxon), "arrangement": int(sampled_arrangements[taxon])}
        for taxon in taxa
    ]
