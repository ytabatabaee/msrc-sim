from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping
import csv
import json

import numpy as np
import yaml

from .analytic import TOPOLOGY_NAMES
from .experiments import _tree_from_config
from .genomic import GenomicInterval, interval_at, ordered_windows, sample_metadata_rows, validate_contiguous_block_ids
from .spatial import _load_or_simulate_history
from .structured_coalescent import simulate_genealogy, simulate_msc_genealogy


SPATIAL_GENEALOGY_FIELDS = [
    "window_id", "chrom", "start", "end", "midpoint", "block_id",
    "topology", "topology_index", "is_rearranged", "rearrangement_id",
    "local_model", "newick",
]
LINKAGE_DIAGNOSTIC_FIELDS = [
    "region", "num_windows", "num_genealogy_blocks", "mean_block_length_bp",
    "median_block_length_bp", "breakpoint_density_per_bp",
    "observed_inside_outside_rate_ratio", "expected_kappa",
]


@dataclass(frozen=True)
class GenealogyBlock:
    block_id: int
    chrom: str
    start: float
    end: float
    is_rearranged: bool
    rearrangement_id: str
    topology_index: int
    topology: str
    newick: str
    local_model: str


def linked_spatial_enabled(config: Mapping[str, Any]) -> bool:
    return bool((config.get("linked_spatial", {}) or {}).get("enabled", False))


def _linked_cfg(config: Mapping[str, Any]) -> Mapping[str, Any]:
    return config.get("linked_spatial", {}) or {}


def _chrom_length(config: Mapping[str, Any]) -> tuple[str, float]:
    cfg = _linked_cfg(config)
    genome = config.get("genome", {}) or {}
    chrom = str(cfg.get("chrom", genome.get("chrom", "chr1")))
    length = float(cfg.get("length_bp", genome.get("length_bp", genome.get("length", 0))))
    if length <= 0.0:
        raise ValueError("linked spatial mode requires linked_spatial.length_bp or genome.length")
    return chrom, length


def _intervals(config: Mapping[str, Any], chrom: str, length: float) -> list[GenomicInterval]:
    cfg = _linked_cfg(config)
    raw = cfg.get("rearrangements")
    if raw is None:
        genome = config.get("genome", {}) or {}
        if "rearrangement_intervals" in genome:
            raw = genome["rearrangement_intervals"]
        elif "rearrangement_interval" in genome:
            raw = [genome["rearrangement_interval"]]
        else:
            raw = []
    intervals = [
        GenomicInterval(
            chrom=str(item.get("chrom", chrom)),
            start=float(item["start"]),
            end=float(item["end"]),
            interval_id=str(item.get("id", item.get("interval_id", f"rearrangement_{i}"))),
        )
        for i, item in enumerate(raw)
    ]
    for interval in intervals:
        if interval.chrom != chrom:
            raise ValueError("v0.8.0 linked spatial mode supports one chromosome per run")
        if interval.end > length:
            raise ValueError("rearrangement interval cannot exceed chromosome length")
    intervals.sort(key=lambda item: item.start)
    for left, right in zip(intervals, intervals[1:]):
        if right.start < left.end:
            raise ValueError("rearrangement intervals must not overlap")
    return intervals


def _breakpoint_rate(config: Mapping[str, Any]) -> float:
    cfg = _linked_cfg(config)
    rec = config.get("recombination", {}) or {}
    rate = float(cfg.get("breakpoint_rate_per_bp", rec.get("genealogy_breakpoint_rate_per_bp", 0.0)))
    if rate < 0.0:
        raise ValueError("breakpoint_rate_per_bp must be nonnegative")
    return rate


def _kappa(config: Mapping[str, Any]) -> float:
    cfg = _linked_cfg(config)
    value = float(cfg.get("kappa", cfg.get("suppression_factor", 1.0)))
    if not (0.0 <= value <= 1.0):
        raise ValueError("linked_spatial.kappa must be between 0 and 1")
    return value


def _forced_boundaries(length: float, intervals: list[GenomicInterval]) -> list[float]:
    points = {0.0, float(length)}
    for interval in intervals:
        points.add(float(interval.start))
        points.add(float(interval.end))
    return sorted(points)


def generate_genealogy_breakpoints(
    length: float,
    intervals: list[GenomicInterval],
    baseline_rate_per_bp: float,
    kappa: float,
    rng: np.random.Generator,
) -> list[float]:
    points = _forced_boundaries(length, intervals)
    extra: list[float] = []
    for start, end in zip(points, points[1:]):
        mid = (start + end) / 2.0
        local_rate = baseline_rate_per_bp * (kappa if interval_at(mid, intervals) else 1.0)
        n = int(rng.poisson(local_rate * (end - start))) if local_rate > 0.0 else 0
        if n:
            extra.extend(float(x) for x in rng.uniform(start, end, size=n))
    return sorted(set(points + extra))


def block_diagnostics(blocks: list[GenealogyBlock], rows: list[Mapping[str, Any]], kappa: float | None = None) -> list[dict[str, Any]]:
    observed_ids = {int(row["block_id"]) for row in rows}
    ordered_rows = sorted(rows, key=lambda row: int(row["window_id"]))
    row_counts = {
        "inside": sum(bool(row["is_rearranged"]) for row in rows),
        "outside": sum(not bool(row["is_rearranged"]) for row in rows),
    }
    out = []
    by_region: dict[str, dict[str, Any]] = {}
    for region, rearranged in (("inside", True), ("outside", False)):
        chunk = [block for block in blocks if block.is_rearranged is rearranged and int(block.block_id) in observed_ids]
        total_length = sum(block.end - block.start for block in chunk)
        breakpoints = 0
        run_blocks: set[int] = set()
        for row in ordered_rows:
            if bool(row["is_rearranged"]) is rearranged:
                run_blocks.add(int(row["block_id"]))
            elif run_blocks:
                breakpoints += max(0, len(run_blocks) - 1)
                run_blocks = set()
        if run_blocks:
            breakpoints += max(0, len(run_blocks) - 1)
        lengths = [block.end - block.start for block in chunk]
        by_region[region] = {
            "region": region,
            "num_windows": int(row_counts[region]),
            "num_genealogy_blocks": int(len(chunk)),
            "mean_block_length_bp": float(total_length / len(chunk)) if chunk else float("nan"),
            "median_block_length_bp": float(np.median(lengths)) if lengths else float("nan"),
            "breakpoint_density_per_bp": float(breakpoints / total_length) if total_length > 0.0 else float("nan"),
        }
    inside_density = float(by_region["inside"]["breakpoint_density_per_bp"])
    outside_density = float(by_region["outside"]["breakpoint_density_per_bp"])
    ratio = inside_density / outside_density if np.isfinite(inside_density) and np.isfinite(outside_density) and outside_density > 0.0 else float("nan")
    for row in by_region.values():
        row["observed_inside_outside_rate_ratio"] = ratio
        row["expected_kappa"] = "" if kappa is None else float(kappa)
        out.append(row)
    return out


def _windows_from_config(config: Mapping[str, Any], chrom: str, length: float):
    cfg = _linked_cfg(config)
    genome = config.get("genome", {}) or {}
    windows_cfg = cfg.get("windows", genome.get("windows", {})) or {}
    count = windows_cfg.get("count", cfg.get("window_count"))
    size = windows_cfg.get("size_bp", cfg.get("window_size_bp"))
    if count is None and size is None:
        loci = genome.get("loci", {}) or {}
        count = int(loci.get("count", 0)) or None
    return ordered_windows(chrom, length, window_size=None if size is None else float(size), count=None if count is None else int(count))


def simulate_linked_spatial(config: Mapping[str, Any]) -> Path:
    if config.get("mode") != "spatial":
        raise ValueError("linked spatial simulation requires mode: spatial")
    rng = np.random.default_rng(int(config.get("seed", 1)))
    tree = _tree_from_config(config)
    history, sampled, history_metadata, _rearrangement = _load_or_simulate_history(config, rng)
    chrom, length = _chrom_length(config)
    intervals = _intervals(config, chrom, length)
    windows = _windows_from_config(config, chrom, length)
    breakpoints = generate_genealogy_breakpoints(length, intervals, _breakpoint_rate(config), _kappa(config), rng)
    rec = config.get("recombination", {}) or {}
    base_rate = float(rec.get("baseline_rate", rec.get("rate", 0.0)))
    inside_fraction = float((config.get("genome", {}) or {}).get("inside_model", {}).get(
        "effective_cross_arrangement_fraction",
        rec.get("effective_cross_arrangement_fraction", rec.get("suppression_factor", 1.0)),
    ))

    raw_blocks: list[GenealogyBlock] = []
    for raw_id, (start, end) in enumerate(zip(breakpoints, breakpoints[1:])):
        mid = (start + end) / 2.0
        interval = interval_at(mid, intervals)
        if interval is None:
            result = simulate_msc_genealogy(raw_id, tree, rng, record_events=False)
            local_model = "msc"
            rearranged = False
            rearrangement_id = ""
        else:
            result = simulate_genealogy(raw_id, tree, history, sampled, base_rate, inside_fraction, rng, record_events=False)
            local_model = "msrc"
            rearranged = True
            rearrangement_id = interval.interval_id
        top = int(result.topology_index)
        raw_blocks.append(GenealogyBlock(raw_id, chrom, float(start), float(end), rearranged, rearrangement_id, top, TOPOLOGY_NAMES[top], result.newick, local_model))

    block_starts = np.asarray([b.start for b in raw_blocks], dtype=float)
    used: dict[int, int] = {}
    observed_blocks: dict[int, GenealogyBlock] = {}
    rows: list[dict[str, Any]] = []
    for window in windows:
        idx = int(np.searchsorted(block_starts, window.midpoint, side="right") - 1)
        idx = min(max(idx, 0), len(raw_blocks) - 1)
        block = raw_blocks[idx]
        if idx not in used:
            used[idx] = len(used)
        block_id = used[idx]
        observed_blocks[block_id] = GenealogyBlock(block_id, block.chrom, block.start, block.end, block.is_rearranged, block.rearrangement_id, block.topology_index, block.topology, block.newick, block.local_model)
        rows.append({
            "window_id": int(window.window_id),
            "chrom": window.chrom,
            "start": float(window.start),
            "end": float(window.end),
            "midpoint": float(window.midpoint),
            "block_id": int(block_id),
            "topology": block.topology,
            "topology_index": int(block.topology_index),
            "is_rearranged": bool(block.is_rearranged),
            "rearrangement_id": block.rearrangement_id,
            "local_model": block.local_model,
            "newick": block.newick,
        })
    validate_contiguous_block_ids(rows)
    diagnostics = block_diagnostics(list(observed_blocks.values()), rows, _kappa(config))

    out = Path((config.get("output", {}) or {}).get("directory", "linked_spatial_output"))
    out.mkdir(parents=True, exist_ok=True)
    if (config.get("output", {}) or {}).get("record_resolved_config", True):
        with (out / "config.resolved.yaml").open("w") as handle:
            yaml.safe_dump(dict(config), handle, sort_keys=False)
    with (out / "spatial_genealogies.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SPATIAL_GENEALOGY_FIELDS)
        writer.writeheader()
        writer.writerows({k: row[k] for k in SPATIAL_GENEALOGY_FIELDS} for row in rows)
    with (out / "spatial_gene_trees.nwk").open("w") as handle:
        for row in rows:
            handle.write(str(row["newick"]) + "\n")
    with (out / "spatial_linkage_diagnostics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LINKAGE_DIAGNOSTIC_FIELDS)
        writer.writeheader()
        writer.writerows(diagnostics)
    metadata = sample_metadata_rows(config, tree.taxa, sampled)
    if metadata:
        fields = sorted({key for row in metadata for key in row})
        with (out / "sample_metadata.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(metadata)
    summary = {
        "version": "0.8.3",
        "mode": "spatial",
        "linked_loci_model": True,
        "spatial_model_note": (
            "This is a piecewise-correlated genealogy process, not a full ARG. "
            "Dense windows within one block are correlated observations, not independent replicates."
        ),
        "seed": int(config.get("seed", 1)),
        "chrom": chrom,
        "chromosome_length_bp": float(length),
        "num_windows": len(rows),
        "num_blocks_observed": len(set(int(row["block_id"]) for row in rows)),
        "kappa": _kappa(config),
        "breakpoint_rate_per_bp": _breakpoint_rate(config),
        "rearrangements": [asdict(interval) for interval in intervals],
        "history_metadata": history_metadata,
        "linkage_diagnostics": diagnostics,
    }
    with (out / "linked_spatial_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    return out
