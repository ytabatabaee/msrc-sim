from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping
import csv
import json

import numpy as np
import yaml

from .analytic import TOPOLOGY_NAMES
from .experiments import _rearrangement_from_config, _tree_from_config
from .history_io import load_frozen_history
from .rearrangement import Rearrangement
from .spatial_statistics import (
    RearrangementInterval,
    classify_region,
    generate_locus_positions,
    sliding_windows,
    spatial_summary,
)
from .structured_coalescent import simulate_genealogy, simulate_msc_genealogy
from .population_process import PopulationFrequencyHistory, resolved_population_process, simulate_population_history


LOCUS_FIELDS = [
    "locus_id", "position", "region", "inside_rearrangement", "local_model",
    "local_effective_cross_arrangement_fraction", "topology_index", "topology",
    "terminal_pattern", "num_switch_events", "mean_or_first_coalescence_time",
]

WINDOW_FIELDS = [
    "window_id", "start_position", "end_position", "center_position", "n_loci",
    "n1", "n2", "n3", "q1", "q2", "q3", "dominant_topology",
    "nearest_msc_topology", "distance_to_nearest_msc_arm", "off_arm_difference",
    "off_arm_p_value", "fraction_inside_rearrangement",
]


def validate_spatial_config(config: Mapping[str, Any]) -> None:
    if config.get("mode") != "spatial":
        raise ValueError("Spatial simulation requires mode: spatial")
    for section in ("species_tree", "rearrangement", "recombination", "genome"):
        if section not in config:
            raise ValueError(f"Missing spatial configuration section: {section}")
    genome = config["genome"]
    interval = genome.get("rearrangement_interval")
    if not interval:
        raise ValueError("genome.rearrangement_interval is required")
    length = int(genome["length"])
    start = int(interval["start"])
    end = int(interval["end"])
    if not (0 <= start < end <= length):
        raise ValueError("Require 0 <= start < end <= genome.length")
    loci = genome.get("loci", {})
    if int(loci.get("count", 0)) <= 0:
        raise ValueError("genome.loci.count must be positive")
    if loci.get("placement", "evenly_spaced") not in {"evenly_spaced", "uniform_random"}:
        raise ValueError("genome.loci.placement must be evenly_spaced or uniform_random")
    if genome.get("inside_model", {}).get("type", "msrc") != "msrc":
        raise ValueError("v0.7.0 inside_model.type must be msrc")
    if genome.get("outside_model", {}).get("type", "msc") != "msc":
        raise ValueError("v0.7.0 outside_model.type must be msc")


def _terminal_pattern(sampled: Mapping[str, int], taxa: tuple[str, ...]) -> str:
    return "".join(str(int(sampled[t])) for t in taxa)


def _load_or_simulate_history(config: Mapping[str, Any], rng: np.random.Generator) -> tuple[PopulationFrequencyHistory, dict[str, int], dict[str, Any], Rearrangement]:
    history_path = config.get("history", {}).get("frozen_history")
    rearrangement = _rearrangement_from_config(config)
    if history_path:
        _, history, sampled, metadata = load_frozen_history(history_path)
        return history, sampled, {"frozen_history": str(history_path), **metadata}, rearrangement
    tree = _tree_from_config(config)
    history = simulate_population_history(tree, rearrangement, rng, dict(config))
    sampled = {t: int(rng.random() < history.terminal_frequency(t)) for t in tree.taxa}
    return history, sampled, {"frozen_history": None, "population_process": resolved_population_process(dict(config))}, rearrangement


def _write_frequency_history(out: Path, history: PopulationFrequencyHistory) -> None:
    rows = [asdict(r) for r in history.records]
    if not rows:
        return
    with (out / "frequency_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def simulate_spatial(config: Mapping[str, Any]) -> Path:
    from .linked_spatial import linked_spatial_enabled, simulate_linked_spatial
    if linked_spatial_enabled(config):
        return simulate_linked_spatial(config)
    validate_spatial_config(config)
    rng = np.random.default_rng(int(config.get("seed", 1)))
    tree = _tree_from_config(config)
    history, sampled, history_metadata, rearrangement = _load_or_simulate_history(config, rng)
    rec = config["recombination"]
    base_rate = float(rec.get("baseline_rate", rec.get("rate")))
    default_fraction = float(rec.get("effective_cross_arrangement_fraction", rec.get("suppression_factor", 1.0)))
    genome = config["genome"]
    loci_cfg = genome.get("loci", {})
    positions = generate_locus_positions(
        int(genome["length"]),
        int(loci_cfg["count"]),
        loci_cfg.get("placement", "evenly_spaced"),
        rng,
    )
    interval = RearrangementInterval(
        int(genome["rearrangement_interval"]["start"]),
        int(genome["rearrangement_interval"]["end"]),
        str(genome["rearrangement_interval"].get("id", rearrangement.rearrangement_id)),
    )
    inside_fraction = float(genome.get("inside_model", {}).get("effective_cross_arrangement_fraction", default_fraction))
    summary_cfg = config.get("spatial_summary", {})
    terminal = _terminal_pattern(sampled, tree.taxa)
    record_gene_trees = bool(config.get("output", {}).get("record_gene_trees", True))
    locus_rows: list[dict[str, Any]] = []
    gene_tree_rows: list[dict[str, Any]] = []

    for locus_id, position in enumerate(positions):
        region, inside = classify_region(int(position), interval)
        if inside:
            result = simulate_genealogy(locus_id, tree, history, sampled, base_rate, inside_fraction, rng, record_events=True)
            local_model = "msrc"
            local_fraction = inside_fraction
            switches = sum(e.event_type == "switch" for e in result.events)
        else:
            result = simulate_msc_genealogy(locus_id, tree, rng, record_events=False)
            local_model = "msc"
            local_fraction = ""
            switches = 0
        mean_time = float(np.mean(result.coalescence_times)) if result.coalescence_times else float("nan")
        locus_rows.append({
            "locus_id": locus_id,
            "position": int(position),
            "region": region,
            "inside_rearrangement": inside,
            "local_model": local_model,
            "local_effective_cross_arrangement_fraction": local_fraction,
            "topology_index": int(result.topology_index),
            "topology": TOPOLOGY_NAMES[int(result.topology_index)],
            "terminal_pattern": terminal,
            "num_switch_events": int(switches),
            "mean_or_first_coalescence_time": mean_time,
        })
        if record_gene_trees:
            gene_tree_rows.append({"locus_id": locus_id, "position": int(position), "newick": result.newick})

    windows = sliding_windows(
        locus_rows,
        int(summary_cfg.get("window_loci", 50)),
        int(summary_cfg.get("step_loci", 10)),
    )
    summary = spatial_summary(
        locus_rows,
        windows,
        interval,
        int(summary_cfg.get("breakpoint_bandwidth_loci", 50)),
    )
    summary.update({
        "version": "0.7.0",
        "mode": "spatial",
        "seed": int(config.get("seed", 1)),
        "genome_length": int(genome["length"]),
        "rearrangement_interval": asdict(interval),
        "terminal_pattern": terminal,
        "sampled_arrangements": sampled,
        "history_metadata": history_metadata,
        "population_process": resolved_population_process(dict(config)),
        "topology_names": TOPOLOGY_NAMES,
    })

    out = Path(config.get("output", {}).get("directory", "spatial_output"))
    out.mkdir(parents=True, exist_ok=True)
    if config.get("output", {}).get("record_resolved_config", True):
        with (out / "config.resolved.yaml").open("w") as handle:
            yaml.safe_dump(dict(config), handle, sort_keys=False)
    if config.get("output", {}).get("record_frequency_history", True):
        _write_frequency_history(out, history)
    with (out / "sampled_arrangements.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", "arrangement", "terminal_frequency_A1"])
        for taxon in tree.taxa:
            writer.writerow([taxon, sampled[taxon], history.terminal_frequency(taxon)])
    with (out / "spatial_loci.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOCUS_FIELDS)
        writer.writeheader()
        writer.writerows(locus_rows)
    if record_gene_trees:
        with (out / "spatial_gene_trees.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["locus_id", "position", "newick"], delimiter="\t")
            writer.writeheader()
            writer.writerows(gene_tree_rows)
    with (out / "spatial_windows.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WINDOW_FIELDS)
        writer.writeheader()
        writer.writerows({k: row[k] for k in WINDOW_FIELDS} for row in windows)
    with (out / "spatial_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    return out
