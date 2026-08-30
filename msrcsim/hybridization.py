from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping
import csv
import json

import numpy as np
import yaml

from .analytic import TOPOLOGY_NAMES
from .ancestry_tracts import AncestryTract, ancestry_label, introgressed_fraction_by_length, simulate_ancestry_tracts
from .model_fitting import msc_probabilities, network_probabilities, off_arm_statistics


TOPOLOGY_LABEL_TO_INDEX = {label: i for i, label in enumerate(TOPOLOGY_NAMES)}


@dataclass(frozen=True)
class HybridizationLocus:
    locus_id: int
    position_bp: float
    ancestry_state: int
    ancestry_label: str
    topology_index: int
    topology_label: str
    conditional_q1: float
    conditional_q2: float
    conditional_q3: float
    tract_id: int
    distance_to_left_tract_boundary: float
    distance_to_right_tract_boundary: float
    model: str = "pulse_hybridization"


def topology_index(value: str | int) -> int:
    if isinstance(value, str):
        if value not in TOPOLOGY_LABEL_TO_INDEX:
            raise ValueError(f"topology must be one of {', '.join(TOPOLOGY_NAMES)}")
        return int(TOPOLOGY_LABEL_TO_INDEX[value])
    value = int(value)
    if value not in (0, 1, 2):
        raise ValueError("topology index must be 0, 1, or 2")
    return value


def _validate_config(config: Mapping[str, Any]) -> None:
    if config.get("mode") != "pulse_hybridization":
        raise ValueError("Hybridization simulation requires mode: pulse_hybridization")
    for section in ("chromosome", "hybridization", "recombination", "windows"):
        if section not in config:
            raise ValueError(f"Missing hybridization configuration section: {section}")
    chrom = config["chromosome"]
    if float(chrom["length_bp"]) <= 0.0:
        raise ValueError("chromosome.length_bp must be positive")
    if int(chrom["num_loci"]) <= 0:
        raise ValueError("chromosome.num_loci must be positive")
    hyb = config["hybridization"]
    gamma = float(hyb["gamma"])
    if not (0.0 <= gamma <= 1.0):
        raise ValueError("hybridization.gamma must be between 0 and 1")
    if float(hyb["generations_since_pulse"]) < 0.0:
        raise ValueError("hybridization.generations_since_pulse must be nonnegative")
    for parent in ("major", "introgressed"):
        topology_index(hyb[parent]["topology"])
        if float(hyb[parent]["internal_branch_length"]) < 0.0:
            raise ValueError(f"hybridization.{parent}.internal_branch_length must be nonnegative")
    if float(config["recombination"]["rate_per_bp_per_generation"]) < 0.0:
        raise ValueError("recombination.rate_per_bp_per_generation must be nonnegative")
    if int(config["windows"]["loci_per_window"]) <= 0 or int(config["windows"]["step_loci"]) <= 0:
        raise ValueError("windows.loci_per_window and windows.step_loci must be positive")


def generate_hybridization_locus_positions(
    chromosome: Mapping[str, Any],
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[int]]:
    length = float(chromosome["length_bp"])
    count = int(chromosome["num_loci"])
    cfg = chromosome.get("locus_positions", {"mode": "evenly_spaced"})
    mode = cfg.get("mode", "evenly_spaced")
    if mode == "evenly_spaced":
        positions = np.linspace(0.0, length, count, endpoint=False)
        locus_ids = list(range(count))
    elif mode == "random_uniform":
        positions = np.sort(rng.uniform(0.0, length, size=count))
        locus_ids = list(range(count))
    elif mode == "file":
        path = Path(cfg["path"])
        column = str(cfg.get("column", "position_bp"))
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
        rows: list[tuple[int, float]] = []
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"position file must contain column {column!r}")
            for row_number, row in enumerate(reader):
                locus_id = int(row["locus_id"]) if "locus_id" in row and row["locus_id"] != "" else row_number
                rows.append((locus_id, float(row[column])))
        if len(rows) != count:
            raise ValueError("chromosome.num_loci must match the number of rows in the position file")
        rows.sort(key=lambda item: item[1])
        locus_ids = [row[0] for row in rows]
        positions = np.asarray([row[1] for row in rows], dtype=float)
    else:
        raise ValueError("locus_positions.mode must be evenly_spaced, random_uniform, or file")
    if np.any(positions < 0.0) or np.any(positions > length):
        raise ValueError("locus positions must lie within the chromosome")
    if len(np.unique(positions)) != count:
        raise ValueError("locus positions must be unique")
    return positions.astype(float), locus_ids


def _tract_for_position(tracts: list[AncestryTract], position: float) -> AncestryTract:
    starts = np.asarray([t.start_bp for t in tracts], dtype=float)
    idx = int(np.searchsorted(starts, position, side="right") - 1)
    idx = min(max(idx, 0), len(tracts) - 1)
    return tracts[idx]


def sample_hybridization_loci(
    positions: np.ndarray,
    locus_ids: list[int],
    tracts: list[AncestryTract],
    q_major: np.ndarray,
    q_introgressed: np.ndarray,
    rng: np.random.Generator,
) -> list[HybridizationLocus]:
    loci: list[HybridizationLocus] = []
    for position, locus_id in zip(positions, locus_ids):
        tract = _tract_for_position(tracts, float(position))
        q = q_major if tract.ancestry_state == 0 else q_introgressed
        top = int(rng.choice(3, p=q))
        loci.append(HybridizationLocus(
            locus_id=int(locus_id),
            position_bp=float(position),
            ancestry_state=int(tract.ancestry_state),
            ancestry_label=ancestry_label(tract.ancestry_state),
            topology_index=top,
            topology_label=TOPOLOGY_NAMES[top],
            conditional_q1=float(q[0]),
            conditional_q2=float(q[1]),
            conditional_q3=float(q[2]),
            tract_id=int(tract.tract_id),
            distance_to_left_tract_boundary=float(position - tract.start_bp),
            distance_to_right_tract_boundary=float(tract.end_bp - position),
        ))
    return loci


def _counts(topologies: list[int]) -> np.ndarray:
    return np.bincount(np.asarray(topologies, dtype=int), minlength=3)[:3]


def hybridization_windows(loci: list[HybridizationLocus], loci_per_window: int, step_loci: int) -> list[dict[str, Any]]:
    if loci_per_window <= 0 or step_loci <= 0:
        raise ValueError("loci_per_window and step_loci must be positive")
    ordered = sorted(loci, key=lambda locus: locus.position_bp)
    out: list[dict[str, Any]] = []
    for window_id, start in enumerate(range(0, len(ordered) - loci_per_window + 1, step_loci)):
        chunk = ordered[start:start + loci_per_window]
        counts = _counts([l.topology_index for l in chunk])
        q = counts / counts.sum()
        model = off_arm_statistics(counts)
        first = float(chunk[0].position_bp)
        last = float(chunk[-1].position_bp)
        out.append({
            "window_id": window_id,
            "start_bp": first,
            "end_bp": last,
            "center_bp": (first + last) / 2.0,
            "num_loci": len(chunk),
            "introgressed_fraction": float(np.mean([l.ancestry_state == 1 for l in chunk])),
            "n1": int(counts[0]), "n2": int(counts[1]), "n3": int(counts[2]),
            "q1": float(q[0]), "q2": float(q[1]), "q3": float(q[2]),
            "dominant_topology": int(np.argmax(q)),
            "distance_to_nearest_msc_arm": float(model["distance_to_nearest_msc_arm"]),
            "off_arm_difference": float(model["off_arm_difference"]),
            "model": "pulse_hybridization",
        })
    return out


def lag_same_probability(loci: list[HybridizationLocus], lags_bp: list[float]) -> list[dict[str, Any]]:
    ordered = sorted(loci, key=lambda locus: locus.position_bp)
    positions = np.asarray([l.position_bp for l in ordered], dtype=float)
    topologies = np.asarray([l.topology_index for l in ordered], dtype=int)
    ancestry = np.asarray([l.ancestry_state for l in ordered], dtype=int)
    out: list[dict[str, Any]] = []
    for lag in lags_bp:
        same_top = []
        same_anc = []
        for i, pos in enumerate(positions):
            j = int(np.searchsorted(positions, pos + float(lag), side="left"))
            if j < len(positions):
                same_top.append(topologies[i] == topologies[j])
                same_anc.append(ancestry[i] == ancestry[j])
        out.append({
            "lag_bp": float(lag),
            "n_pairs": int(len(same_top)),
            "prob_same_topology": float(np.mean(same_top)) if same_top else float("nan"),
            "prob_same_ancestry": float(np.mean(same_anc)) if same_anc else float("nan"),
        })
    return out


def hybridization_summary(
    config: Mapping[str, Any],
    tracts: list[AncestryTract],
    loci: list[HybridizationLocus],
    q_major: np.ndarray,
    q_introgressed: np.ndarray,
    expected_marginal_q: np.ndarray,
    lag_stats: list[dict[str, Any]],
) -> dict[str, Any]:
    hyb = config["hybridization"]
    counts = _counts([l.topology_index for l in loci])
    freqs = counts / counts.sum()
    intro_lengths = [t.length_bp for t in tracts if t.ancestry_state == 1]
    major_lengths = [t.length_bp for t in tracts if t.ancestry_state == 0]
    return {
        "version": "0.7.0",
        "mode": "pulse_hybridization",
        "linked_ancestry_mosaic": True,
        "linked_loci_model": False,
        "spatial_model_note": (
            "The ancestry mosaic is linked along the chromosome, but this release does not simulate "
            "a full ARG or full linked multispecies-network coalescent. Conditional on ancestry state, "
            "local quartet topologies are sampled independently from the corresponding MSC distribution."
        ),
        "seed": int(config.get("seed", 1)),
        "chromosome_length_bp": float(config["chromosome"]["length_bp"]),
        "num_loci": int(config["chromosome"]["num_loci"]),
        "gamma": float(hyb["gamma"]),
        "generations_since_pulse": float(hyb["generations_since_pulse"]),
        "recombination_rate_per_bp_per_generation": float(config["recombination"]["rate_per_bp_per_generation"]),
        "major_topology": int(topology_index(hyb["major"]["topology"])),
        "major_topology_label": TOPOLOGY_NAMES[topology_index(hyb["major"]["topology"])],
        "major_internal_branch_length": float(hyb["major"]["internal_branch_length"]),
        "introgressed_topology": int(topology_index(hyb["introgressed"]["topology"])),
        "introgressed_topology_label": TOPOLOGY_NAMES[topology_index(hyb["introgressed"]["topology"])],
        "introgressed_internal_branch_length": float(hyb["introgressed"]["internal_branch_length"]),
        "donor": hyb.get("donor"),
        "recipient": hyb.get("recipient"),
        "q_major": [float(x) for x in q_major],
        "q_introgressed": [float(x) for x in q_introgressed],
        "expected_marginal_q": [float(x) for x in expected_marginal_q],
        "expected_marginal_q_formula": "network_probabilities from msrcsim.model_fitting",
        "observed_marginal_q": [float(x) for x in freqs],
        "realized_introgressed_fraction_by_length": float(introgressed_fraction_by_length(tracts)),
        "realized_introgressed_fraction_at_loci": float(np.mean([l.ancestry_state == 1 for l in loci])),
        "num_tracts": int(len(tracts)),
        "num_introgressed_tracts": int(len(intro_lengths)),
        "num_major_tracts": int(len(major_lengths)),
        "mean_introgressed_tract_length_bp": float(np.mean(intro_lengths)) if intro_lengths else None,
        "mean_major_tract_length_bp": float(np.mean(major_lengths)) if major_lengths else None,
        "median_introgressed_tract_length_bp": float(np.median(intro_lengths)) if intro_lengths else None,
        "median_major_tract_length_bp": float(np.median(major_lengths)) if major_lengths else None,
        "topology_counts": [int(x) for x in counts],
        "topology_frequencies": [float(x) for x in freqs],
        "num_topology_change_points": int(sum(a.topology_index != b.topology_index for a, b in zip(loci, loci[1:]))),
        "lag_statistics": lag_stats,
    }


TRACT_FIELDS = ["tract_id", "start_bp", "end_bp", "length_bp", "ancestry_state", "ancestry_label"]
LOCUS_FIELDS = [
    "locus_id", "position_bp", "model", "ancestry_state", "ancestry_label",
    "tract_id", "distance_to_left_tract_boundary", "distance_to_right_tract_boundary",
    "topology_index", "topology_label", "conditional_q1", "conditional_q2", "conditional_q3",
]
WINDOW_FIELDS = [
    "window_id", "start_bp", "end_bp", "center_bp", "num_loci", "introgressed_fraction",
    "n1", "n2", "n3", "q1", "q2", "q3", "dominant_topology",
    "distance_to_nearest_msc_arm", "off_arm_difference", "model",
]
LAG_FIELDS = ["lag_bp", "n_pairs", "prob_same_topology", "prob_same_ancestry"]


def _write_outputs(
    out: Path,
    config: Mapping[str, Any],
    tracts: list[AncestryTract],
    loci: list[HybridizationLocus],
    windows: list[dict[str, Any]],
    summary: dict[str, Any],
    lag_stats: list[dict[str, Any]],
) -> None:
    output_cfg = config.get("output", {})
    if output_cfg.get("record_resolved_config", True):
        with (out / "config.resolved.yaml").open("w") as handle:
            yaml.safe_dump(dict(config), handle, sort_keys=False)
    if output_cfg.get("record_tracts", True):
        with (out / "hybridization_tracts.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=TRACT_FIELDS)
            writer.writeheader()
            for tract in tracts:
                writer.writerow({
                    "tract_id": tract.tract_id,
                    "start_bp": tract.start_bp,
                    "end_bp": tract.end_bp,
                    "length_bp": tract.length_bp,
                    "ancestry_state": tract.ancestry_state,
                    "ancestry_label": tract.ancestry_label,
                })
    if output_cfg.get("record_loci", True):
        with (out / "hybridization_loci.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=LOCUS_FIELDS)
            writer.writeheader()
            writer.writerows({k: asdict(locus)[k] for k in LOCUS_FIELDS} for locus in loci)
    if output_cfg.get("record_windows", True):
        with (out / "hybridization_windows.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=WINDOW_FIELDS)
            writer.writeheader()
            writer.writerows({k: row[k] for k in WINDOW_FIELDS} for row in windows)
    with (out / "hybridization_spatial_autocorrelation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LAG_FIELDS)
        writer.writeheader()
        writer.writerows(lag_stats)
    with (out / "hybridization_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)


def simulate_hybridization(config: Mapping[str, Any]) -> Path:
    _validate_config(config)
    rng = np.random.default_rng(int(config.get("seed", 1)))
    hyb = config["hybridization"]
    major_topology = topology_index(hyb["major"]["topology"])
    intro_topology = topology_index(hyb["introgressed"]["topology"])
    gamma = float(hyb["gamma"])
    q_major = msc_probabilities(major_topology, float(hyb["major"]["internal_branch_length"]))
    q_intro = msc_probabilities(intro_topology, float(hyb["introgressed"]["internal_branch_length"]))
    expected = network_probabilities(
        major_topology,
        intro_topology,
        gamma,
        float(hyb["major"]["internal_branch_length"]),
        float(hyb["introgressed"]["internal_branch_length"]),
    ) if major_topology != intro_topology else (1.0 - gamma) * q_major + gamma * q_intro
    positions, locus_ids = generate_hybridization_locus_positions(config["chromosome"], rng)
    tracts = simulate_ancestry_tracts(
        float(config["chromosome"]["length_bp"]),
        gamma,
        float(hyb["generations_since_pulse"]),
        float(config["recombination"]["rate_per_bp_per_generation"]),
        rng,
    )
    loci = sample_hybridization_loci(positions, locus_ids, tracts, q_major, q_intro, rng)
    output_cfg = config.get("output", {})
    if output_cfg.get("record_windows", True) or output_cfg.get("make_plots", True):
        windows = hybridization_windows(loci, int(config["windows"]["loci_per_window"]), int(config["windows"]["step_loci"]))
    else:
        windows = []
    lag_cfg = config.get("spatial_statistics", {})
    lags = [float(x) for x in lag_cfg.get("lags_bp", [10000, 100000, 1000000, 5000000])]
    lag_stats = lag_same_probability(loci, lags)
    summary = hybridization_summary(config, tracts, loci, q_major, q_intro, expected, lag_stats)
    out = Path(config.get("output", {}).get("directory", "hybridization_output"))
    out.mkdir(parents=True, exist_ok=True)
    _write_outputs(out, config, tracts, loci, windows, summary, lag_stats)
    return out
