"""Reproducible, compact validation tables for the neutral diffusion.

This is deliberately a script rather than a CI test: the replicate counts can
be increased for publication-quality Monte Carlo uncertainty.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy
from msrcsim.wf_diffusion import _step


def _tree(ne: int, root_extension: float = 20.0) -> SpeciesTree:
    return SpeciesTree("((A:5,B:5)AB:5,(C:5,D:5)CD:5)ROOT;", ne, root_extension)


def _config(process: str, max_step_cu: float, tolerance: float) -> dict:
    cfg = {"population_process": {"model": process}}
    if process == "wf_diffusion":
        cfg["population_process"]["diffusion"] = {
            "max_step_cu": max_step_cu,
            "tolerance": tolerance,
        }
    return cfg


def _one_history(ne: int, process: str, rng: np.random.Generator,
                 max_step_cu: float, tolerance: float):
    tree = _tree(ne)
    rr = Rearrangement("inv", "inversion", "ROOT", 2, max(1, ne // 2), 0.0)
    history = simulate_population_history(tree, rr, rng, _config(process, max_step_cu, tolerance))
    p = float(history.frequency_at("ROOT", tree.branches["ROOT"].younger_age))
    sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
    result = simulate_genealogy(0, tree, history, sampled, 0.01, 0.1, rng, record_events=False)
    return tree, history, p, result


def run(output: Path, reps: int, seed: int, ne_values: tuple[int, ...],
        tolerances: tuple[float, ...], steps: tuple[float, ...]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for ne in ne_values:
        for process in ("wright_fisher", "moran"):
            rng = np.random.default_rng(seed + ne + (0 if process == "wright_fisher" else 100000))
            values = []
            tops = []
            for _ in range(reps):
                _tree_obj, _history, p, result = _one_history(ne, process, rng, 0.01, 0.001)
                values.append(p)
                tops.append(result.topology_index)
            for top in range(3):
                rows.append({"table": "matched", "process": process, "Ne": ne, "tolerance": "",
                             "max_step_cu": "", "replicates": reps,
                             "mean_terminal_frequency": np.mean(values),
                             "variance_terminal_frequency": np.var(values, ddof=1),
                             "fixation_probability": np.mean(np.asarray(values) == 1.0),
                             "loss_probability": np.mean(np.asarray(values) == 0.0),
                             "segregating_probability": np.mean((np.asarray(values) > 0) & (np.asarray(values) < 1)),
                             "quartet_topology": top, "quartet_probability": np.mean(np.asarray(tops) == top),
                             "mean_history_records": ""})
        for tolerance in tolerances:
            for step in steps:
                rng = np.random.default_rng(seed + ne * 1000 + round(tolerance * 1e6) + round(step * 1e6))
                values = []
                tops = []
                records = []
                for _ in range(reps):
                    _tree_obj, history, p, result = _one_history(ne, "wf_diffusion", rng, step, tolerance)
                    values.append(p)
                    tops.append(result.topology_index)
                    records.append(len(history.records))
                for top in range(3):
                    rows.append({"table": "convergence", "process": "wf_diffusion", "Ne": ne,
                                 "tolerance": tolerance, "max_step_cu": step, "replicates": reps,
                                 "mean_terminal_frequency": np.mean(values),
                                 "variance_terminal_frequency": np.var(values, ddof=1),
                                 "fixation_probability": np.mean(np.asarray(values) == 1.0),
                                 "loss_probability": np.mean(np.asarray(values) == 0.0),
                                 "segregating_probability": np.mean((np.asarray(values) > 0) & (np.asarray(values) < 1)),
                                 "quartet_topology": top, "quartet_probability": np.mean(np.asarray(tops) == top),
                                 "mean_history_records": np.mean(records)})
    fields = list(rows[0])
    with (output / "wf_diffusion_validation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    boundary_rows = []
    for p0 in (0.0005, 0.1, 0.5, 0.9, 0.9995):
        for tolerance in tolerances:
            for step in steps:
                rng = np.random.default_rng(seed + round(p0 * 1e6) + round(tolerance * 1e6) + round(step * 1e6))
                terminal = []
                for _ in range(reps):
                    p = p0
                    for _ in range(1000):
                        if p == 0.0 or p == 1.0:
                            break
                        p = _step(p, step, rng, tolerance)
                    terminal.append(p)
                values = np.asarray(terminal)
                boundary_rows.append({
                    "p0": p0, "tolerance": tolerance, "max_step_cu": step,
                    "replicates": reps, "mean_terminal_frequency": values.mean(),
                    "fixation_probability": np.mean(values == 1.0),
                    "loss_probability": np.mean(values == 0.0),
                    "absorption_probability": np.mean((values == 0.0) | (values == 1.0)),
                })
    with (output / "wf_diffusion_boundary_validation.csv").open("w", newline="") as handle:
        fields = list(boundary_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(boundary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("wf_diffusion_validation"))
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261011)
    args = parser.parse_args()
    run(args.output, args.replicates, args.seed, (20, 100, 500), (0.01, 0.003, 0.001), (0.05, 0.02, 0.01))


if __name__ == "__main__":
    main()
