"""Run the staged neutral WF-diffusion validation in this repository.

The defaults are a compact, resumable screening run. Increase --replicates
for publication-scale uncertainty; every result row records the replicate
count and numerical settings.
"""

from __future__ import annotations

import argparse
import csv
import json
import resource
import time
from pathlib import Path

import numpy as np

from msrcsim.arbitrary import tree_from_config
from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.structured_coalescent import simulate_genealogy
from msrcsim.wf_diffusion import _step
from scripts.validate_wf_diffusion import _config, _tree


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _se(values: np.ndarray) -> float:
    return float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else float("nan")


def analytical(out: Path, reps: int, seed: int) -> None:
    rows = []
    p_values = (0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
    u_values = (0.01, 0.1, 0.5, 1.0, 3.0)
    du = 0.002
    for p0 in p_values:
        for u in u_values:
            n_steps = int(round(u / du))
            rng = np.random.default_rng(seed + int(p0 * 1e6) + int(u * 1000))
            moment_reps = max(reps, 20000) if p0 in (0.001, 0.999) else (max(reps, 5000) if p0 in (0.01, 0.99) else max(reps, 1000))
            values = np.empty(moment_reps)
            het = np.empty(moment_reps)
            for i in range(moment_reps):
                p = p0
                for _ in range(n_steps):
                    p = _step(p, du, rng, 1e-5)
                values[i] = p
                het[i] = 2 * p * (1 - p)
            expected_het = 2 * p0 * (1 - p0) * np.exp(-u)
            rows.extend([
                {"p0": p0, "u": u, "statistic": "mean_frequency", "observed": values.mean(),
                 "expected": p0, "absolute_error": abs(values.mean() - p0), "mc_se": _se(values),
                 "replicates": moment_reps, "du": du, "tolerance": 1e-5},
                {"p0": p0, "u": u, "statistic": "heterozygosity", "observed": het.mean(),
                 "expected": expected_het, "absolute_error": abs(het.mean() - expected_het), "mc_se": _se(het),
                 "replicates": moment_reps, "du": du, "tolerance": 1e-5},
            ])

    for p0 in p_values:
        n = (20000 if p0 in (0.001, 0.999) else
             5000 if p0 in (0.01, 0.99) else max(reps, 1000))
        rng = np.random.default_rng(seed + 700000 + int(p0 * 1e6))
        terminal = np.empty(n)
        absorption_u = np.empty(n)
        for i in range(n):
            p = p0
            u = 0.0
            for _ in range(25000):
                if p == 0.0 or p == 1.0:
                    break
                p = _step(p, du, rng, 1e-5)
                u += du
            terminal[i] = p
            absorption_u[i] = u
        fixed = terminal == 1.0
        rows.append({"p0": p0, "u": "eventual", "statistic": "fixation_probability", "observed": fixed.mean(),
                     "expected": p0, "absolute_error": abs(fixed.mean() - p0),
                     "mc_se": np.sqrt(max(fixed.mean() * (1 - fixed.mean()), 1e-15) / n),
                     "replicates": n, "du": du, "tolerance": 1e-5})
        rows.append({"p0": p0, "u": "eventual", "statistic": "mean_absorption_u", "observed": absorption_u.mean(),
                     "expected": "", "absolute_error": "", "mc_se": _se(absorption_u),
                     "replicates": n, "du": du, "tolerance": 1e-5})
    _write(out / "results/analytical_moments.tsv", rows)


def convergence(out: Path, reps: int, seed: int) -> None:
    rows = []
    steps = (0.05, 0.02, 0.01, 0.005, 0.002)
    tolerances = (1e-2, 1e-3, 1e-4, 1e-5)
    for tol in tolerances:
        for step in steps:
            rng = np.random.default_rng(seed + round(tol * 1e6) + round(step * 1e6))
            values = []
            hetero = []
            absorbed = []
            for _ in range(reps):
                p = 0.5
                for _ in range(int(round(1.0 / step))):
                    p = _step(p, step, rng, tol)
                values.append(p)
                hetero.append(2 * p * (1 - p))
                absorbed.append(p == 0.0 or p == 1.0)
            values = np.asarray(values)
            rows.append({"max_step_cu": step, "tolerance": tol, "p0": 0.5, "u": 1.0,
                         "mean_terminal_frequency": values.mean(), "mean_frequency_se": _se(values),
                         "variance_terminal_frequency": values.var(ddof=1),
                         "heterozygosity": np.mean(hetero), "heterozygosity_se": _se(np.asarray(hetero)),
                         "finite_time_absorption": np.mean(absorbed), "replicates": reps})
    _write(out / "results/convergence.tsv", rows)


def _scenario(process: str, ne: int, rng: np.random.Generator, step: float, tol: float,
              trees_per_history: int) -> tuple[float, tuple[float, float, float], float, float]:
    tree = _tree(ne, root_extension=20.0)
    rr = Rearrangement("inv", "inversion", "ROOT", 2, max(1, ne // 2), 0.0)
    history = simulate_population_history(tree, rr, rng, _config(process, step, tol))
    root_p = history.frequency_at("ROOT", tree.branches["ROOT"].younger_age)
    tops = []
    coal = []
    patterns = []
    for _ in range(trees_per_history):
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        result = simulate_genealogy(0, tree, history, sampled, 0.01, 0.1, rng, record_events=False)
        tops.append(result.topology_index)
        coal.append(float(np.mean(result.coalescence_times)))
        patterns.append("".join(str(sampled[t]) for t in tree.taxa))
    top_fractions = tuple(float(np.mean(np.asarray(tops) == k)) for k in range(3))
    return root_p, top_fractions, float(np.mean(coal)), float(np.mean([p.count("1") for p in patterns]))


def process_comparison(out: Path, reps: int, seed: int) -> None:
    rows = []
    by_process = {"wright_fisher": [], "moran": []}
    for ne in (20, 100, 500):
        for process in ("wright_fisher", "moran", "wf_diffusion"):
            process_offset = {"wright_fisher": 100, "moran": 200, "wf_diffusion": 300}[process]
            rng = np.random.default_rng(seed + ne + process_offset)
            values = []
            coal = []
            pattern = []
            top = {0: [], 1: [], 2: []}
            for _ in range(reps):
                p, top_fractions, c, pat = _scenario(process, ne, rng, 0.01, 0.001, 5)
                values.append(p); coal.append(c); pattern.append(pat)
                for k in range(3):
                    top[k].append(top_fractions[k])
            values = np.asarray(values)
            rows_process = []
            for statistic, arr, expected in (
                ("mean_terminal_frequency", values, values.mean()),
                ("variance_terminal_frequency", values, values.var(ddof=1)),
                ("mean_coalescence_time", np.asarray(coal), np.mean(coal)),
                ("mean_terminal_A1_count", np.asarray(pattern), np.mean(pattern)),
            ):
                rows_process.append({"process": process, "Ne": ne, "statistic": statistic,
                                     "estimate": expected, "mc_se": _se(arr), "replicates": reps,
                                     "trees_per_history": 5})
            for k in range(3):
                arr = np.asarray(top[k], dtype=float)
                rows_process.append({"process": process, "Ne": ne, "statistic": f"quartet_topology_{k}",
                                     "estimate": arr.mean(), "mc_se": _se(arr), "replicates": reps,
                                     "trees_per_history": 5})
            by_process["wright_fisher" if process == "wright_fisher" else "moran"].extend(rows_process)
            rows.extend(rows_process)
    _write(out / "results/wf_comparison.tsv", by_process["wright_fisher"])
    _write(out / "results/moran_comparison.tsv", by_process["moran"])
    _write(out / "results/genealogy_comparison.tsv", rows)


def runtime(out: Path, seed: int) -> None:
    config = json.loads(json.dumps({
        "species_tree": {"newick": "", "default_effective_population_size": 10000},
    }))
    example = Path("examples/arbitrary_30taxa_diffusion.yaml")
    import yaml
    with example.open() as handle:
        config = yaml.safe_load(handle)
    tree = tree_from_config(config)
    rows = []
    for step in (0.05, 0.01, 0.002):
        rng = np.random.default_rng(seed + round(step * 1e6))
        rr = Rearrangement("inv_1", "inversion", "ROOT", 2, 1, 0.0)
        start = time.perf_counter()
        history = simulate_population_history(tree, rr, rng, {"population_process": {"model": "wf_diffusion", "diffusion": {"max_step_cu": step, "tolerance": 0.001}}})
        history_seconds = time.perf_counter() - start
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        for loci in (20, 500):
            start = time.perf_counter()
            results = [simulate_genealogy(i, tree, history, sampled, 0.01, 0.1, rng, record_events=False) for i in range(loci)]
            genealogy_seconds = time.perf_counter() - start
            rows.append({"Ne": 10000, "taxa": len(tree.taxa), "max_step_cu": step, "loci": loci,
                         "history_seconds": history_seconds, "genealogy_seconds": genealogy_seconds,
                         "total_seconds": history_seconds + genealogy_seconds,
                         "history_records": len(history.records), "gene_trees": len(results),
                         "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         "seed": seed})
    _write(out / "results/runtime_benchmark.tsv", rows)


def figures(out: Path) -> None:
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    analytical_rows = list(csv.DictReader((out.parent / "results/analytical_moments.tsv").open(), delimiter="\t"))
    p = 0.5
    rows = [r for r in analytical_rows if float(r["p0"]) == p and r["statistic"] == "heterozygosity"]
    rows.sort(key=lambda r: float(r["u"]))
    x = np.array([float(r["u"]) for r in rows]); y = np.array([float(r["observed"]) for r in rows])
    plt.figure(); plt.plot(x, y, "o", label="simulation"); plt.plot(x, 2*p*(1-p)*np.exp(-x), label="theory"); plt.xlabel("coalescent time u"); plt.ylabel("heterozygosity"); plt.legend(); plt.tight_layout(); plt.savefig(out / "heterozygosity_decay.png", dpi=160); plt.close()
    fix = [r for r in analytical_rows if r["statistic"] == "fixation_probability"]
    fix.sort(key=lambda r: float(r["p0"]))
    fx = np.array([float(r["p0"]) for r in fix]); fy = np.array([float(r["observed"]) for r in fix]); fe = np.array([float(r["mc_se"]) for r in fix])
    plt.figure(); plt.errorbar(fx, fy, yerr=1.96 * fe, fmt="o", label="simulation"); plt.plot(fx, fx, label="theory"); plt.xlabel("initial frequency"); plt.ylabel("eventual fixation probability"); plt.legend(); plt.tight_layout(); plt.savefig(out / "fixation_probability.png", dpi=160); plt.close()
    conv = list(csv.DictReader((out.parent / "results/convergence.tsv").open(), delimiter="\t"))
    for tol in sorted({float(r["tolerance"]) for r in conv}):
        q = [r for r in conv if float(r["tolerance"]) == tol]; q.sort(key=lambda r: float(r["max_step_cu"]))
        plt.plot([float(r["max_step_cu"]) for r in q], [float(r["mean_terminal_frequency"]) for r in q], "o-", label=f"tol={tol:g}")
    plt.axhline(.5, color="k", linestyle="--", label="martingale target"); plt.xlabel("max_step_cu"); plt.ylabel("mean terminal frequency"); plt.legend(); plt.tight_layout(); plt.savefig(out / "convergence.png", dpi=160); plt.close()
    run = list(csv.DictReader((out.parent / "results/runtime_benchmark.tsv").open(), delimiter="\t"))
    for loci in (20, 500):
        q = [r for r in run if int(r["loci"]) == loci]; q.sort(key=lambda r: float(r["max_step_cu"]))
        plt.plot([float(r["max_step_cu"]) for r in q], [float(r["total_seconds"]) for r in q], "o-", label=f"{loci} loci")
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("max_step_cu"); plt.ylabel("seconds"); plt.legend(); plt.tight_layout(); plt.savefig(out / "runtime_precision.png", dpi=160); plt.close()
    genealogy = list(csv.DictReader((out.parent / "results/genealogy_comparison.tsv").open(), delimiter="\t"))
    for process, marker in (("wright_fisher", "o"), ("wf_diffusion", "s"), ("moran", "^")):
        q = [r for r in genealogy if r["process"] == process and r["statistic"] == "quartet_topology_0"]
        q.sort(key=lambda r: int(r["Ne"]))
        if q:
            plt.plot([int(r["Ne"]) for r in q], [float(r["estimate"]) for r in q], marker + "-", label=process)
    plt.xscale("log"); plt.xlabel("Ne"); plt.ylabel("quartet topology 0 probability"); plt.legend(); plt.tight_layout(); plt.savefig(out / "wf_diffusion_quartet_comparison.png", dpi=160); plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("validation/wf_diffusion"))
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261011)
    parser.add_argument("--stage", choices=("all", "analytical", "convergence", "process", "runtime"), default="all")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = {"seed": args.seed, "replicates": args.replicates, "stage": args.stage, "commit": "working-tree"}
    (args.output / "README.md").write_text("Validation outputs are generated by validation/wf_diffusion/scripts/run_validation.py.\n")
    (args.output / "results/metadata.json").parent.mkdir(parents=True, exist_ok=True)
    (args.output / "results/metadata.json").write_text(json.dumps(metadata, indent=2))
    if args.stage in ("all", "analytical"): analytical(args.output, args.replicates, args.seed)
    if args.stage in ("all", "convergence"): convergence(args.output, args.replicates, args.seed)
    if args.stage in ("all", "process"): process_comparison(args.output, args.replicates, args.seed)
    if args.stage in ("all", "runtime"): runtime(args.output, args.seed)
    if args.stage == "all": figures(args.output / "figures")


if __name__ == "__main__":
    main()
