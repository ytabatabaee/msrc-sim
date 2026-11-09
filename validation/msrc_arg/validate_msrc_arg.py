"""Replicate-level validation for the simplified two-background MSRC ARG."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy


class FrozenHistory:
    # The synthetic validation history becomes ordinary MSC above the root
    # after a finite structured interval, avoiding an unbounded structured
    # phase that is not part of the test target.
    structured_origin_age = 100.0

    def __init__(self, p=0.5):
        self.p = float(p)

    def frequency_at(self, branch_id, age):
        return self.p

    def next_frequency_boundary(self, branch_id, age):
        return float("inf")

    def terminal_frequency(self, taxon):
        return self.p


def quartet_topology(tree):
    full = frozenset(tree.samples()); splits = []
    for node in tree.nodes():
        if node == tree.root:
            continue
        side = frozenset(tree.samples(node))
        if len(side) == 2:
            splits.append(tuple(sorted((tuple(sorted(side)), tuple(sorted(full - side))))))
    return min(splits)


def transition_targets():
    rows = []
    history = FrozenHistory(.37)
    for kappa in (0.0, .01, .1, .5, 1.0):
        model = MSRCTwoBackground(history, 500.0, (250.0, 750.0), kappa)
        rec_map = PiecewiseConstantRecombinationMap.constant(1000, 1e-5)
        for state in (0, 1):
            for focal in (100.0, 500.0, 900.0):
                lo, hi = sorted((focal, 500.0)); distance = hi - lo
                if distance == 0:
                    expected = 0.0
                else:
                    expected = history.p if state == 0 else 1 - history.p
                    outside = sum(max(0.0, min(hi, b) - max(lo, a))
                                  for a, b in ((0.0, 250.0), (750.0, 1000.0)))
                    inside = distance - outside
                    expected *= 1e-5 * (outside + kappa * inside)
                rows.append({"kappa": kappa, "state": state, "focal_bp": focal,
                             "expected_cross_background_rate": expected,
                             "target_formula": "p_g * integral(r0 * eta)"})
    return rows


def run_replicates(replicates, seed):
    tree = SpeciesTree("((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;", 20, 100)
    history = FrozenHistory(.5)
    rows = []
    for kappa in (0.0, .01, .1, .5, 1.0):
        model = MSRCTwoBackground(history, 500.0, (250.0, 750.0), kappa)
        for i in range(replicates):
            length = 1000.0
            rec_map = PiecewiseConstantRecombinationMap.constant(length, 1e-5)
            cfg = ARGConfig(length, rec_map, model, True, history,
                            {"1": 1, "2": 1, "3": 0, "4": 0})
            result = simulate_arg(tree, cfg, np.random.default_rng(seed + i + int(kappa * 10000)))
            ts = build_tree_sequence(result).simplify()
            events = [e for e in result.events if e.event_type in {"recombination", "background_transfer"}]
            times = [float(ts.at(x).time(ts.at(x).mrca(0, 1))) for x in (100.0, 900.0)]
            rows.append({"kappa": kappa, "recombination_events": sum(e.event_type == "recombination" for e in events),
                         "background_transfer_events": sum(e.event_type == "background_transfer" for e in events),
                         "inside_recombination": sum(e.event_type == "recombination" and 250 <= e.position < 750 for e in events),
                         "outside_recombination": sum(e.event_type == "recombination" and not (250 <= e.position < 750) for e in events),
                         "marginal_trees": ts.num_trees, "topology": str(quartet_topology(ts.first())),
                         "topology_focal": str(quartet_topology(ts.at(100.0))),
                         "mrca_left": times[0], "mrca_right": times[1]})
    return rows


def single_locus_reference(replicates, seed):
    tree = SpeciesTree("((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;", 20, 100)
    history = FrozenHistory(.5); sampled = {"1": 1, "2": 1, "3": 0, "4": 0}
    output = {}
    for kappa in (0.0, .01, .1, .5, 1.0):
        for focal in (100.0, 900.0):
            lo, hi = sorted((focal, 500.0)); outside = sum(max(0.0, min(hi, b) - max(lo, a)) for a, b in ((0.0, 250.0), (750.0, 1000.0)))
            inside = hi - lo - outside
            rho_eff = 1e-5 * (outside + kappa * inside)
            results = [simulate_genealogy(i, tree, history, sampled, rho_eff, 1.0,
                                          np.random.default_rng(seed + i), record_events=False)
                       for i in range(replicates)]
            counts = np.bincount([r.topology_index for r in results], minlength=3) / replicates
            output[f"{kappa}:{focal}"] = {"rho_eff": rho_eff, "q": counts.tolist(),
                                           "mean_coalescence_time": float(np.mean([r.coalescence_times[-1] for r in results]))}
    return output


def aggregate(rows):
    output = {}
    for kappa in sorted(set(r["kappa"] for r in rows)):
        sub = [r for r in rows if r["kappa"] == kappa]
        left = np.asarray([r["mrca_left"] for r in sub]); right = np.asarray([r["mrca_right"] for r in sub])
        output[str(kappa)] = {
            "replicates": len(sub),
            "mean_recombination_events": float(np.mean([r["recombination_events"] for r in sub])),
            "mean_background_transfer_events": float(np.mean([r["background_transfer_events"] for r in sub])),
            "mean_inside_recombination": float(np.mean([r["inside_recombination"] for r in sub])),
            "mean_outside_recombination": float(np.mean([r["outside_recombination"] for r in sub])),
            "mean_marginal_trees": float(np.mean([r["marginal_trees"] for r in sub])),
            "mrca_correlation": float(np.corrcoef(left, right)[0, 1]),
            "mean_mrca_left": float(np.mean(left)), "mean_mrca_right": float(np.mean(right)),
            "focal_topology_frequencies": {topology: float(np.mean([r["topology_focal"] == topology for r in sub]))
                                           for topology in sorted(set(r["topology_focal"] for r in sub))},
        }
    return output


def write_figure(summary, path):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    x = [float(k) for k in summary]
    plt.figure(figsize=(7, 4))
    plt.plot(x, [summary[k]["mean_background_transfer_events"] for k in summary], marker="o", label="background transfers")
    plt.plot(x, [summary[k]["mean_recombination_events"] for k in summary], marker="o", label="physical recombinations")
    plt.xlabel("Cross-background suppression κ"); plt.ylabel("Events per chromosome"); plt.legend(); plt.tight_layout()
    plt.savefig(path, dpi=160); plt.close(); return True


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--output", type=Path, default=Path("validation/msrc_arg/results/msrc_arg_validation_1000.json"))
    parser.add_argument("--figure", type=Path, default=Path("validation/msrc_arg/results/msrc_suppression_events.png"))
    args = parser.parse_args(argv)
    rows = run_replicates(args.replicates, args.seed)
    summary = {"replicates_per_kappa": args.replicates, "seed": args.seed,
               "transition_rate_targets": transition_targets(), "suppression_summary": aggregate(rows),
               "single_locus_reference": single_locus_reference(args.replicates, args.seed),
               "notes": ["Frozen p(t)=0.5 history; fixed terminal backgrounds 1,1,0,0.",
                         "Neutral marginal tree statistics use tskit.simplify(); event logs retain ancestry-invisible transfers."]}
    args.output.parent.mkdir(parents=True, exist_ok=True); args.figure.parent.mkdir(parents=True, exist_ok=True)
    summary["figure_generated"] = write_figure(summary["suppression_summary"], args.figure)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
