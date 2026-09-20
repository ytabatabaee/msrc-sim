from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import yaml

from run_daughter_sorting_eta import _write_csv, run_forward_sweep


def _mass_conserving_config(config: dict[str, Any]) -> dict[str, Any]:
    copied = yaml.safe_load(yaml.safe_dump(config))
    copied["eta_design"]["primary_split_policy"]["model"] = "mass_conserving"
    return copied


def _summary_rows(rows: list[Mapping[str, Any]], diagnostics: list[Mapping[str, Any]], *, kernel: str) -> list[dict[str, Any]]:
    by_eta_errors: dict[float, list[float]] = {}
    for row in diagnostics:
        eta = float(row["eta"])
        parent = float(row["parent_frequency_A1"])
        target_mean = (float(row["target_frequency_positive"]) + float(row["target_frequency_negative"])) / 2.0
        by_eta_errors.setdefault(eta, []).append(abs(target_mean - parent))

    out = []
    for row in rows:
        eta = float(row["eta"])
        errors = by_eta_errors.get(eta, [])
        out.append({
            "kernel": kernel,
            "split_mode": row["split_mode"],
            "eta": eta,
            "n_histories": int(row["n_histories"]),
            "P_persistence": row["persistence_p"],
            "P_strict_2to2": row["strict_2_2_p"],
            "P_discordant_2to2": row["discordant_2_2_p"],
            "P_discordant_given_persistence": row["p_discordant_given_persistent"],
            "P_1010": row["pattern_1010_p"],
            "P_0101": row["pattern_0101_p"],
            "persistence_ci_low": row["persistence_ci_low"],
            "persistence_ci_high": row["persistence_ci_high"],
            "discordant_ci_low": row["discordant_2_2_ci_low"],
            "discordant_ci_high": row["discordant_2_2_ci_high"],
            "conditional_ci_low": row["conditional_ci_low"],
            "conditional_ci_high": row["conditional_ci_high"],
            "d_without_p_count": row["d_without_p_count"],
            "mean_realized_daughter_difference": row["mean_realized_daughter_difference"],
            "max_target_mass_conservation_error": max(errors) if errors else float("nan"),
        })
    return out


def _plot(rows: list[Mapping[str, Any]], output: Path) -> None:
    main = [row for row in rows if row["split_mode"] == "descendant" and float(row["eta"]) >= 0.0]
    main = sorted(main, key=lambda row: float(row["eta"]))
    etas = np.array([float(row["eta"]) for row in main])

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), constrained_layout=True)
    ax = axes[0]
    for key, label, color in (
        ("P_persistence", "P(persistence)", "#1f77b4"),
        ("P_discordant_2to2", "P(discordant 2:2)", "#2ca02c"),
        ("P_discordant_given_persistence", "P(discordant | persistence)", "#d62728"),
    ):
        ax.plot(etas, [float(row[key]) for row in main], marker="o", linewidth=1.6, label=label, color=color)
    ax.set_xlabel("eta")
    ax.set_ylabel("Probability")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for key, label, color in (
        ("P_1010", "1010", "#9467bd"),
        ("P_0101", "0101", "#ff7f0e"),
    ):
        ax.plot(etas, [float(row[key]) for row in main], marker="o", linewidth=1.6, label=label, color=color)
    ax.set_xlabel("eta")
    ax.set_ylabel("Labeled pattern frequency")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Mass-conserving speciation daughter sorting")
    fig.savefig(output / "mass_conserving_eta_main.pdf", bbox_inches="tight")
    fig.savefig(output / "mass_conserving_eta_main.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_markdown(rows: list[Mapping[str, Any]], output: Path) -> None:
    by_eta = {float(row["eta"]): row for row in rows}
    main_etas = [0.0, 0.5, 1.0, 2.0, 4.0]
    conditional = "\n".join(
        f"- eta = {eta:g}: P(discordant | persistence) = "
        f"{float(by_eta[eta]['P_discordant_given_persistence']):.4f}"
        for eta in main_etas
    )
    plus = by_eta[2.0]
    minus = by_eta[-2.0]
    max_error = max(float(row["max_target_mass_conservation_error"]) for row in rows if math.isfinite(float(row["max_target_mass_conservation_error"])))
    robust = float(by_eta[4.0]["P_discordant_given_persistence"]) > float(by_eta[0.0]["P_discordant_given_persistence"])

    output.joinpath("summary.md").write_text(
        "# Mass-conserving daughter-sorting robustness check\n\n"
        "Kernel: for parental frequency `p`, set `s = tanh(eta)` and "
        "`delta = s * min(p, 1 - p)`, then assign the positive daughter "
        "`p + delta` and the negative daughter `p - delta`.\n\n"
        f"The maximum absolute error in `(p_positive + p_negative) / 2 - p_parent` "
        f"across recorded polymorphic split events was `{max_error:.3g}`.\n\n"
        "Conditional discordant strict 2:2 probabilities:\n\n"
        f"{conditional}\n\n"
        "Sign-symmetry check at eta = ±2:\n\n"
        f"- eta = +2: P(1010) = {float(plus['P_1010']):.4f}, "
        f"P(0101) = {float(plus['P_0101']):.4f}, "
        f"P(discordant | persistence) = {float(plus['P_discordant_given_persistence']):.4f}\n"
        f"- eta = -2: P(1010) = {float(minus['P_1010']):.4f}, "
        f"P(0101) = {float(minus['P_0101']):.4f}, "
        f"P(discordant | persistence) = {float(minus['P_discordant_given_persistence']):.4f}\n\n"
        "Conclusion: the qualitative descendant-split eta result is "
        f"{'robust' if robust else 'not supported in this run'} under the "
        "mass-conserving kernel: stronger daughter sorting increases "
        "discordant strict 2:2 sorting conditional on ancestral persistence.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the mass-conserving daughter-sorting robustness experiment.")
    parser.add_argument("--config", default="experiments/daughter_sorting_eta/config.yaml")
    parser.add_argument("--output", default="results/daughter_sorting_eta/mass_conserving")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    config = _mass_conserving_config(yaml.safe_load(Path(args.config).read_text()))
    master = np.random.SeedSequence(int(config["experiment"]["master_seed"]))
    streams = master.spawn(2)

    main_summary, _main_reps, main_diag, _ = run_forward_sweep(
        config,
        output,
        etas=[0.0, 0.5, 1.0, 2.0, 4.0],
        histories_per_eta=int(config["experiment"]["main_histories_per_eta"]),
        process="wright_fisher",
        split_mode="descendant",
        master_ss=streams[0],
        retain_genealogy_histories=False,
    )
    sign_summary, _sign_reps, sign_diag, _ = run_forward_sweep(
        config,
        output,
        etas=[-2.0],
        histories_per_eta=int(config["experiment"]["sign_symmetry_histories_per_eta"]),
        process="wright_fisher",
        split_mode="descendant",
        master_ss=streams[1],
        retain_genealogy_histories=False,
    )
    for row in sign_summary:
        row["split_mode"] = "sign_symmetry"
    for row in sign_diag:
        row["split_mode"] = "sign_symmetry"

    summary = _summary_rows(main_summary + sign_summary, main_diag + sign_diag, kernel="mass_conserving")
    if any(int(row["d_without_p_count"]) != 0 for row in summary):
        raise RuntimeError("Discordant strict 2:2 without persistence observed")

    fields = list(summary[0])
    _write_csv(output / "mass_conserving_eta_summary.csv", summary, fields)
    _plot(summary, output)
    _write_markdown(summary, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
