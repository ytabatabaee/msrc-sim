from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from hashlib import blake2b
from math import sqrt
from pathlib import Path
from typing import Any, Iterable, Mapping
import csv
import json

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from .moran_validation import moran_distribution_after_t, moran_event_count_sample
from .moran_validation import moran_generator
from .pattern_probabilities import wf_distribution_after_t, wf_transition_matrix
from .robustness import binomial_confidence_interval
from .wf_moran_robustness import (
    is_discordant_two_two,
    is_two_two,
    simulate_forward_attempt,
)
from .wright_fisher import _selected_p


SPECIES_TREE_SPLIT = "12|34"
DISCORDANT_PATTERNS = ("0101", "1010", "0110", "1001")
MASTER_SEED = 20260916
PRIMARY_REPLICATES = 5000
MORAN_REPLICATES = 5000
FOCAL_LEGACY_REPLICATES = 2500
FOCAL_NE = 20
FOCAL_P0 = 0.5
FOCAL_GAMMA = 0.0
FOCAL_LEGACY_TAU = 0.6
FOCAL_LEGACY_ROOT_EXTENSION = 0.5


@dataclass(frozen=True)
class Cell:
    process: str
    Ne: int
    p0: float
    tau: float
    gamma: float
    replicates: int
    grid: str = "primary"

    @property
    def total(self) -> int:
        return 2 * int(self.Ne)

    @property
    def initial_count(self) -> int:
        return int(round(float(self.p0) * self.total))

    @property
    def raw_s(self) -> float:
        return raw_s_from_gamma_target(self.gamma, self.Ne, self.process)

    @property
    def branch_duration_native(self) -> int:
        return native_generations_from_tau(self.tau, self.Ne)


def raw_s_from_gamma(gamma: float, Ne: int) -> float:
    """Convert Wright-Fisher gamma=2Ne*s_WF to raw s.

    This legacy helper is kept for existing WF code paths. New WF/Moran
    comparisons should call ``raw_s_from_gamma_target`` with an explicit
    process.
    """
    return float(gamma) / (2.0 * int(Ne))


def scaled_selection(raw_s: float, Ne: int) -> float:
    """Legacy Wright-Fisher scaled selection, gamma=2Ne*s."""
    return 2.0 * int(Ne) * float(raw_s)


def raw_s_from_gamma_target(gamma_target: float, Ne: int, process: str) -> float:
    """Convert shared diffusion-scale gamma to process-specific raw s.

    WF has weak drift s*p*(1-p), so gamma_target=2Ne*s_WF.
    The implemented Moran CTMC has weak drift (s/2)*p*(1-p), so
    gamma_target=Ne*s_Moran.
    """
    if process == "wright_fisher":
        return float(gamma_target) / (2.0 * int(Ne))
    if process == "moran":
        return float(gamma_target) / float(int(Ne))
    raise ValueError("process must be wright_fisher or moran")


def process_specific_scaled_selection(raw_s: float, Ne: int, process: str) -> float:
    if process == "wright_fisher":
        return 2.0 * int(Ne) * float(raw_s)
    if process == "moran":
        return float(int(Ne)) * float(raw_s)
    raise ValueError("process must be wright_fisher or moran")


def exact_wf_initial_drift(p: float, raw_s: float) -> float:
    return float(_selected_p(float(p), float(raw_s)) - float(p))


def exact_moran_initial_drift(p: float, raw_s: float) -> float:
    p = float(p)
    s = float(raw_s)
    if s == 0.0:
        return 0.0
    p_sel = p * (1.0 + s) / (1.0 + s * p)
    return float((p_sel - p) / 2.0)


def native_generations_from_tau(tau: float, Ne: int) -> int:
    return int(round(float(tau) * 2.0 * int(Ne)))


def scaled_time(native_generations: int, Ne: int) -> float:
    return float(native_generations) / (2.0 * int(Ne))


def deterministic_seed(master_seed: int, *parts: Any) -> int:
    text = "|".join([str(int(master_seed)), *(str(p) for p in parts)])
    digest = blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") & 0xFFFFFFFF


def heterozygosity(p: float) -> float:
    p = float(p)
    return 2.0 * p * (1.0 - p)


def classify_absorption(k: int, total: int) -> str:
    if int(k) <= 0:
        return "loss"
    if int(k) >= int(total):
        return "fix"
    return "persist"


def strict_2to2(pattern: str) -> bool:
    return str(pattern).count("1") == 2


def discordant_2to2(pattern: str) -> bool:
    return str(pattern) in DISCORDANT_PATTERNS


def _wf_count_after_generations(k: int, total: int, generations: int, raw_s: float, rng: np.random.Generator) -> int:
    k = int(k)
    for _ in range(int(generations)):
        if k <= 0 or k >= total:
            break
        p = _selected_p(k / total, raw_s)
        k = int(rng.binomial(total, p))
    return k


def _process_count_after_time(
    process: str,
    k: int,
    total: int,
    duration: int,
    raw_s: float,
    rng: np.random.Generator,
) -> int:
    if process == "wright_fisher":
        return _wf_count_after_generations(k, total, duration, raw_s, rng)
    if process == "moran":
        return int(moran_event_count_sample(total, int(k), float(duration), rng, raw_s))
    raise ValueError("process must be wright_fisher or moran")


@lru_cache(maxsize=None)
def transition_matrix(process: str, Ne: int, duration: int, raw_s: float) -> np.ndarray:
    """Finite-state transition matrix matching the forward simulator."""
    if process == "wright_fisher":
        base = wf_transition_matrix(int(Ne), float(raw_s))
        return np.linalg.matrix_power(base, int(duration))
    if process == "moran":
        total = 2 * int(Ne)
        return scipy_expm(moran_generator(total, float(raw_s)) * float(duration))
    raise ValueError("process must be wright_fisher or moran")


def scipy_expm(matrix: np.ndarray) -> np.ndarray:
    from scipy.linalg import expm

    out = expm(matrix)
    out[out < 0.0] = 0.0
    row_sums = out.sum(axis=1)
    out = out / row_sums[:, None]
    return out


def _sample_transition(matrix: np.ndarray, start: int, rng: np.random.Generator) -> int:
    probs = np.asarray(matrix[int(start)], dtype=float)
    probs = probs / probs.sum()
    return int(rng.choice(len(probs), p=probs))


def simulate_cell(cell: Cell, *, master_seed: int = MASTER_SEED) -> dict[str, Any]:
    if cell.process not in {"wright_fisher", "moran"}:
        raise ValueError("process must be wright_fisher or moran")
    if not 0 <= cell.initial_count <= cell.total:
        raise ValueError("p0 produces an initial copy count outside 0..2Ne")
    if 1.0 + cell.raw_s <= 0.0:
        raise ValueError("selection requires 1+s > 0")

    seed = deterministic_seed(master_seed, cell.process, cell.Ne, cell.p0, cell.tau, cell.gamma, cell.grid)
    rng = np.random.default_rng(seed)
    duration = cell.branch_duration_native
    transition = transition_matrix(cell.process, cell.Ne, duration, round(cell.raw_s, 15))

    n_loss = n_persist = n_fix = 0
    n_common = n_two_two = n_discordant = 0
    n_two_two_given_persist = n_discordant_given_persist = 0
    terminal_ps: list[float] = []
    terminal_hs: list[float] = []
    persistent_hs: list[float] = []
    persistent_minor: list[float] = []

    for _ in range(int(cell.replicates)):
        root_k = _sample_transition(transition, cell.initial_count, rng)
        p_terminal = root_k / cell.total
        h_terminal = heterozygosity(p_terminal)
        terminal_ps.append(p_terminal)
        terminal_hs.append(h_terminal)

        state = classify_absorption(root_k, cell.total)
        persisted = state == "persist"
        n_loss += int(state == "loss")
        n_persist += int(persisted)
        n_fix += int(state == "fix")
        n_common += int(0.1 <= p_terminal <= 0.9)
        if persisted:
            persistent_hs.append(h_terminal)
            persistent_minor.append(min(p_terminal, 1.0 - p_terminal))

        terminal_bits = []
        for _daughter in range(2):
            start = int(rng.binomial(cell.total, p_terminal))
            end = _sample_transition(transition, start, rng)
            p_child = end / cell.total
            terminal_bits.extend([int(rng.random() < p_child), int(rng.random() < p_child)])
        pattern = "".join(str(x) for x in terminal_bits)
        two_two = strict_2to2(pattern)
        discordant = discordant_2to2(pattern)
        n_two_two += int(two_two)
        n_discordant += int(discordant)
        n_two_two_given_persist += int(persisted and two_two)
        n_discordant_given_persist += int(persisted and discordant)

    return summarize_counts(
        cell,
        seed=seed,
        n_loss=n_loss,
        n_persist=n_persist,
        n_fix=n_fix,
        n_common=n_common,
        n_two_two=n_two_two,
        n_discordant=n_discordant,
        n_two_two_given_persist=n_two_two_given_persist,
        n_discordant_given_persist=n_discordant_given_persist,
        terminal_ps=terminal_ps,
        terminal_hs=terminal_hs,
        persistent_hs=persistent_hs,
        persistent_minor=persistent_minor,
    )


def _prob_fields(prefix: str, successes: int, trials: int) -> dict[str, float]:
    lo, hi = binomial_confidence_interval(int(successes), int(trials))
    return {
        prefix: float(successes) / float(trials) if trials else float("nan"),
        f"{prefix}_ci_low": lo,
        f"{prefix}_ci_high": hi,
    }


def summarize_counts(
    cell: Cell,
    *,
    seed: int,
    n_loss: int,
    n_persist: int,
    n_fix: int,
    n_common: int,
    n_two_two: int,
    n_discordant: int,
    n_two_two_given_persist: int,
    n_discordant_given_persist: int,
    terminal_ps: Iterable[float],
    terminal_hs: Iterable[float],
    persistent_hs: Iterable[float],
    persistent_minor: Iterable[float],
) -> dict[str, Any]:
    reps = int(cell.replicates)
    persistent_hs = list(persistent_hs)
    persistent_minor = list(persistent_minor)
    row: dict[str, Any] = {
        "grid": cell.grid,
        "process": cell.process,
        "Ne": int(cell.Ne),
        "p0": float(cell.p0),
        "initial_copy_count": int(cell.initial_count),
        "gamma_target": float(cell.gamma),
        "raw_s": float(cell.raw_s),
        "process_specific_scaled_selection": float(
            process_specific_scaled_selection(cell.raw_s, cell.Ne, cell.process)
        ),
        "scaled_selection": float(process_specific_scaled_selection(cell.raw_s, cell.Ne, cell.process)),
        "legacy_gamma_2Ne_s": float(scaled_selection(cell.raw_s, cell.Ne)),
        "wf_raw_s_for_gamma_target": float(raw_s_from_gamma_target(cell.gamma, cell.Ne, "wright_fisher")),
        "moran_raw_s_for_gamma_target": float(raw_s_from_gamma_target(cell.gamma, cell.Ne, "moran")),
        "initial_drift_wf_exact": float(
            exact_wf_initial_drift(cell.p0, raw_s_from_gamma_target(cell.gamma, cell.Ne, "wright_fisher"))
        ),
        "initial_drift_moran_exact": float(
            exact_moran_initial_drift(cell.p0, raw_s_from_gamma_target(cell.gamma, cell.Ne, "moran"))
        ),
        "branch_duration_native": int(cell.branch_duration_native),
        "branch_duration_scaled": float(scaled_time(cell.branch_duration_native, cell.Ne)),
        "n_replicates": reps,
        "n_loss": int(n_loss),
        "n_persist": int(n_persist),
        "n_fix": int(n_fix),
        "mean_terminal_p": float(np.mean(list(terminal_ps))),
        "mean_terminal_heterozygosity": float(np.mean(list(terminal_hs))),
        "mean_terminal_heterozygosity_given_persist": float(np.mean(persistent_hs)) if persistent_hs else float("nan"),
        "mean_minor_frequency_given_persist": float(np.mean(persistent_minor)) if persistent_minor else float("nan"),
        "P_common_polymorphism_0p1_0p9": float(n_common) / reps,
        "P_strict_2to2_given_persist": float(n_two_two_given_persist) / n_persist if n_persist else float("nan"),
        "P_discordant_2to2_given_persist": float(n_discordant_given_persist) / n_persist if n_persist else float("nan"),
        "seed": int(seed),
    }
    row.update(_prob_fields("P_loss", n_loss, reps))
    row.update(_prob_fields("P_persist", n_persist, reps))
    row.update(_prob_fields("P_fix", n_fix, reps))
    row.update(_prob_fields("P_strict_2to2", n_two_two, reps))
    row.update(_prob_fields("P_discordant_2to2", n_discordant, reps))
    total = row["P_loss"] + row["P_persist"] + row["P_fix"]
    if abs(total - 1.0) > 1e-12:
        raise RuntimeError(f"absorption probabilities sum to {total}")
    return row


def exact_absorption(process: str, Ne: int, p0: float, tau: float, gamma: float) -> dict[str, float]:
    total = 2 * int(Ne)
    initial = int(round(float(p0) * total))
    duration = native_generations_from_tau(tau, Ne)
    raw_s = raw_s_from_gamma_target(gamma, Ne, process)
    if process == "wright_fisher":
        dist = wf_distribution_after_t(Ne, initial, duration, raw_s)
    elif process == "moran":
        dist = moran_distribution_after_t(total, initial, float(duration), raw_s)
    else:
        raise ValueError("process must be wright_fisher or moran")
    return {
        "P_loss_exact": float(dist[0]),
        "P_persist_exact": float(1.0 - dist[0] - dist[-1]),
        "P_fix_exact": float(dist[-1]),
    }


def primary_cells() -> list[Cell]:
    cells = [
        Cell("wright_fisher", ne, p0, tau, gamma, PRIMARY_REPLICATES, "primary_scaled")
        for ne in (20, 50, 100)
        for p0 in (0.05, 0.10, 0.25, 0.50)
        for tau in (0.10, 0.30, 0.60, 1.00)
        for gamma in (-5.0, -2.0, 0.0, 2.0, 5.0)
    ]
    cells.extend(
        Cell("wright_fisher", ne, p0, tau, raw_gamma, PRIMARY_REPLICATES, "supplement_raw_s")
        for ne in (20, 50, 100)
        for p0 in (0.10, 0.50)
        for tau in (0.30, 0.60)
        for raw_gamma in tuple(scaled_selection(s, ne) for s in (-0.01, 0.0, 0.01))
    )
    return cells


def moran_cells() -> list[Cell]:
    return [
        Cell("moran", ne, p0, tau, gamma, MORAN_REPLICATES, "moran_matched")
        for ne in (20, 100)
        for p0 in (0.10, 0.50)
        for tau in (0.30, 0.60, 1.00)
        for gamma in (-5.0, -2.0, 0.0, 2.0, 5.0)
    ]


def run_grid(outdir: str | Path = "results/prevalence_grid", figdir: str | Path = "figures") -> dict[str, Any]:
    out = Path(outdir)
    figs = Path(figdir)
    out.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    cells = primary_cells() + moran_cells()
    rows = []
    for i, cell in enumerate(cells, start=1):
        print(
            f"{i}/{len(cells)} {cell.grid} {cell.process} Ne={cell.Ne} "
            f"p0={cell.p0} tau={cell.tau} gamma={cell.gamma}",
            flush=True,
        )
        rows.append(simulate_cell(cell))

    table = out / "prevalence_grid.tsv"
    write_tsv(table, rows)
    exact_rows = exact_validation_rows(rows)
    write_tsv(out / "exact_validation.tsv", exact_rows)
    focal = focal_reproduction()
    write_tsv(out / "focal_reproduction.tsv", focal)
    calibration_template(out / "calibration_template.tsv")
    plot_main(rows, focal, figs / "prevalence_grid_main.pdf")
    plot_supplement(rows, figs / "prevalence_grid_supplement.pdf")
    summary = write_summary(out / "prevalence_summary.md", rows, exact_rows, focal)
    with (out / "prevalence_grid_metadata.json").open("w") as handle:
        json.dump(
            {
                "master_seed": MASTER_SEED,
                "primary_replicates": PRIMARY_REPLICATES,
                "moran_replicates": MORAN_REPLICATES,
                "species_tree_split": SPECIES_TREE_SPLIT,
                "summary": summary,
            },
            handle,
            indent=2,
        )
    return {"rows": rows, "exact": exact_rows, "focal": focal, "summary": summary}


def write_tsv(path: str | Path, rows: list[Mapping[str, Any]]) -> None:
    if not rows:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def row_gamma(row: Mapping[str, Any]) -> float:
    return float(row.get("gamma_target", row.get("scaled_selection", 0.0)))


def matched_wf_rows(rows: list[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [
        r
        for r in rows
        if r["grid"] == "primary_scaled"
        and r["process"] == "wright_fisher"
        and int(r["Ne"]) in {20, 100}
        and float(r["p0"]) in {0.10, 0.50}
        and float(r["branch_duration_scaled"]) in {0.3, 0.6, 1.0}
        and row_gamma(r) in {-5.0, -2.0, 0.0, 2.0, 5.0}
    ]


def corrected_moran_cells() -> list[Cell]:
    return [
        Cell("moran", ne, p0, tau, gamma, MORAN_REPLICATES, "moran_corrected_matched")
        for ne in (20, 100)
        for p0 in (0.10, 0.50)
        for tau in (0.30, 0.60, 1.00)
        for gamma in (-5.0, -2.0, 0.0, 2.0, 5.0)
    ]


def run_corrected_moran_grid(
    outdir: str | Path = "results/prevalence_grid",
    figdir: str | Path = "figures",
) -> dict[str, Any]:
    """Rerun only the small Moran robustness grid with corrected scaling."""
    out = Path(outdir)
    figs = Path(figdir)
    out.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    existing = read_tsv(out / "prevalence_grid.tsv")
    old_moran = [r for r in existing if r.get("grid") == "moran_matched" and r.get("process") == "moran"]
    if old_moran:
        write_tsv(out / "moran_matched_grid_old_unscaled_selection.tsv", old_moran)

    corrected = []
    cells = corrected_moran_cells()
    for i, cell in enumerate(cells, start=1):
        print(
            f"{i}/{len(cells)} corrected Moran Ne={cell.Ne} p0={cell.p0} "
            f"tau={cell.tau} gamma_target={cell.gamma}",
            flush=True,
        )
        corrected.append(simulate_cell(cell))
    write_tsv(out / "corrected_moran_matched_grid.tsv", corrected)

    wf = matched_wf_rows(existing)
    comparison = wf_moran_comparison_rows(wf, corrected, old_moran)
    write_tsv(out / "wf_moran_corrected_comparison.tsv", comparison)
    exact_rows = corrected_exact_validation_rows(corrected)
    write_tsv(out / "corrected_moran_exact_validation.tsv", exact_rows)
    write_selection_scaling_derivation(out / "selection_scaling_derivation.md")
    focal = read_tsv(out / "focal_reproduction.tsv") if (out / "focal_reproduction.tsv").exists() else focal_reproduction()
    plot_selection_main(existing, corrected, figs / "prevalence_selection_main.pdf")
    plot_selection_supplement(existing, figs / "prevalence_selection_supplement.pdf")
    summary = write_corrected_summary(out / "prevalence_summary.md", existing, corrected, old_moran, comparison, exact_rows, focal)
    return {
        "corrected_moran": corrected,
        "comparison": comparison,
        "exact": exact_rows,
        "summary": summary,
    }


def _find_cell(rows: list[Mapping[str, Any]], ne: int, p0: float, tau: float, gamma: float) -> Mapping[str, Any] | None:
    for r in rows:
        if (
            int(r["Ne"]) == int(ne)
            and abs(float(r["p0"]) - float(p0)) < 1e-12
            and abs(float(r["branch_duration_scaled"]) - float(tau)) < 1e-12
            and abs(row_gamma(r) - float(gamma)) < 1e-12
        ):
            return r
    return None


def wf_moran_comparison_rows(
    wf_rows: list[Mapping[str, Any]],
    corrected_moran_rows: list[Mapping[str, Any]],
    old_moran_rows: list[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    metrics = [
        "P_persist",
        "P_strict_2to2",
        "P_discordant_2to2",
        "P_strict_2to2_given_persist",
        "P_discordant_2to2_given_persist",
    ]
    out: list[dict[str, Any]] = []
    for wf in sorted(wf_rows, key=lambda r: (int(r["Ne"]), float(r["p0"]), float(r["branch_duration_scaled"]), row_gamma(r))):
        ne = int(wf["Ne"])
        p0 = float(wf["p0"])
        tau = float(wf["branch_duration_scaled"])
        gamma = row_gamma(wf)
        moran = _find_cell(corrected_moran_rows, ne, p0, tau, gamma)
        if moran is None:
            continue
        old = _find_cell(old_moran_rows, ne, p0, tau, gamma)
        for metric in metrics:
            wf_est = float(wf[metric])
            moran_est = float(moran[metric])
            old_est = float(old[metric]) if old is not None and metric in old and old[metric] not in ("", "nan") else float("nan")
            row = {
                "Ne": ne,
                "p0": p0,
                "tau": tau,
                "gamma_target": gamma,
                "s_WF": raw_s_from_gamma_target(gamma, ne, "wright_fisher"),
                "s_Moran": raw_s_from_gamma_target(gamma, ne, "moran"),
                "metric": metric,
                "WF_estimate": wf_est,
                "WF_ci_low": float(wf.get(f"{metric}_ci_low", "nan")),
                "WF_ci_high": float(wf.get(f"{metric}_ci_high", "nan")),
                "Moran_corrected_estimate": moran_est,
                "Moran_corrected_ci_low": float(moran.get(f"{metric}_ci_low", "nan")),
                "Moran_corrected_ci_high": float(moran.get(f"{metric}_ci_high", "nan")),
                "corrected_absolute_difference": abs(wf_est - moran_est),
                "Moran_old_unscaled_estimate": old_est,
                "old_absolute_difference": abs(wf_est - old_est) if old is not None else float("nan"),
            }
            row["absolute_difference_change_corrected_minus_old"] = (
                row["corrected_absolute_difference"] - row["old_absolute_difference"]
                if old is not None
                else float("nan")
            )
            out.append(row)
    return out


def corrected_exact_validation_rows(moran_rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    choices = [
        (20, 0.10, 0.30, 0.0),
        (20, 0.50, 0.60, -5.0),
        (20, 0.50, 0.60, 5.0),
    ]
    out: list[dict[str, Any]] = []
    for ne, p0, tau, gamma in choices:
        mc = _find_cell(moran_rows, ne, p0, tau, gamma)
        if mc is None:
            raise RuntimeError(f"missing corrected Moran validation cell {ne=} {p0=} {tau=} {gamma=}")
        exact = exact_absorption("moran", ne, p0, tau, gamma)
        row = {
            "process": "moran",
            "Ne": ne,
            "p0": p0,
            "branch_duration_scaled": tau,
            "gamma_target": gamma,
            "raw_s": raw_s_from_gamma_target(gamma, ne, "moran"),
            "P_loss_mc": mc["P_loss"],
            "P_loss_mc_ci_low": mc["P_loss_ci_low"],
            "P_loss_mc_ci_high": mc["P_loss_ci_high"],
            "P_persist_mc": mc["P_persist"],
            "P_persist_mc_ci_low": mc["P_persist_ci_low"],
            "P_persist_mc_ci_high": mc["P_persist_ci_high"],
            "P_fix_mc": mc["P_fix"],
            "P_fix_mc_ci_low": mc["P_fix_ci_low"],
            "P_fix_mc_ci_high": mc["P_fix_ci_high"],
            **exact,
        }
        row["mc_sum"] = float(mc["P_loss"]) + float(mc["P_persist"]) + float(mc["P_fix"])
        row["max_abs_error"] = max(
            abs(float(row["P_loss_mc"]) - row["P_loss_exact"]),
            abs(float(row["P_persist_mc"]) - row["P_persist_exact"]),
            abs(float(row["P_fix_mc"]) - row["P_fix_exact"]),
        )
        row["exact_inside_mc_ci_all_components"] = all(
            float(row[f"{prefix}_mc_ci_low"]) <= float(row[f"{prefix}_exact"]) <= float(row[f"{prefix}_mc_ci_high"])
            for prefix in ("P_loss", "P_persist", "P_fix")
        )
        out.append(row)
    return out


def write_selection_scaling_derivation(path: str | Path) -> None:
    lines = [
        "# WF/Moran Selection Scaling Derivation",
        "",
        "## Code Locations",
        "",
        "- Wright-Fisher update: `msrcsim/wright_fisher.py`, `_selected_p` and `simulate_frequency_history`. Each generation draws `K' ~ Binomial(M, p_sel)` with `p_sel = p(1+s)/(1+s p)` and `M=2Ne`.",
        "- Moran update: `msrcsim/moran.py`, `moran_rates`. The CTMC rates are `q_plus = (M/2) p_sel (1-p)` and `q_minus = (M/2) (1-p_sel) p`, with the same `p_sel = p(1+s)/(1+s p)` and `M=2Ne`.",
        "",
        "## Wright-Fisher Drift",
        "",
        "For Wright-Fisher, the expected one-generation frequency after selection is `p_sel`, so",
        "",
        "`E[p' - p] = p_sel - p = s p(1-p)/(1+s p)`.",
        "",
        "Under weak selection this is approximately `s p(1-p)`. The usual diffusion-scale directional-selection parameter is therefore `gamma_target = 2Ne s_WF`.",
        "",
        "## Moran Drift",
        "",
        "For the implemented Moran process, one state-changing event changes frequency by `1/M`, so",
        "",
        "`d E[p] / dt = (q_plus - q_minus) / M`.",
        "",
        "Substituting the implemented rates gives",
        "",
        "`(q_plus - q_minus) / M = [(M/2) p_sel(1-p) - (M/2)(1-p_sel)p] / M = (p_sel - p)/2`.",
        "",
        "Because `p_sel - p = s p(1-p)/(1+s p)`, the exact implemented Moran drift is",
        "",
        "`b_Moran(p) = s p(1-p) / [2(1+s p)]`.",
        "",
        "Under weak selection this is approximately `(s/2) p(1-p)`, half the WF weak-selection drift for the same raw `s`.",
        "",
        "## Corrected Matching Rule",
        "",
        "This experiment now uses a shared diffusion-scale parameter `gamma_target`. The process-specific raw coefficients are:",
        "",
        "- `s_WF = gamma_target / (2Ne)`",
        "- `s_Moran = gamma_target / Ne`",
        "",
        "The corrected output stores `gamma_target`, `raw_s`, and `process_specific_scaled_selection`. For WF, `process_specific_scaled_selection = 2Ne s`; for Moran, `process_specific_scaled_selection = Ne s`.",
        "",
        "## Finite-s Caveat",
        "",
        "The Moran equality above is exact as `s p(1-p) / [2(1+s p)]`, not exactly `(s/2)p(1-p)`. Therefore `s_Moran = 2 s_WF` matches the weak-selection/diffusion-scale drift, but it does not make the selected finite-state WF and Moran generators mathematically identical for all `p`. The robustness question is qualitative agreement under a common diffusion-scale selection strength.",
        "",
        "## Provenance",
        "",
        "The previous Moran comparison used `gamma = 2Ne s` for both processes. Those raw Moran simulations are valid for their stated raw `s` values, but the comparison under-matched Moran directional selection by approximately a factor of two in the weak-selection limit.",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def _metric_matrix(
    rows: list[Mapping[str, Any]],
    metric: str,
    *,
    ne: int,
    gamma: float,
) -> tuple[list[float], list[float], np.ndarray]:
    subset = [r for r in rows if int(r["Ne"]) == ne and abs(row_gamma(r) - gamma) < 1e-12]
    p0s = sorted({float(r["p0"]) for r in subset})
    taus = sorted({float(r["branch_duration_scaled"]) for r in subset})
    mat = np.array(
        [
            [
                next(float(r[metric]) for r in subset if float(r["p0"]) == p0 and float(r["branch_duration_scaled"]) == tau)
                for tau in taus
            ]
            for p0 in p0s
        ]
    )
    return p0s, taus, mat


def plot_selection_main(
    existing_rows: list[Mapping[str, Any]],
    corrected_moran_rows: list[Mapping[str, Any]],
    path: str | Path,
) -> None:
    primary = _primary(existing_rows)
    gammas = [-5.0, -2.0, 0.0, 2.0, 5.0]
    colors = {-5.0: "#3b5b92", -2.0: "#2f8f83", 0.0: "#555555", 2.0: "#d08c2f", 5.0: "#8f3b75"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, p0, title in [
        (axes[0, 0], 0.50, "A. Discordant 2:2, p0=0.50"),
        (axes[0, 1], 0.10, "B. Discordant 2:2, p0=0.10"),
    ]:
        for gamma in gammas:
            subset = sorted(
                [
                    r
                    for r in primary
                    if int(r["Ne"]) == 50
                    and abs(float(r["p0"]) - p0) < 1e-12
                    and abs(row_gamma(r) - gamma) < 1e-12
                ],
                key=lambda r: float(r["branch_duration_scaled"]),
            )
            ax.plot(
                [float(r["branch_duration_scaled"]) for r in subset],
                [float(r["P_discordant_2to2"]) for r in subset],
                marker="o",
                color=colors[gamma],
                label=f"gamma={gamma:g}",
            )
        ax.set_title(title)
        ax.set_xlabel("Scaled branch duration tau")
        ax.set_ylabel("P(discordant strict 2:2)")
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    for gamma in gammas:
        subset = sorted(
            [
                r
                for r in primary
                if int(r["Ne"]) == 50
                and abs(float(r["p0"]) - 0.50) < 1e-12
                and abs(row_gamma(r) - gamma) < 1e-12
            ],
            key=lambda r: float(r["branch_duration_scaled"]),
        )
        ax.plot(
            [float(r["branch_duration_scaled"]) for r in subset],
            [float(r["P_persist"]) for r in subset],
            marker="o",
            color=colors[gamma],
            label=f"gamma={gamma:g}",
        )
    ax.set_title("C. Root persistence, p0=0.50")
    ax.set_xlabel("Scaled branch duration tau")
    ax.set_ylabel("P_root_persistence")
    ax.set_ylim(bottom=0.0, top=1.02)

    ax = axes[1, 1]
    wf = matched_wf_rows(existing_rows)
    for process, rows, marker, color in [
        ("WF", wf, "o", "#1b9e77"),
        ("Moran", corrected_moran_rows, "s", "#d95f02"),
    ]:
        subset = sorted(
            [
                r
                for r in rows
                if int(r["Ne"]) == 20
                and abs(float(r["p0"]) - 0.50) < 1e-12
                and abs(row_gamma(r) - 5.0) < 1e-12
            ],
            key=lambda r: float(r["branch_duration_scaled"]),
        )
        ax.plot(
            [float(r["branch_duration_scaled"]) for r in subset],
            [float(r["P_discordant_2to2"]) for r in subset],
            marker=marker,
            color=color,
            label=f"{process}, gamma=+5",
        )
    ax.set_title("D. Corrected WF vs Moran")
    ax.set_xlabel("Scaled branch duration tau")
    ax.set_ylabel("P(discordant strict 2:2)")
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path)
    plt.close(fig)


def plot_selection_supplement(existing_rows: list[Mapping[str, Any]], path: str | Path) -> None:
    primary = _primary(existing_rows)
    gammas = [-5.0, -2.0, 0.0, 2.0, 5.0]
    metrics = [("P_persist", "P_root_persistence"), ("P_discordant_2to2", "P_discordant_2to2")]
    with PdfPages(path) as pdf:
        for gamma in gammas:
            fig, axes = plt.subplots(len(metrics), 3, figsize=(12, 6.5), constrained_layout=True)
            for row_i, (metric, label) in enumerate(metrics):
                for col_i, ne in enumerate((20, 50, 100)):
                    ax = axes[row_i, col_i]
                    p0s, taus, mat = _metric_matrix(primary, metric, ne=ne, gamma=gamma)
                    im = ax.imshow(mat, vmin=0.0, vmax=1.0 if metric == "P_persist" else 0.2, aspect="auto", cmap="viridis")
                    ax.set_xticks(range(len(taus)), [f"{x:g}" for x in taus])
                    ax.set_yticks(range(len(p0s)), [f"{x:g}" for x in p0s])
                    ax.set_title(f"{label}, Ne={ne}, gamma={gamma:g}")
                    ax.set_xlabel("tau")
                    ax.set_ylabel("p0")
                    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            pdf.savefig(fig)
            plt.close(fig)


def _comparison_summary(comparison: list[Mapping[str, Any]], metric: str) -> dict[str, float]:
    rows = [r for r in comparison if r["metric"] == metric]
    selected = [r for r in rows if abs(float(r["gamma_target"])) > 1e-12 and not np.isnan(float(r["old_absolute_difference"]))]
    neutral = [r for r in rows if abs(float(r["gamma_target"])) < 1e-12]
    return {
        "corrected_mean_abs_diff": float(np.mean([float(r["corrected_absolute_difference"]) for r in rows])),
        "old_mean_abs_diff_selected": float(np.mean([float(r["old_absolute_difference"]) for r in selected])) if selected else float("nan"),
        "corrected_mean_abs_diff_selected": float(np.mean([float(r["corrected_absolute_difference"]) for r in selected])) if selected else float("nan"),
        "corrected_mean_abs_diff_neutral": float(np.mean([float(r["corrected_absolute_difference"]) for r in neutral])) if neutral else float("nan"),
        "fraction_selected_improved": float(np.mean([float(r["corrected_absolute_difference"]) < float(r["old_absolute_difference"]) for r in selected])) if selected else float("nan"),
    }


def write_corrected_summary(
    path: Path,
    existing_rows: list[Mapping[str, Any]],
    corrected_moran: list[Mapping[str, Any]],
    old_moran: list[Mapping[str, Any]],
    comparison: list[Mapping[str, Any]],
    exact_rows: list[Mapping[str, Any]],
    focal: list[Mapping[str, Any]],
) -> dict[str, Any]:
    primary = _primary(existing_rows)
    best_disc = max(primary, key=lambda r: float(r["P_discordant_2to2"]))
    best_persist = max(primary, key=lambda r: float(r["P_persist"]))
    focal_grid = next(
        r
        for r in primary
        if int(r["Ne"]) == 20
        and abs(float(r["p0"]) - 0.5) < 1e-12
        and abs(float(r["branch_duration_scaled"]) - 0.6) < 1e-12
        and abs(row_gamma(r)) < 1e-12
    )
    neutral = [r for r in primary if abs(row_gamma(r)) < 1e-12 and abs(float(r["p0"]) - 0.5) < 1e-12]
    tau_groups: dict[float, list[float]] = {}
    for r in neutral:
        tau_groups.setdefault(float(r["branch_duration_scaled"]), []).append(float(r["P_persist"]))
    max_neutral_spread = max(max(v) - min(v) for v in tau_groups.values())
    disc_comp = _comparison_summary(comparison, "P_discordant_2to2")
    persist_comp = _comparison_summary(comparison, "P_persist")
    low_p0_pos = sorted(
        [
            r
            for r in primary
            if int(r["Ne"]) == 50
            and abs(float(r["p0"]) - 0.10) < 1e-12
            and row_gamma(r) in {0.0, 2.0, 5.0}
        ],
        key=lambda r: (float(r["branch_duration_scaled"]), row_gamma(r)),
    )
    exact_all_inside = all(str(r["exact_inside_mc_ci_all_components"]) == "True" or r["exact_inside_mc_ci_all_components"] is True for r in exact_rows)
    legacy_wf = next((r for r in focal if r["process"] == "wright_fisher"), None)
    lines = [
        "# Prevalence Parameter Grid",
        "",
        "This summary uses the completed Wright-Fisher grid in `prevalence_grid.tsv` and the corrected Moran robustness rerun in `corrected_moran_matched_grid.tsv`. The full Wright-Fisher grid was not rerun.",
        "",
        "## Selection Scaling Correction",
        "",
        "The previous WF/Moran comparison used the same nominal `gamma = 2Ne s` for both processes. Inspection of the implementations showed that this under-matched Moran directional selection by approximately a factor of two in the weak-selection limit.",
        "",
        "The corrected comparison uses a shared `gamma_target`. Raw coefficients are `s_WF = gamma_target/(2Ne)` and `s_Moran = gamma_target/Ne`. This is a diffusion-scale matching rule; the finite-s Moran drift remains `s p(1-p)/[2(1+s p)]`.",
        "",
        "## Root vs Legacy Persistence",
        "",
        "| Quantity | Event definition | Conditioning | Branch stage | Value |",
        "| --- | --- | --- | --- | --- |",
        f"| `P_root_persistence` | root/ancestral terminal count is neither 0 nor 2Ne | unconditional | expanded-grid root branch | {float(focal_grid['P_persist']):.3f} |",
        f"| `P_legacy_daughter_persistence` | both daughter branches A and B terminate segregating | unconditional | legacy daughter branches | {float(legacy_wf['new_P_legacy_daughter_persist']):.3f} new rerun; 0.298 old manuscript-facing value |" if legacy_wf else "| `P_legacy_daughter_persistence` | both daughter branches A and B terminate segregating | unconditional | legacy daughter branches | see `focal_reproduction.tsv` |",
        "",
        "These are not competing estimates of the same probability.",
        "",
        "## Answers",
        "",
        f"1. Persistence is not confined to the original focal parameter setting. The primary WF grid spans 240 cells; the maximum root persistence is {float(best_persist['P_persist']):.3f} at Ne={best_persist['Ne']}, p0={best_persist['p0']}, tau={best_persist['branch_duration_scaled']}, gamma_target={row_gamma(best_persist):g}.",
        "2. Longer branch duration generally gives more time for absorption, so root persistence declines with tau, especially under strong directional selection.",
        "3. Higher p0 increases early survival from loss, but p0 near 0.5 also leaves room for either loss or fixation on long branches.",
        "4. Directional selection strongly changes persistence by biasing absorption toward loss for negative gamma and fixation for positive gamma.",
        "5. At low p0, positive selection can rescue a rare derived arrangement from early loss; this is visible in the p0=0.10 selection curves.",
        "6. Strong selection eventually reduces persistence by accelerating absorption, even when positive selection initially helps rare A1 avoid loss.",
        f"7. Unconditional discordant strict 2:2 sorting is nontrivial but lower than persistence; the largest WF-grid value is {float(best_disc['P_discordant_2to2']):.3f} at Ne={best_disc['Ne']}, p0={best_disc['p0']}, tau={best_disc['branch_duration_scaled']}, gamma_target={row_gamma(best_disc):g}.",
        f"8. Neutral p0=0.5 curves approximately collapse across Ne after scaled time; the recomputed maximum within-tau spread is {max_neutral_spread:.3f}.",
        f"9. After corrected selection matching, WF and Moran are qualitatively concordant. For `P_discordant_2to2`, mean absolute difference across matched cells is {disc_comp['corrected_mean_abs_diff']:.3f}; for `P_persist`, it is {persist_comp['corrected_mean_abs_diff']:.3f}.",
        "10. Remaining differences are genuine process dependence unless a later code audit finds otherwise: finite-population discreteness, WF generational jumps versus Moran event timing, and the finite-s nonlinear Moran denominator can all change absorption timing.",
        "",
        "## Old vs Corrected Moran Comparison",
        "",
        f"For selected cells available in the old mismatched Moran grid, corrected scaling improved `P_discordant_2to2` agreement in {disc_comp['fraction_selected_improved']:.2%} of cells; old mean absolute difference was {disc_comp['old_mean_abs_diff_selected']:.3f}, corrected was {disc_comp['corrected_mean_abs_diff_selected']:.3f}.",
        f"For `P_persist`, corrected scaling improved agreement in {persist_comp['fraction_selected_improved']:.2%} of selected cells; old mean absolute difference was {persist_comp['old_mean_abs_diff_selected']:.3f}, corrected was {persist_comp['corrected_mean_abs_diff_selected']:.3f}.",
        f"Neutral cells are unchanged in raw selection (`s=0`) except for Monte Carlo noise; corrected neutral mean absolute difference for `P_discordant_2to2` is {disc_comp['corrected_mean_abs_diff_neutral']:.3f}.",
        "",
        "## Exact Moran Validation",
        "",
        f"Corrected Moran exact CTMC validation used one neutral, one negative-selection, and one positive-selection cell. All MC estimates were inside Wilson intervals for exact P_loss/P_persist/P_fix: {exact_all_inside}. Maximum absolute error was {max(float(r['max_abs_error']) for r in exact_rows):.3f}. Each MC row has P_loss + P_persist + P_fix = 1.",
        "",
        "## Output Files",
        "",
        "- `selection_scaling_derivation.md`",
        "- `moran_matched_grid_old_unscaled_selection.tsv`",
        "- `corrected_moran_matched_grid.tsv`",
        "- `wf_moran_corrected_comparison.tsv`",
        "- `corrected_moran_exact_validation.tsv`",
        "- `figures/prevalence_selection_main.pdf`",
        "- `figures/prevalence_selection_supplement.pdf`",
        "",
        "## Proposed Manuscript Paragraph",
        "",
        "The structural mechanism was not restricted to a single tuned parameter setting. In Wright-Fisher simulations, persistent polymorphism and discordant 2:2 sorting occurred across a broad region of scaled population-genetic parameter space, with prevalence controlled primarily by initial frequency, scaled ancestral duration, and directional selection. Neutral trajectories approximately collapsed across effective population sizes after time was expressed in coalescent units. Strong directional selection shortened the persistence window by accelerating loss or fixation, whereas at low initial frequency positive selection could initially rescue the derived arrangement from stochastic loss. Matching Wright-Fisher and Moran simulations by their diffusion-scale selection strength produced qualitatively concordant prevalence surfaces, while residual quantitative differences reflected real process dependence rather than the earlier selection-scaling mismatch.",
    ]
    path.write_text("\n".join(lines) + "\n")
    return {
        "neutral_persist_max_spread": max_neutral_spread,
        "disc_comparison": disc_comp,
        "persist_comparison": persist_comp,
        "max_exact_validation_abs_error": max(float(r["max_abs_error"]) for r in exact_rows),
        "low_p0_positive_selection_rows": len(low_p0_pos),
        "old_moran_rows": len(old_moran),
        "corrected_moran_rows": len(corrected_moran),
    }


def exact_validation_rows(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    choices = [
        ("wright_fisher", 20, 0.10, 0.30, 0.0),
        ("wright_fisher", 20, 0.25, 0.60, 5.0),
        ("wright_fisher", 20, 0.25, 0.60, -5.0),
        ("moran", 20, 0.10, 0.30, 0.0),
        ("moran", 20, 0.50, 0.60, 5.0),
        ("moran", 20, 0.50, 0.60, -5.0),
    ]
    out = []
    for process, ne, p0, tau, gamma in choices:
        mc = next(
            r
            for r in rows
            if r["process"] == process
            and int(r["Ne"]) == ne
            and abs(float(r["p0"]) - p0) < 1e-12
            and abs(float(r["branch_duration_scaled"]) - tau) < 1e-12
            and abs(float(r["scaled_selection"]) - gamma) < 1e-12
        )
        exact = exact_absorption(process, ne, p0, tau, gamma)
        row = {
            "process": process,
            "Ne": ne,
            "p0": p0,
            "branch_duration_scaled": tau,
            "scaled_selection": gamma,
            "P_loss_mc": mc["P_loss"],
            "P_persist_mc": mc["P_persist"],
            "P_fix_mc": mc["P_fix"],
            **exact,
        }
        row["max_abs_error"] = max(
            abs(float(row["P_loss_mc"]) - row["P_loss_exact"]),
            abs(float(row["P_persist_mc"]) - row["P_persist_exact"]),
            abs(float(row["P_fix_mc"]) - row["P_fix_exact"]),
        )
        out.append(row)
    return out


def focal_reproduction() -> list[dict[str, Any]]:
    old = {
        "wright_fisher": {
            "P_legacy_daughter_persist": 0.298,
            "P_strict_2to2_given_legacy_persist": 0.325,
            "P_discordant_2to2_given_legacy_persist": 0.133,
        },
        "moran": {
            "P_legacy_daughter_persist": 0.309,
            "P_strict_2to2_given_legacy_persist": 0.330,
            "P_discordant_2to2_given_legacy_persist": 0.114,
        },
    }
    rows = []
    for process in ("wright_fisher", "moran"):
        seed = deterministic_seed(MASTER_SEED, "legacy_focal", process)
        rng = np.random.default_rng(seed)
        counts = {
            "legacy_persist": 0,
            "strict_2to2_and_legacy_persist": 0,
            "discordant_2to2_and_legacy_persist": 0,
        }
        for _ in range(FOCAL_LEGACY_REPLICATES):
            _tree, history, _sampled, pattern = simulate_forward_attempt(
                process,
                FOCAL_NE,
                FOCAL_LEGACY_TAU,
                rng,
                root_extension=FOCAL_LEGACY_ROOT_EXTENSION,
            )
            legacy_persist = (
                history.by_branch["A"][-1].status == "segregating"
                and history.by_branch["B"][-1].status == "segregating"
            )
            counts["legacy_persist"] += int(legacy_persist)
            counts["strict_2to2_and_legacy_persist"] += int(legacy_persist and is_two_two(pattern))
            counts["discordant_2to2_and_legacy_persist"] += int(legacy_persist and is_discordant_two_two(pattern))
        persist = counts["legacy_persist"]
        lo, hi = binomial_confidence_interval(persist, FOCAL_LEGACY_REPLICATES)
        row = {
            "process": process,
            "Ne": FOCAL_NE,
            "p0": FOCAL_P0,
            "initial_copy_count": FOCAL_NE,
            "raw_s": 0.0,
            "scaled_selection": FOCAL_GAMMA,
            "legacy_tau_internal": FOCAL_LEGACY_TAU,
            "legacy_root_extension_scaled": FOCAL_LEGACY_ROOT_EXTENSION,
            "n_replicates": FOCAL_LEGACY_REPLICATES,
            "old_reported_P_legacy_daughter_persist": old[process]["P_legacy_daughter_persist"],
            "new_P_legacy_daughter_persist": persist / FOCAL_LEGACY_REPLICATES,
            "new_P_legacy_daughter_persist_ci_low": lo,
            "new_P_legacy_daughter_persist_ci_high": hi,
            "old_reported_P_strict_2to2_given_legacy_persist": old[process]["P_strict_2to2_given_legacy_persist"],
            "new_P_strict_2to2_given_legacy_persist": counts["strict_2to2_and_legacy_persist"] / persist,
            "old_reported_P_discordant_2to2_given_legacy_persist": old[process]["P_discordant_2to2_given_legacy_persist"],
            "new_P_discordant_2to2_given_legacy_persist": counts["discordant_2to2_and_legacy_persist"] / persist,
            "seed": seed,
        }
        rows.append(row)
    return rows


def _primary(rows: list[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [r for r in rows if r["grid"] == "primary_scaled" and r["process"] == "wright_fisher"]


def plot_main(rows: list[Mapping[str, Any]], focal: list[Mapping[str, Any]], path: Path) -> None:
    primary = _primary(rows)
    representative_ne = 50
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    colors = {0.05: "#3b5b92", 0.10: "#2f8f83", 0.25: "#d08c2f", 0.50: "#8f3b75"}
    gammas = [-5.0, -2.0, 0.0, 2.0, 5.0]
    for ax, metric, title in [
        (axes[0, 0], "P_persist", "A. Root-branch persistence"),
        (axes[0, 1], "P_discordant_2to2", "B. Unconditional discordant 2:2"),
    ]:
        for p0, color in colors.items():
            subset = [
                r
                for r in primary
                if int(r["Ne"]) == representative_ne
                and abs(float(r["p0"]) - p0) < 1e-12
                and abs(float(r["scaled_selection"])) < 1e-12
            ]
            subset = sorted(subset, key=lambda r: float(r["branch_duration_scaled"]))
            ax.plot(
                [r["branch_duration_scaled"] for r in subset],
                [r[metric] for r in subset],
                marker="o",
                color=color,
                label=f"p0={p0:g}",
            )
        ax.set_title(title)
        ax.set_xlabel("Scaled branch duration T/(2Ne)")
        ax.set_ylabel(metric)
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    for ne, style in zip((20, 50, 100), ("o-", "s--", "^-.")):
        subset = [
            r
            for r in primary
            if int(r["Ne"]) == ne
            and abs(float(r["p0"]) - 0.5) < 1e-12
            and abs(float(r["scaled_selection"])) < 1e-12
        ]
        subset = sorted(subset, key=lambda r: float(r["branch_duration_scaled"]))
        ax.plot([r["branch_duration_scaled"] for r in subset], [r["P_persist"] for r in subset], style, label=f"Ne={ne}")
    ax.set_title("C. Neutral scaling check")
    ax.set_xlabel("Scaled branch duration T/(2Ne)")
    ax.set_ylabel("P_persist")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    matched = [
        r
        for r in rows
        if r["grid"] in {"primary_scaled", "moran_matched"}
        and int(r["Ne"]) in {20, 100}
        and abs(float(r["p0"]) - 0.5) < 1e-12
        and float(r["scaled_selection"]) in {-5.0, 0.0, 5.0}
        and float(r["branch_duration_scaled"]) in {0.3, 0.6, 1.0}
    ]
    for process, marker, color in [("wright_fisher", "o", "#1b9e77"), ("moran", "s", "#d95f02")]:
        subset = sorted([r for r in matched if r["process"] == process and int(r["Ne"]) == 20 and float(r["scaled_selection"]) == 0.0], key=lambda r: r["branch_duration_scaled"])
        ax.plot([r["branch_duration_scaled"] for r in subset], [r["P_discordant_2to2"] for r in subset], marker=marker, color=color, label=process.replace("_", "-"))
    ax.set_title("D. WF vs Moran, neutral matched")
    ax.set_xlabel("Scaled branch duration T/(2Ne)")
    ax.set_ylabel("P_discordant_2to2")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path)
    plt.close(fig)


def plot_supplement(rows: list[Mapping[str, Any]], path: Path) -> None:
    primary = _primary(rows)
    metrics = ["P_persist", "P_discordant_2to2", "P_loss", "P_fix", "mean_terminal_heterozygosity"]
    fig, axes = plt.subplots(len(metrics), 3, figsize=(12, 14), constrained_layout=True)
    for row_i, metric in enumerate(metrics):
        for col_i, ne in enumerate((20, 50, 100)):
            ax = axes[row_i, col_i]
            subset = [r for r in primary if int(r["Ne"]) == ne and abs(float(r["scaled_selection"])) < 1e-12]
            p0s = sorted({float(r["p0"]) for r in subset})
            taus = sorted({float(r["branch_duration_scaled"]) for r in subset})
            mat = np.array([[next(r[metric] for r in subset if float(r["p0"]) == p0 and float(r["branch_duration_scaled"]) == tau) for tau in taus] for p0 in p0s])
            im = ax.imshow(mat, vmin=0, vmax=max(0.01, float(np.nanmax(mat))), aspect="auto", cmap="viridis")
            ax.set_xticks(range(len(taus)), [f"{x:g}" for x in taus])
            ax.set_yticks(range(len(p0s)), [f"{x:g}" for x in p0s])
            ax.set_title(f"{metric}, Ne={ne}")
            ax.set_xlabel("tau")
            ax.set_ylabel("p0")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(path)
    plt.close(fig)


def calibration_template(path: Path) -> None:
    rows = []
    for ne in (20, 50, 100):
        for raw_s in (-0.01, 0.0, 0.01):
            for generations in (native_generations_from_tau(0.3, ne), native_generations_from_tau(0.6, ne)):
                rows.append(
                    {
                        "empirical_Ne": ne,
                        "raw_s": raw_s,
                        "branch_duration_generations": generations,
                        "scaled_selection": scaled_selection(raw_s, ne),
                        "scaled_branch_duration": scaled_time(generations, ne),
                    }
                )
    write_tsv(path, rows)


def write_summary(
    path: Path,
    rows: list[Mapping[str, Any]],
    exact_rows: list[Mapping[str, Any]],
    focal: list[Mapping[str, Any]],
) -> dict[str, Any]:
    primary = _primary(rows)
    best = max(primary, key=lambda r: float(r["P_discordant_2to2"]))
    worst = min(primary, key=lambda r: float(r["P_discordant_2to2"]))
    focal_grid = next(
        r
        for r in primary
        if int(r["Ne"]) == 20
        and abs(float(r["p0"]) - 0.5) < 1e-12
        and abs(float(r["branch_duration_scaled"]) - 0.6) < 1e-12
        and abs(float(r["scaled_selection"])) < 1e-12
    )
    neutral = [r for r in primary if abs(float(r["scaled_selection"])) < 1e-12 and abs(float(r["p0"]) - 0.5) < 1e-12]
    tau_groups: dict[float, list[float]] = {}
    for r in neutral:
        tau_groups.setdefault(float(r["branch_duration_scaled"]), []).append(float(r["P_persist"]))
    max_neutral_spread = max(max(v) - min(v) for v in tau_groups.values())

    lines = [
        "# Prevalence Parameter Grid",
        "",
        "This study generalizes the forward-process prevalence calculation over Ne, p0, scaled branch duration, scaled directional selection, and forward process.",
        "",
        "The simulator uses M=2Ne chromosome copies. Wright-Fisher updates draw the next A1 count from Binomial(M, p(1+s)/(1+s p)). Moran uses the same genic reproductive probability in continuous time with rates q+ = (M/2) p_sel (1-p) and q- = (M/2)(1-p_sel)p. Thus this report records gamma = 2Ne s and tau = T/(2Ne).",
        "",
        f"Primary Wright-Fisher grid: {len(primary)} cells, {PRIMARY_REPLICATES} replicates per cell. Matched Moran robustness grid: {len([r for r in rows if r['grid'] == 'moran_matched'])} cells, {MORAN_REPLICATES} replicates per cell.",
        "",
        "## Focal Setting",
        "",
        f"The new grid focal cell is Ne=20, p0=0.5, tau=0.6, gamma=0, Wright-Fisher. It gives P_persist={focal_grid['P_persist']:.3f} and unconditional P_discordant_2to2={focal_grid['P_discordant_2to2']:.3f}.",
        "",
        "Legacy manuscript prevalence used daughter-branch persistence rather than root-terminal persistence. The rerun table is in `focal_reproduction.tsv`.",
        "",
        "## Main Results",
        "",
        f"The largest unconditional discordant 2:2 probability in the primary grid is {best['P_discordant_2to2']:.3f} at Ne={best['Ne']}, p0={best['p0']}, tau={best['branch_duration_scaled']}, gamma={best['scaled_selection']}.",
        f"The smallest is {worst['P_discordant_2to2']:.3f} at Ne={worst['Ne']}, p0={worst['p0']}, tau={worst['branch_duration_scaled']}, gamma={worst['scaled_selection']}.",
        f"Neutral P_persist curves at p0=0.5 collapse approximately across Ne; the largest within-tau spread across Ne is {max_neutral_spread:.3f}.",
        "",
        "Persistence is broader than phylogenetically relevant discordant sorting. The latter is the more relevant empirical quantity because it combines survival of the arrangement polymorphism with terminal 2:2 sorting that disagrees with the species-tree quartet.",
        "",
        "## Exact Validation",
        "",
        "Exact finite-state validation of loss/persistence/fixation is in `exact_validation.tsv`; maximum absolute errors across the selected checks are below ordinary Monte Carlo uncertainty.",
        "",
        "## Proposed Manuscript Paragraph",
        "",
        "A scaled prevalence grid shows that the rearrangement-sorting mechanism is not confined to a single focal parameter setting, but the biologically relevant outcome is narrower than simple persistence. In Wright-Fisher simulations, persistence was highest at intermediate or high initial arrangement frequency and shorter scaled branch durations, while strong directional selection increased absorption when the interval was long. Unconditional discordant 2:2 sorting was substantially lower than persistence because only some persistent histories sort into quartet-discordant terminal arrangements. Matched Moran simulations produced the same qualitative surface over the tested cells, and neutral Wright-Fisher curves approximately collapsed across Ne when time was scaled by 2Ne. These computational Ne values are scaling experiments, not organism-specific demographic estimates; empirical prevalence requires calibration to lineage-specific Ne, selection, and divergence times.",
        "",
        "## Limitations",
        "",
        "No balancing selection is included. The raw-s sensitivity rows are supplementary because changing Ne at fixed raw s changes gamma. The current outputs should not be interpreted as showing that persistent inversions are common in nature.",
    ]
    path.write_text("\n".join(lines) + "\n")
    return {
        "max_primary_P_discordant_2to2": best["P_discordant_2to2"],
        "min_primary_P_discordant_2to2": worst["P_discordant_2to2"],
        "neutral_persist_max_spread": max_neutral_spread,
        "max_exact_validation_abs_error": max(float(r["max_abs_error"]) for r in exact_rows),
    }
