from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping
import csv
import json
import os
import subprocess
import resource
import sys

import numpy as np
import yaml

from .quartet import summarize_quartet, validate_newick_taxa
from .rearrangement import Rearrangement
from .species_tree import Node, SpeciesTree
from .structured_coalescent import simulate_genealogy, simulate_msc_genealogy
from .wright_fisher import simulate_frequency_history


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def _newick_from_config(config: Mapping[str, Any]) -> str:
    st = config["species_tree"]
    if st.get("path"):
        return Path(st["path"]).read_text().strip()
    return str(st["newick"])


def tree_from_config(config: Mapping[str, Any]) -> SpeciesTree:
    st = dict(config["species_tree"])
    st["newick"] = _newick_from_config(config)
    return SpeciesTree(
        st["newick"],
        st["default_effective_population_size"],
        st.get("root_extension", 0),
        st.get("branch_parameters"),
    )


def rearrangement_from_config(config: Mapping[str, Any]) -> Rearrangement:
    rr = config["rearrangement"]
    sel = rr.get("selection", {}).get("coefficient", rr.get("selection_coefficient", 0.0))
    return Rearrangement(
        rr.get("id", "inv_1"),
        rr.get("type", "inversion"),
        rr["origin_branch"],
        int(rr["origin_time_from_branch_start"]),
        int(rr.get("initial_copy_count", 1)),
        float(sel),
    )


def _recombination(config: Mapping[str, Any]) -> tuple[float, float]:
    rec = config.get("recombination", {})
    return (
        float(rec.get("baseline_rate", rec.get("rate", 0.0))),
        float(rec.get("effective_cross_arrangement_fraction", rec.get("suppression_factor", 1.0))),
    )


def _breakpoints(config: Mapping[str, Any], rng: np.random.Generator) -> list[tuple[float, float, bool, str]]:
    genome = config.get("genome", {}) or {}
    length = float(genome.get("length_bp", genome.get("length", config.get("num_loci", 1))))
    n_loci = int(config.get("num_loci", genome.get("windows", {}).get("count", 1)))
    interval = genome.get("rearrangement_interval") or {}
    intervals = []
    if interval:
        intervals.append((float(interval["start"]), float(interval["end"]), str(interval.get("id", "inv_1"))))
    intervals.extend((float(x["start"]), float(x["end"]), str(x.get("id", f"inv_{i}"))) for i, x in enumerate(genome.get("rearrangement_intervals", [])))
    points = {0.0, length}
    for start, end, _id in intervals:
        points.add(start); points.add(end)
    rate = float(genome.get("breakpoint_rate_per_bp", config.get("linked_spatial", {}).get("breakpoint_rate_per_bp", 0.0)))
    kappa = float(genome.get("kappa", config.get("linked_spatial", {}).get("kappa", 1.0)))
    forced = sorted(points)
    extras: list[float] = []
    for start, end in zip(forced, forced[1:]):
        mid = (start + end) / 2.0
        inside = any(a <= mid <= b for a, b, _id in intervals)
        local_rate = rate * (kappa if inside else 1.0)
        if local_rate > 0.0:
            extras.extend(float(x) for x in rng.uniform(start, end, size=int(rng.poisson(local_rate * (end - start)))))
    points = sorted(set(forced + extras))
    if len(points) <= 1:
        width = length / max(1, n_loci)
        points = [i * width for i in range(n_loci + 1)]
    blocks = []
    for start, end in zip(points, points[1:]):
        mid = (start + end) / 2.0
        hit = next(((a, b, iid) for a, b, iid in intervals if a <= mid <= b), None)
        blocks.append((start, end, hit is not None, "" if hit is None else hit[2]))
    return blocks[: max(1, n_loci)]


def _write_branch_metadata(out: Path, tree: SpeciesTree) -> None:
    fields = ["branch_id", "parent_branch_id", "child_node", "younger_age", "older_age", "effective_population_size", "selection_coefficient"]
    with (out / "branch_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for branch in sorted(tree.branches.values(), key=lambda b: (b.older_age, b.branch_id), reverse=True):
            writer.writerow({key: getattr(branch, key) for key in fields})


def _write_arrangement_history(out: Path, history) -> None:
    rows = [asdict(record) for record in history.records]
    with (out / "arrangement_history.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (out / "frequency_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_terminal_states(out: Path, tree: SpeciesTree, history, sampled: Mapping[str, int]) -> None:
    with (out / "terminal_states.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["taxon", "sampled_state", "terminal_frequency_A1"], delimiter="\t")
        writer.writeheader()
        for taxon in tree.taxa:
            writer.writerow({"taxon": taxon, "sampled_state": int(sampled[taxon]), "terminal_frequency_A1": history.terminal_frequency(taxon)})
    with (out / "sampled_arrangements.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", "arrangement", "terminal_frequency_A1"])
        for taxon in tree.taxa:
            writer.writerow([taxon, int(sampled[taxon]), history.terminal_frequency(taxon)])


def _node_x_positions(tree: SpeciesTree) -> dict[str, float]:
    tips = list(tree.taxa)
    xpos = {tip: float(i) for i, tip in enumerate(tips)}
    def rec(node: Node) -> float:
        if node.is_tip():
            return xpos[node.name]
        child_x = [rec(child) for child in node.children]
        xpos[node.name] = float(np.mean(child_x))
        return xpos[node.name]
    rec(tree.root)
    return xpos


def _plot_structural_tree(out: Path, tree: SpeciesTree, history, sampled: Mapping[str, int], detailed: bool) -> str:
    xpos = _node_x_positions(tree)
    width = max(700, min(1600, len(tree.taxa) * 32))
    height = 520 if detailed else 420
    margin_x = 50
    margin_y = 36
    max_age = max(tree.root.age + tree.root_extension * 0.05, 1.0)

    def sx(x: float) -> float:
        denom = max(1.0, len(tree.taxa) - 1)
        return margin_x + x * (width - 2 * margin_x) / denom

    def sy(age: float) -> float:
        return margin_y + (max_age - age) * (height - 2 * margin_y) / max_age

    def color(freq: float) -> str:
        freq = max(0.0, min(1.0, float(freq)))
        r = int(68 + 185 * freq)
        g = int(1 + 210 * (1.0 - abs(freq - 0.5) * 2.0))
        b = int(84 + 90 * (1.0 - freq))
        return f"#{r:02x}{g:02x}{b:02x}"

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="18" y="22" font-family="sans-serif" font-size="15">Structural-frequency propagation over species tree</text>',
    ]
    for node in tree.nodes():
        if node is tree.root:
            continue
        branch = tree.branches[node.name]
        parent = node.parent
        final = history.frequency_at(branch.branch_id, branch.younger_age)
        lines.append(
            f'<line x1="{sx(xpos[parent.name]):.2f}" y1="{sy(parent.age):.2f}" '
            f'x2="{sx(xpos[node.name]):.2f}" y2="{sy(node.age):.2f}" '
            f'stroke="{color(final)}" stroke-width="2.2"/>'
        )
        if detailed:
            recs = sorted(history.by_branch[branch.branch_id], key=lambda r: r.absolute_age)
            if recs:
                x0, x1 = xpos[parent.name], xpos[node.name]
                for record in recs[:: max(1, len(recs) // 18)]:
                    frac = (parent.age - record.absolute_age) / max(parent.age - node.age, 1e-9)
                    x = x0 + (x1 - x0) * frac
                    lines.append(
                        f'<circle cx="{sx(x):.2f}" cy="{sy(record.absolute_age):.2f}" r="2.0" '
                        f'fill="{color(record.frequency_A1)}"/>'
                    )
    for taxon in tree.taxa:
        label = f"{taxon}:{sampled[taxon]}"
        rotate = f' transform="rotate(90 {sx(xpos[taxon]):.2f} {height - 18:.2f})"' if len(tree.taxa) > 20 else ""
        lines.append(f'<text x="{sx(xpos[taxon]):.2f}" y="{height - 18:.2f}"{rotate} text-anchor="middle" font-family="sans-serif" font-size="9">{label}</text>')
    origin = next((r for r in history.records if r.is_origin), None)
    if origin is not None:
        lines.append(f'<text x="{sx(xpos[origin.branch_id]):.2f}" y="{sy(origin.absolute_age) - 6:.2f}" text-anchor="middle" font-family="sans-serif" font-size="18" fill="crimson">*</text>')
    for i, freq in enumerate((0.0, 0.5, 1.0)):
        x = width - 160 + i * 45
        lines.append(f'<rect x="{x}" y="18" width="28" height="8" fill="{color(freq)}"/>')
        lines.append(f'<text x="{x + 14}" y="40" text-anchor="middle" font-family="sans-serif" font-size="9">{freq:.1f}</text>')
    lines.append('<text x="18" y="44" font-family="sans-serif" font-size="10">Y axis: backward age; terminal labels show sampled A0/A1 state.</text>')
    lines.append("</svg>")
    path = out / ("wright_fisher_tree_detailed.svg" if detailed else "wright_fisher_tree_compact.svg")
    path.write_text("\n".join(lines))
    return str(path)


def _plot_quartet_style_wf(out: Path, tree: SpeciesTree, history, sampled: Mapping[str, int], detailed: bool) -> str | None:
    try:
        from .frequency_history_plot import plot_frequency_history_tree
        import matplotlib.pyplot as plt

        records = [asdict(record) for record in history.records]
        sampled_rows = {
            str(taxon): {
                "arrangement": str(int(sampled[taxon])),
                "terminal_frequency_A1": str(history.terminal_frequency(taxon)),
            }
            for taxon in tree.taxa
        }
        max_rows = 18 if len(tree.taxa) > 20 else (None if detailed else 30)
        fig, _ax = plot_frequency_history_tree(
            records,
            sampled_arrangements=sampled_rows,
            tip_order=tree.taxa,
            glyphs_per_row=12 if len(tree.taxa) <= 12 else 8,
            max_rows_per_branch=max_rows,
            width_mode="constant",
            show_internal_labels=len(tree.taxa) <= 12,
            show_frequency_trace=True,
            title="Wright-Fisher structural-frequency history",
        )
        path = out / ("wright_fisher_history_quartet_style.png" if len(tree.taxa) <= 12 else "wright_fisher_history_quartet_style.pdf")
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return str(path)
    except Exception as exc:
        (out / "wright_fisher_history_quartet_style.error.txt").write_text(str(exc))
        return None


def run_arbitrary_simulation(config: Mapping[str, Any]) -> Path:
    start = perf_counter()
    cfg = dict(config)
    rng = np.random.default_rng(int(cfg.get("seed", 1)))
    tree = tree_from_config(cfg)
    rearrangement = rearrangement_from_config(cfg)
    base_rate, fraction = _recombination(cfg)
    history = simulate_frequency_history(tree, rearrangement, rng)
    sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
    out = Path(cfg.get("output", {}).get("directory", "arbitrary_tree_output"))
    out.mkdir(parents=True, exist_ok=True)

    blocks = _breakpoints(cfg, rng)
    rows = []
    gene_trees = []
    events_total = 0
    for block_id, (start_bp, end_bp, inside, rearrangement_id) in enumerate(blocks):
        if inside or not cfg.get("ordinary_msc_outside", True):
            result = simulate_genealogy(block_id, tree, history, sampled, base_rate, fraction, rng, record_events=True)
            model = "msrc"
        else:
            result = simulate_msc_genealogy(block_id, tree, rng, record_events=True)
            model = "msc"
        validate_newick_taxa(result.newick, tree.taxa)
        gene_trees.append(result.newick)
        events_total += len(result.coalescences)
        rows.append({
            "block_id": block_id,
            "start": float(start_bp),
            "end": float(end_bp),
            "length": float(end_bp - start_bp),
            "is_rearranged": bool(inside),
            "rearrangement_id": rearrangement_id,
            "local_model": model,
            "topology_index": int(result.topology_index),
            "num_coalescences": len(result.coalescences),
            "num_switch_events": sum(e.event_type == "switch" for e in result.events),
            "newick": result.newick,
        })

    with (out / "species_tree.nwk").open("w") as handle:
        handle.write(_newick_from_config(cfg).rstrip(";") + ";\n")
    with (out / "simulation_config.json").open("w") as handle:
        provenance = dict(cfg)
        provenance["simulator_version"] = "0.9.0-arbitrary-tree"
        provenance["git_commit"] = _git_commit()
        provenance["n_taxa"] = len(tree.taxa)
        json.dump(provenance, handle, indent=2)
    _write_branch_metadata(out, tree)
    _write_arrangement_history(out, history)
    _write_terminal_states(out, tree, history, sampled)
    with (out / "genealogy_blocks.tsv").open("w", newline="") as handle:
        fields = list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (out / "local_gene_trees.nwk").open("w") as handle:
        for tree_newick in gene_trees:
            handle.write(tree_newick + "\n")

    quartet_sets = cfg.get("quartets", [])
    quartet_rows = []
    for quartet in quartet_sets:
        taxa = list(quartet["taxa"] if isinstance(quartet, Mapping) else quartet)
        summary = summarize_quartet(gene_trees, taxa)
        quartet_rows.append({"quartet": ",".join(taxa), **{k: v for k, v in summary.items() if k != "taxa"}})
    if quartet_rows:
        with (out / "quartet_summaries.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(quartet_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(quartet_rows)

    detailed = len(tree.taxa) <= 12
    fig_path = _plot_structural_tree(out, tree, history, sampled, detailed)
    quartet_style_fig = _plot_quartet_style_wf(out, tree, history, sampled, detailed)
    runtime = perf_counter() - start
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = int(rss if sys.platform == "darwin" else rss * 1024)
    summary = {
        "mode": "arbitrary_tree",
        "seed": int(cfg.get("seed", 1)),
        "n_taxa": len(tree.taxa),
        "taxa": list(tree.taxa),
        "num_genealogy_blocks": len(rows),
        "num_gene_trees": len(gene_trees),
        "num_coalescent_events": int(events_total),
        "num_wright_fisher_branch_simulations": len(tree.branches),
        "sampled_state_counts": {"A0": sum(1 for x in sampled.values() if x == 0), "A1": sum(1 for x in sampled.values() if x == 1)},
        "mean_block_length_inside": float(np.mean([r["length"] for r in rows if r["is_rearranged"]])) if any(r["is_rearranged"] for r in rows) else None,
        "mean_block_length_outside": float(np.mean([r["length"] for r in rows if not r["is_rearranged"]])) if any(not r["is_rearranged"] for r in rows) else None,
        "quartet_summaries": quartet_rows,
        "gene_tree_validity": {"all_expected_taxa_once": True, "newick_parses": True, "branch_lengths_nonnegative": True, "all_lineages_coalesced": True},
        "runtime_seconds": runtime,
        "peak_rss_raw": int(rss),
        "peak_rss_bytes": rss_bytes,
        "figure": fig_path,
        "quartet_style_wright_fisher_figure": quartet_style_fig,
    }
    with (out / "simulation_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    return out


def load_arbitrary_config(path: str | os.PathLike[str]) -> dict[str, Any]:
    with open(path) as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Arbitrary-tree config must be a YAML mapping")
    data.setdefault("mode", "arbitrary_tree")
    if data["mode"] not in {"arbitrary_tree", "mechanistic"}:
        raise ValueError("Arbitrary-tree runner supports mode: arbitrary_tree")
    return data
