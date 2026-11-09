from __future__ import annotations

from dataclasses import dataclass, field
import csv
import json
from pathlib import Path
import time

import numpy as np

from ..species_tree import SpeciesTree
from .backgrounds.neutral import NeutralBackground
from .backgrounds.msrc_two_state import MSRCTwoBackground
from .coalescence import union_lineage
from .events import ARGEvent
from .lineages import ARGLineage
from .recombination import recombine
from .recombination_map import PiecewiseConstantRecombinationMap
from .segments import Segment, interval_union, split, total_length


@dataclass(frozen=True)
class ARGConfig:
    sequence_length_bp: float
    recombination_map: PiecewiseConstantRecombinationMap
    background_model: str | object = "neutral"
    record_events: bool = True
    background_history: object | None = None
    sampled_backgrounds: dict[str, int] | None = None


@dataclass
class ARGResult:
    sequence_length: float
    nodes: list[dict]
    edges: list[tuple[float, float, int, int]]
    events: list[ARGEvent]
    runtime_seconds: float = 0.0
    peak_lineages: int = 0
    peak_segments: int = 0

    @property
    def marginal_intervals(self):
        points = {0.0, self.sequence_length}
        points.update(e[0] for e in self.edges)
        points.update(e[1] for e in self.edges)
        return tuple((a, b) for a, b in zip(sorted(points), sorted(points)[1:]) if b > a)


def _add_edge(edges, child, parent, left, right):
    if right - left > 1e-12:
        edges.append((float(left), float(right), int(parent), int(child)))


def simulate_arg(tree: SpeciesTree, config: ARGConfig, rng: np.random.Generator) -> ARGResult:
    """Simulate a neutral, continuous-time Hudson ARG on a species tree.

    Recombination has hazard integral(r(x) dx) over each lineage's current
    ancestral material. Pairwise coalescence has rate 1/(2N) for every pair in
    the same species population, including pairs with disjoint material.
    """
    if isinstance(config.background_model, str):
        if config.background_model != "neutral":
            raise ValueError("unknown background model; pass an MSRC background object")
        background = NeutralBackground()
    else:
        background = config.background_model
    if abs(config.recombination_map.length - config.sequence_length_bp) > 1e-9:
        raise ValueError("recombination map length must equal sequence_length_bp")
    start = time.perf_counter()
    nodes: list[dict] = []
    edges: list[tuple[float, float, int, int]] = []
    lineages: list[ARGLineage] = []
    next_node = 0
    next_lineage = 0
    for taxon in tree.taxa:
        node = next_node
        next_node += 1
        nodes.append({"flags": 1, "time": 0.0, "metadata": {"type": "sample", "taxon": taxon}})
        state = background.initial_state(tree.branch_for_tip(taxon), (0.0, config.sequence_length_bp), taxon=taxon)
        if config.sampled_backgrounds is not None and taxon in config.sampled_backgrounds:
            state = int(config.sampled_backgrounds[taxon])
        lineages.append(ARGLineage(next_lineage, tree.branch_for_tip(taxon),
                                   (Segment(0.0, config.sequence_length_bp, node),), 0.0, state))
        next_lineage += 1
    events: list[ARGEvent] = []
    event_index = 0
    peak = len(lineages)
    peak_segments = sum(len(l.segments) for l in lineages)
    t = 0.0

    def add_node(event_type, t, **metadata):
        nonlocal next_node
        node_id = next_node
        next_node += 1
        nodes.append({"flags": 0, "time": float(t), "metadata": {"type": event_type, **metadata}})
        return node_id

    def refresh_completed_material():
        """Retire intervals with one remaining active carrier."""
        boundaries = {x for lineage in lineages for interval in lineage.active for x in interval}
        if len(boundaries) < 2:
            return
        replacements = []
        for lineage in lineages:
            keep = []
            for left, right in zip(sorted(boundaries), sorted(boundaries)[1:]):
                if not any(a <= left and right <= b for a, b in lineage.active):
                    continue
                carriers = sum(any(a <= left and right <= b for a, b in other.active) for other in lineages)
                if carriers > 1:
                    keep.append((left, right))
            replacements.append(lineage.with_active(keep))
        lineages[:] = [lineage for lineage in replacements if lineage.active]

    while len(lineages) > 1:
        # Boundary transitions are deterministic and must happen before events.
        moved = False
        for lineage in lineages:
            branch = tree.branches[lineage.population]
            if branch.parent_branch_id is not None and t >= branch.older_age - 1e-10:
                old = lineage.population
                lineage.population = branch.parent_branch_id
                moved = True
                if config.record_events:
                    events.append(ARGEvent(event_index, t, "species_boundary", (lineage.lineage_id,), (old,), (lineage.population,)))
                    event_index += 1
        if moved:
            continue

        options = []
        for i, a in enumerate(lineages):
            if hasattr(background, "recombination_hazard"):
                hazard = background.recombination_hazard(a.population, a.background, a.active,
                                                         a.segments, t, config.recombination_map)
            else:
                hazard, _ = config.recombination_map.crossover_hazard_and_region(a.active)
            if hazard > 0:
                options.append(("recombination", (i,), hazard))
            for j in range(i):
                b = lineages[j]
                if a.population != b.population:
                    continue
                rate = background.coalescence_rate(a.population, a.background, b.background,
                                                   tree.branches[a.population].effective_population_size, t)
                if rate < 0:
                    raise ValueError("coalescence rates must be nonnegative")
                if rate > 0:
                    options.append(("coalescence", (j, i), rate))
        total = sum(x[2] for x in options)
        boundaries = [tree.branches[l.population].older_age for l in lineages
                      if tree.branches[l.population].parent_branch_id is not None]
        if hasattr(background, "next_boundary"):
            boundaries.extend(background.next_boundary(l.population, t) for l in lineages)
        next_boundary = min(boundaries + [float("inf")])
        if total <= 0:
            if next_boundary == float("inf"):
                raise RuntimeError("ARG has no possible event and no finite species boundary")
            t = next_boundary
            continue
        wait = float(rng.exponential(1.0 / total))
        if t + wait >= next_boundary:
            t = next_boundary
            continue
        t += wait
        draw = rng.random() * total
        chosen = options[-1]
        for option in options:
            draw -= option[2]
            if draw <= 0:
                chosen = option
                break
        kind, indexes, rate = chosen
        if kind == "recombination":
            i = indexes[0]
            old = lineages[i]
            if hasattr(background, "sample_recombination"):
                sampled = background.sample_recombination(old.population, old.background, old.active,
                                                           old.segments, t, config.recombination_map, rng)
                x = sampled["position"]
                parent_background = sampled.get("parent_background")
                if sampled.get("kind") == "background_transfer":
                    old_background = old.background
                    old.background = background.transfer_state(old.background, parent_background)
                    if config.record_events:
                        events.append(ARGEvent(event_index, t, "background_transfer", (old.lineage_id,),
                                               (old.population,), (old.population,), x, rate, old.material_length,
                                               "marker-mediated transfer", str(old_background), str(old.background),
                                               str(parent_background), False, old_background != old.background))
                        event_index += 1
                    continue
            else:
                old = lineages.pop(i)
                x = config.recombination_map.sample_position(old.active, rng)
                parent_background = None
            if hasattr(background, "sample_recombination"):
                lineages.pop(i)
            left_node, right_node = add_node("recombination_parent", t), add_node("recombination_parent", t)
            old_left, old_right = split(old.segments, x)
            left, right = recombine(old, x, left_node, right_node, t, config.recombination_map,
                                    background, next_lineage, next_lineage + 1, parent_background)
            for s in old_left:
                for active_left, active_right in old.active:
                    _add_edge(edges, s.node, left_node,
                              max(s.left, active_left), min(s.right, active_right))
            for s in old_right:
                for active_left, active_right in old.active:
                    _add_edge(edges, s.node, right_node,
                              max(s.left, active_left), min(s.right, active_right))
            lineages.extend((left, right))
            next_lineage += 2
            if config.record_events:
                events.append(ARGEvent(event_index, t, "recombination", (old.lineage_id, left.lineage_id, right.lineage_id),
                                       (old.population,), (old.population,), x, rate, old.material_length,
                                       "physical crossover", str(old.background),
                                       f"{left.background},{right.background}", str(parent_background),
                                       True, parent_background is not None and parent_background != old.background))
                event_index += 1
        else:
            ia, ib = indexes
            a, b = lineages[ia], lineages[ib]
            overlap = []
            for sa in a.segments:
                for sb in b.segments:
                    left, right = max(sa.left, sb.left), min(sa.right, sb.right)
                    if right > left:
                        overlap.append((left, right))
            overlap = interval_union(overlap)
            parent_node = add_node("coalescence", t, children=[a.lineage_id, b.lineage_id]) if overlap else None
            if parent_node is not None:
                for child in (a, b):
                    for s in child.segments:
                        for left, right in overlap:
                            _add_edge(edges, s.node, parent_node, max(s.left, left), min(s.right, right))
            parent = union_lineage(a, b, parent_node, t, next_lineage)
            next_lineage += 1
            for index in sorted((ia, ib), reverse=True):
                lineages.pop(index)
            lineages.append(parent)
            if config.record_events:
                events.append(ARGEvent(event_index, t, "coalescence", (a.lineage_id, b.lineage_id, parent.lineage_id),
                                       (a.population, b.population), (parent.population,), None, rate,
                                       a.material_length + b.material_length))
                event_index += 1
            refresh_completed_material()
        peak = max(peak, len(lineages))
        peak_segments = max(peak_segments, sum(len(l.segments) for l in lineages))
    return ARGResult(config.sequence_length_bp, nodes, edges, events, time.perf_counter() - start, peak, peak_segments)


def write_event_log(events: list[ARGEvent], path: Path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "time", "event_type", "lineage_ids", "population_before", "population_after", "position", "rate", "material_length", "details",
                         "background_before", "background_after", "parental_backgrounds", "ancestry_changed", "background_changed"])
        for e in events:
            writer.writerow([e.index, e.time, e.event_type, ";".join(map(str, e.lineage_ids)), ";".join(e.population_before), ";".join(e.population_after), e.position, e.rate, e.material_length, e.details,
                             e.background_before, e.background_after, e.parental_backgrounds, e.ancestry_changed, e.background_changed])


def _write_frequency_history(history, path: Path):
    records = getattr(history, "records", [])
    if not records:
        return
    fields = list(records[0].__dataclass_fields__)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({field: getattr(record, field) for field in fields})


def run_arg_config(config: dict) -> Path:
    """Run the experimental CLI mode and write the documented artifacts."""
    from ..experiments import _tree_from_config
    from .output.genealogy_export import export_genealogies
    from .output.tree_sequence import build_tree_sequence
    import yaml

    tree = _tree_from_config(config)
    arg_cfg = config.get("arg", {}) or {}
    length = float(arg_cfg.get("sequence_length_bp", 0))
    if length <= 0:
        raise ValueError("arg.sequence_length_bp must be positive")
    rec = arg_cfg.get("recombination_rate_per_bp", 0.0)
    map_cfg = arg_cfg.get("recombination_map", {}) or {}
    if map_cfg:
        points = tuple(float(x) for x in map_cfg.get("breakpoints", ()))
        rates = tuple(float(x) for x in map_cfg.get("rates", (rec,)))
        rec_map = PiecewiseConstantRecombinationMap(length, points, rates)
    else:
        rec_map = PiecewiseConstantRecombinationMap.constant(length, float(rec))
    rng = np.random.default_rng(int(config.get("seed", 1)))
    model_name = str(arg_cfg.get("background_model", "neutral"))
    history = None
    sampled_backgrounds = arg_cfg.get("sampled_backgrounds")
    if sampled_backgrounds is not None:
        sampled_backgrounds = {str(key): int(value) for key, value in sampled_backgrounds.items()}
    if model_name == "msrc_two_state":
        from ..population_process import simulate_population_history
        from ..rearrangement import Rearrangement
        rr = config.get("rearrangement", {}) or {}
        selection = (rr.get("selection", {}) or {}).get("coefficient", rr.get("selection_coefficient", 0.0))
        rearrangement = Rearrangement(rr.get("id", "inv_1"), rr.get("type", "inversion"), rr["origin_branch"],
                                     int(rr["origin_time_from_branch_start"]), int(rr.get("initial_copy_count", 1)), float(selection))
        history_path = arg_cfg.get("frequency_history_path")
        if history_path:
            from .history import load_frequency_history
            history = load_frequency_history(history_path, arg_cfg.get("structured_origin_age"))
        else:
            history = simulate_population_history(tree, rearrangement, rng, config)
        if sampled_backgrounds is None:
            sampled_backgrounds = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        interval = arg_cfg.get("rearrangement_interval", {}) or {}
        model = MSRCTwoBackground(history, float(arg_cfg["background_anchor_bp"]),
                                  (float(interval["start"]), float(interval["end"])),
                                  float(arg_cfg.get("cross_background_suppression", 1.0)), sampled_backgrounds)
    else:
        model = "neutral"
    result = simulate_arg(tree, ARGConfig(length, rec_map, model,
                                          bool((config.get("output", {}) or {}).get("record_events", True)),
                                          history, sampled_backgrounds), rng)
    out = Path((config.get("output", {}) or {}).get("directory", "arg_output"))
    out.mkdir(parents=True, exist_ok=True)
    with (out / "config.resolved.yaml").open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    ts = build_tree_sequence(result)
    ts.dump(out / "simulated_arg.trees")
    export_genealogies(ts, out / "local_gene_trees.nwk", out / "local_gene_trees.csv",
                       chrom=str(arg_cfg.get("chrom", "chr1")))
    write_event_log(result.events, out / "arg_events.csv")
    if history is not None:
        _write_frequency_history(history, out / "frequency_history.csv")
        (out / "sample_backgrounds.csv").write_text("taxon,background\n" + "\n".join(
            f"{taxon},{sampled_backgrounds[taxon]}" for taxon in tree.taxa) + "\n")
        diagnostic_rows = [e for e in result.events if e.event_type in {"recombination", "background_transfer"}]
        with (out / "msrc_arg_diagnostics.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time", "population", "position", "event_type", "rate", "background_before", "background_after", "parental_backgrounds", "ancestry_changed", "background_changed"])
            for e in diagnostic_rows:
                writer.writerow([e.time, e.population_before[0], e.position, e.event_type, e.rate, e.background_before, e.background_after,
                                 e.parental_backgrounds, e.ancestry_changed, e.background_changed])
    summary = {
        "mode": "arg", "model": "hudson", "background_model": model_name,
        "sequence_length_bp": length, "num_nodes": len(result.nodes), "num_edges": len(result.edges),
        "num_events": len(result.events), "num_marginal_trees": ts.num_trees,
        "runtime_seconds": result.runtime_seconds, "peak_active_lineages": result.peak_lineages,
        "peak_segments": result.peak_segments,
        "recombination_events": sum(e.event_type == "recombination" for e in result.events),
        "coalescence_events": sum(e.event_type == "coalescence" for e in result.events),
        "background_transfer_events": sum(e.event_type == "background_transfer" for e in result.events),
    }
    (out / "arg_summary.json").write_text(json.dumps(summary, indent=2))
    return out
