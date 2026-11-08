from __future__ import annotations

from dataclasses import dataclass, field
import csv
import json
from pathlib import Path
import time

import numpy as np

from ..species_tree import SpeciesTree
from .backgrounds.neutral import NeutralBackground
from .coalescence import union_lineage
from .events import ARGEvent
from .lineages import ARGLineage
from .recombination import recombine
from .recombination_map import PiecewiseConstantRecombinationMap
from .segments import Segment, split, total_length


@dataclass(frozen=True)
class ARGConfig:
    sequence_length_bp: float
    recombination_map: PiecewiseConstantRecombinationMap
    background_model: str = "neutral"
    record_events: bool = True


@dataclass
class ARGResult:
    sequence_length: float
    nodes: list[dict]
    edges: list[tuple[float, float, int, int]]
    events: list[ARGEvent]
    runtime_seconds: float = 0.0
    peak_lineages: int = 0

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
    if config.background_model != "neutral":
        raise ValueError("only background_model: neutral is implemented")
    if abs(config.recombination_map.length - config.sequence_length_bp) > 1e-9:
        raise ValueError("recombination map length must equal sequence_length_bp")
    start = time.perf_counter()
    background = NeutralBackground()
    nodes: list[dict] = []
    edges: list[tuple[float, float, int, int]] = []
    lineages: list[ARGLineage] = []
    next_node = 0
    next_lineage = 0
    for taxon in tree.taxa:
        node = next_node
        next_node += 1
        nodes.append({"flags": 1, "time": 0.0, "metadata": {"type": "sample", "taxon": taxon}})
        state = background.initial_state(tree.branch_for_tip(taxon), (0.0, config.sequence_length_bp))
        lineages.append(ARGLineage(next_lineage, tree.branch_for_tip(taxon),
                                   (Segment(0.0, config.sequence_length_bp, node),), 0.0, state))
        next_lineage += 1
    events: list[ARGEvent] = []
    event_index = 0
    peak = len(lineages)
    t = 0.0

    def add_node(event_type, t, **metadata):
        nonlocal next_node
        node_id = next_node
        next_node += 1
        nodes.append({"flags": 0, "time": float(t), "metadata": {"type": event_type, **metadata}})
        return node_id

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
            hazard = config.recombination_map.integrated_rate([(s.left, s.right) for s in a.segments])
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
            old = lineages.pop(i)
            x = config.recombination_map.sample_position([(s.left, s.right) for s in old.segments], rng)
            left_id, right_id = add_node("recombination_parent", t), add_node("recombination_parent", t)
            old_left, old_right = split(old.segments, x)
            left, right = recombine(old, x, next_lineage, next_lineage + 1, t, config.recombination_map, background)
            for s in old_left:
                _add_edge(edges, s.node, left_id, s.left, s.right)
            for s in old_right:
                _add_edge(edges, s.node, right_id, s.left, s.right)
            lineages.extend((left, right))
            next_lineage += 2
            if config.record_events:
                events.append(ARGEvent(event_index, t, "recombination", (old.lineage_id, left.lineage_id, right.lineage_id),
                                       (old.population,), (old.population,), x, rate, old.material_length))
                event_index += 1
        else:
            ia, ib = indexes
            a, b = lineages[ia], lineages[ib]
            parent_node = add_node("coalescence", t, children=[a.lineage_id, b.lineage_id])
            for child in (a, b):
                for s in child.segments:
                    _add_edge(edges, s.node, parent_node, s.left, s.right)
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
        peak = max(peak, len(lineages))
    return ARGResult(config.sequence_length_bp, nodes, edges, events, time.perf_counter() - start, peak)


def write_event_log(events: list[ARGEvent], path: Path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "time", "event_type", "lineage_ids", "population_before", "population_after", "position", "rate", "material_length", "details"])
        for e in events:
            writer.writerow([e.index, e.time, e.event_type, ";".join(map(str, e.lineage_ids)), ";".join(e.population_before), ";".join(e.population_after), e.position, e.rate, e.material_length, e.details])


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
    result = simulate_arg(tree, ARGConfig(length, rec_map, str(arg_cfg.get("background_model", "neutral")),
                                          bool((config.get("output", {}) or {}).get("record_events", True))),
                          np.random.default_rng(int(config.get("seed", 1))))
    out = Path((config.get("output", {}) or {}).get("directory", "arg_output"))
    out.mkdir(parents=True, exist_ok=True)
    with (out / "config.resolved.yaml").open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    ts = build_tree_sequence(result)
    ts.dump(out / "simulated_arg.trees")
    export_genealogies(ts, out / "local_gene_trees.nwk", out / "local_gene_trees.csv")
    write_event_log(result.events, out / "arg_events.csv")
    summary = {
        "mode": "arg", "model": "hudson", "background_model": "neutral",
        "sequence_length_bp": length, "num_nodes": len(result.nodes), "num_edges": len(result.edges),
        "num_events": len(result.events), "num_marginal_trees": ts.num_trees,
        "runtime_seconds": result.runtime_seconds, "peak_active_lineages": result.peak_lineages,
        "recombination_events": sum(e.event_type == "recombination" for e in result.events),
        "coalescence_events": sum(e.event_type == "coalescence" for e in result.events),
    }
    (out / "arg_summary.json").write_text(json.dumps(summary, indent=2))
    return out
