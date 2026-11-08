from __future__ import annotations

import json


def build_tree_sequence(result):
    try:
        import tskit
    except ImportError as exc:
        raise ImportError("ARG output requires tskit; install msrc-sim with its runtime dependencies") from exc
    tables = tskit.TableCollection(result.sequence_length)
    for node in result.nodes:
        tables.nodes.add_row(flags=node["flags"], time=node["time"], metadata=json.dumps(node["metadata"]).encode())
    for edge in result.edges:
        tables.edges.add_row(edge[0], edge[1], edge[2], edge[3])
    tables.sort()
    # Loading a TreeSequence performs tskit's structural validation. tskit 1.0
    # no longer exposes TreeSequence.validate(); TableCollection.sort/load is
    # the stable validation path across supported versions.
    return tables.tree_sequence()


def write_tree_sequence(result, path):
    ts = build_tree_sequence(result)
    ts.dump(path)
    return ts
