from __future__ import annotations

import csv
import json


def _labels(ts):
    labels = {}
    for node in ts.nodes():
        metadata = node.metadata
        if isinstance(metadata, bytes):
            metadata = metadata.decode()
        try:
            record = json.loads(metadata) if metadata else {}
        except (TypeError, json.JSONDecodeError):
            record = {}
        labels[node.id] = str(record.get("taxon", ""))
    return labels


def export_genealogies(ts, newick_path, csv_path, chrom="chr1"):
    # Unary recombination-provenance nodes are useful in the full ARG tables,
    # but downstream gene-tree tools expect ordinary marginal trees.
    ts = ts.simplify()
    labels = _labels(ts)
    with open(newick_path, "w") as nw, open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chrom", "left", "right", "midpoint", "tree_id", "newick"])
        for tree in ts.trees():
            newick = tree.as_newick(precision=12, node_labels=labels)
            nw.write(newick + "\n")
            writer.writerow([chrom, tree.interval.left, tree.interval.right,
                             (tree.interval.left + tree.interval.right) / 2,
                             tree.index, newick])
