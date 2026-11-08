from __future__ import annotations

import csv


def export_genealogies(ts, newick_path, csv_path):
    with open(newick_path, "w") as nw, open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["left", "right", "newick"])
        for tree in ts.trees():
            newick = tree.as_newick(precision=12)
            nw.write(f"[{tree.interval.left},{tree.interval.right})\t{newick}\n")
            writer.writerow([tree.interval.left, tree.interval.right, newick])
