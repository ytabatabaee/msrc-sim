from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .spatial_plotting import plot_spatial_output
from .spatial import LOCUS_FIELDS, WINDOW_FIELDS
from .spatial_statistics import RearrangementInterval, classify_region, sliding_windows, spatial_summary


TOPOLOGY_TO_INDEX = {"12|34": 0, "13|24": 1, "14|23": 2, "0": 0, "1": 1, "2": 2}


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize an ordered locus topology table")
    parser.add_argument("--input", required=True)
    parser.add_argument("--breakpoints", required=True, help="Comma-separated start,end rearrangement breakpoints")
    parser.add_argument("--window-loci", type=int, default=50)
    parser.add_argument("--step-loci", type=int, default=10)
    parser.add_argument("--output", required=True)
    parser.add_argument("--format", choices=["png", "pdf", "svg"], default="png")
    args = parser.parse_args()
    bps = [int(x) for x in args.breakpoints.split(",")]
    if len(bps) != 2 or bps[0] >= bps[1]:
        raise ValueError("--breakpoints must be start,end with start < end")
    interval = RearrangementInterval(bps[0], bps[1])
    with open(args.input, newline="") as handle:
        raw = sorted(csv.DictReader(handle), key=lambda r: int(r["position"]))
    rows = []
    for locus_id, row in enumerate(raw):
        topo = str(row["topology"])
        idx = TOPOLOGY_TO_INDEX.get(topo)
        if idx is None:
            raise ValueError(f"Unsupported topology value: {topo}")
        region, inside = classify_region(int(row["position"]), interval)
        rows.append({
            "locus_id": int(row.get("locus_id", locus_id)),
            "position": int(row["position"]),
            "region": region,
            "inside_rearrangement": inside,
            "local_model": row.get("local_model", "external"),
            "local_effective_cross_arrangement_fraction": row.get("local_effective_cross_arrangement_fraction", ""),
            "topology_index": idx,
            "topology": ["12|34", "13|24", "14|23"][idx],
            "terminal_pattern": row.get("terminal_pattern", ""),
            "num_switch_events": row.get("num_switch_events", ""),
            "mean_or_first_coalescence_time": row.get("mean_or_first_coalescence_time", ""),
        })
    windows = sliding_windows(rows, args.window_loci, args.step_loci)
    summary = spatial_summary(rows, windows, interval, args.window_loci)
    summary.update({"mode": "spatial_summarize", "rearrangement_interval": {"start": interval.start, "end": interval.end, "interval_id": interval.interval_id}})
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "spatial_loci.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOCUS_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    with (out / "spatial_windows.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WINDOW_FIELDS)
        writer.writeheader()
        writer.writerows({k: row[k] for k in WINDOW_FIELDS} for row in windows)
    with (out / "spatial_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    plot_spatial_output(out, fmt=args.format)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

