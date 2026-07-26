from __future__ import annotations
import argparse, csv
from pathlib import Path
from .model_fitting import compare_models


def _counts(row: dict[str, str]) -> tuple[int, int, int]:
    if all(k in row for k in ("n1", "n2", "n3")):
        return int(row["n1"]), int(row["n2"]), int(row["n3"])
    loci = int(float(row.get("num_loci", row.get("loci_per_replicate", 0))))
    if loci <= 0:
        raise ValueError("input needs n1/n2/n3 or num_loci plus q1/q2/q3")
    raw = [float(row[f"q{i}"]) * loci for i in (1, 2, 3)]
    counts = [int(round(x)) for x in raw]
    counts[0] += loci - sum(counts)
    return tuple(counts)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit MSC and two-tree quartet-mixture models")
    parser.add_argument("--input", required=True, help="replicate_summary.csv or equivalent")
    parser.add_argument("--output", default="model_comparison.csv")
    args = parser.parse_args()
    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle))
    out_rows = []
    for row in rows:
        if str(row.get("accepted", "True")).lower() in {"false", "0"}:
            continue
        result = dict(row)
        result.update(compare_models(_counts(row)))
        out_rows.append(result)
    if not out_rows:
        raise ValueError("no accepted rows found")
    path = Path(args.output)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0]))
        writer.writeheader(); writer.writerows(out_rows)
    print(f"Wrote {path}")

if __name__ == "__main__":
    main()
