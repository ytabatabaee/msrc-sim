from __future__ import annotations

import csv
from pathlib import Path

from ..wright_fisher import FrequencyHistory, FrequencyRecord


def load_frequency_history(path: str | Path, structured_origin_age: float | None = None):
    """Load a previously written full frequency-history CSV/TSV."""
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t" if path.suffix.lower() in {".tsv", ".tab"} else ","))
    if not rows:
        raise ValueError(f"frequency history is empty: {path}")
    records = []
    required = set(FrequencyRecord.__dataclass_fields__)
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"frequency history is missing columns: {sorted(missing)}")
    bools = {"true": True, "false": False, "1": True, "0": False}
    for row in rows:
        values = dict(row)
        for field in ("forward_generation", "population_chromosomes", "effective_population_size"):
            values[field] = int(float(values[field]))
        for field in ("absolute_age", "copy_count_A1", "copy_count_A0", "frequency_A1", "frequency_A0", "selection_coefficient"):
            values[field] = float(values[field])
        for field in ("is_origin", "is_branch_start", "is_branch_end"):
            values[field] = bools[str(values[field]).strip().lower()]
        records.append(FrequencyRecord(**{field: values[field] for field in required}))
    by_branch = {}
    for record in records:
        by_branch.setdefault(record.branch_id, []).append(record)
    origin = structured_origin_age
    if origin is None:
        origin_records = [r for r in records if r.is_origin]
        origin = origin_records[0].absolute_age if origin_records else None
    origin_branch = origin_records[0].branch_id if (origin_records := [r for r in records if r.is_origin]) else None
    return FrequencyHistory(records, by_branch, [], origin, origin_branch)
