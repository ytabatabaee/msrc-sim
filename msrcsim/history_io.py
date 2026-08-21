from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping
import yaml
from .wright_fisher import FrequencyRecord, FrequencyHistory


def save_frozen_history(path: str | Path, config: Mapping[str, Any], history: FrequencyHistory,
                        sampled_arrangements: Mapping[str, int], metadata: Mapping[str, Any] | None = None) -> Path:
    path = Path(path)
    payload = {
        "format_version": 1,
        "config": dict(config),
        "sampled_arrangements": {str(k): int(v) for k, v in sampled_arrangements.items()},
        "metadata": dict(metadata or {}),
        "frequency_records": [asdict(r) for r in history.records],
    }
    with path.open("w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    return path


def load_frozen_history(path: str | Path) -> tuple[dict[str, Any], FrequencyHistory, dict[str, int], dict[str, Any]]:
    with Path(path).open() as handle:
        payload = yaml.safe_load(handle)
    records = [FrequencyRecord(**row) for row in payload["frequency_records"]]
    by: dict[str, list[FrequencyRecord]] = {}
    for rec in records:
        by.setdefault(rec.branch_id, []).append(rec)
    return payload["config"], FrequencyHistory(records, by), {
        str(k): int(v) for k, v in payload["sampled_arrangements"].items()
    }, dict(payload.get("metadata", {}))
