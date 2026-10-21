from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping
import yaml
from .wright_fisher import FrequencyRecord, FrequencyHistory
from .population_process import resolved_population_process, resolved_speciation_split
from .speciation_split import split_events_as_dicts


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
    split_events = split_events_as_dicts(getattr(history, "split_events", []))
    if split_events:
        payload["split_events"] = split_events
    payload["metadata"].setdefault("population_process", resolved_population_process(dict(config)))
    split_provenance = resolved_speciation_split(dict(config))
    if split_provenance is not None:
        payload["metadata"].setdefault("speciation_split", split_provenance)
    if getattr(history, "structured_origin_age", None) is not None:
        payload["metadata"].setdefault("structured_origin_age", float(history.structured_origin_age))
        payload["metadata"].setdefault("structured_origin_branch", getattr(history, "structured_origin_branch", None))
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
    metadata = dict(payload.get("metadata", {}))
    return payload["config"], FrequencyHistory(
        records,
        by,
        list(payload.get("split_events", [])),
        metadata.get("structured_origin_age"),
        metadata.get("structured_origin_branch"),
    ), {
        str(k): int(v) for k, v in payload["sampled_arrangements"].items()
    }, metadata
