from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

@dataclass(frozen=True)
class RunConfig:
    raw: dict
    path: Path


def load_config(path: str|Path) -> RunConfig:
    p=Path(path)
    with p.open() as h: raw=yaml.safe_load(h)
    if raw.get("mode") not in {"conditional","mechanistic"}:
        raise ValueError("mode must be conditional or mechanistic")
    return RunConfig(raw,p)
