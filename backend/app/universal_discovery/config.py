from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib


@dataclass(frozen=True)
class AppConfig:
    raw: dict

    @property
    def output_top_n(self) -> int:
        return int(self.raw["output"]["return_top_n"])

    @property
    def scoring(self) -> dict:
        return self.raw["scoring"]

    @property
    def base_filters(self) -> dict:
        return self.raw["base_filters"]

    @property
    def safety_filters(self) -> dict:
        return self.raw["safety_filters"]

    @property
    def pk_filters(self) -> dict:
        return self.raw["pk_filters"]

    @property
    def docking(self) -> dict:
        return self.raw["docking"]

    @property
    def modules(self) -> dict:
        return self.raw["module"]


def load_config(config_path: Path) -> AppConfig:
    with config_path.open("rb") as fh:
        raw = tomllib.load(fh)
    return AppConfig(raw=raw)
