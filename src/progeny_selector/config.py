"""Twelve-factor configuration read from environment variables (see .env.example).

Responsibility: the only place that touches ``os.environ``. Values are typed,
defaulted, and documented here; nothing else in the package reads the
environment. No secrets are expected: the tool processes local files.

Interface:
    Settings (frozen dataclass); load_settings() -> Settings
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class Settings:
    max_upload_mb: int = 500  # refuse genotype files above this size in the UI
    max_samples_duplicate_check: int = 2000  # pairwise IBS duplicate check is skipped above this
    default_flank_window_cm: float = 5.0
    default_max_marker_coverage_cm: float = 10.0
    log_level: str = "INFO"
    data_dir: str | None = None  # optional directory the UI lists for quick loading


def _env(name: str, default: str) -> str:
    return os.environ.get(f"PROGENY_SELECTOR_{name}", default)


def load_settings() -> Settings:
    return Settings(
        max_upload_mb=int(_env("MAX_UPLOAD_MB", "500")),
        max_samples_duplicate_check=int(_env("MAX_SAMPLES_DUPLICATE_CHECK", "2000")),
        default_flank_window_cm=float(_env("DEFAULT_FLANK_WINDOW_CM", "5")),
        default_max_marker_coverage_cm=float(_env("DEFAULT_MAX_MARKER_COVERAGE_CM", "10")),
        log_level=_env("LOG_LEVEL", "INFO"),
        data_dir=os.environ.get("PROGENY_SELECTOR_DATA_DIR"),
    )
