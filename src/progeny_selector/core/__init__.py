"""Pure, UI-free compute core. Every public function takes arrays/dataclasses and returns arrays/dataclasses."""

from progeny_selector.core.pipeline import AnalysisResult, run_analysis

__all__ = ["AnalysisResult", "run_analysis"]
