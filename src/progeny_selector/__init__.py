"""progeny-selector: marker-assisted backcross progeny ranking and selection.

Public surface: ``load_dataset``, ``read_criteria``, ``run_analysis`` and the
export writers. Everything under ``core`` is pure and UI-free.
"""

from progeny_selector._version import __version__
from progeny_selector.core.pipeline import AnalysisResult, run_analysis
from progeny_selector.io import load_dataset, read_criteria

__all__ = ["AnalysisResult", "__version__", "load_dataset", "read_criteria", "run_analysis"]
