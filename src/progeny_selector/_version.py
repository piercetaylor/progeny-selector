"""Single source of the package version; hatch reads it (pyproject.toml, [tool.hatch.version]).

Keep this module free of imports so the build and ``provenance`` can read it without a cycle.
"""

__version__ = "0.1.0"
