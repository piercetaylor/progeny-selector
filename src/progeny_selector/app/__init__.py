"""Shiny for Python user interface. Every screen is a module with ``ui()`` and ``server()``.

The UI never computes metrics itself; it calls ``progeny_selector.core.run_analysis``
and renders the returned rows. Run locally with ``shiny run src/progeny_selector/app/app.py``
or export as a static site with ``shinylive export src/progeny_selector/app site``.
"""
