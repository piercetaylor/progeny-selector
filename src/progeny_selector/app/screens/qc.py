"""Screen 2, Validate and QC: contract warnings, parent checks, per-individual QC.

Responsibility: show ``result.warnings``, informative-marker counts, and the
``SampleQC`` table (missing, het, hom-donor and non-parental rates versus the
generation expectation, flags). Placeholder in M0; interface fixed.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result)
"""

from __future__ import annotations

from shiny import module, render, ui


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(ui.card_header("Validation and QC"), ui.output_text_verbatim("summary"))


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @render.text
    def summary() -> str:
        result = state.result()
        if result is None:
            return "Load data first."
        lines = [f"informative markers: {result.classification.n_informative}"]
        lines += [f"warning: {w}" for w in result.warnings]
        for q in result.qc:
            if q.flags:
                flags = "|".join(q.flags)
                lines.append(
                    f"{q.sample_id}: het {q.het_rate:.3f} (expected {q.expected_het}), missing {q.missing_rate:.3f}, flags {flags}"
                )
        return "\n".join(lines)


def server(id: str, state) -> None:
    server_(id, state)
