# Accessibility baseline: WCAG 2.2 AA, axe-core in the e2e suite

Status: accepted. Date: 2026-09-25.

## Context and Problem Statement

M3 Phase 5 (`docs/m3-phases.md`) asked for an accessibility pass with no prior baseline in this
repository: no target level, no automated check, and no record of which parts of the UI are the
app's own markup versus a third-party component's. Q3 of that phase fixed the target at WCAG 2.2
Level AA rather than A or AAA.

## Decision Drivers

- The app's screens are simple forms, tables and an SVG summary; AA is reachable without
  redesigning the two components this project does not author, Shiny's `render.DataGrid` and
  Bootstrap's navbar (via bslib).
- A review with no re-runnable check decays; axe-core running inside the existing Playwright e2e
  suite makes the baseline a gate, not a one-time audit.
- Colours are fixed to the Okabe-Ito palette (`docs/adr/0006`, `constants.py`) for colour-vision
  deficiency; contrast and colour-alone status still needed checking against that specific palette.

## Decision Outcome

Axe-core, via the `axe-playwright-python` bridge, runs in `tests/e2e/test_a11y.py` against every
screen and is part of the `e2e` CI job (`CLAUDE.md`, "Gates"). Rule tags `wcag2a`, `wcag2aa`,
`wcag21a`, `wcag21aa` and `wcag22aa` are asserted on; axe's non-WCAG best-practice rules are not.

**Exceptions policy.** A violation inside markup this app does not template or attribute itself
(the DataGrid's own `<table>`, its filter `<input>`s, bslib's fillable card-body scroll containers,
Bootstrap's navbar) is recorded in `tests/e2e/test_a11y.py::KNOWN_A11Y_EXCEPTIONS` by axe rule id
with a one-line reason, not patched. Every exception needs a dated line in this record and an entry
in `docs/accessibility.md`'s exceptions table; a rule id present in one and not the other is a bug
in the review, not a passing state (`test_exceptions_are_documented`).

Exceptions recorded so far, all dated 2026-09-25 (see `docs/accessibility.md` for the fuller
reasons and the affected screens):

- `aria-allowed-attr` — `aria-rowcount`/`aria-multiselectable` on `render.DataGrid`'s own table.
- `label` — `render.DataGrid`'s per-column filter inputs, rendered by the component when
  `filters=True` is passed; this app supplies no template for them.
- `scrollable-region-focusable` — bslib's fillable card-body container, applied to every
  `ui.card()` under `page_navbar(fillable=True)`; the missing `tabindex` is the layout engine's,
  not an attribute a screen module sets.

**The DataGrid keyboard limit stays a documented limit, not a workaround.** `docs/keyboard-walkthrough.md`
already records that the Rank grid needs a mouse: `Tab` does not give its cells keyboard focus, a
property of Shiny's DataGrid component. This review does not add an alternative keyboard selection
path in M3; multi-row selection for Compare and the Selection list already has no keyboard route
and is exercised by mouse in the other e2e tests. Building one would mean re-implementing selection
outside the component the rest of the screen relies on, a larger change than an accessibility
review should make unasked.

**Colours stay Okabe-Ito; only text colours are chosen for contrast, from `constants.py`.** Where
axe's `color-contrast` rule found black text failing against a background colour, the fix changes
only the text colour, and only through a new lookup table in `constants.py`
(`STATE_TEXT_COLORS`), never a literal in `src/progeny_selector/app/**`. No background colour in
`PALETTE_OKABE_ITO`, `STATE_COLORS` or `STATUS_COLORS` changed. The one failing case found: black
on the Okabe-Ito blue used for parent-of-origin state A (`#0072B2`) is 4.05:1, below AA's 4.5:1 for
normal text; `STATE_TEXT_COLORS["A"]` is white there. Black passes on every other state colour and
on all three status colours (`docs/accessibility.md` has the full contrast note).

## Considered Options

1. **axe-core in the e2e suite (chosen).** Reuses the existing Playwright fixtures and screens;
   runs on every commit that touches the `e2e` job, not just at review time.
2. **A one-time manual audit, recorded in a document with no executable check.** Rejected: nothing
   would catch a regression the next time a screen's markup changed.
3. **Lighthouse or `pa11y` in CI.** Both need a Node.js toolchain this repository does not have
   (`docs/accessibility.md`, "Not verified"); axe-core is a Python-native dependency already
   compatible with the pytest-playwright setup.

## Consequences

Good: the baseline is a passing gate (`pytest -m e2e tests/e2e/test_a11y.py`), not a document that
can drift from the code. `STATE_TEXT_COLORS` is one small table, not a redesign, and required no
change to any background colour or to the fixture.

Bad: the exceptions list is a permanent list of things this project is choosing not to fix, in
components it depends on but does not own; a future upgrade of Shiny or bslib could remove one of
these limits without this record noticing, since nothing here re-checks whether an excepted rule id
still actually fires.

## Revisit when

- Shiny's `render.DataGrid` gains keyboard cell focus or accessible filter labels, at which point
  the corresponding exception should be removed and the rule re-asserted on.
- A screen adds a colour not in `PALETTE_OKABE_ITO`, which needs its own contrast check against
  black and white before any text is placed on it.
- A Node.js toolchain becomes available in CI, at which point Lighthouse or `pa11y` could cross-check
  the axe-core results.
