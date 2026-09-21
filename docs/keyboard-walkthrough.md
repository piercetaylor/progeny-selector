# Keyboard walkthrough

From page load to a downloaded `selected.csv`, using only the keyboard, with the
places a mouse is still needed called out. Observed against `shiny` 1.7 in
Chromium (`tests/e2e/test_keyboard_walkthrough.py`); the DataGrid's own keyboard
handling is a third-party component this project does not control.

## Load

Tab order on a fresh page (contract 1.4.0 added the token-profile controls
after this screen was first specced, so the order below is what the browser
does today, not the four-file list an earlier draft named):

1. Genotypes file input
2. samples.csv file input
3. markers.csv file input
4. Token profile select
5. Custom token profile file input
6. "Clear custom token profile" button
7. criteria.yaml file input
8. "Load and analyse" button
9. The criteria text area, "Apply criteria and re-analyse", "Download criteria.yaml"

File inputs are opened and populated with `Enter`/`Space` on the button they
render as, but choosing a file from the native file-picker dialog itself needs
the OS file dialog, which is not keyboard-scriptable from here; the e2e tests
upload files directly through Playwright's file-chooser API instead of a real
dialog. `Enter` on "Load and analyse" runs the analysis.

## Navbar

`Tab` from the page content reaches the tab list; `ArrowRight`/`ArrowLeft` move
the focus between tabs without switching screens, and `Enter` (or `Space`)
activates the focused tab.

## Navigate

The family accordion headers are buttons: `Enter`/`Space` opens or closes a
family's panel. Inside an open panel, the generation choice is a radio group:
arrow keys move the selection between generations, which is a change event
(no extra `Enter` needed).

## Rank

**Needs a mouse.** Tabbing through the page does not give the grid's cells
keyboard focus at all; the `Tab` sequence around the grid skips over it and
returns to the page chrome. A cell must be clicked first. Once a cell has been
clicked (and its row thereby selected), arrow keys move the focused cell
without changing the selection.

What `Space` does to the selection once inside the grid is not documented
here and the e2e test does not assert it: on one run it toggled the focused
row into a growing multi-row selection, and on another (CI) run, with the
same steps, it cleared the selection instead. That is Shiny DataGrid's own
keyboard handling, not something this project implements, and it is not
relied on anywhere in the app: multi-row selection for Compare and Selection
list is exercised through the mouse (click plus Ctrl/Cmd-click) in the other
e2e tests, and there is no documented keyboard path to it.

## Selection list

Tab into the grid, then `Enter` opens the focused cell for editing (only the
`notes` column is editable); typing replaces the text; `Enter` commits the
edit; `Escape` cancels and restores the previous value. Committing a note
re-renders the whole grid on the server (the render function reads
`state.notes()`), so the next keyboard action should wait for that
round-trip rather than being fired immediately after the previous commit.

## Export

Tab to a download button; `Enter` triggers the browser's download exactly as
a click would.
