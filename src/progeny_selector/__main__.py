"""``python -m progeny_selector`` entry point; delegates to the CLI."""

from progeny_selector.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
