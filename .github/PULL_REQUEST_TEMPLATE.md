## Summary

<!-- What changes and why. -->

## Checklist

- [ ] Gates run and pass: `ruff check .`, `ruff format --check .`, `mypy`, `pytest`, `python scripts/check_contract.py`
- [ ] Fixture regenerated (`python scripts/make_fixture.py`) if the generator changed, and `git diff --exit-code -- tests/fixtures` is clean
- [ ] CHANGELOG line added under `[Unreleased]`
- [ ] An ADR in `docs/adr/` for any decision
- [ ] No real genotype data
- [ ] No AI attribution trailer in commits or this description
- [ ] Contract changes start in backcross, then are mirrored here
