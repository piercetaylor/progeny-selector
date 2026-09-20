"""Chromosome lengths per assembly (docs/adr/0015)."""

from __future__ import annotations

from progeny_selector.constants import SOYBEAN_CHROM_LENGTHS_BP, SOYBEAN_CHROM_LENGTHS_BP_WM82A4, SOYBEAN_CHROMOSOMES
from progeny_selector.core.chrom import chrom_length_bp, normalize_chrom
from progeny_selector.model.criteria import ASSEMBLIES, DEFAULT_ASSEMBLY


def test_named_assembly_table():
    # One value per assembly that no other table shares, so a copy-pasted table fails here.
    assert chrom_length_bp("Gm01", None, "Wm82.a2") == 56_831_625
    assert chrom_length_bp("Gm01", None, "Wm82.a1") == 55_915_596
    assert chrom_length_bp("Gm01", None, "Wm82.a4") == 57_932_356
    assert chrom_length_bp("Gm11", None, "Wm82.a1") == 39_172_791
    assert chrom_length_bp("Gm11", None, "Wm82.a2") == 34_766_868
    assert chrom_length_bp("Gm11", None, "Wm82.a4") == 39_643_746
    assert chrom_length_bp("Gm18", None, "Wm82.a1") == 62_308_141  # a1's longest chromosome, 1-based end as published
    assert chrom_length_bp("Gm18", None, "Wm82.a2") == 58_018_743
    assert chrom_length_bp("Gm18", None, "Wm82.a4") == 58_286_271


def test_the_tables_are_distinct():
    values = [tuple(t.values()) for t in SOYBEAN_CHROM_LENGTHS_BP.values()]
    assert len(set(values)) == len(values)
    assert "Wm82.a5" not in SOYBEAN_CHROM_LENGTHS_BP  # left out until the lengths are verified (docs/adr/0015)
    assert "Wm82.a6" not in SOYBEAN_CHROM_LENGTHS_BP


def test_none_and_unknown_assembly_return_the_fallback():
    assert chrom_length_bp("Gm01", 123, "none") == 123
    assert chrom_length_bp("Gm01", None, "none") is None
    assert chrom_length_bp("Gm01", 123, "Wm82.a9") == 123
    assert chrom_length_bp("Gm01", 123, "Wm82.a5") == 123


def test_default_is_wm82a4_and_unchanged():
    assert DEFAULT_ASSEMBLY == "Wm82.a4"
    for chrom, length in SOYBEAN_CHROM_LENGTHS_BP_WM82A4.items():
        assert chrom_length_bp(chrom) == length
        assert chrom_length_bp(chrom, None, DEFAULT_ASSEMBLY) == length
    assert chrom_length_bp("scaffold_1", None) is None
    assert chrom_length_bp("chr6", None) == SOYBEAN_CHROM_LENGTHS_BP_WM82A4["Gm06"]  # names are normalised first


def test_every_assembly_covers_the_twenty_chromosomes():
    assert set(SOYBEAN_CHROM_LENGTHS_BP) == set(ASSEMBLIES) - {"none"}
    for table in SOYBEAN_CHROM_LENGTHS_BP.values():
        assert tuple(table) == SOYBEAN_CHROMOSOMES
        assert all(v > 0 for v in table.values())


def test_normalize_chrom_untouched():
    assert normalize_chrom("Chromosome_06") == "Gm06"
    assert normalize_chrom("scaffold_1") == "scaffold_1"
