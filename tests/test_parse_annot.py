"""Pin the SOFT-format platform annotation parser.

The bug this guards: GEO .annot.gz files open with '^Annotation' and a block of
'!' and '#' metadata. Reading the first non-comment line as the table header
picks up '^Annotation' and the parser then reports 'no gene symbol column',
skipping every cohort. The table fence is the only reliable anchor.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.parse_geo import parse_platform_annot  # noqa: E402

SOFT = """^Annotation
!Annotation_date = Aug 09 2016
!Annotation_platform = GPL571
#ID = ID from Platform data table
#Gene symbol = Entrez Gene symbol
!platform_table_begin
ID\tGene title\tGene symbol\tGene ID
1007_s_at\tmicroRNA 4640///discoidin domain receptor\tMIR4640///DDR1\t100616237///780
1053_at\treplication factor C subunit 2\tRFC2\t5982
117_at\theat shock protein family A member 6\tHSPA6\t3310
empty_at\tno symbol here\t\t
dashes_at\tplaceholder\t---\t
!platform_table_end
"""


@pytest.fixture
def annot_file(tmp_path: Path) -> Path:
    path = tmp_path / "GPLtest.annot.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(SOFT)
    return path


def test_reads_probes_past_the_soft_header(annot_file):
    mapping = parse_platform_annot(annot_file)
    assert mapping.loc["1053_at"] == "RFC2"
    assert mapping.loc["117_at"] == "HSPA6"
    assert "^Annotation" not in mapping.index


def test_multi_symbol_probe_takes_the_first(annot_file):
    mapping = parse_platform_annot(annot_file)
    assert mapping.loc["1007_s_at"] == "MIR4640"


def test_blank_and_placeholder_symbols_are_dropped(annot_file):
    mapping = parse_platform_annot(annot_file)
    assert "empty_at" not in mapping.index
    assert "dashes_at" not in mapping.index
    assert len(mapping) == 3


def test_missing_fence_raises_instead_of_guessing(tmp_path: Path):
    """Break-it case: no fence means we must fail loudly, not return junk."""
    path = tmp_path / "GPLbad.annot.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("^Annotation\n!Annotation_date = x\nID\tGene symbol\np1\tAAA\n")
    with pytest.raises(ValueError, match="platform_table_begin"):
        parse_platform_annot(path)
