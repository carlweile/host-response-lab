"""Pin the GEO directory-naming rule.

The first version of this code special-cased three-digit accessions and turned
GPL571 into GPL571nnn, which 404s. The correct forms below were checked against
the live server on 2026-10-10. These are offline string assertions; no network.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.fetch_geo import platform_annot_url, series_matrix_url  # noqa: E402


@pytest.mark.parametrize("accession,expected_dir", [
    ("GSE63990", "GSE63nnn"),
    ("GSE42026", "GSE42nnn"),
    ("GSE6269", "GSE6nnn"),
    ("GSE60244", "GSE60nnn"),
])
def test_series_directory_grouping(accession, expected_dir):
    url = series_matrix_url(accession)
    assert f"/series/{expected_dir}/{accession}/matrix/" in url
    assert url.endswith(f"{accession}_series_matrix.txt.gz")


@pytest.mark.parametrize("platform,expected_dir", [
    ("GPL571", "GPLnnn"),      # three digits collapse to the bare form
    ("GPL96", "GPLnnn"),
    ("GPL6947", "GPL6nnn"),
    ("GPL10558", "GPL10nnn"),
])
def test_platform_directory_grouping(platform, expected_dir):
    url = platform_annot_url(platform)
    assert f"/platforms/{expected_dir}/{platform}/annot/" in url


def test_three_digit_platform_does_not_repeat_its_own_digits():
    """The exact shape of the bug that was fixed."""
    assert "GPL571nnn" not in platform_annot_url("GPL571")
