"""Download GEO series matrices and platform annotations.

Both live at deterministic URLs on the NCBI FTP mirror, so no API key or
scraping is needed. Files are cached on disk; re-running is free.
"""
from __future__ import annotations

import re
from pathlib import Path

import requests

FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
TIMEOUT = 120
CHUNK = 1 << 16


def _stub(accession: str) -> str:
    """GEO groups directories by replacing the last three digits with 'nnn'.

    GSE63990 -> GSE63nnn, GPL10558 -> GPL10nnn, GPL571 -> GPLnnn.
    Verified against the live server on 2026-10-10; a three-digit accession
    really does collapse to the bare 'nnn' form, which an earlier version of
    this function got wrong and silently turned into 404s.
    """
    prefix = re.match(r"^[A-Za-z]+", accession).group(0)
    digits = accession[len(prefix):]
    return f"{prefix}{digits[:-3]}nnn"


def series_matrix_url(accession: str) -> str:
    return f"{FTP}/series/{_stub(accession)}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def platform_annot_url(platform: str) -> str:
    return f"{FTP}/platforms/{_stub(platform)}/{platform}/annot/{platform}.annot.gz"


def download(url: str, dest: Path, force: bool = False) -> tuple[bool, str]:
    """Fetch url to dest. Returns (ok, message). Never raises on HTTP errors."""
    if dest.exists() and dest.stat().st_size > 0 and not force:
        mb = dest.stat().st_size / 1e6
        return True, f"cached ({mb:.1f} MB)"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with requests.get(url, stream=True, timeout=TIMEOUT) as r:
            if r.status_code != 200:
                return False, f"HTTP {r.status_code}"
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(CHUNK):
                    fh.write(chunk)
        tmp.replace(dest)
    except requests.RequestException as exc:
        tmp.unlink(missing_ok=True)
        return False, f"{type(exc).__name__}: {exc}"
    mb = dest.stat().st_size / 1e6
    return True, f"downloaded ({mb:.1f} MB)"


def matrix_dir_url(accession: str) -> str:
    return f"{FTP}/series/{_stub(accession)}/{accession}/matrix/"


def list_matrix_files(accession: str) -> list[str]:
    """Series run on several platforms have no single matrix file, only one per
    platform (GSE6269-GPL96_series_matrix.txt.gz and friends). List what exists."""
    try:
        r = requests.get(matrix_dir_url(accession), timeout=TIMEOUT)
        if r.status_code != 200:
            return []
        names = re.findall(r'href="([^"]*_series_matrix\.txt\.gz)"', r.text)
        return sorted(set(names))
    except requests.RequestException:
        return []


def fetch_series(
    accession: str, data_dir: Path, force: bool = False, prefer_platform: str | None = None
) -> tuple[Path | None, str]:
    """Download one series matrix, always caching it under the plain accession name.

    When a series exists only as per-platform files we take exactly ONE of them.
    The variants are the same patients from the same study re-run on different
    arrays, so treating them as separate cohorts would split one study across
    train and test and defeat the entire point of this repo.
    """
    dest = data_dir / "series" / f"{accession}_series_matrix.txt.gz"
    msg = f"{prefer_platform} variant required"
    if not prefer_platform:
        ok, msg = download(series_matrix_url(accession), dest, force)
        if ok:
            return dest, msg

    variants = list_matrix_files(accession)
    if not variants:
        return None, msg
    chosen = variants[0]
    if prefer_platform:
        # exact id: "GPL96" must not match "GPL961"; never fall back to another platform
        match = [v for v in variants if f"-{prefer_platform}_" in v]
        if not match:
            return None, f"no {prefer_platform} variant among {variants}"
        chosen = match[0]
    url = f"{matrix_dir_url(accession)}{chosen}"
    ok2, msg2 = download(url, dest, force)
    extra = (f"{msg2}; used {chosen} of {len(variants)} platform variants, "
             f"others excluded as the same study")
    return (dest if ok2 else None), extra


def fetch_platform(platform: str, data_dir: Path, force: bool = False) -> tuple[Path | None, str]:
    dest = data_dir / "platforms" / f"{platform}.annot.gz"
    ok, msg = download(platform_annot_url(platform), dest, force)
    return (dest if ok else None), msg
