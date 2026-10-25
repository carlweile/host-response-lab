"""Parse a GEO series matrix file into an expression table plus sample metadata.

Series matrix format: '!'-prefixed metadata lines, then the expression table
fenced between '!series_matrix_table_begin' and '!series_matrix_table_end'.
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

BEGIN = "!series_matrix_table_begin"
END = "!series_matrix_table_end"


@dataclass
class Series:
    accession: str
    platform: str
    expression: pd.DataFrame                      # probes x samples
    characteristics: dict[str, list[str]]         # sample -> raw characteristic strings
    titles: dict[str, str] = field(default_factory=dict)

    @property
    def samples(self) -> list[str]:
        return list(self.expression.columns)


def _split_meta_line(line: str) -> tuple[str, list[str]]:
    parts = line.rstrip("\n").split("\t")
    key = parts[0].lstrip("!").strip()
    vals = [p.strip().strip('"') for p in parts[1:]]
    return key, vals


def parse_series_matrix(path: Path, accession: str) -> Series:
    meta: dict[str, list[list[str]]] = {}
    table_lines: list[str] = []
    in_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(BEGIN):
                in_table = True
                continue
            if line.startswith(END):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
            elif line.startswith("!"):
                key, vals = _split_meta_line(line)
                meta.setdefault(key, []).append(vals)

    if not table_lines:
        raise ValueError(f"{accession}: no expression table found in {path.name}")

    expr = pd.read_csv(
        io.StringIO("".join(table_lines)),
        sep="\t", index_col=0, low_memory=False,
    )
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = [str(c).strip('"') for c in expr.columns]
    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")

    accessions = meta.get("Sample_geo_accession", [[]])[0]
    if accessions and len(accessions) == expr.shape[1]:
        expr.columns = accessions
    samples = list(expr.columns)

    chars: dict[str, list[str]] = {s: [] for s in samples}
    for row in meta.get("Sample_characteristics_ch1", []):
        for sample, value in zip(samples, row):
            if value:
                chars[sample].append(value)
    # Sample_source_name_ch1 and the title often carry the only usable label.
    for key in ("Sample_source_name_ch1", "Sample_title"):
        for row in meta.get(key, []):
            for sample, value in zip(samples, row):
                if value:
                    chars[sample].append(value)

    titles = {}
    if meta.get("Sample_title"):
        titles = dict(zip(samples, meta["Sample_title"][0]))

    platforms = meta.get("Sample_platform_id", [[]])[0]
    platform = platforms[0] if platforms else ""

    return Series(accession, platform, expr, chars, titles)


PLATFORM_BEGIN = "!platform_table_begin"
PLATFORM_END = "!platform_table_end"


def parse_platform_annot(path: Path) -> pd.Series:
    """Return probe id -> gene symbol from a GEO .annot.gz platform file.

    These files are SOFT format: a '^Annotation' line, '!'-prefixed metadata,
    '#'-prefixed column descriptions, then the real table fenced between
    '!platform_table_begin' and '!platform_table_end'. Reading the first
    non-comment line as the header picks up '^Annotation' instead, so the fence
    is what we key on.

    Where a probe maps to several genes ('MIR4640///DDR1') the first symbol is
    taken. That is arbitrary but applied identically to every platform, which is
    what matters for comparing cohorts.
    """
    rows: list[tuple[str, str]] = []
    header: list[str] | None = None
    sym_idx = id_idx = None
    in_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith(PLATFORM_BEGIN):
                    in_table = True
                continue
            if line.startswith(PLATFORM_END):
                break
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                lower = [h.lower() for h in header]
                id_idx = 0
                for cand in ("gene symbol", "gene_symbol", "symbol"):
                    if cand in lower:
                        sym_idx = lower.index(cand)
                        break
                if sym_idx is None:
                    raise ValueError(f"{path.name}: no gene symbol column in {header}")
                continue
            if len(parts) <= sym_idx:
                continue
            probe, symbol = parts[id_idx].strip(), parts[sym_idx].strip()
            if probe and symbol and symbol not in {"", "-", "---"}:
                rows.append((probe, symbol.split("///")[0].strip()))

    if header is None:
        raise ValueError(f"{path.name}: no '{PLATFORM_BEGIN}' fence found")

    mapping = pd.Series(dict(rows), name="symbol")
    mapping.index.name = "probe"
    return mapping
