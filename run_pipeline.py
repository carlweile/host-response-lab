"""Entry point: fetch -> review labels -> train under both schemes.

    python run_pipeline.py --fetch
    python run_pipeline.py --review
    python run_pipeline.py --train [--rank-normalize | --within-sample-normalize]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

from src import fetch_geo, harmonize, labels as lab, parse_geo, splits, train

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
CONFIG = ROOT / "config" / "cohorts.json"


def load_config() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def cmd_fetch(cfg: dict, force: bool = False) -> int:
    failures = 0
    for entry in cfg["cohorts"]:
        acc = entry["accession"]
        path, msg = fetch_geo.fetch_series(
            acc, DATA, force, prefer_platform=entry.get("prefer_platform")
        )
        print(f"  {acc:<12} series   {msg}")
        if path is None:
            failures += 1
            continue
        try:
            series = parse_geo.parse_series_matrix(path, acc)
        except Exception as exc:
            print(f"  {acc:<12} parse    FAILED: {exc}")
            failures += 1
            continue
        print(f"  {acc:<12} parsed   {series.expression.shape[0]} probes x "
              f"{series.expression.shape[1]} samples, platform {series.platform or '?'}")
        if series.platform:
            _, pmsg = fetch_geo.fetch_platform(series.platform, DATA, force)
            print(f"  {acc:<12} platform {series.platform} {pmsg}")
    return failures


def _load_all(cfg: dict) -> dict[str, tuple[parse_geo.Series, lab.LabelResult]]:
    loaded = {}
    for entry in cfg["cohorts"]:
        acc = entry["accession"]
        path = DATA / "series" / f"{acc}_series_matrix.txt.gz"
        if not path.exists():
            print(f"  {acc}: not downloaded, run --fetch first")
            continue
        try:
            series = parse_geo.parse_series_matrix(path, acc)
        except Exception as exc:
            print(f"  {acc}: parse failed: {exc}")
            continue
        want = entry.get("prefer_platform")
        if want and series.platform != want:
            raise SystemExit(f"{acc}: parsed platform {series.platform}, config requires {want}; "
                             "delete the cached matrix and re-run --fetch")
        result = lab.assign_labels(series, cfg["label_rules"], entry.get("expected_counts"),
                                   entry.get("exclude_samples", {}).get("ids"))
        loaded[acc] = (series, result)
    return loaded


def cmd_review(cfg: dict) -> int:
    loaded = _load_all(cfg)
    RESULTS.mkdir(exist_ok=True)
    chunks, summary = [], []
    for acc, (_series, result) in loaded.items():
        text = lab.review_text(result)
        print(text)
        chunks.append(text)
        summary.append({
            "cohort": acc,
            "bacterial": result.counts.get(lab.BACTERIAL, 0),
            "viral": result.counts.get(lab.VIRAL, 0),
            "excluded": result.counts.get(lab.EXCLUDED, 0),
            "verified": result.verified,
            "verdict": result.verdict,
        })
    if summary:
        df = pd.DataFrame(summary)
        print("\nLabel summary\n" + df.to_string(index=False))
        (RESULTS / "label_review.md").write_text(
            f"# Label review ({date.today().isoformat()})\n\n"
            "Parsed from GEO sample characteristics, checked against published counts\n"
            "where we have them. A cohort marked unverified is still usable but its\n"
            "labels have not been confirmed by an outside number.\n\n```\n"
            + "\n\n".join(chunks) + "\n```\n\n## Summary\n\n"
            + df.to_markdown(index=False) + "\n",
            encoding="utf-8",
        )
        print(f"\nwrote {RESULTS / 'label_review.md'}")
    return 0


def cmd_train(cfg: dict, rank_normalize: bool, within_sample: bool = False) -> int:
    loaded = _load_all(cfg)
    if not loaded:
        print("no cohorts available")
        return 1
    lab.refuse_mismatches([result for _, result in loaded.values()])

    per_cohort, label_map, notes = {}, {}, []
    for acc, (series, result) in loaded.items():
        usable = result.usable
        if len(usable) < 20:
            notes.append(f"{acc}: only {len(usable)} labelled samples, skipped")
            continue
        annot_path = DATA / "platforms" / f"{series.platform}.annot.gz"
        if not annot_path.exists():
            notes.append(f"{acc}: no platform annotation for {series.platform}, skipped")
            continue
        try:
            mapping = parse_geo.parse_platform_annot(annot_path)
        except Exception as exc:
            notes.append(f"{acc}: platform annotation unusable ({exc}), skipped")
            continue
        expr, logged = harmonize.maybe_log2(series.expression)
        genes = harmonize.to_symbols(expr, mapping)
        if genes.empty:
            notes.append(f"{acc}: no probes mapped to symbols, skipped")
            continue
        per_cohort[acc] = genes
        label_map[acc] = usable
        notes.append(f"{acc}: {genes.shape[0]} genes, {len(usable)} labelled samples"
                     f"{', log2 applied' if logged else ''}"
                     f"{'' if result.verified else ', labels UNVERIFIED'}")

    for n in notes:
        print("  " + n)
    if len(per_cohort) < 2:
        print("\nneed at least 2 usable cohorts for leave-one-cohort-out")
        return 1

    X, y, cohort, build_notes = harmonize.build_matrix(per_cohort, label_map)
    for n in build_notes:
        print("  " + n)
    if rank_normalize:
        X = splits.rank_normalize_within_cohort(X, cohort)
        print("  within-cohort rank normalization: ON")
    elif within_sample:
        X = splits.rank_normalize_within_sample(X)
        print("  within-sample rank normalization: ON (one patient at a time)")

    print(f"\n  matrix: {X.shape[0]} samples x {X.shape[1]} genes "
          f"({int(y.sum())} bacterial, {int((1 - y).sum())} viral) "
          f"across {cohort.nunique()} cohorts\n")

    results = train.run_both_schemes(X, y, cohort)
    detail = pd.DataFrame([r.__dict__ for r in results])
    summary = train.summarize(results)

    print("Per-fold\n" + detail.to_string(index=False))
    print("\nSummary\n" + (summary.to_string(index=False) if not summary.empty else "none"))

    RESULTS.mkdir(exist_ok=True)
    tag = "ranknorm" if rank_normalize else "withinsample" if within_sample else "raw"
    detail.to_csv(RESULTS / f"folds_{tag}.csv", index=False)
    if not summary.empty:
        summary.to_csv(RESULTS / f"summary_{tag}.csv", index=False)
    meta = {
        "date": date.today().isoformat(),
        "rank_normalize": rank_normalize,
        "within_sample_normalize": within_sample,
        "cohorts": {k: len(v) for k, v in label_map.items()},
        "n_samples": int(X.shape[0]),
        "n_genes": int(X.shape[1]),
        "notes": notes + build_notes,
    }
    (RESULTS / f"run_{tag}.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"\nwrote {RESULTS / f'summary_{tag}.csv'} and {RESULTS / f'run_{tag}.json'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true", help="download GEO series + platforms")
    ap.add_argument("--review", action="store_true", help="print and save label review")
    ap.add_argument("--train", action="store_true", help="run both validation schemes")
    ap.add_argument("--rank-normalize", action="store_true",
                    help="rank-normalize genes within each cohort before training")
    ap.add_argument("--within-sample-normalize", action="store_true",
                    help="rank each sample's own genes (reads no other patient); "
                         "exclusive with --rank-normalize")
    ap.add_argument("--force", action="store_true", help="re-download cached files")
    args = ap.parse_args(argv)

    if args.rank_normalize and args.within_sample_normalize:
        ap.error("--rank-normalize and --within-sample-normalize are mutually exclusive")

    if not (args.fetch or args.review or args.train):
        ap.print_help()
        return 0

    cfg = load_config()
    rc = 0
    if args.fetch:
        print("== fetch ==")
        rc |= cmd_fetch(cfg, args.force)
    if args.review:
        print("\n== label review ==")
        rc |= cmd_review(cfg)
    if args.train:
        print("\n== train ==")
        rc |= cmd_train(cfg, args.rank_normalize, args.within_sample_normalize)
    return rc


if __name__ == "__main__":
    sys.exit(main())
