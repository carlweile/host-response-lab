# host-response-lab — project rules

## The claim this repo is allowed to make

> Published host-response signatures report AUCs around 0.93–0.95. Here is the
> same kind of model, on public data, scored two ways in one run, so you can see
> how much of that number is cohort-specific.

Anything stronger is out of scope. In particular this repo may **not** claim a
new signature, a better classifier, or clinical usefulness. FDA-cleared products
(Inflammatix TriVerity, MeMed BV) already exist, and a published comparison of 28
signatures across 51 datasets already exists. See README *Prior art*.

## Non-negotiables

1. **No feature selection outside the cross-validation loop.** Choosing genes on
   the full dataset and then cross-validating is the standard way this experiment
   gets accidentally rigged. If you add selection, it goes inside the fold.
2. **Labels are guilty until verified.** `config/cohorts.json` carries published
   class counts per cohort. If our parse disagrees, our parse is wrong. Never
   adjust `expected_counts` to match the code; fix the rules or mark the cohort
   unverified.
3. **Exclusions are evaluated before class rules.** "bacterial and viral
   co-infection" contains the word "bacterial". `tests/test_honest_splits.py`
   pins this and demonstrates the bug if the order flips.
4. **Report the gap, not the best number.** Every result states the random-CV AUC
   and the leave-one-cohort-out AUC together. Reporting only the higher one is
   the failure mode this whole repo exists to measure.
5. **Tests must go red on broken code.** New tests construct the bug they guard
   against and assert the guard fires.

## Conventions

- CPU only. Nothing here needs the GPU; if something seems to, question it first.
- `data/` is git-ignored and fully reproducible from `config/cohorts.json`.
- `results/` is committed, so numbers stay comparable across runs.
- New files target ~230 lines; 500 is the hard ceiling.
- Every work session appends one dated entry to `LABBOOK.md` before it ends:
  question, what was run, what was measured, what we now believe, where it lives.

## Known soft spots

- Co-infections and non-infectious inflammation are excluded rather than
  modelled. This flatters the numbers and is a real limitation, not a detail.
- Within-cohort rank normalization (`--rank-normalize`) assumes you hold a whole
  cohort at once, which a bedside test does not. Report both settings.
- Microarray data dominates. Do not extrapolate to RNA-seq without checking.
