# host-response-lab

**How much of a published diagnostic AUC survives a hospital the model has never seen?**
A small, public, CPU-only reproduction study on blood gene-expression data, built
so that one command computes the optimistic number and the honest number from the
same data and prints the gap between them.

> Constitution: [CORE.md](CORE.md) (question, success criteria, kill criterion).
> Diary: [LABBOOK.md](LABBOOK.md) (every run, dated, with its pre-registered
> rules). Write-up for a general reader: [`reports/2026-11-08-writeup.html`](reports/2026-11-08-writeup.html).

## 1. The question

When a child or an adult arrives in an emergency department with a fever, the
hard decision is *bacterial or viral?* Antibiotics help one and not the other.
A line of research since about 2007 reads the patient's own immune response in
blood (which genes are switched on, the "host response") and trains a classifier
on it. Published papers report areas under the ROC curve (AUC) of **0.93 to
0.95**; two such tests are FDA-cleared.

Those AUCs are almost always measured by splitting one study's patients at
random into training and test sets. Every patient in the test set then shares a
hospital, a lab protocol and a microarray chip with the training set. A model is
free to learn those local quirks and gets credit for them. A real test is used
on a hospital it has never seen.

**So: of a 0.93, how much is biology that travels, and how much is the hospital?**

That gap is measurable from public data, it is reported inconsistently in the
literature, and when this project started (October 2026) we found no public code that
recomputes it in one run. That is the whole contribution: not a better
classifier, not a new gene signature, not a clinical claim. An engineering and
reproducibility piece. (Prior art, including the FDA-cleared products and the
large benchmark that deliberately does *not* measure this, is in section 8.)

## 2. The design in plain words

Five public cohorts from GEO, each from a different study and hospital, all with
bacterial/viral labels, all on microarray chips. After label parsing: **458
patient arrays, 12,040 genes shared by every chip.**

| cohort | study | chip family | age group | bacterial | viral |
|---|---|---|---|---|---|
| GSE63990 | Tsalik 2016 | Affymetrix | mixed, 6–88 y | 70 | 115 |
| GSE42026 | Herberg 2013 | Illumina | children | 18 | 41 |
| GSE60244 | Suarez 2015 | Illumina | adults | 22 | 71 |
| GSE40396 | Hu 2013 | Illumina | children | 8 | 22 |
| GSE6269 | Ramilo 2007 | Affymetrix | mostly children (3 of 91 aged 18–36) | 73 | 18 |

The same two models (an L2 logistic regression and a LightGBM gradient-boosted
tree, both class-weighted, no gene selection anywhere) are trained and scored
two ways on exactly the same matrix:

- **Random cross-validation.** Patients of all five hospitals are shuffled and
  split 5-fold. This is the number papers usually show.
- **Leave-one-cohort-out.** Train on four hospitals, test on the fifth, repeat
  five times, average. This is the number a new hospital would see.

The headline output is not an AUC. It is the **difference** between those two.

Three habits from wet-lab work are enforced in code, because this kind of study
is easy to rig by accident:

1. **Labels are guilty until verified.** The parser's class counts per cohort
   must equal the counts printed in each study's paper, or the pipeline refuses
   to train. Counts are never edited to match the code.
2. **Every rule is written before the run.** Each experiment's thresholds and
   verdict words are written into the LABBOOK first; the run does not get to
   change them. Exceptions are disclosed in the entry concerned. The
   pre-registration commits live in the private working repo, so their timing
   cannot be checked from this repository's history.
3. **Tests go red on broken code.** Each of the 211 tests constructs the bug it
   guards against (a co-infection parsed as "bacterial", a test cohort leaking
   into training, a transform that reads other patients) and asserts the guard
   fires.

## 3. The four results

### Result 1: about 7 to 12 AUC points of the headline are hospital-specific

Raw expression scale, 458 arrays. Random-CV at fold seed 0, then the spread over
10 fold seeds.

| model | random-CV AUC | leave-one-cohort-out AUC | gap | gap over 10 fold seeds |
|---|---|---|---|---|
| logistic L2 | 0.957 | 0.881 | 0.076 | 0.073 ± 0.004 (0.065–0.076) |
| LightGBM | 0.954 | 0.827 | 0.127 | 0.123 ± 0.006 (0.113–0.133) |

In words: both models look like a 0.95 test when the hospitals are mixed.
Tested on an unseen hospital, the logistic model is a 0.88 test and the tree
model a 0.83 test. The gap is positive in all 10 fold seeds and, read cohort by
cohort, positive in **all five hospitals for both models** (smallest single-cohort
gap 0.033 logistic, 0.011 LightGBM). The per-cohort held-out AUC ranges from
0.83 to 0.94 (logistic) and 0.70 to 0.90 (LightGBM).

The more flexible model pays more. LightGBM also carries its own randomness: over
10 model seeds its held-out AUC moves between 0.780 and 0.861 (SD 0.023), so its
gap should be read as "clearly larger than logistic, size uncertain". The
logistic model is deterministic.

Sources: `results/summary_raw.csv`, `results/folds_raw.csv`,
`results/seedsweep_raw_summary.csv`, `results/seedsweep_lgbm_modelseed.csv`;
per-cohort gaps `results/gapdiag_cohorts.csv`.

### Result 2: a bedside-legal transform gives back 1 point (logistic) or 6 points (LightGBM)

Most published pipelines normalise a gene against the other patients of the same
cohort. A bedside test sees one patient at a time, so that is not allowed here.
Instead each patient's 12,040 genes are ranked against *each other* only
("within-sample" ranking; no other patient is read, in training or test).

| model | random-CV AUC | leave-one-cohort-out AUC | gap over 10 fold seeds |
|---|---|---|---|
| logistic L2 | 0.963 | 0.895 | 0.065 ± 0.003 (0.061–0.068) |
| LightGBM | 0.958 | 0.888 | 0.067 ± 0.004 (0.062–0.074) |

The gap does not vanish. For logistic it shrinks by about one point; for
LightGBM by about six, and that six survives LightGBM's own randomness: over 10
paired model seeds the within-sample run beat the raw run every time (mean gain
+0.064, range +0.027 to +0.112; within-sample held-out 0.891 ± 0.006).

So the leak-free answer to the question in section 1 is: **of a random-CV AUC of
about 0.96, about 0.88 survives an unseen hospital on the raw scale and about
0.89 survives with per-patient ranking.** A whole-cohort rank normalisation
*would* make the gap look closed, and was shown on synthetic data to be a
transductive artefact (it reads the test set); its numbers are not quoted.

Sources: `results/summary_withinsample.csv`,
`results/seedsweep_withinsample_summary.csv`,
`results/seedsweep_lgbm_modelseed_withinsample.csv`.

### Result 3: frozen published signatures have no random-CV bonus, which is what fitting predicts

If the gap comes from a model learning hospital quirks, a gene score that is
**not fitted** on these data should show no gap. Two published scores were used
exactly as printed in their papers and supplements, nothing re-fitted: Herberg
2016 (two genes) and Sweeney 2016 (seven genes).

| score | cohorts | random-CV AUC | leave-one-cohort-out AUC | gap (10 fold seeds) | our logistic / LightGBM held out, same cohorts |
|---|---|---|---|---|---|
| Herberg 2-gene | 3 | 0.876 | 0.882 | −0.004 ± 0.005 | 0.846 / 0.847 |
| Sweeney 7-gene | 5 | 0.837 | 0.912 | −0.075 ± 0.002 | 0.881 / 0.827 |
| Sweeney, only its 2 independent cohorts | 2 | — (not computed on 2 cohorts) | 0.901 | — | 0.885 / 0.851 |

The frozen scores have **no** random-CV advantage (the pooled number is the same
or lower, because each cohort puts the fixed score on its own offset). Held out,
they do about as well as or slightly better than the trained logistic model. This fits the reading that
the Result-1 gap comes from fitting; it does not prove it. Limits that matter:
Herberg's gene FAM89A has no probe on two of the five chips, so it is scored on
three cohorts; Sweeney's own paper lists three of our five cohorts among those its
genes were chosen on, so only two are truly independent of it. One of the Sweeney
authors confirmed our reading of their score by e-mail (LABBOOK 2026-10-23).

Sources: `results/frozen_signatures_summary.csv`, `config/signatures.json`.

### Result 4: the drop is not the chip family, and class weighting is not where it lives

The held-out AUC differs by hospital. Does the size of the drop follow something
we can name? Three candidates were pre-registered: chip family (Affymetrix vs
Illumina), age group (children vs adults), class mix.

- **Chip family: no.** On all 20 ordered train-on-one / test-on-one cohort pairs
  (logistic, raw), same-family pairs transferred **worse** than cross-family
  pairs: mean AUC 0.760 vs 0.821, difference −0.061. With every training cohort
  cut to 30 arrays (20 draws) the difference is −0.062, so it is not that the
  Illumina cohorts are small. Within-sample ranking −0.034; LightGBM −0.11 at
  each of three seeds. Verdict by the pre-registered rule: "does not track".
- **Class weighting: not where the drop lives.** Leave-one-cohort-out with class
  weighting switched off moves the mean held-out AUC by +0.0005 (logistic) and
  +0.003 (LightGBM, 10 model seeds; per-seed range −0.029 to +0.032, SD 0.026).
  A first 3-seed read gave +0.0195 and "undecided"; the 10-seed run showed that
  was seed luck.
- **Age: not testable from public sources.** GSE63990 mixes ages 6 to 88 and
  neither GEO (0 of 280 arrays carry an age) nor the paper gives per-patient
  ages, so the children/adults grid cannot be built with these five cohorts.

Before any of this ran, the new code had to reproduce all ten on-file held-out
AUCs to within 0.005 (it did, largest difference 0.0006).

Sources: `results/transfer_matrix.csv`, `results/transfer_sizematched.csv`,
`results/balance_control_10seeds.csv`, `results/run_transfer.json`,
`results/run_balance_10seeds.json`; LABBOOK entries 2026-11-05 and 2026-11-06.

## 4. What the drop is NOT (closed candidates)

Each line was a test with its own control, pre-registered unless its entry says otherwise; the LABBOOK entry has
the rule and the numbers.

| candidate | how it was tested | verdict |
|---|---|---|
| a lucky random split | 10 fold seeds, raw and within-sample | no: gap > 0.06 in every seed, both models |
| LightGBM's own randomness | 10 model seeds, raw and within-sample, paired | no: moves LightGBM's level, not the sign of the gap or the 6-point give-back |
| 7 re-run arrays of GSE63990 | removed; re-run (`results/dedup_sensitivity_raw.csv`) | no: removing them made the gap slightly larger |
| chip family (Affymetrix vs Illumina) | 20 single-cohort transfer pairs + size-matched control | no: same-family pairs transfer worse |
| class weighting | leave-one-cohort-out with weighting off, 10 LightGBM seeds | no: moves the mean by < 0.01 for both models |
| a whole-cohort normalisation "fix" | synthetic positive controls | the fix is the leak: it reads the test set; not quoted |
| sample type (PBMC vs whole blood) | GSE6269 is the only PBMC cohort; the 4 whole-blood cohorts re-run without it, same 12,040 genes (`results/whole_blood_check.csv`) | logistic: no detectable effect (held-out mean 0.857 vs 0.866, gap 0.081 vs 0.076); LightGBM: no verdict, see below |
| age group | GEO fields + paper tables for GSE63990 | cannot be tested with these cohorts |
| a bug in our pipeline | 3 AI code reviews (two by other vendors' models), 211 red-on-bug tests, 4 identical end-to-end re-runs incl. one by the author from a fresh GEO download | no reviewer found a bug that changes a committed number |

The sample-type row carries a limit and a side finding. The limit: GSE6269 is
also the mostly pediatric, 80 %-bacterial, Affymetrix cohort, so with one PBMC cohort
"PBMC" cannot be separated from "that one cohort". The side finding, LightGBM
only: without GSE6269 in the training set the tree model becomes anti-predictive
on the other young-children cohort, GSE40396 (held-out AUC 0.22 to 0.48 over 3
seeds, from 0.84 with it); the other three whole-blood cohorts barely move.
Read as training-set composition, not as sample type.

### Where the cause could still hide

Checked against the per-sample fields in the cached GEO series matrices (08.11.2026).

| group | candidate | status |
|---|---|---|
| tested, ruled out | lucky split · chip family · class weighting · the 7 re-run arrays · sample type (logistic) | see the table above |
| checkable from GEO metadata, not yet tested | blood tube: PAXgene (GSE42026, GSE63990) vs Tempus (GSE40396, GSE60244); GSE6269 is PBMC, no tube | not testable on its own: tube is confounded with cohort (each cohort used one tube type) |
| | lab batch: GSE63990 `batch` (0/1), GSE40396 `sentrix barcode` + `hybridization date` (3 dates) | GSE63990 batch split (run 8 Nov): verdict **undecided**. Batch is fully visible in the values (batch-AUC 1.000). Train batch 0 / test batch 1: AUC 0.849 vs 0.957 for 20 size- and class-matched mixed splits, Δ 0.108, below all 20. Reverse: 0.885 vs 0.924, Δ 0.039, inside the control spread. Per-patient ranking did not close it. Power limit: the small side has 39 arrays, 12 bacterial. **Deciding test (run 8 Nov, 1000 size- and class-matched mixed-split draws): still undecided.** Train batch 0 / test batch 1: AUC 0.849, null mean 0.961, 4 of 1000 draws at or below it (0.4 %). Train batch 1 / test batch 0: AUC 0.885, null mean 0.920, 121 of 1000 at or below it (12.1 %, interval 10.2-14.3 %). So the batch costs AUC clearly in one direction only. GSE40396 hybridization dates: underpowered, not run (bacterial arrays per date 0, 0 and 8; all 8 bacterial arrays sit on one date). The drop stays unattributed |
| | sex: in GEO for 4 of 5 (all but GSE63990) | not tested |
| | age: in GEO for 3 of 5 (GSE40396 months, GSE60244 years, GSE6269 `Age: 10 y`); not for GSE42026, GSE63990 | the children/adults grid still fails on GSE63990 |
| not in the public data | days since symptom onset · illness severity · ward vs ICU · how the bacterial/viral label was adjudicated (GSE6269 lists `illness` and antibiotic `treatment`, the others nothing) | cannot be tested from GEO |

**What remains:** with five cohorts the drop is a *pattern*, not a *cause*. It
stays unattributed. The batch split of the one cohort with two documented array
batches has been run (LABBOOK 2026-11-07 PLAN onward): a
within-hospital batch costs AUC in one direction only, so it is a partial
candidate, not a cause.

## 5. Honest limits

- **Flattered by exclusions.** Co-infections and non-infectious inflammation are
  excluded, not modelled. Several published papers do the same; it makes
  everyone's numbers, including these, look better than an emergency department.
- **Five cohorts, 458 arrays.** Roughly a tenth of the largest published
  comparison. Enough to measure a gap, not enough to name its cause.
- **Microarray only.** Nothing here is checked on RNA-seq.
- **Labels are machine-parsed.** Published counts match for all five cohorts, but
  only two of the five count checks were done by a human (the author); the other
  three are AI-verified against the papers. No author of an original cohort has
  read the label table.
- **Which GSE63990 arrays are the re-runs is inferred**, not proven (GEO has no
  patient ids); a correlation check supports the choice.
- **Not a clinical failure rate.** A gap across archived cohorts is not a
  prospective validation. The FDA-cleared TriVerity test *was* validated
  prospectively (1,222 patients, 22 sites, Nature Medicine 2025) and its numbers
  are not comparable to these.
- **Within-sample ranking still needs the whole gene profile of one patient**,
  which a 2- or 7-gene bedside assay does not have.

## 6. Who did what

This repo was built with an AI coding assistant (Claude Code) directed by the
author, a molecular biologist. The honest split, from the LABBOOK record:

| item | done by | judged by |
|---|---|---|
| question, scope, the frozen five-cohort set, every label-rule decision, every go/no-go | the author | the author |
| pipeline code, tests, all runs, LABBOOK entries | the AI, directed by the author | 3 outside AI code reviews (gpt-6-astra, Gemini 3.1 Pro, a blind Claude agent); every finding fixed with a red test, rejected with a reason checked on the data, or left with its limit stated |
| published class counts | the AI from the papers; GSE40396 and GSE6269 re-checked by the author | author for 2 of 5, AI-only for 3 of 5 |
| every result sentence | AI drafted, author approved | a blind fresh agent re-derives the numbers from the CSVs before the sentence is accepted; pre-registered thresholds decide the verdict word |
| reproducibility | — | the author re-ran the whole pipeline from a fresh GEO download on their own machine: all four headline AUCs identical to the last digit |
| Sweeney score reading | the AI from paper + supplement + authors' R code | confirmed by one of the score's authors by e-mail |

Still open before anyone should rely on this: a human code review, and an author
of one of the original cohorts reading the label table.

## 7. Reproduce it

```
python -m pip install -r requirements.txt   # numpy, pandas, scikit-learn, scipy, lightgbm, requests
python run_pipeline.py --fetch              # download the 5 GEO series into data/ (once, ~120 MB cached)
python run_pipeline.py --review             # label candidates + published-count verification (must say verified 5/5)
python run_pipeline.py --train              # Result 1: both schemes, writes results/ (~1.5 min on a desktop CPU)
python run_pipeline.py --train --within-sample-normalize   # Result 2
python run_frozen_signatures.py             # Result 3
python run_transfer_matrix.py              # Result 4: chip family + 3-seed class weighting (~4 min)
python run_transfer_matrix.py --balance-seeds 10           # Result 4: class weighting, 10 LightGBM seeds (~10 min)
python tools/geo_age_check.py               # Result 4: age fields in GEO
python run_whole_blood_check.py             # sample-type check (~4 min)
python run_batch_split.py                   # GSE63990 batch split
python run_batch_deciding.py                # deciding batch test, 1000 draws (~11 min)
python run_seed_sweep.py                    # 10 fold seeds (~20 min)
python -m pytest -q                         # 211 tests
```

CPU only; no GPU is used anywhere. `data/` is git-ignored and rebuilt by
`--fetch`. `results/` is committed so every number above can be diffed against a
re-run. Each run reports its own reproduction gate before writing anything new.

## 8. Prior art, and the one paragraph that justifies this repo

This field is mature and commercial. Know this before calling the work novel.

| what exists | status |
|---|---|
| **Inflammatix TriVerity**, 29-gene host-response test | FDA-cleared Jan 2025; prospective validation Liesenfeld et al., Nature Medicine 2025 ([doi](https://doi.org/10.1038/s41591-025-03933-y)): 1,222 patients, 22 emergency departments; bacterial AUROC 0.83, viral 0.91 in the adjudicated population |
| **MeMed BV**, 3-protein test | FDA-cleared |
| **Bodkin et al. 2022**, Genome Medicine | 28 signatures, 51 datasets, 4,589 subjects. The largest comparison; see below |
| **Chawla et al. 2022**, Cell Systems | 30 signatures vs a 17,105-profile compendium; measures robustness and cross-reactivity |
| **MetaIntegrator** (Khatri lab) | R package with leave-one-dataset-out built in; the method is not new |
| **bvnGPS**, 16 cohorts incl. these five | has both a held-out and a random-split mode, but as two invocations plus a scoring script, one split not CV |

Bodkin et al. answers a different question. They fitted a fresh lasso logistic
regression **inside each dataset separately** and cross-validated within it,
stating: *"Creating dataset-specific models overcomes batch effects since each
signature is optimized in each dataset."* Sound for their question ("how good is
this gene list when tuned locally"), but it removes the batch effect by
construction instead of measuring it, and the paper reports no discovery-versus-
external drop. The second-largest benchmark measures robustness instead.
Well-resourced groups do cross-cohort validation in primary papers, but their
pipelines are not public.

**Therefore not the contribution here:** a new signature, a better AUC, a
clinical product, or more data. **The contribution that survives scrutiny:** a
small, public, runnable implementation where one command computes the optimistic
number and the cohort-held-out number from the same data, the label parsing
checks itself against published counts, and every verdict was written down
before the run. Say it that way.

## Labels, in detail

| cohort | platform | bacterial | viral | excluded | published count source |
|---|---|---|---|---|---|
| GSE63990 | GPL571 | 70 | 115 | 95 | Tsalik 2016 Sci Transl Med, Table 1 (AI-verified); excluded = 88 non-infectious + 7 re-run arrays |
| GSE42026 | GPL6947 | 18 | 41 | 33 | authors' GEO design text (Herberg 2013), Tsalik 2016 Table 1 (AI-verified) |
| GSE60244 | GPL10558 | 22 | 71 | 65 | Suarez 2015 J Infect Dis, Methods (AI-verified) |
| GSE40396 | GPL10558 | 8 | 22 | 35 | Hu 2013 PNAS, Results; febrile children only (checked by the author 2026-10-23) |
| GSE6269 | GPL96 | 73 | 18 | 6 | Ramilo 2007 Blood, Table 1 footnote (checked by the author 2026-10-23) |

GSE63990 holds 280 arrays of 273 patients in two batches sharing 7 re-run
samples; the last 7 arrays carry exactly the surplus over the published counts
and are excluded (`config/cohorts.json`, guarded by a test). Exclusion rules run
before class rules, so "bacterial and viral co-infection" never becomes
"bacterial".

See `LABBOOK.md` for what was run and measured, in date order, with every
pre-registered rule next to its result.
