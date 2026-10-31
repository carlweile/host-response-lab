# host-response-lab

**How much of a published "bacterial or viral?" blood-test score survives a
hospital it has never seen?** A small, fully re-runnable reproduction study on
public gene-expression data.

## The question

When a patient arrives with a fever, doctors often cannot tell quickly whether a
bacterial or a viral infection causes it, and antibiotics only help against
bacteria. Research groups have published "host-response signatures": patterns
in the patient's own blood gene activity that tell the two apart. Papers report
scores (AUCs) of around 0.93 to 0.95 (literature figure; see *Prior art*, not
re-measured here).

Most of those numbers come from testing a model on patients from the **same
studies** it learned from. This repo asks how much of the score is left when the
model is tested on a **whole hospital cohort it never saw** during training.

This is a reproduction and validation study, not a new test. Commercial tests
already exist (see *Prior art* below).

## Design in one picture

```
5 public patient cohorts (GEO), 458 blood samples, 12,040 shared genes
                │
   ┌────────────┴─────────────┐
   │                          │
random split               leave-one-cohort-out
(patients from every       (train on 4 cohorts,
 cohort mixed into          test on the 5th, then
 train and test)            rotate through all 5)
   │                          │
   └──── same model, same data, same run ────┘
                │
     the GAP between the two scores = the headline
```

- **Two models:** a simple one (logistic regression) and a flexible one
  (LightGBM, gradient-boosted trees).
- **Two published signatures, used exactly as printed** (nothing re-fitted):
  Herberg 2016 (2 genes) and Sweeney 2016 (7 genes).
- **Labels are checked against the papers.** The pipeline counts bacterial and
  viral patients per cohort and refuses to train if the counts disagree with the
  numbers published for that cohort. All 5 cohorts pass (the author read 2 of the
  5 published counts in the papers; 2 were checked by AI only, see *Honest limits*).

## Key numbers, in plain words

**What an AUC means here:** pick one bacterial and one viral patient at random.
An AUC of 0.95 means the model gives the bacterial patient the higher
"bacterial" score 95 times out of 100. 0.50 would be a coin flip.

| model | random split | new hospital (held out) | gap | gap over 10 random seeds |
|---|---|---|---|---|
| logistic regression | 0.957 | 0.881 | 0.076 | 0.073 ± 0.004 |
| LightGBM | 0.954 | 0.827 | 0.127 | 0.123 ± 0.006 |

In words:

1. **A random split looks excellent, a new hospital looks clearly worse.** The
   simple model drops from 0.96 to 0.88; the flexible model drops further, to
   0.83. In the worst case the flexible model reaches only 0.70 on one held-out
   cohort (GSE6269), where the simple model still reaches 0.94.
2. **The drop is not luck of the split.** Repeated with 10 different random
   splits, the gap stays between 0.065 and 0.076 for the simple model and
   between 0.113 and 0.133 for the flexible one. Only the split seed was varied
   here; LightGBM's own model seed adds more noise (see *Technical detail*).
3. **A fair per-patient fix closes part of it.** If each patient's
   genes are ranked against that patient's own other genes (something a
   bedside test could also do, because it needs no other patient), the
   held-out scores rise to 0.895 and 0.888. That closes little of the simple
   model's gap and about half of the flexible model's; both end near 0.07
   (0.065 and 0.067 averaged over 10 random-split seeds; the held-out score itself does not move with the split).
4. **The published signatures get no random-split bonus, as expected.** Herberg's
   gap is −0.004; Sweeney even scores lower on the random split (gap −0.075),
   probably because mixing cohorts puts each cohort's scores on its own offset
   (inferred, not tested; gaps are 10-seed means; see *Technical detail*). A frozen score
   learns nothing from our data, so it cannot learn one hospital's quirks.
   Held out, Herberg scores 0.882 (3 cohorts) and Sweeney 0.912 (5 cohorts;
   0.901 on the 2 cohorts its authors did not use to pick its genes). That is
   slightly better than our trained simple model on the same cohorts (0.846 and
   0.881; no confidence intervals), and 3 of Sweeney's 5 cohorts were used to
   pick its genes. This fits the reading
   that the gap comes from fitting to the data. It does not prove it.

Sources: `results/summary_raw.csv`, `results/seedsweep_raw_summary.csv`,
`results/gapdiag_cohorts.csv` (per-cohort numbers), `results/summary_withinsample.csv`,
`results/frozen_signatures_summary.csv`.

## How to re-run it

Tested with Python 3.12. The GEO downloads are on the order of 100 MB. CPU only, no GPU.

```
python -m pip install -r requirements.txt
python run_pipeline.py --fetch        # download the 5 GEO cohorts into data/ (once)
python run_pipeline.py --review       # label counts vs the published counts
python run_pipeline.py --train        # both schemes, both models -> results/
python run_pipeline.py --train --within-sample-normalize
python run_seed_sweep.py              # the 10-seed repeat (about 20 minutes)
python run_frozen_signatures.py       # the two published signatures
python -m pytest -q                   # 127 tests; they need no download
```

The tests are built to fail on broken code: each one constructs the bug it
guards against (for example, a co-infection counted as bacterial, a gene chosen
on the whole dataset before the split, a score that peeks at other patients).

## Honest limits

- **AI-assisted, and not yet checked by a human expert.** The code, the label
  parsing and most checks were written or run by an AI coding assistant
  (Claude Code), directed by the author, a molecular biologist. Three AI code
  reviews found no bug that changes a number on file, but **no human has
  reviewed the code or the label table yet**. For two of the five cohorts the
  published-count check was done by the AI alone; for a third the author accepted
  the AI's comparison without re-reading the paper. Details in *What the author
  checked vs what AI did*.
- **Small.** 5 cohorts and 458 samples, roughly a tenth of the largest published
  comparison. No confidence intervals on the per-cohort numbers; the smallest
  cohort has 8 bacterial patients.
- **Flattered by exclusions.** Patients with co-infections or non-infectious
  inflammation are left out, as in several published papers. A real emergency
  department would be harder.
- **Microarray only, archived data only.** Results may not transfer to RNA-seq,
  and a gap measured on old cohorts is not a prospective clinical failure rate.
- **One label decision is inferred.** One cohort (GSE63990) has 7 re-run arrays
  of patients already in it; which 7 is inferred from the data, not stated by
  GEO.
- **Not novel as a method.** Leave-one-dataset-out validation is already
  packaged elsewhere (MetaIntegrator). The contribution is a small public
  pipeline where one command produces both numbers from the same data, and the
  labels check themselves against the papers.

## How the work was done

`LABBOOK.md` is the lab diary, one dated entry per session: the question, the
plan (thresholds written down **before** each run), what was measured, and
what we now believe. `CORE.md` holds the lab's fixed question, success and stop
criteria. `CLAUDE.md` holds the rules every contributor (human or AI) follows,
such as "no gene selection outside the cross-validation fold" and "always
report both numbers together".

---

# Technical detail

### Result (v1, 2026-10-21)

Five public microarray cohorts (GSE63990, GSE42026, GSE60244, GSE40396, GSE6269
on GPL96), 458 labelled arrays (the 7 GSE63990 re-run arrays removed; see Labels), 12,040 shared genes,
raw scale.

| data | model | random-CV AUC, fold seed 0 | leave-one-cohort-out AUC | gap, fold seed 0 | gap over 10 fold seeds (mean ± SD, range) |
|---|---|---|---|---|---|
| **458 arrays (main)** | logistic L2 | 0.957 | 0.881 | 0.076 | 0.073 ± 0.004 (0.065-0.076) |
| **458 arrays (main)** | LightGBM | 0.954 | 0.827 | 0.127 | 0.123 ± 0.006 (0.113-0.133) |
| comparison: with all 280 GSE63990 arrays (463) | logistic L2 | 0.953 | 0.882 | 0.071 | 0.071 ± 0.004 (0.065-0.076) |
| comparison: with all 280 GSE63990 arrays (463) | LightGBM | 0.954 | 0.839 | 0.115 | 0.111 ± 0.003 (0.106-0.117) |

The comparison rows keep 5 repeat arrays of GSE63990 patients in the data (an
earlier version of this repo); removing them makes the gap slightly *larger*, not
smaller (for LightGBM the change is inside its own model-seed noise, see below). Sources: `results/summary_raw.csv`, `results/seedsweep_raw_summary.csv`;
comparison `results/summary_raw_all280arrays.csv`,
`results/seedsweep_raw_summary_all280arrays.csv`.

In one sentence: of a random-CV AUC of about 0.95 to 0.96 (10-seed means: logistic 0.953 raw /
0.960 within-sample, LightGBM 0.950 / 0.954), about 0.88 survives a hospital cohort the
model never saw on the raw scale (logistic 0.881, gap 0.073 ± 0.004 over 10 fold
seeds; LightGBM 0.827, gap 0.123 ± 0.006, with a wide model-seed range inferred
from the side check below), and about 0.89 survives when each patient is
normalised against their own genes only (within-sample ranking, legal at the
bedside; logistic 0.895, gap 0.065 ± 0.003; LightGBM 0.888, gap 0.067 ± 0.004,
both over 10 fold seeds, LightGBM's model seed not swept on this scale). Roughly
6.5 to 12 AUC points of the headline are hospital-specific. Sources: raw
`results/summary_raw.csv`, `results/seedsweep_raw_summary.csv`,
`results/seedsweep_raw.csv`; within-sample `results/summary_withinsample.csv`,
`results/seedsweep_withinsample_summary.csv`, `results/seedsweep_withinsample.csv`.

How to read it, and what it does not say:

- **Not a seed accident.** The gap stays above 0.06 for both models in all 10
  random-CV fold seeds. LightGBM's own randomness is much larger: its held-out AUC
  moves between 0.780 and 0.861 across 10 model seeds (SD 0.023), so its gap is
  probably "somewhere around 0.09-0.17" (inferred: random-CV was not re-run per
  model seed). Read the LightGBM gap as "clearly larger than logistic, size
  uncertain". The logistic model is deterministic.
- **The size of the gap depends on how you average.** The headline mixes a pooled
  random-CV AUC with an unweighted mean over five cohorts. Compared cohort by cohort
  (random-CV AUC inside that cohort vs its held-out AUC), the gap is smaller, 0.054
  (logistic) and 0.090 (LightGBM) at fold seed 0, but **positive in all five cohorts
  for both models** (smallest: 0.033 logistic, 0.011 LightGBM).
- **Flattered by the exclusions.** Co-infections and non-infectious
  inflammation are excluded, which makes both numbers look better than a real
  emergency department would.
- **Rank-normalising within a cohort is not reported.** That setting reads the whole
  held-out cohort at once; its random-CV number was shown to be a transductive
  artefact (LABBOOK, 2026-10-19).
- It is not a new signature, not a better classifier, and not a clinical claim.

### v1.1: frozen published signatures (2026-10-21)

Question: does the gap also hit gene lists that were chosen **outside** our data?
Two published scores, used exactly as printed, nothing re-fitted: Herberg 2016
(FAM89A − IFI44L) and Sweeney 2016 (geometric mean of IFI27, JUP, LAX1 minus that of
HK3, TNIP1, GPAA1, CTSB). Gene lists, signs and score rules are pinned to the papers
**and their supplements** in `config/signatures.json`: Herberg eTable 3 (directions),
Supplementary Methods eq. 1 (log2 FAM89A − log2 IFI44L) and eTable 6 (validation
sets); Sweeney fig. S3 and Table S1 (the 7 genes and their directions), plus the
authors' MetaIntegrator `calculateScore` code. The supplements confirmed every gene
and sign; the numbers did not change.

Two limits come first, because they change what can be said:

- **FAM89A has no probe on GPL571 or GPL96**, so Herberg is scored on 3 cohorts
  (GSE42026, GSE60244, GSE40396), never on GSE63990 or GSE6269.
- **Sweeney is not fully outside our data.** Its own Table 1 lists GSE6269, GSE40396
  and GSE42026 as cohorts its genes were *chosen on* (also fig. S3). Only GSE60244
  and GSE63990 are independent of it. Herberg chose its genes on a new cohort, not
  on any of ours (possible patient overlap with GSE42026, the same group's earlier
  cohort, is not stated either way); its eTable 6 scored GSE40396, GSE60244 and a
  *different* part of GSE6269 (the GPL570 arrays, not ours; whether they share
  patients with our GPL96 arrays is not checked). GSE42026 is not named
  as a dataset anywhere in the Herberg 2016 main text or supplement.

| score | cohorts | random-CV AUC (pooled) | leave-one-cohort-out AUC | gap (10 fold seeds) | our logistic / LightGBM, held out, same cohorts |
|---|---|---|---|---|---|
| Herberg 2-gene | 3 | 0.876 | 0.882 | −0.004 ± 0.005 | 0.846 / 0.847 |
| Sweeney 7-gene | 5 | 0.837 | 0.912 | −0.075 ± 0.002 | 0.881 / 0.827 |
| Sweeney 7-gene, only its 2 independent cohorts | 2 | — | 0.901 | — | 0.885 / 0.851 |

What it says:

- **A frozen score has no random-CV bonus, and this is expected.** Nothing is
  trained, so nothing can learn one hospital's quirks. Random-CV mixes cohorts in
  every fold, and each cohort puts the score on its own offset, so the pooled
  number comes out the *same or lower* (gap between −0.100 and +0.010 in every
  variant and fold seed). This fits the reading that the v1 gap comes from
  fitting to the data. It does not prove it.
- **Held out, the frozen scores do about as well as our trained models.** Herberg
  0.88 vs logistic 0.85; Sweeney 0.90 on its independent cohorts vs logistic 0.89.
  These are means over 2-3 cohorts at one fold seed, with no confidence intervals.
  The smallest cohort has 8 bacterial patients, so read this as "about the same as
  logistic, and above LightGBM", not as "better".
- **Robust to the reading of the papers.** Raw instead of log2 intensities, a
  probe choice that looks at no expression value, or the other reading of Sweeney's
  "times the ratio of the number of genes" (factor 3/4 on the bacterial term only)
  moves the held-out means by ≤ 0.016. The direction stays the same in every variant.
- **The authors' own software removes part of the pooling penalty, not all of it.**
  MetaIntegrator's `calculateScore` has no gene-count factor and z-scores the score
  inside each dataset. That reads the whole cohort, so it is not a per-patient
  score, and we report it as a separate row: Sweeney random-CV 0.837 → 0.859, gap
  −0.075 → −0.055, held-out unchanged (0.912, as it must be). Herberg: gap −0.003.
  Our guess at why a gap is left: each cohort has a different bacterial/viral mix,
  so centring on the cohort mean still puts cohorts on different offsets
  (inferred, not tested).

Checks: a pre-registered positive control (the run stops if either score falls
below 0.75 on a cohort its own paper reports it working on; it did not fire). On
GSE40396, the one cohort where Herberg's eTable 6 has the same counts as ours
(8 bacterial / 22 viral), our AUC is 0.892 and theirs is 89.2 %. A
blind fresh AI agent re-built both scores from the papers' definitions without
seeing this code, and got every main number back to 3 decimals. 15 tests fail on
broken code (wrong gene list, a gene list that disagrees with the supplement
tables, flipped sign, a dropped gene, a score that reads other patients, a wrong
ratio factor, z-scoring across instead of within cohorts). Sources:
`results/frozen_signatures_summary.csv`, `results/frozen_signatures_cohorts.csv`,
LABBOOK 2026-10-21 (PLANs c193407, 661ea85).

### What the author checked vs what AI did

This repo was built with an AI coding assistant (Claude Code). The author is a
molecular biologist; the code, the label parsing and most checks were written or
run by AI and then accepted by the author. Who did what, from the record in
`LABBOOK.md`:

| item | done by | judged by |
|---|---|---|
| Question, scope, frozen five-cohort set, every rule change (e.g. dropping afebrile children in GSE40396) | the author decided | the author |
| Pipeline code, tests, all runs | Claude (AI), directed by the author | 3 outside AI code reviews (below); no human code review yet |
| Published counts GSE42026 | Claude (AI), pinned 2026-10-21 to the authors' GEO design text, cross-checked with Tsalik 2016 Table 1 | **AI-verified, not human** |
| Published counts GSE63990 | Claude (AI), 2026-10-21, Tsalik 2016 Table 1 (70 bacterial / 115 viral subjects); the 7 re-run arrays excluded by the author's decision | **AI-verified, not human**; which 7 arrays are the re-runs is inferred (see Labels) |
| Published counts GSE60244 | two blind AI agents extracted them from the papers, Claude compared | accepted by the author 2026-10-19 |
| Published counts GSE40396 (Hu 2013, Results first paragraph) | two blind AI agents extracted them from the paper, Claude compared | accepted by the author 2026-10-19; **checked by the author 2026-10-23: Hu 2013 counts 22 viral / 8 bacterial match** |
| Published counts GSE6269 (Ramilo 2007, Table 1 footnote) | Claude (AI), 2026-10-20 and re-read 2026-10-21 | **checked by the author 2026-10-23: Ramilo 2007 Table 1 counts 73 bacterial / 18 viral match** |
| The one-sentence answer | Claude drafted, the author approved | a blind fresh AI agent re-derived all 10 numbers from the CSVs (same model family) |
| 10-seed sweep, probe-choice control | Claude, pre-registered in the LABBOOK before running | the pre-registered rules |
| Code review | gpt-6-astra, Gemini 3.1 Pro, a blind Claude agent, same code bundle | every finding is fixed (all but one with a test that fails on the old code), rejected with a reason checked on the data, or left with its limit stated: LABBOOK 2026-10-21, `reports/review/2026-10-21/`. No reviewer found a bug that changes a number on file |
| Prior-art sweep | Claude (AI), 2026-10-10 | not re-checked by a human |
| v1.1 signature gene lists, signs, score rules, which cohorts each paper used | Claude (AI), 2026-10-21, from the PMC full texts and the supplements (Herberg eTables 3/6, Sweeney fig. S3 / Table S1) plus MetaIntegrator `calculateScore` (`config/signatures.json`) | **AI-verified, not human** |
| v1.1 scores and AUCs | Claude, pre-registered in the LABBOOK | a blind fresh AI agent re-implemented both scores and matched every main number (same model family) |

Open before anyone should rely on this: a human (ideally an author of one of the
original cohorts) reading the label table and the code.

### Labels

| cohort | platform | bacterial | viral | excluded | check |
|---|---|---|---|---|---|
| GSE63990 | GPL571 | 70 | 115 | 95 | Tsalik 2016 Sci Transl Med, Table 1 (AI-verified); excluded = 88 non-infectious + 7 re-run arrays |
| GSE42026 | GPL6947 | 18 | 41 | 33 | authors' GEO design text (Herberg 2013), Tsalik 2016 Table 1 |
| GSE60244 | GPL10558 | 22 | 71 | 65 | Suarez 2015 J Infect Dis, Methods |
| GSE40396 | GPL10558 | 8 | 22 | 35 | Hu 2013 PNAS, Results (febrile children only; checked by the author 2026-10-23) |
| GSE6269 | GPL96 | 73 | 18 | 6 | Ramilo 2007 Blood, Table 1 footnote (checked by the author 2026-10-23) |

The pipeline refuses to train if any parse disagrees with its published counts.

**GSE63990 re-run arrays (found 2026-10-21).** GEO holds 280 arrays of 273
patients: two batches share 7 re-run samples. The last 7 arrays (sample_274-280,
GSM1562133-GSM1562139) carry exactly the surplus over the published counts
(3 bacterial, 2 viral, 2 non-infectious), so they are excluded as the re-runs
(`config/cohorts.json` `exclude_samples`; a test fails if they come back). Which
arrays are the re-runs is **inferred, not proven**: GEO has no patient ids, and a
cross-batch correlation check supports the choice without pinning it.

### Prior art (checked 2026-10-10, deepened the same day; bvnGPS row corrected 2026-10-30)

This field is mature and commercial. Know this before telling anyone the work is novel.

| What exists | Status |
|---|---|
| **Inflammatix TriVerity** — 29-gene ML host-response test, bacterial + viral + severity scores | **FDA-cleared Jan 2025**, 1,222 patients across 22 sites, CMS add-on payment for FY2026 |
| **MeMed BV** — 3-protein host-response test (TRAIL, IP-10, CRP) | FDA-cleared |
| **Bodkin et al. 2022**, *Genome Medicine* — the big systematic comparison | **28 signatures, 51 datasets, 4,589 subjects.** See the note below: it deliberately does not measure cross-cohort transfer. Code ships only as a supplementary Python file, no repository, no stated licence |
| **Chawla et al. 2022**, *Cell Systems* — benchmarking framework | 30 signatures against a 17,105-profile compendium. Measures **robustness and cross-reactivity**, including signatures firing falsely on ageing. Finds the two are conflicting objectives |
| **MetaIntegrator** (Khatri lab, Stanford) | R package with leave-one-dataset-out built in, so the *method* is packaged. Reportedly removed from CRAN on 2024-10-14 |
| **bvnGPS** (`github.com/Ritchiegit/bvnGPS`) | Gene-pair signatures over 2,680 samples from 16 cohorts, including the five used here, with **both** a held-out-cohort and a random-split mode; understated in the 2026-10-10 check. Not the same measurement: two invocations plus a separate scoring script, one split rather than CV, probably not runnable as shipped |
| `rmaghembe1/host-response-module-transportability` | Nearest-sounding public repo. Scores **fixed, locked modules** across cohorts; its README says it "is not a diagnostic classifier-discovery study". No held-out-cohort training scheme, 0 stars, runnability unverified |

**TriVerity validation (Liesenfeld et al., *Nature Medicine* 2025, [doi:10.1038/s41591-025-03933-y](https://doi.org/10.1038/s41591-025-03933-y)).** SEPSIS-SHIELD, a prospective non-interventional study, enrolled 1,441 adults with suspected infection in 22 emergency departments (21 in the United States, 1 in Europe); 1,222 had valid TriVerity results and 729 of those had a consensus-adjudicated infection status, the main analysis population. In that population the Bacterial score reached an AUROC of 0.83 and the Viral score 0.91; the Severity score (mechanical ventilation, vasopressors or renal replacement within 7 days; 1,120 evaluable) reached 0.78, and with forced adjudication of all 1,222 patients the bacterial and viral AUROCs fall to 0.76 and 0.83. That is a prospective, multi-site validation of one locked, FDA-cleared assay against clinical adjudication, whereas this study is a retrospective reproduction on public cohorts comparing random cross-validation with leave-one-cohort-out; the numbers are not comparable and we make no claim that they are.

#### The one paragraph that justifies this repo

Bodkin et al. is the most thorough comparison in the field, and it answers a
different question than the one asked here. They did not apply one model across
cohorts. They **fitted a fresh lasso logistic regression inside each dataset
separately**, and estimated performance with nested leave-one-out
cross-validation *within* that dataset. Their stated reason, quoted from the
paper: *"Creating dataset-specific models overcomes batch effects since each
signature is optimized in each dataset."*

That is a sound choice for their question, which was "how good is this gene list
when it is tuned locally". It also means the batch effect was removed by
construction rather than measured, and the paper does not report a
discovery-versus-external drop at all.

So the gap this repo measures is not reported by the largest comparison in the
field, and the second-largest benchmark measures robustness and cross-reactivity
instead. Cross-cohort validation does happen in primary papers from
well-resourced groups, for example the 29-mRNA classifier trained across 18
studies, but their pipelines are not public.

**Therefore NOT the contribution here:** a new signature, a better AUC, a
clinical product, or more data than anyone else. Five cohorts and 458 samples is
roughly a tenth of Bodkin's scale.

**The contribution that survives scrutiny:** a small, public, runnable
implementation where one command computes the optimistic number and the
cohort-held-out number from the same data, and the label parsing checks itself
against published counts. An engineering and reproducibility contribution, not a
scientific-novelty one. Say it that way.

## License

MIT, see [LICENSE](LICENSE).
