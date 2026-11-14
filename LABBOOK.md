# LABBOOK — host-response-lab

> **Constitution:** [CORE.md](CORE.md) — the question, what success looks like, the kill criterion, the core ideas, and what a driver agent may propose. Diary: `LABBOOK.md`.

The story of this repo across time. One dated entry per working session, appended
before the session ends. Question drift gets written down the day it happens.

> **About this copy.** This is the lab diary from the private working repo, kept
> as written apart from a privacy scrub: internal tool, session and task-card
> names were removed, and one entry about unsent emails to other researchers
> was dropped. Commit hashes and branch names below refer to that private
> history and do not exist in this repository. Planning briefs under `docs/`
> and the HTML pages under `reports/` that some entries mention were internal
> working notes and are not included; the one exception is
> `reports/2026-11-08-writeup.html`, the write-up for a general reader. Entries "written by Claude" were written
> by an AI coding assistant; "the author" decided what is noted as decided.

---

## 2026-10-10 — repo built, first measurement taken

**Question asked.** Published host-response signatures that separate bacterial
from viral infection report AUCs around 0.93 to 0.95. How much of that survives
being tested on a hospital cohort the model has never seen?

**Why this and not a new classifier.** A prior-art sweep was run *before* any code
(see README). The field is commercially mature: Inflammatix TriVerity is FDA-cleared
since January 2025 (29 genes, bacterial + viral + severity scores, 1,222 patients
across 22 sites), MeMed BV is cleared as a 3-protein test, and Bodkin et al. 2022 in
Genome Medicine already compared 28 published signatures across 51 datasets. So the
question was deliberately narrowed from "build a better signature" to "ship runnable
code that measures the optimistic number and the cohort-held-out number in one run".
A deeper sweep the same evening sharpened this further; see the second entry below.

**What was built.** Fetch from the GEO FTP mirror, SOFT-format parsing, label
assignment with published-count verification, probe-to-symbol harmonization across
platforms, and two validation schemes. 24 tests, each constructing the bug it guards
against. CPU only; the GPU was never used.

**Data actually used.** Five public cohorts, 476 labelled samples, 12,040 gene
symbols shared across all five platforms.

| Cohort | Platform | Bacterial | Viral | Labels |
|---|---|---|---|---|
| GSE63990 | GPL571 | 73 | 117 | **verified against published counts** |
| GSE42026 | GPL6947 | 18 | 41 | **verified against published counts** |
| GSE60244 | GPL10558 | 22 | 71 | unverified |
| GSE40396 | GPL10558 | 8 | 35 | unverified |
| GSE6269 | GPL96 | 73 | 18 | unverified |

The two cohorts with published counts to check against reproduced them **exactly**,
which is the strongest evidence available that the label parsing is not inventing
its own study.

**Measured — raw expression scale.**

| Model | Random CV | Leave-one-cohort-out | Gap |
|---|---|---|---|
| Logistic L2 | 0.954 | 0.885 | **0.069** |
| LightGBM | 0.956 | 0.825 | **0.131** |

**Measured — with within-cohort rank normalization (`--rank-normalize`).**

| Model | Random CV | Leave-one-cohort-out | Gap |
|---|---|---|---|
| Logistic L2 | 0.743 | 0.893 | **−0.150** |
| LightGBM | 0.904 | 0.895 | **0.008** |

**What we now believe.**

1. **The gap is real and reproducible.** Random cross-validation on raw data gives
   0.95, matching the published within-cohort range. Holding out a whole cohort
   costs 0.07 to 0.13 AUC. The one published study that used fully independent
   prospective cohorts reported 0.93 falling to 0.84, so our numbers sit in the
   same place as the literature's own external validation.
2. **The more flexible model pays more.** LightGBM lost nearly twice what the
   regularized linear model lost, and on its worst held-out cohort (GSE6269) it
   fell to **0.625**, close to useless, where logistic regression held **0.944** on
   the same cohort. The booster found cohort-specific structure and leaned on it.
   This is the clearest single result of the day.
3. **Rank normalization appears to erase the gap, and that reading is NOT yet
   earned.** See below.

**Claim-check on finding 3.**
- *Claim one might want to make:* "within-cohort rank normalization fixes
  cross-cohort generalization."
- *What was actually measured:* with rank normalization the LightGBM gap fell to
  0.008 and its worst cohort rose from 0.625 to 0.922.
- *Rival explanation, and the better one:* rank normalization computes each
  patient's percentile **using the other patients in the same cohort**. For a
  held-out cohort that means the transform reads the whole test set before
  predicting. That is transduction, not generalization, and a bedside test seeing
  one patient at a time cannot do it. The apparent fix may be the leak.
- *The control that would settle it:* normalize each sample against itself only
  (within-sample ranking across genes, no other patient involved) and re-measure.
  Until that is run, finding 3 stays **unproven**.
- *Unexplained and flagged:* logistic regression's random-CV AUC **fell** to 0.743
  under rank normalization while its held-out-cohort AUC rose, giving a negative
  gap. A negative gap should not happen if both schemes measure the same thing.
  The leading suspicion is that per-cohort ranking is confounded with each cohort's
  class balance, which ranges from 19% bacterial (GSE40396) to 80% (GSE6269). Not
  yet tested. **Do not quote the rank-normalized numbers to anyone until this is
  understood.**

**Where it lives.** `results/summary_raw.csv`, `results/summary_ranknorm.csv`,
`results/folds_*.csv`, `results/label_review.md`, `results/run_*.json`.

**Two bugs found by running it, both now pinned by tests.** GEO's directory naming
collapses three-digit accessions to the bare `nnn` form, so `GPL571` lives in
`GPLnnn` and the first implementation generated 404s. And `.annot.gz` files are
SOFT format whose real header follows a `!platform_table_begin` fence, so reading
the first non-comment line silently skipped every cohort.

**Next session, in order.**
1. Run the within-sample normalization control. It decides whether finding 3 is a
   result or an artefact, and nothing else should be reported until it is done.
2. Explain or eliminate the negative logistic gap.
3. Find published class counts for GSE60244, GSE40396 and GSE6269 so three more
   cohorts move from unverified to verified.
4. Only then: per-cohort error analysis, and whether the gap tracks platform,
   patient age, or class balance.

---

## 2026-10-10 (later) — deeper prior-art sweep, claim narrowed again

**Question asked.** The author asked directly whether a project like this already
exists. The morning sweep was not deep enough to answer that, so it was redone
against the primary literature and against public code.

**What was found.**

| Source | What it actually does |
|---|---|
| Bodkin et al. 2022, *Genome Medicine* | 28 signatures, 51 datasets, 4,589 subjects. **Fitted a fresh lasso logistic regression inside each dataset** and used nested leave-one-out CV *within* each dataset. Quoted: "Creating dataset-specific models overcomes batch effects since each signature is optimized in each dataset." Does **not** report a discovery-versus-external drop. Code is a supplementary Python file only, no repo, no licence stated |
| Chawla et al. 2022, *Cell Systems* | 30 signatures against a 17,105-profile compendium. Measures **robustness and cross-reactivity**, including false firing on ageing; finds the two conflict |
| MetaIntegrator (Khatri lab) | R package with leave-one-dataset-out built in. Reportedly off CRAN since 2024-10-14 |
| `rmaghembe1/host-response-module-transportability` | Scores **fixed, locked modules** across cohorts. README states it "is not a diagnostic classifier-discovery study". No held-out-cohort training scheme, 0 stars, runnability unverified |
| bvnGPS | Gene-pair signatures, 2,680 samples, 16 cohorts |

**What we now believe.**

1. **The question is not taken.** The largest comparison in the field removed the
   batch effect by construction instead of measuring it, and says so plainly. The
   second-largest benchmark measures a different property. So "train on four
   hospitals, test on the fifth, report the drop" is not what either did.
2. **The method is not novel.** Leave-one-dataset-out is packaged in
   MetaIntegrator and used routinely in multi-cohort meta-analysis. Nothing here
   invents a technique.
3. **Cross-cohort validation does exist in primary papers** from well-resourced
   groups, for example the 29-mRNA classifier trained across 18 studies. Their
   code is not public. So the shortage is runnable code, not insight.
4. **Scale is the honest weakness.** Five cohorts and 476 samples is about a tenth
   of Bodkin's 51 datasets and 4,589 subjects. Any claim must be framed as
   reproducibility and engineering, never as a bigger or better study.

**Correction made.** The morning's README, this labbook and the brain note all
described Bodkin as validating "with leave-one-out validation", which reads as
leave-one-dataset-out. It was leave-one-out *within* each dataset, which is close
to the opposite of what matters here. All copies corrected.

**Where it lives.** README section *Prior art*, including the paragraph that now
justifies the repo by quoting Bodkin's own design choice.

**Next session.** Unchanged, and the control experiment is still first. Add one
item: pull the numbers from Bodkin's Table 2 so this repo's random-CV result has a
published figure to sit beside.

---

## 2026-10-18 — within-sample normalization control

**Question asked.** The 2026-10-10 rank-normalised run seemed to erase the
held-out-cohort gap. But that transform ranks each gene against the OTHER patients
of the same cohort, so on a held-out cohort it reads the whole test set before
predicting. Does the gap survive a transform that reads one patient at a time,
which is the only kind a bedside test could use?

**What was done.** Added `rank_normalize_within_sample` (`src/splits.py`): each
sample's own 12,040 gene values are ranked against each other (percentile, ties
averaged, NaN -> 0.5). No other patient is read, in train or test. Added flag
`--within-sample-normalize` to `run_pipeline.py` (exclusive with
`--rank-normalize`). The model's `StandardScaler` sits inside the pipeline and is
fitted on the training fold only (read in `src/train.py`), so it adds no leak.
New test `tests/test_within_sample_normalize.py` swaps every other row of the
test cohort for noise and asserts the row's output is identical; a second test
builds the leaky cohort-wise transform and asserts the same guard FIRES on it.
Suite: 28 passed. `config/cohorts.json` untouched (no diff).
Command: `python run_pipeline.py --train --within-sample-normalize` (CPU, ~2 min).

**Measured** (mean AUC over 5 folds, one seed, 476 samples, 5 cohorts, 12,040 genes):

| model | transform | random-CV | leave-one-cohort-out | gap |
|---|---|---|---|---|
| LightGBM | raw | 0.956 | 0.825 | **0.131** |
| LightGBM | within-sample | 0.954 | 0.857 | **0.097** |
| logistic L2 | raw | 0.954 | 0.885 | **0.069** |
| logistic L2 | within-sample | 0.961 | 0.900 | **0.061** |

(The cohort-wise rank-normalised run stays in `results/summary_ranknorm.csv`; its
numbers are not quoted outside this repo.)

Held-out AUC per cohort, raw -> within-sample:

| held-out cohort | LightGBM | logistic L2 |
|---|---|---|
| GSE40396 (n=43, 8 bacterial) | 0.868 -> 0.761 | 0.889 -> 0.893 |
| GSE42026 | 0.885 -> 0.818 | 0.812 -> 0.860 |
| GSE60244 | 0.848 -> 0.894 | 0.848 -> 0.866 |
| GSE6269 | 0.625 -> 0.916 | 0.944 -> 0.953 |
| GSE63990 | 0.899 -> 0.895 | 0.930 -> 0.929 |

**What it means in plain words.** Random-CV still sits around 0.95-0.96 under every
transform. Held-out-cohort AUC stays lower. A leak-free transform leaves most of
the gap in place: LightGBM 0.131 -> 0.097, logistic 0.069 -> 0.061. The
LightGBM shrinkage is mostly ONE cohort (GSE6269, 0.625 -> 0.916) while two
others got worse; with 5 folds, one seed and no confidence interval, that is not
evidence of a general improvement. GSE40396 has only 8 positives, so its AUC is
very noisy.

**What we now believe.**
1. **Finding 3 ("rank normalization fixes it") is RETRACTED.** The closure seen
   on 2026-10-10 does not reproduce under a bedside-legal transform, which fits
   the transduction explanation (the cohort-wise transform read the test cohort).
   This run does not by itself prove transduction was the whole cause.
2. The held-out-cohort drop is real under one-patient-at-a-time scoring for both
   models, at 0.06-0.10 AUC. It is not shown to be "cohort-specific structure
   removable by a per-cohort transform".
3. Within-sample ranking keeps only gene ORDER, not absolute level, so this
   control tests one particular bedside-legal transform, not every possible one.
4. Still OPEN (next-step 2): the class-balance rival and the negative logistic gap
   seen in the cohort-wise run. Under the within-sample transform the logistic gap
   is positive (0.061), so the negative gap looks specific to the cohort-wise
   transform, but that is not yet explained.

**Where it lives.** `results/summary_withinsample.csv`, `folds_withinsample.csv`,
`run_withinsample.json`; code in `src/splits.py` + `run_pipeline.py`; test in
`tests/test_within_sample_normalize.py`; branch `within-sample-control`.

**Next session.** Item 2 above (class balance / negative logistic gap), then more
seeds for a spread on the LightGBM gap before quoting any number.

## 2026-10-18 (evening) — CORE.md accepted

**Question.** Does the lab adopt the drafted constitution (branch `docs/core-md`)?
**Done.** The author accepted it as proposed; merged into master with a merge commit (38920b1). The cohort set is FROZEN at the current five GEO cohorts.
**Result (measured).** Nothing run; docs only.
**Belief.** Unchanged. CORE §2 success condition 1 (within-sample control) is already met by the 2026-10-18 within-sample control; conditions 2 (labels for GSE60244, GSE40396, GSE6269 verified against papers) and 4 (one traceable sentence) remain.
**Where.** `CORE.md`; this file.

---

## 2026-10-19 — PLAN, pre-registered before any run: why is the logistic gap negative? (CORE §2 condition 3)

_Night run on branch `night/2026-10-19`, written by Claude while the author sleeps (they
approved the night run at 02:50). This plan is committed BEFORE any experiment
runs and is not edited afterwards. Any change gets a new dated line under it._

**The puzzle.** On 2026-10-10, with cohort-wise rank normalisation
(`--rank-normalize`), logistic regression scored 0.743 in random cross-validation
but 0.893 on held-out cohorts, a gap of −0.150 (`results/summary_ranknorm.csv`).
Holding a whole hospital out should not be easier than a random split, so
something in the measurement is off. The per-fold file points the same way: all
five random-CV folds are low (0.687 to 0.770, `results/folds_ranknorm.csv`), so
it is not one unlucky fold.

**Claim to test (one sentence).** The negative gap is a measurement artefact:
the random-CV AUC is computed over test patients from all five cohorts mixed
together, and under cohort-wise ranking each cohort's scores sit on their own
scale, shifted by that cohort's share of bacterial patients (19 % to 80 %); the
held-out-cohort AUC is computed inside one cohort, where such a shift cancels.

In plain words: AUC asks "pick one bacterial and one viral patient: how often
does the bacterial one get the higher score?". In random CV about three in four
of those pairs mix two hospitals. Ranking inside a mostly-bacterial hospital
pushes all its patients' scores down (and inside a mostly-viral one, up). Then
the mixed-hospital pairs go wrong, even if the model sorts the patients inside
each hospital well.

**What the instrument can see.** The same test patients and the same scores,
sorted into pairs from the SAME cohort and pairs from DIFFERENT cohorts; and
re-runs where every cohort's bacterial share is made equal.

**Experiments.** CPU only. One new script `run_gap_diagnostics.py`, helpers in
`src/` with tests that go red on the broken version.
- **E0 reproduction.** Re-run both schemes (seed 0, same folds) for raw,
  cohort-wise and within-sample transforms, both models, keeping every test
  score. Passes if every fold AUC matches `results/folds_*.csv` within 0.005
  (logistic) or 0.02 (LightGBM; thread count can move it a little). If not:
  report the drift and use the fresh numbers.
- **E1 pair split.** Per random-CV fold: AUC over same-cohort pairs ("within"),
  over different-cohort pairs ("cross"), over all pairs ("pooled", the number on
  file). Pooled must equal the pair-weighted mix of the two (built-in check).
  Also, from the out-of-fold scores: random-CV AUC inside each cohort next to
  that cohort's held-out AUC; mean score per cohort and class; whether any
  logistic fit stopped at its 3,000-iteration cap.
- **E2a class balance.** Drop patients at random until every cohort is 50 %
  bacterial (8+8, 18+18, 22+22, 18+18, 73+73 = 278 patients), redo the
  cohort-wise ranking on that subset, run both schemes. 10 random subsets; the
  same subsets also under raw and within-sample.
- **E2b size control.** Same cohort sizes as E2a, but each cohort keeps its
  ORIGINAL bacterial share. 10 random subsets. This separates "equal balance"
  from "fewer patients".
- **E3 transduction, random-CV side.** Cohort-wise ranking where each test
  patient is ranked only against the TRAINING patients of its own cohort (no
  other test patient is read). Compare with the version on file.
- **E4 transduction, held-out side.** For each held-out cohort, 20 random splits
  into half A and half B. Score B two ways: ranked within the whole cohort (as
  on file), or ranked against half A only (no other B patient is read). Compare
  the AUC on B.

**Pre-registered predictions** (logistic regression, cohort-wise ranking,
unless named otherwise):
- **P1 (E1):** within-cohort-pair AUC ≥ 0.88 (mean of 5 folds). The claim is
  KILLED if it is < 0.85: then the model really is worse in random CV, not just
  the scoring. (Cross-cohort-pair AUC is reported; with pooled at 0.743 and about
  3/4 of pairs being cross-cohort, it must then be about 0.70 or lower.)
- **P1b (E1):** random-CV AUC inside each cohort minus that cohort's held-out
  AUC, averaged over the 5 cohorts, ≥ −0.02. KILLED if ≤ −0.05.
- **P2 (E2a):** with equal class balance, the mean gap over 10 subsets is
  ≥ −0.02. The class-balance part is KILLED if the mean gap is ≤ −0.05.
- **P3 (E2b):** with the original balance at the same smaller sizes, the mean
  gap stays ≤ −0.05. If it is ≥ −0.02 instead, the smaller sample (not the
  balance) closed the gap, and the class-balance part is NOT supported.
- **P4 (E3):** with test patients kept out of the ranking, the random-CV pooled
  AUC stays within ±0.03 of the E0 number (0.743 on file). If it rises to
  ≥ 0.85, transduction on the random-CV side causes the low number and the claim
  is wrong.
- **P5 (E4):** the mean AUC difference on half B (whole-cohort ranking minus
  half-A ranking) is ≤ +0.02. If larger, reading the other test patients lifts
  the held-out number.
- A result between the two numbers of a prediction is reported as "partial".
- LightGBM: direction only (within-cohort-pair AUC ≥ pooled AUC). Raw and
  within-sample transforms: measured and reported, no prediction made.

**Rival explanations and the control for each.**
- Metric mismatch (the claim itself): E1.
- Class balance as the reason for the per-cohort shift: E2a; "it is only fewer
  patients": E2b.
- Transduction (the transform reads test patients): E3 for random CV, E4 for
  held-out cohorts.
- A pipeline bug (rows and labels out of order; a fit that did not finish):
  E0, an index-alignment assert, the iteration-cap check. LightGBM on the same
  matrix already reaches 0.904 in random CV, so a gross row/label mix-up is
  unlikely.
- Other cohort differences (age, platform, severity) that shift scores: what is
  left over if P2 fails; they become the next suspects.

**Who judges.** These thresholds, committed to git before the run (the commit
time is the proof). No blind outside agent is available in this night lane, so
every result is "self-assessed against pre-registered thresholds" until a fresh
agent or the author re-runs `python run_gap_diagnostics.py` from `data/`.

**Budget and fences.** CPU only, about 1,300 model fits, expected 20 to 60
minutes on the desktop. With all 16 threads this can add up to more than one
CPU-hour, so this entry is the plan entry CORE §7 asks for. No change to
`config/cohorts.json`, `expected_counts`, the cohort set, CORE §1–§5 or the
existing results files. Rank-normalised numbers stay inside this repo.

_2026-10-19 03:15, clarification before the real-data run (no threshold changed)._
Averaging, fixed in code before any real-data run: the P2/P3 gap of one subset =
mean random-CV AUC over its 5 folds minus mean held-out AUC over its 5 cohorts,
then the mean over the 10 subsets. P5 = mean over the 20 splits inside each
cohort, then the mean over the 5 cohorts (each cohort counts once, like the
held-out AUC). P1b uses the out-of-fold scores of the 5 random-CV folds, inside
each cohort. Seeds: subsets 0-9, half splits 0-19, random-CV folds as on file
(seed 0). One thing was already seen before the run, on made-up data only: the
positive-control test in `tests/test_gap_diagnostics.py` (4 synthetic cohorts,
one shared signal, NO batch effect) gives pooled random-CV 0.80 vs held-out 0.99
when the bacterial shares are 20/20/80/80 %, and 0.998 vs 0.996 when all are
50 %. So the mechanism CAN produce a negative gap; whether it is what happens in
the real cohorts is what E1-E4 test.

---

## 2026-10-19 — RESULT: a post-hoc explanation for the negative logistic gap (CORE §2 condition 3)

_Wording corrected after the outside check (04:33); every changed sentence is
listed in the "Outside check" entry at the end of this file._

**In short.** The −0.150 was not about held-out cohorts at all. Post-hoc
explanation, consistent with E3 and the exploratory mirror table, and NOT the
cause my PLAN pre-registered (P4 killed that one): the cohort-wise
ranking on file was done on each whole hospital BEFORE random CV split it into
training and test patients. That tie, together with hospitals that are 19 % to
80 % bacterial, pushed each hospital's test patients onto a shifted score scale,
so pairs of patients from two different hospitals were often sorted the wrong
way. Ranking test patients against training patients only, or making every
hospital 50 % bacterial, removes the negative gap.

**Question asked.** Under cohort-wise ranking (`--rank-normalize`), why did
logistic regression score 0.743 in random CV but 0.893 on held-out cohorts, a
gap of −0.150? Pass and kill rules: the PLAN entry above, committed at 03:07
(cc7ba64) before the first run started at 03:15.

**What ran.** `python run_gap_diagnostics.py` (parts E0-E4, about 33 minutes),
then `--parts mirror` (1 minute). Desktop PC, CPU only, branch
`night/2026-10-19`. 476 patients, 5 cohorts, 12,040 genes, the same models and
the same random-CV folds (seed 0) as the files on record. "Score" below means
the model's log-odds for "bacterial": 0 is 50:50, +2 is about 7 to 1 for
bacterial, −2 about 7 to 1 for viral.

**E0, reproduction: passed.** Every fold AUC on file came out again. The
largest difference was 0.0008 (one raw logistic fold; the limit was 0.005).
No logistic fit in any part hit its iteration cap (at most 153 of 3,000
steps), and rows and labels were checked to be in the same order. So the
−0.150 is not a crash, an unfinished fit or a row mix-up.

**E1, where the gap lives.** AUC = "pick one bacterial and one viral patient:
how often does the bacterial one get the higher score". In a random-CV test fold,
23 % of these pairs come from the same hospital and 77 % from two different
hospitals. Random-CV AUC split by pair type (mean of 5 folds), next to the
held-out AUC (mean of 5 cohorts, always inside one hospital):

| transform | model | random-CV, all pairs (on file) | same-hospital pairs | cross-hospital pairs | held-out cohort | gap |
|---|---|---|---|---|---|---|
| raw | logistic | 0.954 | 0.954 | 0.953 | 0.884 | 0.069 |
| raw | LightGBM | 0.956 | 0.949 | 0.958 | 0.825 | 0.131 |
| within-sample | logistic | 0.961 | 0.959 | 0.962 | 0.900 | 0.061 |
| within-sample | LightGBM | 0.954 | 0.951 | 0.954 | 0.857 | 0.097 |
| cohort-wise, as on file | logistic | **0.743** | 0.954 | **0.680** | 0.893 | **−0.150** |
| cohort-wise, as on file | LightGBM | 0.904 | 0.951 | 0.889 | 0.895 | 0.008 |

In plain words: under cohort-wise ranking, logistic regression sorts patients
INSIDE one hospital as well as on raw data (0.954), but for pairs from two
hospitals it gets only 0.680 right. Those mixed pairs are 3 in 4 of all pairs,
so they drag the random-CV number down to 0.743. P1 (same-hospital pairs
≥ 0.88): **supported**, 0.954. P1b (random-CV AUC inside each cohort minus its
held-out AUC, mean of 5 cohorts, ≥ −0.02): **supported**, +0.061.
For raw and within-sample data the two kinds of pairs agree to within 0.01
(raw logistic 0.954 vs 0.953), so the positive gaps on file do NOT carry this
artefact.

**E3, what breaks the mixed pairs.** The cohort-wise ranking on file ranks each
WHOLE cohort before the random split, so every test patient is ranked together
with the training patients of its hospital. E3 ranks each test patient against
the training patients of its own cohort only:

| model | random-CV, all pairs | same-hospital pairs | cross-hospital pairs | held-out (unchanged) | gap |
|---|---|---|---|---|---|
| logistic | 0.915 (was 0.743) | 0.955 | 0.903 (was 0.680) | 0.893 | +0.022 |
| LightGBM | 0.933 (was 0.904) | 0.949 | 0.928 (was 0.889) | 0.895 | +0.037 |

The PLAN said: if this rises to 0.85 or more, transduction on the random-CV side
causes the low number and the claim is wrong. It rose to 0.915, so **P4 KILLED
the cause named in my PLAN claim** (simple ranking arithmetic of class balance).
The "mixed pairs" part of the claim stands (E1); the post-hoc explanation (not
pre-registered) is the ranking of test and training patients together.

**Mirror (exploratory: added after seeing E1 and E3, no threshold was set for
it).** Why does ranking them together hurt? Cohort-wise ranks give each gene
almost the same average in every hospital (a percentile rank averages
(n+1)/2n, which shifts a little with hospital size). For a model that is a
weighted sum, like logistic regression, this almost fixes the average score of
each hospital's training and test patients taken together. Mean logistic score of the SAME fold model,
averaged over the 5 folds (`results/gapdiag_mirror.csv`):

| cohort | bacterial share | training patients | test patients | both together | test patients, E3 ranking |
|---|---|---|---|---|---|
| GSE40396 | 19 % | −1.94 | +1.75 | −1.18 | −1.39 |
| GSE60244 | 24 % | −1.78 | +1.50 | −1.14 | −1.16 |
| GSE42026 | 31 % | −1.53 | +0.01 | −1.16 | −1.83 |
| GSE63990 | 38 % | −1.16 | −0.95 | −1.12 | −1.18 |
| GSE6269 | 80 % | +0.20 | −6.23 | −1.14 | −0.77 |

The model fits each hospital's training patients toward that hospital's
bacterial share: their average score rises from −1.94 (19 % bacterial) to +0.20
(80 %). Because the total is almost fixed (−1.12 to −1.18 in every hospital), the same
hospital's test patients move the other way (+1.75 down to −6.23). Result: the bacterial test patients of GSE6269 average −3.96,
below the viral test patients of GSE40396 (+0.64), so mixed pairs fail. With the
E3 ranking the test column no longer slides (−0.77 to −1.83), and those two
groups sit at +1.88 and −2.68 again. Leave-one-cohort-out never splits a
hospital into training and test patients, so it never meets this effect.

**E2, is unequal class balance needed?** 10 random subsets per design; the
cohort-wise ranking is redone on each subset, still before the split as on file:

| subsets | model | random-CV | held-out | gap: mean (lowest to highest) |
|---|---|---|---|---|
| E2a: every cohort 50 % bacterial (278 patients) | logistic | 0.955 | 0.889 | **+0.066** (+0.029 to +0.120) |
| E2a | LightGBM | 0.931 | 0.880 | +0.051 (+0.036 to +0.071) |
| E2b: same sizes, each cohort's own share | logistic | 0.752 | 0.870 | **−0.118** (−0.144 to −0.076) |
| E2b | LightGBM | 0.860 | 0.867 | −0.007 (−0.032 to +0.035) |

With equal shares the negative gap is gone in all 10 subsets (P2:
**supported**). With the original shares at exactly the same sizes it is back in
all 10 (P3: **supported**). So it is the unequal balance, not the smaller sample.
On the same subsets, raw and within-sample gaps stay positive (logistic +0.04 to
+0.06, LightGBM +0.08 to +0.09; `results/gapdiag_balance_summary.csv`).

**E4, the held-out side.** For each held-out cohort, half B was scored with
ranks from the whole cohort (as on file) or with ranks from the other half A
only (no other B patient read), 20 random halves per cohort. Mean AUC difference
on B: logistic +0.000 (per cohort −0.003 to +0.003), LightGBM −0.002. P5
(≤ +0.02): **supported**. Reading the other test patients does not lift the
held-out number. It still needs other patients of the SAME hospital as a
reference, so it is still not a one-patient-at-a-time (bedside) number.

**What we now believe.**
1. CORE §2 condition 3 is proposed as met, pending the author's decision
   (self-assessed, see below).
   The −0.150 says nothing about cohort generalisation. Post-hoc explanation
   (not pre-registered): it is an artefact of
   random CV meeting a ranking done before the split: unequal bacterial shares
   plus training and test patients tied together.
2. My PLAN claim was half right. "The mixed-hospital pairs fail" held. The cause
   I named was killed by its own pre-registered control (P4); the post-hoc
   explanation, consistent with E3 and the exploratory mirror table, is
   a transduction step, the kind of mistake CORE §4 warns about, here making a
   number look WORSE instead of better.
3. The positive raw and within-sample gaps on file are not touched by this
   artefact: their same-hospital and cross-hospital pairs agree to within 0.01.
4. The 2026-10-10 "LightGBM gap closes to 0.008" was partly the same artefact:
   its random-CV side was pulled down (0.904 instead of 0.933 with E3 ranking).
   This supports the 18.10. retraction of finding 3; it does not reverse it.

**Claim-check (2026-10-19, on belief 1).**
- *Claim (post-hoc, not the pre-registered claim):* the −0.150 is a random-CV artefact: ranking each whole cohort before
  the split ties a cohort's training and test patients together, and with shares
  from 19 % to 80 % bacterial this puts each hospital's test patients on a
  shifted scale, so the cross-hospital pairs fail.
- *Measured:* same-hospital pairs 0.954 vs cross-hospital 0.680 (E1); test
  patients ranked against training patients only: 0.915 (E3); equal shares
  +0.066, own shares at the same sizes −0.118 (E2); held-out side unchanged
  (+0.000, E4); almost fixed train + test average of −1.12 to −1.18 (mirror,
  exploratory).
- *Rivals and their controls:* held-out number inflated by transduction, killed
  by E4. "Only fewer patients", killed by E2b. A bug or unfinished fit, killed
  by E0 and the iteration counts. Ranking arithmetic alone (my PLAN claim),
  killed by E3. The main effect has no rival left. A smaller second effect
  remains: under training-only ranking (E3) cross-hospital pairs still score
  below same-hospital pairs (0.903 vs 0.955), probably the class-balance ranking
  arithmetic; it is not tested separately. The step-by-step wording (the almost
  fixed total) rests on a table I added after seeing data.
- *Judge:* thresholds committed before the run (cc7ba64). No outside judge yet,
  so this is **self-assessed**. It fails if a blind re-run of
  `python run_gap_diagnostics.py` from `data/` gives a different verdict for any
  P1-P5.

**Stays unverified.**
- No outside judge yet (above).
- One random-CV fold seed (0), as in the files on record. E2 varies the patients
  (10 subsets) but not the fold seed.
- Labels of GSE60244, GSE40396 and GSE6269 are still unverified (condition 2),
  and every number here inherits that.
- Why LightGBM is hit less: inferred. A tree model is not a weighted sum, so the
  total is not fixed (its train + test averages differ by hospital in the mirror
  file), but this was not tested further.
- Under E3 ranking a smaller cross-hospital penalty is left (0.903 vs 0.955):
  inside each hospital the class averages still slide with the bacterial share
  (the PLAN's ranking-arithmetic idea, a small effect). Not tested separately.

**Open for the author (nothing decided alone).**
- Whether condition 3 counts as met; then the CORE §6 line "the negative
  logistic gap is unexplained" can be updated (CORE.md was not touched).
- `run_pipeline.py --rank-normalize` still ranks before the split. Option A
  (recommended): leave it; the file stays a documented artefact and its numbers
  stay inside this repo. Option B: rank inside each random-CV fold (as in E3) and
  re-run.
- For condition 4: counting every hospital once (the way the held-out number
  does) gives smaller raw gaps: logistic 0.053 instead of 0.069, LightGBM 0.105
  instead of 0.131 (`results/gapdiag_cohorts.csv`). The gap stays positive
  either way; which number the one-sentence answer uses is a reporting choice.
- Merge `night/2026-10-19` into `master`.

**Where it lives.** Results: `results/gapdiag_folds.csv` (E0/E1),
`gapdiag_cohorts.csv`, `gapdiag_inductive.csv` (E3), `gapdiag_balance_folds.csv`
+ `gapdiag_balance_summary.csv` (E2), `gapdiag_splithalf.csv` (E4),
`gapdiag_mirror.csv`, `run_gapdiag.json` (verdicts, versions, run times).
Code: `src/gap_diagnostics.py`, `src/gap_experiments.py`,
`run_gap_diagnostics.py`. Tests: `tests/test_gap_diagnostics.py` (16 new, each
checked to go red on a broken version; suite 44 passed). Existing results files
and `run_pipeline.py` are unchanged.

**Next session.** CORE §7 order item 3: published class counts for GSE60244,
GSE40396 and GSE6269 (condition 2, a literature lookup). Later: more fold seeds
for a spread on every gap before quoting any number.

## 2026-10-19 (night run) - published class counts, candidates only

Two blind extractions of the published class counts for GSE60244, GSE40396 and GSE6269 are in `docs/published-counts-candidates.md`. They agree on every class count (GSE40396 and GSE6269 differ only in layout and in paper-total 144 vs GEO-total 143). Not yet in `config/cohorts.json`: The author checks against the papers first (CORE §2.2).

## 2026-10-19 — Outside check (04:33)

_Night run, branch `night/2026-10-19`. Written by Claude; the author is asleep and
decided nothing in this entry._

**Who checked.** A fresh agent that did not write the RESULT entry or the code.
It re-ran everything in a separate copy of the repo. **Verdict: "partly".**

**What it reproduced.** Identical numbers for E0-E4 and the mirror table. Its
own re-implementation of E1, E3 and the held-out AUC agrees with ours to 4
decimals. With 2 extra random-CV fold seeds the same pattern appears (our run
had used only seed 0).

**The 5 problems it found.**
1. The headline gave a cause (training and test patients tied together because
   each whole hospital is ranked before the split) that was NOT the
   pre-registered cause. The PLAN named ranking arithmetic, and the PLAN's own
   rule P4 killed it. The tied-together explanation is post-hoc, consistent
   with E3 and the exploratory mirror table.
2. "CORE §2 condition 3 is met" is the author's decision, not mine.
3. "No rival is left standing for the main effect" went too far: under
   training-only ranking (E3) a smaller cross-hospital penalty remains, and the
   checker's regularisation probe agrees (cross-hospital pair AUC rises as C
   drops from 0.05 to 0.00005 but never fully reaches the same-hospital level).
4. "Fixed average score" is only almost fixed: pandas `rank(pct=True)` has a
   mean of (n+1)/2n, which shifts a little with hospital size; the 5 hospitals'
   train + test averages agree to within about 0.03-0.09 log-odds.
5. Test gap: no test guarded E4's key column `auc_b_whole`. Replacing it with
   the half-A AUC makes P5 "supported" for free (difference 0), and all 16 new
   tests still passed.

**Sentences changed in the RESULT entry above** (nothing else in it was
touched; the old wording is in git, commit 8673d5e):
- Heading: "the negative logistic gap is explained" → "a post-hoc explanation
  for the negative logistic gap". A one-line note under the heading points here.
- "In short": added "Post-hoc explanation, consistent with E3 and the
  exploratory mirror table, and NOT the cause my PLAN pre-registered (P4 killed
  that one)". (Problem 1)
- E3 paragraph: "the cause is the ranking of test and training patients
  together" → "the post-hoc explanation (not pre-registered) is the ranking
  ...". (Problem 1)
- Mirror paragraph: "the same average in every hospital" → "almost the same
  average ... (a percentile rank averages (n+1)/2n, which shifts a little with
  hospital size)"; "this fixes the average score" → "this almost fixes the
  average score"; "the total is fixed (about −1.1 in every hospital)" → "the
  total is almost fixed (−1.12 to −1.18 in every hospital)". (Problem 4)
- Belief 1: "CORE §2 condition 3 is met by tests that were run" → "CORE §2
  condition 3 is proposed as met, pending the author's decision"; "It is an
  artefact ..." → "Post-hoc explanation (not pre-registered): it is an
  artefact ...". (Problems 1, 2)
- Belief 2: "the real cause is a transduction step" → "the post-hoc
  explanation, consistent with E3 and the exploratory mirror table, is a
  transduction step". (Problem 1)
- Claim-check: "*Claim:*" → "*Claim (post-hoc, not the pre-registered
  claim):*"; "fixed train + test average" → "almost fixed train + test
  average"; "No rival is left standing for the main effect, but the
  step-by-step wording (the fixed total) ..." → "The main effect has no rival
  left. A smaller second effect remains: under training-only ranking (E3)
  cross-hospital pairs still score below same-hospital pairs (0.903 vs 0.955),
  probably the class-balance ranking arithmetic; it is not tested separately.
  The step-by-step wording (the almost fixed total) ...". (Problems 1, 3, 4)
- Not changed, now out of date: "No outside judge yet" (claim-check *Judge*
  and "Stays unverified"). This entry is the outside judge; its verdict is
  "partly", so the RESULT stays a proposal for the author.
- Not changed (outside this file): the `score_mirror` docstring in
  `src/gap_experiments.py` still says "sum to a fixed total". That is exact
  only for equal-sized cohorts (as in the synthetic test); for the real
  cohorts it is almost fixed.

**Problem 5, fixed.** New test
`test_split_half_whole_column_is_the_whole_cohort_auc_on_half_b` in
`tests/test_gap_diagnostics.py` recomputes the whole-cohort AUC on half B on
its own and checks that the two E4 columns really differ on the synthetic data.
On a scratch copy with `auc_b_whole` replaced by the half-A AUC it fails
(1 failed, 16 passed); on the real code the full suite passes (45 passed).

**Labels (condition 2), one line.** Our own parse is now compared with the two
blind count extractions in `docs/published-counts-candidates.md`, section
"Comparison with our own parse": GSE60244 and GSE6269 (GPL96) match; GSE40396
does NOT (our 35 viral include 13 afebrile virus-positive children the paper
does not count as viral). The rule fix is the author's call; nothing was changed
or re-run.

**What we now believe.** The pair-split finding (E1) and the E2-E4 verdicts
stand and were reproduced by an outside agent. The "why" is a post-hoc
explanation, not a pre-registered result. Condition 3 is proposed as met,
pending the author. All gap numbers still inherit the GSE40396 labels, which are
probably wrong for 13 children.

_Coordinator note, 2026-10-19 04:48: the `score_mirror` docstring listed above
as "not changed" now says "almost fixed" too (commit below). The two "No
outside judge yet" lines in the RESULT entry stay as written, because this
entry already marks them out of date. A second fresh agent checked the three
fix-up commits (79282d9, 6789346, d0601cc): verdict "holds"; its own break-it
test turns only the new E4 test red; it re-counted GSE40396 from the GEO files
(22 febrile virus-positive, 13 afebrile virus-positive, 22 virus-negative,
8 febrile bacterial), so the label mismatch is confirmed._

## 2026-10-19 — PLAN, pre-registered before any run: the seesaw test (CORE §2 condition 3)

**Status.** Designed by two independent planner agents and merged by a critic agent (read-only, nothing ran). **Approved by the author 19.10.2026**, with one condition: the real run waits until their paper check (condition 2) is done and the GSE40396 label rule is fixed. Why a re-test at all: The author decided on 19.10. that the seesaw explanation in the RESULT entry above is post-hoc and must pass its own pre-registered test before condition 3 counts.

**Claim.** Under --rank-normalize, each hospital (cohort) is ranked as a whole before random CV splits it. For any model that adds up one gene at a time, this makes the hospital's total score almost fixed. That part is algebra. It is a MANIPULATION CHECK, not evidence. The testable part, the seesaw: the model holds a hospital's training patients close to their labels, so that hospital's few test patients must take up the whole mismatch. Their push is about k times larger (k = number of folds). That push mis-sorts cross-hospital test pairs, and it is the MAIN reason why random-CV AUC falls below leave-one-cohort-out AUC for logistic regression. The class-balance ranking arithmetic is the fuel, not a separate cause. The seesaw is the step that multiplies it. Supported only if all 3 predictions pass on BOTH cohort sets (all 5 cohorts, and the 4 without GSE40396). One kill in either set = not supported, and condition 3 stays open.

**Budget.** About 50 CPU-minutes, CPU only. Fold seeds 100-104 (the outside checker's extra seeds were 1 and 2, see its report of 19.10.; no clash).

### P1-fold-count-lever (about 25 CPU-minutes)

- **What it tests:** Concentration. A hospital's test patients carry about k times the push. So the EXTRA cross-hospital damage from whole-hospital ranking must grow with the fold count k. Plain class-balance arithmetic, label noise, hospital size and hospital biology do not depend on how the hospital is split, so they predict no change with k.
- **Why it is new:** Every earlier run (E0-E4, the mirror table, the checker's 2 seeds, the checker's C probe) used k=5 only. Nobody has changed the train/test size ratio, and that ratio is the one number the seesaw ties its size to. The transform and the matrix stay exactly as on file. Only the split changes.
- **Measurement:** Logistic L2 from make_models, with C rescaled per k so that the penalty per training patient matches the file: C_k = 0.04 * k/(k-1), which gives 0.08 at k=2, 0.05 at k=5 (the file value) and about 0.044 at k=10. k in {2, 5, 10}. Fold seeds 100-104 (stratified, as in random_splits). Arm W = whole-cohort ranks (rank_normalize_within_cohort, the file matrix). Arm I = training-only ranks (cohort_rank_inductive), on identical folds. Per fold, pair_auc gives same-cohort and cross-cohort pair AUC. Pairs are counted only inside one test fold (one model). Per seed, AUCs are pair-weighted over folds; then the mean is taken over the 5 seeds. D = same-pair AUC minus cross-pair AUC. Judged quantity: excess X(k) = D_W(k) - D_I(k). Mechanism readout: per cohort, the class-adjusted test offset o_g (for each class, the mean test log-odds in cohort g minus that class's mean over all test patients in the fold, averaged over the classes present), W minus I. It is regressed on (cohort bacterial share - pooled share). |S(k)| = absolute slope. Report only: D_I(k) on its own, the per-cohort training-side offset per k, k=5 W versus raw, and a fixed-C arm (C=0.05) at k=2 and 10 on seed 100. Run on both cohort sets. The 4-cohort set keeps the same 12,040 genes (per-cohort ranking does not depend on the other cohorts, so the other 4 cohorts' values do not change).
- **PASS:** In BOTH cohort sets, all of these hold: X(5) > 0; X(2) <= 0.5 * X(5); X(10) >= X(5); |S(2)| < |S(5)| < |S(10)|.
- **KILL:** In EITHER cohort set, any of these: X(5) <= 0; X(2) >= 0.8 * X(5) (the push factor drops from about 5 to about 2, but the excess barely shrinks); X(10) < X(2) (the excess does not grow with k). Anything between pass and kill = 'partial'. Partial counts as not passed.
- **Rival it separates:** Class-balance ranking arithmetic WITHOUT concentration (predicts X flat in k). Regularisation as the driver (C is rescaled per k, and the fixed-C arm shows separately what C adds). 'A different training size just fits differently' (the I arm has the same training sizes at every k, and only X is judged). GSE40396 label noise and hospital size (neither depends on k; the 4-cohort set is judged too).

### P2-mean-fixing-swap (about 12 CPU-minutes)

- **What it tests:** Whether the damage comes from the 'fixed hospital mean' property itself, and not from something special about percentile ranks, or from any transform that reads the whole hospital. Each transform is run in two versions: fitted on the WHOLE hospital (ties test patients to training patients), and fitted on the hospital's TRAINING patients only (no tie). The seesaw predicts: centring each gene per hospital (this fixes the mean exactly, with no ranks) shows a large tie effect. Scaling each gene by the hospital's spread without centring (reads the whole hospital, but leaves the mean free) shows almost none.
- **Why it is new:** Only rank transforms were ever run (raw, cohort rank, within-sample rank, training-only rank). No linear per-hospital transform has been run, and no 'reads the whole hospital but leaves the mean free' arm. Pairing each operation with its own training-only twin also removes E3's confound: E3 cut the tie AND changed the rank grid and reference set. Centring has no grid.
- **Measurement:** Logistic L2, k=5, C=0.05, the same fold seeds 100-104, both cohort sets. StandardScaler stays inside the model pipeline, as now. Arms: center_whole (subtract each gene's mean over the whole hospital, before the split); center_train (subtract each gene's mean over that hospital's training patients in the fold; applied to its train and test rows; reads no test patient); scale_whole (divide each gene by its whole-hospital SD, no centring); scale_train (the same, with the SD from training patients only). Report only: zscore_whole and zscore_train (per-hospital z-score, a common batch correction). D is defined as in P1. Tie effects: Y_c = D(center_whole) - D(center_train); Y_s = D(scale_whole) - D(scale_train); X(5) is reused from P1. Treatment-applied gate, computed before any model is fitted: the median over genes of the spread across hospitals of each gene's whole-hospital mean must be <= 1e-9 after center_whole, and must be >= 0.1 times the raw value after scale_whole. If the gate fails, that arm is void: not judged, and the prediction is not passed.
- **PASS:** In BOTH cohort sets: Y_c > 0; Y_c >= 0.5 * X(5); |Y_s| <= 0.25 * Y_c.
- **KILL:** In EITHER cohort set: Y_c <= 0; or Y_c < 0.25 * X(5) (the effect does not carry over to an exact mean-fixing transform, so it is rank-specific); or Y_s >= 0.5 * Y_c (reading the whole hospital hurts even when the mean is left free: generic transduction, not mean-fixing). Anything between = partial.
- **Rival it separates:** (1) 'A quirk of percentile ranks' (tail squashing, the (n+1)/2n wiggle, grid size): centring is linear and has no ranks. (2) 'Generic transduction: any transform that reads test patients ties them together': scale_whole reads the whole hospital but does not fix its mean. (3) 'E3 worked for another reason (re-ranking, reference size)': each pair compares the same operation with and without the tie, so nothing else changes.

### P3-additive-trees (about 10 CPU-minutes)

- **What it tests:** The fixed total holds for ANY model that adds up one-gene functions, not just a linear one. A LightGBM made of depth-1 trees (stumps) is gene-additive. The default LightGBM (15 leaves, gene interactions) is not. So the seesaw predicts: stumps get hit like logistic regression, and harder than the default LightGBM, even though stumps fit the training data LESS tightly. This is the seesaw's own answer to 'why is LightGBM hit less', and it is not compatible with 'logistic is hit because of its L2 penalty or its straight-line boundary'.
- **Why it is new:** LightGBM was only ever run with its file settings. 'LightGBM is hit less because a tree is not a weighted sum' was inferred afterwards and never tested. No earlier run separated additivity from 'trees are robust' or 'trees fit less'.
- **Measurement:** k=5, fold seeds 100-104, arms W and I as in P1, both cohort sets. Models: logistic (reuse P1's k=5 fits); LightGBM exactly as in make_models; LightGBM stumps (num_leaves=2, max_depth=1, everything else identical). Excess X_m = D_W - D_I per model. Algebra gate (computed in the W arm): spread across hospitals of the (train + test) mean score, divided by the pooled within-class SD of the scores. It must be <= 0.5 times the default LightGBM value for stumps AND for logistic. If not, the code is broken: stop, not interpretable. Fit gate: pinning P = slope of the per-hospital training-mean score on bacterial share, divided by the pooled within-class SD of the training scores. Gate passes if P_stumps >= 0.5 * P_lgbm_default. Fit gate failed = 'not interpretable', which counts as not passed and not as killed. Both gates are committed before the run and never relaxed.
- **PASS:** In BOTH cohort sets, with both gates passed: X_stumps > 0; X_stumps >= 2 * max(X_lgbm_default, 0); X_stumps >= 0.5 * X_logistic.
- **KILL:** In EITHER cohort set, with both gates passed: X_stumps <= X_lgbm_default (additive trees are hit no more than trees with interactions). Anything between = partial.
- **Rival it separates:** 'LightGBM is hit less because trees fit training less, only use order, or are robust to shifts' (predicts stumps hit as little as or less than the default, because they fit less). 'Logistic is hit because of its L2 penalty or linear boundary' (stumps have neither, so it predicts no excess for stumps). 'Any model that memorises training is hit' (predicts the default LightGBM is hit at least as hard as stumps).

### Controls, gates and fences

- Reproduction gate first: new code at k=5, fold seed 0, C=0.05 must reproduce the logistic cohort_rank rows of results/gapdiag_folds.csv (W arm) and results/gapdiag_inductive.csv (I arm), auc_within and auc_cross, to 1e-9. If not: stop, nothing is judged, and the drift is logged.
- Synthetic instrument check before any real data, built on _synthetic() from tests/test_gap_diagnostics.py. Positive control: shares 20/20/80/80 with whole-cohort ranks must give X(10) > X(2) > 0 and Y_c > 0. Negative control: all shares 50 % must give |X(k)| and |Y_c| below 0.25 times the positive-control X(5). If either fails, the real run is not started and the thresholds are not changed.
- Treatment-applied gates (P2) and algebra and fit gates (P3), as written in each prediction. A failed gate makes the arm void or not interpretable. It never counts as a pass.
- Training-only twins (I arm, center_train, scale_train) on identical folds and training sizes, so every judged difference is the train/test tie alone.
- Regularisation: C is rescaled per k in P1, and a fixed-C arm is reported. Stumps in P3 have no L2 penalty.
- No logistic fit may stop at its 3,000-iteration cap. Process CPU time is recorded per part in a run JSON.
- Every prediction runs on all 5 cohorts AND on the 4 without GSE40396, with the same 12,040 genes. Supported only if supported in both. A kill in either set is reported as 'not supported', naming the set.
- Tests in a new tests/test_seesaw.py, each checked to go red on a deliberately broken copy: (a) the fold runner really uses k folds of size about n/k, and C_k follows the formula (red if k is ignored or C stays fixed); (b) cross pairs are counted only inside one fold (red if pooled across folds with per-fold score offsets); (c) center_whole gives exact zero cohort means, and center_train / scale_train output is unchanged when one test row is perturbed (red if switched to whole-cohort statistics); (d) scale_whole gives cohort SD 1 and unequal cohort means; (e) on equal-size synthetic cohorts with whole-cohort ranks, stump cohort totals match to 1e-9 and 15-leaf totals do not (red if the stumps silently use the default settings); (f) adding a constant to one cohort's test scores moves o_g by that constant and leaves the same-cohort pair AUC unchanged; (g) an --evaluate-only verdict function applies the pass/kill/partial rules to a hand-made table (red if any threshold is flipped).
- Pre-registration: a dated PLAN entry in LABBOOK.md with every threshold, the seeds, the C_k formula, the gates and the fallback cut. It is committed on a feature branch (the repo is on master now; CORE section 7 forbids committing to master without the author) BEFORE any real-data run. The commit hash is written into the new runner's output JSON. A timing probe (one fold, one seed) runs only after that commit.
- Fold seeds 100-104 are fixed now. The checker's 2 extra seeds are NOT recorded in LABBOOK (it only says '2 extra random-CV fold seeds'). If the checker's transcript shows that it used any of 100-104, switch to 200-204 BEFORE the commit, never after seeing numbers.
- Pre-registered fallback: if the timing probe projects more than 60 CPU-minutes, fold seeds go from 5 to 3 (100-102) for every part. This is decided now, not after the probe.
- Fences: new files only (src/seesaw.py, run_seesaw_test.py, tests/test_seesaw.py, results/seesaw_*.csv). run_pipeline.py, config/cohorts.json and the existing results files stay untouched. No --fetch, no new cohorts, no RNA-seq, no gene selection anywhere (all 12,040 genes, StandardScaler inside the fold), CPU only.
- Wording fence: every whole-cohort arm is labelled a 'transductive diagnostic of a random-CV artefact'. None of it is called generalisation, leave-one-cohort-out is not judged here, and no rank-normalised AUC number leaves the repo.
- Logistic regression and stumps are the judged models. The default LightGBM appears only as the comparison arm in P3.
- Afterwards a fresh outside agent re-runs everything blind from data/ (the full run plus --evaluate-only) and must reach the same verdicts. Until then every result is 'self-assessed against pre-registered thresholds'.

### Rivals still left

- The seesaw is partly algebra for an additive model on mean-fixed features. These tests can show it is the DOMINANT cause of the deficit. That it exists is already guaranteed by algebra. A pass should be read as 'this algebra is what moves the AUC', not as a discovery.
- Class-balance ranking arithmetic without concentration: the smaller cross-hospital penalty left under training-only ranking (D_I). It is measured at each k but not explained by any prediction. In this framing it is the fuel, and none of the tests removes it.
- Hospital identity (fixed biology, age, severity, platform) is not retested. E2a (every hospital cut to 50/50 removes the effect) already speaks strongly against it as the main cause. Planner A's share-inversion test was dropped for that reason and for its small sizes (GSE40396 about 10 patients).
- Leave-two-cohorts-out (planner B's P3) was not kept: its kill rule could be fired by real differences between two unseen hospitals, which is not the seesaw's rival. E3's re-ranking confound is handled by P2's centring pair instead.
- C and k are not crossed. The lever depends on how tightly training patients are held, and only one rescaling rule is judged.
- Other additive linear models (for example a linear SVM) and stratifying folds by cohort x class are not tested.
- GSE40396 labels (13 afebrile virus-positive children) are still open (condition 2). The 4-cohort set guards the verdict, but the 5-cohort numbers inherit the label error.
- Multiple comparisons: 3 predictions x 2 sets x several sub-rules. A pass needs every sub-rule in both sets, so a false kill is more likely than a false pass. This is deliberate, but there is no formal correction.
- At k=10 a small hospital has about 4 test patients per fold (GSE40396 about 1 bacterial). Pair-weighting over folds and 5 seeds reduce the noise, and the kill rule uses X(10) < X(2), not X(10) < X(5), to guard against a floor effect. The |S(k)| ordering in the pass rule could still fail on noise.

### Risks

- P1 assumes the training-side push stays roughly the same across k. If the k=2 fits (half the data) hold training patients much less tightly even with rescaled C, P1 could come out partial for a reason that is not a rival. The per-k training offset is recorded to show this.
- scale_whole can make hospitals' raw averages differ by a hospital-specific factor. Pairing it with scale_train cancels most of this, but if Y_s is large, read the treatment gate before blaming the seesaw.
- Stumps may hold training patients too loosely. The fit gate turns this into 'not interpretable' instead of a false kill. That is also a way out, so the gate is committed first and never relaxed.
- In the 4-cohort set the range of shares shrinks (about 24-80 % instead of 19-80 %), so every effect may get smaller. The thresholds are relative to guard against this, but 'partial' gets more likely.
- The CPU estimate rests on the 19.10. timings (run_gapdiag.json: E1 about 156 s, E3 about 29 s). Stump and centring timings are not measured yet. The fallback is fixed in advance.
- Default LightGBM can move slightly with the thread count. Fix n_jobs in the runner and record it.
- The designers have seen E1-E4 and the mirror table. The predictions come from the algebra (the k factor, additivity, mean-fixing), not from new numbers, but the blind outside re-run is still essential.

**Verdict rule.** Supported only if all 3 predictions pass on BOTH cohort sets (all 5, and the 4 without GSE40396). One kill in either set = not supported, condition 3 stays open. Partial counts as not passed.

**In plain terms (the author's pick of two).** Each hospital gets one cake of fixed size to share among its patients, and a patient's slice is their 'bacterial' score. If the model gives the patients it learns from slices that fit them well, the few test patients of that hospital get whatever is left: crumbs in a mostly-bacterial hospital, far too much in a mostly-viral one. The test checks three things: the fewer test patients share the leftovers, the worse it must get; any rule that fixes the cake size must cause it, while a rule that leaves the size free must not; and if neither holds, the cake idea is wrong.

**Pages.** Condition 3 in plain words: (`reports/2026-10-19-condition-3.html`). This plan: (`reports/2026-10-19-seesaw-test-plan.html`).

## 2026-10-19 — PLAN, pre-registered before any re-run: GSE40396 label-rule fix + CORE §3 kill check

**Question.** Is the 2026-10-10 raw gap a GSE40396 label artefact? (CORE §3 first kill criterion.)

**Decision (the author, 19.10.).** The GSE40396 numbers in `docs/published-counts-candidates.md` were checked against the paper and are OK: febrile 22 viral (8 HHV-6, 8 adenovirus, 6 enterovirus) and 8 bacterial; afebrile 13 virus-positive and 22 virus-negative controls (Hu et al., PNAS 2013;110:12792-7, PMID 23858444). The 13 afebrile virus-positive children are DROPPED, as in the paper's own viral-vs-bacterial test.

**Rule change (the rule, never the numbers).** Add `"fever: no"` to `label_rules.exclude` in `config/cohorts.json`. Exclusions run first. Checked before the change: the string `fever: no` occurs only in GSE40396 (35 samples = 13 afebrile virus-positive + 22 virus-negative controls, the controls already excluded); no other cohort has a `fever:` field. Expected parse: GSE40396 8 bacterial / 22 viral / 35 excluded; all other cohorts unchanged; 476 → 463 labelled samples. GSE40396 `expected_counts` = {viral 22, bacterial 8}, so `--review` must print "verified".

**Kill rule (fixed now).** Under raw (the 2026-10-10 setting, `run_pipeline.py --train`): KILL if BOTH models' gap (random-CV mean AUC minus leave-one-cohort-out mean AUC) is ≤ 0.01 or negative. Old raw gaps: logistic 0.069, LightGBM 0.131. Rank-normalised and within-sample re-runs are reported, not judged. If KILL fires, the lab stops and the seesaw test does not run.

## Claim-check (2026-10-19, before the run)
Claim:      the raw cohort gap is not produced by the 13 mislabelled GSE40396 children.
Measured:   the same raw pipeline, same seeds, with those 13 removed; one run, 5 folds per scheme.
Rival:      fold noise — removing 13 samples reshuffles the random-CV folds, so the gap can move for that reason alone. Not killed by this run (one seed); it only matters if the new gap lands near the 0.01 line.
Judge:      the threshold above, committed before the run — fails if both gaps are ≤ 0.01.

## 2026-10-19 — RESULT: CORE §3 kill check after the GSE40396 rule fix — NOT killed

**Question.** Is the 2026-10-10 raw gap a GSE40396 label artefact? (Pre-registered in the PLAN entry above, 9ca61f9.)

**Run.** Rule fix committed in 6ad3881 (`"fever: no"` excluded; GSE40396 and GSE60244 now verified against their papers, GSE6269 still unverified). Then `run_pipeline.py --train`, `--train --rank-normalize`, `--train --within-sample-normalize`. 463 labelled samples (GSE40396 43 → 30), 12,040 genes, same seeds as before.

**Measured** (mean AUC over 5 folds per scheme; gap = random-CV minus leave-one-cohort-out):

| setting | model | random-CV old → new | leave-one-cohort-out old → new | gap old → new |
|---|---|---|---|---|
| raw (judged) | logistic L2 | 0.954 → 0.953 | 0.885 → 0.882 | 0.069 → **0.071** |
| raw (judged) | LightGBM | 0.956 → 0.954 | 0.825 → 0.839 | 0.131 → **0.115** |
| within-sample | logistic L2 | 0.961 → 0.959 | 0.900 → 0.900 | 0.061 → 0.059 |
| within-sample | LightGBM | 0.954 → 0.954 | 0.857 → 0.887 | 0.097 → 0.067 |
| cohort rank (transductive) | logistic L2 | 0.743 → 0.751 | 0.893 → 0.887 | −0.150 → −0.136 |
| cohort rank (transductive) | LightGBM | 0.904 → 0.909 | 0.895 → 0.891 | 0.008 → 0.018 |

**Verdict (pre-registered rule).** KILL needed both raw gaps ≤ 0.01. Measured 0.071 and 0.115, so **not killed**; the lab continues and the seesaw test may run.

## Claim-check (2026-10-19, after the run)
Claim:      the raw cohort gap is not produced by the 13 mislabelled GSE40396 children.
Measured:   with the 13 removed, one run (5 folds each scheme), the raw gap is 0.071 (logistic) and 0.115 (LightGBM), against a kill line of 0.01 written down before the run.
Rival:      fold noise from the reshuffled folds — inferred to be too small to matter here: 0.07 is about 7x the kill line, and per-fold AUCs on file move by a few hundredths. Not measured with extra seeds.
Judge:      the pre-registered threshold (9ca61f9). Self-assessed otherwise; no outside agent re-ran this.

**What we now believe.** The 2026-10-10 gap is not a GSE40396 label artefact. The LightGBM gap shrank by 0.016, and its within-sample gap by 0.030: part of LightGBM's held-out loss came from the 13 mislabelled children. The logistic numbers barely moved. The negative logistic gap under cohort ranking is still there (−0.136), so the seesaw test is still needed. Condition 2: 4 of 5 cohorts verified (GSE6269 waits for the author's check of the per-pathogen split, which comes from GEO titles, not the paper).

**Where.** `results/summary_{raw,ranknorm,withinsample}.csv`, `results/run_*.json`, `results/label_review.md`.

## 2026-10-19 — RESULT: seesaw test STOPPED at the synthetic instrument check (positive control failed)

**Question.** Before any real-data run, does the instrument see the seesaw where we built it in? (PLAN 9daaf08, "Synthetic instrument check": the positive control must give X(10) > X(2) > 0 and Y_c > 0; the negative control must stay below 0.25 × the positive X(5). "If either fails, the real run is not started and the thresholds are not changed.")

**Run.** Code e535f1d + bcb52b9 (`src/seesaw.py`, `src/seesaw_verdict.py`, `run_seesaw_test.py`; tests a-g each seen red on a broken copy, 10 mutations). `python run_seesaw_test.py --parts synthetic`: `_synthetic()` from tests/test_gap_diagnostics.py, 4 cohorts × 60 patients, 60 genes, fold seeds 100-104, the same runner and aggregation as the real run. Output in `results/run_seesaw.json` → `synthetic`.

**Measured** (logistic, C_k rescaled; X = excess of whole-cohort ranks over training-only ranks in D = same-pair AUC − cross-pair AUC):

| control | X(2) | X(5) | X(10) | Y_c | pre-registered check |
|---|---|---|---|---|---|
| positive (shares 20/20/80/80) | 0.082 | 0.084 | 0.079 | 0.077 | X(10) > X(2) > 0 and Y_c > 0 → **FAIL** (X(10) < X(2)) |
| negative (all 50 %) | −0.002 | −0.000 | −0.001 | −0.000 | all \|·\| < 0.25 × 0.084 → pass |

Report-only details from the same positive control (seeds 100-104): D_W ≈ 0.22 and D_I ≈ 0.13-0.14 at every k; D_raw = −0.005. |S(k)| = 0.95, 0.89, 0.77 for k = 2, 5, 10 (falling, not rising). Per-cohort offsets: in the W arm the training offset and the test offset have the SAME sign and similar size at every k (for example cohort C0, 20 % bacterial, k=5: o_train +1.11, o_test +1.29). The I arm shows the same pattern, slightly smaller.

**Verdict (pre-registered rule).** The positive control failed, so the real-data run is NOT started and no threshold is changed. Condition 3 stays open.

## Claim-check (2026-10-19)
Claim:      the seesaw instrument cannot confirm its k-lever even where the effect is built in.
Measured:   on synthetic data with a known class-share imbalance, the excess X is flat in k (0.082 / 0.084 / 0.079), and training and test offsets point the same way.
Rival:      a bug in the new runner. Against it (inferred, not proven): the raw arm gives ~0, the negative control gives ~0, the training-only arm is lower than the whole-cohort arm, and tests a-g (k really used, C_k, pairs inside one fold, train-only twins read no test row) each go red on a broken copy. Not yet checked by anyone else.
Judge:      the pre-registered positive-control rule (9daaf08). Self-assessed otherwise.

**What we now believe (inferred).** On the synthetic data, whole-cohort ranking does add cross-cohort damage (X ≈ 0.08), and centring each cohort's genes on the whole cohort does too (Y_c ≈ 0.08). But the damage does not grow when the test share shrinks, and the training patients carry the same cohort shift as the test patients. That is not the seesaw as written ("training held close, the few test patients take k times the push"). It looks more like a shift shared by the whole cohort, which a cross-cohort pair sees and a same-cohort pair does not. This is a reading of one synthetic design, not a tested claim. Next step needs the author: see the options sent to the author.

**Where.** `results/run_seesaw.json` (synthetic block), code `src/seesaw.py`, `src/seesaw_verdict.py`, `run_seesaw_test.py`, `tests/test_seesaw.py`.

## 2026-10-19 — gap diagnostics re-based on the new labels (463 samples)

**Run.** `python run_gap_diagnostics.py` (all parts: E1 137 s, E3 26 s, E4 43 s, E2a 864 s, E2b 892 s, mirror 58 s), so every committed `results/gapdiag_*` file uses the fixed GSE40396 labels.

**Measured.** No verdict word changed against the 476-sample run (E0 reproduced; P1 supported; LightGBM direction as predicted; P4 KILLED; P1b supported; iteration cap 0). The numbers moved slightly: under cohort ranking, logistic same-cohort / cross-cohort / pooled AUC 0.954 / 0.680 / 0.743 → 0.956 / 0.686 / 0.751; inductive pooled 0.916 → 0.919; P2 gap mean 0.066 → 0.069; P3 gap mean −0.118 → −0.113.

**Belief.** Unchanged. These are the 19.10. post-hoc diagnostics on the corrected labels. They are not the seesaw test, which is stopped (entry above).

**Where.** `results/gapdiag_*.csv`, `results/run_gapdiag.json`.

## 2026-10-19 — outside audit of the stopped seesaw test: NO BUG; the k-lever is not supported by its own positive control

**Question.** Is the failed synthetic positive control (entry above) a bug, or a real property of the design? (the author's choice A.)

**Run.** A fresh agent, given only the PLAN (9daaf08) and the code, worked read-only and blind to the RESULT entries. It reviewed the code against every plan definition and wrote its OWN minimal X(k) implementation. It then ran probes of its own choice.

**Measured (by the auditor).**
- **No bug** that could change the verdict. Two harmless differences: S(k) averages over folds and seeds before taking W − I, and the fixed-C/raw arms are separated from the judged arms only by a C filter (correct, but fragile).
- Its own code matches the repo runner to 4 decimals. Positive control X(2/5/10) = 0.0825 / 0.0840 / 0.0792; negative control ≈ 0.
- Probes (seeds 100-102):

| probe | X(2) | X(5) | X(10) |
|---|---|---|---|
| 150 patients per cohort | 0.036 | 0.032 | 0.036 |
| another data seed | 0.081 | 0.087 | 0.091 |
| C = 10 | 0.055 | 0.085 | 0.084 |
| 2,000 genes (more genes than patients, like the real data) | 0.37 | 0.44 | 0.44 |

  With 2,000 genes, X rises from k=2 to k=5 but not from k=5 to k=10. The P1 kill rule would fire there too (X2 ≥ 0.8·X5).
- Offsets: in every setting, the whole-cohort-arm training offsets and test offsets have the SAME sign (low-share cohorts +, high-share cohorts −). The W − I difference on the training side shrinks about as 1/(k−1) (−0.27, −0.06, −0.03), so the fixed cohort total is real algebra. But the training side absorbs it, and the test-side difference stays about +0.21 to +0.27 at every k.

**Verdict.** **The seesaw k-lever is not supported by its own positive control.** No bug was found by an outside check. The pre-registered seesaw test is closed without a real-data run. Condition 3 stays open.

**What we now believe (inferred).** The cross-cohort damage looks like a shared shift: every patient of a cohort, in training and in test alike, is moved in the same direction by an amount set by the cohort's class mix. This hurts only pairs from different cohorts. The seesaw amplification by k does not appear. A new pre-registered test of this reading is DRAFTED in `docs/2026-10-19-shared-shift-plan-DRAFT.md`. It is not approved, and nothing runs until the author approves it.

## 2026-10-19 — PLAN, pre-registered before any run: the shared-shift test (CORE §2 condition 3)

**Status.** Approved by the author 19.10.2026 (answer A), as drafted in dbd09f7. Text below copied unchanged from `docs/2026-10-19-shared-shift-plan-DRAFT.md` (minus its title and draft note). Plain-terms metaphor: not yet picked.

### Claim

Under per-cohort ranking, every patient of a cohort (training AND test) is moved
in the same direction by an amount set by that cohort's class mix (low
bacterial share → scores pushed up, high share → pushed down). A constant shift
per cohort cannot change same-cohort pairs, but it mis-sorts cross-cohort pairs.
That shared shift is the MAIN reason the logistic random-CV AUC under
`--rank-normalize` falls below leave-one-cohort-out.

Rival (the old seesaw): training patients held close, test patients pushed the
other way. The two readings make OPPOSITE predictions about the sign of the
test shift relative to the training shift. H1 uses exactly that.

### Setup (shared by all predictions)

Logistic L2 from make_models, k=5, C=0.05, fold seeds 100-104, both cohort sets
(all 5, and the 4 without GSE40396), 12,040 genes, StandardScaler inside the
fold. Arms: W (whole-cohort ranks), I (training-only ranks), raw. D = same-pair
AUC − cross-pair AUC (pairs inside one fold only, pair-weighted, mean over
seeds), as in `src/seesaw_verdict.pair_weighted`. Code reuse: `src/seesaw.py`
(runner, offsets), new code only for the two new arms below.

### H1 — the training shift predicts the test shift (sign + transfer)

- **Measurement:** per fold and cohort, o_train and o_test (class-adjusted
  offsets, as now). Arm "W+recentre": the W model's test scores of cohort g are
  shifted by −o_train(g), which is estimated from TRAINING scores only. Judged
  quantity R = share of the cross-pair damage that recentring removes:
  R = (cross_recentred − cross_W) / (cross_raw − cross_W).
- **PASS (both sets):** sign(o_test) = sign(o_train) for every cohort in the W
  arm (mean over folds and seeds) AND R ≥ 0.5.
- **KILL (either set):** R ≤ 0 (recentring with the training shift makes cross
  pairs worse, which is the seesaw's sign).
- Between = partial.

### H2 — the class mix of the reference is the lever

- **Measurement:** new inductive arm "I-balanced": each patient is ranked
  against its cohort's TRAINING rows, re-weighted to 50 % bacterial / 50 % viral
  (training labels only; no test row is read). Compared with I (same folds).
- **PASS (both sets):** D(I-balanced) ≤ 0.5 · D(I) AND the absolute slope of
  o_test on (share − pooled share) drops by ≥ 50 % against I.
- **KILL (either set):** D(I-balanced) ≥ 0.9 · D(I).
- Between = partial.
- Why it is new: E2a balanced the cohorts by DROPPING patients (it changed the
  data). H2 keeps every patient and changes only what they are ranked against.

### H3 — no k-lever (the stopped prediction, now as a stated null)

- **Measurement:** X(k) = D_W − D_I at k = 2, 5, 10, C_k rescaled, as in P1.
- **PASS (both sets):** |X(10) − X(2)| ≤ 0.25 · X(5) and X(5) > 0.
- **KILL (either set):** X(10) ≥ 1.5 · X(2) (the seesaw's growth reappears on
  real data).
- Between = partial.

### Gates (before any real-data run, never relaxed)

1. Synthetic instrument check on `_synthetic()` at **60 AND 2,000 genes**:
   positive control (shares 20/20/80/80) must give H1 PASS and H2 PASS. Negative
   control (all 50 %) must give |D_W − D_raw| and |D_I − D_raw| < 0.25 × the
   positive-control value. If this fails: stop, do not change thresholds.
2. Reproduction gate: the reused runner reproduces the seed-0 rows of
   `results/gapdiag_folds.csv` / `gapdiag_inductive.csv` to 1e-9 (already true
   by construction, re-checked).
3. Tests that go red on broken code: I-balanced reads no test row; the reference
   weights sum to equal class mass per cohort; recentring uses training scores
   only (red if it uses test scores); the verdict function on a hand-made table.
4. No logistic fit at its 3,000-iteration cap.

### Verdict rule

Supported only if H1, H2 and H3 all PASS in BOTH cohort sets. A kill in either
set = not supported, and condition 3 stays open. Afterwards a blind outside agent
re-runs from data/ and must reach the same verdicts.

### Budget and fences

About 25-30 CPU-minutes (P1-style k sweep ≈ 15, the rest ≈ 10). New files only
(`src/shared_shift.py`, `run_shared_shift_test.py`, `tests/test_shared_shift.py`,
`results/sharedshift_*`). No cohort, label or config change. No rank-normalised
AUC leaves the repo; every whole-cohort arm is a transductive diagnostic.

### Rivals still left

- Cohort biology/platform that happens to line up with class share (5 cohorts
  only; E2a speaks against it, but nothing here removes it).
- The negative gap could still be partly class-weighting arithmetic in
  `class_weight="balanced"`; not varied here.
- 5 cohorts = 5 points for any slope; the sign rule in H1 is per cohort, but
  the slope in H2 is fragile.
- Designers have seen the synthetic and audit numbers (not real-data seesaw
  numbers, which were never run). The 10-seed sweep is a
  separate flow job and is not used to set any threshold here.

## 2026-10-19 — RESULT: shared-shift test STOPPED at its synthetic gate (fails at 2,000 genes); real run refused

**Question.** Before any real-data run: does the instrument see the shared shift where we built it in, at 60 AND at 2,000 genes? (PLAN 2fd8865, gate 1.)

**Run.** Code 665db45 + fix 61d154a (both repro gates filtered on `DataFrame.transform`, the method, instead of the column, so they matched 0 rows. Neither gate had produced a number before; nothing on file is affected). Guards: 19 new tests, 5 mutations each seen red; 94 tests green. Added before any run: a VOID band when the reference damage is ≤ 0 (cross_raw ≤ cross_W for H1, D(I) ≤ 0 for H2). It did not fire anywhere below. `python run_shared_shift_test.py --parts synthetic repro`. Output: `results/run_sharedshift.json`.

**Measured.**
- Reproduction gate: `run_k5` at fold seed 0 reproduces the logistic W and I rows on file, max difference 1.1e-16 → **pass**.
- Synthetic gate (positive control shares 20/20/80/80, negative control all 50 %, seeds 100-104):

| genes | H1 (R, signs agree) | H2 (D_I → D_I_bal; slope I → I_bal) | H3 (X2 / X5 / X10) | negative control | gate |
|---|---|---|---|---|---|
| 60 | PASS (R = 0.90, yes) | PASS (0.143 → −0.004; 3.0 → 0.5) | PASS (0.082 / 0.084 / 0.079) | pass | pass |
| 2,000 | **PARTIAL** (R = 0.006, yes) | **PARTIAL** (0.122 → −0.097; 1.4 → **6.1**) | PASS (0.41 / 0.46 / 0.46) | pass | **FAIL** |

  At 2,000 genes the whole-cohort arm is hit far harder (cross-pair AUC 0.19, against 0.91 raw). Recentring each cohort by its training offset leaves it at 0.20.

**Verdict (pre-registered rule).** The positive control must give H1 PASS and H2 PASS at both gene counts. At 2,000 genes it gives PARTIAL and PARTIAL, so the gate fails. The real-data run is refused, and no threshold is changed. Condition 3 stays open.

## Claim-check (2026-10-19)
Claim:      the shared-shift instrument does not confirm its own reading in the many-genes regime.
Measured:   on synthetic data with 2,000 genes (more genes than patients, like the real 12,040 × 463), a per-cohort constant correction learned on training scores removes ~1 % of the cross-pair damage. At 60 genes it removes 90 %.
Rival:      a bug specific to 2,000 genes. Against it (inferred): the same code passes at 60 genes, the negative control is clean at both sizes, and the reproduction gate matches the file to 1e-16. Not checked by an outside agent.
Judge:      the pre-registered gate rule (2fd8865). Self-assessed otherwise.

**What we now believe (inferred, not tested).** With few genes, whole-cohort ranking acts like one constant shift per cohort: same sign for training and test, set by the class mix, and removable. In a viral-heavy cohort every patient looks MORE bacterial (positive offset); this direction is checked against our own offset tables. With more genes than patients the damage is much larger and is NOT a constant per-cohort shift: the model seems to use the ranking pattern itself, which separates training from test patients, not just cohorts. Neither the seesaw nor the shared-shift reading explains the regime the real data live in. Condition 3 ("explain the negative logistic gap by a test that was run") is still open. The honest summary for condition 4 may have to be: "rank-normalised random-CV numbers are a transductive artefact whose mechanism we characterised on synthetic data but did not pin down on the real cohorts".

**Where.** `results/run_sharedshift.json`; code `src/shared_shift.py`, `src/shared_shift_verdict.py`, `run_shared_shift_test.py`, `tests/test_shared_shift.py`.

## 2026-10-19 — DECISION (the author): condition 3 rewritten; branch merged into master

**What changed.** CORE §2 condition 3 now reads: "the rank-normalised random-CV number is a transductive artefact: shown on synthetic data, mechanism on the real cohorts unresolved" (dated, marked as the author's decision; the old wording is kept in CORE). The author gave their CORE §7 OK to merge `plan/seesaw-test` into `master` (local repo, no remote).

**Status of CORE §2.** 1 met (within-sample control). 2: 4 of 5 cohorts verified (GSE6269 waits for the author's own check of the per-pathogen split). 3: met in its rewritten form. 4 (the one answer sentence, raw + within-sample numbers only, claim-check first): NOT yet written → next session.

## 2026-10-20 — RESULT: CORE §2 condition 4, the one answer sentence (the author's call)

> **SUPERSEDED 2026-10-22** by the 458-array restatement in the entry "RESULT: within-sample 10-seed sweep" at the bottom of this file. The numbers below are the 463-array, single-run ones and are kept for history only.

**Question (§1).** Of a within-cohort (random-CV) AUC, how much survives a hospital cohort the model never saw?

**Sources, nothing else.** `results/summary_raw.csv` and `results/summary_withinsample.csv` (commit 36cecca, run date 2026-10-19, 463 samples, 12,040 shared genes, 5 cohorts, after the GSE40396 rule fix). Rank-normalised numbers are NOT used (condition 3: transductive artefact). Means over 5 folds / 5 held-out cohorts; no seed sweep is on file, so these are single-run numbers.

## Claim-check (2026-10-20)
Claim:      most of the within-cohort AUC survives a held-out cohort, but a real gap of 0.06-0.12 AUC remains, and the leak-free within-sample transform shrinks it without removing it.
Measured:   mean AUC random-CV vs leave-one-cohort-out, both models, raw and within-sample scales, one run each (files above).
Rival:      (a) the gap is a label-parsing artefact — against it: 4/5 cohorts verified against published counts, kill check 19.10. not triggered (gaps 0.071/0.115 vs pre-registered line 0.01); GSE6269 still unverified (the author's own check pending). (b) single-run noise — not ruled out, the 10-seed sweep never ran; the spread across the 5 held-out cohorts (raw logistic 0.82-0.94, raw LightGBM 0.70-0.89) is the only variance on file. (c) co-infection/non-infectious exclusions flatter both numbers equally (README).
Judge:      none outside the lab yet. The numbers are machine-written CSVs; the sentence below is self-written. A blind outside re-read of this entry against the two CSVs is the open step.

**The sentence.**
Of a random-CV AUC of 0.95 (logistic L2 0.953, LightGBM 0.954), **0.88 survives a held-out cohort on the raw scale** (logistic 0.882, gap 0.071; LightGBM 0.839, gap 0.115) and **0.90 survives under within-sample normalisation** (logistic 0.900, gap 0.059; LightGBM 0.887, gap 0.067) — every number from `results/summary_raw.csv` and `results/summary_withinsample.csv` at commit 36cecca.

**Plain words.** Roughly 7 to 12 AUC points of the headline number are hospital-specific structure, not biology the model carries to the next hospital. Normalising each patient against their own genes only (legal at the bedside) gives back about one to five of those points; the rest stays. The simpler model loses less.

**Status of CORE §2.** 1 met · 2 four of five (GSE6269 open) · 3 met in its rewritten form · **4 written here**. Open: GSE6269 per-pathogen check (the author), seed sweep (would turn "single run" into "reproducible"), an outside re-read of this entry.

**Where.** This entry; branch `cond4-sentence` (no master commit without the author, CORE §7).

## 2026-10-20 — Blind re-read of condition 4

**Question.** Does the condition-4 sentence (entry above, merge 74cbb73) say what the two CSVs say? An outside reader who has not seen the sentence writes its own from the tables; we compare number by number.

**Design.** A fresh sub-agent with no file access got ONLY the §1 question, the plain meaning of the columns, and the full text of `results/summary_raw.csv` and `results/summary_withinsample.csv` (commit 36cecca). It was not shown the LABBOOK sentence. Note the limit: same model family as the author, so this is an independent *read*, not an independent *judge* of whether the experiment is sound.

**Blind agent's sentence, verbatim.** "Of a random-CV AUC of 0.954 (lightgbm) and 0.953 (logistic_l2), 0.839 and 0.882 survive a held-out cohort (raw) and 0.887 and 0.900 survive (within-sample normalised), so the gap is 0.115 (lightgbm) and 0.071 (logistic_l2) raw, shrinking to 0.067 and 0.059 within-sample."

**Number by number.**
| number | LABBOOK sentence | blind read | verdict |
|---|---|---|---|
| random-CV (raw) logistic / LightGBM | 0.953 / 0.954 | 0.953 / 0.954 | match |
| raw held-out logistic / LightGBM | 0.882 / 0.839 | 0.882 / 0.839 | match |
| raw gap logistic / LightGBM | 0.071 / 0.115 | 0.071 / 0.115 | match |
| within-sample held-out logistic / LightGBM | 0.900 / 0.887 | 0.900 / 0.887 | match |
| within-sample gap logistic / LightGBM | 0.059 / 0.067 | 0.059 / 0.067 | match |
Headline roundings (0.95 → 0.88 raw, 0.90 within-sample) also agree with the CSVs (logistic held-out values).

**Claim-check.** Claim: the sentence's numbers are what the files say. Measured: 10 of 10 numbers match a fresh read. Rival: the blind agent just copied the same obvious cells — a weak test for a transcription claim, which is all this checks. Judge: the CSVs themselves; not a judge of the science (single run, no seed sweep, GSE6269 labels still open — unchanged).

**One wording caveat (not a wrong number).** The within-sample gaps (0.059 / 0.067) are measured against the within-sample random-CV AUCs (0.959 logistic / 0.954 LightGBM), not the raw 0.953 / 0.954 the sentence quotes. 0.953 − 0.900 would be 0.053, not 0.059. The sentence is correct but a reader could mis-subtract; the sentence was NOT changed (only a wrong number would change it).

**Verdict: holds.** No number corrected. Status of CORE §2 unchanged (GSE6269 open).

**Where.** This entry, commit `ebb8c8e`; branch `a side branch` (not merged to master; the author's call).

## 2026-10-21 — RESULT: GSE6269 labels verified (CORE §2 condition 2 → 5 of 5, AI-verified)

**Question.** Do our GSE6269 (GPL96) labels match the published per-pathogen counts? Last open item of condition 2.

**Who checked.** Claude (AI), not a human. First check 2026-10-20, re-read 2026-10-21 from the PMC article page (the Europe PMC full-text XML is not open access; the PMC HTML page was read). The author's decision (21.10.): this counts as 5 of 5, labelled "AI-verified, not human".

**What was compared.** Ramilo et al., Blood 2007, PMID 17105821, [doi:10.1182/blood-2006-02-002477](https://doi.org/10.1182/blood-2006-02-002477), PMC1801073:
- Methods, "Patient information": 29 E coli, 50 S aureus, 16 S pneumoniae, 36 influenza A, 7 healthy (all sets, all platforms).
- Table 1 (set 1, training set) footnote: E coli n = 29, S aureus n = 31, S pneumoniae n = 13, influenza A n = 18 → **73 bacterial, 18 viral**.
- Results: "These data, obtained from a total of 91 patients, were generated using Affymetrix U133A GeneChips" (U133A = GPL96, the platform we keep). 73 + 18 = 91.
- GEO "Overall design": 143 samples in the series (144 minus PBMC_Healthy_INF295, no CEL file).

**Measured.** `python run_pipeline.py --review`: GSE6269 parsed bacterial 73, viral 18, excluded 6 (the healthy controls; 7 minus the missing INF295) → "verified against published counts". All five cohorts now print verified (`results/label_review.md`).

**What we now believe.** Condition 2 holds, with the honest label: the published count check for GSE6269 was done by an AI reading the paper. The README "what I checked vs what AI did" section must list it under AI. Kill criterion 1 (label artefact) is not triggered: the parse was not changed, it already matched.

**Where.** `config/cohorts.json` (GSE6269 `expected_counts` + `source_of_expected`, which names who checked), `results/label_review.md`; branch `v1/release`.

## 2026-10-21 — PLAN: 10-seed sweep of the raw cohort gap (run by the author's decision)

**Question.** Is the raw gap (logistic 0.071, LightGBM 0.115; `results/summary_raw.csv`, fold seed 0) a fold-seed accident? CORE §6 calls the gap "reproducible" but every run used one random-CV fold seed.

**Design (fixed before the run).**
- Data: the 463-sample, 12,040-gene raw matrix (after the GSE40396 rule fix; labels now 5/5 verified). Raw scale only, no rank normalisation, no within-sample transform.
- Sweep: random-CV fold seeds **200-209** (0, 1, 2 and 100-104 were used before). LightGBM model seed **fixed at 0**. Both models (logistic L2, LightGBM), same hyperparameters as `src/train.py`.
- Leave-one-cohort-out has **no fold randomness** (one fold per cohort). With a fixed model seed it must be identical in all 10 rows; the script stops if it is not.
- Side check (the card's control): LightGBM LOCO with **model seeds 200-209** at fold seed 0, to show how much LightGBM's own randomness moves the held-out number. Logistic is deterministic, so it has no such check.
- Reproduction gate first: fold seed 0 + model seed 0 must reproduce `results/summary_raw.csv` to 1e-9, else stop.
- Outputs: `results/seedsweep_raw.csv` (20 rows = 10 seeds x 2 models: random-CV AUC, LOCO AUC, gap), `results/seedsweep_raw_summary.csv` (mean/SD/min/max of the gap per model + seeds at or below the kill line), `results/seedsweep_lgbm_modelseed.csv`, `results/run_seedsweep.json`.
- **Kill rule (pre-registered):** if the gap is ≤ 0.01 for BOTH models in a majority (≥ 6) of the 10 seeds, the 0.071 / 0.115 numbers were a seed artefact → stop, write it down, report to the author (CORE §3 first kill criterion, repeated across seeds).
- The RESULT entry reports the full spread, the old single-seed values beside it, and whether each model's gap stays above 0.01 in all 10 seeds.

**Code.** `run_seed_sweep.py`; `src/train.py` gains `seed` (fold seed) and `model_seed` arguments (defaults 0 = old behaviour). Guards: `tests/test_seed_sweep.py`, 6 tests; 4 mutations seen red (seed dropped, model seed ignored, kill rule any() instead of all(), half counted as majority). 100 tests green.

**Cost.** CPU only, € 0, est. 10-30 CPU-min. Free RAM at plan time: 15 GB (line: 4 GB).

## 2026-10-21 — Three outside CODE reviews: every finding, fixed or rejected

**Question.** Does anyone outside the author find a bug in the code that would change a reported number? (v1 step 3, the author's decision.)

**Design.** One identical input for all three reviewers: `reports/review/2026-10-21/bundle.txt` = 29 files (src, run_*, tests, config, requirements) from commit d912d21, code and config only (no data, no results, no LABBOOK/README prose); hashes in `manifest.md`; instructions in `prompt.md` (leakage, labels, splits, AUC, harmonisation, tests that cannot fail; each finding with severity, failure scenario, red test).
- **gpt-6-astra** (OpenAI Responses API, reasoning high), via a separate review script (not included in this copy). **CORE §7 exception:** §7 says never "use a paid API"; the author OK'd this one call for v1 (decision, 21.10.). Only code went out; 41,937 input / 14,956 output tokens. → `astra.md`
- **Gemini** `gemini-3.1-pro-preview` via a separate API script (not included in this copy). → `gemini.json` / `gemini.md`
- **Blind fresh Claude agent**, told to read only `prompt.md` + `bundle.txt` (same model family as the author: an independent read, not an independent family). → `blind-claude.md`

Verdict lines: astra HIGH-FOUND (3 HIGH), Gemini HIGH-FOUND (4 HIGH), blind Claude NO-HIGH. Every reviewer-HIGH was checked against the real data before it was accepted.

**Findings, merged (A = astra, G = Gemini, C = blind Claude; H/M/L = their severity).** "Fixed" = code changed + a test that goes red when the fix is undone (10 mutations, each seen red; 108 tests green).

| # | finding | who | checked on the real data | outcome |
|---|---|---|---|---|
| 1 | Published-count MISMATCH is only a note, training continues | A3 H, G3 H, C1 M | all 5 cohorts verified, no MISMATCH → no number changed | **fixed**: `labels.refuse_mismatches`, hard stop in `cmd_train` and `load_matrix` |
| 2 | Imputation by whole-cohort gene mean before any split (and 0.0 fill) | A2 H, G4 H, C2 M | **0 NaN cells** in the model matrix, all 5 cohorts → the fill never fires | **fixed (on record)**: `build_matrix` writes the imputed-cell count into its notes, so a future NaN is visible |
| 3 | Highest-mean probe per gene chosen over every patient, test patients included | A1 H, G4 H, C6 L | label-free, but an unsupervised choice that reads test X | **control run** (PLAN below): lowest-probe-id rule that reads no value, `to_symbols(collapse="first")` |
| 4 | `prefer_platform` not enforced (combined matrix tried first; fallback to another platform; "GPL96" substring of "GPL961"; cache not checked) | G1 H, A6 M, C5 L | GSE6269 parsed as GPL96, 97 samples (not 143); GEO has no combined matrix for it → no number changed | **fixed**: preferred platform skips the combined file, exact id match, no fallback; `_load_all` refuses a parsed platform that differs from config |
| 5 | "influenza" rule catches *Haemophilus influenzae* (a bacterium) as viral | G2 H | 0 samples mention haemophilus/influenzae → no number changed | **fixed**: "haemophilus", "influenzae" added to the bacterial rules (checked before viral); `--review` unchanged, 5/5 |
| 6 | Negated text ("bacterial culture: negative", "co-infection: no") matched as a diagnosis | A7 M, C1 M | 0 usable samples contain "negative" or "culture"; class counts match the papers in 5/5 | **rejected for v1** (no such string in the data); stays a known limit of substring rules |
| 7 | X / y / cohort row alignment not checked in the training path | A8 M | `build_matrix` builds all three from one index | **fixed**: `run_both_schemes` refuses unequal indexes |
| 8 | Undefined-AUC folds silently dropped; incomplete studies still write headline files | A9 M, C4 L | `summary_raw.csv` folds = 5 for every model and scheme | **fixed in the sweep** (`sweep_rows` refuses incomplete folds; the LightGBM side check refuses a None AUC); headline files already carry the fold count |
| 9 | Logistic `ConvergenceWarning` suppressed | A10 M | `results/gapdiag_folds.csv` raw: max n_iter = 160 of 3,000 | **rejected** (no fit reached the cap) |
| 10 | Repro gates fail open on NaN / missing rows | A11 M | seed-sweep gate gave 1.1e-16 | **fixed in the seed sweep** (non-finite or model-set mismatch = fail); the seesaw / shared-shift gates belong to stopped tests, not v1 claims → left |
| 11 | `--parts probe` overwrites seesaw results (`startswith("p")`) | A12 M | real bug in a stopped test; the files on disk were not overwritten this way | **fixed** (explicit p1/p2/p3); no test (stopped experiment, one-line change) |
| 12 | Log2 decision from one series maximum (> 100) | A5 M, C6 L | series maxima 15.2 / 34,377 / 30,191 / 15.4 / 103,612: far from 100 either way | **rejected** (no decision is near the threshold) |
| 13 | Gene eligibility depends on all-NaN probes in held-out data | A4 M | a per-series deposit property, no label or value use; 0 NaN in the shared genes | **rejected**, documented as how "shared genes" is defined |
| 14 | The gap mixes estimands: random-CV pools cross-cohort pairs; LOCO is an unweighted mean of 5 cohorts | C3 M | `gapdiag_folds.csv` / `gapdiag_cohorts.csv`, raw, fold seed 0 (logistic / LightGBM): random-CV within-cohort-pair AUC 0.954 / 0.951 vs pooled 0.953 / 0.954; size-weighted LOCO 0.896 / 0.847 vs unweighted 0.882 / 0.839; per-cohort paired gap (in-cohort random-CV minus LOCO) mean 0.052 / 0.089, **positive in all 5 cohorts for both models** | **accepted as a wording caveat**: the README reports the paired per-cohort gap beside the headline gap |
| 15 | Smoke test passes a constant-0.5 scorer; LightGBM seed test trivial without lightgbm; NaN-fill test checks only "no NaN" | A13 L, C8 L | — | **fixed**: smoke test requires random-CV > 0.8; `importorskip`; the fill must equal 0.5 |
| 16 | LightGBM `subsample=0.8` is inert without `subsample_freq` | C8 L | true (LightGBM bagging needs freq > 0) | **not changed** (it would change the headline model); the model-seed side check therefore measures `colsample_bytree` randomness only |
| 17 | `cohort_rank_inductive` scales train rank/n vs test rank/(n+1) | C7 L | diagnostic arms of stopped tests only | **rejected for v1** (not on the headline path) |

**Process note.** The seed sweep (PLAN entry above) was already running on the code of d912d21 when these fixes were made; the running process had imported that code. Its `run_seedsweep.json` `commit` field is read at write time and so names a later state. None of the fixes changes the raw numbers it measures (no MISMATCH, 0 NaN cells, GPL96, aligned rows, 5/5 folds).

**What we now believe.** No reviewer found a bug that changes a number on file. Every reviewer-HIGH, checked against the data, either never fires on these five cohorts or is now a hard stop. The one open HIGH is #3, the probe choice, which needs a measurement, not an argument.

## 2026-10-21 — PLAN: probe-choice control (review finding #3)

**Question.** Does the headline gap depend on choosing each gene's probe by its mean over all patients of the series (test patients included)?

**Design.** `run_probe_control.py`: the raw headline setting (fold seed 0, LightGBM model seed 0, 463 samples), both models, twice: `collapse="mean"` (current) and `collapse="first"` (lowest probe id per gene, reads no expression value). Same shared-gene axis by construction (every gene keeps one probe either way). Output `results/probecontrol_raw.csv`, `results/run_probecontrol.json`.

**Rule (fixed now).** If either model's gap moves by more than **0.02** between the two rules, the probe choice materially moves the headline: report both, and the README headline gets the caveat. Otherwise finding #3 is rejected with this number.

**Cost.** CPU, ~2-4 min, € 0.

## 2026-10-21 — RESULT: 10-seed sweep of the raw cohort gap — NOT a seed artefact

**Question.** Is the raw gap (logistic 0.071, LightGBM 0.115 at fold seed 0) a fold-seed accident? (PLAN d912d21, committed before the run.)

**Run.** `python run_seed_sweep.py` (code of d912d21, see the process note in the review entry), 463 samples × 12,040 genes, raw scale, both models. Repro gate: fold seed 0 + model seed 0 reproduces `results/summary_raw.csv`, max difference 1.1e-16 → pass. LOCO identical in all 10 rows (checked by the script). 1,239 s wall, 7,339 CPU-s.

**Measured (gap = mean random-CV AUC − mean LOCO AUC, fold seeds 200-209, model seed 0).**

| model | old (seed 0) | mean | SD | min | max | seeds with gap ≤ 0.01 |
|---|---|---|---|---|---|---|
| logistic L2 | 0.071 | **0.071** | 0.004 | 0.065 | 0.076 | 0 of 10 |
| LightGBM | 0.115 | **0.111** | 0.003 | 0.106 | 0.117 | 0 of 10 |

Random-CV AUC over the 10 seeds: logistic 0.953 (0.947-0.958), LightGBM 0.950 (0.946-0.956). LOCO has no fold randomness: logistic 0.882, LightGBM 0.839 in every row.

**Side check (LightGBM model seed, fold seed 0).** LOCO AUC over model seeds 200-209: mean 0.842, SD **0.012**, range 0.816-0.856 (`results/seedsweep_lgbm_modelseed.csv`). The headline model seed 0 (0.839) sits inside that range. Note (review #16): `subsample=0.8` is inert without `subsample_freq`, so this is `colsample_bytree` randomness.

**Kill rule.** Not fired: both gaps stay above 0.01 in all 10 seeds (the rule needed ≤ 0.01 for both models in ≥ 6).

## Claim-check (2026-10-21)
Claim:      the raw cohort gap is not a fold-seed accident: logistic 0.071 ± 0.004, LightGBM 0.111 ± 0.003 over 10 fold seeds.
Measured:   exactly that, raw scale, one fixed LightGBM model seed, the five frozen cohorts.
Rival:      (a) LightGBM's own randomness: it moves LightGBM's LOCO 3-4× more than the fold seed moves its gap (SD 0.012 vs 0.003), so a LightGBM gap honestly reads "about 0.09-0.13", not "0.115"; the logistic model is deterministic, so its 0.065-0.076 range is the whole story. (b) The headline gap mixes estimands (review #14): the per-cohort paired gap is smaller (0.052 / 0.089 at seed 0) but positive in all 5 cohorts. (c) Exclusions flatter both schemes (README). (d) Probe choice reads test patients (review #3): control moved the gap by +0.007 / −0.003, under its 0.02 line.
Judge:      the pre-registered rule in the PLAN entry; the CSVs. Not judged by a human.

**What we now believe.** The gap is a stable property of these five cohorts under this pipeline, not a lucky split. "Reproducible" in CORE §6 can now be read as measured across 10 fold seeds (and, for LightGBM, with the model-seed spread stated beside it).

**Where.** `results/seedsweep_raw.csv` (20 rows), `results/seedsweep_raw_summary.csv`, `results/seedsweep_lgbm_modelseed.csv`, `results/run_seedsweep.json`.

## 2026-10-21 — RESULT: probe-choice control (review finding #3) — rejected with a number

**Run.** `python run_probe_control.py` (PLAN 1a0ddd7). Both rules give 463 × 12,040.

| model | gap, highest-mean probe (current) | gap, lowest probe id (reads no value) | change |
|---|---|---|---|
| logistic L2 | 0.071 | 0.077 | +0.007 |
| LightGBM | 0.115 | 0.112 | −0.003 |

The "mean" arm reproduces `summary_raw.csv` exactly. Both changes are under the pre-registered 0.02 line → the probe choice does not materially move the headline; finding #3 is **rejected** on this measurement (the code keeps the published convention, highest mean).

**Where.** `results/probecontrol_raw.csv`, `results/run_probecontrol.json`.

## 2026-10-21 — v1 README written (branch v1/release) + claim-check

**Done.** README gains "Result (v1)", "What the author checked vs what AI did" and "Labels"; 476 → 463 fixed; run commands for the sweep, the probe control and the tests. Every number re-checked against `summary_raw.csv`, `summary_withinsample.csv`, `seedsweep_raw.csv`, `seedsweep_lgbm_modelseed.csv`, `gapdiag_cohorts.csv`.

## Claim-check (2026-10-21, before v1 counts as done)
Claim:      of a random-CV AUC of ~0.95, ~0.88 survives a held-out cohort (raw, logistic; LightGBM 0.84) and ~0.90 within-sample; the raw gap is stable over 10 fold seeds and positive in every cohort.
Measured:   raw gaps over 10 fold seeds (both models); LightGBM LOCO over 10 model seeds; within-sample and per-cohort paired gaps at fold seed 0 only.
Rival:      LightGBM model randomness (stated, its gap range is INFERRED); averaging choice (paired per-cohort gap stated); exclusions flatter both (stated); probe choice (control < 0.02); labels (5/5 published-count checks, one AI-only, two with no pinned table).
Judge:      pre-registered rules + three outside AI code reviewers (could and did say HIGH; each HIGH checked on the data). No human judge yet: README says so. Narrowed during this check: within-sample = seed 0 only; LightGBM range marked inferred; "every finding fixed with a test" → "all but one with a test, or rejected/limit stated".

**Where.** `README.md` on `v1/release`. Nothing public, no remote (CORE §7). Next (the author): read the README, decide merge to master; then the outreach draft (v1 brief, "human check before public").

## 2026-10-21 — DECISION (the author): CORE §6 updated; merge waits for the author's own read

CORE §6 "reproducible UNVERIFIED" → "measured over 10 fold seeds" (commit 4e631bd, dated, marked as the author's decision, old wording kept). Merge of `v1/release` into master: NOT yet. The author reads README + the 21.10. entries first (choice A), then says "merge v1".

## 2026-10-21 — Pinning paper tables for GSE63990 and GSE42026 → GSE42026 verified, GSE63990 FAILS (condition 2 back to 4 of 5)

**Question.** The `source_of_expected` for these two cohorts said only "Multi-cohort signature papers describing …" (no table). Which published numbers do our counts actually match? (Open item from the v1 README; the author's go.) Checked by Claude (AI), not a human.

**GSE42026 — verified, now pinned.** GEO Series overall design (deposited by the Herberg 2013 authors, J Infect Dis 208:1664-8, PMID 23901082, [doi:10.1093/infdis/jit348](https://doi.org/10.1093/infdis/jit348)): after excluding co-infections, H1N1/09 19, RSV 22, bacterial 18, controls 33 → viral 41, bacterial 18, matching our parse 18 / 41 / 33. Cross-check: Tsalik 2016 Table 1 lists GSE42026 as 59 subjects (18 + 41; the two class columns of that row look swapped). Small note: the paper's own Results (34 RSV, 11 co-infected) imply 23 RSV, not 22. That is likely one array lost to QC, but not stated.

**GSE63990 — published SUBJECT counts disagree with our parse.** Tsalik et al., Sci Transl Med 2016;8:322ra11, PMID 26791949, [doi:10.1126/scitranslmed.aad6873](https://doi.org/10.1126/scitranslmed.aad6873), Table 1 and text: 273 ill subjects = **70 bacterial, 115 viral, 88 non-infectious**. Methods: "measured in two microarray batches with seven overlapping samples (GSE63990)"; GEO: 280 arrays. Our parse counts **arrays**: 73 / 117 / 90. The old `expected_counts` (73 / 117) matched the arrays, not a published table. That is what CLAUDE.md rule 2 warns about, and it went unnoticed since 2026-10-10.

**Which arrays are the re-runs (inferred, not proven).** GEO has no subject ids (fields: batch, infection_status, "sample_N"). The last 7 arrays, sample_274 … sample_280 (GSM1562133-GSM1562139, all batch 1), carry exactly the surplus: 3 bacterial + 2 viral + 2 non-infectious (73−70, 117−115, 90−88). In a cross-batch correlation check (top-5,000-variance probes, each batch centred), several of the strongest cross-batch matches involve these GSMs, but the pairing is not clean (r 0.80-0.93, with near-ties to other samples). So the identity of the 7 is a strong inference, not a fact.

**Why it matters.** The same patient twice in GSE63990 can sit on both sides of a random-CV split. That inflates random-CV and so the gap; leave-one-cohort-out is unaffected (both copies are in the same cohort). With ≤ 5 labelled duplicates in 463 rows the effect should be small, but it has to be measured, not assumed.

**What was changed now (record only, no headline number changed).** GSE63990 `expected_counts` → null (the cohort is UNVERIFIED, the CLAUDE.md option "mark the cohort unverified"), the published subject numbers and the finding are kept in a new `published_subjects_not_arrays` field. GSE42026 `source_of_expected` pinned. `--review`: 4 of 5 verified. **CORE §2 condition 2 is back to 4 of 5.** Fixing the parse (dropping the 7 candidates, then expected 70 / 115) changes the headline data, so it is the author's decision.

## 2026-10-21 — PLAN: GSE63990 duplicate-array sensitivity (NOT the headline)

**Question.** How much do the 7 candidate re-run arrays move the raw gap?
**Design.** `run_dedup_sensitivity.py`: raw headline setting (fold seed 0, model seed 0), both models, two arms: all 463 arrays vs the same minus the candidates that are labelled (expected 5: 3 bacterial + 2 viral; the 2 non-infectious are not in the matrix). Output `results/dedup_sensitivity_raw.csv`, `results/run_dedup_sensitivity.json`.
**Reading (fixed now).** A gap change ≤ 0.01 for both models = the duplicates do not move the headline; the v1 numbers stand with a stated caveat. > 0.01 for either model = the v1 headline must be re-run on de-duplicated data before it is quoted (the author's call).
**Cost.** CPU, ~3 min.

## 2026-10-21 — RESULT: GSE63990 duplicate-array sensitivity → over the line for LightGBM

**Run.** `python run_dedup_sensitivity.py` (PLAN 62a6819). 5 of the 7 candidates are in the labelled matrix (GSM1562134-1562138: 3 bacterial, 2 viral); the 2 non-infectious ones were never in it. 463 → 458 rows.

| model | gap, all arrays | gap, 5 candidates dropped | change | random-CV | LOCO |
|---|---|---|---|---|---|
| logistic L2 | 0.071 | 0.076 | +0.005 | 0.953 → 0.957 | 0.882 → 0.881 |
| LightGBM | 0.115 | 0.127 | **+0.012** | 0.954 → 0.954 | 0.839 → 0.827 |

**Reading (pre-registered).** LightGBM moves by more than 0.01 → the v1 headline must be re-run on de-duplicated data before it is quoted. The rule fired.

**What we now believe (inferred).** The duplicates did not inflate the gap; removing them makes it slightly larger. The LightGBM change comes from LOCO, not from random-CV. Removing training rows changes the models fitted for the other held-out cohorts, so this is mostly ordinary model variance: LightGBM's LOCO has SD 0.012 across model seeds, the same size as this change. So the change is within LightGBM's own noise, but the rule is the rule. The v1 numbers are **provisional** until the de-duplicated re-run; the README says so in its Result section.

**Decision for the author.** A (rec): accept the 7 candidates as the re-runs (label-surplus match + correlation evidence), add an exclusion rule for them, set GSE63990 expected 70 / 115 from Tsalik Table 1, and re-run raw, within-sample and the 10-seed sweep (~25 CPU-min); v1 numbers update. B: drop all of GSE63990's batch 1 overlap uncertainty by keeping only batch 0 (loses 124 arrays, worse). C: keep the 463 arrays and report the dedup sensitivity as a caveat only (cheapest, but GSE63990 stays unverified and condition 2 stays at 4/5).

**Where.** `results/dedup_sensitivity_raw.csv`, `results/run_dedup_sensitivity.json`, README (Result marked provisional, Labels table corrected).

## 2026-10-21 — PLAN: GSE63990 de-duplicated re-run becomes the v1 headline (the author's decision A)

**Decision (the author).** Accept sample_274-280 (GSM1562133-GSM1562139, batch 1) as the 7 re-run arrays and exclude them; the de-duplicated result is the MAIN number; the 280-array result stays in the README as one comparison row. Brief: `docs/brief-2026-10-21-gse63990-dedup.md`.

**What changes.** `config/cohorts.json` gets a per-cohort `exclude_samples` list (with its reason); `src/labels.py` marks those samples `excluded` before any class rule runs (same principle as rule 3: exclusions first). GSE63990 `expected_counts` → **bacterial 70, viral 115** (Tsalik 2016 Table 1). A test builds a series where the 7 ids come back and asserts the count check fails (MISMATCH) without the exclusion and passes with it; a second test pins the real config list.

**Expected.** `--review`: GSE63990 parses 70 / 115 / 95 excluded (88 non-infectious + 7 re-runs) → verified → **CORE §2 condition 2 back to 5 of 5**. If it does not parse 70 / 115, the cohort stays unverified and nothing below runs. Labelled matrix 463 → **458** rows (the 2 non-infectious re-runs were never in it). At fold seed 0 the raw headline must reproduce the `drop_7_candidates` arm of `results/dedup_sensitivity_raw.csv` (logistic gap 0.076, LightGBM 0.127); a mismatch = stop and find out why.

**Re-run (everything that feeds the README numbers).** `run_pipeline.py --review`, `--train` (raw), `--train --within-sample-normalize`, `--train --rank-normalize` (not reported, kept consistent because the gap diagnostics compare to it), `run_gap_diagnostics.py --parts e1` (the per-cohort paired gap), `run_seed_sweep.py` (10 fold seeds + LightGBM model-seed side check). Before overwriting, the 463-array files the README quotes are copied to `results/*_all280arrays.*` for the comparison row.

**Kill rule (unchanged).** The gap is a seed artefact if it is ≤ 0.01 for BOTH models in ≥ 6 of the 10 fold seeds.

**Cost.** CPU only, ~30-40 min, € 0. 14 GB RAM free at start.

## 2026-10-21 — RESULT: GSE63990 de-duplicated re-run → v1 headline updated, condition 2 at 5 of 5

**Question.** With the 7 GSE63990 re-run arrays excluded (the author's decision A; PLAN 4b003b3), do the labels match the paper, and what are the v1 numbers?

**Done.** `exclude_samples` in `config/cohorts.json` (GSM1562133-GSM1562139, with reason), applied in `src/labels.py` before any class rule; an id that is not in the series is refused. GSE63990 `expected_counts` 70 / 115 from Tsalik 2016 Table 1. `tests/test_dedup_exclusion.py` (4 tests): builds the 280-array case and asserts MISMATCH + training refused without the exclusion; on the old `labels.py` 2 of 4 go red. 112 tests pass. Commit be6e6e2. The 463-array files the README quoted are kept as `results/*_all280arrays.*`.

**Measured.** `--review`: **GSE63990 70 / 115 / 95 → verified; 5 of 5 cohorts verified → CORE §2 condition 2 back to 5 of 5** (GSE63990 and GSE6269/GSE42026 checks are AI-verified, not human). Matrix 463 → **458** × 12,040 (191 bacterial, 267 viral). Fold seed 0 reproduces the `drop_7_candidates` arm of `dedup_sensitivity_raw.csv` exactly (0.076 / 0.127), as the PLAN required. Seed-sweep repro gate 1.1e-16. 916 s wall for the sweep.

| | logistic L2 (458) | LightGBM (458) | logistic (463, comparison) | LightGBM (463, comparison) |
|---|---|---|---|---|
| random-CV / LOCO, raw, seed 0 | 0.957 / 0.881 | 0.954 / 0.827 | 0.953 / 0.882 | 0.954 / 0.839 |
| gap, seed 0 | 0.076 | 0.127 | 0.071 | 0.115 |
| gap, 10 fold seeds (mean ± SD, range) | 0.073 ± 0.004 (0.065-0.076) | 0.123 ± 0.006 (0.113-0.133) | 0.071 ± 0.004 | 0.111 ± 0.003 |
| within-sample random-CV / LOCO (seed 0) | 0.963 / 0.895 | 0.958 / 0.888 | 0.959 / 0.900 | 0.954 / 0.887 |
| per-cohort paired gap (seed 0), smallest cohort | 0.054, min 0.033 | 0.090, min 0.011 | 0.052 | 0.089 |
| LightGBM LOCO over model seeds 200-209 | — | 0.827 ± **0.023** (0.780-0.861) | — | 0.842 ± 0.012 (0.816-0.856) |

**Kill rule.** Not fired: 0 of 10 seeds with gap ≤ 0.01 for either model.

**Surprise.** LightGBM's own model-seed spread doubled (SD 0.012 → 0.023; one seed at 0.780). Removing 5 rows should not do that by itself; it says the LightGBM held-out number is fragile at this sample size. The README now reads the LightGBM gap as "clearly larger than logistic, size uncertain" (inferred range 0.09-0.17).

## Claim-check (2026-10-21, README after the de-duplicated re-run)
Claim:      of ~0.95 random-CV AUC, ~0.88 survives a held-out cohort (logistic; LightGBM 0.83), ~0.89 within-sample; the gap is stable over 10 fold seeds and positive in every cohort; removing the re-runs makes the gap slightly larger.
Measured:   exactly those numbers on 458 arrays (raw over 10 fold seeds; within-sample and per-cohort at fold seed 0 only); the 463-array numbers from the same code before the exclusion.
Rival:      (a) the 7 excluded arrays are not the true re-runs → some duplicates stay and 5 real patients are lost; bounded: the two arms differ by ≤ 0.012 in gap, and the labels now match the paper. (b) LightGBM model randomness (SD 0.023) is larger than the dedup change (+0.012): stated, LightGBM range marked inferred. (c) Exclusions flatter both schemes (stated).
Judge:      the pre-registered expectation (fold seed 0 must reproduce 0.076 / 0.127; it did, it could have failed), the kill rule, the published Table 1 count (could and did say MISMATCH before). No human judge yet. Narrowed: "one per patient" → "the 7 re-run arrays removed"; "larger" now carries "inside LightGBM's model-seed noise".

**What we now believe.** The v1 numbers are no longer provisional: the de-duplicated run is the headline, and the duplicates were not inflating the gap. The LightGBM number is the soft one.

**Where.** `results/summary_raw.csv`, `summary_withinsample.csv`, `seedsweep_raw*.csv`, `seedsweep_lgbm_modelseed.csv`, `gapdiag_cohorts.csv`, `label_review.md`; comparison `results/*_all280arrays.*`; README on `v1/release`. Next (the author): read README + these entries, say "merge v1". Still open: the identity of the 7 re-runs is inferred (tasks.md).

## 2026-10-21 — PLAN: v1.1, FROZEN published signatures scored both ways (the author's GO)

**Brief.** `docs/brief-2026-10-21-v1.1.md`; plan step "v1.1" in `docs/plan-2026-10-20-v1.md`. Branch `v1.1` off master 8f01501, local only.

**Question.** Does the random-CV vs leave-one-cohort-out (LOCO) gap also hit signatures whose genes were chosen OUTSIDE our data? Two published scores, frozen: Herberg 2016 2-transcript DRS (FAM89A − IFI44L) and Sweeney 2016 7-gene bacterial/viral metascore (geomean IFI27, JUP, LAX1 − geomean HK3, TNIP1, GPAA1, CTSB). Gene lists, directions, rules and sources pinned in `config/signatures.json` (AI-verified, not human).

**Two facts found BEFORE any AUC (they change the question).**
1. **Sweeney is not outside our data.** Its Table 1A lists GSE6269 GPL96 (73 / 18, exactly our arms), GSE40396 and GSE42026 as DISCOVERY cohorts. Only GSE60244 and GSE63990 are independent of its gene choice. Herberg used GSE42026, GSE40396, GSE6269, GSE60244 only as validation sets; GSE63990 is new to both.
2. **FAM89A has no probe on GPL571 (GSE63990) or GPL96 (GSE6269)** (annotation search for FAM89A / old name C1orf153 / Gene ID 375061: 0 hits). Herberg can be scored as published on 3 cohorts only (GSE42026, GSE60244, GSE40396). A missing gene is reported, never replaced. Sweeney: all 7 genes on all 5 platforms. GSE60244 has 10 values at 0 in IFI27/JUP (raw values clipped at 0 before log2), which makes a literal geometric mean 0 for those patients.

**What a "frozen" score means here.** No fit at all: each patient's score is computed from that patient's own row, no labels, no other patient, fixed sign (higher = bacterial; never flipped to make AUC > 0.5). "Random CV" = the same 5 stratified folds as the pipeline over all scorable patients (pooled across cohorts), AUC per test fold, mean; "LOCO" = AUC inside each cohort, unweighted mean. The gap is computed the same way as the v1 headline. Seeds: 0 plus 200-209.

**Expectation (written now, may be wrong).** With nothing trained, the v1 gap mechanism (learning local quirks of a cohort that is also in the test fold) cannot act. What is left is the pooling penalty: different cohorts put the same score on different offsets, which hurts a pooled random-CV fold but not a within-cohort AUC. So I expect the frozen gap to be near 0 or NEGATIVE. The informative comparison is then the frozen held-out AUC next to our trained models' held-out AUC on the same cohorts.

**Measured (outputs).** `results/frozen_signatures_cohorts.csv` (per signature × variant × cohort: n, AUC), `results/frozen_signatures_summary.csv` (random-CV mean ± SD over seeds 200-209, seed-0 value, LOCO mean, gap; also LOCO split into "independent of gene choice" vs "discovery" cohorts), `results/run_frozen_signatures.json`. Comparison column: our logistic/LightGBM LOCO AUC on the same cohorts from `results/folds_raw.csv`.

**Variants (pre-registered, all reported).** main = log2 values, highest-mean probe per gene (pipeline default); `linear` = same on 2^x (scale ambiguity in both papers; also removes the zero problem); `firstprobe` = lowest probe id per gene (the probe-choice control, reads no expression value).

**Rival explanations + controls.**
- (a) *Missing probes / wrong probe* make a signature look bad on a platform: coverage table above; Herberg restricted to its 3 covered cohorts; `firstprobe` variant.
- (b) *Gene choice saw the cohort* (Sweeney discovery overlap) flatters its held-out AUC: LOCO reported separately for independent vs discovery cohorts.
- (c) *Scale reading* (log vs linear, zeros): `linear` variant.
- (d) *Implementation bug* (wrong sign, wrong gene, a silent fit): tests that go red (wrong gene list, a scorer that reads other patients, a flipped AUC, a silently dropped gene).

**Kill rule (positive control).** Each paper reports its own score working on some of our cohorts: Sweeney summary AUC 0.97 over its discovery cohorts; Herberg AUC 0.892-0.966 on four published sets (eTable 6, incl. Hu = GSE40396, Ramilo = GSE6269, Suarez = GSE60244). If, in the main variant, Sweeney's AUC is < 0.75 on any of its 3 discovery cohorts, or Herberg's is < 0.75 on GSE40396 or GSE60244, I treat it as an implementation or mapping bug: STOP, no numbers into the README, find the cause first.

**Cost.** CPU, seconds (no training). Free RAM 16.2 GB at start.

## 2026-10-21 — RESULT: v1.1 frozen signatures → no random-CV bonus without fitting; held-out about equal to our logistic model

**Run.** `python run_frozen_signatures.py` (PLAN c193407, code 472f895). 122 tests pass (10 new; 4 deliberate mutations of `src/frozen_signatures.py` each turned at least one red: no published-list check, silently dropped gene, auto-flipped AUC, disabled per-sample guard). Seconds of CPU, 14.6 GB RAM free.

**Measured (main variant; seed 0 for per-cohort, gap over fold seeds 200-209).**

| score | cohorts (n) | random-CV | LOCO | gap, 10 seeds | LOCO, cohorts independent of its genes | our logistic / LightGBM LOCO, same cohorts |
|---|---|---|---|---|---|---|
| Herberg DRS | 3 (182) | 0.876 | 0.882 | −0.004 ± 0.005 (−0.014 … +0.004) | 0.882 (all 3) | 0.846 / 0.847 |
| Sweeney 7-gene | 5 (458) | 0.837 | 0.912 | −0.075 ± 0.002 | 0.901 (GSE63990 0.888, GSE60244 0.914) | 0.881 / 0.827 (independent 2: 0.885 / 0.851) |

Per cohort (Sweeney): GSE6269 0.972, GSE42026 0.940, GSE40396 0.847 (its 3 discovery cohorts); Herberg: GSE42026 0.866, GSE60244 0.889, GSE40396 0.892. Variants: `linear` gaps −0.047 / −0.095, `firstprobe` +0.001 / −0.086; held-out means move ≤ 0.016. Over all variants × seeds the gap lies in −0.100 … +0.010.

**Kill rule (positive control).** Not fired: lowest relevant AUC Sweeney GSE40396 0.847 (line 0.75). Herberg 0.889 / 0.892 on GSE60244 / GSE40396 sit at the low end of its own reported 0.892-0.966 (eTable 6, not read; range from the main text).

**Expectation check.** Pre-written "gap near 0 or NEGATIVE": held. It could have failed (a positive gap would have meant something other than fitting inflates random-CV).

**Blind judge.** A fresh AI agent, barred from the v1.1 code/config/results/LABBOOK/README, re-implemented both scores from the published definitions and matched all 10 per-cohort AUCs, both LOCO means and both pooled random-CV values to 3 decimals (same model family; checks the implementation, not the reading of the papers).

## Claim-check (2026-10-21, v1.1 README sentences)
Claim:      frozen published scores get no random-CV bonus (gap ≤ 0.01, mostly negative), and held out they do about as well as our logistic model and better than LightGBM on the same cohorts.
Measured:   pooled 5-fold mean AUC minus unweighted per-cohort mean AUC, 11 fold seeds × 3 variants; held-out means over 3 (Herberg) / 5 or 2 (Sweeney) cohorts at fold seed 0, no CIs; comparison against v1 `folds_raw.csv` (model seed 0).
Rival:      (a) "no bonus" is partly true by construction (nothing trained) → sentence narrowed to "expected; fits, does not prove, the fitting reading". (b) Sweeney's genes were chosen on 3 of our cohorts → independent-only number reported (0.901). (c) small n (GSE40396: 8 bacterial) and LightGBM seed SD 0.023 → "about the same", not "better"; no CI computed. (d) missing FAM89A → Herberg on 3 cohorts, stated. (e) scale/probe readings → 2 variants, same direction.
Judge:      pre-registered kill rule + pre-written expectation (both could have said no) and a blind re-implementation (could have disagreed). No human judge; the supplements (eTables 3/6, fig. S3) were not read.

**What we now believe.** The v1 gap is a property of training on the data, not of host-response genes as such: two gene lists chosen elsewhere, used frozen, lose nothing to a held-out cohort (by construction) and land at about the held-out AUC our trained models reach (inferred from small means). Pooling cohorts *hurts* a frozen score (Sweeney −0.075), which is the platform-offset problem Sweeney's COCONUT method was built for.

**Where.** `config/signatures.json`, `src/frozen_signatures.py`, `run_frozen_signatures.py`, `tests/test_frozen_signatures.py`, `results/frozen_signatures_*.csv`, `results/run_frozen_signatures.json`, README section "v1.1". Branch `v1.1`, local, NOT merged. Next (the author): read README "v1.1", say "merge v1.1". Open: a human check of the gene lists against the supplements.

## 2026-10-21 — PLAN: v1.1 supplements check — pin gene lists + score rules to the supplements, 2 more variants (the author's decision)

**Brief.** `docs/brief-2026-10-21-v1.1-supplements.md`. Branch `v1.1`, local only. Open point it closes: the gene lists were read from the PMC main texts only.

**Read (before any code or run; Claude, AI, not a human).**
- Herberg supplement `https://pmc.ncbi.nlm.nih.gov/articles/instance/5997174/bin/NIHMS75938-supplement-Supplementary_files.pdf`. **eTable 3** (38-transcript list): IFI44L log2 FC bacterial vs viral **−3.79**, FAM89A **+1.21** (and +1.26 for its 2nd array id) → FAM89A bacterial-up, IFI44L viral-up: as pinned. **Supplementary Methods, eq. (1):** "the DRS was calculated by subtracting the log2 transformed expression value of IFI44L from the log2 transformed FAM89A expression value for every patient. … Higher scores indicate bacterial assignment" → the rule is **log2**; our main reading is right, and the `linear` variant becomes a pure sensitivity (kept, pre-registered). **eTable 6** (external datasets): GSE6269 **GPL570 subset only** (12 bacterial / 10 viral; probes 204439_at, 226448_at) — NOT our GPL96 arrays; GSE40396 8 / 22, **AUC 89.2**; GSE60244 15 / 64 (we have 22 / 71), AUC 89.8. **GSE42026 appears nowhere in the 2016 paper** (main text + supplement). Two config statements were wrong: "GSE42026 ref 13 validation" and "GSE6269 validation" as if it were our arrays. Neither changes a number (the `independent` list means "not used to choose the genes", and Herberg is not scored on GSE6269), but the notes get fixed.
- Sweeney supplement `https://pmc.ncbi.nlm.nih.gov/articles/instance/5348917/bin/NIHMS850538-supplement-supplementary_figures_and_tables.pdf`. **Fig. S3** forest plots name exactly the 7 genes; **Table S1** summary effect sizes (viral vs bacterial, log2, Hedges' g): IFI27 +2.299, JUP +1.008, LAX1 +1.114, HK3 −1.109, TNIP1 −1.023, GPAA1 −1.275, CTSB −1.080 → as pinned. Fig. S3 discovery cohorts include GSE6269 gpl96, GSE40396, GSE42026 (as pinned); Table S2 lists GSE60244 and GSE63990 as direct validation (as pinned). Fig. S1 counts: GSE40396 8 / 22, GSE42026 18 / 41, GSE60244 22 / 71, GSE63990 70 / 115 = ours.
- MetaIntegrator `calculateScore` (CRAN v2.1.3, 2020, `https://github.com/cran/MetaIntegrator/blob/master/R/calculateScore.R`): geomean(pos) − geomean(neg), **no ratio-of-genes factor**; if a dataset's minimum is < 0 the whole dataset is shifted by |min| + 1; then `scale()` = **z-score inside the dataset**. The package postdates the 2016 paper; it is the authors' public implementation, not proof of what the paper ran.

**Expectation.** The gene lists, directions and our main rule are confirmed → the **main variant numbers must reproduce exactly** (same code path); any change = a bug, stop.

**New variants (pre-registered now, reported next to the 3 old ones, never replacing main).**
4. `ratio_one_term` — the other reading of Sweeney's "minus the geometric mean of the bacterial genes, times the ratio of the number of genes in each set": the factor sits on the bacterial term only, ratio in the order the text names the sets, n_viral / n_bacterial = 3/4. Our sign: 0.75 · geomean(bacterial) − geomean(viral). (The reverse order, 4/3, is not run; said here so it is not silently skipped.) Herberg has 1 : 1 genes → factor 1 → identical to main by construction (reported as such).
5. `metaintegrator` — main score, then the `calculateScore` steps: shift a cohort if its min < 0, z-score the score inside each cohort. This reads the whole cohort (like `--rank-normalize`), so it is NOT a per-patient frozen score; it is the authors' software reading. Prediction: LOCO identical to main (a monotone per-cohort transform cannot change a within-cohort AUC); pooled random-CV goes UP because cohort offsets are removed → gap moves toward 0 or above. If this happens it supports, not overturns, the "pooling penalty = platform offsets" reading.

**Headline (fixed now, for the merge rule).** The README v1.1 table's **main-variant** numbers: random-CV, LOCO, 10-seed gap, independent-cohort LOCO, for both scores. Merge rule (brief): none moves by > 0.01 → merge v1.1 into master locally; otherwise no merge, options to the author. Variants 4-5 are sensitivity rows, not headline.

**Tests (must go red on broken code).** A test transcribes the supplement effect signs (eTable 3, Table S1) independently of `PUBLISHED`, derives the up/down lists and asserts they equal `PUBLISHED` and the config; it also builds a wrong list (one gene swapped) and asserts the comparison fires. Plus: the ratio variant weights the bacterial term by 3/4; the MetaIntegrator step gives per-cohort mean 0 / SD 1, shifts a negative cohort, and leaves every within-cohort AUC unchanged. Mutation check: edit `PUBLISHED` → the supplement test must go red.

**Kill rule.** Unchanged (main variant, 0.75 line). Plus: main reproduces to 1e-12 or stop.

**Cost.** CPU, seconds. € 0.

## 2026-10-21 — RESULT: v1.1 supplements check → every gene and sign confirmed, main numbers unchanged (diff 0), merge rule allows the merge

**Run.** `python run_frozen_signatures.py` (PLAN 661ea85). 127 tests pass (5 new). Mutation check: a gene swapped in `PUBLISHED` + config together (which the old config test cannot see) → the new supplement test goes red; ratio inverted (4/3) → red; z-score across instead of within cohorts → red. Seconds of CPU.

**Measured.** Main, `linear`, `firstprobe`: all 6 summary rows and 24 per-cohort AUCs equal the c193407 run, **max abs difference 0**. → No headline number moves; the brief's merge condition (none > 0.01) holds.

| score | variant | random-CV (seed 0) | LOCO | gap, 10 seeds (range) | LOCO independent |
|---|---|---|---|---|---|
| Herberg | main (= ratio_one_term, 1 : 1 genes) | 0.876 | 0.882 | −0.004 ± 0.005 (−0.014 … +0.004) | 0.882 |
| Herberg | metaintegrator | 0.883 | 0.882 | −0.003 ± 0.006 (−0.014 … +0.005) | 0.882 |
| Sweeney | main | 0.837 | 0.912 | −0.075 ± 0.002 | 0.901 |
| Sweeney | ratio_one_term | 0.819 | 0.910 | −0.088 ± 0.002 | 0.898 |
| Sweeney | metaintegrator | 0.859 | 0.912 | −0.055 ± 0.002 (−0.059 … −0.051) | 0.901 |

Per cohort, `ratio_one_term` moves Sweeney by ≤ 0.008 (GSE42026 0.940 → 0.932); `metaintegrator` per-cohort AUCs are identical to main, as predicted. Positive control from the supplement: GSE40396 (same 8 / 22 counts as Herberg eTable 6) our 0.892 vs their 89.2 %; GSE60244 our 0.889 (22 / 71) vs their 89.8 % (15 / 64, a different subset).

**Expectation check.** "Main reproduces exactly": held. "metaintegrator: LOCO identical, random-CV up, gap toward 0 or above": LOCO identical and gap toward 0 held; it did NOT reach 0 (−0.055). Inferred reason, not tested: z-scoring centres each cohort on its own mean, and the cohorts' bacterial shares differ (GSE6269 80 %, GSE63990 38 %), so the class-specific offsets survive.

## Claim-check (2026-10-21, changed README v1.1 sentences)
Claim:      the supplements confirm every gene and sign, the numbers did not change; Herberg chose its genes on none of our cohorts; the ratio reading moves held-out ≤ 0.016; MetaIntegrator's z-scoring shrinks Sweeney's pooling penalty to −0.055 without touching held-out; GSE40396 matches Herberg's own 89.2 %.
Measured:   the supplement PDFs as text (Firecrawl PDF parse; eTable 3, eq. 1, eTable 6, fig. S3, Table S1), the MetaIntegrator v2.1.3 source, and the re-run with diff 0 on the old variants.
Rival:      (a) a PDF-parse error in a table → the FAM89A/IFI44L signs are double-sourced (eTable 3 and eq. 1 text), Sweeney's by fig. S3 and Table S1. (b) "same 8 / 22" ≠ same patients → narrowed to "same counts". (c) "GSE42026 not in the paper" rests on a text search of the PMC main text + supplement → narrowed to "not named as a dataset"; patient overlap stays "not stated". (d) the class-mix explanation for the remaining −0.055 is a guess → marked inferred. (e) MetaIntegrator (2020) may not be what the 2016 paper ran → stated in config and here.
Judge:      the pre-registered "main must reproduce to 1e-12" (could have failed), the mutation-checked tests, Herberg's own published AUC on the same-count cohort (could have disagreed). No human read the supplements; still AI-verified.

**What we now believe.** The v1.1 gene lists and rules were right; the open "supplements not read" point is closed (by AI, not a human). The authors' own software would show a smaller pooling penalty for Sweeney (−0.055 instead of −0.075), still negative, still no random-CV bonus.

**Where.** `config/signatures.json` (supplement links + table/figure per field, Herberg overlap note corrected), `src/frozen_signatures.py`, `run_frozen_signatures.py`, `tests/test_frozen_signatures.py`, `results/frozen_signatures_*.csv`, README "v1.1". Next per brief: merge v1.1 into master locally (the author's OK under the rule), then outreach drafts.

## 2026-10-21 — BLIND CHECK of condition 4

**Question.** Does the 2026-10-20 condition-4 entry (merge 74cbb73) say what `results/summary_raw.csv` and `results/summary_withinsample.csv` say? Second, independent re-read after the 2026-10-20 one.

**Who checked.** A fresh Claude (AI) worker agent, not a human; it did not write the entry. Same model family as the author, so this is an independent *read*, not an independent judge of the science. Nothing was re-run (CPU-free, read-only).

**Order.** The two CSVs were read first; the numbers below were written down before the LABBOOK entry was opened.

**Blind read of the CSVs (mean AUC over 5 folds / 5 held-out cohorts).**
| scale | model | random_cv | leave_one_cohort_out | gap |
|---|---|---|---|---|
| raw | logistic_l2 | 0.953 | 0.882 | 0.071 |
| raw | lightgbm | 0.954 | 0.839 | 0.115 |
| within-sample | logistic_l2 | 0.959 | 0.900 | 0.059 |
| within-sample | lightgbm | 0.954 | 0.887 | 0.067 |

**Criteria.**
| criterion | verdict | detail |
|---|---|---|
| Raw AUCs 0.953/0.882 (logistic), 0.954/0.839 (LightGBM) | PASS | all four quoted in the entry match to 3 decimals |
| Within-sample AUCs 0.959/0.900, 0.954/0.887 | PASS with note | held-out 0.900 and 0.887 are quoted and match. The within-sample random-CV values (0.959 logistic, 0.954 LightGBM) are NOT quoted in the entry: it quotes only the raw 0.953/0.954 as "random-CV AUC of 0.95". Nothing is wrong, but 5 of the 6 values the card lists are never written there as such (0.954 does appear, for raw LightGBM, and equals the within-sample LightGBM value). |
| Headline 0.88 raw and 0.90 within-sample "survive from 0.95" | PASS | 0.882 -> 0.88, 0.900 -> 0.90 (logistic held-out); 0.953 -> 0.95 |
| Gaps quoted (raw 0.071/0.115, within-sample 0.059/0.067) | PASS | match the `gap` column to 3 decimals |
| raw / within-sample labels not swapped | PASS | raw held-out 0.882/0.839 (raw file), within-sample 0.900/0.887 (within-sample file) |
| No rank-normalised number quoted | PASS | the entry only says rank-normalised numbers are NOT used (condition 3); no `summary_ranknorm.csv` value appears |
| Claim-check names the single-run limitation | PASS | "no seed sweep is on file ... single-run numbers" and rival (b) "single-run noise — not ruled out" |
| Extra figures in the entry | PASS | per-cohort spread raw logistic 0.82-0.94 and LightGBM 0.70-0.89 match min_auc/max_auc (0.824-0.941, 0.698-0.894); "7 to 12 points" matches gaps 0.071-0.115 |

**Mismatches verbatim.** None in any number. One omission, same as the first blind re-read's caveat: the entry's wording "gives back about one to five of those points" compares gaps (0.071 -> 0.059, 0.115 -> 0.067), which are measured against the within-sample random-CV AUCs; a reader subtracting from the raw 0.953 would get 0.053, not 0.059.

**Claim-check.** Claim: the entry's numbers and labels equal the CSVs. Measured: 18 numbers compared, 0 mismatches. Rival: the reader copied obvious cells; for a transcription claim that is the whole test. Judge: the CSVs themselves. Not a verdict on the science (single run for these two files; the separate 10-seed sweep result is in its own entry above; GSE6269 labels AI-verified, not human).

**Verdict: PASS.** No number corrected, entry untouched.

**Where.** This entry; branch `a side branch` (from master, not merged; no remote, CORE §7).

## 2026-10-22 — PLAN: within-sample 10-seed sweep + restate the condition-4 sentence at 458 arrays (approved by the author)

**Question.** The within-sample half of the §1 answer sentence rests on fold seed 0 only. Does the within-sample random-CV minus leave-one-cohort-out gap hold over 10 fold seeds?

**Setup.** `run_seed_sweep.py --within-sample`: `rank_normalize_within_sample(X)` once on the 458 x 12,040 matrix (one patient at a time, nothing reads another patient, so applying it before the folds leaks nothing), then the same `sweep_rows` as the raw sweep. Fold seeds 200-209, model seed 0. Repro gate: fold seed 0 must reproduce `results/summary_withinsample.csv` to 1e-12, else stop. LightGBM model-seed side check NOT run (its range is inferred from the raw side check, SD 0.023). CPU only, expected ~15-20 min, EUR 0.

**Claim.** Under leak-free within-sample normalisation about 0.89 of a ~0.96 random-CV AUC survives a held-out cohort, and that is not a fold-seed accident.

**Pre-registered reading.** If the 10-seed within-sample mean gap differs from the seed-0 gap by more than 0.01 for either model, the README/§1 within-sample number is wrong and must be replaced. A within-sample gap <= 0.01 is NOT a kill (CORE §3: it is an answer). Rival: seed 0 was a lucky/unlucky split, or LightGBM's own randomness.

**Where.** Branch `within-sample-seedsweep` (no merge, no remote, CORE §7).

## 2026-10-22 — RESULT: within-sample 10-seed sweep → the within-sample gap is not a seed accident; condition-4 sentence restated at 458 arrays

**Question.** The within-sample half of the §1 sentence rested on fold seed 0 only. Does its gap survive 10 fold seeds? (PLAN above.)

**Design.** Same sweep as the raw one, but each patient's genes are ranked against each other first (`rank_normalize_within_sample`, one patient at a time, no other patient read). Random-CV fold seeds 200-209, model seed 0, 458 arrays x 12,040 genes, 5 cohorts. Repro gate: fold seed 0 vs `summary_withinsample.csv`, max diff 1.1e-16 (limit 1e-12) — passed. Leave-one-cohort-out (LOCO) has no fold randomness, so it is identical in all 10 rows (checked). 584 s wall, CPU only.

**Measured** (`results/seedsweep_withinsample.csv`, `..._summary.csv`, `run_seedsweep_withinsample.json`).
| model | random-CV AUC (10 seeds) | LOCO AUC | gap mean ± SD | gap range | gap at seed 0 | 10-seed mean minus seed 0 |
|---|---|---|---|---|---|---|
| logistic L2 | 0.960 ± 0.003 | 0.895 | 0.065 ± 0.003 | 0.061-0.068 | 0.068 | -0.004 |
| LightGBM | 0.954 ± 0.004 | 0.888 | 0.067 ± 0.004 | 0.062-0.074 | 0.070 | -0.003 |

In words: moving the random split changes the random-CV number by only a few thousandths; the held-out-cohort number does not move with the split at all. Seed 0 was a slightly high draw (about 0.003-0.004 AUC), not a lucky one.

**Pre-registered reading.** Kill = 10-seed mean differs from seed 0 by more than 0.01 for either model. Differences are 0.004 and 0.003 → **not fired**; the within-sample number in the README/§1 stands. Gap <= 0.01 count: 0 of 10 seeds for both models (also not a kill by CORE §3). Not swept: LightGBM's own model seed on the within-sample scale. Its range is INFERRED from the raw side check (LOCO SD 0.023 at 458 arrays, `seedsweep_lgbm_modelseed.csv`); the within-sample LightGBM LOCO 0.888 is therefore one model-seed draw.

**Condition-4 sentence, restated (supersedes the 2026-10-20 one).** Of a random-CV AUC of about 0.96 (logistic L2 0.957 raw / 0.960 within-sample, LightGBM 0.954 / 0.954), **0.88 survives a held-out cohort on the raw scale** (logistic LOCO 0.881, gap 0.073 ± 0.004 over 10 fold seeds; LightGBM 0.827, gap 0.123 ± 0.006, model-seed range inferred wide) and **0.89 survives under within-sample normalisation** (logistic 0.895, gap 0.065 ± 0.003; LightGBM 0.888, gap 0.067 ± 0.004). Sources: raw `results/summary_raw.csv` (seed 0) + `results/seedsweep_raw_summary.csv` + `results/seedsweep_raw.csv` (random-CV means); within-sample `results/summary_withinsample.csv` (seed 0) + `results/seedsweep_withinsample_summary.csv` + `results/seedsweep_withinsample.csv`. All 458 arrays, commit of this entry. Plain words: roughly 6.5 to 12 AUC points of the headline are hospital-specific; the bedside-legal per-patient transform gives back about 1 point for logistic and about 6 for LightGBM (whose raw number is the fragile one).

## Claim-check (2026-10-22)
Claim:      under leak-free within-sample normalisation about 0.89 of a ~0.96 random-CV AUC survives a held-out cohort, and not by fold-seed accident.
Measured:   random-CV minus LOCO AUC, within-sample scale, 10 fold seeds x 2 models, 458 arrays; repro gate 1.1e-16.
Rival:      (a) seed 0 was an odd split — ruled out for fold seeds (mean differs by 0.003-0.004). (b) LightGBM model randomness — NOT tested here; inferred from the raw side check (SD 0.023). (c) only 5 held-out cohorts, so the LOCO value itself rests on 5 numbers, unchanged by this sweep. (d) labels: 3 of 5 cohort checks are AI-only.
Judge:      none outside the lab for this entry; the numbers are machine-written CSVs, the sentence is self-written. The restated sentence still needs a blind re-read against the CSVs (as the earlier blind check did for the old one).

**Next (not part of this card).** Carry the within-sample sweep numbers + this entry into the public copy of this repo (this entry is that carried copy); then the author's own check of 1-2 published count tables; then outreach to other researchers.

**Where.** Branch `within-sample-seedsweep` in the original repo (no merge, no remote, CORE §7). Test: `tests/test_seed_sweep.py::test_within_sample_flag_changes_numbers_and_paths` (goes red when the flag falls back to raw — confirmed by mutation).

## 2026-10-23 — BLIND RE-READ of the restated condition-4 sentence (carry into this copy)

**Who.** A fresh Claude (AI) agent that saw only the sentence and six CSVs (`summary_raw`, `seedsweep_raw_summary`, `seedsweep_raw`, `summary_withinsample`, `seedsweep_withinsample_summary`, `seedsweep_withinsample`), not the reasoning above. Same model family as the author: an independent read, not an independent judge of the science.

**Result.** 21 numbers checked. 19 match: all LOCO AUCs (0.881, 0.827, 0.895, 0.888), all 8 gap means and SDs, n = 10 seeds, model seed 0 everywhere, and "6.5 to 12 points" (gaps 0.065-0.123; 0.123 rounds to 12). **2 mismatches, in the random-CV figures of the 2026-10-22 RESULT sentence above:** it quotes logistic raw 0.957 and LightGBM raw 0.954, which are fold-seed-0 values, while the gaps beside them are 10-seed means (10-seed means: 0.953 and 0.950). The within-sample figures (0.960, 0.954) are 10-seed means. So 0.957 - 0.881 does not equal the quoted 0.073 (seed-0 gap is 0.076). Not checkable from the CSVs: the LightGBM model-seed range and "legal at the bedside".

**What was done.** The 2026-10-22 RESULT sentence is left as written (it is a diary record). The README one-sentence summary in this copy uses the 10-seed means (0.953 / 0.960 and 0.950 / 0.954, "about 0.95 to 0.96"). The headline conclusions (about 0.88 and 0.89 survive) do not change.


## 2026-10-23 — Checked by the author: published class counts for GSE40396 and GSE6269

**Question.** Do the published class counts the label check relies on (`config/cohorts.json` `expected_counts`) really say what we wrote? Until today both were read by Claude (AI) only.

**Done.** Claude opened each paper at the exact passage; the author read the numbers and compared them with ours.
- Hu et al., PNAS 2013 (PMC3732941), Results, first paragraph: 30 febrile children = 8 HHV-6 + 8 adenovirus + 6 enterovirus + 8 bacterial. Ours: viral 22 / bacterial 8.
- Ramilo et al., Blood 2007 (PMC1801073), Table 1 footnote: E. coli 29 + S. aureus 31 + S. pneumoniae 13 = 73 bacterial, influenza A 18. Ours: bacterial 73 / viral 18. The paper's larger totals elsewhere cover all three platforms; we use the GPL96 Table 1 set only.

**Measured.**
- checked by the author (2026-10-23): Hu 2013 PNAS counts 22 viral / 8 bacterial match.
- checked by the author (2026-10-23): Ramilo 2007 Blood Table 1 counts 73 bacterial / 18 viral match.

**Belief.** 2 of the 5 cohort count checks now have a human judge (GSE40396, GSE6269). GSE63990 and GSE42026 are still checked by AI only; GSE60244 was accepted by the author on 2026-10-19 from the AI comparison, not re-read in the paper. No code, config or number changed.

**Where.** README: *Labels are checked*, *Honest limits*, the provenance table and the Labels table.

## 2026-10-23 — Outreach: first author reply confirms the Sweeney score reading

**Done.** Questions sent to the authors of the two frozen signatures and of two cohorts. One of the Sweeney 2016 authors replied the same day.

**Measured (what the reply says, paraphrased).**
- Sweeney score reading: the author confirms that our reading of the score is right. The constant ratio factor was added at a reviewer's request; their own analysis showed it does not matter. This matches our variant test (factor on one term moves held-out AUC 0.912 -> 0.910; on the whole difference it changes no AUC).
- On random CV: in the author's view, random CV leaks information; the classifier has very likely seen samples from the dataset left for validation, one reason classifiers look ready for translation and then fail. Note: this answers the trained-model gap (our headline). It does not directly answer why the FROZEN score lands lower on pooled random-CV (a frozen score cannot leak); that question stays open.
- Pointer: the 7-gene score is part of the 29-gene signature validated prospectively in 1,222 patients across 22 emergency departments (USA + Europe), now FDA-cleared: Liesenfeld et al., Nature Medicine 2025, https://doi.org/10.1038/s41591-025-03933-y (not yet read by us). _[DOI corrected 2026-10-25: this note first carried 10.1038/s41591-025-03956-5, a different Nature Medicine paper; see the 2026-10-24 entry.]_
- Publication: the author has no objection to us publishing this work, in whatever form and at whatever time suits us.

**Belief.** First outside human judge on this study: the Sweeney score is implemented as the authors intended. Labels and code are still unchecked by a human expert (two further author groups asked on 2026-10-23, no reply yet).

**Where.** Label CSVs sent with the questions: `reports/2026-10-23-outreach-labels/`.

## 2026-10-24 — Read Liesenfeld 2025 (TriVerity validation); cited in README Prior art

**Question.** What does the paper behind the FDA-cleared 29-gene test actually report, and how does it relate to this study?

**Done.** Read the open-access full text (PMC12705421). Correction: the DOI in the 2026-10-23 pointer entry (`10.1038/s41591-025-03956-5`) is a different Nature Medicine paper (an immune-dysregulation framework from the SUBSPACE consortium). The Liesenfeld paper is `10.1038/s41591-025-03933-y`; the earlier entry's DOI was corrected in place on 2026-10-25, with a bracketed note saying so.

**Measured (read in the paper, not run by us).**
- SEPSIS-SHIELD: prospective, non-interventional; 1,441 adults enrolled in 22 emergency departments (21 US, 1 Europe), March 2020 to May 2024; 1,222 with valid TriVerity results; 729 consensus-adjudicated (main population); 1,120 evaluable for the severity endpoint.
- Consensus population: Bacterial AUROC 0.83 (80% CI 0.81-0.85), Viral 0.91 (0.89-0.93), Severity 0.78 (0.75-0.81). Forced adjudication: 0.76 bacterial, 0.83 viral.
- The existing README row ("1,222 patients across 22 sites") agrees with the paper; no change to it.

**Belief.** TriVerity's published numbers come from a prospective multi-site validation of one locked assay against clinical adjudication. They are not comparable with our retrospective random-CV vs leave-one-cohort-out numbers, and the README says so.

**Where.** README, *Prior art*, paragraph under the table.

## 2026-10-25 — Repo published privately on GitHub

**Done.** This public-ready copy went to GitHub as a **private** repository under
the author's real-name account, with a **one-commit history** (the local
preparation history was not pushed). Added before the commit: a `LICENSE` file
(MIT), a License line at the end of the README, and this entry. The three
verbatim phrases from the 2026-10-23 author reply were paraphrased.

**Measured.** Nothing new. No number, result file or code path changed.

**Belief.** Unchanged.

**Where.** GitHub, private, read access granted per reviewer on request.
**Gate:** the repo stays private until the two remaining author groups have had
about two weeks to reply (~2026-11-06); the switch to public is the author's
decision alone.

## 2026-10-31 — README bvnGPS row corrected

README bvnGPS row corrected, mirrored to the public copy: the Prior-art row now states the held-out-cohort and random-split modes, the five shared GSE ids, and the two reasons it is still not the same measurement (two invocations plus a separate scoring script; not runnable as shipped). Source of truth: the 2026-10-30 scoop re-check entry above. No other README text touched; no code or results changed.

## 2026-10-30 — PLAN (the LightGBM model-seed check): does the 6-point within-sample give-back for LightGBM survive LightGBM's own model seed?

*(Written and committed BEFORE any code or run. Approved by the author 2026-10-30 15:18.)*

**Question.** The condition-4 sentence says the per-patient (within-sample) transform "gives back about 6 points" for LightGBM (leave-one-cohort-out AUC, LOCO: raw 0.827 -> within-sample 0.888, both at model seed 0). The raw side has a model-seed side check (`results/seedsweep_lgbm_modelseed.csv`, 458 arrays, mean 0.8265, SD 0.0226). The within-sample side does not: 0.888 is one model-seed draw (the within-sample 10-seed sweep said so). Is a 6-point move bigger than LightGBM's own randomness?

**Setup.** `rank_normalize_within_sample` once on 458 x 12,040; LightGBM LOCO only, fold seed 0, model seeds 200-209 (same as the raw side check). Repro gate: model seed 0 must reproduce `results/summary_withinsample.csv` LightGBM LOCO to 1e-9, else stop and write nothing. New option `--lgbm-modelseed` in `run_seed_sweep.py`, valid only with `--within-sample`; defaults still reproduce every existing file.

**Pre-registered readings (copied from the plan, applied exactly as written).**
- (a) If the 10-model-seed mean within-sample LightGBM LOCO differs from 0.888 by more than 0.01, the within-sample LightGBM number in the condition-4 sentence is replaced by mean +- SD.
- (b) If the mean paired difference (within-sample minus raw LOCO, same model seed) is smaller than 2 x the SD of those 10 paired differences, the phrase "gives back about 6 points for LightGBM" is withdrawn and restated as "not distinguishable from LightGBM's own randomness".
- Neither is a lab kill (CORE §3). Note: (a) is measured against the exact file value 0.8876018 (the plan's 0.888 is its rounding); (b) is a strict yardstick for n=10 (needs roughly t > 6.3), so the result may land near the line.

**Claim-check.** Claim: the per-patient transform recovers about 6 AUC points of LightGBM's held-out-cohort loss. Measured: LightGBM LOCO AUC, 10 model seeds x 2 scales, paired by seed, 458 arrays. Rival: the 0.827 -> 0.888 move is LightGBM model-seed noise (raw SD about 0.02), not the transform. Control: the committed raw-scale model-seed run (same seeds, same folds; only the transform differs). Judge outside us: none; the pre-registered thresholds above are the judge, fixed before the run. Limit: `subsample=0.8` is inert without `subsample_freq` (review #16), so model-seed randomness is `colsample_bytree` only.

**Cost.** CPU only, about 10-15 min wall, EUR 0, 0 GPU-hours. Both validation schemes: this run measures LOCO only; random-CV is unchanged from the within-sample 10-seed sweep.

## 2026-10-30 — RESULT (the LightGBM model-seed check): the 6-point give-back survives LightGBM's own model seed; neither pre-registered reading fires

**Question and why.** The condition-4 sentence says the per-patient (within-sample) transform gives back about 6 points of LightGBM's held-out-cohort loss (leave-one-cohort-out AUC, LOCO: raw 0.827 -> within-sample 0.888, both at model seed 0). The raw number had a model-seed side check (SD 0.023); the within-sample number was a single draw. If LightGBM's own randomness is about 0.02, a 6-point move could have been luck. This run gives the within-sample side its own 10 draws and pairs them with the raw ones.

**Design.** `python run_seed_sweep.py --within-sample --lgbm-modelseed`: rank-normalise each patient's own genes (458 arrays x 12,040 genes), then LightGBM LOCO at fold seed 0 for model seeds 200-209, the same seeds as `results/seedsweep_lgbm_modelseed.csv`. Only the transform differs between the two columns. Repro gate (model seed 0 vs `summary_withinsample.csv`): |diff| = 0.0 (tolerance 1e-9), passed. CPU only, 638 s wall, EUR 0, 0 GPU. Existing result files untouched (`git status`: only the two new files). Run commit 60f46b2. Suite: 134 tests green (6 new); 5 hand mutations (skip rank-norm, ignore model seed, wrong gate reference, flipped (b) comparison, overwrite the raw file) each turned a test red.

**Measured.** LOCO scheme only (random-CV unchanged from the within-sample 10-seed sweep; the gap table there stands).

| | mean | SD | min - max |
|---|---|---|---|
| raw, 10 model seeds | 0.8265 | 0.0226 | 0.780 - 0.861 |
| within-sample, 10 model seeds | 0.8909 | 0.0055 | 0.882 - 0.899 |
| paired difference (within-sample minus raw) | +0.0644 | 0.0243 | +0.027 - +0.112 |

Model-seed-0 values for comparison: raw 0.827, within-sample 0.8876.

In plain words: across ten different LightGBM random draws the within-sample number sits at 0.891, almost exactly where the single draw (0.888) was, and it is much steadier (SD 0.006) than on the raw scale (0.023). In every one of the 10 seeds the within-sample run beat the raw run with the same seed; the smallest gain was +0.027, the largest +0.112.

**Pre-registered readings, applied as written.**
- (a) |mean 0.8909 - 0.8876| = 0.0033, not above 0.01: does NOT fire. The within-sample LightGBM number in the condition-4 sentence stays.
- (b) mean paired difference 0.0644 vs 2 x SD of the pairs = 2 x 0.0243 = 0.0485: 0.0644 is not smaller, does NOT fire. "Gives back about 6 points for LightGBM" stays (the 10-seed mean is 6.4).
- Neither fired, so no README edit and no change to the condition-4 sentence. It now has a model-seed control on both sides.

**Descriptive extra (not a reading, not a rescue).** Standard error of the mean difference 0.0077, t about 8.4 on 9 df. (b) was a strict yardstick for n=10 and was cleared, with margin of about 1 SD of the pairs.

**Claim-check.** Claim: the transform recovers about 6 AUC points for LightGBM. Measured: 10 paired model seeds, 458 arrays, LOCO. Rival (model-seed noise): not supported; raw SD 0.023 is smaller than the 0.064 shift and the shift has the same sign in 10/10 seeds. Control: committed raw model-seed run, same seeds and folds. Judge: thresholds fixed before the run in the PLAN entry; no outside human judge. Limits: model-seed randomness is `colsample_bytree` only (`subsample` inert without `subsample_freq`, review #16); fold seed fixed at 0 and the same five cohorts, so this says nothing about other cohorts; within-sample scaling assumes the whole patient profile is available (see Known soft spots); rank-normalised numbers are not quoted outside the repo. The paired SD is driven mostly by the raw side's randomness.

**Implications / options.** (1) Leave everything as is: the claim holds. (2) Optional, only if the author wants it: quote "0.891 +- 0.006 (10 model seeds)" next to the single-seed 0.888 where it is cited. Public copy needs nothing since nothing changed.

**Where.** A worker branch (local only, no remote, CORE §7; merge is the author's). Files: `run_seed_sweep.py` (`--lgbm-modelseed`, `lgbm_loco`, `compare_modelseed`), `tests/test_seed_sweep.py`, `results/seedsweep_lgbm_modelseed_withinsample.csv`, `results/run_seedsweep_lgbm_modelseed_withinsample.json`.
## 2026-10-30 — pre-public scoop re-check (CORE §3 kill 2), the pre-public scoop re-check

**Question asked.** Before the public copy goes out: is there now a public, runnable repository that (1) trains on some host-response cohorts, (2) tests on a held-out cohort, and (3) prints BOTH the random-CV AUC and the held-out-cohort AUC from the same public GEO data in one command? The 2026-10-10 sweep said no. It is 20 days old.

**Rules fixed before searching.** "Same public data" = the code itself pulls public GEO bacterial-vs-viral cohorts (overlap with our five GSE ids is noted, not required). "Both numbers" = one script/command prints a within-dataset or random-split AUC AND a train-on-others/test-on-held-out-cohort AUC. FIRED only if all three are verified by READING code/README, not an abstract. Short of that = "unclear", never rounded to "not fired". Nothing was cloned, installed or run. the author approved the read-only web search 2026-10-30 (plan approval).

**Search strings (one date, 2026-10-30; hit counts in brackets).**
- GitHub code, per id: `GSE63990`, `GSE42026`, `GSE60244`, `GSE40396`, `GSE6269` (all returned hits, mostly refinebio/accession lists; relevant: `ambaran3/LRTI_CV`, `katrinakalantar/NMLHC`). Phrase strings in code search (`leave-one-dataset-out bacterial viral`, `leave-one-cohort-out bacterial viral`, `host response cross-cohort AUC`, `bacterial viral GEO classifier held-out cohort`): [0 each].
- GitHub repos: `host-response bacterial viral` [2], `host-response cross-cohort` [3], `leave-one-dataset-out` [20, all unrelated], `leave-one-cohort-out transcriptomic` [1, long-COVID], `MetaIntegrator` [2], `bacterial-vs-viral transcriptomic classifier` [0], `bacterial viral infection gene expression classifier` [0], `infection signature GEO validation AUC` [0]. (Quoted multi-word repo searches return 0 on GitHub; the hyphenated forms work.)
- Europe PMC 2025-01 to 2026-10: the five GSE ids AND bacterial AND viral [12]; abstract `"leave-one-dataset-out"` AND bacterial AND viral [0]; host response AND bacterial AND viral AND (held-out OR external validation OR leave-one) AND (github OR code availability), 2025-09 to 2026-10 [20]; cross-cohort / leave-one-study/cohort/dataset-out AND bacterial AND viral AND transcript* [187, nearly all off-topic]; preprints only (SRC:PPR) 2026-09-01 to 2026-10-30, host response AND bacterial AND viral [1: a review of RNA host-response tests], same window with the GSE ids [0]; forward citations of Bodkin 2022 (PMID 35184750) [38, 14 from 2025-26, titles read].
- PubMed (MCP) 2025-2026: three natural-language/GSE queries [0 each; PubMed does not index GSE ids as text, so Europe PMC carried this]. The bioRxiv MCP has no keyword search (date/category only), so medRxiv/bioRxiv was covered through Europe PMC preprints instead (limit).
- Web (standard search): two broad queries (GitHub / Zenodo / Hugging Face / Kaggle / GitLab plus the concepts). They returned papers only; no repo on those hosts.
- Papers with Code: `paperswithcode.com` now redirects to `huggingface.co/papers/trending` (HTTP 302). Source gone; covered by GitHub + web. Listed as a coverage gap.

**Candidate table.**

| Candidate | URL | What it computes | Runnable | Both numbers, one command | Same public data | Verdict |
|---|---|---|---|---|---|---|
| bvnGPS (Ritchiegit), last commit 2023-02-01 | github.com/Ritchiegit/bvnGPS | 3-class (bacterial/viral/non-infected) gene-pair-score model on ~16 GEO cohorts incl. GSE40396, 42026, 60244, 63990, 6269. `main.py --test_mode cohort --test_cohort_index N` holds one cohort out; `--test_mode random` does a random split; AUC comes from a separate `predict_and_roc_0719.py` | unknown, probably not as shipped: no requirements file, no `data/` folder in the repo, hard-coded `data/coconut_...` paths, R download script with a fixed mirror | NO: two invocations plus a separate scoring script; "random" is one split, not CV | yes (GEO, same five series) | does NOT meet (fails "one command"). **Closest hit; understated on 10.10.** |
| High-AUROC-sepsis-decision-failure (Zhenghongwei11); paper PLOS ONE 2026-09-03 | github.com/Zhenghongwei11/High-AUROC-sepsis-decision-failure | sepsis vs control, GSE65682 discovery with 5-fold CV inside it, then external cohorts GSE95233 / GSE154918 / GSE28750; compares scalers incl. sample-wise rank normalisation (our within-sample control). README: "one-command reproduction" from GEO download | yes by README (not run) | yes (internal CV AUROC + external-cohort AUROC) | GEO yes, but the task is sepsis, NOT bacterial-vs-viral, and other cohorts | does NOT meet (condition 3, different task). **Same design pattern; new since 10.10.** |
| ambaran3/LRTI_CV (Falsey et al., Nat Commun 2025-11) | github.com/ambaran3/LRTI_CV | 4-gene LRTI signature; nested-CV scripts plus one validation script per cohort (GSE60244, GSE63990, ...) | partly (R scripts, data folders checked in) | no: separate scripts | partly (RNA-seq, own cohort) | does not |
| rmaghembe1/host-response-module-transportability | github.com/rmaghembe1/host-response-module-transportability | fixed locked modules projected onto external cohorts; PLOS ONE 2026-09-10; README: "not a diagnostic classifier-discovery study" | unknown | no (no trained classifier, no random-CV) | yes (GEO) | does not. Changed since 10.10.: published, v1.1.0 2026-09-12, same design |
| katrinakalantar/NMLHC | github.com/katrinakalantar/NMLHC | noisy-label logistic regression; trains GSE60244, tests GSE42026 | no (MATLAB via R, local paths); last push 2019 | no | yes | does not |
| ashm97/Developing-an-infection-state-predictive-model (Myall 2020) | github.com/ashm97/Developing-an-infection-state-predictive-model | batch-corrected merged data, 60/20/20 split | partly | no held-out cohort | yes | does not |
| MetaIntegrator | cran.r-project.org/package=MetaIntegrator | leave-one-dataset-out meta-analysis | CRAN still archived (2024-10-14); `cran/MetaIntegrator` mirror last push 2023-05-22; no khatrilab GitHub repo found by API | no one-command random-CV + held-out report | yes | does not. UNCHANGED since 10.10. |
| Bodkin 2022 (Genome Med) | doi 10.1186/s13073-022-01025-x | within-dataset lasso / LOO, no held-out cohort | availability statement reads "data available from the authors on request"; no repo | no | yes | does not. UNCHANGED |
| Chawla 2022 (Cell Syst) | not found by any code search | robustness / cross-reactivity | n/a | no | yes | does not. UNCHANGED (search gap, see below) |
| 29-mRNA classifier (Mayhew 2020) / TriVerity / MeMed BV | papers / commercial | cross-study training, prospective validation | no public code found | no | n/a | does not. UNCHANGED |
| UC two-transcript model (Front Immunol 2025); febrile-children ANN/RF (Front Pediatr 2025) | PMC12491250; PMC12364931 | single random split, no code link | no | no | partly | do not |
| Other repos (keratitis, sepsis LODO, long-COVID LOSO) | see search strings | cross-cohort designs on other diseases | partly | partly | no (other diseases) | do not (condition 3) |

**Recall check (the control).** Items from 10.10. found again by generic strings, without naming them: `rmaghembe1/...` (repo search `host-response bacterial viral`), bvnGPS (web search), Bodkin 2022 (web search), 29-mRNA / Mayhew 2020 (web search). Found only when searched by name: MetaIntegrator (the `leave-one-dataset-out` strings did not surface it). **SEARCH GAPS (not found again by any generic string):** Chawla 2022, TriVerity and MeMed BV (the last two are commercial, no code expected). Reading: generic strings recover public code reasonably (4 of 4 code-bearing items, MetaIntegrator only by name), but GitHub code-search phrase queries return 0, so repos whose README does not say `host-response` or `leave-one-...` could be missed.

**What changed since 10.10.** (1) bvnGPS was described then as "gene-pair signatures, 2,680 samples, 16 cohorts". Read properly, it has a `cohort` held-out mode and a `random` mode and uses our five GSE ids. A near-miss, not a scoop, but not what the 10.10. line said. (2) A sepsis repo with the "internal CV vs external cohort in one command, plus sample-wise rank normalisation" design exists (Sept 2026). Different disease. (3) rmaghembe1 is now published. (4) Nothing relevant in the Sept-Oct 2026 preprints.

**Verdict line.** **CORE §3 scooped: not fired** (no row meets all three conditions). Two near-misses go to the author instead of being decided here: bvnGPS fails only on "one command / runnable as shipped" (could change if someone wraps it); the sepsis repo fails only on disease. The README claim "no public runnable code computes both numbers from the same public data" still holds for bacterial-vs-viral, but a reader should be told that the same design exists for sepsis. No README change made (plan rule).

**Limits.** One date; search engines and GitHub code search are incomplete ("not found" is not "does not exist"); private or unindexed code is invisible; nothing was run, so "runnable" is judged by reading; Papers with Code is gone; bioRxiv has no keyword API, so preprints rely on Europe PMC's index; PDF-only supplements were not read.

**Claim-check.** Claim = README Prior-art line still holds on 2026-10-30. Measured = table above. Rival = the search missed it. Control = recall check: 4 of 4 code-bearing 10.10. items re-found (one only by name), 3 non-code items not found (listed as gaps). Judge = the repos' own READMEs and entry-point code read by us; no outside judge yet.

**Where.** This entry only; a worker branch (local commit, no merge, no remote, CORE §7).

## 2026-10-30 — Back-port of the public-copy corrections into the private working repo

Process only. The 2026-10-23 corrections made in this public copy (blind re-read entry, README 10-seed-means sentence, "fold seed 0" column label) were carried back into the private working repo, and a second blind re-read checked every number of that sentence against the six CSVs (0 mismatches). No numbers changed.

## 2026-11-01 — This repository switched to public

**Question.** May this work leave the lab? The author looked at this copy and decided on 01.11.2026 to make it public.

**Done.** Read-only pre-check of the whole history before the switch: one commit identity (the real-name noreply address); 0 hits for private e-mail addresses, phone numbers, local paths and secret shapes. Visibility switched; checked as a logged-out visitor (repo page, API `private=false`, raw README).

**Measured.** Nothing scientific; no code, results or config touched.

**Belief.** Unchanged. CORE §7 now allows naming the results outside the repo, always with both validation schemes side by side.

## 2026-11-05 — Re-run of the headline table: numbers identical

`python run_pipeline.py --train` on cached data, CPU only, 1 min 19 s: `results/summary_raw.csv` and all 20 fold rows byte-identical (logistic 0.9566 / 0.8807, gap 0.0759; LightGBM 0.9543 / 0.8272, gap 0.1271). Only the `date` field in `results/run_raw.json` moved. Third identical end-to-end run.

## 2026-11-05 — PLAN (not run): CORE §7 item (4), does the per-cohort drop track platform, age or class mix?

**Status of CORE first.** All four §2 conditions are met (1: 18.10./22.10.; 2: 21.10. + author checks 23.10.; 3: rewritten by the author 19.10.; 4: 20.10., restated 22.10.). The open CORE work is §7's proposal list; (1) and (3) are done, (2) is parked by the 19.10. decision, (4) per-cohort error analysis was never started. This plan is for (4). Interpretation flagged to the author, who may redirect.

**Question.** Held-out AUC differs by cohort (logistic 0.83–0.94, LightGBM 0.70–0.90, `results/folds_raw.csv`). Does the drop line up with platform family (Affymetrix: GSE63990, GSE6269; Illumina: GSE42026, GSE60244, GSE40396), age group (children: GSE42026, GSE40396, GSE6269; adults: GSE63990, GSE60244 — to be pinned to the papers in `config/cohorts.json` before the run), or class mix (24–80 % bacterial)?

**Design.** Gate 0: new runner reproduces the five on-file held-out AUCs to 0.005. Arm A: 20 ordered train-one/test-one cohort pairs, logistic L2 on raw and within-sample scale, LightGBM model seeds 0–2 secondary; pairs tagged same/cross family and same/cross age (2×2 cells: 2 / 6 / 6 / 6). Control: every training cohort cut to n = 30, 20 stratified draws, logistic. Arm B: leave-one-cohort-out with `class_weight=None` instead of `balanced`, raw, both models; plus the descriptive rank order of bacterial share vs per-cohort gap (n = 5, no verdict).

**Pre-registered rules.** R-platform: same-family minus cross-family mean AUC (logistic, raw) ≥ 0.03 AND size-matched control same sign ≥ 0.02 → "tracks platform family"; ≤ 0.01 → "does not"; between → undecided. R-age: same thresholds. R-both: both fire → report 2×2 cell means, "confounded with five cohorts", no cause named. R-scale: platform fires on raw and the difference halves on within-sample → "mostly a scale offset". R-balance: mean held-out AUC moves ≥ 0.02 when weighting is off → "class weighting matters here"; < 0.01 → "not where the drop lives". README untouched in this round whatever fires.

**Rivals + controls.** Training size → size-matched arm. Platform = age in disguise → 2×2 grid + R-both. Seed luck → logistic deterministic, LightGBM 3 seeds. Label quality per cohort → NOT controllable, stated limit. n = 5 → effect sizes only, "pattern" never "cause". One-cohort training ≠ four-cohort training → Gate 0 anchors; matrix read as "which training cohort carries over".

**Judge.** These rules, committed before any run; blind re-read of the RESULT sentence vs the CSVs by a fresh agent (same model family); tests red on broken code (test cohort leaking into training, size cap not applied, wrong tag, class-weight switch ignored).

**Runtime.** Calibrated on tonight's `--train` (1 min 19 s for 20 fits): ≈ 10–15 CPU-minutes in all. € 0, 0 GPU. Outputs: `results/transfer_matrix.csv`, `results/transfer_sizematched.csv`, `results/balance_control.csv`, `results/run_transfer.json`; code `src/transfer.py`, `run_transfer_matrix.py`, `tests/test_transfer.py`.

**Where.** This entry (its HTML page was an internal working page and is not included). Branch `night/2026-11-05` in the private repo.


## 2026-11-05 — PRE-RUN CLARIFICATIONS for the transfer-matrix plan (the transfer-matrix run; committed BEFORE any code or run)

No threshold of the 2026-11-05 PLAN entry changes. This entry only fixes how they are read, and records the Step-0 age pin.

**Step 0, age pin (`config/cohorts.json`, new keys `platform_family`, `age_group`, `age_source`).** Affymetrix: GSE63990 (GPL571), GSE6269 (GPL96). Illumina: GSE42026 (GPL6947), GSE60244, GSE40396 (GPL10558). Age groups: GSE42026 children (Herberg, PMID 23901082), GSE40396 children (Hu, PMID 23858444), GSE6269 children (Ramilo, PMID 17105821), GSE60244 adults (Suarez, PMID 25637350, median age 61). **GSE63990 is NOT pinned**: Tsalik (PMID 26791949) Table 1 gives ages 6-88 y (mean 45) for the derivation cohort, minors and adults mixed, no under-18 / 18+ split. The PLAN's guess of "adults" for it was wrong. Consequence, by the plan's own rule: **R-age = "not testable"; the 2x2 grid and R-both are skipped.** Only R-platform, R-scale and R-balance are read. Age is shown descriptively for the 4 pinned cohorts only (no verdict word).

**Readings fixed now.**
(a) "Mean pairwise AUC" = unweighted mean over the ordered pairs of a cell (same-family 8 pairs vs cross-family 12); in the size-matched arm each pair's AUC is first averaged over its 20 draws.
(b) "tracks" = main difference >= 0.03 AND size-matched difference >= 0.02, same positive sign, edges inclusive. "does not" = main difference <= 0.01 (negative counts). Anything else = "undecided".
(c) R-scale is evaluated only if R-platform says "tracks": within-sample difference <= half the raw difference (logistic, same 20 pairs).
(d) R-balance per model: delta = mean leave-one-cohort-out AUC without weighting minus with weighting over the 5 held-out cohorts. |delta| >= 0.02 for either model = "class weighting matters" (direction + per-cohort deltas reported); |delta| < 0.01 for BOTH = "not where the drop lives"; otherwise "undecided".
(e) Size-matched draw: n_bacterial = round(30 x training cohort's bacterial share), rest viral, rng seeds 0-19 fixed in code, test cohort whole. GSE40396 has exactly 30 arrays, so its 20 draws are identical (stated in the result).
(f) LightGBM and within-sample runs are secondary: ranges over 3 model seeds, no verdict word.
(g) Gate 0 covers all TEN on-file held-out AUCs of `results/folds_raw.csv` (logistic 5 + LightGBM seed 0, 5), tolerance 0.005 each. If one misses, nothing is written to `results/`.
(h) Class mix: descriptive rank order of bacterial share vs per-cohort gap plus Spearman rho on n = 5, labelled descriptive.
(i) All pairs use the one shared 12,040-gene set of the 458-array matrix, not a per-pair intersection.


## 2026-11-05 — RESULT: per-cohort drop vs platform / age / class mix (CORE §7 item 4; the transfer-matrix run; a worker branch)

**Question (recap).** Held-out AUC differs by cohort. Does the drop line up with platform family, age group or class mix? Rules were pre-registered in the 2026-11-05 PLAN entry (c7202b2) and the pre-run clarifications (fb2b99b), both committed before any run; no threshold was changed.

**Design.** Gate 0 first: the new runner reproduced all ten on-file leave-one-cohort-out AUCs (`results/folds_raw.csv`), max abs difference 6.4e-04 (tolerance 0.005). Arm A: 20 ordered train-one/test-one pairs (8 same-family, 12 cross-family), logistic L2 raw and within-sample, LightGBM raw at model seeds 0-2. Size-matched control: every training cohort cut to 30 arrays, 20 stratified draws, logistic raw. Arm B: leave-one-cohort-out with `class_weight=None` vs `balanced`, both models. **Validation schemes named:** single-cohort transfer pairs (Arm A + control) and leave-one-cohort-out (Gate 0 + Arm B); no random-CV number is produced or quoted in this round.

**RESULT.** On the 20 ordered pairs (logistic L2, raw, class-weighted) same-platform-family pairs scored LOWER than cross-family pairs: mean AUC 0.760 vs 0.821, difference -0.061 (size-matched control -0.062). **R-platform = "does not"** (rule: difference <= 0.01, negative counts) - the sign is opposite to "the drop tracks platform family". **R-age = "not testable"**: GSE63990 (Tsalik, Table 1: ages 6-88 y) cannot be pinned to children or adults, so the 2x2 grid and R-both were skipped. **R-scale** not evaluated (platform did not fire). **R-balance = "undecided"**: switching class weighting off moved the mean held-out AUC by +0.0005 (logistic) and +0.0195 (LightGBM, mean over 3 model seeds; per seed +0.031 / +0.026 / +0.002), so neither the >= 0.02 nor the < 0.01-for-both band is met. Class mix (descriptive, no verdict): Spearman rho = -0.6 between bacterial share and the per-cohort random-CV minus held-out gap, n = 5.

**What the numbers mean in plain words.** If platform were the cause, training on a cohort from the same chip family should transfer better. It transfers worse, and it does so with equal training sizes, so training-set size is not the explanation. Secondary: within-sample ranking, difference -0.034; LightGBM seeds -0.114 / -0.111 / -0.108 (same sign, larger). The class-weight switch does not move the logistic model at all; for LightGBM it moves the result by about 0.02 on average, mostly through GSE6269 (+0.080, the 80 % bacterial cohort), and the size depends on the model seed - this is why the verdict stays "undecided" and not "matters".

**Disclosures.** (1) Reading (d) "mean over the 5 cohorts" for R-balance was implemented as "per-cohort mean over the model seeds, then mean over cohorts"; this was in the code committed before the run but NOT written in the pre-run clarifications entry. Seed 0 alone (+0.031) would have crossed 0.02. (2) Three bugs were fixed between the first full run and this one, none by looking at a verdict: `matrix.transform` resolved to the pandas method (first run printed NaN, an empty arm, so its output was discarded and a NaN guard added), the same slip in the class-mix step (NaN Spearman, guard added), and a `cell_mean` indexing error. Each time results/ was deleted and the whole run repeated; the numbers above are from the last, clean run. (3) GSE40396 has exactly 30 arrays, so its 20 size-matched draws are identical and its size-matched rows equal its main rows. (4) `class_weight=None` is written as the string "None" in `balance_control.csv`; read it with `keep_default_na=False` or pandas turns it into NaN.

**Rivals.** Training size: controlled (size-matched arm). Platform = age in disguise: not testable (see R-age). Seed luck: logistic deterministic; LightGBM 3 seeds, same sign. **Not controllable (stated limit): per-cohort label quality.** n = 5 cohorts means a *pattern*, never a *cause*. The finding is a negative for platform family in this design, not a statement about why the drop exists; the drop stays unattributed.

**Guard tests.** `tests/test_transfer.py` (22 tests, suite 156 green). Mutations checked red on scratch copies: (1) leakage assert removed -> `test_leaking_test_cohort_into_training_raises` fails; (2) size cap 30 -> 31 -> `test_size_cap_applied_in_every_row_and_stratified` fails; (3) same/cross tag comparison inverted -> `test_real_config_tags_are_pinned` + `test_cell_mean_averages_draws_first_then_pairs` fail; (4) class_weight override dropped -> 3 class-weight tests fail.

**Claim-check.** Claim: the per-cohort drop does not line up with platform family in this single-cohort-transfer design. Measured: 20 pairs + size-matched control (the claim matches what was measured, but only for logistic raw; LightGBM and within-sample are secondary). Rival with control: training size. Judge: pre-registered rules in a prior commit + blind re-read below, not us. Not claimed: any cause, anything about age, anything about class mix beyond a descriptive n = 5 correlation. README untouched.

**Where.** `results/transfer_matrix.csv`, `transfer_sizematched.csv`, `balance_control.csv`, `run_transfer.json`; code `src/transfer.py`, `src/transfer_verdict.py`, `run_transfer_matrix.py`. A worker branch, local only (no remote); not merged to master. Age pin re-checked on 05.11. against Tsalik Table 1 (PMC4905578): derivation cohort mean 45 y, range 6-88 y, confirmed.

**Blind re-read:** a fresh subagent (same model family) was given only the RESULT sentences, the rules and the CSV paths and recomputed everything: all numbers MATCH (same 0.7603 / cross 0.8212 / diff -0.0609; size-matched -0.0617; balance +0.00048 logistic, +0.01949 LightGBM, per seed +0.0310 / +0.0259 / +0.0016; Spearman -0.6, p = 0.28; within-sample -0.0344; LightGBM per-seed -0.1144 / -0.1106 / -0.1077) and both verdict words follow from the rules. BLIND VERDICT: all match. (Same model family, so not an independent judge of the design, only of sentence-vs-CSV.)

**Belief.** In single-cohort transfer the per-cohort drop does not track platform family (same-family pairs transfer worse, also at equal training size); age cannot be tested with these five cohorts while GSE63990 has no per-sample ages; class weighting is at most a ~0.02 LightGBM effect that moves with the model seed, nothing for logistic. The drop stays unattributed (five cohorts = pattern, not cause). Merged to master 05.11.2026 after the author's OK. Next, as planned follow-ups: Arm B with 10 LightGBM seeds (settle R-balance, EUR 0, CPU) and a read-only check whether GEO carries per-sample ages for GSE63990 (would make R-age testable).


## 2026-11-05 — PLAN (not run): Arm B with 10 LightGBM seeds, settle R-balance (the 10-seed Arm B run; committed BEFORE any code or run)

**Question.** the transfer-matrix run left R-balance "undecided": switching class weighting off moved the mean leave-one-cohort-out (LOCO) AUC by +0.0005 (logistic) and +0.0195 (LightGBM, mean over model seeds 0-2; per seed +0.031 / +0.026 / +0.002). The LightGBM mean sits 0.0005 under the 0.02 line while the per-seed spread is 0.03, so the 3-seed verdict is seed luck either way. Do ten seeds settle it, or show it cannot be settled at this n? The 3-seed verdict in the 2026-11-05 RESULT entry stays on record; this is a NEW entry, not an edit.

**Design (fixed now).** LightGBM model seeds 0-9 (the original 0-2 included). Delta per model = mean over the 5 held-out cohorts of [per-cohort mean over the seeds of AUC(class_weight=None) minus AUC(balanced)], exactly reading (d) of the pre-run clarifications. Logistic is deterministic: re-run once, as the Gate-0 anchor and for the logistic delta. Gate 0 (ten on-file LOCO AUCs, tolerance 0.005) runs first and refuses to write on a miss. Arm A is NOT re-run (R-platform did not depend on seeds, all three about -0.11, same sign).

**Thresholds, unchanged.** |delta| >= 0.02 for either model = "class weighting matters here"; |delta| < 0.01 for BOTH models = "not where the drop lives"; anything else = "undecided". Edges inclusive for the 0.02 line. No threshold change after seeing data.

**Also reported next to the verdict (descriptive, no verdict word):** the 10 per-seed deltas (LightGBM, mean over cohorts per seed), their min, max and sd (ddof=1), and the per-cohort deltas. If the 10-seed mean lands in the undecided band again, the answer is "cannot be settled at this n", stated as such.

**Rivals + controls.** Seed luck: ten seeds instead of three, spread reported. Cohort luck (n = 5): not controllable, stated limit; effect sizes only, "pattern" never "cause". Judge: these rules committed before any code, a regression anchor (seeds 0-2 of the new run must reproduce the committed `results/balance_control.csv` to 1e-6), and a blind re-read by a fresh subagent.

**Runtime / cost.** CPU only, about 5-10 minutes, EUR 0. New outputs only: `results/balance_control_10seeds.csv`, `results/run_balance_10seeds.json`; the transfer-matrix run result files stay byte-identical. Code: `--balance-seeds N` option in `run_transfer_matrix.py` (default 3 keeps today's behaviour), seed list in `src/transfer.py`. A worker branch, local only (no remote). No HTML report.


## 2026-11-05 — RESULT: Arm B with 10 LightGBM seeds (the 10-seed Arm B run; a worker branch; PLAN c146725)

**Question (recap).** the transfer-matrix run left R-balance "undecided": dropping class weighting moved the mean leave-one-cohort-out (LOCO) AUC by +0.0005 (logistic) and +0.0195 (LightGBM, 3 seeds), 0.0005 under the 0.02 line, with a per-seed spread of 0.03. Ten seeds were meant to settle it or show it cannot be settled.

**Design.** Exactly as pre-registered in the PLAN entry (committed before any code): Gate 0 first (PASS, max abs diff 6.4e-04, tolerance 0.005); LightGBM model seeds 0-9; per cohort the mean over seeds, then the mean over the 5 held-out cohorts (reading d); thresholds unchanged (>= 0.02 either model = "matters"; < 0.01 both = "not where the drop lives"; else "undecided"). Logistic re-run once (deterministic). Arm A not re-run. Files: `results/balance_control_10seeds.csv`, `results/run_balance_10seeds.json`; the transfer-matrix run result files are byte-identical to master.

**RESULT.** **R-balance = "not where the drop lives".** Logistic delta (no weighting minus weighting) +0.0005; LightGBM 10-seed delta +0.0033; both under 0.01. Per-seed LightGBM deltas, seeds 0-9: +0.0310, +0.0259, +0.0016, -0.0179, -0.0286, +0.0217, -0.0261, +0.0206, -0.0265, +0.0316. Min -0.0286, max +0.0316, sd 0.0257 (ddof 1). Seeds 0-2 alone reproduce the old +0.0195 (regression anchor, 1e-6).

**What the numbers mean in plain words.** Whether class weighting is on or off changes the LightGBM score by about a third of a percentage point on average, and the sign flips from seed to seed (5 seeds positive, 5 negative-or-near-zero). The 3-seed "+0.0195" was a lucky draw of seeds 0-2, which happen to be mostly positive. Weighting does not systematically move the held-out AUC. Caveat: a single seed can still move by +-0.03, so any single LightGBM run is itself noisy at that size; that noise is the model's seed variance, not an effect of weighting.

**Claim-check.** Claim: class weighting is not where the per-cohort drop lives. Measured: LOCO AUC with and without weighting, 5 cohorts, LightGBM 10 seeds + logistic. Rival + control: seed luck, controlled by 10 seeds and the spread reported; cohort luck (n = 5) not controllable, stated limit, a pattern not a cause. Judge: the pre-registered rule plus a blind re-read by a fresh subagent given only the sentences, the rule and the CSV path: all 5 sentences MATCH (it noted that the verdict depends on using all 10 seeds: the 3-seed subset alone is "undecided").

**Belief.** R-balance is closed at "not where the drop lives"; the drop stays unattributed. No HTML report (skipped, per plan). README and CORE.md untouched. Local only (no remote); not merged.


## 2026-11-06 — RESULT: does GEO carry per-sample ages for GSE63990? (the GEO age check; a worker branch)

**Question.** the transfer-matrix run could not test R-age ("does the per-cohort drop track age group?") because GSE63990 (Tsalik 2016) is a mixed cohort with no children/adult split in the paper's Table 1. Does GEO's series matrix carry a per-sample age that would let us split it? This is a count of metadata fields, not a test, so no pre-registration; nothing from the pipeline was run and no result changed.

**Measured** (`tools/geo_age_check.py` -> `results/geo_age_check.csv`, cached matrices only, nothing downloaded). **GSE63990: 280 samples in the matrix, 0 carry an age.** Its only characteristics are `batch`, `infection_status`, and tissue/sample-id text; no key contains "age". Under-18 count: not available. Cross-check cohorts: GSE42026 (92 samples) also has no age field. GSE40396 (65): age in months (key `age (month)`), 0.2-2.75 y, all under 18. GSE6269 (97): key `Age`, mixed units ("10 y", "4m"; converted), 0.04-36 y, median 0.9 y, 94 under 18 and 3 aged 18 or over. GSE60244 (158): key `age`, 23-92 y, all adults. So the `age_group` pins in `config/cohorts.json` hold for three cohorts; GSE6269 "children" contains 3 samples of 18+ (a small impurity, not changed here; the config is untouched). Per the plan, no GSM pages or SRA were scraped.

**Belief.** **R-age testable by splitting GSE63990: no, GEO carries no per-sample age for GSE63990.** Next option (not filed): fetch the per-patient table in the Tsalik 2016 paper's supplement (or the authors' data), if one exists; that would be a different source than GEO's matrix.

**Where.** `tools/geo_age_check.py`, `results/geo_age_check.csv`, `tests/test_geo_age_check.py` (4 tests; suite 168 passed). `config/cohorts.json`, README, CORE.md, `src/` and all pipeline results untouched. Local only (no remote).

## 2026-11-07 — Re-run of the headline table: numbers identical again

Same command as on 05.11., 1 min 33 s, md5-checked before and after: `results/summary_raw.csv` and `results/folds_raw.csv` byte-identical. Only the `date` field in `results/run_raw.json` moved (to 2026-11-07). Fourth identical end-to-end run.

## 2026-11-07 — R-age closed from public sources (adds the paper check to the GEO age check; nothing run)

**Question.** the GEO age check found no per-sample age in GEO for GSE63990 and named the paper's supplement as the next place to look. Does Tsalik 2016 (PMID 26791949, PMC4905578, doi 10.1126/scitranslmed.aad6873) give per-subject ages that could be mapped to the 280 arrays?

**Checked (read-only, PubMed Central full text + the PMC page's table captions).** Table 1: cohort-level age only ("45 (6-88)" for the derivation cohort). Table S2 is "Summary of clinical features for the derivation cohort" (aggregate). Table S4, "Subjects with discordant predictions compared to clinical assignments", is the only per-subject table and covers the 35 discordant subjects, not 273. No GSM ids appear anywhere in the text. The supplement PDF (`NIHMS791803-supplement-01.pdf`, 313 KB) could not be fetched by script tonight (PMC and Europe PMC both answered non-browser requests with HTML/XML, not the file); its table captions were read from the article page. Related fact from the Methods, relevant to the plan below: the 317 arrays "were measured in two microarray batches with seven overlapping samples", with "substantial batch differences" corrected in the authors' analysis.

**Belief.** **R-age is not testable from public sources** for GSE63990: GEO has no age field, the paper has no per-subject age table for the full cohort. Closed unless the author opens the supplement PDF in a browser and finds otherwise (an impossibility for 280 arrays even then, since S4 holds 35 subjects). `config/cohorts.json` stays "mixed / NOT PINNED".

**Where.** This entry only. The GEO age check branch (a worker branch, commit 674d018, its own LABBOOK entry 06.11.) is still unmerged into master; the author's merge.


## 2026-11-07 — PLAN (not run): does GSE63990's own array batch already cost AUC? (CORE §7 item (4), continued; committed BEFORE any code or run)

**Status of CORE first.** All four §2 conditions met (05.11. check, unchanged). §7 list: (1) done, (2) parked 19.10., (3) done, (4) transfer matrix done 05.11. (R-platform "does not track", R-balance "not where the drop lives", R-age closed above). The drop stays unattributed. This plan continues (4) on the one cohort with a documented internal structure that is not patient biology: GSE63990's two microarray batches.

**Question.** Trained on batch 0, tested on batch 1 of the same hospital cohort (and the reverse): how much AUC is lost compared with a random split of the same shape, and does the loss shrink under within-sample normalisation?

**Data.** The 185 labelled GSE63990 arrays already in the pipeline; `batch` from the cached series matrix. Counted tonight with the pipeline's own rules and exclusion list (counts only, no AUC looked at): batch 0 = 58 bacterial / 88 viral (146), batch 1 = 12 / 27 (39); total 70 / 115 = Table 1.

**Design.** Gate 0: the runner rebuilds exactly these 185 arrays and reproduces the on-file held-out-GSE63990 AUCs (logistic 0.9257, LightGBM seed 0 0.8927) to 0.005; on a miss nothing is written. Gate B: logistic 5-fold CV predicting *batch* from expression (labels unused); batch-AUC ≥ 0.80 → proceed, else stop with verdict "batch not visible in GEO's values; test uninformative". Arm A: direction 1 train batch 0 / test batch 1; direction 2 the reverse; logistic L2 raw + within-sample primary, LightGBM seeds 0-2 raw secondary (ranges, no verdict). Control: per direction 20 random draws (rng seeds 0-19) with the same training n and class counts as the training batch and the same test n and class counts as the held-out batch, both sides mixing batches; Δ = control mean − batch-out AUC.

**Pre-registered rules.** R-batch (logistic raw): "a within-hospital batch costs AUC" = Δ ≥ 0.03 in BOTH directions AND batch-out AUC below all 20 control draws in both; "no batch cost" = Δ ≤ 0.01 in both; else "undecided"; edges inclusive. R-scale, only if R-batch fires: within-sample Δ ≤ half the raw Δ in both directions → "mostly a scale offset". Descriptive, no verdict word: Δ next to the overall logistic gap 0.076 and GSE63990's own held-out loss (~0.03); 1000-resample bootstrap interval per batch-out AUC; the Gate-B batch-AUC; LightGBM ranges. Power limit stated now: batch 1 has 12 bacterial arrays, one AUC on 39 arrays has SE ≈ 0.07, so direction 1 sees only a large loss and "no batch cost" is weak evidence there. No threshold change after data.

**Rivals + controls.** Training size → control n. Class mix (40 % vs 31 %) → control class counts. GEO already batch-corrected → Gate B. Batch = sub-population (four sites, 48 Duke students in the viral group) → NOT controllable, stated: "batch or sub-population", pattern not cause. Seed luck → logistic deterministic, LightGBM secondary. Test noise → control spread + bootstrap + "below all 20 draws". n = 1 cohort → says nothing about the other four.

**Judge.** These rules before code; tests red on broken code (batch leaks into training, control not class-matched, Gate B skipped, rule edges inverted); blind re-read of the RESULT sentences vs the CSVs by a fresh subagent.

**Runtime / cost.** ≈ 120 fits on ≤ 185 arrays (tonight's calibration: 20 fits on 458 arrays = 1 min 33 s) → ≈ 5-10 CPU-min, € 0, 0 GPU. New files only: `results/batch_split.csv`, `results/batch_control.csv`, `results/run_batch_split.json`; code `src/batch_split.py`, `src/batch_split_verdict.py`, `run_batch_split.py`, `tests/test_batch_split.py`. Existing results byte-identical.

**Morning options (in the report).** A GO as written (recommended) · B logistic only · C Bodkin Table 2 pull instead (last listed §7 item, literature only) · D stop the lab at the current belief. Plus metaphor A/B. Also: merge the GEO age check.

**Where.** This entry (its HTML page was an internal working page and is not included). Branch `night/2026-11-07` in the private repo.

## 2026-11-08 — PRE-RUN CLARIFICATIONS for the batch-split run (the batch-split run; committed BEFORE any code or run)

Reads of the 07.11. PLAN above; **no threshold, rule or reading from it changes.** This entry only fixes what the PLAN leaves open. A worker branch, cut from `night/2026-11-07` (so the PLAN is an ancestor of the code). Approval = the author's option A (logistic + LightGBM secondary).

(a) **Gate 0** = two held-out-GSE63990 fits (train on the other four cohorts, test on GSE63990; logistic and LightGBM seed 0, raw) must match `results/folds_raw.csv` (0.9257 / 0.8927) within 0.005, AND the 185 arrays must count exactly 58 bacterial / 88 viral (batch 0) and 12 / 27 (batch 1) with exactly one `batch:` string per array. Any miss = nothing is written to `results/`.
(b) **Gate B** number = mean of the 5 per-fold AUCs (the repo convention in `train.summarize`), `StratifiedKFold(5, shuffle=True, random_state=0)` stratified on batch, same logistic settings as `src/train.py` (`transfer.make_model`); the pooled out-of-fold AUC is shown beside it (descriptive). ≥ 0.80 inclusive proceeds. Below it only `run_batch_split.json` and the Gate-B rows of `batch_split.csv` are written (no Arm A, no control file), verdict text "batch not visible in GEO's values; test uninformative" plus the number.
(c) **Control draw d (0-19)** = one random partition of the 185 arrays: the 39-array side is drawn with `default_rng(d)` at 12 bacterial / 27 viral, the other side is the remaining 146 (58 / 88). Direction 1 trains on the 146 side, direction 2 on the 39 side, so train/test n and class counts equal the batch split's by construction, the sides are disjoint and cover all 185. Each draw is asserted to contain both batches on both sides. Both directions and both transforms use the same 20 partitions. Within-sample Δ is defined like raw Δ (control mean − batch-out AUC); always computed, only READ if R-batch fires.
(d) **Edges.** "Below all 20 draws" = strictly below the minimum of the 20 control AUCs (a tie is not below). Δ edges inclusive with EPS 1e-12 as in `transfer_verdict.py`. "No batch cost" = Δ ≤ 0.01 in both directions (negative Δ counts). Every other combination = "undecided".
(e) **Bootstrap.** 1000 resamples of the test arrays (rng seed 0), model NOT refitted, 95 % percentile interval, resamples with one class skipped and counted. Descriptive only.
(f) **LightGBM** = Arm A only, model seeds 0-2, raw; no LightGBM control draws (the PLAN has none); min/max ranges only, no verdict word.
(g) **Comparators** are read from files, not typed in: overall logistic gap from `results/summary_raw.csv`, GSE63990's own held-out loss from `results/gapdiag_cohorts.csv` (`rcv_auc_in_cohort − loco_auc`, raw, logistic).

**Disclosures fixed now (for the RESULT entry).** (1) The pipeline's highest-mean probe choice and per-cohort NaN imputation read BOTH batches of GSE63990 (label-free; identical to the on-file pipeline Gate 0 reproduces; affects batch split and control equally). (2) The 7 re-run arrays were excluded by an INFERRED identity (2026-10-21); a subject still in both batches would inflate cross-batch AUC, biasing toward "no batch cost". (3) Batch can be a sub-population (four sites, 48 students): wording stays "batch or sub-population", pattern not cause. (4) Power: 12 bacterial arrays in batch 1, one AUC on 39 arrays has SE ≈ 0.07, so "no batch cost" in direction 1 is weak evidence. (5) One cohort says nothing about the other four. README.md and CORE.md stay untouched.

## 2026-11-08 — RESULT: GSE63990's own array batch, split both ways (the batch-split run; PLAN 2026-11-07 + clarifications 2026-11-08)

**Question, and why.** The model loses AUC when the test cohort is a whole different hospital (overall logistic gap 0.0759). Four explanations were already tried and none attributed the drop (platform family, class balance, age; see the 05.11. and 07.11. entries). GSE63990 is the one cohort with a documented internal split that is *not* patient biology: two microarray batches (Tsalik 2016, Methods). If a model taught on batch 0 already does worse on batch 1 of the same hospital (and the reverse), then part of the drop can live in how samples were handled, not in the hospital.

**Design, before the numbers.** The 185 labelled arrays of GSE63990: batch 0 = 58 bacterial / 88 viral, batch 1 = 12 / 27. Gate 0 rebuilt them and had to reproduce the on-file held-out AUCs. Gate B asked: can a model tell the batches apart from expression alone? Arm A trained on one batch and tested on the other, both directions (logistic L2, raw and within-sample; LightGBM seeds 0-2 raw). The control repeated exactly the same train/test sizes and class counts, but with 20 random splits that mix both batches; Δ = control mean − batch-out AUC. Rules were fixed before any code (commits 4531829, 341fe29); nothing was changed after seeing data.

**Numbers in plain words.**
- Gate 0 passed: the refit differs from the file by 0.0005 (logistic) and 0.0000 (LightGBM), limit 0.005. Counts 58/88/12/27 exact.
- Gate B: batch is perfectly recognisable (batch-AUC 1.000 in all 5 folds, bar 0.80). So the two batches differ systematically in the GEO values. That alone does not say the difference hurts infection prediction.
- Direction 1 (train batch 0, test batch 1), logistic raw: AUC 0.849 (bootstrap 95 % interval 0.70-0.96). The 20 mixed control splits average 0.957 (min 0.892, max 1.000, sd 0.037). Δ = 0.108; the batch-out AUC is below all 20 draws.
- Direction 2 (train batch 1, test batch 0), logistic raw: AUC 0.885 (0.81-0.94). Control mean 0.924 (min 0.858, max 0.948, sd 0.021). Δ = 0.039, but the batch-out AUC is NOT below all 20 draws (it sits inside the control spread).
- **R-batch verdict: "undecided"** (direction 1 clears both parts of the rule, direction 2 clears Δ ≥ 0.03 but fails "below all 20 draws"). R-scale: not evaluated, because R-batch did not fire.
- Descriptive, no verdict word: within-sample normalisation gives direction 1 0.858 vs control 0.958 (Δ 0.100) and direction 2 0.770 vs control 0.930 (Δ 0.160), so it did not close the gap. LightGBM raw (seeds 0-2): direction 1 0.719-0.725, direction 2 0.899-0.904. Next to these: overall logistic gap 0.0759, GSE63990's own held-out loss 0.0329 (both read from files).

**Claim-check.** Claim tested: a within-hospital array batch already costs held-out AUC. Measured: logistic AUC across the two documented batches, both directions, against 20 size- and class-matched mixed splits of the same arrays. Rival explanations and their controls: training size and class mix (control matches n and class counts by construction); GEO already batch-corrected (Gate B says no, batch is fully visible); test noise on 39 arrays (control spread, bootstrap, "below all 20" rule); LightGBM seed luck (three seeds, ranges only). Not controllable: batch vs sub-population (four sites, 48 students in the viral group). Judge who is not us: the blind re-read below, same model family, so it checks sentence-versus-CSV, not the design. Kill outcomes: neither fired (Gate B passed; not "no batch cost").

**Disclosures.** (1) The probe choice (highest mean) and per-cohort NaN imputation read both batches (label-free, identical to the on-file pipeline Gate 0 reproduces; affects batch split and control equally). (2) The 7 re-run arrays were excluded by an inferred identity; a subject left in both batches would inflate cross-batch AUC, biasing toward "no batch cost". (3) Batch can be a sub-population: the wording is "batch or sub-population", a pattern, not a cause. (4) Power: one AUC on 39 arrays has SE ≈ 0.07 and batch 1 has only 12 bacterial arrays; direction 1 rests on that small side. (5) One cohort says nothing about the other four. (6) Runner bug, fixed and the WHOLE run repeated: a `c.transform` pandas-method clash in the comparator read crashed the first run before anything was written; the fix is a column lookup by name and cannot touch a number or a verdict. (7) `tests/` hold red-then-green proofs for the four named guards (batch leak, control not class-matched, Gate B bypass, rule edges); the control's small-side class check is redundant with its big-side check (the sides are complements), shown by removing both.

**What we now believe.** Batch is visible in the values and, in one direction, the model taught on one batch scores clearly worse on the other than a mixed split would (Δ 0.108, larger than the whole-cohort loss of 0.033 and than the overall gap 0.076); in the other direction the loss (Δ 0.039) is within what random splits can produce. The pre-registered word is "undecided", and that is what we claim: the batch-or-sub-population candidate is not excluded and not established; the cross-hospital drop stays unattributed. Within-sample normalisation did not help here. This does not close CORE §7 item (4) either way; a larger-test-side design (not available in GSE63990) would be needed to decide.

**Where.** `results/batch_split.csv`, `results/batch_control.csv`, `results/run_batch_split.json`; code `src/batch_split.py`, `src/batch_split_verdict.py`, `run_batch_split.py`, `tests/test_batch_split.py`. A worker branch (local only, no remote). Existing `results/` files byte-identical (md5 before/after, 49 files). README.md and CORE.md untouched.

**Blind re-read (fresh subagent, read-only, given only the RESULT sentences, the rules and the three file paths).** All 9 sentences MATCH the CSV/JSON; it recomputed the verdict independently and also got "undecided" (direction 2 fails "below all 20 draws"; "no batch cost" fails because Δ > 0.01). No anomalies except that `n_train`/`n_test` are stored as floats ("146.0"). Same model family: it checks sentence-versus-file, not the design.


## 2026-11-08 — WRITE-UP for a general reader (the author's GO, 20:05; docs only)

**Question.** Can the lab's state as of 07.11. (all four CORE §2 conditions met; R-platform, R-balance and R-age closed; drop unattributed) be told as one reading story for a non-specialist reader without claiming more than CLAUDE.md allows?

**Done (docs only, no experiment, no code change, `results/` untouched).**
(1) `README.md` rewritten as the reading story: question → why it matters → design in plain words → four results with their key numbers explained → what the drop is NOT (closed candidates, one table) → honest limits → who did what (author vs AI) → reproduce → prior art. The old README's Labels table, prior-art table and author-vs-AI split are kept in condensed form; every number keeps its `results/` source line.
(2) `reports/2026-11-08-writeup.html`, the same story as an HTML page with one dumbbell chart (random-CV vs leave-one-cohort-out, 4 rows).
(3) A privacy check of the tracked files (no private e-mail address, no credential-like file).
(4) This entry.

**The four results as told (numbers re-read from the CSVs tonight, all MATCH the files):** R1 raw gap logistic 0.957→0.881 (0.076; 10 seeds 0.073 ± 0.004), LightGBM 0.954→0.827 (0.127; 0.123 ± 0.006), LightGBM model-seed LOCO 0.780–0.861. R2 per-patient ranks logistic 0.963→0.895 (0.065 ± 0.003), LightGBM 0.958→0.888 (0.067 ± 0.004), paired model-seed gain +0.064 in 10/10 seeds, LOCO 0.891 ± 0.006. R3 frozen Herberg 0.876/0.882 (−0.004 ± 0.005), Sweeney 0.837/0.912 (−0.075 ± 0.002), independent-cohort LOCO 0.901. R4 same- vs cross-family 0.760 vs 0.821 (−0.061; size-matched −0.062), balance delta +0.0005 / +0.0033 (10 seeds, sd 0.026), age not testable (0/280 GEO ages).

## Claim-check (2026-11-08, write-up)
Claim:      the write-up states only what CORE/CLAUDE.md allow — the gap between two validation schemes on public data, four pre-registered results, closed candidates, limits — and every number is traceable to a committed `results/` file.
Measured:   every quoted number re-read from `summary_*.csv`, `seedsweep_*_summary.csv`, `seedsweep_lgbm_modelseed*.csv`, `frozen_signatures_summary.csv`, `folds_raw.csv`, `run_transfer.json`, `run_balance_10seeds.json` tonight; `pytest --co` = 168. No new measurement; the write-up adds no claim the LABBOOK does not already carry.
Rival:      "story inflation" — the story reads as a better classifier or a clinical finding. Control: section 1 and 8 of both documents say in so many words what is NOT claimed; the blind cold-read below checks each sentence against CLAUDE.md's allowed-claim rule.
Judge:      a blind fresh subagent (same model family) given only README.md, CLAUDE.md lines 1-30 and the `results/` files, told to recompute every number and flag any overstep. Fails if any number mismatches or any sentence oversteps. Result (2026-11-08, same night): every printed number with a data file reproduces to the printed precision (about 60 numbers checked, 0 mismatches), 168 tests collected, no sentence claims a signature/classifier/clinical use; 4 low-severity wording notes (Sweeney 2-cohort row has no random-CV cell → footnoted; 'about as well' → 'as well as or slightly better'; two source files added to the Sources lines) all applied before the commit. Same model family, so a sentence-vs-CSV judge, not a design judge. No human reader yet.

**Belief.** Unchanged (docs only).


## 2026-11-08 — PLAN (not run): whole-blood check — does the one PBMC cohort change the held-out number? (the author; committed BEFORE any run)

**Why.** A read-only pass over the cached `data/series/*.gz` metadata (this session, re-verified field by field) shows GSE6269 is the only PBMC cohort (`Sample_source_name` = `PBMC_…`, overall design "PBMC"); GSE40396, GSE42026, GSE60244 and GSE63990 say whole blood (series text; PAXgene tubes for GSE42026 + GSE63990, Tempus for GSE40396 + GSE60244). Sample type is therefore a candidate for part of the drop that was never tested. GSE6269 is also the 80 %-bacterial, pediatric, GPL96 (Affymetrix) cohort and LightGBM's worst held-out (0.703), so whatever moves here is **confounded with "that one cohort"** — stated now, not after.

**Question.** Trained and tested on the 4 whole-blood cohorts only, how do the leave-one-cohort-out (LOCO) AUCs compare with the 5-cohort run, and does the random-CV minus LOCO gap change?

**Design.** `run_whole_blood_check.py`. Gate 0: re-fit the 5-cohort LOCO rows with the same code path and reproduce all 10 on-file AUCs (`results/folds_raw.csv`) to 0.005, else nothing is written. Primary: the SAME 12,040-gene matrix with GSE6269 removed from train and test (so only the training data changes, not the gene set); logistic L2 raw, 5-fold random-CV (fold seed 0) + 4-way LOCO; LightGBM model seeds 0–2 secondary. Secondary, descriptive only: the natural 4-cohort gene intersection (will be larger than 12,040 once GPL96 is gone), logistic LOCO only, to show the gene-set effect separately.

**Pre-registered readings (logistic, raw, fixed genes).** Δ_train = mean LOCO over the 4 whole-blood cohorts when trained on 3 whole-blood cohorts, minus the mean of the SAME 4 held-out rows in the 5-cohort run (where GSE6269 was in the training set). Δ_train ≥ +0.02 → "the PBMC cohort hurts whole-blood transfer"; ≤ −0.02 → "the PBMC cohort helps"; between → "no detectable effect at n = 4". Gap: 4-cohort random-CV minus 4-cohort LOCO next to the 5-cohort 0.076; a shrink ≥ 0.02 is reported as "part of the 5-cohort gap sat with GSE6269" — descriptive, no causal word, because the 5-cohort LOCO mean also contained GSE6269's own held-out row (0.938 logistic, 0.703 LightGBM), which is dropped here by construction. LightGBM: ranges over 3 seeds, no verdict word. No threshold change after data.

**Rivals + controls.** PBMC vs. "the smallest-chip / most-bacterial / pediatric / Affymetrix cohort": NOT separable with one PBMC cohort — the result is about *removing GSE6269*, PBMC is one of its properties. Gene-set change: controlled by fixing the 12,040 genes (primary) and showing the natural intersection separately. Training-set size (458 → 367 arrays): not controlled; stated. Seed luck: logistic deterministic; LightGBM 3 seeds. n = 4 held-out cohorts: effect sizes only.

**Judge.** These readings, committed before the run; Gate 0; tests red on broken code (`tests/test_whole_blood_check.py`: PBMC cohort leaking into the whole-blood matrix, gene set not fixed, Gate 0 tolerance); blind re-read of the RESULT sentences vs the CSV by a fresh subagent.

**Runtime / cost.** ≈ 30 fits on ≤ 458 arrays, ≈ 3–4 CPU-min, € 0, 0 GPU. New files only: `results/whole_blood_check.csv`, `results/run_whole_blood_check.json`. Existing results byte-identical. Branch `docs/2026-11-08-writeup` (the author's instruction: same branch).

_**Filing disclosure (08.11., 21:50).** This PLAN text was written before the run, but the append landed in the wrong folder (a private notes file, committed there at 21:41, i.e. AFTER the run at ~21:20); it was moved here afterwards. What WAS committed before the run is the runner `run_whole_blood_check.py` (431c30b, 21:17), which fixes the design (Gate 0, fixed gene set, GSE6269 removed from train and test, logistic primary, LightGBM seeds 0–2) but not the ±0.02 reading. Treat the Δ_train reading as pre-written, not pre-committed._

## 2026-11-08 — RESULT: whole-blood check (PLAN above; runner 431c30b, see filing disclosure; run on branch `docs/2026-11-08-writeup`)

**Question (recap).** GSE6269 is the only PBMC cohort. With it removed from training AND testing, how do the leave-one-cohort-out (LOCO) AUCs of the 4 whole-blood cohorts compare with the same 4 held-out rows of the 5-cohort run, and does the gap change?

**Run.** `python run_whole_blood_check.py`, CPU, 2 min 18 s. Gate 0 PASS: all 10 on-file LOCO AUCs reproduced, max |diff| 1.1e-16. Existing `results/folds_raw.csv` and `summary_raw.csv` md5-identical before and after. New files only: `results/whole_blood_check.csv`, `results/run_whole_blood_check.json`.

**Measured (logistic L2, raw, the fixed 12,040 genes, 367 arrays).** LOCO per whole-blood cohort, 4-cohort run vs 5-cohort run: GSE40396 0.858 vs 0.864; GSE42026 0.835 vs 0.832; GSE60244 0.819 vs 0.844; GSE63990 0.914 vs 0.926. Mean 0.857 vs 0.866, **Δ_train = −0.010**. 4-cohort random-CV 0.938, LOCO 0.857, **gap 0.081** (5-cohort gap 0.076). Side fact: the natural 4-cohort gene intersection is ALSO 12,040 genes — GPL96 (GSE6269) was never the chip that limited the shared set — so the secondary natural-intersection rows equal the primary rows exactly.

**Measured (LightGBM, raw, secondary, seeds 0–2).** Mean LOCO 0.711 / 0.765 / 0.740 vs 0.858 in the 5-cohort run (Δ ≈ −0.15 … −0.09). Almost all of that is ONE cohort: held-out GSE40396 (30 arrays, 8 bacterial) falls from 0.835 to 0.216 / 0.477 / 0.324 once GSE6269 leaves the training set; the other three move little (GSE42026 0.896 → 0.908 / 0.904 / 0.916; GSE60244 0.809 → 0.853 / 0.809 / 0.868; GSE63990 0.893 → 0.866 / 0.868 / 0.850). 4-cohort random-CV 0.922, LOCO 0.711 (seed 0), gap 0.212.

**Pre-registered reading, applied as written.** Logistic Δ_train = −0.010 lies inside (−0.02, +0.02) → **"no detectable effect at n = 4"**. The gap did not shrink (0.076 → 0.081), so the descriptive line "part of the 5-cohort gap sat with GSE6269" is NOT written. LightGBM carries no verdict word; its numbers are reported as ranges.

**What the numbers mean in plain words.** For the logistic model it makes almost no difference whether the PBMC cohort is in the training set: the four whole-blood hospitals are predicted about as well either way, and the gap between mixed and held-out stays about 8 points. So "the PBMC cohort is where the drop lives" is not supported for the model the headline rests on. For LightGBM the picture is different and stranger: without the pediatric, bacterial-heavy GSE6269 in training, the tree model becomes anti-predictive on the smallest cohort (Hu 2013, 30 children, 8 bacterial; AUC 0.22–0.48 means it ranks bacterial cases as *less* bacterial). A plausible reading (inferred, not tested): GSE6269 supplies most of the training signal for young children with bacterial infections, and the tree model leans on it for GSE40396, the other young-children cohort. That is a statement about training-set composition, not about PBMC.

**Claim-check.** Claim: removing the PBMC cohort does not detectably change the logistic held-out AUC on the whole-blood cohorts (Δ −0.010, n = 4), and the random-CV minus LOCO gap does not shrink. Measured: exactly that (LOCO on 4 cohorts, fixed genes, Gate 0). Rival: PBMC is confounded with "that one cohort" (pediatric, 80 % bacterial, Affymetrix GPL96, smallest held-out LightGBM AUC) — NOT separable with one PBMC cohort, so this is a result about *removing GSE6269*; training-set size (458 → 367) not controlled, stated. Judge: readings committed before the run (431c30b); 4 new guard tests (PBMC leaking into the whole-blood matrix, gene set not fixed, Gate 0 misses a moved number / missing rows, committed CSV ≠ JSON); blind re-read by a fresh subagent given only the two new result files + `folds_raw.csv`: all 5 numbered sentences MATCH, no GSE6269 row in the CSV, verdict word follows from the rule — BLIND VERDICT: all match (same model family; sentence-vs-CSV judge only). It added one reading note, adopted here: the "5-cohort gap 0.076" is the all-five-cohort gap; computed on the same 4 whole-blood held-out rows only it would be 0.090 (random-CV 0.957 − 0.866), so the gap did not shrink under either reading. Suite: 172 tests.

**Belief.** Sample type (PBMC vs whole blood) is not where the logistic drop lives, within what one PBMC cohort can show. New open pattern, LightGBM only: GSE40396 depends on GSE6269 being in the training set (anti-predictive without it). The drop stays unattributed. **Correction to the record:** the in-flight board block of 06.11. (removed 08.11.) said "GEO carries no per-sample age for the 5 cohorts"; that was wrong — per the GEO age check entry (06.11.) GEO carries per-sample ages for GSE40396 (months), GSE60244 (years) and GSE6269 (`Age: 10 y, when sample taken`), and NOT for GSE63990 and GSE42026. README, CORE and the 07.11. entry only ever said "no age for GSE63990", which stands.

**Where.** `results/whole_blood_check.csv`, `results/run_whole_blood_check.json`, `run_whole_blood_check.py`, `tests/test_whole_blood_check.py`; README section 4 (new "where the cause could still hide" table + the sample-type row); `reports/2026-11-08-writeup.html` (same artifact link, updated). Branch `docs/2026-11-08-writeup`, local only.


## 2026-11-08 — causes table filled from the batch-split run (docs only; the author)

**Done.** README §4 "Where the cause could still hide" + `reports/2026-11-08-writeup.html`: lab-batch row now carries the batch-split run result with the verdict word exactly "undecided" (Gate B 1.000; dir 1 Δ 0.108 below all 20 controls; dir 2 Δ 0.039 inside the control spread; power limit 39 arrays / 12 bacterial), numbers re-read from `results/run_batch_split.json` on a worker branch. Tube row marked "not testable on its own (tube confounded with cohort)". Follow-up planned: **the deciding batch test** (200-draw null both directions in GSE63990 + GSE40396 hybridization-date batches behind a power gate; plan first, nothing runs before the author's approval).


## 2026-11-08 — PLAN (not run): deciding test of the batch cause, 1000-draw null + GSE40396 date batches behind a power gate (the deciding batch test, approved by the author 08.11. 22:54; committed BEFORE any code or run)

**Question.** the batch-split run left "does GSE63990's own array batch cost AUC?" at "undecided" (direction 1 below all 20 control draws, direction 2 not). With a much larger null, does the batch-out AUC sit in the tail of the size- and class-matched mixed-split null in BOTH directions?

**Findings from planning (counts only; no model fitted, no new AUC looked at).**
(1) Labelled GSE40396 arrays per hybridization date (the pipeline's own label rules): 2010-08-03 = 0 bacterial / 2 viral; 2011-07-22 = 0 / 8; 2012-04-24 = 8 / 12. All 8 bacterial arrays sit on one date. The plan's power gate (every held-out date: >= 5 bacterial AND >= 10 arrays) fails and cannot pass on this cohort (3 x 5 = 15 bacterial needed, 8 exist). Arm B therefore ends with "GSE40396 dates: underpowered, not run" plus the counted fact that date is confounded with infection class. (2) From the 20 draws on file: direction 1 (AUC 0.849) is below all 20; direction 2 (AUC 0.885) has 1 of 20 draws at or below it, right on the 5 % line, so the Monte-Carlo error of the tail share at 200 draws is about +-1.5 points.

**Deviation from the plan, with reason: 1000 draws (rng seeds 0-999), not 200.** At 200 draws the verdict near the 5 % line would partly be decided by Monte-Carlo noise, not by the data. Seeds 0-199 are a prefix, so the plan's 200-draw count is ALSO reported (descriptive, no verdict word). Output: `results/batch_null1000.csv`. CPU only, EUR 0, at most ~15 CPU-min.

**Pre-registered rule R-batch-1000** (logistic L2, raw, per direction; obs = batch-out AUC already on file in `results/batch_split.csv`; null = N size- and class-matched mixed splits of the same 185 arrays, built by `batch_split.partition` / `check_partition` / `fit_scores` unchanged; EPS 1e-12). Count form: `k_le` = null draws with AUC <= obs (a tie counts as at-or-below); `k_lt` = null draws with AUC < obs.
- **"batch costs AUC"** = `k_le / N <= 0.05` in BOTH directions.
- **"no batch cost"** = `k_lt / N >= 0.25` in BOTH directions.
- **"undecided"** otherwise. Gate B first (batch-AUC >= 0.80, else the batch-split run stop text). Both tie rules lean toward "undecided".
- Within-sample (same null, secondary): R-scale as in the batch-split run (within-sample Delta <= half the raw Delta in both directions, Delta = null mean - obs), read only if "batch costs AUC" fires.
- Descriptive only, no verdict word: both tail shares per direction, 95 % Wilson interval of each tail share, the 200-draw-prefix counts (if they would read differently, the RESULT says so).
- Arm B (GSE40396, leave-one-date-out) runs only if the power gate passes; same rule. No threshold change after data.

**Gates (`--gates-only` writes nothing).** 0a = the batch-split run's (held-out GSE63990 AUCs within 0.005 of `folds_raw.csv`, counts 58/88/12/27, one `batch:` string per array); 0b = seeds 0-19 of the new null reproduce `results/batch_control.csv` (both transforms, both directions) to 1e-6; 0c = the four batch-out AUCs reproduce `results/batch_split.csv` to 1e-6; Gate B = 1.000 again, bar 0.80. Any miss: nothing written, stop.

**Rivals + controls.** Training size / class mix: each draw matches the held-out side's n and class counts, asserted per draw. Batch = sub-population (four sites): NOT controllable; wording stays "batch or sub-population". GEO already batch-corrected: Gate B. Test noise (39 arrays, 12 bacterial small side): the null has the same test size, so the percentile is the control, plus the Wilson interval. The 7 re-run arrays: excluded as on file. Probe choice / NaN imputation read both batches: label-free, equal in test and null (inherited the batch-split run disclosure).
**Honest limit.** The rule is written after seeing the batch-split run's result on the same 185 arrays (both batch-out AUCs and 20 of the 1000 draws are known). It sharpens the null; it is not an independent replication and adds no new patients. The 5 % / 25 % cut-offs are the plan's, fixed now. "Undecided" is a normal outcome. One cohort says nothing about the other four.

**Judge.** Rules committed before code; Gates 0a/0b/0c; tests red on broken code (null not class-matched, date leaking into training, power gate bypassed, rule edges incl. N = 200, inverted tail, seed shifted by one); blind re-read of the RESULT sentences by a fresh subagent that recomputes k_le / k_lt from the CSV itself (same model family: sentence-vs-file, not design).

**Claim-check block (before the run).** Claim: if the rule fires "batch costs AUC", the within-hospital batch lowers held-out AUC beyond what size/class-matched mixing explains, for GSE63990 only. Measured by: batch-out AUC vs the 1000-draw null. Rival: batch = sub-population; test noise; Monte-Carlo noise (1000 draws + Wilson interval). Judge not us: the registered rule + independent recomputation. Files: new only (`results/batch_null1000.csv`, `results/batch_dates_gse40396.csv`, `results/run_batch_deciding.json`); existing results byte-identical; no README / CORE edit; local only, never pushed, public mirror untouched.


## 2026-11-08 — RESULT: deciding batch test — R-batch-1000 = "undecided" (the deciding batch test; PLAN entry 0e529bd, code cbd5d60, a worker branch)

**Question (recap).** The cross-hospital drop (logistic random-CV 0.957 vs leave-one-cohort-out 0.881) is still unattributed. the batch-split run asked whether GSE63990's own two array batches cost AUC and answered "undecided" with only 20 control draws. Here the control is made 50 times larger so the answer is not decided by luck of the draw.

**Design before numbers.** Train on one batch of GSE63990, test on the other (both directions). To know what AUC is "normal" for a test side of that size and class mix, we shuffle all 185 arrays into new train/test splits of exactly the same sizes and class counts (batches mixed on purpose), 1000 times (seeds 0-999), and score each. If the real batch-out AUC sits in the extreme low tail of that shuffle cloud in BOTH directions, the batch matters beyond size and class mix. Rule written before code: "batch costs AUC" = at most 5 % of draws at or below the observed AUC in both directions; "no batch cost" = at least 25 % strictly below in both; else "undecided". Deviation from the plan (stated in the PLAN): 1000 draws instead of 200. Gate 0a (the batch-split run's) passed, max diff 0.0005; Gate 0b: draws 0-19 reproduce `batch_control.csv`, max diff 1.1e-16; Gate 0c: the four batch-out AUCs reproduce `batch_split.csv` exactly; Gate B (can expression predict the batch?) = 1.000, bar 0.80.

**Measured (logistic L2, raw).**
- Direction 1 (train batch 0, test batch 1): observed AUC 0.849. Null mean 0.961. **4 of 1000 draws are at or below it (0.4 %, 95 % Wilson interval 0.16-1.0 %).** Delta = 0.112.
- Direction 2 (train batch 1, test batch 0): observed AUC 0.885. Null mean 0.920. **121 of 1000 draws are at or below it (12.1 %, interval 10.2-14.3 %); 120 strictly below.** Delta = 0.035.
- First 200 draws only (the plan's literal count): direction 1 = 0 of 200; direction 2 = 25 of 200 (12.5 %). Same word as at 1000.
- Arm B (GSE40396, hybridization dates): labelled bacterial/viral arrays 2010-08-03 = 0/2, 2011-07-22 = 0/8, 2012-04-24 = 8/12. Power gate (>= 5 bacterial and >= 10 arrays per date) fails on two dates and cannot pass on this cohort. **"GSE40396 dates: underpowered, not run."** No AUC was looked at. Counted fact: date is confounded with infection class there (all 8 bacterial arrays sit on one date).

**Pre-registered reading, applied as written: "undecided".** Direction 1 meets the "costs AUC" bar (0.4 % <= 5 %); direction 2 does not (12.1 %); neither direction pair meets "no batch cost". R-scale: not evaluated.

**In plain words.** In one direction the real batch split is clearly worse than almost any random split of the same size (4 in 1000). In the other direction it is only a little worse than typical and about one random split in eight does as badly or worse. With 1000 draws the interval for that 12 % is tight (10-14 %), so "undecided" is now a real finding, not noise from too few draws: the batch effect is not symmetric. Only one direction shows it, and the weak direction is the one where the test side is the large batch (148 arrays), so a large effect would be harder to hide there. A plausible but UNTESTED reading: batch 1 (39 arrays, 12 bacterial) is a small or different sub-population, which hurts when it is the thing being predicted from but not when it is predicted.

**Claim-check.** Claim: batch costs AUC beyond size/class mix. Measured: fires in one direction, not the other. Rival 1, batch = sub-population (four sites): not controllable, wording stays "batch or sub-population". Rival 2, test noise: the null has the same test size; Wilson interval shown. Rival 3, Monte-Carlo noise: 1000 draws, 200-draw prefix agrees. Honest limit: the rule was written after seeing the batch-split run's result on the same 185 arrays (not an independent replication; adds no patients); one cohort says nothing about the other four. Judge not us: rule committed first (0e529bd), Gates 0a/0b/0c, tests red on broken code (4 mutations of tie rule, power gate, 5 % edge, 25 % edge each made tests fail), and a blind fresh subagent recomputed k_le / k_lt / null means / date counts from the CSVs: all 7 sentences MATCH (same model family: sentence-versus-file check, not a design check). Suite 211 passed (191 + 20 new); all 55 pre-existing `results/` files md5-identical before and after.

**Belief.** Unchanged in substance: within-hospital batch is a partial, one-directional candidate for the drop, not an established cause, and the drop stays unattributed. The README causes-table line is NOT edited here (plan rule): the lab-batch row can move from "undecided (20 draws)" to "undecided (1000 draws; costs AUC in one direction only)" after the author reads this.

**Where.** `results/batch_null1000.csv`, `results/batch_dates_gse40396.csv`, `results/run_batch_deciding.json`, `run_batch_deciding.py`, `src/batch_deciding.py`, `src/batch_deciding_verdict.py`, `tests/test_batch_deciding.py`. A worker branch, local only (repo has no remote; never pushed; public mirror untouched). Runtime 11 CPU-min, EUR 0.

## 2026-11-09 — Deciding batch test and causes-table update merged

Process only. The deciding batch test (above) and a docs update were merged into the private working repo's main branch; full suite 211 passed. The docs update (README causes table and `reports/2026-11-08-writeup.html`) now carries the deciding-test numbers: direction 1 AUC 0.849 vs null mean 0.961, 4 of 1000 draws at or below; direction 2 AUC 0.885 vs null mean 0.920, 121 of 1000 (12.1 %, interval 10.2-14.3 %); verdict still undecided; GSE40396 dates underpowered, not run. **Belief.** Unchanged: the drop stays unattributed; within-hospital batch is a partial, one-directional candidate.


## 2026-11-09 — Results of 30.10.–09.11. carried into this public copy

**Question.** This copy stood at the 31.10. state: no Result 4, no whole-blood check, no batch tests, no write-up. Can the newer work come over without private material and without any number changing on the way?

**Done.** Code, tests, config and result files for the LightGBM model-seed check (30.10.), the transfer matrix, the 10-seed Arm B run, the GEO age check, the whole-blood check, the batch-split run and the deciding batch test copied from the private working repo; README replaced by its reading-story version; `reports/2026-11-08-writeup.html` added; the diary entries above appended. Same privacy scrub as for the first copy (names, session and task-card ids, private links and folders removed; pronouns for the author → they). Card ids in three result JSONs (`card` label field) became plain run names; no number in any result file changed. Three internal HTML plan/result pages were left out; their content is in the entries above. CORE.md unchanged: its §7 already reads as the private one does now.

**Measured.** Every number new to this copy was re-computed from the committed result files before the carry-over (chip-family difference −0.061, size-matched −0.062, within-sample −0.034, LightGBM −0.114 / −0.111 / −0.108; class-weighting +0.0005 / +0.0033 (SD 0.026); whole-blood logistic 0.857 vs 0.866, gap 0.081 vs 0.076; batch split 0.849 / 0.885; deciding test 4 / 1000 and 121 / 1000): all match. One stale count fixed: README and write-up said 172 / 168 tests; the suite now has 211 (211 passed in this copy). A blind fresh agent (same model family, so a text-versus-file judge, not a design judge) then read the whole diff: no leak; every Result 4 and causes-table number reproduced from `results/`; four sentences softened on its advice ("we found no public code", pre-registration exceptions and the private location of the pre-registration commits stated, "independent" dropped from the AI code reviews, "help one and not the other"), and the Reproduce list now names every runner whose numbers are quoted. Side note from that read, not yet checked: GSE6269 (tagged "children") has 3 arrays with a GEO age of 18 or more.

**Belief.** Unchanged: about 7 to 12 AUC points of the headline are hospital-specific; the per-cohort drop is a pattern without a named cause; within-hospital batch is a partial, one-directional candidate.

**In plain terms.** The public edition of a lab notebook gets its new chapters: names blacked out, every number checked against the raw tables, printed only after a proofread.

**Where.** This copy: README.md, LABBOOK.md, `reports/2026-11-08-writeup.html`, `results/`, `src/`, `tests/`, `tools/`.

## 2026-11-10 — GSE6269 age tag checked (wording only, no number changes)

**Question.** GSE6269 is tagged "children", but the blind read of the last carry-over saw 3 arrays with a GEO age of 18 or more. Is the tag wrong, and does any reported number move?

**Measured.** GPL96 series matrix, 97 arrays (91 labelled + 6 healthy). 88 of the 91 labelled arrays are under 18; 3 are adults: GSM173201 (S. aureus, 18 y), GSM173221 (influenza A, 24 y), GSM173268 (influenza A, 36 y). The paper agrees independently: Ramilo 2007 (PMC1801073) Table 1 footnote gives S. aureus 3 months to 18 years and influenza A 3 weeks to 36 years, and its patient table lists one influenza patient at 36 y.

**Effect on numbers: none.** The age tag feeds only the age test of the transfer matrix, which was already "not testable" (GSE63990 has no ages). Every AUC uses the class labels, not age.

**Claim-check.** Claim: "children" is imprecise for GSE6269 (3 of 91 adults); no reported number changes. Rival: a GEO data-entry error, ruled out by the paper's own table. Judge: the paper; otherwise self-assessed. "Children" is the paper's own cohort-level word ("pediatric patients"), so the tag is called imprecise, not wrong.

**Changed.** README and `reports/2026-11-08-writeup.html`: GSE6269 age group "mostly children (3 of 91 aged 18–36)" and "the mostly pediatric … cohort"; `config/cohorts.json` `age_source` note lists the 3 sample ids (no value read by code changed). Tests: 211 passed.

**Belief.** Unchanged.

## 2026-11-10 — PLAN, pre-registered before any read: Bodkin 2022 side-by-side

**Question.** Is our random-CV AUC per cohort in line with what the largest published comparison (Bodkin et al. 2022, Genome Medicine, 28 signatures, 51 datasets) measured on the same GEO series? This is the last open CORE s7 "MAY propose" item. Literature only: no run, no GPU, EUR 0. Committed BEFORE the supplement is opened for any number of our five series.

**Why Table 2 is not used.** Table 2 is one row per signature (weighted mean AUC, sens/spec/DOR), no dataset column. Per-dataset numbers are in Additional file 6 = Table S9 (XLSX): median AUC and IQR (Q1-Q3) across the signatures, per dataset and per task. Min and max are not published as numbers (box plots in Fig S3 only, no GSE labels), so the original "above max / inside min-max" rule cannot be applied. This replacement (the author's approved option A) is recorded openly.

**Same series** = GEO accession match only. GSE6269 counts as "in Bodkin" if either of its two Bodkin rows is present; row `-2` (GPL96, our platform) is the one the rule uses.

**Our number** per series = `rcv_auc_in_cohort` (random-CV predictions of the pooled model scored inside each cohort) from `results/gapdiag_cohorts.csv` (transform=raw, model=logistic_l2), with `loco_auc` beside it. Sensitivity lines, never the verdict: the pooled headline random-CV mean AUC (0.957, `summary_raw.csv`) and Bodkin's viral-vs-non-viral task.

**Task column.** Bodkin's bacterial vs non-bacterial (positive class = bacterial, like ours). Viral vs non-viral is copied as a sensitivity line only.

**Rule (A, verbatim).** n = number of our 5 series found in Table S9. If n < 3 -> NOT-COMPARABLE. OPTIMISTIC-FLAG if our random-CV AUC is above Bodkin's Q3 by more than 0.03 for EVERY overlapping series. CONSISTENT if it sits inside [Q1, Q3] for MORE THAN HALF of the overlapping series (3 of 5 when n = 5). Anything else (including values below Q1) -> NOT-COMPARABLE, reported per series, nothing claimed. If the primary and the sensitivity readings give different outcome words, the RESULT says so and claims no more than the weaker word.

**Known weakening, stated up front.** Bodkin's negative class includes healthy / SIRS subjects (ours is bacterial vs viral only); their models are per-dataset lasso refits with nested leave-one-out; they use 29 signatures (28 plus one composite); Q3 is lower than the max, so the flag is easier to trigger than the card intended.

**Kill outcome.** OPTIMISTIC-FLAG = recommend (not start) a pipeline check before further outside use (CORE s3 kill 1, "the measurement was a bug").

**Judge not us.** The published table, plus a blind subagent re-read of every copied cell (gets the bodkin_* columns and the xlsx URL, not our numbers, the rule or the verdict).

## 2026-11-10 — RESULT: Bodkin 2022 side-by-side = NOT-COMPARABLE (PLAN commit `b11104b`)

**Question (recap).** Is our random-CV AUC per cohort in line with what Bodkin et al. 2022 measured on the same GEO series? Rule A was written in the PLAN commit before any number was read; it was not edited.

**Download log.** PMC8858657 Additional file 6 (Table S9), 18,337 bytes, sheet `Table S9` 76x19, downloaded 10.11.2026 by the author by hand into `data/papers/bodkin_s9.xlsx` (git-ignored), sha256 `3cfd0eaee695dfb81dc01324746114063d9fa5249345c105de076fac1aa12da2`.

**Design.** All 5 series are in Table S9 (GSE6269 via row -2, GPL96; row -1 listed only in the note, not used). Bodkin publishes median and IQR across signatures per dataset, not min/max. Their design: per-dataset lasso refits, nested leave-one-out, negative class includes healthy/SIRS, ours is bacterial vs viral only.

**Measured** (bacterial task; ours = raw logistic_l2 `rcv_auc_in_cohort` / `loco_auc` from `gapdiag_cohorts.csv`):

| series | Table S9 row | Bodkin median [Q1-Q3] | ours random-CV | random-CV minus Q3 | ours held-out (LOCO) |
|---|---|---|---|---|---|
| GSE63990 | 40 | 0.899 [0.844-0.909] | 0.959 | +0.050 | 0.926 |
| GSE42026 | 33 | 0.780 [0.650-0.837] | 0.893 | +0.056 | 0.832 |
| GSE60244 | 35 | 0.762 [0.682-0.793] | 0.894 | +0.101 | 0.843 |
| GSE40396 | 31 | 0.868 [0.820-0.910] | 0.938 | +0.027 | 0.864 |
| GSE6269 | 39 | 0.845 [0.818-0.876] | 0.988 | +0.112 | 0.938 |

**Rule applied mechanically.** Primary: n=5; flag needs every series above Q3 by more than 0.03, but GSE40396 is +0.027; none sits inside [Q1,Q3] -> **NOT-COMPARABLE**. Sensitivities: viral-vs-non-viral columns -> NOT-COMPARABLE; held-out AUC vs the bacterial Q3 -> NOT-COMPARABLE; the pooled headline 0.957 against each series' Q3 -> OPTIMISTIC-FLAG (by +0.047 to +0.164). Readings disagree, so per the PLAN the claim is no stronger than the weaker word: NOT-COMPARABLE. The primary margin for GSE40396 is 0.003 under the line; no rule change made.

**Claim-check.** Claim: our optimistic number is not an artefact of our pipeline. Measured: nothing supports that claim; ours sits above Bodkin's upper quartile for all five series, but their numbers come from a different design (rival: not comparable), so no verdict either way. Judge not us: published Table S9 plus two independent re-reads of every copied cell (blind subagent: all cells match; own accession-keyed lookup: all match). Outcome word: NOT-COMPARABLE.

**Belief.** Unchanged. Recommendation only (not started): the consistent direction (ours higher in random-CV for 5 of 5) does not prove a bug, since Bodkin's lasso signatures are not tuned for our narrower task; if the author wants more, the pipeline check of CORE s3 kill 1 is the next task, not part of this one.

**Where.** `results/bodkin_table2_comparison.csv`, README section 8, `tests/test_bodkin_comparison.py`; suite 213 passed (211 + 2 new).

## 2026-11-12 — PLAN, pre-registered before any code: label-shuffle null of the raw random-CV headline

**Question.** Does the pooled random-CV number (logistic_l2, raw, 458 arrays x 12,040 genes, fold seed 0) contain something that is not biology? Two nulls, 50 shuffles each (seeds 0-49), through the same code path as the headline (`splits.random_splits`, `train._fit_score`, `make_models()["logistic_l2"]`): (1) GLOBAL shuffle of bacterial/viral over all 458 arrays; (2) WITHIN-COHORT shuffle (labels permuted only inside each cohort, so each cohort keeps its class mix). CPU only, EUR 0. This check is not on the CORE s7 list; the author approved it as a one-off item on 12.11. No CORE edit.

**Why.** The Bodkin side-by-side (10.11.) ended NOT-COMPARABLE; our random-CV AUC sat above Bodkin's Q3 for 5 of 5 series. The RESULT named this pipeline check as the next task.

**Reading rule (verbatim from the task, thresholds unchanged).**
- LEAK: global-shuffle mean AUC > 0.55, or any single shuffle > 0.65. Then the headline is not quoted anywhere until the cause is found (CORE s3 kill 1, "the measurement was a bug").
- CLEAN: global mean within 0.45-0.55 and no single shuffle > 0.65.
- SHORTCUT SIZE (reported, not a kill): within-cohort mean AUC as a number; SMALL if < 0.60, LARGE if >= 0.70, else MEDIUM.

**Three gaps closed here, no threshold changed.** (a) One shuffle's AUC = the mean of its 5 fold AUCs (same definition as `summary_raw.csv`); p95 = numpy linear 95th percentile over the 50. (b) A global mean below 0.45 (with no shuffle above 0.65) is neither LEAK nor CLEAN: word ODD-LOW, no leak claim, report and stop. CLEAN includes 0.45 and 0.55. (c) The RESULT quotes the on-file logistic held-out AUC (0.8807) beside the random-CV null (CLAUDE.md rule 4); no held-out null is run.

**Planning-probe facts (read-only, before this entry).** Logistic random-CV at fold seed 0 = 0.9566006614665998 in 3.7 s vs 0.9566006614666 on file (diff 2e-16), so Gate 0 at 1e-12 should pass. Class mix (viral/bacterial): GSE63990 115/70, GSE42026 41/18, GSE60244 71/22, GSE40396 22/8, GSE6269 18/73. A score that knows only the cohort already reaches AUC 0.7004 on these labels (labels only, no model); this is the same for every within-cohort shuffle. So LARGE (>= 0.70) is reachable only if the model recovers nearly all of the cohort-only shortcut. Printed as a reference line, not a rule.

**Gate 0.** With labels unshuffled, logistic random-CV and held-out reproduce `results/summary_raw.csv` to 1e-12, else exit 1 before any shuffle, no output files.

**Claim-check.** Claim: the random-CV headline reflects host-response signal, not a pipeline leak or cohort shortcut. Measured: random-CV AUC under both shuffles, 50 seeds each, same code path. Rival 1 (label leak via preprocessing): control = global shuffle. Rival 2 (cohort identity + class mix): control = within-cohort shuffle. Judge not us: a fresh sub-agent recomputes verdict words from the CSV. Limit: shuffling is per array, so it cannot see a same-patient-twice leak (the 7 known GSE63990 re-runs are excluded; others would be unknown).

**Outputs.** `results/label_null_raw.csv` (100 rows), `results/run_label_null.json`. Run on a worker branch.

## 2026-11-12 — RESULT: label-shuffle null of the raw random-CV headline = CLEAN, cohort shortcut MEDIUM (PLAN commit `81252f2`)

**Question.** Does the pooled random-CV AUC (logistic_l2, raw, 458 arrays x 12,040 genes, fold seed 0) contain something that is not biology? Rival 1: a label leak through the pipeline. Rival 2: cohort identity plus per-cohort class mix.

**Design.** Same code path as the headline. Gate 0 first: unshuffled random-CV differs from `summary_raw.csv` by 2.2e-16 and held-out by 0.0 (limit 1e-12): passed. Then 50 shuffles (seeds 0-49) per null, folds re-stratified on the shuffled labels. One shuffle's AUC = mean of its 5 fold AUCs. Reading rule: the PLAN entry above, unchanged.

**Measured** (`results/label_null_raw.csv`, 100 rows; `results/run_label_null.json`):

| null | mean | p95 | max | min | sd |
|---|---|---|---|---|---|
| global shuffle | 0.502 | 0.543 | 0.574 | 0.428 | 0.031 |
| within-cohort shuffle | 0.668 | 0.699 | 0.701 | 0.632 | 0.020 |

Reference line (labels only, no model): a score that knows only each array's cohort reaches AUC 0.7004. Both-schemes line (CLAUDE.md rule 4): headline random-CV 0.957 and on-file logistic held-out (leave-one-cohort-out) 0.881.

**Verdict by the pre-registered rule.** Global mean 0.502 is inside 0.45-0.55 and no global shuffle is above 0.65 (max 0.574): **CLEAN**. Within-cohort mean 0.668 (not < 0.60, not >= 0.70): **MEDIUM**. In plain words: with the labels scrambled over all arrays the model scores like a coin flip, so no label leak shows up. With labels scrambled only inside each hospital, the model still scores 0.668, because each hospital has its own mix of viral and bacterial; that is about 95 % of the most a cohort-only score could get (0.700). So a real part of the 0.957 headline can come from knowing the hospital, not from the host response. This matches the direction of the held-out gap (0.957 vs 0.881).

**Claim-check.** Claim: the headline is not a pipeline label leak. Supported only for the label-leak route the global shuffle can see. Rival 2 is NOT ruled out; it is measured (MEDIUM), and the null does not say how much of 0.957 is shortcut, only how much a shortcut alone can reach (about 0.67). Judge not us: a fresh sub-agent given only the CSV and the rule recomputed every number and both words: all MATCH (it also noted 40 of 50 within-cohort shuffles are above 0.65; the 0.65 cap applies to the global null only). Limits: shuffle is per array, so a same-patient-twice leak is invisible (the 7 known GSE63990 re-runs are already excluded); one fold seed (0), one model, random-CV only; no held-out null run.

**Belief.** The headline random-CV number is not a label-leak bug, so it may stay quoted under CLAUDE.md:5-7. Cohort shortcut is a real, medium-size contributor and stays part of the gap story. CORE unchanged.

**Disclosures.** `commit` in the JSON reads `81252f2-dirty`: the run started before the code was committed (the same files were committed right after as `6c98d2f`, unchanged). No run was repeated or patched by hand. The matrix was cached in the job tmp folder (not in the repo). Tests: red-first checked by breaking the code three ways (no-op shuffle, shuffle crossing cohorts, Gate 0 bypassed); suite 234 passed (215 + 19 new).

**Where.** `src/label_null.py`, `run_label_null.py`, `tests/test_label_null.py`, `results/label_null_raw.csv`, `results/run_label_null.json`.

## 2026-11-13 — PLAN, pre-registered before any code: does the headline lean on the cohort shortcut? within-cohort-pair AUC vs pooled

**Question.** The 12.11. label-shuffle null said the within-cohort shuffle reaches 0.668 (cohort-only ceiling 0.700): "a real part of the 0.957 headline can come from knowing the hospital". The audit table (row 14) said the opposite way: within-cohort-pair AUC about equal to pooled. The two were never put side by side. A pair of patients from the SAME cohort cannot be told apart by knowing the cohort, so if the within-pair AUC equals the pooled AUC, the headline does not lean on the shortcut; the null showed only how far a shortcut COULD reach. Raw random-CV, 458 arrays, 5 folds, fold seeds 0-9, logistic_l2 and LightGBM. CPU only, EUR 0, no network. Not on the CORE s7 list; the author approved it as a one-off (13.11., option B: both extras below).

**Reading rule (verbatim from the task, thresholds unchanged).**
- NOT-USED: mean over 10 seeds of (pooled minus within-pair) is <= 0.01 for logistic_l2. Then the 12.11. belief sentence is corrected: "a cohort shortcut could reach about 0.67, but the model's ranking inside each cohort is as good as its pooled ranking, so the headline does not lean on it".
- USED: mean (pooled minus within-pair) >= 0.02 for logistic_l2. Then the 12.11. sentence stands, now with a size.
- In between: UNCLEAR, both numbers reported, sentence softened to "not shown either way".
- LightGBM is reported beside it, direction only (no band word).

**Expected numbers (today's file, NOT row 14).** Row 14 (logistic 0.954 within / 0.953 pooled, LightGBM 0.951 / 0.954) comes from the old 463-array `gapdiag_folds.csv`; master has carried the 458-array re-run since 21.10. Read-only probe of today's file, seed 0, pair-weighted over 5 folds: logistic pooled 0.9565, within 0.9542, cross 0.9571, D = +0.002; LightGBM pooled 0.9543, within 0.9496, cross 0.9557, D = +0.005. So the likely word is NOT-USED. If it comes out otherwise, that is the finding.

**Gaps closed here, no threshold changed.**
(a) Aggregation: per test fold `pair_auc` (pairs counted only inside one test fold). Over the 5 folds pair-weighted: pooled weighted by all pairs of the fold, within by its same-cohort pairs, cross by its cross-cohort pairs. Per seed D = pooled minus within from these. The verdict uses the plain mean of D over seeds 0-9.
(b) Fit path: the E1 path (`gap_experiments._fit`, `predict_proba`, `pair_auc`), the only one that reproduces `gapdiag_folds.csv` (the `train._fit_score` path differs on logistic fold 4 by 5e-4). LightGBM model seed 0. Fold seeds via `splits.random_splits(y, seed=s)`.
(c) Gate 0: seed 0, both models, every fold: pooled / within / cross reproduce the raw random_cv rows of `gapdiag_folds.csv` to 1e-9 and the integer pair counts are equal. Any failure: exit 1, no files.
(d) Instrument check (positive control 1): logistic, 10 within-cohort shuffles (`label_null.shuffle_within_cohort` seeds 0-9 plus `check_shuffle`), folds re-stratified on the shuffled labels at fold seed 0, same aggregator. Expected within-pair mean about 0.50, pooled about 0.67. STOP if the 10-shuffle mean within-pair AUC is outside 0.45-0.55 (edges pass): write only the JSON with verdict "NO VERDICT (instrument)", no main read.
(e) Positive control 2 (the author's extra 1): plant a KNOWN hospital offset into the real out-of-fold logistic scores. Offset of a patient = k x sd x z, where sd = standard deviation of that seed's out-of-fold logits, z = the patient's cohort bacterial share, standardised over patients (mean 0, sd 1), so high-bacterial hospitals get higher scores (a shortcut by construction). k in {0.1, 0.25, 0.5, 1.0}. Within-pair AUC must stay unchanged (to 1e-9) because an offset is constant inside a cohort. PASS = mean D at k = 1.0 is >= 0.02 (the method CAN see a shortcut at this signal strength). I also report the smallest k whose mean D reaches 0.02. A FAIL does not change the verdict word, but then NOT-USED is written as "not detected by a method whose sensitivity here is unproven".
(f) Interval (the author's extra 2): patient-level bootstrap of the 10-seed mean D (logistic and LightGBM). Each replicate draws a bootstrap count per patient (multinomial over 458 patients), shared by all 10 seeds; within each test fold every pair is weighted by the product of its two patients' counts; the same pair-weighted aggregation. B = 2000, numpy seed 0, percentile 2.5 / 97.5. Descriptive only: the verdict word is fixed by the rule on the point mean. If the interval crosses 0.01 or 0.02 it is reported as such.
(g) Limits written in advance. D is a NET number: a shortcut gain on cross-cohort pairs could be cancelled by hospital offsets that hurt cross-cohort pairs (platform / batch shifts; the 0.076 held-out gap shows such shifts exist). NOT-USED therefore means "the pooled number is not above the within-hospital number", which is what the headline claim needs, not "the model carries no hospital information". Ten fold seeds average away fold assignment, not the patient sample (same 458 arrays); within pairs are only about 23 % of all pairs. The instrument check works in a no-biology regime (D about 0.17); control 2 is what speaks to D about 0.02 on top of a real 0.95 signal.

**Kill line.** USED keeps the 12.11. sentence and sizes it. NOT-USED retracts it: the null measured a ceiling, not a contribution. Neither is a CORE s3 kill; both change the gap story's wording.

**Judge not us.** A fresh sub-agent gets only the CSV and the reading-rule text and returns the logistic 10-seed mean of D, the verdict word and the LightGBM direction; MATCH / MISMATCH is pasted in the RESULT.

**Outputs.** `results/shortcut_pairs_raw.csv` (20 rows), `results/run_shortcut_pairs.json`, `results/shortcut_pairs_planted.csv`. Run on a worker branch.

**In plain terms.** Imagine judging which of two students is the stronger reader. If you only know the school, and one school is famous for strong readers, you can often guess right without hearing either student read. To check whether our model is secretly guessing from the school (the hospital), we compare only pairs from the SAME school, where the school name cannot help, and see if the model is still just as good. If it is, the high score comes from the reading, not from the school's reputation.

## 2026-11-13 — PLAN ADDENDUM (before the real run, no result read): control 2 redesigned, reading rule untouched (PLAN commit `65a9d0b`)

**What happened.** While testing the code I probed control 2 on the real logistic scores, fold seed 0 only (the planted offsets, NOT the main read; no verdict quantity was looked at beyond that seed's known D = +0.0022). Result: every planted hospital offset LOWERS pooled AUC (k 0.1: D -0.0003; k 0.25: -0.006; k 0.5: -0.019; k 1.0: -0.042) and a label-aware search for the best possible per-hospital constant offset (start at zero, coordinate grid, objective = pooled AUC) lifts D only to +0.0035. Reason: pooled AUC is already 0.957, so about 0.043 of headroom is left, and a constant per hospital can fix few of the remaining errors but breaks many correct pairs. So control 2 as written in the PLAN (PASS = mean D at k = 1.0 is >= 0.02) cannot pass on this model, for a reason that has nothing to do with whether the method works.

**Change (control 2 only; the reading rule, thresholds, gate, instrument check and bootstrap are unchanged).** Control 2 now reports three things, all descriptive:
(1) **Invariance + sensitivity.** Planted offsets (same k grid, plus the k = 0 base) leave the within-pair AUC unchanged (to 1e-9) while pooled moves. PASS = |mean D at k = 0.5 minus mean D at k = 0| >= 0.02 (the metric registers a planted hospital offset of that size).
(2) **Headroom bound.** For each fold seed, the best label-aware per-hospital offset (start at zero, 41-point grid on -1.0..+1.0 x sd of the out-of-fold logits, 4 coordinate sweeps, objective = pooled AUC through the same aggregator) and the D it reaches: D_max(seed), mean over seeds. It is an upper bound for what any hospital-only additive shortcut can add to this model's pooled AUC.
(3) **Reading the verdict with the bound.** If mean D_max < 0.02, the USED line cannot be reached by a hospital offset on this model; then NOT-USED is written as "no sign of it, and little room for it" instead of as a strong refutation, and the RESULT says that the control the author asked for ("show D can pass 0.02") is impossible here and why. If mean D_max >= 0.02, the USED line is reachable and NOT-USED stands as written. The verdict WORD is still set only by the rule on the point mean.

## 2026-11-13 — RESULT: the headline does NOT lean on the cohort shortcut = NOT-USED (PLAN `65a9d0b`, ADDENDUM `892e79e`, run on code commit `ffed5fd`)

**Question, in plain words.** The 12.11. null said a model that only knows the hospital already scores about 0.67, so "part of the 0.957 may be hospital mix, not biology". The audit table pointed the other way. Test: compare the model's AUC on ALL (bacterial, viral) pairs with its AUC on pairs where both patients come from the SAME hospital. A same-hospital pair cannot be separated by knowing the hospital. If the two AUCs are equal, the headline does not rely on the shortcut.

**Design before numbers.** Raw random-CV, 458 arrays, 5 folds, fold seeds 0-9, logistic_l2 and LightGBM, E1 fit path, pairs counted only inside one test fold, pair-weighted over folds. D = pooled minus within-pair. Rule fixed in the PLAN: NOT-USED if the 10-seed mean D (logistic) is <= 0.01, USED if >= 0.02, else UNCLEAR.

**Gates.** Gate 0 (seed 0, both models, every fold, vs `gapdiag_folds.csv`): max diff 1.1e-16, 0 pair-count mismatches: passed. Instrument check (10 within-hospital label shuffles, logistic): within-pair AUC mean 0.520 (sd 0.060, per shuffle 0.42-0.60), inside the 0.45-0.55 window; pooled mean 0.671 (matches the 12.11. null to about 5e-4 per shuffle). So with nothing but the hospital to go on, pooled reads 0.67 and the within-pair AUC reads about chance (0.52): the metric separates the two.

**Main numbers (10-seed means).**

| model | pooled | within | cross | D = pooled - within | sd of D over seeds |
|---|---|---|---|---|---|
| logistic_l2 | 0.9561 | 0.9547 | 0.9566 | **+0.0014** | 0.0034 |
| LightGBM | 0.9521 | 0.9483 | 0.9533 | +0.0038 | 0.0050 |

Per-seed logistic D runs from -0.0038 to +0.0078; every seed is under 0.01. Patient-level bootstrap interval for the logistic mean D (B = 2000): -0.0079 to +0.0126 (LightGBM -0.0061 to +0.0146). The upper end touches the 0.01 line but stays under 0.02; the verdict word is set by the point mean, as pre-registered.

**Verdict: NOT-USED** (mean D 0.0014 <= 0.01). LightGBM direction: pooled slightly above within (+0.004), small. Judge not us: a fresh sub-agent that saw only the CSV and the rule text got logistic mean D +0.001414, LightGBM +0.003767, NOT-USED, "above zero": **MATCH** with the JSON; its recomputed D column matched on all 20 rows.

**The 12.11. belief sentence: CORRECTED** (per the rule): "a cohort shortcut could reach about 0.67, but the model's ranking inside each cohort is as good as its pooled ranking, so the headline does not lean on it". The 12.11. null measured a ceiling, not a contribution.

**Control 2 (the author's extra 1) and what it showed.** As pre-registered (planted offset, pass at D >= 0.02 for k = 1.0) the control could not pass; found before the real run and recorded in the ADDENDUM: pooled is 0.957 already, so every planted hospital offset LOWERS pooled. Redesigned control, results: the within-pair AUC does not move at all (shift 0.0) while pooled moves, and mean D shifts by 0.024 at k = 0.5 (D: k=0 +0.0014, k=0.1 -0.0010, k=0.25 -0.0080, k=0.5 -0.0223, k=1.0 -0.0441): sensitivity PASS (the metric registers a planted hospital offset of that size, in the downward direction). Headroom bound: the best label-aware per-hospital constant offset lifts D to only +0.0027 on average over seeds, far below 0.02. So (a) "show D can pass 0.02" is impossible on this model by any hospital-only offset, which is a property of the headroom, not a failure of the method; (b) per the ADDENDUM, NOT-USED is read as **no sign of it, and little room for it**, not as a strong refutation: a hospital-only additive shortcut cannot add more than about 0.003 on top of a model already at 0.957.

**Both schemes (CLAUDE.md rule 4).** Random-CV 0.957 (logistic, `summary_raw.csv`) next to leave-one-cohort-out 0.881: the 0.076 gap stays. Rival 1 (within-hospital pairs are an easier subset for other reasons): in-cohort random-CV AUC / held-out AUC from `gapdiag_cohorts.csv`: GSE63990 0.959 / 0.926, GSE42026 0.893 / 0.832, GSE60244 0.894 / 0.843, GSE40396 0.938 / 0.864, GSE6269 0.988 / 0.938. The model ranks well inside each hospital (0.89-0.99); it is higher than held-out in every cohort, which is the gap, not the shortcut.

**Claim-check.** Claim: the 0.957 random-CV headline does not depend on the cohort/class-mix shortcut. Measured: pooled vs same-cohort-pair AUC, 10 fold seeds, raw. Rival 2 (seed 0 lucky): 10 seeds, sd 0.003. Judge not us: blind sub-agent, MATCH. Limits: D is NET (offsets that help some cross-cohort pairs and hurt others cancel), so NOT-USED means "pooled is not above within", not "no hospital information in the model"; ten fold seeds do not vary the 458 patients (the bootstrap interval is the patient-sample error range and reaches +0.013); within pairs are about 23 % of all pairs; this does NOT say the model transfers to new hospitals (0.881).

**Disclosures.** Row 14's numbers (0.954 / 0.953) came from the old 463-array file; the planning probe used today's 458-array file (seed 0 D = +0.002, reproduced). The E1 fit path differs from the headline `train._fit_score` path by about 5e-4 on one logistic fold (seed-0 pooled 0.95647 vs headline 0.95660); not chased, cannot move a verdict. The control-2 redesign was made after a probe on fold seed 0 and before the main run (ADDENDUM). Red-first checked: breaking the code four ways (within/cross masks swapped, plain fold mean, `<` instead of `<=`, Gate 0 skipped) each turned tests red. Suite: 255 passed (234 + 21 new).

**Where.** `results/shortcut_pairs_raw.csv` (20 rows), `results/shortcut_pairs_planted.csv`, `results/run_shortcut_pairs.json`, `src/shortcut_pairs.py`, `src/shortcut_extras.py`, `run_shortcut_pairs.py`, `tests/test_shortcut_pairs.py`.

**Belief.** Updated: the headline's 0.957 does not come from knowing the hospital. The gap to the held-out 0.881 still stands and is the real story.


## 2026-11-14 — Results of 10.–13.11. carried into this public copy

**Question.** This copy stood at the 10.11. state. Can the Bodkin side-by-side, the label-shuffle null and the same-hospital-pair check come over without private material and without any number changing?

**Done.** Code, tests and result files for the three checks copied from the private working repo; the PLAN / ADDENDUM / RESULT entries above appended; README: a short section "Two checks on the random-CV headline itself", the Bodkin side-by-side paragraph in section 8, two Reproduce lines, test count 211 → 253. Same privacy scrub as before (names, task-card ids, branch names and private-repo notes removed; author → "the author"). An internal German status page was left out. The 08.11. write-up page is a dated snapshot and was not changed.

**Measured.** Re-computed from the committed files in this copy: label shuffle global mean 0.502 (max 0.574), within-hospital 0.668; same-hospital pairs mean D logistic +0.0014, LightGBM +0.0038; headroom bound 0.0027. All match the RESULT entries. `pytest`: 253 passed.

**Belief.** The 0.957 headline does not come from knowing the hospital; the 0.076 gap to a new hospital (0.957 vs 0.881) stands and is the result.

**Where.** README.md sections 3 and 8, LABBOOK.md, `results/`, `src/label_null.py`, `src/shortcut_pairs.py`, `src/shortcut_extras.py`, `run_label_null.py`, `run_shortcut_pairs.py`, `tests/`.
