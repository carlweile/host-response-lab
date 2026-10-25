# host-response-lab — CORE

_Constitution of this lab. Changes only by the author's decision. Diary = LABBOOK.md._
_Public-copy note: line references like `README.md:36` or `LABBOOK.md:82-84` below point to those files as they were when this was drafted in the private working repo; they may not match the files in this copy. §6 is a dated snapshot; current numbers are in README.md._
_Drafted 2026-10-11 from README.md, LABBOOK.md, CLAUDE.md, config/cohorts.json, run_pipeline.py, results/summary_raw.csv, results/summary_ranknorm.csv, results/folds_raw.csv, results/run_raw.json, results/label_review.md. No docs/ folder exists. Red-teamed 2026-10-11 (every citation re-checked against the repo)._

## 1. Question
Verbatim (README.md:3-5): "when a published host-response signature says 'this infection is viral, not bacterial' with an AUC of 0.93, how much of that number survives being tested on a hospital cohort the model has never seen?"

Gloss: AUC = a 0.5–1.0 score of how well a model separates the two classes (1.0 = perfect, 0.5 = coin flip). "Cohort" = one hospital's patient set, one GEO series.

## 2. Success looks like
**ACCEPTED 2026-10-18** by the author (LABBOOK, 2026-10-18 evening entry). Drafted because the repo stated no "done when"; README.md:36 only says "The headline output is not an AUC. It is the *difference* between two AUCs.").

Note: the side-by-side table itself ALREADY exists (`results/summary_raw.csv` has random-CV, leave-one-cohort-out and gap per model). So "a table exists" is not success. Success is the table being **trustworthy**. Done when ALL four hold:
1. **Within-sample normalisation control run** (each patient ranked against their own genes only, no other patient involved — LABBOOK.md:82-84). Artifact: a new `results/summary_withinsample.csv` (now on file) in the same format as `summary_raw.csv`, produced by one new `run_pipeline.py --train` flag, plus a LABBOOK entry stating the gap under this leak-free setting for both models. The gap may survive or vanish — EITHER is a result; what matters is that it was measured without the test set leaking.
2. **Labels verified for all five cohorts.** `expected_counts` in `config/cohorts.json` filled for GSE60244, GSE40396, GSE6269 (today null; config/cohorts.json:25,31,38), each with a `source_of_expected` naming a paper (DOI or PMID + table/figure), and `run_pipeline.py --review` prints "verified" for 5 of 5. Counts invented or copied from our own parse do not count (CLAUDE.md:19-22).
3. **REWRITTEN 19.10.2026 by the author's decision (answer A):** the rank-normalised random-CV number is a transductive artefact: shown on synthetic data, mechanism on the real cohorts unresolved. (Evidence: LABBOOK 19.10. entries "seesaw test STOPPED", "outside audit", "shared-shift test STOPPED". Both pre-registered mechanism tests stopped at their own synthetic gates; no real-data mechanism run.) _Old wording, kept for the record: "The negative logistic gap under rank normalisation (−0.150, results/summary_ranknorm.csv:4-5) is explained in the LABBOOK by a test that was actually run (LABBOOK.md:85-91 names class balance as the suspect)."_
4. **One sentence answers §1** in the LABBOOK: "of a random-CV AUC of X, Y survives a held-out cohort (raw) and Z survives (within-sample normalised)", with every number traceable to a committed `results/` file.

Why this and not an AUC target: the claim is "the drop is cohort-specific structure, not a bug in our pipeline". A control that could shrink the drop by cheating (reading the test set) must be ruled out, and the labels must be checked, before any number means anything. We do not pick an AUC target on purpose.

What it is NOT: not a better classifier, not a new signature, not clinical use (CLAUDE.md:9-12). "An engineering and reproducibility contribution, not a scientific-novelty one. Say it that way." (README.md:80-81).

## 3. Kill criterion
**ACCEPTED 2026-10-18** by the author (same entry). Stop the lab if either happens; both are cheaper than §2:
- **The measurement was a bug.** Published counts for any of the 3 still-unverified cohorts disagree with our parse, and fixing the parse (never the counts, CLAUDE.md:21-22) changes the gap's sign or removes it for both models. Then the 2026-10-10 numbers were a label artefact and there is nothing to report but the code.
- **Scooped.** A public, runnable repository is found that already computes the random-CV number and the held-out-cohort number from the same public data in one command. The README's whole justification is that no such code exists (README.md:77-81); if it does, the reproduction contribution is gone. Check: repeat the 2026-10-10 prior-art sweep (LABBOOK.md:119-127) and find one.

NOT a kill: the gap vanishing under within-sample normalisation. That is an answer to §1 ("the drop is removable with a bedside-legal transform"), and it ends the lab as finished, not failed. Likewise "no cohort effect left" is a finding to write down, not a reason to hide the run.

## 4. Core ideas
- **Report the gap, not the best number** — every result shows random-CV and held-out-cohort AUC side by side (CLAUDE.md:26-28). Without it the lab becomes the failure mode it measures.
- **Labels are guilty until verified** — `config/cohorts.json` carries published counts; a mismatch means our parse is wrong, never the counts (CLAUDE.md:19-22). Label parsing is "the single easiest place to silently poison this kind of study" (README.md:28-29).
- **Exclusions are checked before class rules** — "bacterial and viral co-infection" contains the word "bacterial"; exclusions run first (CLAUDE.md:23-25). Flip the order and co-infections silently become bacterial cases.
- **No feature selection outside the fold** — gene choice happens inside cross-validation or not at all (CLAUDE.md:16-18). Otherwise the experiment is rigged by accident.
- **A held-out cohort is seen one patient at a time** — any transform that reads the whole test cohort is transduction (learning from the test set), not generalisation; "a bedside test seeing one patient at a time cannot do it" (LABBOOK.md:78-81; CLAUDE.md:45-46). Without this rule the "fix" is the leak.
- **Tests go red on broken code** — each test constructs the bug it guards against (CLAUDE.md:29-30).

## 5. Not goals
- Not a new signature, better AUC, clinical product, or bigger dataset (README.md:73-75). The cohort set is frozen at these five unless the author says otherwise (fence accepted 2026-10-18).
- Not modelling co-infections or non-infectious inflammation; they are excluded, and the README says this flatters the numbers (README.md:105-107).
- Not claiming the method is novel — leave-one-dataset-out is packaged in MetaIntegrator (LABBOOK.md:135-137).
- Not extrapolating to RNA-seq; the data are microarray (CLAUDE.md:47).
- Not a prospective clinical failure rate (README.md:110-111).

## 6. Current belief (as of 2026-10-11)
- Data: 5 GEO cohorts, 476 labelled samples, 12,040 shared gene symbols (results/run_raw.json; LABBOOK.md:28-29). Labels verified for GSE63990 and GSE42026 only (results/label_review.md:10,30); GSE60244, GSE40396, GSE6269 are unverified (results/label_review.md:50,70,90).
- Raw scale, ONE run on 2026-10-10 (results/run_raw.json): logistic L2 random-CV 0.954 vs leave-one-cohort-out 0.885, gap 0.069; LightGBM 0.956 vs 0.825, gap 0.131 (results/summary_raw.csv:2-5). LABBOOK.md:59 calls this "real and reproducible", but no second run or seed sweep is on file — treat "reproducible" as UNVERIFIED until a repeat run is committed. **UPDATED 21.10.2026 by the author's decision: "reproducible" is now MEASURED over 10 fold seeds** (seeds 200-209, raw, 463 samples after the GSE40396 rule fix; LABBOOK 2026-10-21 RESULT, `results/seedsweep_raw.csv`): gap logistic 0.071 ± 0.004 (0.065-0.076), LightGBM 0.111 ± 0.003 (0.106-0.117); LightGBM's own model seed moves its held-out AUC 0.816-0.856. _Old wording above kept for the record._
- The flexible model pays more: LightGBM fell to 0.625 on held-out GSE6269 where logistic held 0.944 (results/folds_raw.csv:20,10; LABBOOK.md:64-68).
- Rank-normalised results: UNVERIFIED — do not quote. The apparent gap closure may be transduction leakage, and the negative logistic gap is unexplained (LABBOOK.md:69-91).
- Prior art: the question is not taken (Bodkin 2022 removed the batch effect by construction), but the method is not novel and the scale is about a tenth of Bodkin's (LABBOOK.md:131-143).

## 7. Working rules for anyone extending this lab (human or AI agent)
MAY propose:
- CPU-only runs of `run_pipeline.py` with new normalisation flags, in this order: (1) within-sample normalisation control, (2) explain/eliminate the negative logistic gap, (3) published counts for GSE60244, GSE40396, GSE6269 (literature lookup, no code), (4) per-cohort error analysis (LABBOOK.md:102-109). Also: pull Bodkin Table 2 numbers for comparison (LABBOOK.md:153-155). One proposal = one of these, not several.
- Budget: € 0 and 0 GPU-hours per proposal; nothing here needs the GPU, "if something seems to, question it first" (CLAUDE.md:34). Cap: no single run over 1 CPU-hour without a LABBOOK plan entry first.
- Network use only to fetch the five listed GEO series (`--fetch`, already cached in `data/`) and to read published papers for counts.
- Work on a feature branch; never commit to the default branch without the author.
- Every session appends a dated LABBOOK entry (CLAUDE.md:38-40); every result lists both schemes.

NEVER:
- Change §1–§5, or claim more than CLAUDE.md:5-7 allows.
- Quote rank-normalised numbers outside the repo (LABBOOK.md:90-91) or cite any §6 line marked UNVERIFIED.
- Edit `expected_counts` to match the code (CLAUDE.md:21-22), or fill them from anything but a named published source.
- Add feature selection outside the fold, or a transform that reads other test-cohort patients, and call it generalisation.
- Add or swap cohorts, add RNA-seq data, or model co-infections as a side quest.
- Buy anything, serve or ship a trained model, or present it as diagnostic.
- Report a number without its two-scheme partner (random-CV and leave-one-cohort-out together, CLAUDE.md rule 4).

Status (2026-10-21): §2 and §3 were accepted on 2026-10-18. All four §2 conditions are met per the LABBOOK (condition 2: `results/label_review.md` prints "verified" for 5 of 5 cohorts; three of the five count checks were done by an AI assistant and not yet by a human, see README). Still open: a human read of the label table and the code.
