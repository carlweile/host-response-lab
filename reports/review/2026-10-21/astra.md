# External review — 2026-10-21

Model: `gpt-6-astra` (Responses API, reasoning high) · 2026-10-21 03:22 · status `completed` · in 41937 / out 14956 tokens · 5.6 min

Prompt: `reports/review/2026-10-21/prompt.md` (sha256 `68f54f8b777101df`)

Inputs, exactly as sent:

| file | role | chars | sha256 (first 16) |
|---|---|---|---|
| `reports/review/2026-10-21/bundle.txt` | code_bundle | 144,516 | `d27e828c216be272` |

---
1. **Severity: HIGH — Test patients select the probe used as a feature**
   - **Where:** `src/harmonize.py::to_symbols`; called before splitting by `run_pipeline.py::cmd_train` and `run_gap_diagnostics.py::load_matrix`.
   - **What:** The highest-mean probe for each gene is selected using every patient in the series, including test patients and subsequently excluded samples. This is data-dependent feature selection outside CV and makes even the within-sample setting depend on other patients.
   - **Failure scenario:** Two probes map to gene A: training values are `[10, 10]` and `[9, 9]`, while the test values are `0` and `100`. The second probe wins; changing only its test value to `0` selects the first probe and changes the training features.
   - **Suggested test:** Hold training expression and annotations fixed, perturb only held-out expression, and assert that the selected probes and transformed training features remain unchanged. The current implementation fails.

2. **Severity: HIGH — Imputation is transductive in both validation schemes**
   - **Where:** `src/harmonize.py::build_matrix`.
   - **What:** Missing values are filled from whole-cohort means before any split. Random-CV training rows can therefore depend on test rows, while LOCO test patients depend on other patients in the held-out cohort; the subsequent within-sample ranking does not undo this dependence.
   - **Failure scenario:** A gene has values `[NaN, 2, 4]` in a held-out cohort. The first patient's value becomes `3`; changing only the third patient to `100` changes it to `51`, potentially changing that first patient's gene ranks and prediction.
   - **Suggested test:** Exercise `build_matrix` followed by within-sample normalization on a fixture containing missing values. Change another held-out patient and assert that the target patient's representation is unchanged; also assert that changing random-CV test values cannot alter imputed training values.

3. **Severity: HIGH — Published-count mismatches do not block training**
   - **Where:** `src/labels.py::_verify`; `run_pipeline.py::cmd_review`, `cmd_train`; `run_gap_diagnostics.py::load_matrix`.
   - **What:** A count mismatch only produces `verified=False`; every training path still uses those labels. `cmd_train` even reduces the distinction between a known mismatch and absent verification to “labels UNVERIFIED,” violating the rule that a mismatch means the parse is wrong.
   - **Failure scenario:** A cohort parses as 80 bacterial and 110 viral instead of the configured 73 and 117. Provided enough samples remain, the pipeline fits models and writes ordinary headline results.
   - **Suggested test:** Supply a mismatched `LabelResult` with at least 20 usable samples and otherwise valid expression/annotation. Assert that training returns failure before calling any estimator or replacing result files; test the diagnostic loader too.

4. **Severity: MEDIUM — Held-out expression controls gene eligibility**
   - **Where:** `src/parse_geo.py::parse_series_matrix`; `src/harmonize.py::build_matrix`.
   - **What:** The shared feature axis is not determined solely by fixed platform annotations: the parser first drops all-NaN probes, and the resulting gene indexes are intersected across cohorts. Thus held-out expression missingness can remove a gene from the training feature set.
   - **Failure scenario:** Gene G is annotated on every platform and measured in all training cohorts. Making its last finite measurement in the held-out cohort missing drops G from the entire experiment, including training.
   - **Suggested test:** Keep annotations and training data fixed, change only held-out missingness, and assert that training feature eligibility remains fixed. This should distinguish a legitimate predefined platform intersection from expression-dependent filtering.

5. **Severity: MEDIUM — Log-scale detection reads the entire series**
   - **Where:** `src/harmonize.py::maybe_log2`; its pre-CV callers.
   - **What:** A maximum over all samples determines whether every value in the series is log-transformed. One held-out or excluded sample can therefore change the preprocessing of training patients; the threshold also assumes that every linear-scale series exceeds 100.
   - **Failure scenario:** All training values are below 100. Changing a single held-out maximum from `100` to `101` changes training preprocessing from identity to `log2(x + 1)`, potentially changing the raw-model AUC.
   - **Suggested test:** Perturb only a held-out value across the threshold and assert unchanged training preprocessing. Separately test a known linear-scale fixture whose maximum is below 100, using an explicit scale declaration rather than the heuristic as the oracle.

6. **Severity: MEDIUM — GPL96 is a preference, not an enforced restriction**
   - **Where:** `src/fetch_geo.py::fetch_series`; `run_pipeline.py::_load_all`; `src/parse_geo.py::parse_series_matrix`.
   - **What:** An existing cache or successful unqualified download is accepted without checking `prefer_platform`; if the preferred variant is absent, another platform is silently selected. Variant matching uses substring matching, and the loader never validates the parsed platform against the configuration.
   - **Failure scenario:** A previous fetch cached a non-GPL96 GSE6269 matrix. Adding or retaining `"prefer_platform": "GPL96"` does not replace or reject it, so training uses the wrong patient/platform subset.
   - **Suggested test:** Seed the cache with a matrix declaring another platform and request GPL96. Require rejection or a verified replacement; also test missing preferred variants and exact matching against similarly prefixed platform IDs.

7. **Severity: MEDIUM — Label rules interpret negated tests as diagnoses**
   - **Where:** `src/labels.py::_classify`, `assign_labels`; `config/cohorts.json::label_rules`.
   - **What:** Unrestricted substring searches over concatenated fields do not distinguish diagnoses from negative test results, field names, or negated exclusions. Correct exclusion ordering does not solve this.
   - **Failure scenario:** `"diagnosis: viral infection | bacterial culture: negative"` becomes bacterial. `"diagnosis: viral infection | co-infection: no"` is excluded despite explicitly denying co-infection.
   - **Suggested test:** Add these cases using the actual configuration, requiring viral labels. Include representative negative pathogen tests and negated exclusion fields; aggregate count agreement alone cannot detect compensating misclassifications.

8. **Severity: MEDIUM — The production training path does not check row alignment**
   - **Where:** `src/train.py::run_both_schemes`, `_fit_score`.
   - **What:** Training and splitting use positional indexing without checking that `X`, `y`, and `cohort` have identical indexes in identical order. The alignment assertion tested in `both_schemes_with_scores` does not protect the main pipeline or seed sweep.
   - **Failure scenario:** A cached tuple contains `y` or `cohort` in a different row order from `X`. The code assigns labels or cohort memberships to the wrong patients; the cohort-overlap assertion can still pass because it checks the same misordered cohort vector.
   - **Suggested test:** Call `run_both_schemes` with `y.iloc[::-1]`, then separately with a permuted cohort series, and require an alignment error before fitting. Both currently proceed.

9. **Severity: MEDIUM — Partial studies still produce ordinary headline files**
   - **Where:** `run_pipeline.py::cmd_train`; `src/train.py::make_models`, `summarize`.
   - **What:** Training accepts any two surviving cohorts, silently omits LightGBM on `ImportError`, and computes the gap after dropping undefined AUCs. Although notes and fold counts disclose some omissions, there is no completeness gate for the prescribed five-cohort, two-model comparison.
   - **Failure scenario:** One LOCO test cohort has only one class. Its AUC is removed, and the headline becomes five random-CV folds minus four LOCO cohorts; alternatively, a missing annotation removes a cohort and overwrites the same `summary_raw.csv`.
   - **Suggested test:** Supply one undefined LOCO fold and require the headline gap to be invalid or explicitly marked incomplete. Integration-test missing configured cohorts and a failed LightGBM import, requiring refusal to publish a complete-study result.

10. **Severity: MEDIUM — Non-converged logistic fits are reported without warning**
    - **Where:** `src/train.py` module-level warning filter; `_fit_score`.
    - **What:** Suppressing all `UserWarning`s also suppresses scikit-learn's `ConvergenceWarning`. The main scorer neither records `n_iter_` nor marks capped fits, so an unfinished optimization produces a normal AUC and empty note.
    - **Failure scenario:** A high-dimensional fold reaches `max_iter=3000`; its approximate solution changes predictions, but the main results contain no indication that the intended logistic fit did not converge.
    - **Suggested test:** Force a demonstrably non-converged fit with `max_iter=1`. Require a surfaced warning, failure, or explicit non-convergence status in the fold result rather than an unqualified score.

11. **Severity: MEDIUM — Reproduction gates fail open on missing or NaN comparisons**
    - **Where:** `run_seed_sweep.py::repro_gate`, `main`; `run_seesaw_test.py::repro_gate`; `run_shared_shift_test.py::repro_gate`; `run_gap_diagnostics.py::_with_file_auc`, `evaluate`.
    - **What:** These checks do not consistently require finite comparisons and complete reference coverage. `NaN > tolerance` is false, `max(0.0, NaN)` can retain zero, and pandas maxima skip missing comparisons.
    - **Failure scenario:** Reference AUC columns are NaN while fold keys still match. The seesaw/shared-shift gates can return `ok=True` with zero worst difference, and the seed-sweep gate can proceed with a NaN difference; E0 can ignore unmatched reference rows.
    - **Suggested test:** Replace reference AUCs with NaN and remove individual reference rows. Require every gate to reject missing/nonfinite comparisons and verify the exact expected model/scheme/fold coverage.

12. **Severity: MEDIUM — The timing probe overwrites real diagnostic results**
    - **Where:** `run_seesaw_test.py::main`.
    - **What:** `todo = [p for p in real if p.startswith("p")]` includes `"probe"`. After timing finishes, a second `run_parts` call receives `{"probe"}`, runs no experimental jobs, and writes empty frames over the primary folds and offsets files.
    - **Failure scenario:** Run the documented `--parts probe` command after a completed experiment. Existing CSV results are erased, while the loaded metadata can retain an old verdict alongside newly updated run metadata.
    - **Suggested test:** Create nonempty result files, run `main(["--parts", "probe"])` with timing/model work stubbed, and assert that the primary results and judged-result metadata are unchanged. Restrict experimental parts by explicit membership, not a prefix.

13. **Severity: LOW — The main scoring smoke test accepts a completely broken scorer**
    - **Where:** `tests/test_pipeline_smoke.py::test_both_schemes_run_and_random_cv_looks_better`.
    - **What:** The guard requires only valid-range scores and `random_cv >= LOCO`, so equality at chance passes. Replacing the production scorer with a constant `0.5` therefore satisfies the test advertised as protecting the measurement machinery.
    - **Failure scenario:** `_fit_score` always returns `(0.5, "")` without fitting anything. The smoke test still passes, as do the basic fold-size and LOCO seed-invariance checks.
    - **Suggested test:** Add a mutation test requiring the existing guard to reject that constant-score implementation—it currently does not. Strengthen the production-path test with independently calculated per-fold AUCs or a deterministic fixture that must score above chance.

HIGH-FOUND
