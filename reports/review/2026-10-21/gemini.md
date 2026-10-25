# External review - Gemini

Model: gemini-3.1-pro-preview | 2026-10-21 03:16:49 | prompt sha256 e1c2a1c862c28eb8 (= prompt.md + bundle.txt, see manifest.md)

---

Here are the findings from the review, ordered by severity.

### 1. Platform preference is completely ignored for multi-platform series
- **Severity**: HIGH
- **Where**: `src/fetch_geo.py` in `fetch_series`
- **What**: The code attempts to download the combined series matrix (`series_matrix_url(accession)`) first and immediately returns it if successful, completely bypassing the `prefer_platform` argument.
- **Failure scenario**: For GSE6269, GEO provides a combined matrix containing all samples across its three platforms. The pipeline successfully downloads this combined matrix instead of the specified `GPL96` variant. As a result, the same patients are included multiple times (once for each platform) in the single cohort, artificially inflating the sample size and causing exact patient duplicates across the train/test splits.
- **Suggested test**: Mock `download` to successfully return a dummy combined matrix, provide `prefer_platform="GPL96"`, and assert that the function skips the combined matrix and fetches the platform-specific one instead.

### 2. Viral substring rule erroneously catches bacterial pathogens
- **Severity**: HIGH
- **Where**: `config/cohorts.json` in `label_rules`
- **What**: The substring `"influenza"` is used as a rule to classify samples as viral. However, this incorrectly catches *Haemophilus influenzae*, which is a very common respiratory bacterium. 
- **Failure scenario**: A sample with the characteristic `"pathogen: haemophilus influenzae"` will not match any of the bacterial rules. It will then fall through to the viral rules, match `"influenza"`, and be incorrectly classified as viral. 
- **Suggested test**: Add a test in `test_honest_splits.py` passing `"infection: haemophilus influenzae"` to `assign_labels` and assert that it maps to `"bacterial"`. It will currently fail and return `"viral"`.

### 3. The published-count label verification is silently bypassed
- **Severity**: HIGH
- **Where**: `run_pipeline.py` in `cmd_train`
- **What**: When `assign_labels` detects a mismatch against the published expected counts, it flags the cohort as `verified = False` and attaches a `"MISMATCH"` verdict. However, `cmd_train` merely logs this as `"labels UNVERIFIED"` in the notes and proceeds to use the broken labels in the dataset anyway.
- **Failure scenario**: If a label parsing rule is broken or misconfigured (causing parsed counts to heavily disagree with published counts), the pipeline will log a warning but still execute the whole experiment, directly violating the rule that "a mismatch means the parse is wrong... tests must go red". 
- **Suggested test**: Write a test that feeds `cmd_train` a dummy cohort configuration with mismatched expected counts and assert that it drops the cohort or raises an exception. 

### 4. Leakage of test-set statistics during Harmonisation
- **Severity**: HIGH
- **Where**: `src/harmonize.py` in `to_symbols` and `build_matrix`
- **What**: The dataset applies feature extraction and imputation on the *entire dataset* before the cross-validation loops begin. 
- **Failure scenario**: 
  1. `to_symbols` collapses multiple probes into one gene by picking the probe with the highest mean. By calculating this mean over the whole cohort, test-set expression levels influence the feature selection for the training set.
  2. `X.fillna(col.mean())` imputes missing gene values using the mean of the entire cohort, bleeding test-set expression values directly into the training features.
- **Suggested test**: Extract the selected probes and imputed values for a training split, change the expression values of the test-split rows, re-run harmonisation, and assert that the training-split features remain mathematically identical.

HIGH-FOUND
