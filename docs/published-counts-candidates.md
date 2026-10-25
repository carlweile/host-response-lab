# Published class counts: CANDIDATES (19.10.2026)

_Historical record. These candidates have since been checked and the counts are now in `config/cohorts.json` with their sources (see README, Labels)._

> **CANDIDATES ONLY, not in config/cohorts.json.** The author checks each number
> against the paper, then fills `expected_counts` + `source_of_expected`
> (CORE §2.2).

Two blind agents read each paper on their own, without seeing each other. Below,
their results sit side by side. "Agree" means both found the same numbers.
The automatic flag only compares the raw lists; the "my reading" line says
whether the difference is real or only how the list was written down.

---

## GSE60244 (adults, lower respiratory tract infection) — AGREE

| | Extraction 1 | Extraction 2 |
|---|---|---|
| Paper | Suarez et al., J Infect Dis 2015;212:213-22, PMID 25637350, doi:10.1093/infdis/jiv047 | same |
| Where | Abstract; Results "Patient Demographic Characteristics and Etiologic Diagnosis"; Methods "Microarray Data and Statistical Analysis" | same, plus GEO "condition" field |
| Viral | 71 | 71 |
| Bacterial | 22 | 22 |
| Bacterial + viral co-infection | 25 | 25 |
| Healthy controls | 40 | 40 |
| Total | 158 | 158 |
| Unit | patients = samples (one blood array per person) | same |
| Confidence | high | high |

Quote (both): "158 samples [22 bacterial LRTI, 71 viral LRTI, 25 bacterial-viral
coinfections, and 40 healthy controls]".

Notes: 6 of the 71 viral are virus+virus; 3 of the 22 bacterial are
bacteria+bacteria (they stay in their groups). An erratum exists (PMID 26598315)
that neither agent could read; a possible clue is that GEO says 97 %
sensitivity and the paper says 95 %. This is not a class count.

**My reading:** agree on every number.

---

## GSE40396 (children, fever) — flag says DISAGREE, numbers really AGREE

| | Extraction 1 | Extraction 2 |
|---|---|---|
| Paper | Hu et al., PNAS 2013;110:12792-7, PMID 23858444, doi:10.1073/pnas.1302968110 | same |
| Where | Results first paragraph; Methods "Microarray Data Analysis"; Abstract | Results first subsection; Methods > Subjects; GEO fields |
| Febrile viral total | 22 (8 HHV-6 + 8 adenovirus + 6 enterovirus; agent's own sum) | 22 (same split) |
| Febrile bacterial | 8 | 8 |
| Afebrile virus-positive | 2 HHV-6 + 3 adeno + 8 rhino (listed apart) | 13 (same split) |
| Afebrile virus-negative controls | 22 | 22 |
| Co-infection | 0 (excluded by design) | 0 |
| Total | 65 | 65 |
| Unit | children = arrays (70 run, 5 outliers removed) | same |
| Confidence | high | high |

Quote (both): "...30 febrile children [8 ... HHV-6, 8 ... adenovirus, 6 ...
enterovirus, and 8 children with acute bacterial infection]. ... 35 afebrile
children (2 ... HHV-6, 3 ... adenovirus, 8 ... rhinovirus, and 22 virus-negative
control children)."

Notes: the flag is "disagree" only because agent 1 listed the 13 afebrile
virus-positive children as three lines and agent 2 as one line. Same numbers.
Watch out: GEO summary says "51 patients" and the abstract implies 52; both are
wrong against the 65 in Results and in GEO. The 13 afebrile virus-positive
children are neither "viral infection" nor "control"; the paper's own
viral-vs-bacterial test uses only the 30 febrile children (22 vs 8).

**My reading:** agree. Decision for the author: how to label the 13.

---

## GSE6269 (mostly children, PBMC, three platforms) — DISAGREE on the total only

| | Extraction 1 | Extraction 2 |
|---|---|---|
| Paper | Ramilo et al., Blood 2007;109:2066-77, PMID 17105821, doi:10.1182/blood-2006-02-002477 | same |
| Where | Methods "Patient information"; Tables 1-3; Fig. 6 legend | Methods; Tables 1-3; Microarray assay; Discussion |
| Influenza (patients) | 36 (really 30 A + 6 B) | 36 (same remark) |
| E. coli | 29 | 29 |
| S. aureus | 50 | 50 |
| S. pneumoniae | 16 patients (22 arrays) | 16 patients |
| Bacterial total | 95 patients | 95 patients / 101 arrays |
| Healthy controls | 7 | 7 in paper, 6 in GEO |
| Arrays per platform | U133A 91 infected, U133plus2 22, Illumina Hu6 24 | same |
| Total | **144** (arrays in the paper) | **143** (arrays in GEO) |
| Unit | mixed: 131 patients, 137 infected arrays | same |
| Confidence | medium | high |

Quote (both): "Altogether 144 microarray analyses ... including 137 on samples
collected from 131 patients" and GEO: "Only 143 are included in this series
because the CHP/CEL file was not available for sample PBMC_Healthy_INF295."

Notes: 144 vs 143 is not a real conflict: 144 = the paper, 143 = GEO (one
healthy control missing). Both agents report both facts. Other flags: 6
S. pneumoniae patients (96, 265, 268, 277, 287, 289) are on two platforms, so
the same person appears twice (matters for a held-out test). The paper says
"influenza A" but 6 are influenza B. GEO sample INF057 says "S. pneumoniae" in
its name but "E. coli" as pathogen (paper says E. coli). The paper says 98
bacterial patients in one place, 95 elsewhere.

**My reading:** agree on every class count; choose whether
`expected_counts` uses patients (paper) or arrays in GEO (143).

---

## How to check (about 10 minutes)

1. **GSE60244** — open https://doi.org/10.1093/infdis/jiv047 . Read the
   abstract Results sentence and the Methods line "158 samples [...]". Compare:
   22 / 71 / 25 / 40, total 158.
2. **GSE40396** — open https://doi.org/10.1073/pnas.1302968110 . Read the
   first paragraph of Results. Compare: febrile 8 HHV-6, 8 adeno, 6 entero,
   8 bacterial; afebrile 2 / 3 / 8 virus-positive, 22 controls; total 65.
   Then decide how to label the 13 afebrile virus-positive children.
3. **GSE6269** — open https://doi.org/10.1182/blood-2006-02-002477 . Read
   "Patient information" (first sentence) and the footnotes of Tables 1 and 2.
   Compare: 29 E. coli, 50 S. aureus, 16 S. pneumoniae, 36 influenza, 7
   healthy. Then open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6269
   and check the "Overall design" line says 143.
4. Write the numbers you confirmed into `config/cohorts.json`
   (`expected_counts` + `source_of_expected`).

---
---

## Comparison with our own parse (results/label_review.md, run of 2026-10-10)

_Added 19.10.2026 (night run), after the two blind extractions above. Nothing
was re-run and nothing was changed: not the parser, not `config/cohorts.json`.
Our numbers are copied from `results/label_review.md`; the per-group splits in
the right-hand column were re-counted by hand from the GEO series files in
`data/series/` (read only)._

| Cohort | Our parse (bacterial / viral / excluded) | Blind extractions (paper) | Match? |
|---|---|---|---|
| GSE60244 | 22 / 71 / 65 | 22 bacterial, 71 viral; 25 co-infections + 40 healthy = 65 excluded | **yes** |
| GSE6269, GPL96 only (`prefer_platform`) | 73 / 18 / 6 | U133A (= GPL96): 91 infected arrays. GEO titles on GPL96: 18 influenza A; 29 E. coli + 31 S. aureus (19 MRSA + 12 MSSA) + 13 S. pneumoniae = 73; 6 healthy | **yes** |
| GSE40396 | 8 / 35 / 22 | febrile: 8 bacterial, 22 viral; afebrile: 13 virus-positive, 22 virus-negative controls | **NO** |

**GSE60244.** GEO "condition" field: BACTERIA 22, VIRUS 71, COINFECTION 25,
Healthy Control 40. Our rule puts the 25 co-infections and 40 controls in
"excluded", as the paper's own bacterial-vs-viral comparison does.

**GSE6269.** On GPL96 our counts equal the GEO sample titles and add up to the
paper's 91 U133A infected arrays (73 + 18). Two small notes that do not change
any class: GEO has 6 healthy GPL96 arrays (the paper says 7 healthy overall;
one CEL file is missing in GEO), and sample INF057 is named S. pneumoniae but
listed as E. coli; both are bacterial, so our label is the same either way.
Not re-checked by me: the per-pathogen split on U133A inside the paper's own
tables (only the 91 total is in the extractions above).

**GSE40396: our parse does not match the paper's design.** Our 35 "viral" =
22 FEBRILE virus-positive children (8 HHV-6, 8 adenovirus, 6 enterovirus) + 13
AFEBRILE virus-positive children (2 HHV-6, 3 adenovirus, 8 rhinovirus; "fever:
no" in GEO). Our rule matches the virus name and ignores the fever field. The
paper does not count the 13 afebrile children as viral infections; its
viral-vs-bacterial classifier used only the 30 febrile children (22 viral vs 8
bacterial). Our 8 bacterial and 22 excluded (virus-negative controls) match.

**What this means.**
- By `CLAUDE.md` rule 2 (labels are guilty until verified) our GSE40396 rule is
  probably wrong: 13 of its 35 viral labels are children without an infection
  episode in the paper's sense.
- The fix is to the RULE (for example: in GSE40396, viral only with
  "fever: yes"), never to the counts. Which rule to use is the author's call after
  they read the paper (step 2 of "How to check" above). Nothing was changed.
- CORE §3 kill criterion: if fixing the parse (never the counts) changes the
  gap's sign or removes it for both models, "the 2026-10-10 numbers were a
  label artefact and there is nothing to report but the code", so the lab
  stops. The re-run after a rule fix is a real test, not a formality.
- Every number on file (including tonight's gap diagnostics) inherits the
  current GSE40396 labels until then.
