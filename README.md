# GSE278446 NeuroNovelty

## Status

This package is the **authoritative corrected analysis pipeline** for the GSE278446 NeuroNovelty study.

**Pipeline version:** `GSE278446_NeuroNovelty_FINAL_1.0`  
**Analysis output version:** `GSE278446_v4_final`

The previous `GSE278446_v3` outputs are treated as **legacy/intermediate results** and must not be mixed with the final run.

---

## 1. Critical source-level finding

The original GEO source contains **28 GSM records**.

An exact source-level duplicate was independently confirmed:

- Retained source representative: `GSM8546664`
- Excluded duplicate: `GSM8546673`

The duplicate was verified across:

1. matrix content
2. barcode content
3. feature content

using decompressed-content SHA-256 fingerprints.

### Analytical rule

| Layer | GSM count |
|---|---:|
| Original GEO source inventory | 28 |
| Excluded source duplicate | 1 |
| Independent biological analysis | 27 |

The raw source files are **not deleted or modified**.

`GSM8546673` is excluded from downstream biological inference.

This distinction must be stated explicitly in the manuscript.

---

## 2. Common temporal window

The definitive reversal analysis uses the shared timepoints:

`0, 24, 72 hpi`

After exclusion of the source duplicate, the complete independent trajectories are:

- Exp1_R2
- Exp2_R1
- Exp2_R2

Therefore the strict validation criterion is:

**3/3 available independent trajectories**, not 4/4.

A gene is a strict reversal candidate only when the same 24-h peak/dip pattern is present in all three trajectories and both early and late relative-change magnitudes pass the predefined 20% threshold.

---

## 3. Biological unit of inference

Single cells are used for:

- QC
- descriptive cell-level program scoring
- localization/exploratory analyses

Cells are **not treated as independent biological replicates** for temporal inference.

The biological sample/GSM remains the primary experimental unit.

---

## 4. Pipeline stages

| Stage | Purpose |
|---|---|
| 00 | Raw-file smoke test |
| 01 | Representative matrix load/dimension test |
| 02 | Representative QC + marker availability |
| 03 | Complete source dataset inventory |
| 04 | Sample identity/design audit |
| 05A | Source-level duplicate verification |
| 05 | Full raw matrix SHA-256/content audit |
| 06 | Single-cell preprocessing + sample gene means |
| 07 | Program temporal discovery |
| 08 | Gene temporal discovery |
| 09 | Program replicate validation |
| 10 | Common 0/24/72-hpi reversal discovery |
| 11 | Strict 3/3 independent trajectory validation |
| 12 | Temporal sampling sensitivity |
| 13 | Marker/program overlap audit |
| 14 | Cell-level reversal localization |
| 15 | Gene characterization |
| 16 | Gene prioritization |
| 17 | Enrichment input preparation |
| 18 | Functional enrichment |
| 19 | Pathway reproducibility interpretation |
| 20 | NCBI gene-symbol/alias audit |
| 21 | Broad literature context audit |
| 22 | Independent literature cross-check |
| 23 | Targeted literature verification |
| 24 | Paper-level verification |
| 25 | Independent literature recheck |
| 26 | Candidate biological-context classification |
| 27 | Candidate evidence audit |
| 28 | Ranking sensitivity/robustness |
| 29 | Final candidate synthesis |
| 30 | Final paper-level recheck |
| 31 | Final reproducibility/provenance audit |

---

## 5. How to run

Place this folder at:

```text
Research/
└── Final Tests/
    └── Final Pipeline/
        ├── RUN_ALL.py
        ├── config.py
        └── ...
```

and keep the raw source at:

```text
Research/
└── GSE278446_RAW/
```

Install dependencies:

```text
pip install -r requirements.txt
```

Then run:

```text
python RUN_ALL.py
```

The runner automatically infers the Research root from the folder structure. An environment variable can override it:

```text
GSE278446_RESEARCH_ROOT
```

The final outputs are written to:

```text
Research/
└── NeuroNovelty/
    └── temporal_results/
        └── GSE278446_v4_final/
```

---

## 6. Important reproducibility rule

Do **not** mix outputs from:

- `GSE278446_v3`
- earlier experimental pipelines
- older enrichment runs
- older literature audits

with the final `GSE278446_v4_final` output.

The final run creates a fresh versioned result directory specifically to prevent accidental mixing.

---

## 7. Enrichment

Functional enrichment uses **Danio rerio (`drerio`)**.

The pipeline uses the g:Profiler API and records the API/data-version metadata used during the run.

The tested gene universe is used as the enrichment background rather than the entire genome.

If enrichment returns no valid terms, the pipeline records that result as unavailable. It does **not** fabricate pathway interpretation.

---

## 8. Literature and novelty rule

The literature stages are evidence audits, not automatic novelty detectors.

The following are deliberately prohibited:

- treating zero PubMed hits as proof of novelty
- calling a gene novel because no direct paper was retrieved
- inferring cell type from arbitrary pathway names
- converting pathway-list evidence into a gene-specific score without explicit mapping
- treating single cells as independent biological replicates
- silently excluding the source duplicate

Final candidates are therefore presented as **review/prioritization candidates**, not automatically as novel discoveries.

---

## 9. Final audit gate

Stage 31 must finish with:

```text
PASS — pipeline output, provenance, sample-independence, and code audit complete.
```

It also creates:

```text
pipeline_audit/
├── PIPELINE_OUTPUT_FILE_AUDIT.csv
├── PIPELINE_LOGICAL_CHECKS.csv
├── PIPELINE_CODE_HASHES.csv
└── PIPELINE_MANIFEST.json
```

The manifest records the final pipeline version, analysis version, sample/exclusion information and SHA-256 hashes of the authoritative Python files.

---

## 10. Manuscript handoff

Do not write the manuscript from screenshots or from memory.

Use only:

```text
GSE278446_v4_final/
```

after Stage 31 passes.

Recommended manuscript evidence chain:

```text
Raw GEO source
      ↓
Source inventory
      ↓
Source duplicate verification
      ↓
Raw SHA-256 audit
      ↓
QC / preprocessing
      ↓
Sample-level temporal analysis
      ↓
3/3 independent trajectory validation
      ↓
Sensitivity / robustness
      ↓
Literature evidence
      ↓
Final candidate review
```

The final manuscript must preserve both positive and negative results, including unavailable enrichment or literature evidence when applicable.
