# Final Pipeline Audit Summary

## Package reviewed

Input package: `Final Tests (2).zip`

The package contained 33 executable pipeline stages plus configuration/runner files, README, requirements and old execution artifacts.

## Critical defects found and corrected

### 1. Stage 05A → Stage 05 contract failure
**Problem:** Stage 05 expected `05A_source_duplicate_verification_summary.csv`, but Stage 05A only produced a text report.

**Correction:** Stage 05A now writes:
- `05A_original_geo_source_verification.txt`
- `05A_source_duplicate_verification_summary.csv`
- `05A_source_duplicate_file_comparison.csv`

Stage 05 validates this machine-readable record before proceeding.

### 2. Duplicate detection based only on compressed bytes
**Problem:** Stage 05 used only compressed-file SHA-256 for duplicate detection.

**Correction:** Stage 05 now records:
- compressed SHA-256
- decompressed content SHA-256

Duplicate classification uses decompressed content SHA-256.

### 3. Cell-level analysis reintroduced the excluded duplicate
**Problem:** Stage 14 iterated over all `SAMPLE_METADATA`, which included GSM8546673.

**Correction:** Stage 14 now iterates only over `ANALYSIS_GSMS` and explicitly fails if the excluded GSM enters the cell-level output.

### 4. Stage 28 boolean parsing was unsafe
**Problem:** `astype(bool)` can interpret the string `"False"` as `True`.

**Correction:** Explicit string-to-boolean parsing is now used.

### 5. Stage 26/27/29 schema mismatch
**Problem:** Stage 29 expected fields that Stage 26/27 did not produce:
- `biological_context_class`
- `strict_4of4_evidence`

The final stage was therefore structurally incapable of completing reliably.

**Correction:** The final chain consistently uses:
- `biological_context`
- `strict_3of3_available_evidence`

Stage 26 now evaluates the TOP-20 candidates rather than a single candidate.

### 6. Stale-result mixing
**Problem:** Re-running the pipeline reused `GSE278446_v3`, allowing old outputs to coexist with new outputs.

**Correction:** The authoritative pipeline now writes to:
`GSE278446_v4_final`

### 7. Hard-coded Research path in RUN_ALL
**Problem:** RUN_ALL contained a user-specific absolute path.

**Correction:** The runner infers the Research root from the package location and supports `GSE278446_RESEARCH_ROOT` as an override.

### 8. Stage 31 audit was incomplete
**Problem:** The old final audit checked several files but did not comprehensively verify the duplicate handoff, cell-level exclusion, final schema, or pipeline manifest.

**Correction:** Stage 31 now checks:
- source count
- analysis count
- duplicate decision
- inventory
- identity audit
- SHA-256 duplicate status
- absence of the excluded GSM from gene/cell analysis
- strict 3/3 validation
- final candidate schema
- novelty flag
- required outputs
- code hashes
- pipeline manifest

## Scientific safeguards preserved

The package deliberately does NOT:
- call zero literature hits proof of novelty
- fabricate pathway interpretation
- treat cells as independent biological replicates
- silently remove the duplicate from the raw source
- call the three available trajectories four independent trajectories
- convert pathway-list evidence into a gene-specific score without explicit mapping

## Structural testing performed

The corrected pipeline was syntax-compiled successfully.

The following key stages were also exercised on a synthetic structural dataset:
- Stage 00
- Stage 03
- Stage 04
- Stage 05A
- Stage 05
- Stage 26
- Stage 27
- Stage 28
- Stage 29

The Stage-05 duplicate handoff successfully passed when the synthetic data contained exactly the expected duplicate pair and correctly failed when unrelated duplicate matrices were introduced.

## Important limitation

This uploaded ZIP did not contain the full 27-sample/28-source raw dataset or the complete final result directory. Therefore this package is a **corrected authoritative pipeline**, not a fabricated claim that the entire biological analysis has already been rerun successfully on the real dataset.

The real final result directory becomes authoritative only after the user runs this package against the actual `GSE278446_RAW` source and Stage 31 passes.
