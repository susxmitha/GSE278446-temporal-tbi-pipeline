# Manuscript Handoff — GSE278446 NeuroNovelty

## What is safe to state after the final run

The computational pipeline is designed around a source-level audit of the GEO supplementary data.

The source inventory contains 28 GSM records. One record, GSM8546673, was independently verified as an exact source-level duplicate of GSM8546664 across matrix, barcode and feature content. The duplicate remains in the raw source for provenance but is excluded from downstream biological inference, leaving 27 analysis samples.

The shared 0/24/72 hpi temporal window contains three complete independent trajectories after duplicate exclusion: Exp1_R2, Exp2_R1 and Exp2_R2.

The strict reversal analysis requires the same 24-h peak/dip pattern and at least 20% relative change in both temporal phases for all three available independent trajectories.

## What must not be claimed automatically

Do not write:

- "GSM8546673 was corrupted" — the pipeline establishes duplication, not corruption.
- "The gene is novel" solely because literature searches returned zero hits.
- "The pathway proves the candidate belongs to a particular cell type" without explicit mapping evidence.
- "n cells were independent biological replicates."
- "4/4 replicate trajectories validated the reversal."

## Results to extract for the manuscript

After Stage 31 passes, use the v4_final result directory as the single source of truth.

### Dataset / provenance
- `03_dataset_inventory.csv`
- `04_sample_identity_audit.csv`
- `05A_original_geo_source_verification.txt`
- `05A_source_duplicate_verification_summary.csv`
- `05_full_sha256_matrix_verification.csv`

### QC / preprocessing
- `GSE278446_cell_qc_summary.csv`
- `GSE278446_sample_program_summary.csv`
- `GSE278446_gene_level_sample_data.csv`

### Temporal analysis
- `PROGRAM_TEMPORAL_TRAJECTORIES.csv`
- `COMMON_WINDOW_REVERSAL_ALL_GENES.csv`
- `STRONG_REPRODUCIBLE_REVERSALS.csv`
- `ALL_AVAILABLE_TRAJECTORY_VALIDATION.csv`
- `STRICT_3_OF_3_AVAILABLE_REVERSALS.csv`

### Robustness
- `PHASE_SHIFT_SENSITIVITY_SUMMARY.csv`
- `SENSITIVITY_TOP20_STABILITY.csv`
- `SENSITIVITY_GENE_RANK_STABILITY.csv`

### Functional interpretation
- `ENRICHMENT_RUN_SUMMARY.csv`
- `MASTER_PATHWAY_INTERPRETATION.csv`
- `PATHWAY_REPRODUCIBILITY_SUMMARY.csv`

### Literature
- `BROAD_LITERATURE_CONTEXT_AUDIT.csv`
- `TARGETED_LITERATURE_VERIFICATION.csv`
- `TOP20_PAPER_LEVEL_SUMMARY.csv`
- `LITERATURE_RECHECK_CONCORDANCE.csv`
- `FINAL_5_PAPER_RECHECK.csv`

### Final candidates
- `FINAL_CANDIDATE_EVIDENCE_AUDIT.csv`
- `FINAL_CANDIDATE_RANKING.csv`
- `FINAL_20_RANKED_CANDIDATES.csv`
- `FINAL_5_STRONGEST_CANDIDATES.csv`

## Recommended manuscript wording for the duplicate

"During source-level integrity auditing, one GEO sample record (GSM8546673) was found to contain matrix, barcode, and feature content identical to GSM8546664. The duplicate record was retained in the provenance audit but excluded from downstream biological inference, resulting in 27 independent analysis samples from the original 28 GEO records."

This wording should only be used after confirming the final Stage-05A output contains the corresponding decision.

## Recommended manuscript wording for the common window

"Because the experimental sampling grids were not identical across experiments, reproducible reversal analysis was restricted to the shared 0-, 24-, and 72-hpi timepoints. After exclusion of the source-level duplicate, three complete independent trajectories were available for strict validation."

## Final caution

Stage 31 is a reproducibility and provenance gate, not a substitute for biological interpretation or expert review. The final five candidates remain candidates for manual scientific review unless the manuscript contains additional experimental or literature evidence.
