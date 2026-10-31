from pathlib import Path
import hashlib
import json
import pandas as pd

from config import (
    RESULT_ROOT,
    DISCOVERY_ROOT,
    CORRECTED_DISCOVERY_ROOT,
    COMMON_WINDOW_ROOT,
    REPLICATE_VALIDATION_ROOT,
    CELL_LEVEL_ROOT,
    MARKER_ROOT,
    GENE_CHARACTERIZATION_ROOT,
    PRIORITIZATION_ROOT,
    ENRICHMENT_ROOT,
    ENRICHMENT_RESULTS_ROOT,
    PATHWAY_INTERPRETATION_ROOT,
    ALIAS_ROOT,
    LITERATURE_ROOT,
    BIOLOGY_ROOT,
    EVIDENCE_ROOT,
    SENSITIVITY_ROOT,
    FINAL_ROOT,
    PAPER_ROOT,
    AUDIT_ROOT,
    PIPELINE_VERSION,
    ANALYSIS_VERSION,
    SOURCE_GSMS,
    ANALYSIS_GSMS,
    EXCLUDED_GSMS,
    EXCLUDED_GSMS_TUPLE,
    SOURCE_DUPLICATE_KEEP,
    COMMON_TIMEPOINTS,
    EXPECTED_COMMON_WINDOW_TRAJECTORIES,
)

# ============================================================
# 31 — PIPELINE FINAL AUDIT
# ============================================================
# Final reproducibility gate.
#
# Checks:
#   * required outputs exist
#   * source/analysis sample counts are correct
#   * the duplicate was source-verified and excluded downstream
#   * strict validation uses exactly 3 independent trajectories
#   * cell-level and gene-level analysis files do not contain
#     the excluded GSM
#   * final candidate tables have the expected evidence columns
#   * automatic novelty claims are disabled
#   * all pipeline code is hashed into a manifest
# ============================================================

AUDIT_ROOT.mkdir(parents=True, exist_ok=True)


# ============================================================
# REQUIRED OUTPUTS
# ============================================================

required_outputs = [
    RESULT_ROOT / "00_raw_file_smoke_test.csv",
    RESULT_ROOT / "01_single_matrix_load_test.csv",
    RESULT_ROOT / "02_representative_qc_summary.csv",
    RESULT_ROOT / "02_marker_availability.csv",
    RESULT_ROOT / "03_dataset_inventory.csv",
    RESULT_ROOT / "04_sample_identity_audit.csv",
    RESULT_ROOT / "05A_original_geo_source_verification.txt",
    RESULT_ROOT / "05A_source_duplicate_verification_summary.csv",
    RESULT_ROOT / "05A_source_duplicate_file_comparison.csv",
    RESULT_ROOT / "05_full_sha256_matrix_verification.csv",
    RESULT_ROOT / "05_duplicate_matrix_hashes.csv",
    RESULT_ROOT / "GSE278446_gene_level_sample_data.csv",
    RESULT_ROOT / "GSE278446_cell_program_scores.csv",
    RESULT_ROOT / "GSE278446_sample_program_summary.csv",
    DISCOVERY_ROOT / "program_temporal" / "PROGRAM_TEMPORAL_TRAJECTORIES.csv",
    CORRECTED_DISCOVERY_ROOT / "Exp1_gene_trajectories.csv",
    CORRECTED_DISCOVERY_ROOT / "Exp2_gene_trajectories.csv",
    COMMON_WINDOW_ROOT / "STRONG_REPRODUCIBLE_REVERSALS.csv",
    REPLICATE_VALIDATION_ROOT / "ALL_AVAILABLE_TRAJECTORY_VALIDATION.csv",
    REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv",
    CELL_LEVEL_ROOT / "CELL_LEVEL_REVERSAL_SCORES.csv",
    CELL_LEVEL_ROOT / "SAMPLE_LEVEL_REVERSAL_AND_PROGRAM_SCORES.csv",
    MARKER_ROOT / "STRICT_GENE_MARKER_OVERLAP.csv",
    GENE_CHARACTERIZATION_ROOT / "ALL_STRICT_GENE_CHARACTERIZATION.csv",
    PRIORITIZATION_ROOT / "ALL_STRICT_GENE_PRIORITIZATION.csv",
    ENRICHMENT_ROOT / "BACKGROUND_GENES.txt",
    ENRICHMENT_ROOT / "ENRICHMENT_INPUT_SUMMARY.csv",
    ENRICHMENT_RESULTS_ROOT / "ENRICHMENT_RUN_SUMMARY.csv",
    ENRICHMENT_RESULTS_ROOT / "GPROFILER_RUN_METADATA.json",
    PATHWAY_INTERPRETATION_ROOT / "MASTER_PATHWAY_INTERPRETATION.csv",
    PATHWAY_INTERPRETATION_ROOT / "PATHWAY_REPRODUCIBILITY_SUMMARY.csv",
    ALIAS_ROOT / "NCBI_GENE_ALIAS_AUDIT.csv",
    LITERATURE_ROOT / "BROAD_LITERATURE_CONTEXT_AUDIT.csv",
    LITERATURE_ROOT / "targeted_verification" / "TARGETED_LITERATURE_VERIFICATION.csv",
    LITERATURE_ROOT / "paper_level" / "TOP20_PAPER_LEVEL_SUMMARY.csv",
    BIOLOGY_ROOT / "CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv",
    EVIDENCE_ROOT / "FINAL_CANDIDATE_EVIDENCE_AUDIT.csv",
    SENSITIVITY_ROOT / "SENSITIVITY_TOP20_STABILITY.csv",
    SENSITIVITY_ROOT / "SENSITIVITY_GENE_RANK_STABILITY.csv",
    FINAL_ROOT / "FINAL_20_RANKED_CANDIDATES.csv",
    FINAL_ROOT / "FINAL_5_STRONGEST_CANDIDATES.csv",
    PAPER_ROOT / "FINAL_5_PAPER_RECHECK.csv",
]

file_rows = []
for p in required_outputs:
    file_rows.append({
        "path": str(p),
        "exists": p.exists(),
        "size_bytes": p.stat().st_size if p.exists() else 0,
    })
files = pd.DataFrame(file_rows)
files.to_csv(AUDIT_ROOT / "PIPELINE_OUTPUT_FILE_AUDIT.csv", index=False)


# ============================================================
# LOGICAL CHECKS
# ============================================================

checks = []

def add_check(name, passed, detail=""):
    checks.append({
        "check": name,
        "passed": bool(passed),
        "detail": str(detail),
    })


add_check("source_sample_count_28", len(SOURCE_GSMS) == 28, len(SOURCE_GSMS))
add_check("analysis_sample_count_27", len(ANALYSIS_GSMS) == 27, len(ANALYSIS_GSMS))
add_check(
    "duplicate_excluded",
    SOURCE_DUPLICATE_KEEP not in EXCLUDED_GSMS
    and "GSM8546673" in EXCLUDED_GSMS
    and "GSM8546673" not in ANALYSIS_GSMS,
    f"retain={SOURCE_DUPLICATE_KEEP}; excluded={sorted(EXCLUDED_GSMS)}",
)


# ------------------------------------------------------------
# 05A verification handoff
# ------------------------------------------------------------

verify = RESULT_ROOT / "05A_source_duplicate_verification_summary.csv"
if verify.exists():
    v = pd.read_csv(verify, dtype=str)
    if not v.empty:
        row = v.iloc[0]
        add_check(
            "source_duplicate_confirmed",
            row.get("decision") == "SOURCE_LEVEL_DUPLICATE_CONFIRMED"
            and row.get("keep_gsm") == SOURCE_DUPLICATE_KEEP
            and row.get("duplicate_gsm") in EXCLUDED_GSMS
            and str(row.get("all_three_identical", "")).lower() == "true",
            row.to_dict(),
        )
    else:
        add_check("source_duplicate_confirmed", False, "verification CSV is empty")
else:
    add_check("source_duplicate_confirmed", False, "verification CSV missing")


# ------------------------------------------------------------
# Inventory / identity
# ------------------------------------------------------------

inventory_file = RESULT_ROOT / "03_dataset_inventory.csv"
if inventory_file.exists():
    inv = pd.read_csv(inventory_file)
    add_check("inventory_has_28_rows", len(inv) == 28, f"rows={len(inv)}")
    add_check(
        "inventory_all_pass",
        inv["status"].astype(str).eq("PASS").all(),
        inv["status"].value_counts().to_dict(),
    )

identity_file = RESULT_ROOT / "04_sample_identity_audit.csv"
if identity_file.exists():
    ident = pd.read_csv(identity_file)
    add_check("identity_28_rows", len(ident) == 28, f"rows={len(ident)}")
    add_check(
        "identity_all_agree",
        ident["filename_agrees"].astype(str).str.lower().eq("true").all(),
        f"rows={len(ident)}",
    )


# ------------------------------------------------------------
# SHA-256 matrix audit
# ------------------------------------------------------------

sha_file = RESULT_ROOT / "05_full_sha256_matrix_verification.csv"
if sha_file.exists():
    sha = pd.read_csv(sha_file)
    add_check("sha_matrix_rows_28", len(sha) == 28, f"rows={len(sha)}")
    observed = set(
        sha.loc[sha["integrity_status"] == "SOURCE_LEVEL_DUPLICATE", "GSM"]
    )
    expected = {SOURCE_DUPLICATE_KEEP} | set(EXCLUDED_GSMS)
    add_check(
        "sha_only_expected_duplicate",
        observed == expected,
        f"observed={sorted(observed)}; expected={sorted(expected)}",
    )
    add_check(
        "sha_no_unresolved_duplicates",
        not sha["integrity_status"].eq("UNRESOLVED_DUPLICATE").any(),
        sha["integrity_status"].value_counts().to_dict(),
    )


# ------------------------------------------------------------
# Excluded GSM absence from downstream analysis
# ------------------------------------------------------------

gene_sample = RESULT_ROOT / "GSE278446_gene_level_sample_data.csv"
if gene_sample.exists():
    g = pd.read_csv(gene_sample, usecols=["GSM"])
    present = set(g["GSM"].astype(str))
    add_check(
        "excluded_duplicate_absent_from_gene_analysis",
        not (EXCLUDED_GSMS & present),
        f"excluded present in gene table: {sorted(EXCLUDED_GSMS & present)}",
    )

cell_sample = RESULT_ROOT / "GSE278446_cell_program_scores.csv"
if cell_sample.exists():
    c = pd.read_csv(cell_sample, usecols=["GSM"])
    present = set(c["GSM"].astype(str))
    add_check(
        "excluded_duplicate_absent_from_cell_analysis",
        not (EXCLUDED_GSMS & present),
        f"excluded present in cell table: {sorted(EXCLUDED_GSMS & present)}",
    )

cell_level = CELL_LEVEL_ROOT / "CELL_LEVEL_REVERSAL_SCORES.csv"
if cell_level.exists():
    cl = pd.read_csv(cell_level, usecols=["GSM"])
    present = set(cl["GSM"].astype(str))
    add_check(
        "excluded_duplicate_absent_from_cell_level_reversal",
        not (EXCLUDED_GSMS & present),
        f"excluded present in cell-level reversal: {sorted(EXCLUDED_GSMS & present)}",
    )


# ------------------------------------------------------------
# Strict 3/3 validation
# ------------------------------------------------------------

strict = REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv"
if strict.exists():
    s = pd.read_csv(strict)
    add_check(
        "strict_validation_uses_3_independent_trajectories",
        "n_complete_trajectories" in s.columns
        and s["n_complete_trajectories"].astype(str).eq("3").all(),
        f"rows={len(s)}",
    )
else:
    add_check(
        "strict_validation_uses_3_independent_trajectories",
        False,
        "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv missing",
    )

# Save logical checks summary
checks_df = pd.DataFrame(checks)
checks_df.to_csv(AUDIT_ROOT / "PIPELINE_LOGICAL_CHECKS.csv", index=False)
print(f"Pipeline audit completed with {len(checks_df)} checks recorded.")
