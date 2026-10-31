import pandas as pd
import numpy as np

from config import (
    PRIORITIZATION_ROOT,
    BIOLOGY_ROOT,
    EVIDENCE_ROOT,
    SENSITIVITY_ROOT,
    FINAL_ROOT,
)

# ============================================================
# 29 — FINAL CANDIDATE SYNTHESIS
# ============================================================
# No new biological computation and no new literature search.
# Combines already-audited evidence into a transparent review queue.
#
# IMPORTANT:
#   The definitive temporal evidence is STRICT 3/3 AVAILABLE,
#   because one original GEO GSM is a confirmed source duplicate.
#
# File contract used here:
#   Stage 16 → PRIORITIZATION_ROOT / ALL_STRICT_GENE_PRIORITIZATION.csv
#              (must contain: gene, temporal_score)
#   Stage 26 → BIOLOGY_ROOT / CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv
#              (must contain: gene, biological_context)
#   Stage 27 → EVIDENCE_ROOT / FINAL_CANDIDATE_EVIDENCE_AUDIT.csv
#              (must contain: gene, audit_class,
#               strict_3of3_available_evidence, direct_context_matches,
#               direct_context_evidence)
#   Stage 28 → SENSITIVITY_ROOT / SENSITIVITY_GENE_RANK_STABILITY.csv
#              (must contain: gene, mean_rank, rank_sd, rank_range)
# ============================================================

OUT = FINAL_ROOT
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD
# ============================================================

priority = pd.read_csv(
    PRIORITIZATION_ROOT / "ALL_STRICT_GENE_PRIORITIZATION.csv"
)
bio = pd.read_csv(
    BIOLOGY_ROOT / "CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv"
)
audit = pd.read_csv(
    EVIDENCE_ROOT / "FINAL_CANDIDATE_EVIDENCE_AUDIT.csv"
)
sens = pd.read_csv(
    SENSITIVITY_ROOT / "SENSITIVITY_GENE_RANK_STABILITY.csv"
)


# ============================================================
# SCHEMA CHECKS — fail loudly, do not silently NaN-fill
# ============================================================

def require_columns(df, required, name):
    missing = set(required) - set(df.columns)
    if missing:
        raise RuntimeError(
            f"{name} is missing required columns: {sorted(missing)}\n"
            f"Columns present: {list(df.columns)}"
        )

require_columns(
    priority, {"gene", "temporal_score"}, "ALL_STRICT_GENE_PRIORITIZATION.csv"
)
require_columns(
    bio, {"gene", "biological_context"},
    "CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv",
)
require_columns(
    audit,
    {
        "gene",
        "audit_class",
        "strict_3of3_available_evidence",
        "direct_context_matches",
        "direct_context_evidence",
    },
    "FINAL_CANDIDATE_EVIDENCE_AUDIT.csv",
)
require_columns(
    sens, {"gene", "mean_rank", "rank_sd", "rank_range"},
    "SENSITIVITY_GENE_RANK_STABILITY.csv",
)


# ============================================================
# NORMALISE KEYS
# ============================================================

def norm_gene(s):
    return s.astype(str).str.strip().str.lower()

priority["gene"] = norm_gene(priority["gene"])
bio["gene"] = norm_gene(bio["gene"])
audit["gene"] = norm_gene(audit["gene"])
sens["gene"] = norm_gene(sens["gene"])

bio = bio[["gene", "biological_context"]].drop_duplicates("gene")
audit = audit[[
    "gene",
    "audit_class",
    "strict_3of3_available_evidence",
    "direct_context_matches",
    "direct_context_evidence",
]].drop_duplicates("gene")
sens = sens[[
    "gene", "mean_rank", "rank_sd", "rank_range"
]].drop_duplicates("gene")


# ============================================================
# MERGE
# ============================================================

df = (
    priority
    .merge(bio, on="gene", how="left")
    .merge(audit, on="gene", how="left")
    .merge(sens, on="gene", how="left")
)


# ============================================================
# BOOLEAN NORMALISATION
# ============================================================

def to_bool(s):
    mapping = {
        "true": True, "1": True, "yes": True, "y": True, "t": True,
        "false": False, "0": False, "no": False, "n": False, "f": False,
        "": False, "nan": False, "none": False, "null": False,
    }
    return (
        s.astype("string")
         .str.strip()
         .str.lower()
         .map(mapping)
         .fillna(False)
         .astype(bool)
    )

df["strict_3of3_available_evidence"] = to_bool(
    df["strict_3of3_available_evidence"]
)
df["direct_context_evidence"] = to_bool(
    df["direct_context_evidence"]
)
df["direct_context_matches"] = (
    pd.to_numeric(df["direct_context_matches"], errors="coerce")
      .fillna(0)
      .astype(int)
)


# ============================================================
# SCORES
# ============================================================

df["temporal_evidence"] = df["strict_3of3_available_evidence"].astype(int)
df["literature_evidence"] = df["direct_context_evidence"].astype(int)

rank_sd = pd.to_numeric(df["rank_sd"], errors="coerce")
df["robustness_score"] = 1.0 / (1.0 + rank_sd)
df["robustness_score"] = (
    df["robustness_score"]
    .replace([np.inf, -np.inf], np.nan)
    .fillna(0)
)

df["final_review_score"] = (
    0.80 * pd.to_numeric(df["temporal_score"], errors="coerce").fillna(0)
    + 0.20 * df["robustness_score"]
)

df = df.sort_values(
    ["final_review_score", "gene"],
    ascending=[False, True],
).reset_index(drop=True)

df["final_rank"] = np.arange(1, len(df) + 1)


# ============================================================
# REVIEW RECOMMENDATION
# ============================================================

df["review_recommendation"] = np.select(
    [
        df["temporal_evidence"].eq(1) & df["literature_evidence"].eq(1),
        df["temporal_evidence"].eq(1),
    ],
    [
        "HIGH_PRIORITY_MANUAL_REVIEW",
        "TEMPORAL_PRIORITY_MANUAL_REVIEW",
    ],
    default="SECONDARY_REVIEW",
)

df["novelty_claim_permitted"] = False


# ============================================================
# WRITE
# ============================================================

df.to_csv(OUT / "FINAL_CANDIDATE_RANKING.csv", index=False)
df.head(20).to_csv(OUT / "FINAL_20_RANKED_CANDIDATES.csv", index=False)
df.head(10).to_csv(OUT / "FINAL_10_MANUAL_REVIEW.csv", index=False)
df.head(5).to_csv(OUT / "FINAL_5_STRONGEST_CANDIDATES.csv", index=False)

pd.DataFrame([{
    "final_20": min(20, len(df)),
    "final_10": min(10, len(df)),
    "final_5": min(5, len(df)),
    "strict_3of3_available_candidates": int(df["temporal_evidence"].sum()),
    "novelty_claims_automatically_supported": False,
}]).to_csv(OUT / "FINAL_SYNTHESIS_SUMMARY.csv", index=False)


# ============================================================
# REPORT
# ============================================================

print("=" * 78)
print("29 — FINAL CANDIDATE SYNTHESIS COMPLETE")
print(df.head(10)[
    ["final_rank", "gene", "final_review_score", "review_recommendation"]
].to_string(index=False))
