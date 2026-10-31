import pandas as pd
import numpy as np

from config import (
    BIOLOGY_ROOT,
    REPLICATE_VALIDATION_ROOT,
    LITERATURE_ROOT,
    EVIDENCE_ROOT,
)

# ============================================================
# 27 — CANDIDATE EVIDENCE AUDIT
# ============================================================
# Checks that every TOP-20 candidate has traceable temporal
# evidence and, where available, paper-level literature evidence.
#
# Missing evidence is never converted into positive evidence.
# ============================================================

OUT = EVIDENCE_ROOT
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD
# ============================================================

candidate_file = (
    BIOLOGY_ROOT / "CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv"
)
strict_file = (
    REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv"
)
paper_file = (
    LITERATURE_ROOT / "paper_level" / "TOP20_PAPER_LEVEL_SUMMARY.csv"
)

if not candidate_file.exists():
    raise FileNotFoundError(
        f"Stage-26 output not found:\n{candidate_file}"
    )
if not strict_file.exists():
    raise FileNotFoundError(
        f"Stage-11 strict output not found:\n{strict_file}"
    )

cand = pd.read_csv(candidate_file)
strict = pd.read_csv(strict_file)
paper = pd.read_csv(paper_file) if paper_file.exists() else pd.DataFrame()


# ============================================================
# NORMALISE
# ============================================================

if "gene" not in cand.columns or "gene" not in strict.columns:
    raise RuntimeError(
        "Candidate or strict validation file lacks the 'gene' column."
    )

cand["gene"] = cand["gene"].astype(str).str.strip().str.lower()
strict["gene"] = strict["gene"].astype(str).str.strip().str.lower()

strict_genes = set(strict["gene"])
cand["strict_3of3_available_evidence"] = cand["gene"].isin(strict_genes)


def parse_bool_series(s):
    """Robustly convert a Series of mixed bool / string values to bool."""
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)

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


# ============================================================
# LITERATURE EVIDENCE
# ============================================================

if not paper.empty:

    required = {
        "gene",
        "papers_reviewed",
        "direct_context_matches",
        "direct_context_evidence",
    }
    missing = required - set(paper.columns)
    if missing:
        raise RuntimeError(
            f"Paper-level summary missing columns: {sorted(missing)}"
        )

    paper["gene"] = paper["gene"].astype(str).str.strip().str.lower()
    paper = paper.drop_duplicates("gene")

    cand = cand.merge(
        paper[[
            "gene",
            "papers_reviewed",
            "direct_context_matches",
            "direct_context_evidence",
        ]],
        on="gene",
        how="left",
    )
else:
    cand["papers_reviewed"] = 0
    cand["direct_context_matches"] = 0
    cand["direct_context_evidence"] = False


cand["papers_reviewed"] = (
    pd.to_numeric(cand["papers_reviewed"], errors="coerce")
      .fillna(0)
      .astype(int)
)
cand["direct_context_matches"] = (
    pd.to_numeric(cand["direct_context_matches"], errors="coerce")
      .fillna(0)
      .astype(int)
)
cand["direct_context_evidence"] = parse_bool_series(
    cand["direct_context_evidence"]
)


# ============================================================
# AUDIT CLASSIFICATION
# ============================================================

cand["evidence_completeness"] = (
    cand["strict_3of3_available_evidence"].astype(int)
    + cand["direct_context_evidence"].astype(int)
)

cand["audit_class"] = np.select(
    [
        cand["strict_3of3_available_evidence"]
        & cand["direct_context_evidence"],
        cand["strict_3of3_available_evidence"],
    ],
    [
        "temporal_plus_literature_traceable",
        "temporal_traceable_only",
    ],
    default="insufficient_traceable_evidence",
)

cand["novelty_claim_permitted"] = False


# ============================================================
# WRITE
# ============================================================

cand.to_csv(OUT / "FINAL_CANDIDATE_EVIDENCE_AUDIT.csv", index=False)

cand[
    cand["audit_class"] == "temporal_plus_literature_traceable"
].to_csv(OUT / "HIGH_EVIDENCE_CANDIDATES.csv", index=False)

pd.DataFrame([{
    "candidate_count": len(cand),
    "strict_3of3_count": int(
        cand["strict_3of3_available_evidence"].sum()
    ),
    "direct_literature_count": int(
        cand["direct_context_evidence"].sum()
    ),
    "temporal_plus_literature_count": int(
        (cand["audit_class"] == "temporal_plus_literature_traceable").sum()
    ),
    "novelty_claims_automatically_supported": False,
}]).to_csv(OUT / "EVIDENCE_AUDIT_SUMMARY.csv", index=False)

print("=" * 78)
print("27 — CANDIDATE EVIDENCE AUDIT COMPLETE")
print(f"Candidates audited: {len(cand)}")
print(f"Output: {OUT}")
