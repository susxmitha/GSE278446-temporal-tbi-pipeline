"""
26 — CANDIDATE BIOLOGICAL CONTEXT CLASSIFICATION

This stage classifies the TOP-20 candidates using only explicit
pathway evidence produced by Stage 19. It deliberately refuses to
infer cell type or biological context from pathway names alone.

Contract with Stage 19:
  MASTER_PATHWAY_INTERPRETATION.csv columns include:
    input_set, source, term_id, term_name, description, p_value,
    intersection_size, term_size, query_size, effective_domain_size,
    precision, recall, significant, intersection_genes, p_value_numeric

  There is NO "status" column in this file. The old version of this
  script filtered on "status" == "SUCCESS", which silently produced
  an empty table and mislabelled every candidate as UNRESOLVED.

If no pathway evidence is available, the candidate remains:
    biological_context = UNRESOLVED

No biological interpretation is fabricated.
"""

from pathlib import Path
import sys
import traceback
import pandas as pd

from config import (
    PRIORITIZATION_ROOT,
    PATHWAY_INTERPRETATION_ROOT,
    BIOLOGY_ROOT,
)


MASTER_PATHWAY_FILE = (
    PATHWAY_INTERPRETATION_ROOT / "MASTER_PATHWAY_INTERPRETATION.csv"
)
PRIORITY_FILE = PRIORITIZATION_ROOT / "TOP_20_PRIORITY_GENES.csv"
OUTPUT_ROOT = BIOLOGY_ROOT
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def safe_read_csv(path):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception as exc:
        print(f"WARNING: Could not read {path}: {exc}")
        return pd.DataFrame()


def load_candidates():
    if not PRIORITY_FILE.exists():
        raise FileNotFoundError(
            f"TOP-20 priority file not found:\n{PRIORITY_FILE}"
        )
    df = pd.read_csv(PRIORITY_FILE)
    if "gene" not in df.columns:
        raise RuntimeError(
            "TOP_20_PRIORITY_GENES.csv lacks the 'gene' column."
        )
    df["gene"] = df["gene"].astype(str).str.strip().str.lower()
    df = (
        df[df["gene"].ne("") & df["gene"].ne("nan")]
        .drop_duplicates("gene")
        .copy()
    )
    if df.empty:
        raise RuntimeError(
            "TOP_20_PRIORITY_GENES.csv contains no usable genes."
        )
    return df


def pathway_table_status(pathway_df):
    """
    Return (available, reason). 'available' means Stage 19 produced
    at least one pathway row with a usable term_id.
    """
    if pathway_df.empty:
        return False, "MASTER_PATHWAY_INTERPRETATION.csv is empty or missing"

    if "term_id" not in pathway_df.columns:
        return False, "MASTER_PATHWAY_INTERPRETATION.csv lacks 'term_id'"

    usable = pathway_df[
        pathway_df["term_id"].astype(str).str.strip().ne("")
    ]
    if usable.empty:
        return False, "No non-empty term_id rows in pathway table"

    return True, "PATHWAY_TABLE_USABLE"


def build_supporting_pathways(pathway_df, max_chars=5000):
    """
    Concatenate observed pathway term names. This is list-level
    evidence about the enrichment run. It is NOT a gene-specific
    assignment and must not be interpreted as one.
    """
    if pathway_df.empty or "term_name" not in pathway_df.columns:
        return ""

    names = [
        str(v).strip()
        for v in pathway_df["term_name"].dropna()
        if str(v).strip()
    ]
    seen = list(dict.fromkeys(names))
    text = "; ".join(seen)
    return text[:max_chars]


# ============================================================
# CLASSIFICATION
# ============================================================

def classify_candidates(candidates, pathway_df):

    pathway_available, pathway_reason = pathway_table_status(pathway_df)
    supporting = (
        build_supporting_pathways(pathway_df)
        if pathway_available else ""
    )

    rows = []
    for rank, gene in enumerate(candidates["gene"], start=1):

        # Pathway evidence in this pipeline is list-level. There is
        # no verified gene-to-term mapping available at this stage,
        # so we do not assign a cell-type / biological context from
        # pathway names alone.
        if not pathway_available:
            context = "UNRESOLVED"
            status = "INSUFFICIENT_PATHWAY_EVIDENCE"
            evidence = "UNAVAILABLE"
        else:
            context = "UNRESOLVED"
            status = (
                "PATHWAY_EVIDENCE_AVAILABLE_BUT_NO_EXPLICIT_GENE_MAPPING"
            )
            evidence = "AVAILABLE"

        rows.append({
            "rank_in_top20": rank,
            "gene": gene,
            "biological_context": context,
            "biological_context_class": context,
            "classification_status": status,
            "pathway_evidence": evidence,
            "supporting_pathways": supporting,
            "pathway_table_reason": pathway_reason,
        })

    return pd.DataFrame(rows)


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 90)
    print("26 — CANDIDATE BIOLOGICAL CONTEXT CLASSIFICATION")
    print("=" * 90)

    candidates = load_candidates()
    pathway_df = safe_read_csv(MASTER_PATHWAY_FILE)

    print(f"TOP-20 candidates loaded : {len(candidates)}")
    print(f"Pathway rows loaded      : {len(pathway_df)}")

    result = classify_candidates(candidates, pathway_df)

    result.to_csv(
        OUTPUT_ROOT / "CANDIDATE_BIOLOGICAL_CONTEXT_CLASSIFICATION.csv",
        index=False,
    )

    # Emit a transparent copy of the pathway evidence used, using
    # only columns Stage 19 actually writes.
    evidence_columns = [
        "input_set", "source", "term_id", "term_name",
        "p_value", "intersection_size",
    ]

    if pathway_df.empty:
        evidence_df = pd.DataFrame(columns=evidence_columns)
    else:
        evidence_df = pathway_df.copy()
        for col in evidence_columns:
            if col not in evidence_df.columns:
                evidence_df[col] = ""
        evidence_df = evidence_df[evidence_columns]

    evidence_df.to_csv(
        OUTPUT_ROOT / "CANDIDATE_PATHWAY_EVIDENCE.csv",
        index=False,
    )

    resolved_context_count = int(
        (result["biological_context"] != "UNRESOLVED").sum()
    )

    pd.DataFrame([{
        "candidate_count": len(result),
        "pathway_rows_available": len(pathway_df),
        "resolved_context_count": resolved_context_count,
        "automatic_novelty_claim_supported": False,
    }]).to_csv(
        OUTPUT_ROOT / "BIOLOGICAL_CONTEXT_SUMMARY.csv",
        index=False,
    )

    print(f"Candidates evaluated     : {len(result)}")
    print(f"Resolved biological ctx  : {resolved_context_count}")
    print(f"Output                   : {OUTPUT_ROOT}")
    print("No biological context was fabricated.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
