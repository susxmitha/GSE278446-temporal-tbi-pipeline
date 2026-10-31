import time
import requests
import pandas as pd
from config import *

# ============================================================
# 25 — INDEPENDENT LITERATURE RECHECK
# ============================================================
# Purpose:
#   Independently re-query PubMed for the TOP-20 candidates
#   using alternative wording and compare the results against
#   Stage 23 targeted literature verification.
#
# IMPORTANT:
#   This stage is a literature-search robustness check.
#   It does NOT assign biological novelty.
#
#   Zero literature hits are reported as zero retrieved evidence,
#   NOT as proof of novelty.
#
# Improvements over the previous version:
#   1. Records independent PubMed IDs.
#   2. Records search errors separately from true zero hits.
#   3. Compares hit/no-hit status.
#   4. Calculates PMID overlap between Stage 23 and Stage 25.
#   5. Preserves the exact query used.
#   6. Produces a transparent concordance summary.
# ============================================================

INPUT = PRIORITIZATION_ROOT / "TOP_20_PRIORITY_GENES.csv"

TARGET_INPUT = (
    LITERATURE_ROOT
    / "targeted_verification"
    / "TARGETED_LITERATURE_VERIFICATION.csv"
)

OUT = LITERATURE_ROOT / "independent_recheck"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 90)
print("25 — INDEPENDENT LITERATURE RECHECK")
print("=" * 90)


# ------------------------------------------------------------
# LOAD INPUTS
# ------------------------------------------------------------

if not INPUT.exists():
    raise FileNotFoundError(
        f"TOP-20 prioritization file not found:\n{INPUT}"
    )

if not TARGET_INPUT.exists():
    raise FileNotFoundError(
        f"Stage 23 output not found:\n{TARGET_INPUT}"
    )

df = pd.read_csv(INPUT)
target = pd.read_csv(TARGET_INPUT)

if "gene" not in df.columns:
    raise RuntimeError(
        "TOP_20_PRIORITY_GENES.csv does not contain 'gene'."
    )

if "gene" not in target.columns:
    raise RuntimeError(
        "TARGETED_LITERATURE_VERIFICATION.csv does not contain 'gene'."
    )

if "targeted_pmids" not in target.columns:
    raise RuntimeError(
        "Stage 23 output does not contain 'targeted_pmids'."
    )

df["gene"] = (
    df["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
)

target["gene"] = (
    target["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
)

print(f"Genes loaded: {len(df)}")


# ------------------------------------------------------------
# PUBMED SEARCH FUNCTION
# ------------------------------------------------------------

def search_pubmed(query, retmax=50):

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "tool": NCBI_TOOL
    }

    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL

    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params,
        timeout=30
    )

    response.raise_for_status()

    result = response.json()["esearchresult"]

    count = int(
        result.get("count", 0)
    )

    ids = result.get(
        "idlist",
        []
    )

    return count, ids


# ------------------------------------------------------------
# ALTERNATIVE QUERY
# ------------------------------------------------------------

def build_query(gene):

    return (
        f'"{gene}"[Title/Abstract] '
        f'AND '
        f'(zebrafish[Title/Abstract] OR '
        f'"Danio rerio"[Title/Abstract]) '
        f'AND '
        f'(brain[Title/Abstract] OR '
        f'telencephalon[Title/Abstract] OR '
        f'forebrain[Title/Abstract]) '
        f'AND '
        f'(injury[Title/Abstract] OR '
        f'lesion[Title/Abstract] OR '
        f'trauma[Title/Abstract] OR '
        f'regeneration[Title/Abstract])'
    )


# ------------------------------------------------------------
# RUN INDEPENDENT SEARCH
# ------------------------------------------------------------

rows = []

genes = (
    df["gene"]
    .drop_duplicates()
    .tolist()
)

for i, gene in enumerate(genes, start=1):

    query = build_query(gene)

    print()
    print(
        f"[{i:02d}/{len(genes)}] {gene}"
    )

    print(
        f"    Query: {query}"
    )

    try:

        count, ids = search_pubmed(
            query,
            retmax=50
        )

        ids = list(
            dict.fromkeys(ids)
        )

        status = "SUCCESS"

        print(
            f"    PubMed count: {count}"
        )

        print(
            f"    Retrieved PMIDs: {len(ids)}"
        )

    except Exception as exc:

        count = None
        ids = []
        status = "ERROR"

        print(
            f"    ERROR: {exc}"
        )

    rows.append({
        "gene": gene,
        "independent_query": query,
        "independent_context_count": count,
        "independent_retrieved_pmids": len(ids),
        "independent_pmids": ";".join(ids),
        "search_status": status
    })

    # Small delay for API courtesy.
    if i < len(genes):
        time.sleep(0.25)


independent = pd.DataFrame(rows)


# ------------------------------------------------------------
# SAVE RAW INDEPENDENT SEARCH
# ------------------------------------------------------------

independent_output = (
    OUT /
    "INDEPENDENT_LITERATURE_RECHECK.csv"
)

independent.to_csv(
    independent_output,
    index=False
)


# ------------------------------------------------------------
# PREPARE STAGE 23 DATA
# ------------------------------------------------------------

target_small = target[
    [
        "gene",
        "targeted_unique_pmids",
        "targeted_pmids"
    ]
].copy()

target_small = (
    target_small
    .drop_duplicates(
        subset=["gene"]
    )
)


# ------------------------------------------------------------
# PMID OVERLAP FUNCTION
# ------------------------------------------------------------

def parse_pmids(value):

    if pd.isna(value):
        return set()

    text = str(value).strip()

    if not text:
        return set()

    if text.lower() in {
        "nan",
        "none"
    }:
        return set()

    return {
        x.strip()
        for x in text.split(";")
        if x.strip()
    }


# ------------------------------------------------------------
# MERGE STAGE 23 + STAGE 25
# ------------------------------------------------------------

merged = target_small.merge(
    independent,
    on="gene",
    how="outer"
)


# ------------------------------------------------------------
# CONCORDANCE ANALYSIS
# ------------------------------------------------------------

overlap_counts = []
concordance_status = []
target_status = []
independent_status = []

for _, row in merged.iterrows():

    targeted_ids = parse_pmids(
        row.get("targeted_pmids", "")
    )

    independent_ids = parse_pmids(
        row.get("independent_pmids", "")
    )

    overlap = (
        targeted_ids
        .intersection(
            independent_ids
        )
    )

    overlap_counts.append(
        len(overlap)
    )

    targeted_n = len(
        targeted_ids
    )

    independent_n = len(
        independent_ids
    )

    target_status.append(
        "hit"
        if targeted_n > 0
        else "zero_retrieved_hits"
    )

    if independent_n > 0:

        independent_status.append(
            "hit"
        )

    elif row.get("search_status") == "ERROR":

        independent_status.append(
            "search_error"
        )

    else:

        independent_status.append(
            "zero_retrieved_hits"
        )

    # --------------------------------------------------------
    # SEARCH CONCORDANCE
    # --------------------------------------------------------

    if row.get("search_status") == "ERROR":

        status = "NOT_ASSESSED_SEARCH_ERROR"

    elif targeted_n > 0 and independent_n > 0:

        if len(overlap) > 0:
            status = "CONCORDANT_WITH_PAPER_OVERLAP"
        else:
            status = "CONCORDANT_HIT_STATUS_NO_PAPER_OVERLAP"

    elif targeted_n == 0 and independent_n == 0:

        status = "CONCORDANT_ZERO_HIT"

    else:

        status = "DISCORDANT_HIT_STATUS"

    concordance_status.append(
        status
    )


merged["targeted_hit_status"] = (
    target_status
)

merged["independent_hit_status"] = (
    independent_status
)

merged["shared_pmids"] = [
    len(
        parse_pmids(
            row.get("targeted_pmids", "")
        )
        &
        parse_pmids(
            row.get("independent_pmids", "")
        )
    )
    for _, row in merged.iterrows()
]

merged["search_concordance"] = (
    concordance_status
)


# ------------------------------------------------------------
# OVERLAP PERCENTAGE
# ------------------------------------------------------------

def calculate_overlap(row):

    targeted = parse_pmids(
        row.get("targeted_pmids", "")
    )

    independent_ids = parse_pmids(
        row.get("independent_pmids", "")
    )

    if not targeted or not independent_ids:
        return 0.0

    union = targeted.union(
        independent_ids
    )

    if not union:
        return 0.0

    return (
        len(
            targeted.intersection(
                independent_ids
            )
        )
        / len(union)
        * 100.0
    )


merged["pmid_overlap_percent"] = (
    merged.apply(
        calculate_overlap,
        axis=1
    )
)


# ------------------------------------------------------------
# SAVE CONCORDANCE TABLE
# ------------------------------------------------------------

concordance_output = (
    OUT /
    "LITERATURE_RECHECK_CONCORDANCE.csv"
)

merged.to_csv(
    concordance_output,
    index=False
)


# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

successful = merged[
    merged["search_status"] == "SUCCESS"
]

summary = {
    "genes_evaluated": len(merged),

    "successful_searches": len(
        successful
    ),

    "search_errors": int(
        (merged["search_status"] == "ERROR")
        .sum()
    ),

    "targeted_stage23_hits": int(
        (
            pd.to_numeric(
                merged["targeted_unique_pmids"],
                errors="coerce"
            ).fillna(0)
            > 0
        ).sum()
    ),

    "independent_stage25_hits": int(
        (
            pd.to_numeric(
                merged["independent_context_count"],
                errors="coerce"
            ).fillna(0)
            > 0
        ).sum()
    ),

    "concordant_zero_hit": int(
        (
            merged["search_concordance"]
            == "CONCORDANT_ZERO_HIT"
        ).sum()
    ),

    "concordant_with_paper_overlap": int(
        (
            merged["search_concordance"]
            == "CONCORDANT_WITH_PAPER_OVERLAP"
        ).sum()
    ),

    "concordant_hit_without_overlap": int(
        (
            merged["search_concordance"]
            == "CONCORDANT_HIT_STATUS_NO_PAPER_OVERLAP"
        ).sum()
    ),

    "discordant_hit_status": int(
        (
            merged["search_concordance"]
            == "DISCORDANT_HIT_STATUS"
        ).sum()
    )
}

summary_df = pd.DataFrame(
    [summary]
)

summary_output = (
    OUT /
    "LITERATURE_RECHECK_SUMMARY.csv"
)

summary_df.to_csv(
    summary_output,
    index=False
)


# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

print()
print("=" * 90)
print("25 — INDEPENDENT LITERATURE RECHECK COMPLETE")
print("=" * 90)

print(
    f"Genes evaluated                 : "
    f"{summary['genes_evaluated']}"
)

print(
    f"Successful searches             : "
    f"{summary['successful_searches']}"
)

print(
    f"Search errors                   : "
    f"{summary['search_errors']}"
)

print(
    f"Stage-23 candidates with hits   : "
    f"{summary['targeted_stage23_hits']}"
)

print(
    f"Stage-25 candidates with hits   : "
    f"{summary['independent_stage25_hits']}"
)

print(
    f"Concordant zero-hit candidates  : "
    f"{summary['concordant_zero_hit']}"
)

print(
    f"Concordant + PMID overlap       : "
    f"{summary['concordant_with_paper_overlap']}"
)

print(
    f"Concordant hit/no overlap       : "
    f"{summary['concordant_hit_without_overlap']}"
)

print(
    f"Discordant hit status           : "
    f"{summary['discordant_hit_status']}"
)

print()
print("Outputs:")
print(independent_output)
print(concordance_output)
print(summary_output)

print()
print("IMPORTANT:")
print(
    "This stage evaluates literature-search robustness only."
)
print(
    "Literature-search absence is NOT interpreted as proof "
    "of biological novelty."
)
