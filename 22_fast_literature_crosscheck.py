import time
import requests
import pandas as pd
from config import *

# ============================================================
# 22 — FAST INDEPENDENT LITERATURE CROSS-CHECK
# ============================================================
#
# Purpose:
#   Independently cross-check the Stage-21 PubMed literature
#   audit for ALL TOP-50 priority genes.
#
# Important:
#   This stage does NOT claim novelty.
#   Zero PubMed hits are reported only as absence of
#   retrieved evidence under the stated query.
#
# Improvements over the earlier version:
#   1. Audits ALL 50 genes, not only the first 20.
#   2. Requires zebrafish context for brain/TBI queries.
#   3. Records exact queries used.
#   4. Records returned PMIDs.
#   5. Records query errors explicitly.
#   6. Compares Stage 21 and Stage 22 results.
#   7. Does not fabricate or infer novelty.
# ============================================================

INPUT = LITERATURE_ROOT / "BROAD_LITERATURE_CONTEXT_AUDIT.csv"

OUT = LITERATURE_ROOT / "fast_crosscheck"
OUT.mkdir(parents=True, exist_ok=True)

if not INPUT.exists():
    raise FileNotFoundError(
        f"Stage-21 literature audit not found:\n{INPUT}"
    )

df = pd.read_csv(INPUT)

if "gene" not in df.columns:
    raise RuntimeError(
        "Required column 'gene' is missing from Stage-21 literature audit."
    )

# ------------------------------------------------------------
# Clean and validate gene list
# ------------------------------------------------------------

genes = (
    df["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
)

genes = genes[genes.notna()]
genes = genes[genes != ""]
genes = genes.drop_duplicates().tolist()

if len(genes) == 0:
    raise RuntimeError("No genes were found in the Stage-21 audit.")

print("=" * 90)
print("22 — FAST INDEPENDENT LITERATURE CROSS-CHECK")
print("=" * 90)

print(f"Input file : {INPUT}")
print(f"Genes loaded: {len(genes)}")

if len(genes) != 50:
    print(
        f"WARNING: Expected 50 TOP-priority genes, "
        f"but found {len(genes)}."
    )

# ------------------------------------------------------------
# PubMed query function
# ------------------------------------------------------------

def pubmed_count(query):

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 20,
        "tool": NCBI_TOOL,
    }

    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL

    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params,
        timeout=30,
    )

    response.raise_for_status()

    result = response.json()["esearchresult"]

    return (
        int(result.get("count", 0)),
        result.get("idlist", []),
    )


# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

rows = []

for index, gene in enumerate(genes, start=1):

    print("-" * 90)
    print(f"[{index:02d}/{len(genes)}] {gene}")

    # --------------------------------------------------------
    # Query 1:
    # Brain / telencephalon + zebrafish
    # --------------------------------------------------------

    q_brain = (
        f'"{gene}"[Title/Abstract] '
        f'AND ("Danio rerio"[Title/Abstract] OR zebrafish[Title/Abstract]) '
        f'AND (brain[Title/Abstract] OR telencephalon[Title/Abstract])'
    )

    # --------------------------------------------------------
    # Query 2:
    # TBI / brain injury + zebrafish
    # --------------------------------------------------------

    q_tbi = (
        f'"{gene}"[Title/Abstract] '
        f'AND ("Danio rerio"[Title/Abstract] OR zebrafish[Title/Abstract]) '
        f'AND ('
        f'"traumatic brain injury"[Title/Abstract] '
        f'OR "brain injury"[Title/Abstract] '
        f'OR trauma[Title/Abstract]'
        f')'
    )

    # --------------------------------------------------------
    # Query 3:
    # Regeneration + zebrafish
    # --------------------------------------------------------

    q_regeneration = (
        f'"{gene}"[Title/Abstract] '
        f'AND ("Danio rerio"[Title/Abstract] OR zebrafish[Title/Abstract]) '
        f'AND ('
        f'regeneration[Title/Abstract] '
        f'OR regenerative[Title/Abstract]'
        f')'
    )

    record = {
        "gene": gene,
        "stage22_brain_zebrafish_query": q_brain,
        "stage22_tbi_zebrafish_query": q_tbi,
        "stage22_regeneration_zebrafish_query": q_regeneration,
    }

    # --------------------------------------------------------
    # Brain / zebrafish
    # --------------------------------------------------------

    try:

        brain_count, brain_pmids = pubmed_count(q_brain)

        record["crosscheck_brain_zebrafish_count"] = brain_count
        record["crosscheck_brain_zebrafish_pmids"] = (
            ";".join(brain_pmids[:10])
        )
        record["crosscheck_brain_zebrafish_error"] = ""

        print(
            f"    brain + zebrafish: {brain_count} hits"
        )

    except Exception as exc:

        record["crosscheck_brain_zebrafish_count"] = None
        record["crosscheck_brain_zebrafish_pmids"] = ""
        record["crosscheck_brain_zebrafish_error"] = str(exc)

        print(
            f"    brain + zebrafish: ERROR — {exc}"
        )

    time.sleep(NCBI_DELAY_SECONDS)

    # --------------------------------------------------------
    # TBI / zebrafish
    # --------------------------------------------------------

    try:

        tbi_count, tbi_pmids = pubmed_count(q_tbi)

        record["crosscheck_tbi_zebrafish_count"] = tbi_count
        record["crosscheck_tbi_zebrafish_pmids"] = (
            ";".join(tbi_pmids[:10])
        )
        record["crosscheck_tbi_zebrafish_error"] = ""

        print(
            f"    TBI + zebrafish: {tbi_count} hits"
        )

    except Exception as exc:

        record["crosscheck_tbi_zebrafish_count"] = None
        record["crosscheck_tbi_zebrafish_pmids"] = ""
        record["crosscheck_tbi_zebrafish_error"] = str(exc)

        print(
            f"    TBI + zebrafish: ERROR — {exc}"
        )

    time.sleep(NCBI_DELAY_SECONDS)

    # --------------------------------------------------------
    # Regeneration / zebrafish
    # --------------------------------------------------------

    try:

        regen_count, regen_pmids = pubmed_count(q_regeneration)

        record["crosscheck_regeneration_zebrafish_count"] = regen_count
        record["crosscheck_regeneration_zebrafish_pmids"] = (
            ";".join(regen_pmids[:10])
        )
        record["crosscheck_regeneration_zebrafish_error"] = ""

        print(
            f"    regeneration + zebrafish: {regen_count} hits"
        )

    except Exception as exc:

        record["crosscheck_regeneration_zebrafish_count"] = None
        record["crosscheck_regeneration_zebrafish_pmids"] = ""
        record["crosscheck_regeneration_zebrafish_error"] = str(exc)

        print(
            f"    regeneration + zebrafish: ERROR — {exc}"
        )

    time.sleep(NCBI_DELAY_SECONDS)

    # --------------------------------------------------------
    # Independent evidence status
    # --------------------------------------------------------

    counts = [
        record.get("crosscheck_brain_zebrafish_count"),
        record.get("crosscheck_tbi_zebrafish_count"),
        record.get("crosscheck_regeneration_zebrafish_count"),
    ]

    numeric_counts = [
        x for x in counts
        if isinstance(x, (int, float))
    ]

    if not numeric_counts:

        status = "QUERY_ERROR"

    elif any(x > 0 for x in numeric_counts):

        status = "RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE"

    else:

        status = "NO_RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE"

    record["stage22_context_status"] = status

    rows.append(record)


# ------------------------------------------------------------
# Create Stage-22 dataframe
# ------------------------------------------------------------

crosscheck = pd.DataFrame(rows)

crosscheck_file = (
    OUT / "FAST_LITERATURE_CROSSCHECK.csv"
)

crosscheck.to_csv(
    crosscheck_file,
    index=False
)


# ------------------------------------------------------------
# Merge with Stage 21
# ------------------------------------------------------------

stage21 = df.copy()

stage21["gene"] = (
    stage21["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
)

merged = stage21.merge(
    crosscheck,
    on="gene",
    how="left",
    suffixes=("_stage21", "_stage22")
)


# ------------------------------------------------------------
# Compare Stage 21 vs Stage 22
# ------------------------------------------------------------

stage21_direct = pd.to_numeric(
    merged.get("direct_tbi_zebrafish_count"),
    errors="coerce"
)

stage22_direct = pd.to_numeric(
    merged.get("crosscheck_tbi_zebrafish_count"),
    errors="coerce"
)

merged["direct_tbi_count_difference"] = (
    stage22_direct - stage21_direct
)

merged["direct_tbi_count_agreement"] = (
    stage21_direct == stage22_direct
)

merged["stage21_stage22_status"] = merged.apply(
    lambda row: (
        "AGREE"
        if pd.notna(row.get("direct_tbi_zebrafish_count"))
        and pd.notna(row.get("crosscheck_tbi_zebrafish_count"))
        and row["direct_tbi_zebrafish_count"]
        == row["crosscheck_tbi_zebrafish_count"]
        else "DISAGREEMENT_OR_UNAVAILABLE"
    ),
    axis=1
)

merged_file = (
    OUT / "FAST_LITERATURE_CROSSCHECK_MERGED.csv"
)

merged.to_csv(
    merged_file,
    index=False
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_rows = []

summary_rows.append({
    "metric": "genes_audited",
    "value": len(crosscheck)
})

summary_rows.append({
    "metric": "brain_zebrafish_queries_successful",
    "value": crosscheck[
        "crosscheck_brain_zebrafish_error"
    ].eq("").sum()
})

summary_rows.append({
    "metric": "tbi_zebrafish_queries_successful",
    "value": crosscheck[
        "crosscheck_tbi_zebrafish_error"
    ].eq("").sum()
})

summary_rows.append({
    "metric": "regeneration_zebrafish_queries_successful",
    "value": crosscheck[
        "crosscheck_regeneration_zebrafish_error"
    ].eq("").sum()
})

summary_rows.append({
    "metric": "genes_with_retrieved_zebrafish_context",
    "value": (
        crosscheck["stage22_context_status"]
        == "RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE"
    ).sum()
})

summary_rows.append({
    "metric": "genes_with_no_retrieved_zebrafish_context",
    "value": (
        crosscheck["stage22_context_status"]
        == "NO_RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE"
    ).sum()
})

summary_rows.append({
    "metric": "genes_with_query_errors",
    "value": (
        crosscheck["stage22_context_status"]
        == "QUERY_ERROR"
    ).sum()
})

summary = pd.DataFrame(summary_rows)

summary_file = (
    OUT / "FAST_LITERATURE_CROSSCHECK_SUMMARY.csv"
)

summary.to_csv(
    summary_file,
    index=False
)


# ------------------------------------------------------------
# Final console report
# ------------------------------------------------------------

print()
print("=" * 90)
print("22 — FAST INDEPENDENT LITERATURE CROSS-CHECK COMPLETE")
print("=" * 90)

print(f"Genes audited                     : {len(crosscheck)}")

print(
    "Genes with retrieved context       : "
    f"{(crosscheck['stage22_context_status'] == 'RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE').sum()}"
)

print(
    "Genes with no retrieved context    : "
    f"{(crosscheck['stage22_context_status'] == 'NO_RETRIEVED_ZEBRAFISH_CONTEXT_EVIDENCE').sum()}"
)

print(
    "Genes with query errors            : "
    f"{(crosscheck['stage22_context_status'] == 'QUERY_ERROR').sum()}"
)

print()
print("Outputs:")
print(f"  Raw cross-check : {crosscheck_file}")
print(f"  Merged audit    : {merged_file}")
print(f"  Summary         : {summary_file}")

print()
print(
    "IMPORTANT:"
)
print(
    "A zero-hit PubMed result is NOT interpreted as proof of novelty."
)
print(
    "Stage 22 provides an independent literature retrieval cross-check only."
)

if len(crosscheck) == 50 and (
    crosscheck[
        [
            "crosscheck_brain_zebrafish_error",
            "crosscheck_tbi_zebrafish_error",
            "crosscheck_regeneration_zebrafish_error",
        ]
    ]
    .fillna("")
    .eq("")
    .all()
    .all()
):

    print()
    print(
        "PASS — all 50 TOP-priority genes completed all three "
        "independent PubMed cross-check queries."
    )

else:

    print()
    print(
        "REVIEW REQUIRED — not all queries completed successfully."
    )
