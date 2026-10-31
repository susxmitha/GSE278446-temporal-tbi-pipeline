from pathlib import Path
import time
import requests
import pandas as pd
from config import *

# ============================================================
# 21 — BROAD LITERATURE CONTEXT AUDIT
# ============================================================
#
# Purpose:
#   Audit the TOP-50 prioritized genes against PubMed.
#
# IMPORTANT:
#   No retrieved paper is NOT treated as proof of novelty.
#
# The audit distinguishes:
#   1. direct TBI / injury literature
#   2. regeneration literature
#   3. general zebrafish literature
#   4. no retrieved evidence
#   5. API/query failure
#
# No biological conclusion is fabricated.
# ============================================================


# ============================================================
# PATHS
# ============================================================

INPUT = (
    PRIORITIZATION_ROOT
    / "TOP_50_PRIORITY_GENES.csv"
)

OUT = LITERATURE_ROOT

OUT.mkdir(
    parents=True,
    exist_ok=True
)

ALIAS = (
    ALIAS_ROOT
    / "NCBI_GENE_ALIAS_AUDIT.csv"
)

OUTPUT_FILE = (
    OUT
    / "BROAD_LITERATURE_CONTEXT_AUDIT.csv"
)


# ============================================================
# HEADER
# ============================================================

print("=" * 90)
print("21 — BROAD LITERATURE CONTEXT AUDIT")
print("=" * 90)

print(
    f"Research root : {RESEARCH_ROOT}"
)

print(
    f"Input file    : {INPUT}"
)

print(
    f"Alias file    : {ALIAS}"
)

print(
    f"Output root   : {OUT}"
)

print("")


# ============================================================
# INPUT VALIDATION
# ============================================================

if not INPUT.exists():

    raise FileNotFoundError(
        "TOP-50 prioritization file was not found:\n"
        f"{INPUT}\n\n"
        "Stage 20/previous prioritization stages must be "
        "completed before Stage 21."
    )


df = pd.read_csv(
    INPUT
)


if "gene" not in df.columns:

    raise RuntimeError(
        "TOP_50_PRIORITY_GENES.csv does not contain "
        "the required 'gene' column.\n"
        f"Columns found: {list(df.columns)}"
    )


# ============================================================
# CLEAN TOP-50 GENES
# ============================================================

genes = (
    df["gene"]
    .astype(str)
    .str.strip()
)

genes = [
    g
    for g in genes
    if g
    and g.lower() not in {
        "nan",
        "none",
        "null",
    }
]

# Preserve original order while removing duplicates.
genes = list(
    dict.fromkeys(
        genes
    )
)

print(
    f"Genes loaded: {len(genes)}"
)

if len(genes) == 0:

    raise RuntimeError(
        "No valid genes were found in "
        "TOP_50_PRIORITY_GENES.csv."
    )


# ============================================================
# ALIAS TABLE
# ============================================================

if ALIAS.exists():

    alias = pd.read_csv(
        ALIAS
    )

    print(
        f"Alias audit loaded: {len(alias)} rows"
    )

else:

    alias = pd.DataFrame()

    print(
        "WARNING: NCBI alias audit file was not found."
    )

    print(
        "Queries will use gene symbols only."
    )


# ============================================================
# ALIAS LOOKUP
# ============================================================

def aliases_for(gene):

    if alias.empty:

        return []

    if "gene_query" not in alias.columns:

        return []

    x = alias[
        alias[
            "gene_query"
        ]
        .astype(str)
        .str.strip()
        .str.lower()
        == gene.lower()
    ]

    if x.empty:

        return []

    raw = str(
        x.iloc[0].get(
            "aliases",
            ""
        )
    )

    if raw.lower() in {
        "",
        "nan",
        "none",
    }:

        return []

    return [
        a.strip()
        for a in
        raw
        .replace(
            "|",
            ";"
        )
        .split(";")
        if a.strip()
    ]


# ============================================================
# PUBMED QUERY
# ============================================================

def esearch_count(query):

    params = {
        "db":
            "pubmed",

        "term":
            query,

        "retmode":
            "json",

        "retmax":
            20,

        "tool":
            NCBI_TOOL,
    }

    if NCBI_EMAIL:

        params[
            "email"
        ] = NCBI_EMAIL

    if NCBI_API_KEY:

        params[
            "api_key"
        ] = NCBI_API_KEY

    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params,
        timeout=30,
    )

    response.raise_for_status()

    payload = response.json()

    result = payload[
        "esearchresult"
    ]

    count = int(
        result.get(
            "count",
            0
        )
    )

    ids = result.get(
        "idlist",
        []
    )

    return count, ids


# ============================================================
# AUDIT
# ============================================================

rows = []

total = len(
    genes
)

print("")
print("-" * 90)
print("BEGINNING PUBMED AUDIT")
print("-" * 90)
print("")


for index, gene_raw in enumerate(
    genes,
    start=1
):

    gene = gene_raw.strip().lower()

    aliases = aliases_for(
        gene
    )

    terms = [
        gene
    ] + aliases

    terms = sorted(
        set(
            t.strip()
            for t in terms
            if t.strip()
        )
    )

    name_clause = (
        " OR ".join(
            f'"{x}"[Title/Abstract]'
            for x in terms
        )
    )

    queries = {

        "direct_tbi_zebrafish":

            f'({name_clause}) '
            f'AND "Danio rerio"[Title/Abstract] '
            f'AND (brain OR telencephalon) '
            f'AND (injury OR trauma OR '
            f'"traumatic brain injury")',

        "regeneration_zebrafish":

            f'({name_clause}) '
            f'AND "Danio rerio"[Title/Abstract] '
            f'AND (regeneration OR regenerative)',

        "general_zebrafish":

            f'({name_clause}) '
            f'AND "Danio rerio"[Title/Abstract]',
    }

    rec = {

        "gene":
            gene,

        "aliases_used":
            ";".join(aliases),

        "n_aliases":
            len(aliases),
    }

    errors = []

    print(
        f"[{index:02d}/{total:02d}] "
        f"{gene}"
    )

    for label, query in queries.items():

        try:

            count, ids = esearch_count(
                query
            )

            rec[
                f"{label}_count"
            ] = count

            rec[
                f"{label}_pmids"
            ] = ";".join(
                ids[:10]
            )

            print(
                f"    {label}: "
                f"{count} hits"
            )

        except Exception as exc:

            rec[
                f"{label}_count"
            ] = None

            rec[
                f"{label}_pmids"
            ] = ""

            rec[
                f"{label}_error"
            ] = str(exc)

            errors.append(
                label
            )

            print(
                f"    {label}: "
                f"ERROR — {exc}"
            )

        # Respect configured NCBI delay.
        time.sleep(
            NCBI_DELAY_SECONDS
        )

    # ========================================================
    # CLASSIFY RESULT
    # ========================================================

    direct = rec.get(
        "direct_tbi_zebrafish_count"
    )

    regen = rec.get(
        "regeneration_zebrafish_count"
    )

    general = rec.get(
        "general_zebrafish_count"
    )

    # --------------------------------------------------------
    # API FAILURE TAKES PRIORITY OVER "NO EVIDENCE"
    # --------------------------------------------------------

    if errors:

        if len(errors) == 3:

            status = (
                "NCBI_QUERY_FAILED"
            )

        else:

            status = (
                "PARTIAL_NCBI_QUERY_FAILURE"
            )

    elif (
        isinstance(
            direct,
            (int, float)
        )
        and direct > 0
    ):

        status = (
            "direct_context_literature_found"
        )

    elif (
        isinstance(
            regen,
            (int, float)
        )
        and regen > 0
    ):

        status = (
            "related_regeneration_literature_found"
        )

    elif (
        isinstance(
            general,
            (int, float)
        )
        and general > 0
    ):

        status = (
            "gene_in_zebrafish_literature"
        )

    else:

        status = (
            "no_retrieved_context_evidence"
        )

    rec[
        "literature_context_status"
    ] = status

    rec[
        "query_error_count"
    ] = len(errors)

    rows.append(
        rec
    )


# ============================================================
# BUILD OUTPUT
# ============================================================

audit = pd.DataFrame(
    rows
)


# ============================================================
# SAVE
# ============================================================

audit.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# VERIFY OUTPUT
# ============================================================

if not OUTPUT_FILE.exists():

    raise RuntimeError(
        "Stage 21 completed internally but "
        "the expected output file was not created."
    )


if len(audit) != len(genes):

    raise RuntimeError(
        "Stage 21 output row count does not match "
        "the number of genes audited."
    )


# ============================================================
# SUMMARY
# ============================================================

status_counts = (
    audit[
        "literature_context_status"
    ]
    .value_counts()
)


print("")
print("=" * 90)
print("21 — BROAD LITERATURE CONTEXT AUDIT COMPLETE")
print("=" * 90)

print(
    f"Genes audited : {len(audit)}"
)

print(
    f"Output        : {OUTPUT_FILE}"
)

print("")
print("STATUS SUMMARY")
print("-" * 90)

for status, count in status_counts.items():

    print(
        f"{status}: {count}"
    )


# ============================================================
# QUERY FAILURE SUMMARY
# ============================================================

failed = audit[
    audit[
        "query_error_count"
    ]
    > 0
]

print("")
print(
    f"Genes with NCBI query errors: "
    f"{len(failed)}"
)


# ============================================================
# FINAL STATUS
# ============================================================

if len(failed) == 0:

    print("")
    print(
        "PASS — all TOP-50 genes completed "
        "all three PubMed queries."
    )

    print(
        "A zero-hit result is reported only as "
        "'no_retrieved_context_evidence'; "
        "it is NOT interpreted as proof of novelty."
    )

else:

    print("")
    print(
        "WARNING — one or more genes had "
        "NCBI query failures."
    )

    print(
        "Those records are explicitly marked and "
        "must not be interpreted as zero-literature evidence."
    )

print("")
print(
    "Stage 21 finished successfully."
)
