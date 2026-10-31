import time
import requests
import pandas as pd
from config import *

# ============================================================
# 23 — TARGETED LITERATURE VERIFICATION
# ============================================================
# Purpose:
#   Independently verify the TOP-20 prioritized genes using
#   multiple direct-context PubMed searches.
#
# Important:
#   - A retrieved paper does NOT prove biological relevance.
#   - A zero-hit result does NOT prove novelty.
#   - Search failures are recorded, never silently treated
#     as zero hits.
# ============================================================

INPUT = PRIORITIZATION_ROOT / "TOP_20_PRIORITY_GENES.csv"
OUT = LITERATURE_ROOT / "targeted_verification"
OUT.mkdir(parents=True, exist_ok=True)

PUBMED_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


# ============================================================
# PUBMED SEARCH
# ============================================================

def search_pubmed(query, retmax=10):

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
        PUBMED_URL,
        params=params,
        timeout=30
    )

    response.raise_for_status()

    data = response.json()

    result = data.get("esearchresult", {})

    count = int(result.get("count", 0))
    ids = result.get("idlist", [])

    return count, ids


# ============================================================
# LOAD TOP-20
# ============================================================

if not INPUT.exists():

    raise FileNotFoundError(
        f"TOP-20 prioritization file not found:\n{INPUT}"
    )

df = pd.read_csv(INPUT)

if "gene" not in df.columns:

    raise RuntimeError(
        "TOP_20_PRIORITY_GENES.csv does not contain a 'gene' column."
    )

genes = (
    df["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
    .drop_duplicates()
    .tolist()
)

print("=" * 90)
print("23 — TARGETED LITERATURE VERIFICATION")
print("=" * 90)

print(f"Input : {INPUT}")
print(f"Genes : {len(genes)}")
print()

if len(genes) == 0:
    raise RuntimeError("No genes found in TOP_20_PRIORITY_GENES.csv.")


# ============================================================
# SEARCH DEFINITIONS
# ============================================================

SEARCH_TYPES = {

    "brain_zebrafish":
        '"{gene}"[Title/Abstract] AND '
        '(zebrafish[Title/Abstract] OR "Danio rerio"[Title/Abstract]) AND '
        '(brain[Title/Abstract] OR telencephalon[Title/Abstract])',

    "brain_injury_zebrafish":
        '"{gene}"[Title/Abstract] AND '
        '(zebrafish[Title/Abstract] OR "Danio rerio"[Title/Abstract]) AND '
        '("brain injury"[Title/Abstract] OR '
        '"traumatic brain injury"[Title/Abstract] OR '
        'trauma[Title/Abstract])',

    "regeneration_zebrafish":
        '"{gene}"[Title/Abstract] AND '
        '(zebrafish[Title/Abstract] OR "Danio rerio"[Title/Abstract]) AND '
        '(regeneration[Title/Abstract] OR regenerative[Title/Abstract])',

    "tbi_general":
        '"{gene}"[Title/Abstract] AND '
        '("traumatic brain injury"[Title/Abstract] OR '
        '"brain injury"[Title/Abstract] OR '
        'trauma[Title/Abstract])'
}


# ============================================================
# RUN AUDIT
# ============================================================

rows = []

total_requests = len(genes) * len(SEARCH_TYPES)

request_number = 0


for gene_index, gene in enumerate(genes, start=1):

    print("-" * 90)
    print(f"[{gene_index:02d}/{len(genes)}] {gene}")

    record = {
        "gene": gene
    }

    all_pmids = set()

    for search_name, template in SEARCH_TYPES.items():

        request_number += 1

        query = template.format(gene=gene)

        try:

            count, ids = search_pubmed(
                query,
                retmax=10
            )

            record[f"{search_name}_count"] = count
            record[f"{search_name}_pmids"] = ";".join(ids[:10])

            all_pmids.update(ids)

            print(
                f"    {search_name}: "
                f"{count} hits"
            )

            record[f"{search_name}_status"] = "SUCCESS"

        except Exception as exc:

            record[f"{search_name}_count"] = None
            record[f"{search_name}_pmids"] = ""
            record[f"{search_name}_status"] = "ERROR"
            record[f"{search_name}_error"] = str(exc)

            print(
                f"    {search_name}: ERROR — {exc}"
            )

        # Small delay to avoid unnecessary API pressure.
        # Much faster than the previous version.
        time.sleep(0.25)

    # ========================================================
    # UNIQUE RETRIEVED PAPERS
    # ========================================================

    record["targeted_unique_pmids"] = len(all_pmids)

    record["targeted_pmids"] = ";".join(
        sorted(all_pmids)[:30]
    )

    # ========================================================
    # EVIDENCE FLAGS
    # ========================================================

    counts = []

    for search_name in SEARCH_TYPES:

        value = record.get(
            f"{search_name}_count"
        )

        if isinstance(value, (int, float)):

            counts.append(value)

    if any(x > 0 for x in counts):

        record["direct_paper_evidence"] = True

        record["literature_status"] = (
            "RETRIEVED_CONTEXT_EVIDENCE"
        )

    elif any(
        record.get(f"{s}_status") == "ERROR"
        for s in SEARCH_TYPES
    ):

        record["direct_paper_evidence"] = False

        record["literature_status"] = (
            "QUERY_ERROR_REQUIRES_REVIEW"
        )

    else:

        record["direct_paper_evidence"] = False

        record["literature_status"] = (
            "NO_RETRIEVED_CONTEXT_EVIDENCE"
        )

    rows.append(record)


# ============================================================
# SAVE RESULTS
# ============================================================

out = pd.DataFrame(rows)

output_file = (
    OUT / "TARGETED_LITERATURE_VERIFICATION.csv"
)

out.to_csv(
    output_file,
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

successful_genes = (
    out["literature_status"]
    .eq("RETRIEVED_CONTEXT_EVIDENCE")
    .sum()
)

zero_hit_genes = (
    out["literature_status"]
    .eq("NO_RETRIEVED_CONTEXT_EVIDENCE")
    .sum()
)

error_genes = (
    out["literature_status"]
    .eq("QUERY_ERROR_REQUIRES_REVIEW")
    .sum()
)


print()
print("=" * 90)
print("23 — TARGETED LITERATURE VERIFICATION COMPLETE")
print("=" * 90)

print(f"Genes audited          : {len(out)}")
print(f"Genes with evidence    : {successful_genes}")
print(f"Genes with zero hits   : {zero_hit_genes}")
print(f"Genes requiring review : {error_genes}")

print()
print("Output:")
print(output_file)

print()
print(
    "IMPORTANT:"
)
print(
    "Zero retrieved literature is NOT interpreted as biological novelty."
)

if error_genes == 0:

    print()
    print(
        "PASS — all targeted PubMed searches completed without errors."
    )

else:

    print()
    print(
        "REVIEW REQUIRED — one or more PubMed searches failed."
    )
