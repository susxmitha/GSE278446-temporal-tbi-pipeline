import time
import requests
import xml.etree.ElementTree as ET
import pandas as pd
from config import *

# ============================================================
# 24 — PAPER-LEVEL VERIFICATION OF TOP CANDIDATES
# ============================================================
# Purpose:
#   Retrieve the actual PubMed records associated with Stage 23
#   and inspect title/abstract text for direct biological-context
#   evidence.
#
# IMPORTANT:
#   This stage does NOT claim novelty.
#   Zero evidence means only that no qualifying evidence was
#   retrieved from the reviewed PubMed records.
#
# Optimization:
#   Instead of making one EFetch request per gene, all unique
#   PMIDs are retrieved in batches.
# ============================================================

INPUT = (
    LITERATURE_ROOT
    / "targeted_verification"
    / "TARGETED_LITERATURE_VERIFICATION.csv"
)

OUT = LITERATURE_ROOT / "paper_level"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 90)
print("24 — PAPER-LEVEL VERIFICATION OF TOP CANDIDATES")
print("=" * 90)

# ------------------------------------------------------------
# LOAD STAGE 23
# ------------------------------------------------------------

if not INPUT.exists():
    raise FileNotFoundError(
        f"Stage 23 output not found:\n{INPUT}"
    )

cand = pd.read_csv(INPUT)

if "gene" not in cand.columns:
    raise RuntimeError(
        "Stage 23 output does not contain the required 'gene' column."
    )

if "targeted_pmids" not in cand.columns:
    raise RuntimeError(
        "Stage 23 output does not contain the required "
        "'targeted_pmids' column."
    )

cand["gene"] = (
    cand["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
)

print(f"Genes loaded: {len(cand)}")

# ------------------------------------------------------------
# COLLECT UNIQUE PUBMED IDS
# ------------------------------------------------------------

gene_to_pmids = {}

all_pmids = []

for _, row in cand.iterrows():

    gene = row["gene"]

    raw = str(row.get("targeted_pmids", ""))

    if raw.lower() in {"nan", "none", ""}:
        ids = []
    else:
        ids = [
            x.strip()
            for x in raw.split(";")
            if x.strip()
        ]

    # Remove duplicates while preserving order
    ids = list(dict.fromkeys(ids))

    gene_to_pmids[gene] = ids

    all_pmids.extend(ids)

all_pmids = list(dict.fromkeys(all_pmids))

print(f"Unique PubMed records to retrieve: {len(all_pmids)}")

if not all_pmids:

    print()
    print("WARNING: No PubMed IDs were supplied by Stage 23.")
    print("No paper-level evidence can be evaluated.")

    empty_papers = pd.DataFrame(
        columns=[
            "gene",
            "pmid",
            "title",
            "abstract",
            "journal",
            "year",
            "direct_context_match",
            "context_reason"
        ]
    )

    empty_summary = pd.DataFrame(
        columns=[
            "gene",
            "papers_reviewed",
            "direct_context_matches",
            "direct_context_evidence",
            "literature_evidence_status"
        ]
    )

    empty_papers.to_csv(
        OUT / "TOP20_PAPER_LEVEL_EVIDENCE.csv",
        index=False
    )

    empty_summary.to_csv(
        OUT / "TOP20_PAPER_LEVEL_SUMMARY.csv",
        index=False
    )

    print()
    print("24 — PAPER-LEVEL VERIFICATION COMPLETE")
    print(f"Output: {OUT}")
    raise SystemExit(0)


# ------------------------------------------------------------
# FETCH PUBMED RECORDS IN BATCHES
# ------------------------------------------------------------

def fetch_batch(ids):

    if not ids:
        return []

    params = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "xml",
        "tool": NCBI_TOOL
    }

    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL

    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params=params,
        timeout=90
    )

    response.raise_for_status()

    root = ET.fromstring(response.text)

    records = []

    for article in root.findall(".//PubmedArticle"):

        pmid_node = article.find(".//PMID")

        if pmid_node is not None:
            pmid = "".join(pmid_node.itertext()).strip()
        else:
            pmid = ""

        title_node = article.find(".//ArticleTitle")

        if title_node is not None:
            title = " ".join(
                "".join(title_node.itertext()).split()
            )
        else:
            title = ""

        abstract_parts = []

        for node in article.findall(
            ".//Abstract/AbstractText"
        ):
            abstract_parts.append(
                " ".join(
                    "".join(node.itertext()).split()
                )
            )

        abstract = " ".join(abstract_parts)

        journal = article.findtext(
            ".//Journal/Title",
            default=""
        )

        year = ""

        # Prefer explicit publication year
        year_node = article.find(
            ".//PubDate/Year"
        )

        if year_node is not None:
            year = "".join(
                year_node.itertext()
            ).strip()

        # Fallback for electronic publication dates
        if not year:
            medline_date = article.find(
                ".//PubDate/MedlineDate"
            )

            if medline_date is not None:
                year = "".join(
                    medline_date.itertext()
                ).strip()[:4]

        records.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year
        })

    return records


# ------------------------------------------------------------
# RETRIEVE ALL RECORDS
# ------------------------------------------------------------

BATCH_SIZE = 100

paper_records = {}

batches = [
    all_pmids[i:i + BATCH_SIZE]
    for i in range(
        0,
        len(all_pmids),
        BATCH_SIZE
    )
]

print()
print(
    f"Retrieving PubMed records in {len(batches)} batch(es)..."
)

for i, batch in enumerate(batches, start=1):

    print(
        f"  Batch {i}/{len(batches)} "
        f"({len(batch)} PMIDs)"
    )

    try:

        records = fetch_batch(batch)

        for record in records:

            pmid = record["pmid"]

            if pmid:
                paper_records[pmid] = record

    except Exception as exc:

        print(
            f"  WARNING: batch {i} failed: {exc}"
        )

    # Only a short pause between batches
    if i < len(batches):
        time.sleep(0.25)


print(
    f"PubMed records successfully retrieved: "
    f"{len(paper_records)}"
)


# ------------------------------------------------------------
# PAPER-LEVEL CONTEXT CLASSIFICATION
# ------------------------------------------------------------

def classify_context(title, abstract):

    text = (
        str(title) +
        " " +
        str(abstract)
    ).lower()

    has_zebrafish = (
        "zebrafish" in text
        or "danio rerio" in text
    )

    has_brain = (
        "brain" in text
        or "telencephalon" in text
    )

    has_injury = (
        "injury" in text
        or "trauma" in text
        or "traumatic brain injury" in text
    )

    has_regeneration = (
        "regeneration" in text
        or "regenerative" in text
    )

    # Strongest context:
    # zebrafish + brain + injury/regeneration
    if (
        has_zebrafish
        and has_brain
        and (
            has_injury
            or has_regeneration
        )
    ):
        return True, "direct_zebrafish_brain_injury_or_regeneration_context"

    if (
        has_zebrafish
        and has_brain
    ):
        return False, "zebrafish_brain_context_without_injury_or_regeneration"

    if (
        has_zebrafish
        and (
            has_injury
            or has_regeneration
        )
    ):
        return False, "zebrafish_injury_or_regeneration_without_brain_context"

    return False, "no_direct_target_context_detected"


# ------------------------------------------------------------
# BUILD PAPER-LEVEL TABLE
# ------------------------------------------------------------

rows = []

for gene, ids in gene_to_pmids.items():

    for pmid in ids:

        record = paper_records.get(pmid)

        if record is None:

            rows.append({
                "gene": gene,
                "pmid": pmid,
                "title": "",
                "abstract": "",
                "journal": "",
                "year": "",
                "direct_context_match": False,
                "context_reason": "pubmed_record_not_retrieved"
            })

            continue

        direct, reason = classify_context(
            record["title"],
            record["abstract"]
        )

        rows.append({
            "gene": gene,
            "pmid": record["pmid"],
            "title": record["title"],
            "abstract": record["abstract"],
            "journal": record["journal"],
            "year": record["year"],
            "direct_context_match": direct,
            "context_reason": reason
        })


papers = pd.DataFrame(rows)

# ------------------------------------------------------------
# SAVE PAPER-LEVEL EVIDENCE
# ------------------------------------------------------------

paper_output = (
    OUT /
    "TOP20_PAPER_LEVEL_EVIDENCE.csv"
)

papers.to_csv(
    paper_output,
    index=False
)


# ------------------------------------------------------------
# CREATE GENE-LEVEL SUMMARY
# ------------------------------------------------------------

summary = []

for gene in cand["gene"].drop_duplicates():

    g = papers[
        papers["gene"] == gene
    ]

    reviewed = len(g)

    direct_matches = int(
        g["direct_context_match"]
        .fillna(False)
        .sum()
    )

    if direct_matches > 0:

        status = (
            "direct_context_literature_found"
        )

    elif reviewed > 0:

        status = (
            "no_direct_context_evidence_in_reviewed_papers"
        )

    else:

        status = (
            "no_paper_records_available"
        )

    summary.append({
        "gene": gene,
        "papers_reviewed": reviewed,
        "direct_context_matches": direct_matches,
        "direct_context_evidence": (
            direct_matches > 0
        ),
        "literature_evidence_status": status
    })


summary_df = pd.DataFrame(summary)

summary_output = (
    OUT /
    "TOP20_PAPER_LEVEL_SUMMARY.csv"
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
print("24 — PAPER-LEVEL VERIFICATION COMPLETE")
print("=" * 90)

print(
    f"Genes evaluated              : "
    f"{len(summary_df)}"
)

print(
    f"Unique PubMed records found : "
    f"{len(paper_records)}"
)

print(
    f"Paper-level rows generated  : "
    f"{len(papers)}"
)

print(
    f"Genes with direct context   : "
    f"{int(summary_df['direct_context_evidence'].sum())}"
)

print()
print("Outputs:")
print(paper_output)
print(summary_output)

print()
print(
    "IMPORTANT:"
)
print(
    "Paper-level evidence is based on the retrieved "
    "PubMed title/abstract text."
)
print(
    "Absence of a direct-context match is NOT interpreted "
    "as proof of biological novelty."
)
