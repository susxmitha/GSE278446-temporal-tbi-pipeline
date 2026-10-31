import time
import requests
import pandas as pd
from config import *

# ============================================================
# 20 — NCBI GENE SYMBOL / ALIAS AUDIT
# ============================================================
# This script verifies naming. It does NOT invent ortholog mappings.
# Orthology is reported only if the external record explicitly
# exposes a corresponding field; otherwise it remains unresolved.
#
# Every row emits the same columns, so downstream stages never see
# NaN-filled schema drift when one gene fails.
# ============================================================

INPUT = PRIORITIZATION_ROOT / "TOP_50_PRIORITY_GENES.csv"
OUT = ALIAS_ROOT
OUT.mkdir(parents=True, exist_ok=True)

if not INPUT.exists():
    raise FileNotFoundError(
        f"Stage 16 prioritization output not found:\n{INPUT}"
    )

df = pd.read_csv(INPUT)
if "gene" not in df.columns:
    raise RuntimeError(
        f"TOP_50_PRIORITY_GENES.csv lacks 'gene' column. "
        f"Columns present: {list(df.columns)}"
    )

genes = (
    df["gene"]
    .astype(str)
    .str.strip()
    .str.lower()
    .drop_duplicates()
    .tolist()
)

# ------------------------------------------------------------
# COLUMN SCHEMA — every row uses exactly these keys.
# ------------------------------------------------------------

COLUMNS = [
    "gene_query",
    "ncbi_gene_id",
    "official_symbol",
    "description",
    "organism",
    "aliases",
    "ncbi_resolved",
    "ortholog_status",
    "error",
]


def empty_row(gene, error=""):
    return {
        "gene_query": gene,
        "ncbi_gene_id": "",
        "official_symbol": "",
        "description": "",
        "organism": "",
        "aliases": "",
        "ncbi_resolved": False,
        "ortholog_status": "not_assessed_without_explicit_mapping",
        "error": error,
    }


# ------------------------------------------------------------
# NCBI HELPERS
# ------------------------------------------------------------

def _add_ncbi_auth(params):
    if NCBI_EMAIL and NCBI_EMAIL.strip():
        params["email"] = NCBI_EMAIL.strip()
    if NCBI_API_KEY and NCBI_API_KEY.strip():
        params["api_key"] = NCBI_API_KEY.strip()
    return params


def ncbi_esearch(term):
    params = {
        "db": "gene",
        "term": f'("{term}"[All Fields]) AND "Danio rerio"[Organism]',
        "retmode": "json",
        "retmax": 5,
        "tool": NCBI_TOOL,
    }
    _add_ncbi_auth(params)
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params,
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def ncbi_summary(ids):
    if not ids:
        return {}
    params = {
        "db": "gene",
        "id": ",".join(ids),
        "retmode": "json",
        "tool": NCBI_TOOL,
    }
    _add_ncbi_auth(params)
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params=params,
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("result", {})


# ------------------------------------------------------------
# MAIN LOOP
# ------------------------------------------------------------

rows = []

for i, gene in enumerate(genes, start=1):
    print(f"[{i:02d}/{len(genes)}] {gene}")

    try:
        ids = ncbi_esearch(gene)
        summ = ncbi_summary(ids)
        first = next(
            (summ[k] for k in ids if k in summ),
            {},
        )

        row = empty_row(gene, error="")
        row.update({
            "ncbi_gene_id": first.get("uid", ""),
            "official_symbol": first.get("name", ""),
            "description": first.get("description", ""),
            "organism": first.get("organism", ""),
            "aliases": first.get("otheraliases", ""),
            "ncbi_resolved": bool(ids),
        })
        rows.append(row)

        if not ids:
            print(f"    no NCBI match")
        else:
            print(f"    id={row['ncbi_gene_id']}  symbol={row['official_symbol']}")

    except Exception as exc:
        rows.append(empty_row(gene, error=str(exc)))
        print(f"    ERROR: {exc}")

    time.sleep(NCBI_DELAY_SECONDS)

out = pd.DataFrame(rows, columns=COLUMNS)
out.to_csv(OUT / "NCBI_GENE_ALIAS_AUDIT.csv", index=False)

print("=" * 78)
print("20 — NCBI ALIAS AUDIT COMPLETE")
print(f"Genes audited : {len(out)}")
print(f"Resolved      : {int(out['ncbi_resolved'].sum())}")
print(f"Unresolved    : {int((~out['ncbi_resolved']).sum())}")
print(f"Output        : {OUT}")
