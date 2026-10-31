"""
18 — FUNCTIONAL ENRICHMENT

GSE278446 NeuroNovelty pipeline.

Reads only the authoritative Stage-17 enrichment inputs from
ENRICHMENT_ROOT (defined in config.py). Writes enrichment results
into ENRICHMENT_RESULTS_ROOT.
"""

from pathlib import Path
import json
import re
import sys
import time
import traceback

import pandas as pd
import requests
from config import (
    RESULT_ROOT,
    ENRICHMENT_ROOT,
    ENRICHMENT_RESULTS_ROOT,
    GP_API,
    GP_API_STABLE,
    GP_ORGANISM,
)

# ============================================================
# OUTPUT
# ============================================================

OUTPUT_ROOT = ENRICHMENT_RESULTS_ROOT
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

print("=" * 90)
print("18 — FUNCTIONAL ENRICHMENT")
print("=" * 90)
print(f"Result root   : {RESULT_ROOT}")
print(f"Input root    : {ENRICHMENT_ROOT}")
print(f"Output root   : {OUTPUT_ROOT}")

# ============================================================
# G:PROFILER
# ============================================================

ORGANISM = GP_ORGANISM

SOURCES = [
    "GO:BP",
    "GO:CC",
    "GO:MF",
    "KEGG",
    "REAC",
    "WP",
]

REQUEST_TIMEOUT = 120
P_VALUE_THRESHOLD = 0.05

EXPECTED_GENESETS = [
    "STRICT_3OF3_AVAILABLE",
    "TOP50",
    "REPRODUCIBLE_24H_PEAK",
    "REPRODUCIBLE_24H_DIP",
]

OUTPUT_COLUMNS = [
    "gene_set",
    "source",
    "term_id",
    "term_name",
    "description",
    "p_value",
    "intersection_size",
    "term_size",
    "query_size",
    "effective_domain_size",
    "precision",
    "recall",
    "significant",
    "intersection_genes",
    "status",
]

# ============================================================
# HELPERS
# ============================================================

def clean_text(value):
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').strip("'")

AUTHORITATIVE_INPUTS = {
    "STRICT_3OF3_AVAILABLE": ENRICHMENT_ROOT / "STRICT_3OF3_AVAILABLE.csv",
    "TOP50": ENRICHMENT_ROOT / "TOP50.csv",
    "REPRODUCIBLE_24H_PEAK": ENRICHMENT_ROOT / "REPRODUCIBLE_24H_PEAK.csv",
    "REPRODUCIBLE_24H_DIP": ENRICHMENT_ROOT / "REPRODUCIBLE_24H_DIP.csv",
}

def find_input_file(gene_set):
    path = AUTHORITATIVE_INPUTS.get(gene_set)
    if path is None or not path.exists():
        return None
    return path

def identify_gene_column(df):
    if df.empty:
        return None
    preferred = [
        "gene", "gene_id", "geneid", "gene_name", "gene_symbol",
        "symbol", "feature", "feature_id", "ensembl", "ensembl_id",
        "candidate",
    ]
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for name in preferred:
        if name in lower_map:
            return lower_map[name]
    best = None
    best_score = -1
    for column in df.columns:
        values = df[column].dropna().astype(str).str.strip()
        if values.empty:
            continue
        score = 0
        for value in values.head(500):
            if not value:
                continue
            if re.fullmatch(r"[-+]?\d+(\.\d+)?", value):
                continue
            if len(value) <= 100:
                score += 1
        if score > best_score:
            best_score = score
            best = column
    return best

def read_gene_list(path):
    suffix = path.suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep=None, engine="python", header=None)

    if df.empty:
        return []

    column = identify_gene_column(df)
    if column is None:
        return []

    genes = []
    for value in df[column]:
        value = clean_text(value)
        if not value:
            continue
        if value.lower() in {"gene", "gene_id", "gene_name", "symbol",
                             "nan", "none", "null"}:
            continue
        if ";" in value:
            pieces = value.split(";")
        elif "," in value and " " not in value:
            pieces = value.split(",")
        else:
            pieces = [value]
        for piece in pieces:
            piece = piece.strip()
            if piece:
                genes.append(piece)
    return list(dict.fromkeys(genes))

def find_background():
    """
    Prefer the CSV form of the background gene list. Stage 17 writes
    both BACKGROUND_GENES.csv and BACKGROUND_GENES.txt from the same
    set of symbols; the CSV is quoted by pandas and is not vulnerable
    to the delimiter-collision problem that can corrupt the .txt form
    when a gene symbol itself contains a comma.
    """
    csv_path = ENRICHMENT_ROOT / "BACKGROUND_GENES.csv"
    txt_path = ENRICHMENT_ROOT / "BACKGROUND_GENES.txt"

    candidate = csv_path if csv_path.exists() else txt_path

    if not candidate.exists():
        print(
            f"WARNING: No background file found in {ENRICHMENT_ROOT} "
            f"(looked for BACKGROUND_GENES.csv and .txt)"
        )
        return []

    try:
        genes = read_gene_list(candidate)
        print(f"Background genes found: {len(genes):,}  ({candidate.name})")
        return genes
    except Exception as exc:
        print(f"WARNING: Background read failed from {candidate}: {exc}")
        return []

def run_gprofiler(genes, background):
    payload = {
        "organism": ORGANISM,
        "query": genes,
        "sources": SOURCES,
        "user_threshold": P_VALUE_THRESHOLD,
        "all_results": True,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
        "output": "json",
    }
    if background:
        payload["domain_scope"] = "custom"
        payload["background"] = background

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "GSE278446-NeuroNovelty-Pipeline/1.0",
    }

    endpoints = [GP_API, GP_API_STABLE]
    last_error = None

    for endpoint in endpoints:
        for attempt in range(1, 4):
            try:
                print(f"g:Profiler POST {endpoint} attempt {attempt}/3 ...")
                response = requests.post(
                    endpoint,
                    json=payload,
                    headers=headers,
                    timeout=REQUEST_TIMEOUT,
                )
                if response.status_code == 200:
                    return response.json()
                last_error = RuntimeError(
                    f"HTTP {response.status_code} from {endpoint}: "
                    f"{response.text[:1000]}"
                )
                print(f"WARNING: {last_error}")
            except Exception as exc:
                last_error = exc
                print(f"WARNING: {endpoint} attempt {attempt} failed: {exc}")
            if attempt < 3:
                time.sleep(attempt * 3)

    raise RuntimeError(f"g:Profiler failed on all endpoints: {last_error}")

def make_status_row(gene_set, status):
    return pd.DataFrame(
        [{
            "gene_set": gene_set,
            "source": "",
            "term_id": "",
            "term_name": "",
            "description": "",
            "p_value": "",
            "intersection_size": 0,
            "term_size": 0,
            "query_size": 0,
            "effective_domain_size": 0,
            "precision": 0,
            "recall": 0,
            "significant": False,
            "intersection_genes": "",
            "status": status,
        }],
        columns=OUTPUT_COLUMNS,
    )

def flatten_intersections(intersections):
    values = []
    for item in intersections or []:
        if isinstance(item, (list, tuple)):
            values.extend(str(x) for x in item if x is not None)
        elif item is not None:
            values.append(str(item))
    return values

def parse_result(data, gene_set):
    result = data.get("result", [])
    if not result:
        return make_status_row(gene_set, "NO_ENRICHMENT_RESULTS")

    rows = []
    for item in result:
        intersections = item.get("intersections", []) or []
        if not isinstance(intersections, list):
            intersections = [intersections]
        rows.append({
            "gene_set": gene_set,
            "source": item.get("source", ""),
            "term_id": item.get("native", ""),
            "term_name": item.get("name", ""),
            "description": item.get("description", ""),
            "p_value": item.get("p_value", ""),
            "intersection_size": item.get("intersection_size", 0),
            "term_size": item.get("term_size", 0),
            "query_size": item.get("query_size", 0),
            "effective_domain_size": item.get("effective_domain_size", 0),
            "precision": item.get("precision", 0),
            "recall": item.get("recall", 0),
            "significant": item.get("significant", False),
            "intersection_genes": ";".join(
                flatten_intersections(intersections)
            ),
            "status": "SUCCESS",
        })
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

def fetch_gprofiler_metadata():
    metadata = {
        "organism": ORGANISM,
        "primary_endpoint": GP_API,
        "fallback_endpoint": GP_API_STABLE,
        "retrieved_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
    }
    try:
        response = requests.get(
            "https://biit.cs.ut.ee/gprofiler/api/util/data_versions/",
            params={"organism": ORGANISM},
            timeout=30,
        )
        response.raise_for_status()
        metadata["data_versions_response"] = response.json()
    except Exception as exc:
        metadata["data_versions_error"] = str(exc)
    return metadata

# ============================================================
# MAIN
# ============================================================

def main():

    background = find_background()

    metadata = fetch_gprofiler_metadata()
    with open(
        OUTPUT_ROOT / "GPROFILER_RUN_METADATA.json",
        "w",
        encoding="utf-8",
    ) as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)

    print("")

    summary = []

    for gene_set in EXPECTED_GENESETS:

        print("-" * 80)
        print(f"Running {gene_set}")
        print("-" * 80)

        input_file = find_input_file(gene_set)

        if input_file is None:
            print(f"WARNING: No input file found for {gene_set}")
            output_file = OUTPUT_ROOT / f"{gene_set}_enrichment.csv"
            df = make_status_row(gene_set, "INPUT_NOT_FOUND")
            df.to_csv(output_file, index=False)
            summary.append({
                "gene_set": gene_set,
                "input_file": "",
                "genes": 0,
                "status": "INPUT_NOT_FOUND",
                "rows": 0,
            })
            continue

        print(f"Input file: {input_file}")

        try:
            genes = read_gene_list(input_file)
        except Exception as exc:
            print(f"WARNING: Failed to read input: {exc}")
            genes = []

        print(f"Genes supplied: {len(genes):,}")

        output_file = OUTPUT_ROOT / f"{gene_set}_enrichment.csv"

        if not genes:
            df = make_status_row(gene_set, "EMPTY_GENE_LIST")
            df.to_csv(output_file, index=False)
            summary.append({
                "gene_set": gene_set,
                "input_file": str(input_file),
                "genes": 0,
                "status": "EMPTY_GENE_LIST",
                "rows": 0,
            })
            continue

        try:
            data = run_gprofiler(genes, background)

            raw_file = OUTPUT_ROOT / f"{gene_set}_gprofiler_raw.json"
            with open(raw_file, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2, ensure_ascii=False)

            df = parse_result(data, gene_set)
            df.to_csv(output_file, index=False)

            successful_rows = int((df["status"] == "SUCCESS").sum())

            print(f"Enrichment rows: {successful_rows:,}")

            summary.append({
                "gene_set": gene_set,
                "input_file": str(input_file),
                "genes": len(genes),
                "status": "SUCCESS",
                "rows": successful_rows,
            })

        except Exception as exc:
            print(f"WARNING: {gene_set} enrichment failed:")
            print(exc)

            df = make_status_row(gene_set, "GPROFILER_FAILED")
            df.to_csv(output_file, index=False)

            summary.append({
                "gene_set": gene_set,
                "input_file": str(input_file),
                "genes": len(genes),
                "status": "GPROFILER_FAILED",
                "rows": 0,
            })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(
        OUTPUT_ROOT / "ENRICHMENT_RUN_SUMMARY.csv",
        index=False,
    )

    print("")
    print("=" * 90)
    print("18 — FUNCTIONAL ENRICHMENT COMPLETE")
    print("=" * 90)
    print(
        f"Summary: {OUTPUT_ROOT / 'ENRICHMENT_RUN_SUMMARY.csv'}"
    )
    print("")
    print(summary_df.to_string(index=False))

if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
