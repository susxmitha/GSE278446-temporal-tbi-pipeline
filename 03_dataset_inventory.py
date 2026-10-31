import gzip
import pandas as pd
from config import *

# ============================================================
# 03 — COMPLETE DATASET INVENTORY
# ============================================================
# 28 samples are audited at the source level.
# 27 samples are included in downstream biological analysis.
#
# GSM8546673 is retained in the source audit but excluded from
# biological analysis because it is an exact source-level
# duplicate of GSM8546664.
# ============================================================


def locate_files(gsm):
    result = {}

    for p in RAW_ROOT.rglob("*"):
        if not p.is_file():
            continue

        if gsm.lower() not in p.name.lower():
            continue

        n = p.name.lower()

        if "matrix" in n and n.endswith((".mtx", ".mtx.gz")):
            result["matrix"] = p

        elif "barcodes" in n and n.endswith(
            (".tsv.gz", ".tsv", ".txt.gz", ".txt")
        ):
            result["barcodes"] = p

        elif "features" in n and n.endswith(
            (".tsv.gz", ".tsv", ".txt.gz", ".txt")
        ):
            result["features"] = p

    return result


def count_lines(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open

    with opener(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:
        return sum(1 for _ in fh)


def matrix_shape_from_header(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open

    with opener(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        line = fh.readline()

        while line.startswith("%"):
            line = fh.readline()

        parts = line.split()

        if len(parts) < 3:
            raise ValueError(
                f"Invalid Matrix Market header in {path}"
            )

        return (
            int(parts[0]),
            int(parts[1]),
            int(parts[2])
        )


# ============================================================
# SOURCE-LEVEL INVENTORY
# ============================================================

SOURCE_GSMS = tuple(SAMPLE_METADATA.keys())

ANALYSIS_GSMS = tuple(
    gsm
    for gsm in SOURCE_GSMS
    if gsm not in EXCLUDED_GSMS
)


rows = []

for gsm in SOURCE_GSMS:

    meta = SAMPLE_METADATA[gsm]

    f = locate_files(gsm)

    missing = [
        k
        for k in ("matrix", "barcodes", "features")
        if k not in f
    ]

    included = gsm in ANALYSIS_GSMS

    exclusion_reason = ""

    if gsm in EXCLUDED_GSMS:
        exclusion_reason = (
            "Exact source-level duplicate of GSM8546664"
        )

    if missing:

        rows.append({
            "GSM": gsm,
            **meta,
            "included_in_analysis": included,
            "exclusion_reason": exclusion_reason,
            "matrix_rows": "",
            "matrix_cols": "",
            "feature_rows": "",
            "barcode_rows": "",
            "nonzero_entries": "",
            "matrix_file": "",
            "features_file": "",
            "barcodes_file": "",
            "status": "MISSING",
        })

        continue

    feature_count = count_lines(f["features"])
    barcode_count = count_lines(f["barcodes"])

    matrix_rows, matrix_cols, nnz = (
        matrix_shape_from_header(f["matrix"])
    )

    status = "PASS"

    if feature_count != matrix_rows:
        status = "FEATURE_COUNT_MISMATCH"

    if barcode_count != matrix_cols:
        status = "BARCODE_COUNT_MISMATCH"

    rows.append({
        "GSM": gsm,
        **meta,
        "included_in_analysis": included,
        "exclusion_reason": exclusion_reason,
        "matrix_rows": matrix_rows,
        "matrix_cols": matrix_cols,
        "feature_rows": feature_count,
        "barcode_rows": barcode_count,
        "nonzero_entries": nnz,
        "matrix_file": str(f["matrix"]),
        "features_file": str(f["features"]),
        "barcodes_file": str(f["barcodes"]),
        "status": status,
    })


inventory = pd.DataFrame(rows)

out = RESULT_ROOT / "03_dataset_inventory.csv"

inventory.to_csv(
    out,
    index=False
)


# ============================================================
# REPORT
# ============================================================

print("=" * 78)
print("03 — COMPLETE DATASET INVENTORY")
print("=" * 78)

print()
print(inventory.to_string(index=False))

print()
print("-" * 78)

if not (inventory["status"] == "PASS").all():

    raise SystemExit(
        "DATASET INVENTORY FAILED — "
        "fix missing/mismatched raw files before continuing."
    )


source_count = len(SOURCE_GSMS)
analysis_count = len(ANALYSIS_GSMS)
excluded_count = len(EXCLUDED_GSMS)

print(
    f"PASS — {source_count} original GEO samples "
    f"are complete and dimensionally consistent."
)

print(
    f"Analysis samples: {analysis_count}"
)

print(
    f"Excluded source duplicates: {excluded_count}"
)

print(
    f"Output: {out}"
)
