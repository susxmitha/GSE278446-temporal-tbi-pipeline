from config import *
import gzip
from pathlib import Path
import pandas as pd
import scipy.io
import scipy.sparse as sp
from config import *

# ============================================================
# 01 — SINGLE MATRIX LOAD / DIMENSION TEST
# ============================================================
# Loads one representative GSM and verifies:
# matrix rows == feature rows; matrix columns == barcode rows.
# No filtering or biological scoring is performed.
# ============================================================

TEST_GSM = "GSM8546663"

def locate_files(gsm):
    result = {}
    for p in RAW_ROOT.rglob("*"):
        if not p.is_file() or gsm.lower() not in str(p).lower():
            continue
        n = p.name.lower()
        if "matrix" in n and n.endswith((".mtx", ".mtx.gz")):
            result["matrix"] = p
        elif "barcodes" in n and n.endswith((".tsv.gz", ".tsv", ".txt.gz", ".txt")):
            result["barcodes"] = p
        elif "features" in n and n.endswith((".tsv.gz", ".tsv", ".txt.gz", ".txt")):
            result["features"] = p
    missing = [k for k in ("matrix", "barcodes", "features") if k not in result]
    if missing:
        raise FileNotFoundError(f"{gsm}: missing {missing}")
    return result

def read_lines(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        return [x.rstrip("\n\r") for x in fh]

f = locate_files(TEST_GSM)
barcodes = read_lines(f["barcodes"])
features = pd.read_csv(
    f["features"], sep="\t", header=None,
    compression="gzip" if str(f["features"]).lower().endswith(".gz") else None
)
matrix = scipy.io.mmread(f["matrix"]).tocsr()

if matrix.shape[0] != len(features):
    raise RuntimeError(
        f"Feature/matrix mismatch: matrix rows={matrix.shape[0]}, features={len(features)}"
    )
if matrix.shape[1] != len(barcodes):
    raise RuntimeError(
        f"Barcode/matrix mismatch: matrix cols={matrix.shape[1]}, barcodes={len(barcodes)}"
    )
if matrix.nnz == 0:
    raise RuntimeError("Expression matrix contains zero non-zero entries.")

# Cell Ranger feature files normally contain gene_id, gene_name, feature_type.
# Verify the feature type column if present.
feature_type_summary = {}
if features.shape[1] >= 3:
    feature_type_summary = features.iloc[:, 2].astype(str).value_counts().to_dict()

report = pd.DataFrame([{
    "GSM": TEST_GSM,
    "matrix_rows_genes": matrix.shape[0],
    "matrix_cols_cells": matrix.shape[1],
    "features_rows": len(features),
    "barcodes_rows": len(barcodes),
    "nonzero_entries": matrix.nnz,
    "matrix_orientation_verified": True,
    "feature_types": "; ".join(f"{k}={v}" for k, v in feature_type_summary.items()),
}])

out = RESULT_ROOT / "01_single_matrix_load_test.csv"
report.to_csv(out, index=False)

print("=" * 78)
print("01 — SINGLE MATRIX LOAD TEST")
print("=" * 78)
print(report.to_string(index=False))
print(f"\nPASS — {TEST_GSM} matrix, features and barcodes are dimensionally consistent.")
print(f"Output: {out}")
