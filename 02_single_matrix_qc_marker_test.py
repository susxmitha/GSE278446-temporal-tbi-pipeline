import gzip
import numpy as np
import pandas as pd
import scipy.io
from config import *

# ============================================================
# 02 — SINGLE MATRIX QC + MARKER AVAILABILITY TEST
# ============================================================
# Representative-sample QC only. It does not decide final QC
# thresholds; those are applied consistently in script 06.
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
    return result

f = locate_files(TEST_GSM)
if len(f) != 3:
    raise FileNotFoundError(f"Could not locate all three files for {TEST_GSM}: {f}")

features = pd.read_csv(
    f["features"], sep="\t", header=None,
    compression="gzip" if str(f["features"]).lower().endswith(".gz") else None
)
matrix = scipy.io.mmread(f["matrix"]).tocsr()

if features.shape[1] >= 2:
    gene_names = features.iloc[:, 1].astype(str).str.strip().str.lower().tolist()
else:
    gene_names = features.iloc[:, 0].astype(str).str.strip().str.lower().tolist()

total_counts = np.asarray(matrix.sum(axis=0)).ravel()
genes_detected = np.asarray((matrix > 0).sum(axis=0)).ravel()

mt_mask = np.array([g.upper().startswith("MT-") or g.startswith("mt:") for g in gene_names])
if mt_mask.sum():
    mito_counts = np.asarray(matrix[mt_mask, :].sum(axis=0)).ravel()
    pct_mito = 100 * mito_counts / np.maximum(total_counts, 1)
else:
    pct_mito = np.zeros_like(total_counts, dtype=float)

qc = pd.DataFrame({
    "GSM": TEST_GSM,
    "total_counts": total_counts,
    "genes_detected": genes_detected,
    "pct_mito": pct_mito,
})
summary = qc[["total_counts", "genes_detected", "pct_mito"]].describe(
    percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
).T
summary.to_csv(RESULT_ROOT / "02_representative_qc_summary.csv")

marker_rows = []
gene_set = set(gene_names)
for program, markers in MARKER_PROGRAMS.items():
    found = sorted(set(m.lower() for m in markers) & gene_set)
    marker_rows.append({
        "program": program,
        "markers_requested": len(set(m.lower() for m in markers)),
        "markers_found": len(found),
        "marker_fraction": len(found) / max(len(set(m.lower() for m in markers)), 1),
        "found_genes": ";".join(found),
    })
marker_df = pd.DataFrame(marker_rows)
marker_df.to_csv(RESULT_ROOT / "02_marker_availability.csv", index=False)

print("=" * 78)
print("02 — QC + MARKER TEST")
print("=" * 78)
print(summary)
print("\nMarker availability:")
print(marker_df.to_string(index=False))
print("\nPASS — representative QC metrics and marker availability were calculated.")
