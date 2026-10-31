from pathlib import Path
import gzip
import re
import pandas as pd
from config import *

# ============================================================
# 00 — RAW FILE SMOKE TEST
# ============================================================
# Purpose:
#   1. Verify ALL original GEO GSM records are physically present.
#   2. Verify each GSM has exactly one matrix, barcode and feature file.
#   3. Preserve the duplicate GSM in the source audit.
#   4. Mark GSM8546673 as excluded from biological analysis.
#
# IMPORTANT:
# SOURCE_GSMS  = original GEO source inventory (28)
# ANALYSIS_GSMS = downstream biological analysis set (27)
# ============================================================

def find_gsm_files():
    found = {gsm: {} for gsm in SOURCE_GSMS}

    for path in RAW_ROOT.rglob("*"):
        if not path.is_file():
            continue

        m = re.search(r"(GSM\d+)", path.name, flags=re.I)
        if not m:
            continue

        gsm = m.group(1).upper()

        if gsm not in found:
            continue

        name = path.name.lower()

        if "matrix" in name and name.endswith((".mtx", ".mtx.gz")):
            found[gsm].setdefault("matrix", []).append(path)

        elif "barcodes" in name and name.endswith(
            (".tsv", ".tsv.gz", ".txt", ".txt.gz")
        ):
            found[gsm].setdefault("barcodes", []).append(path)

        elif "features" in name and name.endswith(
            (".tsv", ".tsv.gz", ".txt", ".txt.gz")
        ):
            found[gsm].setdefault("features", []).append(path)

    return found


if not RAW_ROOT.exists():
    raise FileNotFoundError(
        f"RAW_ROOT does not exist:\n{RAW_ROOT}"
    )


found = find_gsm_files()

rows = []
errors = []

for gsm in SOURCE_GSMS:

    rec = {
        "GSM": gsm,
        "source_status": (
            "EXCLUDED_SOURCE_DUPLICATE"
            if gsm in EXCLUDED_GSMS
            else "ANALYSIS_INCLUDED"
        ),
    }

    for kind in ("matrix", "barcodes", "features"):

        paths = found[gsm].get(kind, [])

        rec[f"{kind}_count"] = len(paths)
        rec[f"{kind}_path"] = (
            str(paths[0])
            if len(paths) == 1
            else ""
        )

        if len(paths) != 1:
            errors.append(
                f"{gsm}: expected exactly 1 {kind} file, "
                f"found {len(paths)}"
            )
            continue

        path = paths[0]

        try:
            if path.name.lower().endswith(".gz"):
                with gzip.open(path, "rb") as fh:
                    chunk = fh.read(1024)
            else:
                with open(path, "rb") as fh:
                    chunk = fh.read(1024)

            if not chunk:
                raise ValueError("empty file")

            rec[f"{kind}_readable"] = True

        except Exception as exc:

            rec[f"{kind}_readable"] = False

            errors.append(
                f"{gsm}: {kind} unreadable: {exc}"
            )

    rows.append(rec)


manifest = pd.DataFrame(rows)

out = RESULT_ROOT / "00_raw_file_smoke_test.csv"
manifest.to_csv(out, index=False)


print("=" * 78)
print("00 — RAW FILE SMOKE TEST")
print("=" * 78)

print(f"Original GEO GSMs : {len(SOURCE_GSMS)}")
print(f"Analysis GSMs     : {len(ANALYSIS_GSMS)}")
print(f"Excluded GSMs     : {len(EXCLUDED_GSMS)}")
print(f"Output            : {out}")
print(f"Errors            : {len(errors)}")

if EXCLUDED_GSMS:
    print("\nExcluded source duplicate(s):")
    for gsm in sorted(EXCLUDED_GSMS):
        print(f"  - {gsm}")

if errors:

    print("\nERRORS:")
    for e in errors:
        print(f"  - {e}")

    raise SystemExit(
        "RAW FILE SMOKE TEST FAILED"
    )

print()
print(
    "PASS — all original GEO GSMs have one readable "
    "matrix, barcode and feature file."
)

print(
    "The source-level duplicate is retained for audit "
    "but excluded from downstream biological analysis."
)
