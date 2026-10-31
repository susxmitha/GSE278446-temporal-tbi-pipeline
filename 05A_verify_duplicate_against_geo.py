from pathlib import Path
import hashlib
import gzip
import pandas as pd
from config import (
    RESEARCH_ROOT,
    RAW_ROOT,
    RESULT_ROOT,
    SOURCE_DUPLICATE_KEEP,
    SOURCE_DUPLICATE_PAIRS,
    EXCLUDED_GSMS,
)

# ============================================================
# 05A — ORIGINAL GEO SOURCE DUPLICATE VERIFICATION
# ============================================================
# Purpose:
#   Independently verify the known source-level duplicate pair
#   against the local extraction of the original GEO RAW source.
#
# Comparison is performed on DECOMPRESSED CONTENT, not merely
# compressed-byte identity, for:
#   1. matrix
#   2. barcodes
#   3. features
#
# This stage NEVER modifies or deletes raw data.
# ============================================================

DUPLICATE_GSM = next(
    pair[1]
    for pair in SOURCE_DUPLICATE_PAIRS
    if pair[0] == SOURCE_DUPLICATE_KEEP
)
OUTPUT_FILE = RESULT_ROOT / "05A_original_geo_source_verification.txt"
SUMMARY_FILE = RESULT_ROOT / "05A_source_duplicate_verification_summary.csv"

def sha256_content(path):
    h = hashlib.sha256()
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

def find_gsm_files(gsm):
    result = {"matrix": [], "barcodes": [], "features": []}
    if not RAW_ROOT.exists():
        return result

    for path in RAW_ROOT.rglob("*"):
        if not path.is_file() or not path.name.startswith(gsm + "_"):
            continue
        name = path.name.lower()
        if "matrix" in name and name.endswith((".mtx", ".mtx.gz")):
            result["matrix"].append(path)
        elif "barcodes" in name and name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            result["barcodes"].append(path)
        elif "features" in name and name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            result["features"].append(path)
    return {k: sorted(v) for k, v in result.items()}

def main():
    print("=" * 78)
    print("05A — ORIGINAL GEO SOURCE VERIFICATION")
    print("=" * 78)
    print(f"\nResearch root:\n{RESEARCH_ROOT}")
    print(f"\nRAW source:\n{RAW_ROOT}")

    if not RAW_ROOT.exists():
        raise FileNotFoundError(f"GEO RAW directory not found:\n{RAW_ROOT}")

    keep_files = find_gsm_files(SOURCE_DUPLICATE_KEEP)
    duplicate_files = find_gsm_files(DUPLICATE_GSM)

    lines = [
        "=" * 78,
        "05A — ORIGINAL GEO SOURCE VERIFICATION",
        "=" * 78,
        "",
        f"RAW source: {RAW_ROOT}",
        "",
        f"Representative GSM : {SOURCE_DUPLICATE_KEEP}",
        f"Suspected duplicate : {DUPLICATE_GSM}",
        "",
    ]

    comparison = []
    all_identical = True
    files_complete = True

    for gsm, files in [
        (SOURCE_DUPLICATE_KEEP, keep_files),
        (DUPLICATE_GSM, duplicate_files),
    ]:
        print(f"\n{gsm}:")
        lines.append(f"{gsm}:")
        for kind in ("matrix", "barcodes", "features"):
            paths = files[kind]
            if len(paths) != 1:
                files_complete = False
                all_identical = False
                print(f"  {kind:<10} expected exactly 1 file; found {len(paths)}")
                lines.append(f"  {kind:<10} INVALID_FILE_COUNT={len(paths)}")
                for p in paths:
                    lines.append(f"             {p}")
            else:
                print(f"  {kind:<10} {paths[0].name}")
                lines.append(f"  {kind:<10} {paths[0].name}")

    lines += ["", "CONTENT COMPARISON", "-" * 78]

    for kind in ("matrix", "barcodes", "features"):
        a = keep_files[kind]
        b = duplicate_files[kind]

        if len(a) != 1 or len(b) != 1:
            status = "NOT_VERIFIED"
            hash_a = ""
            hash_b = ""
            identical = False
        else:
            hash_a = sha256_content(a[0])
            hash_b = sha256_content(b[0])
            identical = hash_a == hash_b
            status = "IDENTICAL" if identical else "DIFFERENT"

        if not identical:
            all_identical = False

        comparison.append({
            "file_type": kind,
            "keep_gsm": SOURCE_DUPLICATE_KEEP,
            "duplicate_gsm": DUPLICATE_GSM,
            "keep_content_sha256": hash_a,
            "duplicate_content_sha256": hash_b,
            "identical": identical,
            "status": status,
        })

        print(f"{kind.upper():<12} {status}")
        lines.append(f"{kind.upper():<12} {status}")
        if hash_a:
            lines.append(f"  {SOURCE_DUPLICATE_KEEP} SHA256    : {hash_a}")
            lines.append(f"  {DUPLICATE_GSM} SHA256: {hash_b}")

    confirmed = files_complete and all_identical

    if confirmed:
        decision = "SOURCE_LEVEL_DUPLICATE_CONFIRMED"
        decision_text = (
            f"CONFIRMED: {DUPLICATE_GSM} is an exact source-level duplicate "
            f"of {SOURCE_DUPLICATE_KEEP} across matrix, barcode and feature files."
        )
        print(f"\n{decision_text}")
        print(f"{SOURCE_DUPLICATE_KEEP} = RETAIN")
        print(f"{DUPLICATE_GSM} = EXCLUDE FROM BIOLOGICAL ANALYSIS")
        lines += [
            "",
            "-" * 78,
            decision_text,
            "",
            f"RETAIN  : {SOURCE_DUPLICATE_KEEP}",
            f"EXCLUDE : {DUPLICATE_GSM}",
            "",
            "The duplicate remains in the raw source for provenance and audit.",
            "It is excluded from downstream biological inference.",
        ]
    else:
        decision = "SOURCE_LEVEL_DUPLICATE_NOT_CONFIRMED"
        decision_text = (
            f"WARNING: {DUPLICATE_GSM} could not be confirmed as an exact "
            f"source-level duplicate of {SOURCE_DUPLICATE_KEEP}."
        )
        print(f"\n{decision_text}")
        lines += [
            "",
            "-" * 78,
            decision_text,
            "DO NOT exclude the duplicate until the discrepancy is resolved.",
        ]

    OUTPUT_FILE.write_text("\n".join(lines), encoding="utf-8")

    summary = pd.DataFrame([{
        "decision": decision,
        "keep_gsm": SOURCE_DUPLICATE_KEEP,
        "duplicate_gsm": DUPLICATE_GSM,
        "matrix_identical": next(r["identical"] for r in comparison if r["file_type"] == "matrix"),
        "barcodes_identical": next(r["identical"] for r in comparison if r["file_type"] == "barcodes"),
        "features_identical": next(r["identical"] for r in comparison if r["file_type"] == "features"),
        "all_three_identical": confirmed,
        "excluded_gsm": DUPLICATE_GSM if confirmed else "",
        "raw_data_modified": False,
        "verification_method": "decompressed_content_sha256_against_local_GEO_RAW_source",
    }])
    summary.to_csv(SUMMARY_FILE, index=False)

    comparison_df = pd.DataFrame(comparison)
    comparison_df.to_csv(
        RESULT_ROOT / "05A_source_duplicate_file_comparison.csv",
        index=False,
    )

    print(f"\nAudit report written to:\n{OUTPUT_FILE}")
    print(f"Machine-readable verification written to:\n{SUMMARY_FILE}")

    if not confirmed:
        raise RuntimeError("SOURCE-LEVEL DUPLICATE VERIFICATION FAILED")

if __name__ == "__main__":
    main()
