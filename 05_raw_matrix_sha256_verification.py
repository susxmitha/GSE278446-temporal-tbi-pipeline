import gzip
import hashlib
import pandas as pd
from config import *

# ============================================================
# 05 — FULL RAW MATRIX SHA-256 VERIFICATION
# ============================================================
# Verifies every original GEO GSM has exactly one raw matrix.
#
# Two fingerprints are recorded:
#   compressed_sha256 — exact local file bytes
#   content_sha256    — decompressed matrix content
#
# Duplicate detection uses content_sha256 so gzip metadata cannot
# hide a biologically identical matrix.
#
# The known duplicate pair is accepted ONLY after Stage 05A has
# produced a machine-readable source-verification record.
#
# No raw files are modified or deleted.
# ============================================================

SOURCE_VERIFY = RESULT_ROOT / "05A_source_duplicate_verification_summary.csv"

def matrix_files_for(gsm):
    out = []
    for p in RAW_ROOT.rglob("*"):
        if not p.is_file():
            continue
        name = p.name.lower()
        if gsm.lower() in p.name.lower() and name.endswith((".mtx", ".mtx.gz")):
            out.append(p)
    return sorted(out)

def sha256_file(path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

def sha256_content(path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

def require_source_verification():
    if not SOURCE_VERIFY.exists():
        raise RuntimeError(
            "Stage 05A verification record is missing. "
            "Run 05A_verify_duplicate_against_geo.py first."
        )

    verify = pd.read_csv(SOURCE_VERIFY, dtype=str)
    required = {"decision", "keep_gsm", "duplicate_gsm", "all_three_identical"}
    missing = required - set(verify.columns)
    if missing:
        raise RuntimeError(
            f"Stage 05A verification record is incomplete; missing columns: {sorted(missing)}"
        )

    if verify.empty:
        raise RuntimeError("Stage 05A verification record is empty.")

    row = verify.iloc[0]
    if str(row["decision"]) != "SOURCE_LEVEL_DUPLICATE_CONFIRMED":
        raise RuntimeError(
            "Stage 05A did not confirm the source-level duplicate."
        )

    if str(row["keep_gsm"]) != SOURCE_DUPLICATE_KEEP:
        raise RuntimeError("Stage 05A keep GSM does not match config.py.")

    duplicate_gsm = str(row["duplicate_gsm"])
    if duplicate_gsm not in EXCLUDED_GSMS:
        raise RuntimeError(
            "Stage 05A duplicate GSM is not present in EXCLUDED_GSMS."
        )

    if str(row["all_three_identical"]).strip().lower() != "true":
        raise RuntimeError(
            "Stage 05A did not confirm all three source files as identical."
        )

    return duplicate_gsm

duplicate_gsm = require_source_verification()

rows = []
for gsm in SOURCE_GSMS:
    paths = matrix_files_for(gsm)
    if len(paths) != 1:
        raise RuntimeError(
            f"{gsm}: expected exactly one matrix, found {len(paths)}"
        )

    p = paths[0]
    rows.append({
        "GSM": gsm,
        "matrix_file": p.name,
        "path": str(p),
        "size_bytes": p.stat().st_size,
        "compressed_sha256": sha256_file(p),
        "content_sha256": sha256_content(p),
    })

df = pd.DataFrame(rows)

duplicate_groups = (
    df.groupby("content_sha256")["GSM"]
      .apply(list)
)
duplicate_groups = duplicate_groups[duplicate_groups.apply(len) > 1]

df["integrity_status"] = "UNIQUE_MATRIX"

for digest, gsm_list in duplicate_groups.items():
    gsm_set = set(gsm_list)
    expected_pair_matches = any(
        set(pair).issubset(gsm_set) for pair in SOURCE_DUPLICATE_PAIRS
    )
    status = "SOURCE_LEVEL_DUPLICATE" if expected_pair_matches else "UNRESOLVED_DUPLICATE"
    df.loc[df["GSM"].isin(gsm_list), "integrity_status"] = status

out = RESULT_ROOT / "05_full_sha256_matrix_verification.csv"
df.to_csv(out, index=False)

duplicate_report = RESULT_ROOT / "05_duplicate_matrix_hashes.csv"
if len(duplicate_groups):
    dup_rows = []
    for digest, gsm_list in duplicate_groups.items():
        for gsm in gsm_list:
            row = df[df["GSM"] == gsm].iloc[0]
            dup_rows.append({
                "content_sha256": digest,
                "GSM": gsm,
                "matrix_file": row["matrix_file"],
                "size_bytes": row["size_bytes"],
                "integrity_status": row["integrity_status"],
            })
    pd.DataFrame(dup_rows).to_csv(duplicate_report, index=False)
else:
    pd.DataFrame(
        columns=["content_sha256", "GSM", "matrix_file", "size_bytes", "integrity_status"]
    ).to_csv(duplicate_report, index=False)

unresolved = df[df["integrity_status"] == "UNRESOLVED_DUPLICATE"]
source_duplicates = df[df["integrity_status"] == "SOURCE_LEVEL_DUPLICATE"]

print("=" * 78)
print("05 — FULL SHA-256 MATRIX VERIFICATION")
print("=" * 78)
print()
print(df[
    ["GSM", "matrix_file", "size_bytes", "compressed_sha256",
     "content_sha256", "integrity_status"]
].to_string(index=False))

print()
print("=" * 78)

if len(unresolved):
    print("FAIL — unresolved duplicate matrix content detected.")
    print(unresolved[["GSM", "content_sha256"]].to_string(index=False))
    raise RuntimeError("Unresolved duplicate matrix content requires source verification.")

expected_duplicate_gsms = {SOURCE_DUPLICATE_KEEP, duplicate_gsm}
observed_duplicate_gsms = set(source_duplicates["GSM"])
if observed_duplicate_gsms != expected_duplicate_gsms:
    raise RuntimeError(
        "Source-level duplicate set does not match the verified duplicate pair: "
        f"expected {sorted(expected_duplicate_gsms)}, observed {sorted(observed_duplicate_gsms)}"
    )

print("PASS WITH EXCEPTION — SOURCE-LEVEL DUPLICATE CONFIRMED.")
print(
    f"Verified duplicate pair: {SOURCE_DUPLICATE_KEEP} / {duplicate_gsm}"
)
print(
    "The duplicate is retained in raw provenance but excluded from biological inference."
)
print("No raw files were modified or deleted.")
print(f"\nOutput: {out}")
print(f"Duplicate report: {duplicate_report}")
