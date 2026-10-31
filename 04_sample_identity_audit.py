import re
import pandas as pd
from config import *

# ============================================================
# 04 — SAMPLE IDENTITY / DESIGN AUDIT
# ============================================================
# Audits ALL 28 original GEO samples.
#
# The audit verifies that raw filenames agree with the
# predefined experimental design in config.py.
#
# The source duplicate GSM8546673 remains visible here because
# this is an identity/source audit.
#
# It is excluded only from downstream biological analysis.
# ============================================================


def locate_files(gsm):

    result = {}

    for p in RAW_ROOT.rglob("*"):

        if not p.is_file():
            continue

        if gsm.lower() not in p.name.lower():
            continue

        n = p.name.lower()

        if "matrix" in n and n.endswith(
            (".mtx", ".mtx.gz")
        ):
            result["matrix"] = p

    return result


# ============================================================
# SOURCE + ANALYSIS SAMPLE SETS
# ============================================================

SOURCE_GSMS = tuple(SAMPLE_METADATA.keys())

ANALYSIS_GSMS = tuple(
    gsm
    for gsm in SOURCE_GSMS
    if gsm not in EXCLUDED_GSMS
)


# ============================================================
# AUDIT
# ============================================================

rows = []

for gsm in SOURCE_GSMS:

    f = locate_files(gsm)

    meta = SAMPLE_METADATA[gsm]

    included = gsm in ANALYSIS_GSMS

    exclusion_reason = ""

    if gsm in EXCLUDED_GSMS:
        exclusion_reason = (
            "Exact source-level duplicate of GSM8546664"
        )

    if "matrix" not in f:

        rows.append({
            "GSM": gsm,
            **meta,
            "included_in_analysis": included,
            "exclusion_reason": exclusion_reason,
            "filename": "",
            "inferred_time_hpi": None,
            "inferred_experiment": None,
            "inferred_replicate": None,
            "inferred_condition": None,
            "filename_agrees": False,
        })

        continue

    name = f["matrix"].name.lower()

    # --------------------------------------------------------
    # Infer time
    # --------------------------------------------------------

    inferred_time = None

    if "undamaged" in name:
        inferred_time = 0

    match = re.search(
        r"_(\d+)hpi_",
        name
    )

    if match:
        inferred_time = int(match.group(1))

    # --------------------------------------------------------
    # Infer experiment
    # --------------------------------------------------------

    if "_exp1_" in name:
        inferred_exp = "Exp1"

    elif "_exp2_" in name:
        inferred_exp = "Exp2"

    else:
        inferred_exp = None

    # --------------------------------------------------------
    # Infer replicate
    # --------------------------------------------------------

    if "_rep1_" in name:
        inferred_rep = 1

    elif "_rep2_" in name:
        inferred_rep = 2

    else:
        inferred_rep = None

    # --------------------------------------------------------
    # Infer condition
    # --------------------------------------------------------

    if "undamaged" in name:
        inferred_condition = "undamaged"

    else:
        inferred_condition = "injured"

    # --------------------------------------------------------
    # Compare against config
    # --------------------------------------------------------

    agree = (
        inferred_time == meta["time_hpi"]
        and inferred_exp == meta["experiment"]
        and inferred_rep == meta["replicate"]
        and inferred_condition == meta["condition"]
    )

    rows.append({
        "GSM": gsm,
        **meta,
        "included_in_analysis": included,
        "exclusion_reason": exclusion_reason,
        "filename": f["matrix"].name,
        "inferred_time_hpi": inferred_time,
        "inferred_experiment": inferred_exp,
        "inferred_replicate": inferred_rep,
        "inferred_condition": inferred_condition,
        "filename_agrees": agree,
    })


audit = pd.DataFrame(rows)

out = RESULT_ROOT / "04_sample_identity_audit.csv"

audit.to_csv(
    out,
    index=False
)


# ============================================================
# REPORT
# ============================================================

print("=" * 78)
print("04 — SAMPLE IDENTITY / DESIGN AUDIT")
print("=" * 78)

print()
print(audit.to_string(index=False))

print()
print("-" * 78)

if not audit["filename_agrees"].all():

    bad = audit[
        ~audit["filename_agrees"]
    ]

    print("FAILED SAMPLE IDENTITIES:")
    print(bad.to_string(index=False))

    raise SystemExit(
        "SAMPLE IDENTITY AUDIT FAILED — "
        "raw filenames do not match config.py design."
    )


source_count = len(SOURCE_GSMS)
analysis_count = len(ANALYSIS_GSMS)
excluded_count = len(EXCLUDED_GSMS)

print(
    f"PASS — all {source_count} original GEO sample identities "
    f"agree with the predefined experimental design."
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
