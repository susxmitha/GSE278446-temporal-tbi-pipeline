import pandas as pd
import numpy as np
from config import *

# ============================================================
# 11 — STRICT 3/3 AVAILABLE COMMON-WINDOW VALIDATION
# ============================================================
# GSM8546673 is excluded because GEO source verification shows it
# is a duplicate of GSM8546664. The 0/24/72 hpi common window
# therefore contains 3 complete independent trajectories:
# Exp1_R2, Exp2_R1, Exp2_R2.
#
# A strict candidate must show the same 24 h peak/dip pattern and
# pass the 20% early/late magnitude threshold in ALL 3 trajectories.
# ============================================================

gene_dir = RESULT_ROOT / "gene_means"
OUT = REPLICATE_VALIDATION_ROOT
OUT.mkdir(parents=True, exist_ok=True)

records = []
for gsm in ANALYSIS_GSMS:
    meta = SAMPLE_METADATA[gsm]
    if meta["time_hpi"] not in COMMON_TIMEPOINTS:
        continue
    path = gene_dir / f"{gsm}_gene_means.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    x = pd.read_csv(path, usecols=["gene", "mean_expression"])
    x["gene"] = x["gene"].astype(str).str.strip().str.lower()
    x["GSM"] = gsm
    x["time_hpi"] = meta["time_hpi"]
    x["experiment"] = meta["experiment"]
    x["replicate"] = meta["replicate"]
    records.append(x)

data = pd.concat(records, ignore_index=True)

wide = data.pivot_table(
    index="gene",
    columns=["experiment", "replicate", "time_hpi"],
    values="mean_expression",
    aggfunc="mean",
)
wide.columns = [f"{e}_R{r}_{int(t)}" for e, r, t in wide.columns]
wide = wide.reset_index()

trajectory_keys = [
    ("Exp1", 2),
    ("Exp2", 1),
    ("Exp2", 2),
]

required = [
    f"{e}_R{r}_{t}"
    for e, r in trajectory_keys
    for t in COMMON_TIMEPOINTS
]
missing_cols = [c for c in required if c not in wide.columns]
if missing_cols:
    raise RuntimeError(f"Missing common-window trajectory columns: {missing_cols}")

wide = wide.dropna(subset=required).copy()

def pattern(v0, v24, v72):
    early, late = v24 - v0, v72 - v24
    if early > 0 and late < 0:
        return "peak_at_24"
    if early < 0 and late > 0:
        return "dip_at_24"
    return "none"

for exp, rep in trajectory_keys:
    prefix = f"{exp}_R{rep}"
    wide[f"{prefix}_pattern"] = [
        pattern(a, b, c)
        for a, b, c in zip(
            wide[f"{prefix}_0"],
            wide[f"{prefix}_24"],
            wide[f"{prefix}_72"],
        )
    ]
    wide[f"{prefix}_relative_early"] = (
        (wide[f"{prefix}_24"] - wide[f"{prefix}_0"])
        / wide[f"{prefix}_0"].abs().clip(lower=1e-8)
    )
    wide[f"{prefix}_relative_late"] = (
        (wide[f"{prefix}_72"] - wide[f"{prefix}_24"])
        / wide[f"{prefix}_24"].abs().clip(lower=1e-8)
    )
    wide[f"{prefix}_magnitude_pass"] = (
        wide[f"{prefix}_relative_early"].abs() >= MIN_RELATIVE_CHANGE
    ) & (
        wide[f"{prefix}_relative_late"].abs() >= MIN_RELATIVE_CHANGE
    )

pattern_cols = [f"{e}_R{r}_pattern" for e, r in trajectory_keys]
magnitude_cols = [f"{e}_R{r}_magnitude_pass" for e, r in trajectory_keys]

wide["n_complete_trajectories"] = len(trajectory_keys)
wide["supporting_trajectories"] = wide[pattern_cols].isin(
    ["peak_at_24", "dip_at_24"]
).sum(axis=1)
wide["all_same_pattern"] = (
    wide[pattern_cols].nunique(axis=1).eq(1)
    & wide[pattern_cols[0]].isin(["peak_at_24", "dip_at_24"])
)
wide["magnitude_pass_all_available"] = wide[magnitude_cols].all(axis=1)
wide["strict_3_of_3_available"] = (
    wide["all_same_pattern"] & wide["magnitude_pass_all_available"]
)
wide["consensus_pattern"] = wide[pattern_cols].apply(
    lambda r: r.iloc[0] if r.nunique() == 1 else "mixed",
    axis=1,
)
wide["replicate_consensus"] = np.select(
    [
        wide["strict_3_of_3_available"] & wide["consensus_pattern"].eq("peak_at_24"),
        wide["strict_3_of_3_available"] & wide["consensus_pattern"].eq("dip_at_24"),
        wide["supporting_trajectories"] >= 2,
    ],
    ["3/3_peak", "3/3_dip", "2/3_reversal"],
    default="mixed",
)

wide.to_csv(OUT / "ALL_AVAILABLE_TRAJECTORY_VALIDATION.csv", index=False)
wide[wide["strict_3_of_3_available"]].to_csv(
    OUT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv", index=False
)
wide[wide["supporting_trajectories"] >= 2].to_csv(
    OUT / "2_OF_3_REVERSALS.csv", index=False
)

print("=" * 78)
print("11 — STRICT 3/3 AVAILABLE TRAJECTORY VALIDATION")
print(f"Complete independent trajectories: {len(trajectory_keys)}")
print(f"Strict 3/3 candidates: {int(wide['strict_3_of_3_available'].sum())}")
print(f"Output: {OUT}")
