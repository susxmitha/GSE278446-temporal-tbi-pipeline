import pandas as pd
import numpy as np
from config import *

# ============================================================
# 10 — COMMON-WINDOW (0/24/72 hpi) REVERSAL DISCOVERY
# ============================================================
# Experiment-level discovery. Script 11 independently validates
# the same phenomenon at all four biological replicates.
# ============================================================

INPUT = RESULT_ROOT / "GSE278446_gene_level_sample_data.csv"
OUT = COMMON_WINDOW_ROOT
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
df = df[df["time_hpi"].isin(COMMON_TIMEPOINTS)].copy()
df["gene"] = df["gene"].astype(str).str.strip().str.lower()

# Mean of the available independent biological samples within each experiment/timepoint.
tp = (
    df.groupby(["gene", "experiment", "time_hpi"], as_index=False)
      .agg(
          mean_expression=("mean_expression", "mean"),
          n_biological_samples=("GSM", "nunique"),
      )
)
wide = tp.pivot_table(
    index="gene",
    columns=["experiment", "time_hpi"],
    values="mean_expression",
    aggfunc="mean",
)
wide.columns = [f"{e}_{int(t)}" for e, t in wide.columns]
wide = wide.reset_index()

nwide = tp.pivot_table(
    index="gene",
    columns=["experiment", "time_hpi"],
    values="n_biological_samples",
    aggfunc="max",
)
nwide.columns = [f"{e}_{int(t)}_n" for e, t in nwide.columns]
nwide = nwide.reset_index()
wide = wide.merge(nwide, on="gene", how="left")

required = [f"{e}_{t}" for e in EXPECTED_EXPERIMENTS for t in COMMON_TIMEPOINTS]
wide = wide.dropna(subset=required).copy()

for exp in EXPECTED_EXPERIMENTS:
    wide[f"{exp}_early_change"] = wide[f"{exp}_24"] - wide[f"{exp}_0"]
    wide[f"{exp}_late_change"] = wide[f"{exp}_72"] - wide[f"{exp}_24"]
    wide[f"{exp}_relative_early"] = wide[f"{exp}_early_change"] / wide[f"{exp}_0"].abs().clip(lower=1e-8)
    wide[f"{exp}_relative_late"] = wide[f"{exp}_late_change"] / wide[f"{exp}_24"].abs().clip(lower=1e-8)

def reversal(early, late):
    if early > 0 and late < 0:
        return "peak_at_24"
    if early < 0 and late > 0:
        return "dip_at_24"
    return "none"

for exp in EXPECTED_EXPERIMENTS:
    wide[f"{exp}_pattern"] = [
        reversal(a, b) for a, b in zip(wide[f"{exp}_early_change"], wide[f"{exp}_late_change"])
    ]
    wide[f"{exp}_strong"] = (
        wide[f"{exp}_relative_early"].abs() >= MIN_RELATIVE_CHANGE
    ) & (
        wide[f"{exp}_relative_late"].abs() >= MIN_RELATIVE_CHANGE
    )

wide["pattern_agreement"] = (
    (wide["Exp1_pattern"] == wide["Exp2_pattern"])
    & wide["Exp1_pattern"].isin(["peak_at_24", "dip_at_24"])
)
wide["strong_reproducible"] = (
    wide["pattern_agreement"] & wide["Exp1_strong"] & wide["Exp2_strong"]
)
wide["reproducible_strength"] = np.minimum(
    wide["Exp1_early_change"].abs() + wide["Exp1_late_change"].abs(),
    wide["Exp2_early_change"].abs() + wide["Exp2_late_change"].abs(),
)

wide.to_csv(OUT / "COMMON_WINDOW_REVERSAL_ALL_GENES.csv", index=False)
strong = wide[wide["strong_reproducible"]].sort_values("reproducible_strength", ascending=False)
strong.to_csv(OUT / "STRONG_REPRODUCIBLE_REVERSALS.csv", index=False)

print("=" * 78)
print("10 — COMMON-WINDOW REVERSAL")
print(f"Complete trajectories: {len(wide)}")
print(f"Strong reproducible reversals: {len(strong)}")
print(f"Output: {OUT}")
