import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from config import *

# ============================================================
# 07 — PROGRAM-LEVEL TEMPORAL DISCOVERY
# ============================================================

INPUT = RESULT_ROOT / "GSE278446_sample_program_summary.csv"
OUT = DISCOVERY_ROOT / "program_temporal"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
programs = [p for p in MARKER_PROGRAMS if f"{p}_mean" in df.columns]
df["time_hpi"] = pd.to_numeric(df["time_hpi"], errors="coerce")

# Experiment-specific trajectory means across the two biological replicates.
trajectory_rows = []
for (exp, time), g in df.groupby(["experiment", "time_hpi"]):
    row = {"experiment": exp, "time_hpi": time, "n_replicates": g["replicate"].nunique()}
    for p in programs:
        row[f"{p}_mean"] = g[f"{p}_mean"].mean()
        row[f"{p}_sd"] = g[f"{p}_mean"].std(ddof=1)
    trajectory_rows.append(row)
trajectory = pd.DataFrame(trajectory_rows).sort_values(["experiment", "time_hpi"])
trajectory.to_csv(OUT / "PROGRAM_TEMPORAL_TRAJECTORIES.csv", index=False)

# Peak time per experiment.
peak_rows = []
for exp, g in trajectory.groupby("experiment"):
    for p in programs:
        c = f"{p}_mean"
        if g[c].notna().sum() == 0:
            continue
        idx = g[c].idxmax()
        peak_rows.append({
            "experiment": exp,
            "program": p,
            "peak_time_hpi": float(g.loc[idx, "time_hpi"]),
            "peak_score": float(g.loc[idx, c]),
        })
peaks = pd.DataFrame(peak_rows)
peaks.to_csv(OUT / "PROGRAM_PEAK_ORDER.csv", index=False)

# Pairwise temporal correlations are descriptive because the number of
# timepoints is small. BH correction is applied within each experiment.
pair_rows = []
for exp, g in trajectory.groupby("experiment"):
    g = g.sort_values("time_hpi")
    for i, a in enumerate(programs):
        for b in programs[i+1:]:
            valid = g[[f"{a}_mean", f"{b}_mean"]].dropna()
            if len(valid) < 4 or valid.iloc[:, 0].nunique() < 2 or valid.iloc[:, 1].nunique() < 2:
                rho, p = np.nan, np.nan
            else:
                rho, p = spearmanr(valid.iloc[:, 0], valid.iloc[:, 1])
            pair_rows.append({
                "experiment": exp, "program_A": a, "program_B": b,
                "n_timepoints": len(valid), "spearman_rho": rho, "p_value": p,
            })
pairs = pd.DataFrame(pair_rows)
for exp in pairs["experiment"].unique():
    mask = pairs["experiment"].eq(exp) & pairs["p_value"].notna()
    if mask.any():
        pairs.loc[mask, "q_value"] = multipletests(
            pairs.loc[mask, "p_value"], method="fdr_bh"
        )[1]
pairs.to_csv(OUT / "PAIRWISE_PROGRAM_TEMPORAL_RELATIONSHIPS.csv", index=False)

# Early/late descriptive contrast.
early_times = {0, 5, 12, 18, 24}
early_late_rows = []
for exp, g in df.groupby("experiment"):
    for p in programs:
        c = f"{p}_mean"
        early = g[g["time_hpi"].isin(early_times)][c].dropna()
        late = g[g["time_hpi"] >= 48][c].dropna()
        if len(early) >= 2 and len(late) >= 2:
            stat, pval = mannwhitneyu(early, late, alternative="two-sided")
        else:
            stat, pval = np.nan, np.nan
        early_late_rows.append({
            "experiment": exp, "program": p, "n_early": len(early), "n_late": len(late),
            "early_mean": early.mean() if len(early) else np.nan,
            "late_mean": late.mean() if len(late) else np.nan,
            "difference_late_minus_early": (late.mean()-early.mean()) if len(early) and len(late) else np.nan,
            "mannwhitney_u": stat, "p_value": pval,
        })
el = pd.DataFrame(early_late_rows)
for exp in el["experiment"].unique():
    mask = el["experiment"].eq(exp) & el["p_value"].notna()
    if mask.any():
        el.loc[mask, "q_value"] = multipletests(el.loc[mask, "p_value"], method="fdr_bh")[1]
el.to_csv(OUT / "PROGRAM_EARLY_LATE_CONTRAST.csv", index=False)

print("=" * 78)
print("07 — PROGRAM TEMPORAL DISCOVERY COMPLETE")
print(f"Programs analysed: {len(programs)}")
print(f"Outputs: {OUT}")
