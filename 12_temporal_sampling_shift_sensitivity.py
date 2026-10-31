import numpy as np
import pandas as pd
from config import *

# ============================================================
# 12 — TEMPORAL SAMPLING / PHASE-SHIFT SENSITIVITY
# ============================================================
# Tests whether apparent timing differences depend on the unequal
# temporal grids of Exp1 and Exp2.
# ============================================================

e1 = pd.read_csv(CORRECTED_DISCOVERY_ROOT / "Exp1_gene_trajectories.csv")
e2 = pd.read_csv(CORRECTED_DISCOVERY_ROOT / "Exp2_gene_trajectories.csv")
OUT = CORRECTED_DISCOVERY_ROOT / "phase_shift_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

keep = ["gene", "max_time_hpi", "max_change", "pattern"]
a = e1[keep].rename(columns={c:f"{c}_Exp1" for c in keep if c!="gene"})
b = e2[keep].rename(columns={c:f"{c}_Exp2" for c in keep if c!="gene"})
df = a.merge(b, on="gene", how="inner")

# Full-grid peak shift.
df["full_grid_peak_shift_h"] = df["max_time_hpi_Exp2"] - df["max_time_hpi_Exp1"]
df["full_grid_abs_peak_shift_h"] = df["full_grid_peak_shift_h"].abs()

# Common-window peaks: use 0/24/72 trajectories only.
gene_data = pd.read_csv(RESULT_ROOT / "GSE278446_gene_level_sample_data.csv")
common = gene_data[gene_data["time_hpi"].isin(COMMON_TIMEPOINTS)]
tp = common.groupby(["gene","experiment","time_hpi"], as_index=False)["mean_expression"].mean()
p = tp.pivot_table(index="gene", columns=["experiment","time_hpi"], values="mean_expression")
for e in EXPECTED_EXPERIMENTS:
    cols = [f"{e}_{t}" for t in COMMON_TIMEPOINTS]
    # pivot column labels may be tuples; normalize below
# Rebuild explicitly.
rows=[]
for gene,g in tp.groupby("gene"):
    q=g.pivot(index="time_hpi",columns="experiment",values="mean_expression").reindex(COMMON_TIMEPOINTS)
    if {"Exp1","Exp2"}.issubset(q.columns) and not q[["Exp1","Exp2"]].isna().any().any():
        for exp in EXPECTED_EXPERIMENTS:
            vals=q[exp].to_numpy(float)
            peak_time=COMMON_TIMEPOINTS[int(np.argmax(vals))]
            rows.append({"gene":gene, f"common_peak_time_{exp}":peak_time})
common_peaks=pd.DataFrame(rows).groupby("gene").first().reset_index()
df=df.merge(common_peaks,on="gene",how="left")
df["common_window_peak_shift_h"] = df["common_peak_time_Exp2"] - df["common_peak_time_Exp1"]
df["common_window_abs_peak_shift_h"] = df["common_window_peak_shift_h"].abs()
df["peak_shift_changed_by_sampling"] = (
    df["full_grid_peak_shift_h"].fillna(0) != df["common_window_peak_shift_h"].fillna(0)
)
df["pattern_agreement"] = df.apply(
    lambda r: r["pattern_Exp1"] == r["pattern_Exp2"] and r["pattern_Exp1"] not in ["constant","other"],
    axis=1
)

df.to_csv(OUT / "ALL_PHASE_SHIFT_SENSITIVITY.csv", index=False)

summary = pd.DataFrame([{
    "genes_compared": len(df),
    "genes_with_changed_peak_shift_under_common_window": int(df["peak_shift_changed_by_sampling"].sum()),
    "fraction_changed": float(df["peak_shift_changed_by_sampling"].mean()) if len(df) else np.nan,
    "median_full_grid_abs_shift_h": df["full_grid_abs_peak_shift_h"].median(),
    "median_common_window_abs_shift_h": df["common_window_abs_peak_shift_h"].median(),
    "pattern_agreement_fraction": df["pattern_agreement"].mean() if len(df) else np.nan,
}])
summary.to_csv(OUT / "PHASE_SHIFT_SENSITIVITY_SUMMARY.csv", index=False)

print("=" * 78)
print("12 — TEMPORAL SAMPLING SENSITIVITY COMPLETE")
print(summary.to_string(index=False))
print(f"Output: {OUT}")
