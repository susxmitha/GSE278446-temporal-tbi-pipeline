import numpy as np
import pandas as pd
from scipy.stats import spearmanr, linregress
from statsmodels.stats.multitest import multipletests
from config import *

# ============================================================
# 08 — GENE-LEVEL TEMPORAL DISCOVERY
# ============================================================
# This is exploratory discovery across the full sampling grids.
# The definitive 0/24/72 reversal is handled independently by 10/11.
# ============================================================

INPUT = RESULT_ROOT / "GSE278446_gene_level_sample_data.csv"
OUT = CORRECTED_DISCOVERY_ROOT
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
required = ["gene", "GSM", "time_hpi", "experiment", "replicate", "mean_expression"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df["time_hpi"] = pd.to_numeric(df["time_hpi"], errors="coerce")
df["mean_expression"] = pd.to_numeric(df["mean_expression"], errors="coerce")
df = df.dropna(subset=["gene", "time_hpi", "mean_expression", "experiment"]).copy()
df["gene"] = df["gene"].astype(str).str.strip().str.lower()

# Average biological replicates only after retaining replicate identity.
tp = (
    df.groupby(["gene", "experiment", "time_hpi"], as_index=False)["mean_expression"]
      .mean()
)

def analyze(exp):
    sub = tp[tp["experiment"] == exp]
    rows = []
    total_genes = sub["gene"].nunique()
    print(f"    {exp}: analyzing {total_genes:,} genes ...", flush=True)
    for i, (gene, g) in enumerate(sub.groupby("gene", sort=False), start=1):
        if i % 5000 == 0:
            print(f"    {exp}: {i:,}/{total_genes:,}", flush=True)
        s = g.groupby("time_hpi")["mean_expression"].mean().sort_index()
        if len(s) < MIN_TEMPORAL_POINTS:
            continue
        t = s.index.to_numpy(dtype=float)
        y = s.to_numpy(dtype=float)
        if np.ptp(y) <= 1e-12:
            rho = rho_p = slope = slope_p = np.nan
            max_time = float(t[0])
            max_value = float(y[0])
            max_change = 0.0
        else:
            rho, rho_p = spearmanr(t, y) if np.ptp(t) > 0 else (np.nan, np.nan)
            lr = linregress(t, y)
            slope, slope_p = float(lr.slope), float(lr.pvalue)
            baseline = float(y[0] if 0 not in s.index else s.loc[0])
            change = y - baseline
            idx = int(np.argmax(np.abs(change)))
            max_time, max_value, max_change = float(t[idx]), float(y[idx]), float(change[idx])

        baseline = float(s.loc[0]) if 0 in s.index else float(y[0])
        early = float(y[t <= 24].mean()) if np.any(t <= 24) else np.nan
        late = float(y[t >= 72].mean()) if np.any(t >= 72) else np.nan
        early_change = early - baseline if np.isfinite(early) else np.nan
        late_change = late - baseline if np.isfinite(late) else np.nan
        end_change = float(y[-1] - baseline)
        dynamic_range = float(np.max(y) - np.min(y))
        max_deviation = float(np.max(np.abs(y - baseline)))
        recovery = 1 - abs(end_change) / max_deviation if max_deviation > 0 else 0.0
        recovery = float(np.clip(recovery, -1, 1))

        if max_time <= 24 and recovery >= 0.50:
            pattern = "early_transient"
        elif max_time >= 72 and np.isfinite(early_change) and np.isfinite(late_change) and abs(late_change) > 1.25 * abs(early_change):
            pattern = "late_emergent"
        elif max_deviation > 0 and abs(end_change) >= 0.70 * max_deviation:
            pattern = "persistent"
        elif len(y) >= 5 and 0 < np.argmax(np.abs(y-baseline)) < len(y)-1:
            pattern = "midcourse_peak"
        else:
            pattern = "other"

        rows.append({
            "gene": gene, "experiment": exp, "n_timepoints": len(t),
            "baseline": baseline, "max_time_hpi": max_time, "max_value": max_value,
            "max_change": max_change, "relative_change": max_change / max(abs(baseline), 1e-8),
            "dynamic_range": dynamic_range, "spearman_rho": rho, "spearman_p": rho_p,
            "linear_slope": slope, "linear_p": slope_p,
            "early_mean": early, "late_mean": late,
            "early_change": early_change, "late_change": late_change,
            "transient_score": (abs(early_change)-abs(late_change)) if np.isfinite(early_change) and np.isfinite(late_change) else np.nan,
            "late_score": (abs(late_change)-abs(early_change)) if np.isfinite(early_change) and np.isfinite(late_change) else np.nan,
            "recovery_fraction": recovery, "pattern": pattern,
        })
    return pd.DataFrame(rows)

exp1 = analyze("Exp1")
exp2 = analyze("Exp2")
exp1.to_csv(OUT / "Exp1_gene_trajectories.csv", index=False)
exp2.to_csv(OUT / "Exp2_gene_trajectories.csv", index=False)

for result, exp in [(exp1, "Exp1"), (exp2, "Exp2")]:
    for col in ["spearman_p", "linear_p"]:
        mask = result[col].notna()
        result.loc[mask, f"{col}_q"] = multipletests(result.loc[mask, col], method="fdr_bh")[1]
    result.to_csv(OUT / f"{exp}_gene_trajectories.csv", index=False)

# Independent common-time validation.
shared = tp[tp["time_hpi"].isin(COMMON_TIMEPOINTS)]
rows = []
for gene, g in shared.groupby("gene"):
    p = g.pivot(index="time_hpi", columns="experiment", values="mean_expression").reindex(COMMON_TIMEPOINTS)
    if not {"Exp1", "Exp2"}.issubset(p.columns) or p[["Exp1", "Exp2"]].isna().any().any():
        continue
    x, y = p["Exp1"].to_numpy(float), p["Exp2"].to_numpy(float)
    rho, pval = spearmanr(x, y) if np.ptp(x) and np.ptp(y) else (np.nan, np.nan)
    dx, dy = np.diff(x), np.diff(y)
    rows.append({
        "gene": gene, "rho_shared_times": rho, "p_shared_times": pval,
        "direction_agreement": float(np.mean(np.sign(dx) == np.sign(dy))),
        "Exp1_24h_change": dx[0], "Exp1_72h_change": x[2]-x[0],
        "Exp2_24h_change": dy[0], "Exp2_72h_change": y[2]-y[0],
    })
validation = pd.DataFrame(rows)
if len(validation):
    mask = validation["p_shared_times"].notna()
    validation.loc[mask, "q_shared_times"] = multipletests(validation.loc[mask, "p_shared_times"], method="fdr_bh")[1]
validation.to_csv(OUT / "shared_timepoint_validation.csv", index=False)

# Merge the two full trajectories.
e1 = exp1.add_suffix("_Exp1").rename(columns={"gene_Exp1": "gene"})
e2 = exp2.add_suffix("_Exp2").rename(columns={"gene_Exp2": "gene"})
merged = e1.merge(e2, on="gene", how="inner").merge(validation, on="gene", how="left")

early_patterns = {"early_transient", "midcourse_peak"}
late_patterns = {"late_emergent", "persistent"}
def compatible(a, b):
    return bool(a == b or (a in early_patterns and b in early_patterns) or (a in late_patterns and b in late_patterns))
merged["pattern_agreement"] = [
    compatible(a, b) for a, b in zip(merged["pattern_Exp1"], merged["pattern_Exp2"])
]
merged["trajectory_strength"] = (
    merged[["spearman_rho_Exp1", "spearman_rho_Exp2"]].abs().mean(axis=1)
)
merged["effect_strength"] = np.minimum(
    merged["dynamic_range_Exp1"].fillna(0), merged["dynamic_range_Exp2"].fillna(0)
)
merged["candidate_score"] = (
    merged["trajectory_strength"].fillna(0)
    * merged["effect_strength"].clip(lower=0).map(np.log1p)
    * merged["direction_agreement"].fillna(0)
)
merged.to_csv(OUT / "ALL_CROSS_EXPERIMENT_GENE_RESULTS.csv", index=False)

candidates = merged[
    merged["pattern_agreement"]
    & (merged["direction_agreement"].fillna(0) >= 0.5)
    & (merged["trajectory_strength"].fillna(0) >= 0.3)
    & (merged["effect_strength"].fillna(0) > 0)
].sort_values("candidate_score", ascending=False)
candidates.to_csv(OUT / "EXPLORATORY_TEMPORAL_CANDIDATES.csv", index=False)

print("=" * 78)
print("08 — GENE TEMPORAL DISCOVERY COMPLETE")
print(f"Exp1 genes: {len(exp1)}")
print(f"Exp2 genes: {len(exp2)}")
print(f"Exploratory candidates: {len(candidates)}")
print(f"Output: {OUT}")
