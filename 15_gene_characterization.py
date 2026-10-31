import numpy as np
import pandas as pd
from config import *

# ============================================================
# 15 — STRICT-CANDIDATE GENE CHARACTERIZATION
# ============================================================

OUT = GENE_CHARACTERIZATION_ROOT
OUT.mkdir(parents=True, exist_ok=True)
strict = pd.read_csv(REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv")
e1 = pd.read_csv(CORRECTED_DISCOVERY_ROOT / "Exp1_gene_trajectories.csv")
e2 = pd.read_csv(CORRECTED_DISCOVERY_ROOT / "Exp2_gene_trajectories.csv")

strict["gene"]=strict["gene"].astype(str).str.strip().str.lower()
e1["gene"]=e1["gene"].astype(str).str.strip().str.lower()
e2["gene"]=e2["gene"].astype(str).str.strip().str.lower()

# Keep only one row per gene in each table.
e1=e1.drop_duplicates("gene")
e2=e2.drop_duplicates("gene")

a=e1.add_suffix("_Exp1").rename(columns={"gene_Exp1":"gene"})
b=e2.add_suffix("_Exp2").rename(columns={"gene_Exp2":"gene"})
merged=strict[["gene","replicate_consensus","supporting_trajectories","strict_3_of_3_available"]].merge(a,on="gene",how="inner").merge(b,on="gene",how="inner")

merged["peak_time_difference"]=(
    pd.to_numeric(merged["max_time_hpi_Exp1"],errors="coerce")
    - pd.to_numeric(merged["max_time_hpi_Exp2"],errors="coerce")
).abs()

merged["mean_dynamic_range"]=merged[["dynamic_range_Exp1","dynamic_range_Exp2"]].mean(axis=1)
merged["mean_transient_score"]=merged[["transient_score_Exp1","transient_score_Exp2"]].mean(axis=1)
merged["mean_recovery"]=merged[["recovery_fraction_Exp1","recovery_fraction_Exp2"]].mean(axis=1)
merged["rho_agreement"]=np.sign(merged["spearman_rho_Exp1"])==np.sign(merged["spearman_rho_Exp2"])
merged["pattern_agreement"]=merged["pattern_Exp1"]==merged["pattern_Exp2"]
merged["temporal_effect_consistency"]=(
    np.minimum(merged["max_change_Exp1"].abs(),merged["max_change_Exp2"].abs())
    / np.maximum(merged["max_change_Exp1"].abs(),merged["max_change_Exp2"].abs()).clip(lower=1e-12)
)

def characterize(r):
    p1,p2=str(r["pattern_Exp1"]),str(r["pattern_Exp2"])
    if p1=="early_transient" and p2=="early_transient": return "reproducible_early_transient"
    if p1=="midcourse_peak" and p2=="midcourse_peak": return "reproducible_midcourse_peak"
    if p1=="late_emergent" and p2=="late_emergent": return "reproducible_late"
    if p1=="persistent" and p2=="persistent": return "reproducible_persistent"
    if r["pattern_agreement"]: return "same_pattern_other"
    return "mixed_pattern"
merged["characterization"]=merged.apply(characterize,axis=1)

def pct_rank(s,ascending=True):
    return s.rank(pct=True,ascending=ascending,method="average").fillna(0)
merged["range_rank"]=pct_rank(merged["mean_dynamic_range"])
merged["recovery_rank"]=pct_rank(merged["mean_recovery"])
merged["effect_consistency_rank"]=pct_rank(merged["temporal_effect_consistency"])
merged["peak_agreement_rank"]=1/(1+merged["peak_time_difference"].fillna(np.inf))
merged["temporal_priority_score"]=(
    30*merged["pattern_agreement"].astype(int)
    +20*merged["range_rank"]
    +20*merged["effect_consistency_rank"]
    +15*merged["recovery_rank"]
    +15*merged["peak_agreement_rank"]
)

merged=merged.sort_values(["temporal_priority_score","gene"],ascending=[False,True])
merged.to_csv(OUT/"ALL_STRICT_GENE_CHARACTERIZATION.csv",index=False)
merged.head(50).to_csv(OUT/"TOP_50_PRIORITY_GENES.csv",index=False)
merged.head(20).to_csv(OUT/"TOP_20_PRIORITY_GENES.csv",index=False)

print("="*78)
print("15 — GENE CHARACTERIZATION COMPLETE")
print(f"Strict genes: {len(merged)}")
print(f"Output: {OUT}")
