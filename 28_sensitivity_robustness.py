import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from config import *

# ============================================================
# 28 — TEMPORAL RANKING SENSITIVITY / ROBUSTNESS
# ============================================================
# Pathway enrichment is retained as list-level evidence only.
# It is NOT converted into a gene-specific score. Robustness is
# therefore tested by perturbing temporal-score weights.
# ============================================================

IN = PRIORITIZATION_ROOT / "ALL_STRICT_GENE_PRIORITIZATION.csv"
OUT = SENSITIVITY_ROOT
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN)
required = [
    "gene", "pattern_agreement", "mean_dynamic_range",
    "temporal_effect_consistency", "mean_recovery",
    "peak_agreement_rank",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

def pct(s):
    return pd.to_numeric(s, errors="coerce").rank(
        pct=True, method="average"
    ).fillna(0)

def parse_bool_series(s):
    mapping = {
        "true": True, "1": True, "yes": True, "y": True,
        "false": False, "0": False, "no": False, "n": False,
    }
    return (
        s.astype(str)
         .str.strip()
         .str.lower()
         .map(mapping)
         .fillna(False)
         .astype(bool)
    )

base = pd.DataFrame({
    "gene": df["gene"].astype(str).str.lower(),
    "pattern": parse_bool_series(df["pattern_agreement"]).astype(float),
    "range": pct(df["mean_dynamic_range"]),
    "consistency": pct(df["temporal_effect_consistency"]),
    "recovery": pct(df["mean_recovery"]),
    "peak_agreement": pd.to_numeric(
        df["peak_agreement_rank"], errors="coerce"
    ).fillna(0),
})

scenarios = {
    "S0_predeclared": {"pattern": .30, "range": .20, "consistency": .20, "recovery": .15, "peak_agreement": .15},
    "S1_range_heavy": {"pattern": .30, "range": .30, "consistency": .20, "recovery": .10, "peak_agreement": .10},
    "S2_consistency_heavy": {"pattern": .30, "range": .15, "consistency": .30, "recovery": .15, "peak_agreement": .10},
    "S3_recovery_heavy": {"pattern": .30, "range": .15, "consistency": .15, "recovery": .30, "peak_agreement": .10},
}

rows = []
for name, weights in scenarios.items():
    score = sum(base[k] * w for k, w in weights.items())
    tmp = pd.DataFrame({"gene": base["gene"], "scenario": name, "score": score})
    tmp["rank"] = tmp["score"].rank(ascending=False, method="min")
    rows.append(tmp)

res = pd.concat(rows, ignore_index=True)
res.to_csv(OUT / "SENSITIVITY_ALL_RESULTS.csv", index=False)

k = min(20, len(base))
ref = set(res[res["scenario"] == "S0_predeclared"].nsmallest(k, "rank")["gene"])
overlap = []
for name in scenarios:
    cur = set(res[res["scenario"] == name].nsmallest(k, "rank")["gene"])
    inter = len(ref & cur)
    overlap.append({
        "scenario": name,
        "top20_overlap_with_S0": inter,
        "top20_overlap_fraction": inter / max(1, k),
    })
pd.DataFrame(overlap).to_csv(
    OUT / "SENSITIVITY_TOP20_STABILITY.csv", index=False
)

pivot = res.pivot(index="gene", columns="scenario", values="rank")
rank_rows = []
for gene, row in pivot.iterrows():
    vals = row.dropna()
    rank_rows.append({
        "gene": gene,
        "mean_rank": vals.mean(),
        "rank_sd": vals.std(ddof=1) if len(vals) > 1 else 0.0,
        "rank_range": vals.max() - vals.min(),
    })
pd.DataFrame(rank_rows).sort_values(["mean_rank", "rank_sd"]).to_csv(
    OUT / "SENSITIVITY_GENE_RANK_STABILITY.csv", index=False
)

base_score = res[res["scenario"] == "S0_predeclared"].set_index("gene")["score"]
corr = []
for name in scenarios:
    s = res[res["scenario"] == name].set_index("gene")["score"]
    common = base_score.index.intersection(s.index)
    rho, _ = spearmanr(
        base_score.loc[common], s.loc[common]
    ) if len(common) >= 3 else (np.nan, np.nan)
    corr.append({"scenario": name, "spearman_score_correlation_with_S0": rho})
pd.DataFrame(corr).to_csv(
    OUT / "SENSITIVITY_SCORE_CORRELATIONS.csv", index=False
)

print("=" * 78)
print("28 — TEMPORAL SENSITIVITY ANALYSIS COMPLETE")
print(f"Output: {OUT}")
