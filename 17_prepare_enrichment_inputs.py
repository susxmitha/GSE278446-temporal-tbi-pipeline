import pandas as pd
from config import *

# ============================================================
# 17 — ENRICHMENT INPUT PREPARATION
# ============================================================
# All gene lists are derived from fixed pipeline outputs.
# The enrichment background is the tested gene universe, not the
# entire zebrafish genome.
# ============================================================

OUT=ENRICHMENT_ROOT
OUT.mkdir(parents=True,exist_ok=True)

all_tested=pd.read_csv(REPLICATE_VALIDATION_ROOT/"ALL_AVAILABLE_TRAJECTORY_VALIDATION.csv")
strict=pd.read_csv(REPLICATE_VALIDATION_ROOT/"STRICT_3_OF_3_AVAILABLE_REVERSALS.csv")
top=pd.read_csv(PRIORITIZATION_ROOT/"TOP_50_PRIORITY_GENES.csv")

def clean(x):
    return sorted(set(
        x.astype(str).str.strip().str.upper()
        .replace({"NAN":"","NONE":"","GENE":""})
        .loc[lambda s:s.ne("")]
    ))

strict_genes=clean(strict["gene"])
background=clean(all_tested["gene"])
top_genes=clean(top["gene"])

# Early/24h peak lists from the common-window reversal table.
common=pd.read_csv(COMMON_WINDOW_ROOT/"COMMON_WINDOW_REVERSAL_ALL_GENES.csv")
peak24=clean(common.loc[
    (common["Exp1_pattern"]=="peak_at_24") &
    (common["Exp2_pattern"]=="peak_at_24") &
    common["strong_reproducible"],"gene"
])
dip24=clean(common.loc[
    (common["Exp1_pattern"]=="dip_at_24") &
    (common["Exp2_pattern"]=="dip_at_24") &
    common["strong_reproducible"],"gene"
])

lists={
    "STRICT_3OF3_AVAILABLE":strict_genes,
    "TOP50":top_genes,
    "REPRODUCIBLE_24H_PEAK":peak24,
    "REPRODUCIBLE_24H_DIP":dip24,
}
for name,genes in lists.items():
    pd.DataFrame({"gene":genes}).to_csv(OUT/f"{name}.csv",index=False)
    (OUT/f"{name}.txt").write_text("\n".join(genes)+"\n",encoding="utf-8")

pd.DataFrame({"gene":background}).to_csv(OUT/"BACKGROUND_GENES.csv",index=False)
(OUT/"BACKGROUND_GENES.txt").write_text("\n".join(background)+"\n",encoding="utf-8")

summary=pd.DataFrame([
    {"list":name,"n_genes":len(genes),"background_n":len(background)}
    for name,genes in lists.items()
])
summary.to_csv(OUT/"ENRICHMENT_INPUT_SUMMARY.csv",index=False)

print("="*78)
print("17 — ENRICHMENT INPUTS READY")
print(summary.to_string(index=False))
print(f"Background genes: {len(background)}")
print(f"Output: {OUT}")
