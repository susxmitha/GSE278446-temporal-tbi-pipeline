import numpy as np
import pandas as pd
from config import *

# ============================================================
# 16 — GENE PRIORITIZATION
# ============================================================
# No circular reuse of "in existing TOP50". The ranking is based
# on temporal evidence plus independent pathway-level evidence.
# ============================================================

OUT=PRIORITIZATION_ROOT
OUT.mkdir(parents=True,exist_ok=True)

genes=pd.read_csv(GENE_CHARACTERIZATION_ROOT/"ALL_STRICT_GENE_CHARACTERIZATION.csv")
pathway_file=PATHWAY_INTERPRETATION_ROOT/"MASTER_PATHWAY_INTERPRETATION.csv"
pathways=pd.read_csv(pathway_file) if pathway_file.exists() else pd.DataFrame()

# Temporal evidence is already characterized; normalize only within
# this fixed strict candidate universe.
def pct(s):
    s=pd.to_numeric(s,errors="coerce")
    return s.rank(pct=True,method="average").fillna(0)

genes["temporal_score"] = (
    0.30*genes["pattern_agreement"].astype(float)
    +0.20*pct(genes["mean_dynamic_range"])
    +0.20*pct(genes["temporal_effect_consistency"])
    +0.15*pct(genes["mean_recovery"])
    +0.15*genes["peak_agreement_rank"].fillna(0)
)

# Pathway evidence is list-level, not a fabricated gene-pathway assignment.
# Therefore it is reported separately and NOT used as a gene-specific score
# unless the enrichment service actually returns an intersection containing
# that gene.
# Stage-19 pathway enrichment is list-level evidence. The returned
# intersection identifiers are not assumed to be the original gene
# symbols, so no gene-specific pathway score is fabricated.
genes["pathway_gene_evidence"] = 0
genes["pathway_terms_supporting_gene"] = 0
genes["pathway_evidence_rank"] = 0.0

genes["final_priority_score"] = genes["temporal_score"]

genes["priority_class"]=np.select(
    [
        (genes["pattern_agreement"]) & (genes["temporal_score"]>=0.80),
        (genes["pattern_agreement"]) & (genes["temporal_score"]>=0.60),
        genes["pattern_agreement"],
    ],
    ["HIGH_TEMPORAL_PRIORITY","STRONG_TEMPORAL_PRIORITY","REPRODUCIBLE"],
    default="REVIEW",
)

genes=genes.sort_values(["final_priority_score","gene"],ascending=[False,True]).reset_index(drop=True)
genes["rank"]=np.arange(1,len(genes)+1)
genes.to_csv(OUT/"ALL_STRICT_GENE_PRIORITIZATION.csv",index=False)
genes.head(50).to_csv(OUT/"TOP_50_PRIORITY_GENES.csv",index=False)
genes.head(20).to_csv(OUT/"TOP_20_PRIORITY_GENES.csv",index=False)

print("="*78)
print("16 — GENE PRIORITIZATION COMPLETE")
print(f"Genes ranked: {len(genes)}")
print(f"Top candidate: {genes.iloc[0]['gene'] if len(genes) else 'NONE'}")
print(f"Output: {OUT}")
