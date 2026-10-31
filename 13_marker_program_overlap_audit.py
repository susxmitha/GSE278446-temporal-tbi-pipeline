import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from config import *

# ============================================================
# 13 — MARKER / PROGRAM OVERLAP AUDIT
# ============================================================
# Marker overlap is descriptive. Hypergeometric tests use the
# expressed-gene universe actually present in this dataset.
# ============================================================

strict_file = REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv"
gene_file = RESULT_ROOT / "GSE278446_gene_level_sample_data.csv"
OUT = MARKER_ROOT
OUT.mkdir(parents=True, exist_ok=True)

strict = pd.read_csv(strict_file)
strict_genes = set(strict["gene"].astype(str).str.strip().str.lower())
all_genes = set(pd.read_csv(gene_file, usecols=["gene"])["gene"].astype(str).str.strip().str.lower().unique())

rows=[]
for program, markers in MARKER_PROGRAMS.items():
    marker_set = set(m.lower() for m in markers)
    marker_bg = marker_set & all_genes
    overlap = strict_genes & marker_bg
    N = len(all_genes)
    K = len(marker_bg)
    n = len(strict_genes & all_genes)
    k = len(overlap)
    pval = hypergeom.sf(k-1, N, K, n) if N and K and n else np.nan
    rows.append({
        "program": program,
        "background_genes": N,
        "marker_genes_in_background": K,
        "strict_genes_in_background": n,
        "overlap_count": k,
        "overlap_fraction_of_strict": k/max(n,1),
        "hypergeometric_p": pval,
        "overlap_genes": ";".join(sorted(overlap)),
    })

out_df=pd.DataFrame(rows)
mask=out_df["hypergeometric_p"].notna()
out_df.loc[mask,"hypergeometric_q"]=multipletests(out_df.loc[mask,"hypergeometric_p"],method="fdr_bh")[1]
out_df.to_csv(OUT / "STRICT_GENE_MARKER_OVERLAP.csv",index=False)

print("=" * 78)
print("13 — MARKER / PROGRAM OVERLAP AUDIT")
print(out_df.to_string(index=False))
print(f"Output: {OUT}")
