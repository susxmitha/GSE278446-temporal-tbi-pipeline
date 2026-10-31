from config import *
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from config import *

# ============================================================
# 14 — CELL-LEVEL REVERSAL LOCALIZATION
# ============================================================
# Primary protection against pseudo-replication:
#   * cell scores are retained descriptively;
#   * inferential temporal association is performed on GSM-level
#     sample means, not treating thousands of cells as independent.
# Normalization is CP10K followed by log1p.
# ============================================================

STRICT = pd.read_csv(REPLICATE_VALIDATION_ROOT / "STRICT_3_OF_3_AVAILABLE_REVERSALS.csv")
strict_genes = set(STRICT["gene"].astype(str).str.strip().str.lower())

def locate_files(gsm):
    result={}
    for p in RAW_ROOT.rglob("*"):
        if not p.is_file() or gsm.lower() not in str(p).lower(): continue
        n=p.name.lower()
        if "matrix" in n and n.endswith((".mtx",".mtx.gz")): result["matrix"]=p
        elif "barcodes" in n and n.endswith((".tsv",".tsv.gz",".txt",".txt.gz")): result["barcodes"]=p
        elif "features" in n and n.endswith((".tsv",".tsv.gz",".txt",".txt.gz")): result["features"]=p
    return result

def read_tsv(path):
    return pd.read_csv(path,sep="\t",header=None,
                        compression="gzip" if str(path).lower().endswith(".gz") else None,dtype=str)

cell_tables=[]
for gsm in ANALYSIS_GSMS:
    meta = SAMPLE_METADATA[gsm]
    if meta["time_hpi"] not in COMMON_TIMEPOINTS: continue
    f=locate_files(gsm)
    feat=read_tsv(f["features"])
    mat=scipy.io.mmread(f["matrix"]).tocsr()
    if feat.shape[1]>=3:
        mask=feat.iloc[:,2].astype(str).str.lower().eq("gene expression").to_numpy()
    else:
        mask=np.ones(len(feat),dtype=bool)
    feat=feat.loc[mask].reset_index(drop=True); mat=mat[mask,:].tocsr()
    genes=(feat.iloc[:,1] if feat.shape[1]>=2 else feat.iloc[:,0]).astype(str).str.strip().str.lower().to_numpy()
    total=np.asarray(mat.sum(axis=0)).ravel().astype(float)
    detected=np.asarray((mat>0).sum(axis=0)).ravel()
    mito=np.array([g.startswith("mt-") or g.startswith("mt:") or g.startswith("mt.") for g in genes])
    mito_counts=np.asarray(mat[mito,:].sum(axis=0)).ravel() if mito.any() else np.zeros(mat.shape[1])
    pct_mito=100*mito_counts/np.maximum(total,1)
    keep=(total>=MIN_TOTAL_COUNTS)&(total<=MAX_TOTAL_COUNTS)&(detected>=MIN_GENES_DETECTED)&(detected<=MAX_GENES_DETECTED)&(pct_mito<=MAX_MITO_PERCENT)
    mat=mat[:,keep]
    scale=NORMALIZATION_TARGET/np.maximum(total[keep],1)
    norm=mat.multiply(scale[np.newaxis,:]).tocsr()
    norm.data=np.log1p(norm.data)

    # Aggregate duplicate gene rows, then score strict genes.
    codes, unique=pd.factorize(genes,sort=False)
    group=sp.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(unique),len(codes)))
    agg=(group@norm).tocsr()
    idx={g:i for i,g in enumerate(unique)}
    strict_idx=[idx[g] for g in strict_genes if g in idx]
    if not strict_idx: raise RuntimeError(f"{gsm}: no strict genes found")
    reversal=np.asarray(agg[strict_idx,:].mean(axis=0)).ravel()

    rec=pd.DataFrame({
        "cell_index":np.arange(len(reversal)),
        "GSM":gsm,"experiment":meta["experiment"],"replicate":meta["replicate"],
        "time_hpi":meta["time_hpi"],"condition":meta["condition"],
        "reversal_score":reversal,
    })
    for program,markers in MARKER_PROGRAMS.items():
        idxp=[idx[m.lower()] for m in markers if m.lower() in idx]
        rec[f"{program}_score"]=np.asarray(agg[idxp,:].mean(axis=0)).ravel() if idxp else np.nan
    cell_tables.append(rec)

cells=pd.concat(cell_tables,ignore_index=True)
if set(cells["GSM"].astype(str)) & set(EXCLUDED_GSMS):
    raise RuntimeError("Excluded source duplicate unexpectedly entered cell-level analysis.")
cells.to_csv(CELL_LEVEL_ROOT/"CELL_LEVEL_REVERSAL_SCORES.csv",index=False)

# GSM-level summaries are the primary analysis unit.
sample_rows=[]
for keys,g in cells.groupby(["GSM","experiment","replicate","time_hpi"],sort=False):
    row=dict(zip(["GSM","experiment","replicate","time_hpi"],keys))
    row["n_cells"]=len(g)
    row["reversal_mean"]=g["reversal_score"].mean()
    row["reversal_median"]=g["reversal_score"].median()
    for p in MARKER_PROGRAMS:
        row[f"{p}_mean"]=g[f"{p}_score"].mean()
    sample_rows.append(row)
sample=pd.DataFrame(sample_rows)
sample.to_csv(CELL_LEVEL_ROOT/"SAMPLE_LEVEL_REVERSAL_AND_PROGRAM_SCORES.csv",index=False)

# Temporal summary.
temporal=sample.groupby(["experiment","time_hpi"],as_index=False).agg(
    reversal_mean=("reversal_mean","mean"),
    reversal_sd=("reversal_mean","std"),
    replicate_n=("replicate","nunique"),
)
temporal.to_csv(CELL_LEVEL_ROOT/"TEMPORAL_REVERSAL_SUMMARY.csv",index=False)

# Program associations: sample-level and cell-level are both reported.
rows=[]
for p in MARKER_PROGRAMS:
    c=f"{p}_score"
    x=sample["reversal_mean"].to_numpy(float); y=sample[f"{p}_mean"].to_numpy(float)
    valid=np.isfinite(x)&np.isfinite(y)
    rho,pv=spearmanr(x[valid],y[valid]) if valid.sum()>=5 and np.ptp(x[valid]) and np.ptp(y[valid]) else (np.nan,np.nan)
    xc=cells["reversal_score"].to_numpy(float); yc=cells[c].to_numpy(float)
    v=np.isfinite(xc)&np.isfinite(yc)
    crho,cpv=spearmanr(xc[v],yc[v]) if v.sum()>=20 and np.ptp(xc[v]) and np.ptp(yc[v]) else (np.nan,np.nan)
    rows.append({"program":p,"sample_level_rho":rho,"sample_level_p":pv,"cell_level_rho_exploratory":crho,"cell_level_p_exploratory":cpv})
assoc=pd.DataFrame(rows)
for col in ["sample_level_p","cell_level_p_exploratory"]:
    m=assoc[col].notna()
    if m.any(): assoc.loc[m,col.replace("_p","_q")]=multipletests(assoc.loc[m,col],method="fdr_bh")[1]
assoc.to_csv(CELL_LEVEL_ROOT/"PROGRAM_ASSOCIATION_RESULTS.csv",index=False)

print("="*78)
print("14 — CELL-LEVEL REVERSAL VALIDATION COMPLETE")
print(f"Cells analysed: {len(cells):,}")
print(f"Samples analysed: {len(sample)}")
print(f"Output: {CELL_LEVEL_ROOT}")
