import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from config import *

# ============================================================
# 09 — PROGRAM-LEVEL REPLICATE VALIDATION
# ============================================================

INPUT = RESULT_ROOT / "GSE278446_sample_program_summary.csv"
OUT = DISCOVERY_ROOT / "program_replicate_validation"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
rows = []

# Within each experiment: correlation between replicate trajectories.
for exp in EXPECTED_EXPERIMENTS:
    sub = df[df["experiment"] == exp].copy()
    for program in MARKER_PROGRAMS:
        c = f"{program}_mean"
        piv = sub.pivot(index="time_hpi", columns="replicate", values=c)
        if not {1, 2}.issubset(piv.columns):
            continue
        piv = piv[[1, 2]].dropna()
        if len(piv) >= 4 and piv[1].nunique() > 1 and piv[2].nunique() > 1:
            rho, p = spearmanr(piv[1], piv[2])
        else:
            rho, p = np.nan, np.nan
        rows.append({
            "experiment": exp, "program": program,
            "n_timepoints": len(piv),
            "replicate_spearman_rho": rho,
            "replicate_p": p,
        })

rep = pd.DataFrame(rows)
rep.to_csv(OUT / "PROGRAM_REPLICATE_CORRELATIONS.csv", index=False)

# Common-window direction agreement across replicate trajectories.
common_rows = []
for program in MARKER_PROGRAMS:
    for exp in EXPECTED_EXPERIMENTS:
        sub = df[(df["experiment"] == exp) & df["time_hpi"].isin(COMMON_TIMEPOINTS)]
        piv = sub.pivot(index="time_hpi", columns="replicate", values=f"{program}_mean").reindex(COMMON_TIMEPOINTS)
        if not {1,2}.issubset(piv.columns) or piv.isna().any().any():
            continue
        d1 = np.diff(piv[1].to_numpy(float))
        d2 = np.diff(piv[2].to_numpy(float))
        common_rows.append({
            "experiment": exp, "program": program,
            "early_direction_agreement": bool(np.sign(d1[0]) == np.sign(d2[0])),
            "late_direction_agreement": bool(np.sign(d1[1]) == np.sign(d2[1])),
            "both_phase_direction_agreement": bool(np.all(np.sign(d1) == np.sign(d2))),
        })
pd.DataFrame(common_rows).to_csv(OUT / "PROGRAM_COMMON_WINDOW_REPLICATE_AGREEMENT.csv", index=False)

print("=" * 78)
print("09 — PROGRAM REPLICATE VALIDATION COMPLETE")
print(f"Output: {OUT}")
