import time
import requests
import pandas as pd
from config import *

# ============================================================
# 30 — FINAL PAPER-LEVEL RECHECK
# ============================================================
# Rechecks only the final 5 candidates and produces a claim-safe
# table for manuscript preparation.
# ============================================================

INPUT=FINAL_ROOT/"FINAL_5_STRONGEST_CANDIDATES.csv"
OUT=PAPER_ROOT
OUT.mkdir(parents=True,exist_ok=True)

def count(q):
    p={"db":"pubmed","term":q,"retmode":"json","retmax":20,"tool":NCBI_TOOL}
    if NCBI_EMAIL:p["email"]=NCBI_EMAIL
    if NCBI_API_KEY:p["api_key"]=NCBI_API_KEY
    r=requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",params=p,timeout=30)
    r.raise_for_status()
    return int(r.json()["esearchresult"].get("count",0))

df=pd.read_csv(INPUT)
rows=[]
for gene in df["gene"].astype(str).str.lower():
    q_direct=f'"{gene}"[Title/Abstract] AND zebrafish[Title/Abstract] AND (telencephalon[Title/Abstract] OR brain[Title/Abstract]) AND (injury[Title/Abstract] OR trauma[Title/Abstract] OR regeneration[Title/Abstract])'
    q_exact=f'"{gene}"[Title/Abstract] AND ("Danio rerio"[Title/Abstract]) AND ("traumatic brain injury"[Title/Abstract] OR brain injury[Title/Abstract])'
    try:a=count(q_direct)
    except Exception:a=None
    time.sleep(NCBI_DELAY_SECONDS)
    try:b=count(q_exact)
    except Exception:b=None
    rows.append({
        "gene":gene,
        "final_direct_context_count":a,
        "final_exact_tbi_count":b,
        "novelty_statement_supported":False,
        "required_manual_full_text_review":True,
    })
    time.sleep(NCBI_DELAY_SECONDS)
out=pd.DataFrame(rows)
out.to_csv(OUT/"FINAL_5_PAPER_RECHECK.csv",index=False)

print("="*78)
print("30 — FINAL PAPER RECHECK COMPLETE")
print(out.to_string(index=False))
