from pathlib import Path
import pandas as pd
import traceback

from config import (
    RESEARCH_ROOT,
    COMMON_WINDOW_ROOT,
    ENRICHMENT_RESULTS_ROOT,
    PATHWAY_INTERPRETATION_ROOT,
)

# ============================================================
# 19 — PATHWAY REPRODUCIBILITY INTERPRETATION
# ============================================================
# Reads the four Stage-18 enrichment result files and reports
# which pathway terms are reproduced across input gene sets.
#
# Contract with Stage 18:
#   Stage 18 writes columns:
#     gene_set, source, term_id, term_name, description, p_value,
#     intersection_size, term_size, query_size, effective_domain_size,
#     precision, recall, significant, intersection_genes, status
#   Stage 18 writes files named:
#     {gene_set}_enrichment.csv
#
# This stage does NOT fabricate biological interpretation.
# It only counts observed g:Profiler term IDs.
# ============================================================


# ============================================================
# PATHS
# ============================================================

COMMON_ROOT = COMMON_WINDOW_ROOT
ENRICHMENT_ROOT = ENRICHMENT_RESULTS_ROOT
OUTPUT_ROOT = PATHWAY_INTERPRETATION_ROOT
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


# ============================================================
# EXACT STAGE-18 INPUT FILENAMES
# ============================================================

INPUT_FILES = {
    "STRICT_3OF3_AVAILABLE":
        ENRICHMENT_ROOT / "STRICT_3OF3_AVAILABLE_enrichment.csv",

    "TOP50":
        ENRICHMENT_ROOT / "TOP50_enrichment.csv",

    "REPRODUCIBLE_24H_PEAK":
        ENRICHMENT_ROOT / "REPRODUCIBLE_24H_PEAK_enrichment.csv",

    "REPRODUCIBLE_24H_DIP":
        ENRICHMENT_ROOT / "REPRODUCIBLE_24H_DIP_enrichment.csv",
}


# ============================================================
# OUTPUT FILES
# ============================================================

MASTER_FILE = OUTPUT_ROOT / "MASTER_PATHWAY_INTERPRETATION.csv"
REPRODUCIBILITY_FILE = OUTPUT_ROOT / "PATHWAY_REPRODUCIBILITY_SUMMARY.csv"
STATUS_FILE = OUTPUT_ROOT / "PATHWAY_INTERPRETATION_STATUS.csv"


# ============================================================
# REQUIRED COLUMNS FROM STAGE 18
# ============================================================
# NOTE: Stage 18 emits "gene_set", not "input_set".
# This stage renames it internally to "input_set" for clarity,
# but the required-column check must use the real Stage-18 name.
# ============================================================

REQUIRED_COLUMNS = [
    "gene_set",
    "source",
    "term_id",
    "term_name",
    "description",
    "p_value",
    "intersection_size",
    "term_size",
    "query_size",
    "effective_domain_size",
    "precision",
    "recall",
    "significant",
    "intersection_genes",
]


# ============================================================
# HEADER
# ============================================================

print("=" * 90)
print("19 — PATHWAY REPRODUCIBILITY INTERPRETATION")
print("=" * 90)

print(f"Research root    : {RESEARCH_ROOT}")
print(f"Common root      : {COMMON_ROOT}")
print(f"Enrichment root  : {ENRICHMENT_ROOT}")
print(f"Output directory : {OUTPUT_ROOT}")


# ============================================================
# READ ONE STAGE-18 FILE
# ============================================================

def read_enrichment_file(gene_set, path):

    print("")
    print("-" * 80)
    print(f"Reading: {gene_set}")
    print(f"File   : {path}")

    if not path.exists():
        raise FileNotFoundError(
            f"Required Stage-18 enrichment file not found:\n{path}"
        )

    df = pd.read_csv(path, dtype=str)

    print(f"Rows loaded: {len(df):,}")
    print(f"Columns    : {list(df.columns)}")

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"{path.name} is missing required columns: {missing}"
        )

    df = df[REQUIRED_COLUMNS].copy()

    # Force the authoritative set name so a stale label in the
    # file cannot misclassify the row.
    df["gene_set"] = gene_set

    # Drop records that never received a pathway ID.
    df["term_id"] = (
        df["term_id"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    df = df[df["term_id"] != ""].copy()

    print(f"Valid pathway records: {len(df):,}")

    return df


# ============================================================
# MAIN
# ============================================================

def main():

    frames = []
    status_rows = []

    for gene_set, path in INPUT_FILES.items():

        try:
            df = read_enrichment_file(gene_set, path)
            frames.append(df)
            status_rows.append({
                "gene_set": gene_set,
                "file": str(path),
                "rows_loaded": len(df),
                "status": "LOADED",
            })
        except Exception as exc:
            print("")
            print(f"ERROR: {gene_set}")
            print(str(exc))
            status_rows.append({
                "gene_set": gene_set,
                "file": str(path),
                "rows_loaded": 0,
                "status": "ERROR",
            })

    # --------------------------------------------------------
    # Combine
    # --------------------------------------------------------

    if frames:
        master = pd.concat(frames, ignore_index=True)
    else:
        master = pd.DataFrame(columns=REQUIRED_COLUMNS)

    # Rename to the canonical internal name used downstream.
    master = master.rename(columns={"gene_set": "input_set"})

    # --------------------------------------------------------
    # Numeric p-value and sort
    # --------------------------------------------------------

    if not master.empty:
        master["p_value_numeric"] = pd.to_numeric(
            master["p_value"], errors="coerce"
        )
        master = master.sort_values(
            ["input_set", "p_value_numeric"],
            ascending=[True, True],
            na_position="last",
        )

    master.to_csv(MASTER_FILE, index=False)

    # --------------------------------------------------------
    # Pathway reproducibility
    # --------------------------------------------------------

    if not master.empty:

        occurrence_base = (
            master[["term_id", "term_name", "source", "input_set"]]
            .drop_duplicates(subset=["term_id", "input_set"])
        )

        occurrence = (
            occurrence_base
            .groupby(
                ["term_id", "term_name", "source"],
                dropna=False,
                as_index=False,
            )
            .agg(
                input_sets=(
                    "input_set",
                    lambda x: ";".join(
                        sorted(
                            set(
                                str(v)
                                for v in x
                                if str(v).strip()
                            )
                        )
                    ),
                )
            )
        )

        occurrence["n_input_sets"] = (
            occurrence["input_sets"]
            .str.split(";")
            .apply(
                lambda values: len([v for v in values if v])
            )
        )

        occurrence = occurrence.sort_values(
            ["n_input_sets", "source", "term_id"],
            ascending=[False, True, True],
        )

    else:
        occurrence = pd.DataFrame(
            columns=[
                "term_id",
                "term_name",
                "source",
                "input_sets",
                "n_input_sets",
            ]
        )

    occurrence.to_csv(REPRODUCIBILITY_FILE, index=False)

    # --------------------------------------------------------
    # Status
    # --------------------------------------------------------

    status_df = pd.DataFrame(status_rows)
    status_df.to_csv(STATUS_FILE, index=False)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    actual_rows = len(master)

    unique_terms = (
        int(master["term_id"].nunique()) if not master.empty else 0
    )

    if not occurrence.empty:
        reproducible_terms = int(
            (occurrence["n_input_sets"] >= 2).sum()
        )
        three_of_four_terms = int(
            (occurrence["n_input_sets"] >= 3).sum()
        )
        four_of_four_terms = int(
            (occurrence["n_input_sets"] == 4).sum()
        )
    else:
        reproducible_terms = 0
        three_of_four_terms = 0
        four_of_four_terms = 0

    print("")
    print("=" * 90)
    print("19 — PATHWAY REPRODUCIBILITY COMPLETE")
    print("=" * 90)

    print(f"Enrichment files loaded : {len(frames)}")
    print(f"Actual enrichment rows  : {actual_rows:,}")
    print(f"Unique pathway terms    : {unique_terms:,}")
    print(f"Terms in >=2 gene sets  : {reproducible_terms:,}")
    print(f"Terms in >=3 gene sets  : {three_of_four_terms:,}")
    print(f"Terms in all 4 sets     : {four_of_four_terms:,}")
    print("")
    print(f"Master output          : {MASTER_FILE}")
    print(f"Reproducibility output : {REPRODUCIBILITY_FILE}")
    print(f"Status output          : {STATUS_FILE}")

    if actual_rows == 0:
        print("")
        print("WARNING: No actual enrichment results were available.")
    else:
        print("")
        print("PASS — all available Stage-18 enrichment results were "
              "loaded successfully.")
        print("Pathway reproducibility was calculated from observed "
              "g:Profiler term IDs.")
        print("No biological interpretation was fabricated.")


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
