from config import *

import gzip
import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


# ============================================================
# 06 — SINGLE-CELL PREPROCESSING + SAMPLE GENE MEANS
# ============================================================
#
# Canonical preprocessing:
#
# 1. Gene-expression features only when feature_type is available.
# 2. Cell QC from raw counts.
# 3. Library-size normalization to 10,000 counts/cell.
# 4. Gene means = mean normalized expression across QC-passing cells.
# 5. Cell program scores = mean log1p(normalized expression)
#    across available markers.
# 6. Duplicate gene symbols are aggregated, never overwritten.
#
# IMPORTANT DATA-PROVENANCE NOTE
# --------------------------------
# GSM8546664 and GSM8546673 have been independently confirmed
# as byte-for-byte identical in the original GSE278446_RAW.tar.
#
# They are therefore retained in the raw/sample inventory, but
# explicitly marked as SOURCE_LEVEL_DUPLICATES.
#
# They MUST NOT be interpreted as two independent biological
# observations during downstream statistical inference.
#
# This script does NOT remove, modify, relabel, or overwrite
# either sample.
# ============================================================


# ============================================================
# CONFIRMED SOURCE-LEVEL DUPLICATE PAIR
# ============================================================

SOURCE_LEVEL_DUPLICATE_PAIR = {
    "GSM8546664",
    "GSM8546673",
}


# ============================================================
# SAMPLE INDEPENDENCE STATUS
# ============================================================

def biological_independence_status(gsm):

    if gsm in SOURCE_LEVEL_DUPLICATE_PAIR:

        return "SOURCE_LEVEL_DUPLICATE"

    return "INDEPENDENT_RAW_SAMPLE"


# ============================================================
# LOCATE RAW FILES
# ============================================================

def locate_files(gsm):

    result = {}

    for p in RAW_ROOT.rglob("*"):

        if not p.is_file():
            continue

        if gsm.lower() not in str(p).lower():
            continue

        n = p.name.lower()

        if (
            "matrix" in n
            and n.endswith((".mtx", ".mtx.gz"))
        ):

            result["matrix"] = p

        elif (
            "barcodes" in n
            and n.endswith(
                (".tsv.gz", ".tsv", ".txt.gz", ".txt")
            )
        ):

            result["barcodes"] = p

        elif (
            "features" in n
            and n.endswith(
                (".tsv.gz", ".tsv", ".txt.gz", ".txt")
            )
        ):

            result["features"] = p

    if set(result) != {
        "matrix",
        "barcodes",
        "features"
    }:

        raise FileNotFoundError(
            f"{gsm}: incomplete raw files: {result}"
        )

    return result


# ============================================================
# READ TSV
# ============================================================

def read_tsv(path):

    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        compression=(
            "gzip"
            if str(path).lower().endswith(".gz")
            else None
        ),
        dtype=str,
    )


# ============================================================
# LOAD + PREPROCESS ONE SAMPLE
# ============================================================

def load_sample(gsm):

    f = locate_files(gsm)

    features = read_tsv(
        f["features"]
    )

    barcodes = read_tsv(
        f["barcodes"]
    )

    matrix = scipy.io.mmread(
        f["matrix"]
    ).tocsr()


    # --------------------------------------------------------
    # DIMENSION CHECKS
    # --------------------------------------------------------

    if matrix.shape[1] != len(barcodes):

        raise RuntimeError(
            f"{gsm}: matrix/barcode mismatch"
        )


    if matrix.shape[0] != len(features):

        raise RuntimeError(
            f"{gsm}: matrix/feature mismatch"
        )


    # --------------------------------------------------------
    # FEATURE TYPE FILTER
    # --------------------------------------------------------

    if features.shape[1] >= 3:

        gene_mask = (
            features.iloc[:, 2]
            .astype(str)
            .str.lower()
            .eq("gene expression")
            .to_numpy()
        )

    else:

        gene_mask = np.ones(
            len(features),
            dtype=bool
        )


    features = (
        features
        .loc[gene_mask]
        .reset_index(drop=True)
    )

    matrix = (
        matrix[
            gene_mask,
            :
        ]
        .tocsr()
    )


    # --------------------------------------------------------
    # GENE NAMES
    # --------------------------------------------------------

    gene_names = (

        features.iloc[:, 1]
        if features.shape[1] >= 2
        else features.iloc[:, 0]

    ).astype(str).str.strip().str.lower().to_numpy()


    # --------------------------------------------------------
    # FALLBACK TO GENE ID FOR BLANK NAMES
    # --------------------------------------------------------

    if features.shape[1] >= 1:

        gene_ids = (
            features.iloc[:, 0]
            .astype(str)
            .str.strip()
            .str.lower()
            .to_numpy()
        )

        blank = np.isin(
            gene_names,
            ["", "nan", "none"]
        )

        gene_names[blank] = gene_ids[blank]


    # --------------------------------------------------------
    # RAW CELL QC METRICS
    # --------------------------------------------------------

    total_counts = (
        np.asarray(
            matrix.sum(axis=0)
        )
        .ravel()
        .astype(float)
    )


    genes_detected = (
        np.asarray(
            (matrix > 0).sum(axis=0)
        )
        .ravel()
        .astype(int)
    )


    # --------------------------------------------------------
    # MITOCHONDRIAL CONTENT
    # --------------------------------------------------------

    mt_mask = np.array(

        [
            (
                g.startswith("mt-")
                or g.startswith("mt:")
                or g.startswith("mt.")
            )

            for g in gene_names
        ],

        dtype=bool
    )


    if mt_mask.any():

        mito_counts = (
            np.asarray(
                matrix[
                    mt_mask,
                    :
                ].sum(axis=0)
            )
            .ravel()
        )

    else:

        mito_counts = np.zeros(
            matrix.shape[1]
        )


    pct_mito = (
        100.0
        * mito_counts
        / np.maximum(
            total_counts,
            1.0
        )
    )


    # --------------------------------------------------------
    # CELL QC FILTER
    # --------------------------------------------------------

    keep = (

        (total_counts >= MIN_TOTAL_COUNTS)

        &

        (total_counts <= MAX_TOTAL_COUNTS)

        &

        (genes_detected >= MIN_GENES_DETECTED)

        &

        (genes_detected <= MAX_GENES_DETECTED)

        &

        (pct_mito <= MAX_MITO_PERCENT)

    )


    if keep.sum() == 0:

        raise RuntimeError(
            f"{gsm}: zero cells pass QC."
        )


    # --------------------------------------------------------
    # QC-PASSING COUNTS
    # --------------------------------------------------------

    counts = (
        matrix[
            :,
            keep
        ]
        .tocsr()
    )


    kept_total = (
        total_counts[
            keep
        ]
    )


    # --------------------------------------------------------
    # LIBRARY-SIZE NORMALIZATION
    # --------------------------------------------------------

    scale = (
        NORMALIZATION_TARGET
        /
        np.maximum(
            kept_total,
            1.0
        )
    )


    normalized = (
        counts
        .multiply(
            scale[
                np.newaxis,
                :
            ]
        )
        .tocsr()
    )


    # --------------------------------------------------------
    # AGGREGATE DUPLICATE GENE SYMBOLS
    # --------------------------------------------------------

    codes, unique_genes = (
        pd.factorize(
            gene_names,
            sort=False
        )
    )


    group = sp.csr_matrix(

        (
            np.ones(
                len(codes),
                dtype=np.float32
            ),

            (
                codes,
                np.arange(
                    len(codes)
                )
            )
        ),

        shape=(
            len(unique_genes),
            len(codes)
        )
    )


    gene_norm = (
        group
        @ normalized
    ).tocsr()


    gene_means = (
        np.asarray(
            gene_norm.mean(
                axis=1
            )
        )
        .ravel()
    )


    # ========================================================
    # GENE-LEVEL SAMPLE TABLE
    # ========================================================

    gene_table = pd.DataFrame({

        "gene":
            pd.Series(
                unique_genes,
                dtype=str
            ),

        "mean_expression":
            gene_means.astype(float),

        "mean_expression_log1p":
            np.log1p(
                gene_means.astype(float)
            ),

        "GSM":
            gsm,

        "condition":
            SAMPLE_METADATA[gsm]["condition"],

        "time_hpi":
            SAMPLE_METADATA[gsm]["time_hpi"],

        "experiment":
            SAMPLE_METADATA[gsm]["experiment"],

        "replicate":
            SAMPLE_METADATA[gsm]["replicate"],

        "biological_independence_status":
            biological_independence_status(
                gsm
            ),

    })


    # ========================================================
    # CELL PROGRAM SCORES
    # ========================================================

    log_norm = (
        normalized.copy()
    )


    log_norm.data = np.log1p(
        log_norm.data
    )


    gene_index = {
        g: i
        for i, g
        in enumerate(unique_genes)
    }


    gene_log = (
        group
        @ log_norm
    ).tocsr()


    # ========================================================
    # CELL TABLE
    # ========================================================

    cell_table = pd.DataFrame({

        "cell_index":
            np.arange(
                keep.sum()
            ),

        "GSM":
            gsm,

        "condition":
            SAMPLE_METADATA[gsm]["condition"],

        "time_hpi":
            SAMPLE_METADATA[gsm]["time_hpi"],

        "experiment":
            SAMPLE_METADATA[gsm]["experiment"],

        "replicate":
            SAMPLE_METADATA[gsm]["replicate"],

        "biological_independence_status":
            biological_independence_status(
                gsm
            ),

        "total_counts":
            total_counts[keep],

        "genes_detected":
            genes_detected[keep],

        "pct_mito":
            pct_mito[keep],

    })


    # --------------------------------------------------------
    # PROGRAM SCORES
    # --------------------------------------------------------

    for program, markers in (
        MARKER_PROGRAMS.items()
    ):

        idx = [

            gene_index[m.lower()]

            for m in markers

            if m.lower()
            in gene_index

        ]


        if idx:

            score = np.asarray(

                gene_log[
                    idx,
                    :
                ].mean(
                    axis=0
                )

            ).ravel()

            cell_table[
                f"{program}_score"
            ] = score

        else:

            cell_table[
                f"{program}_score"
            ] = np.nan


    # ========================================================
    # QC SUMMARY
    # ========================================================

    qc = {

        "raw_cells":
            int(
                matrix.shape[1]
            ),

        "qc_cells":
            int(
                keep.sum()
            ),

        "raw_genes":
            int(
                len(gene_names)
            ),

        "gene_expression_features":
            int(
                len(features)
            ),

        "unique_genes":
            int(
                len(unique_genes)
            ),

        "mean_total_counts":
            float(
                np.mean(
                    total_counts
                )
            ),

        "median_total_counts":
            float(
                np.median(
                    total_counts
                )
            ),

        "mean_genes_detected":
            float(
                np.mean(
                    genes_detected
                )
            ),

        "median_genes_detected":
            float(
                np.median(
                    genes_detected
                )
            ),

        "mean_pct_mito":
            float(
                np.mean(
                    pct_mito
                )
            ),

        "qc_fraction":
            float(
                keep.mean()
            ),

        "biological_independence_status":
            biological_independence_status(
                gsm
            ),

    }


    return (
        gene_table,
        cell_table,
        qc
    )


# ============================================================
# MAIN PROCESSING
# ============================================================

print("=" * 78)
print("06 — SINGLE-CELL PREPROCESSING + SAMPLE GENE MEANS")
print("=" * 78)

print()
print(
    "IMPORTANT: GSM8546664 and GSM8546673 are a "
    "confirmed source-level duplicate pair."
)

print(
    "They will be retained but flagged as "
    "SOURCE_LEVEL_DUPLICATE."
)

print(
    "They must not be treated as independent biological "
    "observations in downstream inference."
)

print()


gene_tables = []
cell_tables = []
qc_rows = []


# ============================================================
# PROCESS ALL EXPECTED GSMs
# ============================================================

for n, gsm in enumerate(
    ANALYSIS_GSMS,
    1
):

    print(
        f"[{n:02d}/{len(ANALYSIS_GSMS)}] {gsm}"
    )


    gene_df, cell_df, qc = (
        load_sample(
            gsm
        )
    )


    gene_out = (
        RESULT_ROOT
        / "gene_means"
    )


    cell_out = (
        RESULT_ROOT
        / "cell_program_scores"
    )


    gene_out.mkdir(
        exist_ok=True
    )


    cell_out.mkdir(
        exist_ok=True
    )


    gene_df.to_csv(

        gene_out
        /
        f"{gsm}_gene_means.csv",

        index=False

    )


    cell_df.to_csv(

        cell_out
        /
        f"{gsm}_cell_program_scores.csv",

        index=False

    )


    gene_tables.append(
        gene_df
    )

    cell_tables.append(
        cell_df
    )


    qc_rows.append({

        "GSM":
            gsm,

        **SAMPLE_METADATA[gsm],

        **qc

    })


# ============================================================
# COMBINE RESULTS
# ============================================================

gene_all = pd.concat(
    gene_tables,
    ignore_index=True
)


cell_all = pd.concat(
    cell_tables,
    ignore_index=True
)


qc_df = pd.DataFrame(
    qc_rows
)


# ============================================================
# SAVE MAIN OUTPUTS
# ============================================================

gene_all.to_csv(

    RESULT_ROOT
    /
    "GSE278446_gene_level_sample_data.csv",

    index=False

)


cell_all.to_csv(

    RESULT_ROOT
    /
    "GSE278446_cell_program_scores.csv",

    index=False

)


qc_df.to_csv(

    RESULT_ROOT
    /
    "GSE278446_cell_qc_summary.csv",

    index=False

)


# ============================================================
# SAMPLE-LEVEL PROGRAM SUMMARY
# ============================================================
#
# IMPORTANT:
# The sample is the unit for temporal inference.
#
# Individual cells are NOT treated as independent biological
# replicates.
# ============================================================

program_cols = [
    f"{p}_score"
    for p in MARKER_PROGRAMS
]


agg = []


for keys, g in cell_all.groupby(

    [
        "GSM",
        "condition",
        "time_hpi",
        "experiment",
        "replicate",
        "biological_independence_status"
    ],

    sort=False

):

    row = dict(

        zip(

            [
                "GSM",
                "condition",
                "time_hpi",
                "experiment",
                "replicate",
                "biological_independence_status"
            ],

            keys

        )

    )


    row["cells"] = len(g)


    for p in MARKER_PROGRAMS:

        c = f"{p}_score"


        row[
            f"{p}_mean"
        ] = g[c].mean()


        row[
            f"{p}_median"
        ] = g[c].median()


    agg.append(
        row
    )


sample_program = (

    pd.DataFrame(
        agg
    )

    .sort_values(
        [
            "experiment",
            "time_hpi",
            "replicate"
        ]
    )

)


sample_program.to_csv(

    RESULT_ROOT
    /
    "GSE278446_sample_program_summary.csv",

    index=False

)


# ============================================================
# DUPLICATE SAMPLE REPORT
# ============================================================

duplicate_report = pd.DataFrame({

    "GSM": [
        "GSM8546664",
        "GSM8546673"
    ],

    "condition": [
        SAMPLE_METADATA[
            "GSM8546664"
        ]["condition"],

        SAMPLE_METADATA[
            "GSM8546673"
        ]["condition"]
    ],

    "time_hpi": [
        SAMPLE_METADATA[
            "GSM8546664"
        ]["time_hpi"],

        SAMPLE_METADATA[
            "GSM8546673"
        ]["time_hpi"]
    ],

    "experiment": [
        SAMPLE_METADATA[
            "GSM8546664"
        ]["experiment"],

        SAMPLE_METADATA[
            "GSM8546673"
        ]["experiment"]
    ],

    "replicate": [
        SAMPLE_METADATA[
            "GSM8546664"
        ]["replicate"],

        SAMPLE_METADATA[
            "GSM8546673"
        ]["replicate"]
    ],

    "biological_independence_status": [
        "SOURCE_LEVEL_DUPLICATE",
        "SOURCE_LEVEL_DUPLICATE"
    ],

    "source_verification": [
        "Confirmed against original GSE278446_RAW.tar",
        "Confirmed against original GSE278446_RAW.tar"
    ]

})


duplicate_report.to_csv(

    RESULT_ROOT
    /
    "06_source_level_duplicate_samples.csv",

    index=False

)


# ============================================================
# FINAL OUTPUT
# ============================================================

print()
print("=" * 78)
print("06 — PREPROCESSING COMPLETE")
print("=" * 78)

print(
    f"Samples processed: {len(qc_df)}"
)

print(
    f"Total QC-passing cells: {len(cell_all):,}"
)

print(
    f"Unique gene/sample rows: {len(gene_all):,}"
)

print()

print(
    "SOURCE-LEVEL DUPLICATE SAMPLES:"
)

print(
    "  GSM8546664 (retained source representative)"
)

print(
    "  GSM8546673 (excluded source duplicate)"
)

print()

print(
    "These samples were retained for provenance and "
    "preprocessing but flagged as non-independent."
)

print()

print(
    f"Results: {RESULT_ROOT}"
)
