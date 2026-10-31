from pathlib import Path
import os

# ============================================================
# GSE278446 / NeuroNovelty — SINGLE SOURCE OF TRUTH
# ============================================================

_DEFAULT_RESEARCH_ROOT = Path(__file__).resolve().parents[2]
RESEARCH_ROOT = Path(
    os.environ.get(
        "GSE278446_RESEARCH_ROOT",
        str(_DEFAULT_RESEARCH_ROOT),
    )
).expanduser().resolve()

PROJECT_ROOT = RESEARCH_ROOT / "NeuroNovelty"
RAW_ROOT = RESEARCH_ROOT / "GSE278446_RAW"

ANALYSIS_VERSION = "GSE278446_v4_final"
PIPELINE_VERSION = "GSE278446_NeuroNovelty_FINAL_1.0"
RESULT_ROOT = PROJECT_ROOT / "temporal_results" / ANALYSIS_VERSION

DISCOVERY_ROOT = RESULT_ROOT / "discovery"
CORRECTED_DISCOVERY_ROOT = RESULT_ROOT / "corrected_discovery"
COMMON_WINDOW_ROOT = RESULT_ROOT / "common_window"
REPLICATE_VALIDATION_ROOT = COMMON_WINDOW_ROOT / "replicate_validation"
CELL_LEVEL_ROOT = COMMON_WINDOW_ROOT / "cell_level_reversal"
MARKER_ROOT = COMMON_WINDOW_ROOT / "marker_program_audit"
GENE_CHARACTERIZATION_ROOT = COMMON_WINDOW_ROOT / "gene_characterization"
PRIORITIZATION_ROOT = COMMON_WINDOW_ROOT / "gene_prioritization"
ENRICHMENT_ROOT = COMMON_WINDOW_ROOT / "pathway_enrichment"
ENRICHMENT_RESULTS_ROOT = ENRICHMENT_ROOT / "ENRICHMENT_RESULTS"
PATHWAY_INTERPRETATION_ROOT = COMMON_WINDOW_ROOT / "pathway_interpretation"
ALIAS_ROOT = COMMON_WINDOW_ROOT / "alias_audit"
LITERATURE_ROOT = COMMON_WINDOW_ROOT / "literature_audit"
BIOLOGY_ROOT = COMMON_WINDOW_ROOT / "biological_context"
EVIDENCE_ROOT = COMMON_WINDOW_ROOT / "evidence_audit"
SENSITIVITY_ROOT = COMMON_WINDOW_ROOT / "sensitivity_analysis"
FINAL_ROOT = COMMON_WINDOW_ROOT / "final_candidates"
PAPER_ROOT = COMMON_WINDOW_ROOT / "paper_recheck"
AUDIT_ROOT = COMMON_WINDOW_ROOT / "pipeline_audit"
MANUSCRIPT_ROOT = PROJECT_ROOT / "manuscript_handoff"

for _p in [
    PROJECT_ROOT, RESULT_ROOT, DISCOVERY_ROOT, CORRECTED_DISCOVERY_ROOT,
    COMMON_WINDOW_ROOT, REPLICATE_VALIDATION_ROOT, CELL_LEVEL_ROOT,
    MARKER_ROOT, GENE_CHARACTERIZATION_ROOT, PRIORITIZATION_ROOT,
    ENRICHMENT_ROOT, ENRICHMENT_RESULTS_ROOT, PATHWAY_INTERPRETATION_ROOT,
    ALIAS_ROOT, LITERATURE_ROOT, BIOLOGY_ROOT, EVIDENCE_ROOT,
    SENSITIVITY_ROOT, FINAL_ROOT, PAPER_ROOT, AUDIT_ROOT, MANUSCRIPT_ROOT,
]:
    _p.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Original GEO design: 28 GSM records.
# ------------------------------------------------------------
SAMPLE_METADATA = {
    "GSM8546663": {"time_hpi": 0,   "experiment": "Exp1", "replicate": 1, "condition": "undamaged"},
    "GSM8546664": {"time_hpi": 0,   "experiment": "Exp1", "replicate": 2, "condition": "undamaged"},
    "GSM8546665": {"time_hpi": 5,   "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546666": {"time_hpi": 5,   "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546667": {"time_hpi": 12,  "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546668": {"time_hpi": 12,  "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546669": {"time_hpi": 18,  "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546670": {"time_hpi": 18,  "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546671": {"time_hpi": 48,  "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546672": {"time_hpi": 48,  "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546673": {"time_hpi": 72,  "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546674": {"time_hpi": 72,  "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546675": {"time_hpi": 0,   "experiment": "Exp2", "replicate": 1, "condition": "undamaged"},
    "GSM8546676": {"time_hpi": 0,   "experiment": "Exp2", "replicate": 2, "condition": "undamaged"},
    "GSM8546677": {"time_hpi": 72,  "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546678": {"time_hpi": 72,  "experiment": "Exp2", "replicate": 2, "condition": "injured"},
    "GSM8546679": {"time_hpi": 96,  "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546680": {"time_hpi": 96,  "experiment": "Exp2", "replicate": 2, "condition": "injured"},
    "GSM8546681": {"time_hpi": 120, "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546682": {"time_hpi": 120, "experiment": "Exp2", "replicate": 2, "condition": "injured"},
    "GSM8546683": {"time_hpi": 144, "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546684": {"time_hpi": 144, "experiment": "Exp2", "replicate": 2, "condition": "injured"},
    "GSM8546685": {"time_hpi": 168, "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546686": {"time_hpi": 168, "experiment": "Exp2", "replicate": 2, "condition": "injured"},
    "GSM8546687": {"time_hpi": 24,  "experiment": "Exp1", "replicate": 1, "condition": "injured"},
    "GSM8546688": {"time_hpi": 24,  "experiment": "Exp1", "replicate": 2, "condition": "injured"},
    "GSM8546689": {"time_hpi": 24,  "experiment": "Exp2", "replicate": 1, "condition": "injured"},
    "GSM8546690": {"time_hpi": 24,  "experiment": "Exp2", "replicate": 2, "condition": "injured"},
}

SOURCE_GSMS = tuple(SAMPLE_METADATA.keys())

# Confirmed source-level duplicate: GSM8546673 == GSM8546664.
EXCLUDED_GSMS = {"GSM8546673"}
EXCLUDED_GSMS_TUPLE = tuple(sorted(EXCLUDED_GSMS))
SOURCE_DUPLICATE_KEEP = "GSM8546664"
SOURCE_DUPLICATE_PAIRS = (("GSM8546664", "GSM8546673"),)
ANALYSIS_GSMS = tuple(g for g in SOURCE_GSMS if g not in EXCLUDED_GSMS)

EXPECTED_GSMS = ANALYSIS_GSMS

COMMON_TIMEPOINTS = (0, 24, 72)
EXPECTED_EXPERIMENTS = ("Exp1", "Exp2")
EXPECTED_REPLICATES = (1, 2)

EXPECTED_COMMON_WINDOW_TRAJECTORIES = 3

# ------------------------------------------------------------
# QC / normalization
# ------------------------------------------------------------
MIN_TOTAL_COUNTS = 1000
MAX_TOTAL_COUNTS = 25000
MIN_GENES_DETECTED = 200
MAX_GENES_DETECTED = 10000
MAX_MITO_PERCENT = 30.0
NORMALIZATION_TARGET = 10000.0

# ------------------------------------------------------------
# Temporal criteria
# ------------------------------------------------------------
MIN_RELATIVE_CHANGE = 0.20
MIN_TEMPORAL_POINTS = 4
FDR_THRESHOLD = 0.05

RANDOM_SEED = 42

MARKER_PROGRAMS = {
    "neuronal": ["elavl3", "elavl4", "snap25a", "snap25b", "tubb3", "rbfox3"],
    "excitatory": ["slc17a6a", "slc17a6b", "slc17a7"],
    "inhibitory": ["gad1a", "gad1b", "gad2a", "gad2b", "slc6a1a", "slc6a1b"],
    "radial_glia": ["fabp7a", "fabp7b", "slc1a3b", "gfap", "s100b"],
    "progenitor": ["sox2", "neurog1", "neurod1", "ascl1a", "ascl1b", "pax6"],
    "oligodendrocyte": ["olig2", "sox10", "mbpa", "plp1"],
    "microglia": ["mpeg1.1", "aif1", "csf1ra", "apoea"],
}

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
NCBI_TOOL = "GSE278446_NeuroNovelty"
NCBI_DELAY_SECONDS = 0.34

# g:Profiler: Danio rerio. Primary + stable fallback endpoints.
GP_API = "https://biit.cs.ut.ee/gprofiler_beta/api/gost/profile/"
GP_API_STABLE = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GP_ORGANISM = "drerio"

def sample_rows():
    return [
        {
            "GSM": gsm,
            **meta,
            "analysis_status": (
                "EXCLUDED_SOURCE_DUPLICATE"
                if gsm in EXCLUDED_GSMS
                else "ANALYSIS_INCLUDED"
            ),
        }
        for gsm, meta in SAMPLE_METADATA.items()
    ]
