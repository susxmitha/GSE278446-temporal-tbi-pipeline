from pathlib import Path
import os
import subprocess
import sys
import time
from datetime import datetime

# ============================================================
# GSE278446 NeuroNovelty — AUTHORITATIVE PIPELINE RUNNER
# ============================================================
# Each stage runs as an isolated child process with a hard
# timeout. A stage that hangs is killed; the runner stops and
# tells you which stage and after how long.
#
# Environment overrides:
#   GSE278446_RESEARCH_ROOT
#       Root containing GSE278446_RAW/ and NeuroNovelty/
#
#   GSE278446_STAGE_TIMEOUT_SECONDS
#       Per-stage timeout. Default 1800 (30 minutes).
#
#   GSE278446_START_AT
#       Numeric prefix of the first stage to run. E.g. 18
#       skips stages 00 through 17.
#
#   GSE278446_STOP_AFTER
#       Numeric prefix of the last stage to run. E.g. 18 runs
#       only stage 18.
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESEARCH_ROOT = SCRIPT_DIR.parents[1]
RESEARCH_ROOT = os.environ.get(
    "GSE278446_RESEARCH_ROOT",
    str(DEFAULT_RESEARCH_ROOT),
)

STAGE_TIMEOUT_SECONDS = int(
    os.environ.get("GSE278446_STAGE_TIMEOUT_SECONDS", "1800")
)

START_AT = os.environ.get("GSE278446_START_AT", "").strip()
STOP_AFTER = os.environ.get("GSE278446_STOP_AFTER", "").strip()


TESTS = [
    "00_raw_file_smoke_test.py",
    "01_single_matrix_load_test.py",
    "02_single_matrix_qc_marker_test.py",
    "03_dataset_inventory.py",
    "04_sample_identity_audit.py",
    "05A_verify_duplicate_against_geo.py",
    "05_raw_matrix_sha256_verification.py",
    "06_single_cell_preprocessing_and_gene_means.py",
    "07_program_temporal_discovery.py",
    "08_gene_temporal_discovery.py",
    "09_program_replicate_validation.py",
    "10_common_window_reversal.py",
    "11_strict_replicate_reversal_validation.py",
    "12_temporal_sampling_shift_sensitivity.py",
    "13_marker_program_overlap_audit.py",
    "14_cell_level_reversal_validation.py",
    "15_gene_characterization.py",
    "16_gene_prioritization.py",
    "17_prepare_enrichment_inputs.py",
    "18_run_functional_enrichment.py",
    "19_pathway_reproducibility_interpretation.py",
    "20_ncbi_alias_ortholog_audit.py",
    "21_literature_novelty_audit.py",
    "22_fast_literature_crosscheck.py",
    "23_targeted_literature_verification.py",
    "24_targeted_paper_level_verification.py",
    "25_independent_literature_recheck.py",
    "26_candidate_biological_context_classification.py",
    "27_candidate_evidence_audit.py",
    "28_sensitivity_robustness.py",
    "29_final_candidate_selection.py",
    "30_final_paper_recheck.py",
    "31_pipeline_final_audit.py",
]

RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR = SCRIPT_DIR / "stage1_run_logs" / RUN_ID
LOG_DIR.mkdir(parents=True, exist_ok=True)

MASTER_LOG = LOG_DIR / "RUN_ALL.log"


def write_log(text):
    print(text, flush=True)
    with MASTER_LOG.open("a", encoding="utf-8") as f:
        f.write(text + "\n")


def stage_number(name):
    """
    Return the numeric prefix of a stage filename as a string.
    '05A_verify...' -> '05A'. '18_run_...' -> '18'.
    """
    return name.split("_", 1)[0]


def should_run(name):
    """
    Apply START_AT and STOP_AFTER filters.
    Both are compared lexicographically on the stage prefix.
    '05A' sorts after '05', which is what we want for the
    duplicate-verification script.
    """
    num = stage_number(name)

    if START_AT and num < START_AT:
        return False

    if STOP_AFTER and num > STOP_AFTER:
        return False

    return True


def run_stage(filename, index, total):
    """
    Run one stage as a subprocess with a hard timeout.
    Returns True on success, False on failure or timeout.
    """
    stage_path = SCRIPT_DIR / filename
    stage_log = LOG_DIR / f"{index:02d}_{filename}.log"

    if not stage_path.exists():
        write_log(f"FAILED: file not found: {stage_path}")
        return False

    write_log("")
    write_log("-" * 72)
    write_log(f"[{index:02d}/{total}] RUNNING: {filename}")
    write_log(f"        timeout : {STAGE_TIMEOUT_SECONDS}s")
    write_log(f"        log     : {stage_log}")
    write_log("-" * 72)

    env = os.environ.copy()
    env["GSE278446_RESEARCH_ROOT"] = str(RESEARCH_ROOT)
    # Force line-buffered stdout/stderr from the child so we see
    # progress in real time.
    env["PYTHONUNBUFFERED"] = "1"

    started = time.monotonic()

    try:
        with stage_log.open("w", encoding="utf-8") as log_fh:
            proc = subprocess.Popen(
                [sys.executable, str(stage_path)],
                cwd=str(SCRIPT_DIR),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            # Stream output line by line with a timeout watchdog.
            try:
                for line in proc.stdout:
                    print(line, end="", flush=True)
                    log_fh.write(line)

                    elapsed = time.monotonic() - started
                    if elapsed > STAGE_TIMEOUT_SECONDS:
                        proc.kill()
                        raise TimeoutError(
                            f"Stage exceeded {STAGE_TIMEOUT_SECONDS}s"
                        )

                proc.wait(timeout=STAGE_TIMEOUT_SECONDS)

            except TimeoutError:
                try:
                    proc.kill()
                    proc.wait(timeout=10)
                except Exception:
                    pass
                raise

    except TimeoutError as exc:
        elapsed = time.monotonic() - started
        write_log("")
        write_log(f"TIMEOUT: {filename} after {elapsed:.1f}s")
        write_log(f"        Full output: {stage_log}")
        write_log("Pipeline STOPPED.")
        return False

    except KeyboardInterrupt:
        try:
            proc.kill()
        except Exception:
            pass
        write_log("")
        write_log(f"ABORTED by user during {filename}")
        return False

    except Exception as exc:
        elapsed = time.monotonic() - started
        write_log("")
        write_log(f"FAILED: {filename} after {elapsed:.1f}s — {exc}")
        write_log(f"        Full output: {stage_log}")
        write_log("Pipeline STOPPED.")
        return False

    elapsed = time.monotonic() - started

    if proc.returncode != 0:
        write_log("")
        write_log(
            f"FAILED: {filename} exited with code {proc.returncode} "
            f"after {elapsed:.1f}s"
        )
        write_log(f"        Full output: {stage_log}")
        write_log("Pipeline STOPPED.")
        return False

    write_log(f"COMPLETED: {filename} | elapsed: {elapsed:.1f}s")
    return True


def main():
    os.chdir(SCRIPT_DIR)

    root = Path(RESEARCH_ROOT).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"RESEARCH_ROOT does not exist:\n{root}")

    raw = root / "GSE278446_RAW"
    if not raw.exists():
        raise FileNotFoundError(
            f"GSE278446_RAW directory not found:\n{raw}"
        )

    os.environ["GSE278446_RESEARCH_ROOT"] = str(root)

    if str(SCRIPT_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPT_DIR))

    selected = [
        name for name in TESTS if should_run(name)
    ]

    write_log("=" * 72)
    write_log("GSE278446 — FINAL REPRODUCIBLE PIPELINE RUN")
    write_log("=" * 72)
    write_log(f"Started         : {datetime.now().isoformat(timespec='seconds')}")
    write_log(f"Research root   : {root}")
    write_log(f"Script dir      : {SCRIPT_DIR}")
    write_log(f"Stages in list  : {len(TESTS)}")
    write_log(f"Stages to run   : {len(selected)}")
    write_log(f"Stage timeout   : {STAGE_TIMEOUT_SECONDS}s")
    if START_AT:
        write_log(f"Start at        : {START_AT}")
    if STOP_AFTER:
        write_log(f"Stop after      : {STOP_AFTER}")
    write_log(f"Run log dir     : {LOG_DIR}")
    write_log(f"Master log      : {MASTER_LOG}")
    write_log("")

    if not selected:
        write_log("No stages matched the START_AT / STOP_AFTER filter.")
        return

    overall_start = time.monotonic()

    for i, filename in enumerate(selected, start=1):
        ok = run_stage(filename, i, len(selected))
        if not ok:
            write_log("")
            write_log("=" * 72)
            write_log("PIPELINE STOPPED")
            write_log("=" * 72)
            sys.exit(1)

    total_elapsed = time.monotonic() - overall_start

    write_log("")
    write_log("=" * 72)
    write_log(f"ALL {len(selected)} STAGES COMPLETED SUCCESSFULLY")
    write_log("=" * 72)
    write_log(f"Finished      : {datetime.now().isoformat(timespec='seconds')}")
    write_log(f"Total elapsed : {total_elapsed:.1f}s")
    write_log(f"Master log    : {MASTER_LOG}")
    write_log(f"Per-stage logs: {LOG_DIR}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print()
        print("Aborted by user.")
        sys.exit(130)
