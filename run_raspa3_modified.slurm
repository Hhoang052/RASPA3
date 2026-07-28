#!/usr/bin/env bash
#
# Submit with:
#   sbatch run_raspa3_modified.slurm /absolute/path/to/simulation-folder
#
# The simulation folder must contain simulation.json.
#
#SBATCH --job-name=raspa3-modified
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=raspa3-modified-%j.out
#SBATCH --error=raspa3-modified-%j.err
#SBATCH --signal=B:TERM@120

set -Eeuo pipefail


# =============================================================================
# CONFIGURATION
#
# These defaults match:
#   source tree:  ~/RASPA3-MODIF
#   Conda env:    ~/miniconda3/envs/raspa3-modified
#
# Override without editing this file if necessary:
#   RASPA_MODIFIED_ROOT=/other/path \
#   RASPA_ENV=/other/conda/env \
#   sbatch run_raspa3_modified.slurm /path/to/simulation
# =============================================================================

RASPA_MODIFIED_ROOT="${RASPA_MODIFIED_ROOT:-$HOME/RASPA3-MODIF}"
RASPA_ENV="${RASPA_ENV:-$HOME/miniconda3/envs/raspa3-modified}"
RASPA_EXE="${RASPA_EXE:-$RASPA_MODIFIED_ROOT/build/app/raspa3}"

WORK_DIR="${1:-${SLURM_SUBMIT_DIR:-$PWD}}"
ALLOCATED_CPUS="${SLURM_CPUS_PER_TASK:-1}"


# =============================================================================
# VALIDATION
# =============================================================================

if [[ ! -d "$WORK_DIR" ]]; then
    echo "ERROR: simulation folder does not exist: $WORK_DIR" >&2
    exit 2
fi

WORK_DIR="$(cd "$WORK_DIR" && pwd -P)"

if [[ ! -f "$WORK_DIR/simulation.json" ]]; then
    echo "ERROR: simulation.json not found in: $WORK_DIR" >&2
    exit 2
fi

if [[ ! -x "$RASPA_EXE" ]]; then
    echo "ERROR: modified RASPA3 executable not found or not executable:" >&2
    echo "       $RASPA_EXE" >&2
    echo "Build it first with:" >&2
    echo "       ninja -C \"$RASPA_MODIFIED_ROOT/build\" -j\"\$(nproc)\"" >&2
    exit 2
fi

if [[ ! -d "$RASPA_ENV/lib" ]]; then
    echo "ERROR: Conda environment not found: $RASPA_ENV" >&2
    echo "Set RASPA_ENV to the environment used to build RASPA3-MODIF." >&2
    exit 2
fi


# =============================================================================
# RUNTIME ENVIRONMENT
# =============================================================================

# Loading the environment this way works in non-interactive Slurm shells and
# avoids relying on `conda activate` being initialized in .bashrc.
export PATH="$RASPA_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$RASPA_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# RASPA_DIR must be an installation prefix containing share/raspa3.
# Local files in WORK_DIR can still be used when this directory is unavailable.
if [[ -d "$RASPA_ENV/share/raspa3" ]]; then
    export RASPA_DIR="$RASPA_ENV"
elif [[ -n "${RASPA_DIR:-}" && -d "$RASPA_DIR/share/raspa3" ]]; then
    export RASPA_DIR
else
    echo "WARNING: no valid RASPA_DIR/share/raspa3 directory was found."
    echo "Only force fields, molecules, and frameworks available in the simulation folder can be used."
    echo
fi

export OMP_NUM_THREADS="$ALLOCATED_CPUS"
export OMP_PLACES=cores
export OMP_PROC_BIND=close


# =============================================================================
# START SUMMARY
# =============================================================================

START_TIME_TEXT="$(date --iso-8601=seconds)"
START_TIME_SECONDS="$(date +%s)"

echo "============================================================"
echo "RASPA3 MODIFIED JOB STARTING"
echo "============================================================"
echo "Start time        : $START_TIME_TEXT"
echo "Host              : $(hostname)"
echo "Slurm job ID      : ${SLURM_JOB_ID:-not-running-under-slurm}"
echo "Allocated CPUs    : $ALLOCATED_CPUS"
echo "Simulation folder : $WORK_DIR"
echo "RASPA executable  : $RASPA_EXE"
echo "Conda environment : $RASPA_ENV"
echo "RASPA_DIR         : ${RASPA_DIR:-not-set}"
echo "============================================================"
echo


# =============================================================================
# THREAD CHECK
# =============================================================================

if (( ALLOCATED_CPUS > 1 )); then
    if ! grep -qi '"NumberOfThreads"[[:space:]]*:' \
        "$WORK_DIR/simulation.json"
    then
        echo "WARNING:"
        echo "Slurm allocated $ALLOCATED_CPUS CPU cores."
        echo "However, simulation.json does not contain NumberOfThreads."
        echo "RASPA3 may not use all allocated CPU cores."
        echo "No input file will be modified."
        echo
    fi
fi


# =============================================================================
# RUN RASPA3 MODIFIED
#
# RASPA3 reads simulation.json from the current working directory.
# Results produced by RASPA3 will be written inside the simulation folder.
# =============================================================================

cd "$WORK_DIR"

set +e

srun \
    --ntasks=1 \
    --cpus-per-task="$ALLOCATED_CPUS" \
    --cpu-bind=cores \
    --kill-on-bad-exit=1 \
    "$RASPA_EXE"

RASPA_EXIT_CODE=$?

set -e


# =============================================================================
# FINAL SUMMARY
# =============================================================================

END_TIME_TEXT="$(date --iso-8601=seconds)"
END_TIME_SECONDS="$(date +%s)"
ELAPSED_SECONDS=$((END_TIME_SECONDS - START_TIME_SECONDS))

if [[ "$RASPA_EXIT_CODE" -eq 0 ]]; then
    FINAL_STATUS="SUCCESS"
else
    FINAL_STATUS="FAILED"
fi

echo
echo "============================================================"
echo "RASPA3 MODIFIED JOB FINISHED"
echo "============================================================"
echo "End time          : $END_TIME_TEXT"
echo "Elapsed seconds   : $ELAPSED_SECONDS"
echo "RASPA exit code   : $RASPA_EXIT_CODE"
echo "Result            : $FINAL_STATUS"
echo "Simulation folder : $WORK_DIR"
echo "RASPA executable  : $RASPA_EXE"
echo "============================================================"

exit "$RASPA_EXIT_CODE"
