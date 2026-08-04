#!/usr/bin/env bash
#
# Submit:
#   sbatch run_raspa3_modified.slurm /absolute/path/to/simulation-folder
#
# The simulation folder must contain:
#   simulation.json
#

#SBATCH --job-name=raspa3-modified
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
# Without this the partition default (DefMemPerNode=UNLIMITED) grants the job
# all node memory, so only one job runs at a time. Calibrate per simulation:
#   sacct -j <JOBID> --format=JobID,MaxRSS
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=raspa3-modified-%j.out
#SBATCH --error=raspa3-modified-%j.err
#SBATCH --signal=B:TERM@120

set -Eeuo pipefail


# =============================================================================
# CONFIGURATION
# =============================================================================

# RASPA3 source directory
RASPA_ROOT="${RASPA_ROOT:-$HOME/RASPA3}"

# Conda environment used to build RASPA3
RASPA_ENV="${RASPA_ENV:-$HOME/miniconda3/envs/raspa3-modified}"

# RASPA3 executable
RASPA_EXE="${RASPA_EXE:-$RASPA_ROOT/build/app/raspa3}"

# Simulation directory passed as the first argument
WORK_DIR="${1:-${SLURM_SUBMIT_DIR:-$PWD}}"

# Number of CPUs allocated by Slurm
ALLOCATED_CPUS="${SLURM_CPUS_PER_TASK:-1}"


# =============================================================================
# VALIDATION
# =============================================================================

if [[ ! -d "$WORK_DIR" ]]; then
    echo "ERROR: Simulation directory does not exist:" >&2
    echo "       $WORK_DIR" >&2
    exit 2
fi

WORK_DIR="$(cd "$WORK_DIR" && pwd -P)"

if [[ ! -f "$WORK_DIR/simulation.json" ]]; then
    echo "ERROR: simulation.json was not found in:" >&2
    echo "       $WORK_DIR" >&2
    exit 2
fi

if [[ ! -x "$RASPA_EXE" ]]; then
    echo "ERROR: RASPA3 executable was not found or is not executable:" >&2
    echo "       $RASPA_EXE" >&2
    echo >&2
    echo "Build RASPA3 first with:" >&2
    echo "       cmake --build \"$RASPA_ROOT/build\" --parallel 16" >&2
    exit 2
fi

if [[ ! -d "$RASPA_ENV/bin" || ! -d "$RASPA_ENV/lib" ]]; then
    echo "ERROR: Invalid Conda environment:" >&2
    echo "       $RASPA_ENV" >&2
    exit 2
fi


# =============================================================================
# CONDA RUNTIME ENVIRONMENT
#
# This method works in non-interactive Slurm shells and does not require
# `conda activate` to be initialized in .bashrc.
# =============================================================================

export CONDA_PREFIX="$RASPA_ENV"
export PATH="$RASPA_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$RASPA_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"


# =============================================================================
# RASPA DATA DIRECTORY
#
# RASPA_DIR must point directly to a directory containing:
#   forcefields/
#   molecules/
#
# Installed data:
#   $RASPA_ENV/share/raspa3
#
# Source-tree data:
#   $RASPA_ROOT/data
# =============================================================================

if [[ -n "${RASPA_DIR:-}" && -d "$RASPA_DIR" ]]; then
    export RASPA_DIR

elif [[ -d "$RASPA_ENV/share/raspa3" ]]; then
    export RASPA_DIR="$RASPA_ENV/share/raspa3"

elif [[ -d "$RASPA_ROOT/data/forcefields" &&
        -d "$RASPA_ROOT/data/molecules" ]]; then
    export RASPA_DIR="$RASPA_ROOT/data"

else
    unset RASPA_DIR || true

    echo "WARNING: No valid RASPA3 data directory was found."
    echo "Only force fields, molecules, and frameworks available in the"
    echo "simulation directory can be used."
    echo
fi


# =============================================================================
# DETECT THREADING MODE
#
# Serial:
#   "ThreadingType": "Serial"
#
# OpenMP:
#   "ThreadingType": "OpenMP"
# =============================================================================

THREADING_TYPE="Serial"

if grep -Eqi \
    '"ThreadingType"[[:space:]]*:[[:space:]]*"OpenMP"' \
    "$WORK_DIR/simulation.json"
then
    THREADING_TYPE="OpenMP"
fi

if [[ "$THREADING_TYPE" == "OpenMP" ]]; then
    RASPA_THREADS="$ALLOCATED_CPUS"
else
    RASPA_THREADS=1

    if (( ALLOCATED_CPUS > 1 )); then
        echo "WARNING:"
        echo "simulation.json is configured for Serial execution,"
        echo "but Slurm allocated $ALLOCATED_CPUS CPU cores."
        echo "RASPA3 will use only one CPU core."
        echo
    fi
fi


# =============================================================================
# THREAD ENVIRONMENT
# =============================================================================

# OpenMP threads used by RASPA3
export OMP_NUM_THREADS="$RASPA_THREADS"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_WAIT_POLICY=PASSIVE

# Keep BLAS libraries single-threaded to avoid nested parallelism
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1


# =============================================================================
# CHECK RUNTIME LIBRARIES
# =============================================================================

MISSING_LIBRARIES="$(
    ldd "$RASPA_EXE" 2>/dev/null |
    awk '/not found/ { print }'
)"

if [[ -n "$MISSING_LIBRARIES" ]]; then
    echo "ERROR: The RASPA3 executable has missing runtime libraries:" >&2
    echo "$MISSING_LIBRARIES" >&2
    echo >&2
    echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >&2
    exit 2
fi


# =============================================================================
# START SUMMARY
# =============================================================================

START_TIME_TEXT="$(date --iso-8601=seconds)"
START_TIME_SECONDS="$(date +%s)"

echo "============================================================"
echo "RASPA3 MODIFIED JOB STARTING"
echo "============================================================"
echo "Start time         : $START_TIME_TEXT"
echo "Host               : $(hostname)"
echo "Slurm job ID       : ${SLURM_JOB_ID:-not-running-under-slurm}"
echo "Allocated CPUs     : $ALLOCATED_CPUS"
echo "Threading type     : $THREADING_TYPE"
echo "OMP threads        : $OMP_NUM_THREADS"
echo "Simulation folder  : $WORK_DIR"
echo "RASPA executable   : $RASPA_EXE"
echo "Conda environment  : $RASPA_ENV"
echo "RASPA_DIR          : ${RASPA_DIR:-not-set}"
echo "============================================================"
echo


# =============================================================================
# RUN RASPA3
#
# RASPA3 reads simulation.json from the current working directory.
# Output files are written inside the simulation directory.
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
echo "End time           : $END_TIME_TEXT"
echo "Elapsed seconds    : $ELAPSED_SECONDS"
echo "RASPA exit code    : $RASPA_EXIT_CODE"
echo "Result             : $FINAL_STATUS"
echo "Simulation folder  : $WORK_DIR"
echo "RASPA executable   : $RASPA_EXE"
echo "============================================================"

exit "$RASPA_EXIT_CODE"