# Running RASPA3-MODIF with Slurm

This guide explains how to run the modified RASPA3 executable on a Slurm
cluster without accidentally using the original RASPA3 installation.

## Expected server layout

The supplied Slurm script uses the following default paths:

```text
Source and build directory:  ~/RASPA3-MODIF
Modified executable:         ~/RASPA3-MODIF/build/app/raspa3
Conda environment:           ~/miniconda3/envs/raspa3-modified
```

The simulation itself must be stored in a separate directory containing at
least:

```text
simulation.json
```

Depending on the simulation, the directory may also contain custom molecule
JSON files, a force-field JSON file, CIF files, restart files, or other input
data.

## 1. Build the modified executable

Activate the environment used for the modified source tree:

```bash
conda activate raspa3-modified
cd ~/RASPA3-MODIF
ninja -C build -j"$(nproc)"
```

Confirm that the executable exists:

```bash
ls -l ~/RASPA3-MODIF/build/app/raspa3
```

Test it:

```bash
~/RASPA3-MODIF/build/app/raspa3 --help
```

Using this absolute path ensures that the modified executable is selected
instead of an original RASPA3 executable found through `PATH`.

## 2. Install the Slurm run script

Copy `run_raspa3_modified.slurm` into the modified source directory:

```bash
cp /path/to/run_raspa3_modified.slurm ~/RASPA3-MODIF/
chmod +x ~/RASPA3-MODIF/run_raspa3_modified.slurm
```

The script uses these defaults:

```bash
RASPA_MODIFIED_ROOT="$HOME/RASPA3-MODIF"
RASPA_ENV="$HOME/miniconda3/envs/raspa3-modified"
RASPA_EXE="$RASPA_MODIFIED_ROOT/build/app/raspa3"
```

If the server uses different paths, either edit these values in the script or
override them when submitting the job.

## 3. Configure Slurm resources

The default resource settings are:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
```

Edit the CPU count and time limit according to the simulation and cluster
rules. A partition or account may also be required, for example:

```bash
#SBATCH --partition=compute
#SBATCH --account=my-project
```

The script starts one RASPA3 process. The CPU cores assigned with
`--cpus-per-task` are available to that process.

## 4. Configure RASPA3 threads

When more than one CPU is requested, set `NumberOfThreads` in
`simulation.json` to the same value as `--cpus-per-task`:

```json
{
  "NumberOfThreads": 8,
  "SimulationType": "MonteCarlo"
}
```

For example:

```text
#SBATCH --cpus-per-task=8
"NumberOfThreads": 8
```

The Slurm script prints a warning when multiple CPUs are allocated but
`NumberOfThreads` is missing. It does not modify `simulation.json`.

## 5. Submit a simulation

Pass the simulation directory as the first argument:

```bash
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm \
    /absolute/path/to/simulation-directory
```

For example:

```bash
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm \
    "$HOME/simulations/methane"
```

The submitted directory must contain:

```text
$HOME/simulations/methane/simulation.json
```

RASPA3 reads `simulation.json` from its current working directory. The script
therefore changes to the supplied simulation directory before starting
RASPA3. Simulation results are written into that directory.

If no simulation directory is supplied, the script uses the directory from
which `sbatch` was called:

```bash
cd "$HOME/simulations/methane"
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm
```

## 6. Override paths without editing the script

If the modified source tree or Conda environment is stored elsewhere:

```bash
RASPA_MODIFIED_ROOT=/work/user/RASPA3-MODIF \
RASPA_ENV=/work/user/conda/envs/raspa3-modified \
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm \
    /absolute/path/to/simulation-directory
```

The executable can also be overridden directly:

```bash
RASPA_EXE=/work/user/RASPA3-MODIF/build/app/raspa3 \
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm \
    /absolute/path/to/simulation-directory
```

## 7. Monitor the job

Show queued and running jobs:

```bash
squeue -u "$USER"
```

Show detailed information for a job:

```bash
scontrol show job JOB_ID
```

Follow standard output:

```bash
tail -f raspa3-modified-JOB_ID.out
```

Follow standard error:

```bash
tail -f raspa3-modified-JOB_ID.err
```

Cancel a job:

```bash
scancel JOB_ID
```

## 8. Output and exit status

At startup, the script reports:

- The exact modified executable.
- The Conda environment.
- The simulation directory.
- The allocated CPU count.
- The Slurm job ID and host.
- The value of `RASPA_DIR`.

At completion, it reports:

- `SUCCESS` or `FAILED`.
- The RASPA3 exit code.
- The elapsed time.
- The simulation directory and executable.

An exit code of `0` means RASPA3 completed successfully. A nonzero exit code
means the simulation or job step failed.

## 9. Shared libraries and `RASPA_DIR`

The script adds the modified Conda environment to `PATH` and
`LD_LIBRARY_PATH`. This works in a non-interactive Slurm shell without
requiring `conda activate`.

When this directory exists:

```text
~/miniconda3/envs/raspa3-modified/share/raspa3
```

the script sets:

```bash
RASPA_DIR="$HOME/miniconda3/envs/raspa3-modified"
```

`RASPA_DIR` must be an installation prefix containing `share/raspa3`; it
should not normally be set directly to `~/RASPA3-MODIF`.

## 10. Troubleshooting

### Modified executable not found

Check the build:

```bash
ninja -C ~/RASPA3-MODIF/build -j"$(nproc)"
find ~/RASPA3-MODIF/build -type f -executable -name raspa3
```

### Missing shared library

Inspect linked libraries:

```bash
ldd ~/RASPA3-MODIF/build/app/raspa3 | grep "not found"
```

Confirm the environment:

```bash
ls ~/miniconda3/envs/raspa3-modified/lib
```

### `simulation.json` not found

Supply the directory containing the file, not the path to the file:

```bash
sbatch ~/RASPA3-MODIF/run_raspa3_modified.slurm \
    /path/to/directory-containing-simulation-json
```

### Confirm that the modified version is used

The Slurm output should show:

```text
RASPA executable  : /home/USER/RASPA3-MODIF/build/app/raspa3
```

Do not use only the command `raspa3` when both original and modified versions
are installed, because shell `PATH` ordering may select the original version.
