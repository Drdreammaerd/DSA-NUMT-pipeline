#!/bin/bash
set -e

# ==============================================================================
# run_local.sh
# 
# Executes the DSA-NUMT pipeline LOCALLY using the Docker container.
# This mode uses your local CPU cores and does NOT submit to an LSF cluster.
# ==============================================================================

# Get absolute path to the pipeline script directory
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(pwd)"

# Determine number of cores to use (defaults to all)
CORES=${CORES:-all}

echo "================================================================="
echo "[*] Starting DSA-NUMT Pipeline in LOCAL mode"
echo "[*] Pipeline Dir: $PIPELINE_DIR"
echo "[*] Working Dir:  $WORK_DIR"
echo "[*] Cores:        $CORES"
echo "================================================================="

# The Docker container contains Snakemake and all dependencies.
# We mount the pipeline directory to /opt/numt-dsa-pipeline inside the container (read-only)
# and mount the current working directory so inputs/outputs can be accessed.
# Any additional mounts (e.g. for reference genomes on other drives) must be added as -v flags below.

docker run --rm \
    -v "${PIPELINE_DIR}:/opt/numt-dsa-pipeline:ro" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${WORK_DIR}" \
    dreammaerd/numt_dsa:v1 \
    snakemake -s /opt/numt-dsa-pipeline/Snakefile \
    --cores "${CORES}" \
    "$@"
