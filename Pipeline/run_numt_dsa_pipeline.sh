#!/bin/bash
# ============================================================
# NUMT-DSA Pipeline — Single Command Runner (V1)
# ============================================================
# Usage:
#   bash run_numt_dsa_pipeline.sh --configfile config.yaml --subscription
# ============================================================

set -euo pipefail

SNAKEMAKE_ARGS=()
CONFIG_FILE="$(pwd)/config.yaml"

# Default LSF parameters
LSF_QUEUE="general"
LSF_GROUP="compute-jin810"

while [[ $# -gt 0 ]]; do
    case $1 in
        --configfile)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --subscription)
            LSF_QUEUE="subscription"
            LSF_GROUP="compute-jin810-t3"
            shift 1
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift 1
            ;;
    esac
done

if [[ "$CONFIG_FILE" != /* ]]; then
    CONFIG_FILE="$(pwd)/${CONFIG_FILE}"
fi
CONFIG_DIR="$(cd "$(dirname "$CONFIG_FILE")" && pwd)"

# Resolve path helper
resolve_path() {
    local val="$1"
    if [[ "$val" == /* ]]; then echo "$val"; else echo "${CONFIG_DIR}/${val}"; fi
}

# Parse configs
SAMPLE_SHEET_RAW="$(grep 'sample_sheet:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)"
SAMPLE_SHEET=$(resolve_path "$SAMPLE_SHEET_RAW")

TARGET_FASTA=$(grep 'target_fasta:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)

CHAIN_DIR_RAW="$(grep 'chain_dir:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)"
CHAIN_OUT_DIR=$(resolve_path "$CHAIN_DIR_RAW")

OUTPUT_BASE_RAW="$(grep 'output_base:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)"
OUTPUT_BASE=$(resolve_path "$OUTPUT_BASE_RAW")

# Docker Images
DOCKER_IMAGE="dreammaerd/numt_dsa:v1.2"
NFLO_IMAGE="dreammaerd/nf-lo:v3"

mkdir -p "$OUTPUT_BASE/logs"
mkdir -p "$CHAIN_OUT_DIR"

echo "============================================================"
echo " NUMT-DSA Pipeline (Dockerized V1)"
echo "============================================================"
echo " Config:    ${CONFIG_FILE}"
echo " Log dir:   ${OUTPUT_BASE}/logs"
echo " LSF Queue: ${LSF_QUEUE} (${LSF_GROUP})"
echo "============================================================"

# Required volumes
MOUNT_SCRATCH1="${SCRATCH1:-/scratch1/fs1/jin810}"
MOUNT_STORAGE1="${STORAGE1:-/storage1/fs1/jin810/Active}"
MOUNT_STORAGE2="${STORAGE2:-/storage2/fs1/epigenome/Active}"
VOLUMES_STR="${MOUNT_SCRATCH1}:${MOUNT_SCRATCH1} ${MOUNT_STORAGE1}:${MOUNT_STORAGE1} ${MOUNT_STORAGE2}:${MOUNT_STORAGE2} ${HOME}:${HOME}"
export LSF_DOCKER_VOLUMES="${VOLUMES_STR}"

# Fix XDG
export XDG_CACHE_HOME="/tmp/${USER}/.cache"
export XDG_CONFIG_HOME="/tmp/${USER}/.config"
export MPLCONFIGDIR="/tmp/${USER}/.matplotlib"
mkdir -p "${XDG_CACHE_HOME}" "${XDG_CONFIG_HOME}" "${MPLCONFIGDIR}"

JOB_IDS=()

echo "[1] Submitting nf-LO Chain Generation Jobs..."
while IFS=$'\t' read -r donor haplotype fasta fraction; do
    # Skip header
    if [[ "$donor" == "DonorID" || -z "$donor" ]]; then continue; fi

    sample="${donor}_${haplotype}"
    out_chain="${CHAIN_OUT_DIR}/${sample}_hg38"
    
    if [ -f "${out_chain}/chainnet/liftover.chain" ]; then
        echo "  - Skipping ${sample} (chain already exists)"
        continue
    fi

    # Submit nf-LO
    BSUB_OUT=$(bsub -J "nfLO_${sample}" \
         -o "${OUTPUT_BASE}/logs/nfLO_${sample}.log" -e "${OUTPUT_BASE}/logs/nfLO_${sample}.err" \
         -G "$LSF_GROUP" -q "$LSF_QUEUE" -n 5 \
         -R "select[mem>200000] rusage[mem=200000] span[hosts=1]" \
         -a "docker(${NFLO_IMAGE})" \
         bash -c "
export PATH=\"/venv/bin/:\$PATH\"
set -euo pipefail
mkdir -p \"${out_chain}\"
nextflow run /app/nf-LO/main.nf \
    --source \"${fasta}\" \
    --target \"${TARGET_FASTA}\" \
    --aligner minimap2 \
    --distance near \
    --outdir \"${out_chain}\" \
    --liftover_name liftover \
    --max_cpus 5 \
    --max_memory 200.GB \
    -process.maxForks 1 \
    -w \"${out_chain}/work\"
")
    
    JOB_ID=$(echo "$BSUB_OUT" | grep -oE '[0-9]+' | head -n 1)
    JOB_IDS+=("$JOB_ID")
    echo "  - Submitted ${sample} -> Job ID: ${JOB_ID}"
done < "$SAMPLE_SHEET"

echo ""
echo "[2] Submitting NUMT-DSA Orchestrator..."

DEP_ARGS=()
if [ ${#JOB_IDS[@]} -gt 0 ]; then
    COND_STR=""
    for id in "${JOB_IDS[@]}"; do
        if [ -z "$COND_STR" ]; then COND_STR="done($id)"
        else COND_STR="$COND_STR && done($id)"
        fi
    done
    DEP_ARGS=("-w" "$COND_STR")
fi

rm -rf "${OUTPUT_BASE}/.snakemake/locks" 2>/dev/null || true

bsub -J "DSA_NUMT_Orchestrator" \
     -o "${OUTPUT_BASE}/DSA_NUMT_Pipeline_%J.log" -e "${OUTPUT_BASE}/DSA_NUMT_Pipeline_%J.err" \
     -G "$LSF_GROUP" -q "$LSF_QUEUE" -n 1 \
     -R "span[hosts=1] rusage[mem=8000]" \
     -a "docker(${DOCKER_IMAGE})" \
     ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} \
     snakemake -s "/opt/numt-dsa-pipeline/Snakefile" \
         --directory "${CONFIG_DIR}" \
         --configfile "${CONFIG_FILE}" \
         --cluster-generic-submit-cmd "LSF_DOCKER_VOLUMES='${VOLUMES_STR}' LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
             -G ${LSF_GROUP} \
             -q ${LSF_QUEUE} \
             -R 'rusage[mem={resources.mem_mb}]' \
             -a 'docker({resources.docker})' \
             -J {rule} \
             -o ${OUTPUT_BASE}/logs/{rule}_%J.out \
             -e ${OUTPUT_BASE}/logs/{rule}_%J.err" \
         --executor cluster-generic \
         --cluster-generic-status-cmd "python3 /opt/numt-dsa-pipeline/helpers/lsf_status.py" \
         --jobs 100 \
         --latency-wait 120 \
         --retries 3 \
         --rerun-incomplete \
         --keep-going \
         --envvars XDG_CACHE_HOME XDG_CONFIG_HOME MPLCONFIGDIR \
         "${SNAKEMAKE_ARGS[@]:+${SNAKEMAKE_ARGS[@]}}"

echo ""
echo "============================================================"
echo " Submitted! Pipeline is running autonomously."
echo "============================================================"
