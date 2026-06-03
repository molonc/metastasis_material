#!/bin/bash
# ============================================================
# Run DICE with multiple parameter combinations on a single
# input TSV file. Runs k-means clustering and evaluation for
# each combination.
#
# Usage:
#   ./run_DICE_batch.sh -i <input.tsv> [-o <output_dir>] [-k <max_k>]
#
# Example:
#   ./run_DICE_batch.sh -i SA919_total_merged_filtered_states_longformat_rand_0.2.tsv
#   ./run_DICE_batch.sh -i data.tsv -o results -k 10
# ============================================================

set -euo pipefail

# ---- Defaults -----------------------------------------------
CONDA_ENV="/Users/hoatran/anaconda3/envs/DICE"
DICE_BIN="${CONDA_ENV}/bin/dice"
FASTME_BIN="${CONDA_ENV}/bin/fastme"
PYTHON="${CONDA_ENV}/bin/python"
OUTPUT_DIR="DICE_output"
MAX_K=15
INPUT_FILE=""

# DICE parameter grid
METHODS=("NJ" "balME")
DIST_TYPES=("root")
# DIST_TYPES=("root" "log")
USE_NNI=("no" "yes")       # NNI only applies to ME methods
USE_BREAKPOINT=("no" "yes") # -b flag: DICE-bar (breakpoint profiles)

# ---- Parse arguments ----------------------------------------
usage() {
    echo "Usage: $0 -i <input.tsv> [-o <output_dir>] [-k <max_k>]"
    echo ""
    echo "Options:"
    echo "  -i  Input TSV file (total CN format, required)"
    echo "  -o  Output directory (default: DICE_output)"
    echo "  -k  Max k for k-means silhouette selection (default: 15)"
    echo ""
    echo "Runs DICE with all combinations of:"
    echo "  Methods:      ${METHODS[*]}"
    echo "  Distances:    ${DIST_TYPES[*]}"
    echo "  NNI:          off/on (ME methods only)"
    echo "  Breakpoint:   off/on (DICE-star vs DICE-bar)"
    exit 1
}

while getopts "i:o:k:h" opt; do
    case ${opt} in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        k) MAX_K="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT_FILE}" ]; then
    echo "Error: -i <input.tsv> is required"
    usage
fi

if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: input file not found: ${INPUT_FILE}"
    exit 1
fi

# ---- Setup --------------------------------------------------
BASENAME=$(basename "${INPUT_FILE}" .tsv)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "${OUTPUT_DIR}"

LOG_FILE="${OUTPUT_DIR}/batch_run_log.txt"
SUMMARY_FILE="${OUTPUT_DIR}/batch_summary.csv"

echo "============================================================" | tee "${LOG_FILE}"
echo "DICE Batch Run" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "  Input file:   ${INPUT_FILE}" | tee -a "${LOG_FILE}"
echo "  Output dir:   ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "  Max k:        ${MAX_K}" | tee -a "${LOG_FILE}"
echo "  Methods:      ${METHODS[*]}" | tee -a "${LOG_FILE}"
echo "  Distances:    ${DIST_TYPES[*]}" | tee -a "${LOG_FILE}"
echo "  Start time:   $(date)" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

# ---- Step 1: Run k-means once (shared across all DICE runs) -
echo "Step 1: Running k-means clustering (k=2..${MAX_K}) ..." | tee -a "${LOG_FILE}"

KMEANS_DIR="${OUTPUT_DIR}/kmeans"
mkdir -p "${KMEANS_DIR}"

${PYTHON} -c "
import sys
sys.path.insert(0, '${SCRIPT_DIR}')
from generate_and_cluster import read_total_cn_matrix, run_kmeans_silhouette
import pandas as pd
import numpy as np

cell_names, matrix = read_total_cn_matrix('${INPUT_FILE}')
best_k, best_labels = run_kmeans_silhouette(matrix, k_range=range(2, ${MAX_K}+1))

# Save labels
df = pd.DataFrame({'cell': cell_names, 'kmeans_cluster': best_labels})
df.to_csv('${KMEANS_DIR}/kmeans_labels.tsv', sep='\t', index=False)

# Save gt.tsv (k-means as ground truth)
df['true_clone'] = 'cluster_' + df['kmeans_cluster'].astype(str)
df[['cell', 'true_clone', 'kmeans_cluster']].to_csv('${KMEANS_DIR}/gt.tsv', sep='\t', index=False)

print(f'Optimal k = {best_k}')
print(f'Saved kmeans_labels.tsv and gt.tsv to ${KMEANS_DIR}')
" 2>&1 | tee -a "${LOG_FILE}"

echo "" | tee -a "${LOG_FILE}"

# ---- Step 2: Run DICE with each parameter combination -------
echo "Step 2: Running DICE with parameter grid ..." | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

# Write summary CSV header
echo "run_label,method,dist_type,nni,breakpoint,n_cells,n_bins,time_seconds,tree_file,mean_f1,ari,nrfd,matching_accuracy" > "${SUMMARY_FILE}"

RUN_COUNT=0

for METHOD in "${METHODS[@]}"; do
    for DIST_TYPE in "${DIST_TYPES[@]}"; do
        for NNI in "${USE_NNI[@]}"; do
            for BP in "${USE_BREAKPOINT[@]}"; do

            # NNI only applies to ME methods (balME, olsME)
            if [[ "${NNI}" == "yes" && "${METHOD}" != "balME" && "${METHOD}" != "olsME" ]]; then
                continue
            fi

            # Build run label
            NNI_SUFFIX=""
            if [[ "${NNI}" == "yes" ]]; then
                NNI_SUFFIX="_NNI"
            fi
            BP_SUFFIX=""
            if [[ "${BP}" == "yes" ]]; then
                BP_SUFFIX="_bar"
            else
                BP_SUFFIX="_star"
            fi
            RUN_LABEL="${METHOD}_${DIST_TYPE}${NNI_SUFFIX}${BP_SUFFIX}_${BASENAME}"
            RUN_DIR="${OUTPUT_DIR}/${RUN_LABEL}"
            mkdir -p "${RUN_DIR}"

            RUN_COUNT=$((RUN_COUNT + 1))
            echo "------------------------------------------------------------" | tee -a "${LOG_FILE}"
            echo "Run ${RUN_COUNT}: ${RUN_LABEL}" | tee -a "${LOG_FILE}"
            echo "  Method: ${METHOD}, Distance: ${DIST_TYPE}, NNI: ${NNI}, Breakpoint: ${BP}" | tee -a "${LOG_FILE}"
            echo "  Output: ${RUN_DIR}" | tee -a "${LOG_FILE}"
            echo "  Start:  $(date)" | tee -a "${LOG_FILE}"

            # Copy k-means labels into run directory
            cp "${KMEANS_DIR}/kmeans_labels.tsv" "${RUN_DIR}/"
            cp "${KMEANS_DIR}/gt.tsv" "${RUN_DIR}/"

            # Build DICE command
            DICE_CMD="${DICE_BIN} -i ${INPUT_FILE} -o ${RUN_DIR} -d ${DIST_TYPE} -m ${METHOD} -t -f ${FASTME_BIN}"
            if [[ "${NNI}" == "yes" ]]; then
                DICE_CMD="${DICE_CMD} -n"
            fi
            if [[ "${BP}" == "yes" ]]; then
                DICE_CMD="${DICE_CMD} -b"
            fi

            echo "  Command: ${DICE_CMD}" | tee -a "${LOG_FILE}"

            # Run DICE and capture time
            START_TIME=$(date +%s)
            ${DICE_CMD} >> "${LOG_FILE}" 2>&1
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))

            echo "  DICE finished in ${ELAPSED}s" | tee -a "${LOG_FILE}"

            # Find tree file
            TREE_FILE=$(find "${RUN_DIR}" -name "*_tree.nwk" -size +0 | head -1)

            if [ -z "${TREE_FILE}" ]; then
                echo "  WARNING: No tree file found!" | tee -a "${LOG_FILE}"
                echo "${RUN_LABEL},${METHOD},${DIST_TYPE},${NNI},${BP},,${ELAPSED},,,,," >> "${SUMMARY_FILE}"
                continue
            fi

            TREE_SIZE=$(stat -f%z "${TREE_FILE}" 2>/dev/null || stat --printf="%s" "${TREE_FILE}" 2>/dev/null)
            echo "  Tree: ${TREE_FILE} (${TREE_SIZE} bytes)" | tee -a "${LOG_FILE}"

            # Run evaluation
            echo "  Running evaluation ..." | tee -a "${LOG_FILE}"
            EVAL_OUTPUT=$(${PYTHON} "${SCRIPT_DIR}/evaluate_tree_f1.py" \
                "${TREE_FILE}" \
                "${RUN_DIR}/gt.tsv" \
                "${RUN_DIR}" 2>&1)
            echo "${EVAL_OUTPUT}" >> "${LOG_FILE}"

            # Extract metrics from eval_summary.csv
            if [ -f "${RUN_DIR}/eval_summary.csv" ]; then
                METRICS=$(tail -1 "${RUN_DIR}/eval_summary.csv")
                N_LEAVES=$(echo "${METRICS}" | cut -d',' -f3)
                N_CLADES=$(echo "${METRICS}" | cut -d',' -f4)
                MEAN_F1=$(echo "${METRICS}" | cut -d',' -f6)
                ARI=$(echo "${METRICS}" | cut -d',' -f7)
                NRFD=$(echo "${METRICS}" | cut -d',' -f8)
                MATCH_ACC=$(echo "${METRICS}" | cut -d',' -f9)

                echo "  Results: F1=${MEAN_F1}, ARI=${ARI}, NRFD=${NRFD}, Accuracy=${MATCH_ACC}" | tee -a "${LOG_FILE}"
                echo "${RUN_LABEL},${METHOD},${DIST_TYPE},${NNI},${BP},${N_LEAVES},${N_CLADES},${ELAPSED},${TREE_FILE},${MEAN_F1},${ARI},${NRFD},${MATCH_ACC}" >> "${SUMMARY_FILE}"
            else
                echo "  WARNING: eval_summary.csv not found" | tee -a "${LOG_FILE}"
                echo "${RUN_LABEL},${METHOD},${DIST_TYPE},${NNI},${BP},,${ELAPSED},${TREE_FILE},,,," >> "${SUMMARY_FILE}"
            fi

            echo "  End:    $(date)" | tee -a "${LOG_FILE}"
            echo "" | tee -a "${LOG_FILE}"

            done
        done
    done
done

# ---- Step 3: Print final summary ----------------------------
echo "" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "Batch Run Complete" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "  Total runs:    ${RUN_COUNT}" | tee -a "${LOG_FILE}"
echo "  End time:      $(date)" | tee -a "${LOG_FILE}"
echo "  Summary CSV:   ${SUMMARY_FILE}" | tee -a "${LOG_FILE}"
echo "  Log file:      ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

echo "Summary table:"
column -t -s',' "${SUMMARY_FILE}"
