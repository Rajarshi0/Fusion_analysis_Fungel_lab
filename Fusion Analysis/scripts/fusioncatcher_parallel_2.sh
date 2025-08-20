#!/bin/bash
set -euo pipefail

#================================================================================#
#                           CONFIGURATION & INPUTS                               #
#================================================================================#

# --- INPUTS ---
if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 <R1_FASTQ_GZ> <R2_FASTQ_GZ>"
    exit 1
fi
R1_FILE="$1"
R2_FILE="$2"

# --- MAIN CONFIGURATION ---
THREADS=16
FUSIONCATCHER_DB="../fusioncatcher/data/current"
OUTPUT_BASE_DIR="fusioncatcher_analysis_output"

#================================================================================#
#                           PIPELINE SETUP                                     #
#================================================================================#

# --- Activate Conda Environment ---
eval "$(conda shell.bash hook)"
conda activate fusioncatcher

# --- Extract Sample Name ---
# This logic handles filenames with or without lane info like _L001
#SAMPLE_NAME=$(basename "$R1_FILE" | sed -E 's/(_L[0-9]+)?_R[12](_001)?\.fastq\.gz//')
SAMPLE_NAME=$(basename "$R1_FILE" | sed -E 's/(_L[0-9]+)?_R[12](_001)?(\.fastp)?\.fastq\.gz//')


echo "================================================================="
echo "?? Processing Sample: ${SAMPLE_NAME}"
echo "================================================================="

# --- Define Output Directories ---
TRIMMED_DIR="${OUTPUT_BASE_DIR}/trimmed/${SAMPLE_NAME}"
INPUT_DIR_FC="${OUTPUT_BASE_DIR}/input/${SAMPLE_NAME}"
RESULT_DIR="${OUTPUT_BASE_DIR}/results/${SAMPLE_NAME}"
FINAL_LIST="${RESULT_DIR}/final-list_candidate-fusion-genes.txt"

# --- Skip if already completed successfully ---
if [[ -s "$FINAL_LIST" ]]; then
    echo "? Skipping ${SAMPLE_NAME}: Final output already exists."
    exit 0
fi

# Create clean directories for the run
mkdir -p "$TRIMMED_DIR" "$INPUT_DIR_FC" "$RESULT_DIR"
rm -f "$INPUT_DIR_FC"/* # Clean symlinks from previous runs

#================================================================================#
#                                 ANALYSIS                                     #
#================================================================================#

# --- Step 1: Trim Raw Files with fastp ---
echo "?? Step 1: Trimming raw files with fastp..."
TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"
TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"

fastp \
    -i "$R1_FILE" -I "$R2_FILE" \
    -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
    --thread "$THREADS" \
    --detect_adapter_for_pe \
    --length_required 30 \
    --qualified_quality_phred 20 \
    --html "${TRIMMED_DIR}/fastp_report.html"

# --- Step 2: Prepare Input for FusionCatcher ---
echo "?? Step 2: Creating symlinks to new trimmed files..."
# FusionCatcher requires specific input filenames (e.g., sample_1.fastq.gz)
ln -sf "$(realpath "$TRIMMED_R1")" "$INPUT_DIR_FC/sample_1.fastq.gz"
ln -sf "$(realpath "$TRIMMED_R2")" "$INPUT_DIR_FC/sample_2.fastq.gz"

# --- Step 3: Run FusionCatcher ---
echo "?? Step 3: Running FusionCatcher for ${SAMPLE_NAME}..."
echo "================================================================="

START_TIME=$(date +%s)

if fusioncatcher \
    -d "$FUSIONCATCHER_DB" \
    -i "$INPUT_DIR_FC" \
    -o "$RESULT_DIR" \
    -p "$THREADS" ; then
    
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    
    echo ""
    echo "================================================================="
    echo "? FusionCatcher completed for ${SAMPLE_NAME}"
    echo "??  Runtime: $((RUNTIME / 3600))h $((RUNTIME % 3600 / 60))m $((RUNTIME % 60))s"
    
else
    echo ""
    echo "? FusionCatcher failed for ${SAMPLE_NAME}"
    echo "?? Check logs in: $RESULT_DIR/"
    exit 1
fi

echo "================================================================="
echo "? Processing fully completed for: ${SAMPLE_NAME}"
echo "================================================================="
