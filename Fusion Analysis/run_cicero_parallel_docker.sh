#!/bin/bash

# ====== RESOURCE CONFIGURATION ======
THREADS_PER_JOB=16
MAX_PARALLEL_JOBS=2         # Adjust based on your system's capacity
LOG_DIR="cicero_logs"

# ====== PATH CONFIGURATION (FINAL) ======
# The base path that will be mounted into the Docker container
HOST_BASE_PATH="/Fungel/rajarshi"

# Your current project directory
PROJECT_DIR="${HOST_BASE_PATH}/02_processed_data_08-1-2025"

# Base directory where Cicero will write its output
CICERO_BASE_OUTPUT_DIR="${PROJECT_DIR}/cicero_output"

# Memory management
MEMORY_CHECK_INTERVAL=30
MIN_FREE_MEMORY_GB=100

mkdir -p "${PROJECT_DIR}/${LOG_DIR}"
mkdir -p "$CICERO_BASE_OUTPUT_DIR"

# ====== REFERENCE FILE PATHS (VERIFY THESE ARE CORRECT) ======
REF_FASTA="${HOST_BASE_PATH}/cicero/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"
REF_FLAT="${HOST_BASE_PATH}/cicero/refFlat_hg38.txt"
CICERO_REF_DIR="${HOST_BASE_PATH}/cicero/reference"

# ====== HELPER FUNCTIONS (Memory & Job Management) ======
check_memory() {
    local free_mem_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    local free_mem_gb=$((free_mem_kb / 1024 / 1024))
    echo $free_mem_gb
}

wait_for_memory() {
    while true; do
        local free_mem=$(check_memory)
        if [ $free_mem -ge $MIN_FREE_MEMORY_GB ]; then
            echo "[$(date)] Memory OK: ${free_mem}GB available"
            break
        else
            echo "[$(date)] Low memory: ${free_mem}GB available, waiting..."
            sleep $MEMORY_CHECK_INTERVAL
        fi
    done
}

get_running_jobs() {
    jobs -p | wc -l
}

echo "================================================================="
echo "🔬 Parallel Cicero Analysis (with Auto-Indexing)"
echo "================================================================="
echo "  - Project Directory: $PROJECT_DIR"
echo "  - Max parallel jobs: $MAX_PARALLEL_JOBS"
echo "================================================================="

# ====== FIND SAMPLES FROM FASTQ FILES ======
cd "$PROJECT_DIR" || { echo "Error: Could not change to $PROJECT_DIR"; exit 1; }

SAMPLES=()
for r1_file in *_R1.fastp.fastq.gz; do
    SAMPLES+=("$(basename "$r1_file" _R1.fastp.fastq.gz)")
done

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "❌ No sample FASTQ files (*_R1.fastp.fastq.gz) found in $PROJECT_DIR!"
    exit 1
fi

echo "Found ${#SAMPLES[@]} samples to process:"
printf "  - %s\n" "${SAMPLES[@]}"
echo "================================================================="

# ====== RUN IN PARALLEL BATCHES ======
total_samples=${#SAMPLES[@]}

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    
    # 1. Define the unusual directory name and filename prefix
    filename_prefix="${sample}_R1.fastp.fastq.gz"
    
    # 2. Build the full path using the prefix for both the directory and the filename
    BAM_FILE="${PROJECT_DIR}/fusion_analysis_output/aligned/${filename_prefix}/${filename_prefix}_Aligned.sortedByCoord.out.bam"
    JUNCTIONS_FILE="${PROJECT_DIR}/fusion_analysis_output/aligned/${filename_prefix}/${filename_prefix}_Chimeric.out.junction"
    
    SAMPLE_CICERO_OUT_DIR="${CICERO_BASE_OUTPUT_DIR}/${sample}"
    LOGFILE="${PROJECT_DIR}/${LOG_DIR}/${sample}.log"

    # Verify the main BAM file exists
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "⚠️ WARNING: Cannot find BAM file for sample $sample. Skipping."
        echo "   - Looked for: $BAM_FILE"
        continue
    fi

    # 3. Check for the BAM index file (.bai) and create it if it doesn't exist.
    BAM_INDEX_FILE="${BAM_FILE}.bai"
    if [[ ! -f "$BAM_INDEX_FILE" ]]; then
        echo "  -> Index (.bai) not found for $sample. Creating it now with samtools..."
        samtools index "$BAM_FILE"
        echo "  -> Index created successfully."
    fi

    mkdir -p "$SAMPLE_CICERO_OUT_DIR"
    
    # Wait for an available job slot
    while [ $(get_running_jobs) -ge $MAX_PARALLEL_JOBS ]; do
        echo "[$(date)] Max jobs running ($(get_running_jobs)/$MAX_PARALLEL_JOBS), waiting..."
        sleep 15
        wait -n 2>/dev/null || true
    done
    
    wait_for_memory
    
    echo "🚀 Launching [$((i+1))/$total_samples]: $sample"
    
    (
      echo "Starting RNApeg for $sample"
      singularity run --rm -v "${HOST_BASE_PATH}:${HOST_BASE_PATH}" ghcr.io/stjude/rnapeg:latest \
        -b "$BAM_FILE" \
        -f "$REF_FASTA" \
        -r "$REF_FLAT" \
        -O "$SAMPLE_CICERO_OUT_DIR"

      echo "RNApeg finished. Starting Cicero."
      singularity run --containall --bind "${HOST_BASE_PATH}:${HOST_BASE_PATH}" ghcr.io/stjude/cicero:latest \
        -n "$THREADS_PER_JOB" \
        -b "$BAM_FILE" \
        -g GRCh38_no_alt \
        -r "$CICERO_REF_DIR" \
        -j "$JUNCTIONS_FILE" \
        -o "$SAMPLE_CICERO_OUT_DIR"
      echo "✅ Cicero finished for $sample"
    ) > "$LOGFILE" 2>&1 &
    
    sleep 5
done

echo "[$(date)] All jobs launched. Waiting for completion..."
wait

echo "================================================================="
echo "🎉 All Cicero jobs completed!"
echo "================================================================="
