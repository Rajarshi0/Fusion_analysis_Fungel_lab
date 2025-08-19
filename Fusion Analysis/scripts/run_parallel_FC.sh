#!/bin/bash

# ====== RESOURCE CONFIGURATION ======
THREADS_PER_JOB=16
MAX_PARALLEL_JOBS=6          # Your system can handle much more with 1TB RAM!
SCRIPT="fusioncatcher_parallel_2.sh"
LOGDIR="fusioncatcher_logs"

# Memory management
MEMORY_CHECK_INTERVAL=30     # Check memory every 30 seconds
MIN_FREE_MEMORY_GB=200        # Keep 200GB free as buffer with your 1TB system

mkdir -p "$LOGDIR"

# Function to check available memory in GB
check_memory() {
    local free_mem_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    local free_mem_gb=$((free_mem_kb / 1024 / 1024))
    echo $free_mem_gb
}

# Function to wait for memory to be available
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

# Function to get current running job count
get_running_jobs() {
    jobs -r | wc -l
}

echo "================================================================="
echo "?? FusionCatcher Batch Processing"
echo "================================================================="
echo "Configuration:"
echo "  - Max parallel jobs: $MAX_PARALLEL_JOBS"
echo "  - Threads per job: $THREADS_PER_JOB"
echo "  - Total system load: $((MAX_PARALLEL_JOBS * THREADS_PER_JOB)) threads max"
echo "  - Minimum free memory required: ${MIN_FREE_MEMORY_GB}GB"
echo "================================================================="

# ====== FIND FASTQ PAIRS ======
# Your files have pattern: *_L001_R1_001.fastq.gz
#R1_FILES=(*_L001_R1_001.fastq.gz)
R1_FILES=(*_R1.fastp.fastq.gz)
SAMPLES=()

for r1 in "${R1_FILES[@]}"; do
    # Extract sample name: 40907403471-P486-PB-NCGM-2732_L001_R1_001.fastq.gz -> 40907403471-P486-PB-NCGM-2732
    sample_prefix="${r1%_R1.fastp.fastq.gz}"
    r2="${sample_prefix}_R2.fastp.fastq.gz"
    if [[ -f "$r2" ]]; then
        SAMPLES+=("$sample_prefix")
        echo "Found pair: $sample_prefix"
    else
        echo "Missing R2 for $sample_prefix, skipping..."
    fi
done

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "? No valid sample pairs found!"
    exit 1
fi

echo "Found ${#SAMPLES[@]} sample pairs to process:"
for sample in "${SAMPLES[@]}"; do
    echo "  - $sample"
done
echo "================================================================="

# ====== RUN IN PARALLEL BATCHES WITH MEMORY MANAGEMENT ======
total_samples=${#SAMPLES[@]}
completed=0

for sample in "${SAMPLES[@]}"; do
    R1="${sample}_R1.fastp.fastq.gz"
    R2="${sample}_R2.fastp.fastq.gz"
    LOGFILE="$LOGDIR/${sample##*/}.log"
    
    # Wait for available job slot
    while [ $(get_running_jobs) -ge $MAX_PARALLEL_JOBS ]; do
        echo "[$(date)] Max jobs running ($(get_running_jobs)/$MAX_PARALLEL_JOBS), waiting..."
        sleep 10
        # Check for completed jobs
        wait -n 2>/dev/null || true
    done
    
    # Wait for sufficient memory
    wait_for_memory
    
    echo "[$(date)] Launching [$((completed+1))/$total_samples]: $sample"
    nohup bash "$SCRIPT" "$R1" "$R2" > "$LOGFILE" 2>&1 &
    
    # Give the job a moment to start and allocate initial memory
    sleep 5
done

# Wait for all remaining jobs to complete
echo "[$(date)] All jobs launched. Waiting for completion..."
while [ $(get_running_jobs) -gt 0 ]; do
    running=$(get_running_jobs)
    free_mem=$(check_memory)
    echo "[$(date)] Jobs still running: $running, Free memory: ${free_mem}GB"
    sleep 30
done

echo "================================================================="
echo "? All FusionCatcher jobs completed!"
echo "================================================================="

# Summary report
echo "?? Processing Summary:"
echo "  - Total samples processed: $total_samples"
echo "  - Log files location: $LOGDIR/"
echo ""
echo "?? Quick status check:"
for sample in "${SAMPLES[@]}"; do
    result_file="fusioncatcher_analysis_output/results/${sample}/final-list_candidate-fusion-genes.txt"
    if [[ -s "$result_file" ]]; then
        fusion_count=$(wc -l < "$result_file" 2>/dev/null || echo "0")
        echo "  ? ${sample}: $fusion_count candidates"
    else
        echo "  ? ${sample}: FAILED or incomplete"
    fi
done
