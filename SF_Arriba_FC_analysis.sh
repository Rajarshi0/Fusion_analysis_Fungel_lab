echo "--- [${SAMPLE_NAME}] Starting analysis for tools: ${TOOLS_TO_RUN} ---"#!/bin/bash
set -euo pipefail

# ===================================================================================
# --- Self-Cleaning Function to Remove Old/Stopped Processes ---
# ===================================================================================
cleanup_old_runs() {
    # Get the name and Process ID (PID) of the currently running script
    local script_name
    script_name=$(basename "$0")
    local current_pid=$$

    # Find all PIDs matching the script name, but exclude the current script's PID
    # xargs --no-run-if-empty ensures kill is not run if no other processes are found
    local old_pids
    old_pids=$(pgrep -f "$script_name" | grep -v "^${current_pid}$")

    if [ -n "$old_pids" ]; then
        echo "🧹 Found leftover processes from previous runs. Cleaning up..."
        # The '2>/dev/null || true' part silences the "No such process" error
        # and ensures the script doesn't exit even if kill fails.
        echo "$old_pids" | xargs --no-run-if-empty kill -9 2>/dev/null || true
        echo "   -> Cleanup complete."
    fi
}

# Run the cleanup function at the very beginning of the script
cleanup_old_runs
# ===================================================================================


# --- Default Configuration ---
THREADS=16
OUTPUT_BASE_DIR="fusion_analysis_output"
MAX_JOBS=2
FORCE_RUN="false"
RUN_TOOLS="all"

# --- Function to Display Usage ---
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "This script runs a modular fusion detection pipeline with robust, step-by-step resume."
    echo "It automatically cleans up the output of failed steps and old processes before retrying."
    echo
    echo "--------------------------------------------------------------------------------"
    echo "Input Mode (choose ONE):"
    echo "  -I, --input-dir <dir>    Batch Mode: Directory with FASTQ pairs."
    echo "  -1, --r1 <file>          Single Sample: Path to the R1 FASTQ file."
    echo "  -2, --r2 <file>          Single Sample: Path to the R2 FASTQ file."
    echo "  -s, --sample <name>      Single Sample: A name for the sample (required)."
    echo "--------------------------------------------------------------------------------"
    echo
    echo "Reference Arguments (use --ref-bundle-dir for ease-of-use):"
    echo "  -R, --ref-bundle-dir <dir> Optional: Directory containing all reference files."
    echo "  --ctat-lib <dir>           Path to the CTAT Genome Library directory."
    echo "  --star-index <dir>         Path to the STAR index directory."
    echo "  --arriba-path <path>       Path to the Arriba executable."
    echo "  --arriba-blacklist <file>  Path to the Arriba blacklist TSV file."
    echo "  --fusioncatcher-db <dir>   Path to the FusionCatcher database directory."
    echo "--------------------------------------------------------------------------------"
    echo
    echo "Analysis Control:"
    echo "  --run-tools <tool>         (Optional) Which fusion callers to run. Options: starfusion, arriba, fusioncatcher, all, or comma-separated combinations (e.g., starfusion,arriba). (Default: ${RUN_TOOLS})"
    echo
    echo "Optional Arguments:"
    echo "  -o, --output-dir <dir>     Base directory for all output. (Default: ${OUTPUT_BASE_DIR})"
    echo "  -t, --threads <int>        Number of threads per job. (Default: ${THREADS})"
    echo "  -j, --max-jobs <int>       Max parallel jobs in batch mode. (Default: ${MAX_JOBS})"
    echo "  -f, --force                Force re-analysis of all steps by removing all previous output."
    echo "  -h, --help                 Display this help message."
    exit 1
}

# --- Worker Function (The actual analysis for one sample) ---
run_single_sample() {
    local SAMPLE_NAME="$1"; local R1_FILE="$2"; local R2_FILE="$3"; local OUTPUT_BASE_DIR="$4"
    local THREADS="$5"; local CTAT_LIB_DIR="$6"; local STAR_INDEX_DIR="$7"; local ARRIBA_PATH="$8"
    local ARRIBA_BLACKLIST="$9"; local FUSIONCATCHER_DB="${10}"; local TOOLS_TO_RUN="${11}"; local FORCE="${12}"
    set -euo pipefail

    # Helper function to check if a tool should run
    should_run_tool() {
        local tool="$1"
        local tools_to_run="$2"
        
        if [[ "$tools_to_run" == "all" ]]; then
            return 0
        elif [[ "$tools_to_run" == *"$tool"* ]]; then
            return 0
        else
            return 1
        fi
    }
    
    # --- Define all output and status paths ---
    local GENOME_FASTA="${CTAT_LIB_DIR}/ref_genome.fa"; local ANNOTATION_GTF="${CTAT_LIB_DIR}/ref_annot.gtf"
    local STATUS_DIR="${OUTPUT_BASE_DIR}/status/${SAMPLE_NAME}"
    local TRIMMED_DIR="${OUTPUT_BASE_DIR}/trimmed/${SAMPLE_NAME}"
    local ALIGNED_DIR="${OUTPUT_BASE_DIR}/aligned/${SAMPLE_NAME}"
    local STAR_FUSION_OUT_DIR="${OUTPUT_BASE_DIR}/star_fusion/${SAMPLE_NAME}"
    local ARRIBA_OUT_DIR="${OUTPUT_BASE_DIR}/arriba/${SAMPLE_NAME}"
    local FUSIONCATCHER_OUT_DIR="${OUTPUT_BASE_DIR}/fusioncatcher/${SAMPLE_NAME}"
    local FUSIONCATCHER_INPUT_DIR="${OUTPUT_BASE_DIR}/fusioncatcher_input/${SAMPLE_NAME}"
    
    # Define paths for completion marker files
    local fastp_complete="${STATUS_DIR}/fastp.complete"; local star1_complete="${STATUS_DIR}/star1.complete"
    local star2_complete="${STATUS_DIR}/star2.complete"; local starfusion_complete="${STATUS_DIR}/starfusion.complete"
    local arriba_complete="${STATUS_DIR}/arriba.complete"; local fusioncatcher_complete="${STATUS_DIR}/fusioncatcher.complete"

    if [[ "$FORCE" == "true" ]]; then
        echo "--- [${SAMPLE_NAME}] --force specified. Cleaning all previous output. ---"
        rm -rf "$TRIMMED_DIR" "$ALIGNED_DIR" "$STAR_FUSION_OUT_DIR" "$ARRIBA_OUT_DIR" "$FUSIONCATCHER_OUT_DIR" "$FUSIONCATCHER_INPUT_DIR" "$STATUS_DIR"
    fi
    
    mkdir -p "$STATUS_DIR"; eval "$(conda shell.bash hook)"

    # --- CHECKPOINT 1: fastp Trimming ---
    if [[ ! -f "$fastp_complete" ]]; then
        echo "--- [${SAMPLE_NAME}] Step 1: Cleaning and running fastp... ---"
        rm -rf "$TRIMMED_DIR" && mkdir -p "$TRIMMED_DIR"
        local TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"; local TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"
        conda activate fusion
        fastp \
              -i "$R1_FILE" -I "$R2_FILE" \
              -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
              --thread "$THREADS" \
              --detect_adapter_for_pe \
              --length_required 30 \
              --qualified_quality_phred 20 \
              --html "${TRIMMED_DIR}/fastp_report.html" \
              --json "${TRIMMED_DIR}/fastp_report.json" \
              --report_title "${SAMPLE_NAME} Fastp Report"
        touch "$fastp_complete"
    else
        echo "--- [${SAMPLE_NAME}] Step 1: fastp trimming already complete. Skipping. ---"
    fi

    # --- CHECKPOINT 2.1: First STAR Alignment ---
    if should_run_tool "starfusion" "$TOOLS_TO_RUN"; then
        if [[ ! -f "$star1_complete" ]]; then
            echo "--- [${SAMPLE_NAME}] Step 2.1: Cleaning and running First STAR Alignment... ---"
            rm -rf "$ALIGNED_DIR" && mkdir -p "$ALIGNED_DIR"
            local TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"; local TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"
            conda activate fusion
            STAR --genomeDir "${STAR_INDEX_DIR}" \
                    --runThreadN "${THREADS}" \
                    --readFilesIn "${TRIMMED_R1}" "${TRIMMED_R2}" \
                    --outSAMtype BAM SortedByCoordinate \
                    --outSAMunmapped Within \
                    --outSAMattributes Standard \
                    --readFilesCommand zcat \
                    --outFileNamePrefix "${ALIGNED_DIR}/${SAMPLE_NAME}_" \
                    --outReadsUnmapped None \
                    --twopassMode Basic \
                    --outSAMstrandField intronMotif \
                    --chimSegmentMin 12 \
                    --chimJunctionOverhangMin 8 \
                    --chimOutJunctionFormat 1 \
                    --alignSJDBoverhangMin 10 \
                    --alignMatesGapMax 100000 \
                    --alignIntronMax 100000 \
                    --alignSJstitchMismatchNmax 5 -1 5 5 \
                    --outSAMattrRGline ID:GRPundef \
                    --chimMultimapScoreRange 3 \
                    --chimScoreJunctionNonGTAG -4 \
                    --chimMultimapNmax 20 \
                    --chimNonchimScoreDropMin 10 \
                    --peOverlapNbasesMin 12 \
                    --peOverlapMMp 0.1 \
                    --alignInsertionFlush Right \
                    --alignSplicedMateMapLminOverLmate 0 \
                    --alignSplicedMateMapLmin 30
            touch "$star1_complete"
        else
            echo "--- [${SAMPLE_NAME}] Step 2.1: First STAR Alignment already complete. Skipping. ---"
        fi
    fi

    # --- CHECKPOINT 2.2: Second STAR Alignment (for Arriba) ---
    if should_run_tool "arriba" "$TOOLS_TO_RUN"; then
        if [[ ! -f "$star2_complete" ]]; then
            echo "--- [${SAMPLE_NAME}] Step 2.2: Running Second STAR Alignment (optimized for Arriba)... ---"
            local TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"; local TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"
            conda activate fusion
            STAR \
                  --genomeDir "${STAR_INDEX_DIR}" \
                  --runThreadN "${THREADS}" \
                  --readFilesIn "${TRIMMED_R1}" "${TRIMMED_R2}" \
                  --readFilesCommand zcat \
                  --outFileNamePrefix "${ALIGNED_DIR}/${SAMPLE_NAME}_" \
                  --outSAMtype BAM SortedByCoordinate \
                  --outSAMunmapped Within \
                  --outSAMattributes Standard \
                  --twopassMode Basic \
                  --outSAMstrandField intronMotif \
                  --chimSegmentMin 12 \
                  --chimJunctionOverhangMin 8 \
                  --chimOutJunctionFormat 1 \
                  --alignSJDBoverhangMin 10 \
                  --alignMatesGapMax 100000 \
                  --alignIntronMax 100000 \
                  --alignSJstitchMismatchNmax 5 -1 5 5 \
                  --outSAMattrRGline ID:GRPundef \
                  --chimMultimapScoreRange 3 \
                  --chimScoreJunctionNonGTAG -4 \
                  --chimMultimapNmax 20 \
                  --chimNonchimScoreDropMin 10 \
                  --peOverlapNbasesMin 12 \
                  --peOverlapMMp 0.1 \
                  --alignInsertionFlush Right \
                  --alignSplicedMateMapLminOverLmate 0 \
                  --alignSplicedMateMapLmin 30 \
                  --chimOutType WithinBAM SoftClip
            touch "$star2_complete"
        else
            echo "--- [${SAMPLE_NAME}] Step 2.2: Second STAR Alignment already complete. Skipping. ---"
        fi
    fi

    # --- CHECKPOINT 3: Fusion Detection ---
    if should_run_tool "starfusion" "$TOOLS_TO_RUN"; then
        if [[ ! -f "$starfusion_complete" ]]; then
            echo "--- [${SAMPLE_NAME}]   -> Cleaning and running STAR-Fusion..."
            rm -rf "$STAR_FUSION_OUT_DIR" && mkdir -p "$STAR_FUSION_OUT_DIR"
            local CHIMERIC_JUNCTION="${ALIGNED_DIR}/${SAMPLE_NAME}_Chimeric.out.junction"
            conda activate fusion
            STAR-Fusion --genome_lib_dir "${CTAT_LIB_DIR}" \
                    -J "${CHIMERIC_JUNCTION}" \
                    --output_dir "${STAR_FUSION_OUT_DIR}" --CPU "$THREADS"
            touch "$starfusion_complete"
        else
            echo "--- [${SAMPLE_NAME}]   -> STAR-Fusion already complete. Skipping. ---"
        fi
    fi

    if should_run_tool "arriba" "$TOOLS_TO_RUN"; then
        if [[ ! -f "$arriba_complete" ]]; then
            echo "--- [${SAMPLE_NAME}]   -> Cleaning and running Arriba..."
            rm -rf "$ARRIBA_OUT_DIR" && mkdir -p "$ARRIBA_OUT_DIR"
            local BAM_FILE="${ALIGNED_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"
            conda activate fusion
            "$ARRIBA_PATH" \
                -x "${BAM_FILE}" \
                -a "${GENOME_FASTA}" \
                -g "${ANNOTATION_GTF}" \
                -b "${ARRIBA_BLACKLIST}" \
                -o "${ARRIBA_OUT_DIR}/fusions.tsv" \
                -O "${ARRIBA_OUT_DIR}/fusions.discarded.tsv"
            touch "$arriba_complete"
        else
            echo "--- [${SAMPLE_NAME}]   -> Arriba already complete. Skipping. ---"
        fi
    fi

    if should_run_tool "fusioncatcher" "$TOOLS_TO_RUN"; then
        if [[ ! -f "$fusioncatcher_complete" ]]; then
            echo "--- [${SAMPLE_NAME}]   -> Cleaning and running FusionCatcher..."
            rm -rf "$FUSIONCATCHER_OUT_DIR" "$FUSIONCATCHER_INPUT_DIR" && mkdir -p "$FUSIONCATCHER_OUT_DIR" "$FUSIONCATCHER_INPUT_DIR"
            local TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"; local TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"
            
            # Create symlinks for FusionCatcher input format
            ln -sf "$(realpath "$TRIMMED_R1")" "$FUSIONCATCHER_INPUT_DIR/sample_1.fastq.gz"
            ln -sf "$(realpath "$TRIMMED_R2")" "$FUSIONCATCHER_INPUT_DIR/sample_2.fastq.gz"
            
            conda activate fusioncatcher
            fusioncatcher \
                -d "$FUSIONCATCHER_DB" \
                -i "$FUSIONCATCHER_INPUT_DIR" \
                -o "$FUSIONCATCHER_OUT_DIR" \
                -p "$THREADS"
            touch "$fusioncatcher_complete"
        else
            echo "--- [${SAMPLE_NAME}]   -> FusionCatcher already complete. Skipping. ---"
        fi
    fi
    echo "--- [${SAMPLE_NAME}] Analysis complete! ---"
}
export -f run_single_sample

# --- Main Script Logic ---
INPUT_DIR=""; R1_FILE=""; R2_FILE=""; SAMPLE_NAME=""; REF_BUNDLE_DIR=""
CTAT_LIB_DIR=""; STAR_INDEX_DIR=""; ARRIBA_PATH=""; ARRIBA_BLACKLIST=""; FUSIONCATCHER_DB=""
# local GENOME_BUILD="not_specified" # Commented out as requested

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -I|--input-dir) INPUT_DIR="$2"; shift ;;
        -1|--r1) R1_FILE="$2"; shift ;;
        -2|--r2) R2_FILE="$2"; shift ;;
        -s|--sample) SAMPLE_NAME="$2"; shift ;;
        -o|--output-dir) OUTPUT_BASE_DIR="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        -j|--max-jobs) MAX_JOBS="$2"; shift ;;
        -f|--force) FORCE_RUN="true" ;;
        -R|--ref-bundle-dir) REF_BUNDLE_DIR="$2"; shift ;;
        # --genome-build) GENOME_BUILD="$2"; shift ;; # Commented out as requested
        --run-tools)
            # Validate tool combinations
            IFS=',' read -ra TOOLS <<< "$2"
            for tool in "${TOOLS[@]}"; do
                if [[ "$tool" != "starfusion" && "$tool" != "arriba" && "$tool" != "fusioncatcher" && "$tool" != "all" ]]; then
                    echo "Error: Invalid tool '$tool' in --run-tools. Valid options: starfusion, arriba, fusioncatcher, all, or comma-separated combinations." >&2
                    exit 1
                fi
            done
            RUN_TOOLS="$2"; shift ;;
        --ctat-lib) CTAT_LIB_DIR="$2"; shift ;;
        --star-index) STAR_INDEX_DIR="$2"; shift ;;
        --arriba-path) ARRIBA_PATH="$2"; shift ;;
        --arriba-blacklist) ARRIBA_BLACKLIST="$2"; shift ;;
        --fusioncatcher-db) FUSIONCATCHER_DB="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

if [[ -n "$REF_BUNDLE_DIR" ]]; then
    CTAT_LIB_DIR=${CTAT_LIB_DIR:-$(find "$REF_BUNDLE_DIR" -maxdepth 2 -name "ctat_genome_lib_build_dir" -type d | head -n 1)}
    STAR_INDEX_DIR=${STAR_INDEX_DIR:-$(find "$REF_BUNDLE_DIR" -maxdepth 2 -name "star_index" -type d | head -n 1)}
    ARRIBA_PATH=${ARRIBA_PATH:-$(find "$REF_BUNDLE_DIR" -name "arriba" -type f -executable | head -n 1)}
    ARRIBA_BLACKLIST=${ARRIBA_BLACKLIST:-$(find "$REF_BUNDLE_DIR" -name "blacklist*.tsv.gz" -type f | head -n 1)}
    FUSIONCATCHER_DB=${FUSIONCATCHER_DB:-$(find -L "$REF_BUNDLE_DIR" -type d -path '*/fusioncatcher/data
/current')}
fi

# Updated validation to handle tool-specific requirements
validate_tools_and_paths() {
    local tools_to_run="$1"
    local need_star_arriba=false
    local need_fusioncatcher=false
    
    # Determine what tools are needed
    if [[ "$tools_to_run" == "all" ]]; then
        need_star_arriba=true
        need_fusioncatcher=true
    else
        if should_run_tool "starfusion" "$tools_to_run" || should_run_tool "arriba" "$tools_to_run"; then
            need_star_arriba=true
        fi
        if should_run_tool "fusioncatcher" "$tools_to_run"; then
            need_fusioncatcher=true
        fi
    fi
    
    # Validate paths based on needed tools
    if [[ "$need_star_arriba" == "true" ]]; then
        if [[ -z "$OUTPUT_BASE_DIR" || -z "$CTAT_LIB_DIR" || -z "$STAR_INDEX_DIR" ]]; then
            echo "Error: Missing required paths for STAR-Fusion/Arriba (--ctat-lib, --star-index)." >&2
            usage
        fi
        if should_run_tool "arriba" "$tools_to_run" && [[ -z "$ARRIBA_PATH" || -z "$ARRIBA_BLACKLIST" ]]; then
            echo "Error: Missing required paths for Arriba (--arriba-path, --arriba-blacklist)." >&2
            usage
        fi
    fi
    
    if [[ "$need_fusioncatcher" == "true" && -z "$FUSIONCATCHER_DB" ]]; then
        echo "Error: Missing required path for FusionCatcher (--fusioncatcher-db)." >&2
        usage
    fi
    
    if [[ -z "$OUTPUT_BASE_DIR" ]]; then
        echo "Error: Missing required --output-dir." >&2
        usage
    fi
}

validate_tools_and_paths "$RUN_TOOLS"

if [[ -n "$INPUT_DIR" ]]; then # BATCH MODE
    echo "🔄 Batch Mode Detected. Tools to run: ${RUN_TOOLS}"
    for r1 in "$INPUT_DIR"/*_R1*.fastq.gz; do
        sample=$(basename "$r1" | sed -e 's/_R1.*.fastq.gz//' -e 's/_R1.*.fq.gz//')
        
        # Quick check for final completion markers to skip fully completed samples early
        sample_complete=false
        status_dir="${OUTPUT_BASE_DIR}/status/${sample}"
        starfusion_complete="${status_dir}/starfusion.complete"
        arriba_complete="${status_dir}/arriba.complete"
        fusioncatcher_complete="${status_dir}/fusioncatcher.complete"
        
        # Check if all requested tools are complete
        if [[ "$RUN_TOOLS" == "all" ]]; then
            if [ -f "$starfusion_complete" ] && [ -f "$arriba_complete" ] && [ -f "$fusioncatcher_complete" ]; then
                sample_complete=true
            fi
        else
            sample_complete=true
            if should_run_tool "starfusion" "$RUN_TOOLS" && [ ! -f "$starfusion_complete" ]; then
                sample_complete=false
            fi
            if should_run_tool "arriba" "$RUN_TOOLS" && [ ! -f "$arriba_complete" ]; then
                sample_complete=false
            fi
            if should_run_tool "fusioncatcher" "$RUN_TOOLS" && [ ! -f "$fusioncatcher_complete" ]; then
                sample_complete=false
            fi
        fi

        if [[ "$FORCE_RUN" != "true" && "$sample_complete" == "true" ]]; then
            echo "🔄 Skipping sample: $sample (final completion markers exist)."
            continue
        fi
        
        while (($(jobs -r | wc -l) >= MAX_JOBS)); do sleep 5; done
        r2=$(find "$INPUT_DIR" -name "${sample}_R2*.fastq.gz" -o -name "${sample}_R2*.fq.gz" | head -n 1)
        if [[ -f "$r1" && -f "$r2" ]]; then
            echo "🔄 Launching job for sample: $sample"
            log_file="${OUTPUT_BASE_DIR}/logs/${sample}.log"; mkdir -p "$(dirname "$log_file")"
            run_single_sample "$sample" "$r1" "$r2" "$OUTPUT_BASE_DIR" "$THREADS" "$CTAT_LIB_DIR" "$STAR_INDEX_DIR" "$ARRIBA_PATH" "$ARRIBA_BLACKLIST" "$FUSIONCATCHER_DB" "$RUN_TOOLS" "$FORCE_RUN" > "$log_file" 2>&1 &
        else
             echo "🔄 Warning: Skipping $sample as its R2 pair was not found."
        fi
    done
    echo "🔄 All jobs launched. Waiting for completion..."
    wait
elif [[ -n "$R1_FILE" ]]; then # SINGLE SAMPLE MODE
    if [[ -z "$SAMPLE_NAME" || -z "$R2_FILE" ]]; then echo "Error: --r2 and --sample are required." >&2; usage; fi
    echo "🔄 Single Sample Mode Detected. Tools to run: ${RUN_TOOLS}"
    run_single_sample "$SAMPLE_NAME" "$R1_FILE" "$R2_FILE" "$OUTPUT_BASE_DIR" "$THREADS" "$CTAT_LIB_DIR" "$STAR_INDEX_DIR" "$ARRIBA_PATH" "$ARRIBA_BLACKLIST" "$FUSIONCATCHER_DB" "$RUN_TOOLS" "$FORCE_RUN"
else
    echo "Error: No input specified. Use --input-dir or --r1/--r2." >&2; usage
fi

echo "✅ All processing complete!"
