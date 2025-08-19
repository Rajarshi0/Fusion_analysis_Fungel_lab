#!/bin/bash
set -euo pipefail

# --- INPUTS ---
# Usage: ./run_fusion_analysis.sh /path/to/reads_R1.fastq.gz /path/to/reads_R2.fastq.gz
R1_FILE="$1"
R2_FILE="$2"

# --- CONFIGURATION ---
THREADS=16
OUTPUT_BASE_DIR="fusion_analysis_output"

# Tool and Reference Paths
CTAT_LIB_DIR="../existing_scripts/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"
STAR_INDEX_DIR="../existing_scripts/star_index"
ARRIBA_PATH="../arriba_v1.2.0/arriba"


# --- DERIVED VARIABLES ---
SAMPLE_NAME=$(basename "$R1_FILE" _R1_001.fastq.gz)
GENOME_FASTA="${CTAT_LIB_DIR}/ref_genome.fa"
ANNOTATION_GTF="${CTAT_LIB_DIR}/ref_annot.gtf"
ARRIBA_BLACKLIST="../arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"

# Output Directories
TRIMMED_DIR="${OUTPUT_BASE_DIR}/trimmed/${SAMPLE_NAME}"
ALIGNED_DIR="${OUTPUT_BASE_DIR}/aligned/${SAMPLE_NAME}"
STAR_FUSION_OUT_DIR="${OUTPUT_BASE_DIR}/star_fusion/${SAMPLE_NAME}"
ARRIBA_OUT_DIR="${OUTPUT_BASE_DIR}/arriba/${SAMPLE_NAME}"
BAM_FILE="${ALIGNED_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"

echo "================================================================="
echo "? Processing Sample: ${SAMPLE_NAME}"
echo "================================================================="

mkdir -p "$TRIMMED_DIR" "$ALIGNED_DIR" "$STAR_FUSION_OUT_DIR" "$ARRIBA_OUT_DIR"

# --- Step 1: Trimming with fastp ---
echo "?? Step 1: Trimming reads with fastp..."
eval "$(conda shell.bash hook)"
conda activate fusion

TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fq.gz"
TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fq.gz"

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

# --- Step 2: STAR Alignment (Two Runs as Requested) ---
echo "?? Step 2.1: First STAR Alignment Run..."
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

echo "?? Step 2.2: Second STAR Alignment Run..."
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

# --- Step 3: Fusion Detection ---
echo "?? Step 3: Detecting fusions with STAR-Fusion and Arriba..."

# Run STAR-Fusion
STAR-Fusion --genome_lib_dir "${CTAT_LIB_DIR}" \
        -J "${ALIGNED_DIR}/${SAMPLE_NAME}_Chimeric.out.junction" \
        --output_dir "${STAR_FUSION_OUT_DIR}" --CPU "$THREADS"

# Run Arriba
"$ARRIBA_PATH" \
    -x "${BAM_FILE}" \
    -a "${GENOME_FASTA}" \
    -g "${ANNOTATION_GTF}" \
    -b "${ARRIBA_BLACKLIST}" \
    -o "${ARRIBA_OUT_DIR}/fusions.tsv" \
    -O "${ARRIBA_OUT_DIR}/fusions.discarded.tsv"

echo "? Analysis complete for ${SAMPLE_NAME}!"
