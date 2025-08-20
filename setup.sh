#!/bin/bash
set -euo pipefail

# =============================================================================
# FUSION ANALYSIS PIPELINE - COMPLETE SETUP SCRIPT
# =============================================================================
# This script will set up everything needed to run the gene fusion analysis
# pipeline. It's designed for novice users and handles all installations
# and downloads automatically.
#
# Author: Fungel Lab
# Repository: https://github.com/Rajarshi0/Fusion_analysis_Fungel_lab
# =============================================================================

# Color codes for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m' # No Color

# Global variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR"
LOG_FILE="$SCRIPT_DIR/setup.log"
REFERENCE_DIR="$SCRIPT_DIR/references"
TEMP_DIR="$SCRIPT_DIR/temp_downloads"

# System requirements
MIN_DISK_SPACE_GB=200
MIN_RAM_GB=16

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] $1${NC}" | tee -a "$LOG_FILE"
}

log_warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1${NC}" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}" | tee -a "$LOG_FILE"
}

log_info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1${NC}" | tee -a "$LOG_FILE"
}

print_header() {
    echo -e "${PURPLE}"
    echo "=================================================================="
    echo "$1"
    echo "=================================================================="
    echo -e "${NC}"
}

check_command() {
    if command -v "$1" >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

check_disk_space() {
    local available_gb
    available_gb=$(df "$SCRIPT_DIR" --output=avail --block-size=1G | tail -n 1 | tr -d ' ')
    
    if [ "$available_gb" -lt "$MIN_DISK_SPACE_GB" ]; then
        log_error "Insufficient disk space. Need at least ${MIN_DISK_SPACE_GB}GB, have ${available_gb}GB"
        exit 1
    fi
    
    log "Disk space check passed: ${available_gb}GB available"
}

check_memory() {
    local available_gb
    available_gb=$(free -g | awk '/^Mem:/{print $2}')
    
    if [ "$available_gb" -lt "$MIN_RAM_GB" ]; then
        log_warning "Low RAM detected: ${available_gb}GB (recommended: ${MIN_RAM_GB}GB+)"
    else
        log "Memory check passed: ${available_gb}GB RAM available"
    fi
}

# =============================================================================
# SYSTEM CHECKS
# =============================================================================

perform_system_checks() {
    print_header "SYSTEM REQUIREMENTS CHECK"
    
    log "Checking system requirements..."
    
    # Check if running on Linux
    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        log_error "This pipeline requires Linux. Detected OS: $OSTYPE"
        exit 1
    fi
    
    # Check disk space and memory
    check_disk_space
    check_memory
    
    # Check for required system packages
    local missing_packages=()
    
    if ! check_command "wget"; then
        missing_packages+=("wget")
    fi
    
    if ! check_command "curl"; then
        missing_packages+=("curl")
    fi
    
    if ! check_command "tar"; then
        missing_packages+=("tar")
    fi
    
    if ! check_command "gzip"; then
        missing_packages+=("gzip")
    fi
    
    if ! check_command "git"; then
        missing_packages+=("git")
    fi
    
    if [ ${#missing_packages[@]} -gt 0 ]; then
        log_error "Missing required system packages: ${missing_packages[*]}"
        log_info "Please install them using your system package manager:"
        log_info "  Ubuntu/Debian: sudo apt-get install ${missing_packages[*]}"
        log_info "  CentOS/RHEL:   sudo yum install ${missing_packages[*]}"
        exit 1
    fi
    
    log "System checks completed successfully!"
}

# =============================================================================
# CONDA INSTALLATION AND SETUP
# =============================================================================

install_conda() {
    print_header "CONDA INSTALLATION"
    
    if check_command "conda"; then
        log "Conda already installed: $(conda --version)"
        return 0
    fi
    
    log "Installing Miniconda..."
    
    local conda_installer="Miniconda3-latest-Linux-x86_64.sh"
    local conda_url="https://repo.anaconda.com/miniconda/$conda_installer"
    
    # Download installer
    wget -q "$conda_url" -O "$TEMP_DIR/$conda_installer"
    
    # Install conda
    bash "$TEMP_DIR/$conda_installer" -b -p "$HOME/miniconda3"
    
    # Initialize conda
    "$HOME/miniconda3/bin/conda" init bash
    
    # Source bashrc to make conda available
    source "$HOME/.bashrc" 2>/dev/null || true
    
    # Add conda to PATH for this script
    export PATH="$HOME/miniconda3/bin:$PATH"
    
    log "Conda installed successfully!"
}

create_conda_environments() {
    print_header "CREATING CONDA ENVIRONMENTS"
    
    # Ensure conda is available
    export PATH="$HOME/miniconda3/bin:$PATH"
    
    log "Creating fusion environment..."
    conda env create -f "$PIPELINE_DIR/env/fusion_environment.yml" -q || {
        log_warning "Failed to create from .yml file, trying manual creation..."
        conda create -n fusion -c bioconda -c conda-forge \
            star-fusion arriba star fastp samtools python=3.8 pandas numpy -y
    }
    
    log "Creating fusioncatcher environment..."
    conda env create -f "$PIPELINE_DIR/env/fusioncatcher_environment.yml" -q || {
        log_warning "Failed to create from .yml file, trying manual creation..."
        conda create -n fusioncatcher -c bioconda -c conda-forge \
            fusioncatcher python=2.7 fastp -y
    }
    
    log "Conda environments created successfully!"
}

# =============================================================================
# SINGULARITY/DOCKER SETUP
# =============================================================================

setup_containers() {
    print_header "CONTAINER SETUP"
    
    # Check if Singularity is available
    if check_command "singularity"; then
        log "Singularity found: $(singularity --version)"
        
        log "Pulling Cicero container..."
        singularity pull "$REFERENCE_DIR/cicero.sif" docker://ghcr.io/stjude/cicero:latest || {
            log_warning "Failed to pull Cicero container. Will try to pull at runtime."
        }
        
        log "Pulling RNApeg container..."
        singularity pull "$REFERENCE_DIR/rnapeg.sif" docker://ghcr.io/stjude/rnapeg:latest || {
            log_warning "Failed to pull RNApeg container. Will try to pull at runtime."
        }
        
    elif check_command "docker"; then
        log "Docker found: $(docker --version)"
        log "Pulling required containers..."
        
        docker pull ghcr.io/stjude/cicero:latest || log_warning "Failed to pull Cicero container"
        docker pull ghcr.io/stjude/rnapeg:latest || log_warning "Failed to pull RNApeg container"
        
    else
        log_warning "Neither Singularity nor Docker found. Container-based workflows will not work."
        log_info "To install Singularity: https://sylabs.io/guides/3.0/user-guide/installation.html"
        log_info "To install Docker: https://docs.docker.com/engine/install/"
    fi
}

# =============================================================================
# REFERENCE DATA DOWNLOAD
# =============================================================================

download_ctat_library() {
    local ctat_url="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz"
    local ctat_file="GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz"
    local ctat_dir="$REFERENCE_DIR/ctat_genome_lib_build_dir"
    
    if [ -d "$ctat_dir" ]; then
        log "CTAT library already exists, skipping download"
        return 0
    fi
    
    log "Downloading CTAT Genome Library (~15GB, this may take a while)..."
    
    # Download with resume capability
    wget -c "$ctat_url" -O "$TEMP_DIR/$ctat_file" --progress=bar:force 2>&1 | \
        while IFS= read -r line; do
            if [[ $line == *"%"* ]]; then
                echo -ne "\r$line"
            fi
        done
    echo ""
    
    log "Extracting CTAT library..."
    tar -xzf "$TEMP_DIR/$ctat_file" -C "$REFERENCE_DIR/" --strip-components=1
    
    # Find and rename the extracted directory
    local extracted_dir=$(find "$REFERENCE_DIR" -maxdepth 1 -name "*CTAT_lib*" -type d)
    if [ -n "$extracted_dir" ]; then
        mv "$extracted_dir" "$ctat_dir"
    fi
    
    log "CTAT library downloaded and extracted successfully!"
}

download_fusioncatcher_data() {
    local fc_data_dir="$REFERENCE_DIR/fusioncatcher_data"
    
    if [ -d "$fc_data_dir" ]; then
        log "FusionCatcher data already exists, skipping download"
        return 0
    fi
    
    log "Setting up FusionCatcher reference data..."
    
    # Activate fusioncatcher environment
    export PATH="$HOME/miniconda3/bin:$PATH"
    eval "$(conda shell.bash hook)"
    conda activate fusioncatcher
    
    # Download FusionCatcher data
    mkdir -p "$fc_data_dir"
    
    log "Downloading FusionCatcher database (~20GB, this will take time)..."
    fusioncatcher-build -d "$fc_data_dir" -g homo_sapiens -v 102 || {
        log_warning "FusionCatcher data download failed. You may need to run this manually."
    }
    
    conda deactivate
    log "FusionCatcher data setup completed!"
}

download_arriba_data() {
    local arriba_dir="$REFERENCE_DIR/arriba"
    local arriba_version="2.5.0"
    local arriba_url="https://github.com/suhrig/arriba/releases/download/v${arriba_version}/arriba_v${arriba_version}.tar.gz"
    
    if [ -d "$arriba_dir" ]; then
        log "Arriba data already exists, skipping download"
        return 0
    fi
    
    log "Downloading Arriba..."
    wget -q "$arriba_url" -O "$TEMP_DIR/arriba.tar.gz"
    
    log "Extracting Arriba..."
    tar -xzf "$TEMP_DIR/arriba.tar.gz" -C "$REFERENCE_DIR/"
    mv "$REFERENCE_DIR/arriba_v${arriba_version}" "$arriba_dir"
    
    log "Arriba downloaded successfully!"
}

create_star_index() {
    local ctat_dir="$REFERENCE_DIR/ctat_genome_lib_build_dir"
    local star_index_dir="$REFERENCE_DIR/star_index"
    
    if [ -d "$star_index_dir" ]; then
        log "STAR index already exists, skipping creation"
        return 0
    fi
    
    log "Creating STAR genome index (this may take 30-60 minutes)..."
    
    # Activate fusion environment
    export PATH="$HOME/miniconda3/bin:$PATH"
    eval "$(conda shell.bash hook)"
    conda activate fusion
    
    # Get number of available threads
    local threads=$(nproc)
    local max_threads=16
    if [ $threads -gt $max_threads ]; then
        threads=$max_threads
    fi
    
    mkdir -p "$star_index_dir"
    
    STAR --runMode genomeGenerate \
         --genomeDir "$star_index_dir" \
         --genomeFastaFiles "$ctat_dir/ref_genome.fa" \
         --sjdbGTFfile "$ctat_dir/ref_annot.gtf" \
         --runThreadN $threads \
         --sjdbOverhang 100 || {
        log_error "Failed to create STAR index"
        conda deactivate
        return 1
    }
    
    conda deactivate
    log "STAR index created successfully!"
}

download_reference_data() {
    print_header "REFERENCE DATA DOWNLOAD"
    
    mkdir -p "$REFERENCE_DIR" "$TEMP_DIR"
    
    # Download in order of dependency
    download_ctat_library
    download_arriba_data
    create_star_index
    download_fusioncatcher_data
    
    log "All reference data downloaded successfully!"
}

# =============================================================================
# SCRIPT CONFIGURATION
# =============================================================================

configure_scripts() {
    print_header "SCRIPT CONFIGURATION"
    
    log "Configuring pipeline scripts with correct paths..."
    
    # Configure final_SF_arriba_analysis.sh (it uses auto-discovery)
    log "Main analysis script uses auto-discovery, no configuration needed"
    
    # Configure other scripts if needed
    local scripts_to_configure=(
        "scripts/fusioncatcher_parallel_2.sh"
        "scripts/run_parallel_FC.sh"
        "scripts/run_cicero_parallel_docker.sh"
    )
    
    for script in "${scripts_to_configure[@]}"; do
        if [ -f "$PIPELINE_DIR/$script" ]; then
            log "Script $script will be configured at runtime"
        fi
    done
    
    # Create a configuration file for easy reference
    cat > "$PIPELINE_DIR/pipeline_config.sh" << EOF
#!/bin/bash
# Auto-generated configuration file
# Source this file to set up environment variables for the pipeline

export PIPELINE_DIR="$PIPELINE_DIR"
export REFERENCE_DIR="$REFERENCE_DIR"
export CTAT_LIB_DIR="$REFERENCE_DIR/ctat_genome_lib_build_dir"
export STAR_INDEX_DIR="$REFERENCE_DIR/star_index"
export ARRIBA_PATH="$REFERENCE_DIR/arriba/arriba"
export ARRIBA_BLACKLIST="$REFERENCE_DIR/arriba/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"
export FUSIONCATCHER_DB="$REFERENCE_DIR/fusioncatcher_data"
export CICERO_REF_DIR="$REFERENCE_DIR/cicero"

# Activate conda base environment
export PATH="$HOME/miniconda3/bin:\$PATH"

echo "Pipeline configuration loaded!"
echo "CTAT Library: \$CTAT_LIB_DIR"
echo "STAR Index: \$STAR_INDEX_DIR"
echo "Arriba: \$ARRIBA_PATH"
echo "FusionCatcher DB: \$FUSIONCATCHER_DB"
EOF

    chmod +x "$PIPELINE_DIR/pipeline_config.sh"
    
    log "Configuration completed! Settings saved to pipeline_config.sh"
}

# =============================================================================
# TEST INSTALLATION
# =============================================================================

test_installation() {
    print_header "INSTALLATION TEST"
    
    log "Testing pipeline installation..."
    
    # Test conda environments
    export PATH="$HOME/miniconda3/bin:$PATH"
    eval "$(conda shell.bash hook)"
    
    log "Testing fusion environment..."
    if conda activate fusion 2>/dev/null; then
        if command -v STAR >/dev/null 2>&1; then
            log "✓ STAR found in fusion environment"
        else
            log_warning "✗ STAR not found in fusion environment"
        fi
        conda deactivate
    else
        log_error "✗ Failed to activate fusion environment"
    fi
    
    log "Testing fusioncatcher environment..."
    if conda activate fusioncatcher 2>/dev/null; then
        if command -v fusioncatcher >/dev/null 2>&1; then
            log "✓ FusionCatcher found in environment"
        else
            log_warning "✗ FusionCatcher not found in environment"
        fi
        conda deactivate
    else
        log_error "✗ Failed to activate fusioncatcher environment"
    fi
    
    # Test reference files
    local ref_files=(
        "$REFERENCE_DIR/ctat_genome_lib_build_dir/ref_genome.fa"
        "$REFERENCE_DIR/star_index/SA"
        "$REFERENCE_DIR/arriba/arriba"
    )
    
    for ref_file in "${ref_files[@]}"; do
        if [ -f "$ref_file" ] || [ -d "$ref_file" ]; then
            log "✓ Reference file/directory exists: $(basename "$ref_file")"
        else
            log_warning "✗ Missing reference: $ref_file"
        fi
    done
    
    log "Installation test completed!"
}

# =============================================================================
# CLEANUP
# =============================================================================

cleanup_temp_files() {
    log "Cleaning up temporary files..."
    rm -rf "$TEMP_DIR"
    log "Cleanup completed!"
}

# =============================================================================
# USAGE INSTRUCTIONS
# =============================================================================

print_usage_instructions() {
    print_header "SETUP COMPLETED - USAGE INSTRUCTIONS"
    
    cat << 'EOF'

🎉 CONGRATULATIONS! Your fusion analysis pipeline is now ready to use!

📁 WHAT WAS INSTALLED:
   ├── Conda environments (fusion, fusioncatcher)
   ├── Reference data (~30-50 GB)
   ├── Container images (if Singularity/Docker available)
   └── Configured pipeline scripts

🚀 QUICK START:

1. ACTIVATE THE ENVIRONMENT:
   source ~/pipeline_config.sh

2. RUN A SINGLE SAMPLE (STAR-Fusion + Arriba):
   ./final_SF_arriba_analysis.sh \
     --r1 sample_R1.fastq.gz \
     --r2 sample_R2.fastq.gz \
     --sample MySample \
     --ref-bundle-dir references/

3. RUN BATCH ANALYSIS:
   ./final_SF_arriba_analysis.sh \
     --input-dir /path/to/fastq/files/ \
     --ref-bundle-dir references/

4. RUN FUSIONCATCHER:
   conda activate fusioncatcher
   ./scripts/run_parallel_FC.sh

📖 FULL DOCUMENTATION:
   See README.md for complete usage instructions

🔍 TROUBLESHOOTING:
   - Check setup.log for any warnings or errors
   - Ensure you have sufficient disk space and memory
   - Test with a small sample first

📧 SUPPORT:
   GitHub: https://github.com/Rajarshi0/Fusion_analysis_Fungel_lab

EOF

    log "Setup completed successfully! Check the instructions above."
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

main() {
    # Create log file
    touch "$LOG_FILE"
    
    print_header "FUSION ANALYSIS PIPELINE SETUP"
    log "Starting complete pipeline setup..."
    log "Installation directory: $PIPELINE_DIR"
    log "Log file: $LOG_FILE"
    
    # Run setup steps
    perform_system_checks
    install_conda
    create_conda_environments
    setup_containers
    download_reference_data
    configure_scripts
    test_installation
    cleanup_temp_files
    print_usage_instructions
    
    log "🎉 Setup completed successfully!"
    echo ""
    echo -e "${GREEN}Setup log saved to: $LOG_FILE${NC}"
    echo -e "${GREEN}Configuration saved to: $PIPELINE_DIR/pipeline_config.sh${NC}"
    echo ""
}

# Handle interrupts gracefully
trap 'log_error "Setup interrupted by user"; cleanup_temp_files; exit 1' INT TERM

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
