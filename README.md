# 🧬 High-Throughput Gene Fusion Analysis Pipeline – Fungel Lab

Welcome! This repository contains a complete suite of scripts for performing high-throughput gene fusion analysis on RNA-seq data. It includes parallel execution wrappers for multiple fusion callers (**STAR-Fusion**, **Arriba**, and **FusionCatcher**), Conda environment files for reproducible setup, and Python scripts for post-processing results into a standardized format.

This is created using the help of AI like Gemini and Claude...

---

## 📌 Table of Contents
- [Features](#-features)
- [System Requirements](#-system-requirements)
- [Setup & Installation](#-setup--installation)
- [A Beginner's Guide to Parallel Execution & Resources](#-a-beginners-guide-to-parallel-execution--resources)
- [Analysis Workflows](#-analysis-workflows)
- [How to Run an Analysis](#-how-to-run-an-analysis)
- [Script Descriptions](#-script-descriptions)
- [Output Structure](#-output-structure)
- [Post-Processing & Format Conversion](#-post-processing--format-conversion)
- [Troubleshooting](#-troubleshooting)
- [Contact](#-contact)

---

## ✨ Features

### **Multiple Fusion Calling Pipelines**
- **STAR-Fusion + Arriba**: Combined pipeline with optimized STAR alignment parameters
- **FusionCatcher**: Independent pipeline with Python 2.7 environment


### **Advanced Execution Management** 
- **Modular Checkpointing**: Resume analysis from any failed step without losing progress
- **Smart Resource Management**: Automatic memory monitoring and job throttling
- **Parallel Processing**: Sophisticated wrapper scripts manage concurrent jobs efficiently
- **Self-Cleaning**: Automatic cleanup of old/failed processes before starting new runs

### **Reproducible Environments**
- Pre-configured Conda environments with explicit package lists
- Docker/Singularity integration for containerized workflows
- Comprehensive dependency management

### **Automated Post-Processing**
- Python scripts for converting outputs to the standardized Cicero format
- Batch processing capabilities with automatic file discovery
- Quality control and error handling

---

## 🛠 System Requirements

### **Operating System**
- Linux (tested on CentOS/Ubuntu)

### **Software Dependencies**
- [Conda/Miniconda](https://docs.conda.io/en/latest/miniconda.html) for environment management


### **Hardware Recommendations**
- **CPU**: Multi-core processor (16+ cores recommended)
- **RAM**: 64GB minimum, 128GB+ recommended for large datasets
- **Storage**: High-speed SSD with sufficient space for reference data and outputs

### **Reference Data Requirements**
- STAR genome indices
- CTAT Genome Library
- FusionCatcher database
- Arriba blacklist and reference files


---

## ⚙️ Setup & Installation

Majorly, **setup.sh** should handle the complete setup and download of libraries and the database.

### 1. **Clone the Repository**
```bash
git clone https://github.com/Rajarshi0/Fusion_analysis_Fungel_lab.git
cd Fusion_analysis_Fungel_lab
```

### 2. **Create Conda Environments**

**For STAR-Fusion, Arriba, and Python 3 post-processing:**
```bash
conda env create -f env/fusion_environment.yml
```

**For FusionCatcher (requires Python 2.7):**
```bash
conda env create -f env/fusioncatcher_environment.yml
```

### 3. **Prepare Reference Data**
Download and configure the following reference datasets:
- **CTAT Genome Library**: For STAR-Fusion
  ```
  ## Site link: https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
  wget -c https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
  tar -xvzf GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
  ```
- **STAR Index**: Pre-built genome index
   ```
  ## refer STAR github: 
  https://github.com/alexdobin/STAR/
  ```
- **FusionCatcher Database**: Species-specific fusion database
  ```
  ## refer fusioincatcher github for buiding the database:
  https://github.com/ndaniel/fusioncatcher
  ```
- **Arriba Files**: Blacklist and reference genome
  
   ```
  ## Arriba release site: https://github.com/suhrig/arriba/releases
  #v2.5.0 contains blacklisted, you can use it if it supports. 
  #I have used the #v1.2.0
  
  wget -c https://github.com/suhrig/arriba/releases/download/v1.2.0/arriba_v1.2.0.tar.gz

  wget -c https://github.com/suhrig/arriba/releases/download/v2.5.0/arriba_v2.5.0.tar.gz
   ```
  
---

## 🧠 A Beginner's Guide to Parallel Execution & Resources

### Understanding Your System's Power

Think of your computer's CPU as a kitchen:
- **Cores** = Individual chefs working simultaneously
- **Threads** = Hands each chef has (usually 2 per core)
- **Total Threads** = Your system's maximum processing power

**Check your system's capacity:**
```bash
# Quick check - shows total threads available
nproc

# Detailed breakdown
lscpu
```

### Configuring Resource Usage

The pipeline uses two key settings:
- `THREADS_PER_JOB`: Threads allocated to each sample analysis
- `MAX_PARALLEL_JOBS`: Number of samples processed simultaneously

**Golden Rule:**
```
THREADS_PER_JOB × MAX_PARALLEL_JOBS ≤ Total System Threads
```

**Example Configuration:**
- System has 64 threads (`nproc` returns 64)
- Option 1: `THREADS_PER_JOB=16`, `MAX_PARALLEL_JOBS=4` (intensive per sample)
- Option 2: `THREADS_PER_JOB=8`, `MAX_PARALLEL_JOBS=8` (more samples simultaneously)

### Memory Management
Set `MIN_FREE_MEMORY_GB` to 10-20% of your total RAM as a safety buffer. The scripts will automatically pause new jobs when memory runs low.

---

## 🔬 Analysis Workflows

### **Workflow 1: STAR-Fusion**
```
Raw FASTQs → fastp trimming → STAR alignment → STAR-Fusion → Results
```

### **Workflow 2: Arriba**
```
Raw FASTQs → fastp trimming → STAR alignment → Arriba → Results
```

### **Workflow 3: FusionCatcher**
```
Raw FASTQs → fastp trimming → FusionCatcher → Results
```


---

## 🚀 How to Run an Analysis

### **Option A: Single Sample Mode**

**STAR-Fusion + Arriba:**
```bash
bash final_SF_arriba_analysis.sh \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --sample MySample \
  --ref-bundle-dir /path/to/references \
  --threads 16
```

### **Option B: Batch Processing Mode**

**STAR-Fusion + Arriba Batch:**
```bash
bash final_SF_arriba_analysis.sh \
  -I /path/to/fastq/directory \
  -R /path/to/references \
  --threads 16 \
  --max-jobs 4
```

**FusionCatcher Batch:**
```bash
# Place script in directory with FASTQ files
./scripts/run_parallel_FC.sh
```

### **Command Line Options (final_SF_arriba_analysis.sh)**

#### **Input Modes**
```bash
# Single Sample Mode
-1, --r1 <file>          R1 FASTQ file
-2, --r2 <file>          R2 FASTQ file  
-s, --sample <name>      Sample name

# Batch Mode
-I, --input-dir <dir>    Directory with FASTQ pairs
```

#### **Reference Configuration**
```bash
-R, --ref-bundle-dir <dir>    Directory with all reference files
--ctat-lib <dir>              CTAT Genome Library path
--star-index <dir>            STAR index directory  
--arriba-path <path>          Arriba executable path
--arriba-blacklist <file>     Arriba blacklist file
```

#### **Analysis Control**
```bash
-o, --output-dir <dir>        Output directory (default: fusion_analysis_output)
-t, --threads <int>           Threads per job (default: 16)
-j, --max-jobs <int>          Max parallel jobs (default: 2)
--run-tools <tool>            Tools to run: starfusion, arriba, both (default: both)
-f, --force                   Force reanalysis by removing previous outputs
```

### **Advanced Usage Examples**

**Run only STAR-Fusion on a batch:**
```bash
bash final_SF_arriba_analysis.sh \
  --input-dir fastq_files/ \
  --ref-bundle-dir /data/references/ \
  --run-tools starfusion \
  --max-jobs 6
```

**Run only multi-tools on a batch:**
```bash
bash final_SF_arriba_analysis.sh \
  --input-dir fastq_files/ \
  --ref-bundle-dir /data/references/ \
  --run-tools starfusion,arriba \
  --max-jobs 6
```

**Resume failed analysis:**
```bash
# The script automatically detects completed steps and resumes from failures
bash final_SF_arriba_analysis.sh --input-dir fastq_files/ --ref-bundle-dir /data/references/
```

**Force complete reanalysis:**
```bash
bash final_SF_arriba_analysis.sh \
  --input-dir fastq_files/ \
  --ref-bundle-dir /data/references/ \
  --force
```

---

## 📂 Output Structure

```
fusion_analysis_output/
├── status/                    # Checkpoint markers
│   └── sample_name/
│       ├── fastp.complete
│       ├── star1.complete
│       ├── star2.complete
│       ├── starfusion.complete
│       └── arriba.complete
├── trimmed/                   # fastp outputs
│   └── sample_name/
├── aligned/                   # STAR alignment outputs
│   └── sample_name/
├── star_fusion/               # STAR-Fusion results
│   └── sample_name/
│       └── star-fusion.fusion_predictions.tsv
├── arriba/                    # Arriba results
│   └── sample_name/
│       ├── fusions.tsv
│       └── fusions.discarded.tsv
└── logs/                      # Execution logs
    └── sample_name.log
```

---

## 🔄 Post-Processing & Format Conversion

### **Convert STAR-Fusion to Cicero Format**
```bash
# Edit paths in the script
python scripts/starfusion2cicero.py
```

### **Convert Arriba to Cicero Format**
```bash
# Edit paths in the script  
python scripts/automated_convert_arriba2cicero.py
```

### **Configuration for Conversion Scripts**
Edit these variables in the Python scripts:
```python
# Directory containing sample folders
SEARCH_DIRECTORY = "fusion_analysis_output/star_fusion"  # or "arriba"

# Output directory for converted files
OUTPUT_DIRECTORY = "cicero_converted_outputs"
```

---

## 🔧 Script Descriptions

### **Main Analysis Scripts**

**`final_SF_arriba_analysis.sh`**
- **Primary analysis pipeline** with advanced features
- Modular checkpointing system
- Single sample and batch processing modes
- Automatic resource management
- Self-cleaning functionality

**`scripts/fusioncatcher_parallel_2.sh`**
- FusionCatcher analysis for single samples
- Integrated with Python 2.7 environment
- fastp trimming + FusionCatcher execution

### **Parallel Execution Managers**

**`scripts/run_parallel_FC.sh`**  
- Batch manager for FusionCatcher
- Memory monitoring and job throttling
- Automatic FASTQ pair discovery

**`scripts/run_cicero_parallel_docker.sh`**
- Containerized Cicero/RNApeg pipeline
- Singularity-based execution
- Automatic BAM indexing

### **Format Conversion Tools**

**`scripts/starfusion2cicero.py`**
- Converts STAR-Fusion outputs to Cicero format
- Batch processing with automatic file discovery

**`scripts/automated_convert_arriba2cicero.py`**
- Converts Arriba outputs to Cicero format  
- Handles aggregation of multiple fusion calls per gene pair

---

## 🛠 Troubleshooting

### **Common Issues**

**Memory Problems:**
- Reduce `MAX_PARALLEL_JOBS` 
- Increase `MIN_FREE_MEMORY_GB` buffer
- Monitor system resources with `htop`

**Reference Path Errors:**
- Verify all reference files exist and are readable
- Check file permissions
- Use absolute paths when possible

**Conda Environment Issues:**
```bash
# Recreate environments if corrupted
conda env remove -n fusion
conda env create -f env/fusion_environment.yml

# If facing issues with the YML files, use the other env files:
```

**Incomplete Analysis:**
- Check log files in the `logs/` directory
- Use `--force` flag to restart from beginning
- Verify input FASTQ file integrity

### **Performance Optimization**

**For High-Memory Systems (>256GB RAM):**
- Increase `MAX_PARALLEL_JOBS` to 8-12
- Set `MIN_FREE_MEMORY_GB` to 50-100

**For Many-Core Systems (>32 cores):**
- Increase `THREADS_PER_JOB` to 20-24
- Balance with `MAX_PARALLEL_JOBS` accordingly

---

## 📧 Contact

**Developed in the Fungel Lab by Rajarshi**

For questions, issues, or contributions:
- Open an issue on this GitHub repository
- Check the documentation for common solutions
- Review log files for detailed error messages

---

## 🚀 Quick Start Checklist

- [ ] Clone repository
- [ ] Create Conda environments  
- [ ] Download and configure reference data
- [ ] Edit script paths to match your system
- [ ] Test with a single sample first
- [ ] Scale up to batch processing
- [ ] Convert outputs to standardized format

---

**Happy Fusion Hunting! 🔬✨**
