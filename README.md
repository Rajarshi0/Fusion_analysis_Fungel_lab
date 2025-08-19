# 🧬 High-Throughput Gene Fusion Analysis Pipeline – Fungel Lab

Welcome! This repository contains a complete suite of scripts for performing high-throughput gene fusion analysis on RNA-seq data. It includes parallel execution wrappers for multiple fusion callers (**STAR-Fusion**, **Arriba**, **FusionCatcher**, and **Cicero/RNApeg**), Conda environment files for reproducible setup, and Python scripts for post-processing results into a standardized format.

---
## 📌 Table of Contents
- [Features](#-features)
- [System Requirements](#-system-requirements)
- [Setup & Installation](#-setup--installation)
- [A Beginner's Guide to Parallel Execution & Resources](#-a-beginners-guide-to-parallel-execution--resources)
- [Analysis Workflow](#-analysis-workflow)
- [How to Run an Analysis (Current Usage)](#-how-to-run-an-analysis-current-usage)
- [Script Descriptions](#-script-descriptions)
- [Planned Improvements](#-planned-improvements)
- [Contact](#-contact)

---
## ✨ Features
- **Multiple Fusion Calling Pipelines**: End-to-end scripts for:
  - STAR-Fusion + Arriba
  - FusionCatcher
  - Cicero + RNApeg
- **Parallel Execution**: Sophisticated wrapper scripts (`run_parallel.sh`, etc.) manage concurrent jobs, system threads, and available memory to process large sample batches efficiently.
- **Reproducible Environments**: Conda `.yml` and explicit package lists (`fusion_environment.yml`, `fusioncatcher_environment.yml`) ensure a consistent and correct software setup.
- **Automated Post-Processing**: Python scripts (`starfusion2cicero.py`, `automated_convert_arriba2cicero.py`) automatically find and convert raw fusion calls into the standardized Cicero format.
- **Containerized Workflow**: Integration with **Docker** and **Singularity/Apptainer** for the Cicero/RNApeg pipeline, ensuring portability and reproducibility.

---
## 🛠 System Requirements
- **Operating System**: Linux
- **Software**:
    - [Conda / Miniconda](https://docs.conda.io/en/latest/miniconda.html) for environment management.
    - [Docker](https://www.docker.com/) (for the RNApeg container).
    - [Singularity / Apptainer](https://apptainer.org/) (for the Cicero container).
- **Reference Data**: Pre-built reference genomes, annotations, and databases are required for each fusion caller (e.g., STAR indexes, CTAT library, FusionCatcher DB).

---
## ⚙️ Setup & Installation

1.  **Clone the Repository**
    ```bash
    git clone [https://github.com/Rajarshi0/Fusion_analysis_Fungel_lab.git](https://github.com/Rajarshi0/Fusion_analysis_Fungel_lab.git)
    cd Fusion_analysis_Fungel_lab
    ```

2.  **Create Conda Environments**
    Two distinct environments are needed. Create them using the provided files:

    * **For STAR-Fusion, Arriba, and Python 3 post-processing:**
        ```bash
        conda env create -f fusion_environment.yml
        ```

    * **For FusionCatcher (which requires Python 2.7):**
        ```bash
        conda env create -f fusioncatcher_environment.yml
        ```

3.  **Prepare Reference & Container Files**
    - Download or build all necessary reference genomes, aligner indexes, and fusion databases.
    - Pull or build the required Docker/Singularity images for the Cicero/RNApeg pipeline.
    - **Crucially, you must edit the configuration paths** at the top of all `.sh` scripts to point to the correct locations of these resources on your system.

---
## 🧠 A Beginner's Guide to Parallel Execution & Resources

To use this pipeline effectively, you need to tell the scripts how to use your server's power. Let's break it down.

### What are Cores and Threads? (An Analogy)
Imagine your computer's brain (the CPU) is a large kitchen.
-   **Cores** are like individual chefs in the kitchen. A kitchen with 8 cores has 8 chefs who can all work on different recipes at the same time.
-   **Threads** are like the number of hands each chef has. If each chef has 2 hands (threads), they can chop vegetables and stir a pot simultaneously.

So, a CPU with 8 cores and 2 threads per core has a total of **16 threads** (8 chefs * 2 hands each). This is the total processing power you can use.

### How to Check Your Server's Power
To see how many "hands" your server's kitchen has, open a terminal and run one of these simple commands:
```bash
# Shows the total number of threads. It's the quickest way!
nproc
```
or for more detail:
```bash
# Shows a detailed breakdown, including cores and threads per core.
lscpu
```

### How to Configure the Scripts
Our parallel runner scripts have two important settings you need to adjust:
-   `THREADS_PER_JOB`: How many "hands" (threads) to give to **one sample's analysis**. A good value for alignment is usually between 8 and 16.
-   `MAX_PARALLEL_JOBS`: How many **different samples (recipes) to work on at the same time**.

**The Golden Rule:**
To keep your server happy and running smoothly, the total number of threads you use should not be more than what the server has.

`THREADS_PER_JOB` × `MAX_PARALLEL_JOBS` ≤ `Total Threads from nproc`

**Example:**
Your server's `nproc` command returns `64`.
-   **Safe Choice:** You could set `THREADS_PER_JOB=16` and `MAX_PARALLEL_JOBS=4`. (16 * 4 = 64). This gives a lot of power to 4 samples at once.
-   **Alternative Choice:** You could set `THREADS_PER_JOB=8` and `MAX_PARALLEL_JOBS=8`. (8 * 8 = 64). This processes more samples at once, but with slightly less power for each.

Finally, the scripts also have a `MIN_FREE_MEMORY_GB` setting. This is a safety net. If the server's memory (RAM) gets too low, the script will wait before starting new jobs. Set this to about 10-20% of your total RAM to be safe.

---
##  workflow Analysis Workflow
This repository supports three main analysis workflows that can be run independently.

#### Workflow 1: STAR-Fusion & Arriba
`Raw FASTQs` → `run_parallel.sh` → `STAR_f_arriba_combined.sh` (Trimming & Alignment) → `STAR-Fusion & Arriba Calls`

#### Workflow 2: FusionCatcher
`Raw FASTQs` → `run_parallel_FC.sh` → `fusioncatcher_parallel_2.sh` (Trimming & Calling) → `FusionCatcher Calls`

#### Workflow 3: Cicero & RNApeg
`Aligned BAMs` → `run_cicero_parallel_docker.sh` (Junction Processing & Calling) → `Cicero Calls`

#### Post-Processing (Optional)
`Raw Fusion Calls` → `starfusion2cicero.py` / `automated_convert_arriba2cicero.py` → `Cicero-Formatted TSV`

---
## 🚀 How to Run an Analysis (Current Usage)

**Important**: Before running, you must **manually edit the variables inside the scripts** to set paths and resource limits.

### 1️⃣ Running the STAR-Fusion / Arriba Pipeline
1.  **Edit `STAR_f_arriba_combined.sh`**: Open this file and set the correct paths for `CTAT_LIB_DIR`, `STAR_INDEX_DIR`, etc.
2.  **Edit `run_parallel.sh`**:
    -   Set `SCRIPT="STAR_f_arriba_combined.sh"`.
    -   Configure `THREADS_PER_JOB` and `MAX_PARALLEL_JOBS` according to your server's resources.
3.  **Execute**: Place `run_parallel.sh` in the directory containing your FASTQ files and run it:
    ```bash
    bash run_parallel.sh
    ```

### 2️⃣ Running the FusionCatcher Pipeline
1.  **Edit `fusioncatcher_parallel_2.sh`**: Set the correct path for `FUSIONCATCHER_DB`.
2.  **Edit `run_parallel_FC.sh`**:
    -   Set `SCRIPT="fusioncatcher_parallel_2.sh"`.
    -   Configure `THREADS_PER_JOB` and `MAX_PARALLEL_JOBS`.
3.  **Execute**: Place `run_parallel_FC.sh` in the directory with your FASTQ files and run it:
    ```bash
    bash run_parallel_FC.sh
    ```

### 3️⃣ Running the Cicero / RNApeg Pipeline
1.  **Edit `run_cicero_parallel_docker.sh`**: This script requires careful path setup. Update `HOST_BASE_PATH`, `PROJECT_DIR`, and all reference file paths.
2.  **Execute**: This script finds BAM files from previous alignment runs. Run it from your main analysis directory:
    ```bash
    bash run_cicero_parallel_docker.sh
    ```

### 4️⃣ Converting Outputs to Cicero Format
1.  **Activate Environment**: `conda activate fusion`
2.  **Edit Python Script**: Open either `starfusion2cicero.py` or `automated_convert_arriba2cicero.py` and set the `SEARCH_DIRECTORY` and `OUTPUT_DIRECTORY`.
3.  **Execute**:
    ```bash
    # For STAR-Fusion results
    python starfusion2cicero.py
    ```

---
## 📂 Script Descriptions

### Master Scripts (Orchestrators)
-   **`run_parallel.sh` / `run_parallel_FC.sh`**: These are the main entry points for batch processing. They act as job managers by:
    1.  Automatically discovering paired-end FASTQ samples in a directory.
    2.  Launching a dedicated analysis script for each sample.
    3.  Enforcing limits on the number of concurrent jobs (`MAX_PARALLEL_JOBS`).
    4.  Monitoring system memory and pausing new jobs if it runs low.
    5.  Logging the output of each job separately.

-   **`run_cicero_parallel_docker.sh`**: The master script for the containerized workflow. It finds input BAM files and orchestrates the execution of RNApeg and Cicero using Singularity.

### Core Analysis Scripts (Per-Sample)
-   **`STAR_f_arriba_combined.sh`**: A comprehensive script that processes one sample. It performs `fastp` trimming, STAR alignment, and then runs both STAR-Fusion and Arriba on the resulting BAM file.
-   **`fusioncatcher_parallel_2.sh`**: A focused script that processes one sample by running `fastp` trimming followed by the FusionCatcher tool.

### Post-Processing Scripts
-   **`starfusion2cicero.py` / `automated_convert_arriba2cicero.py`**: Python scripts that recursively search a directory for raw fusion caller outputs and convert them into a uniform, Cicero-compliant TSV format.
-   **`*.yml` / `*.txt`**: Conda environment files for ensuring a reproducible software stack across different systems.

---
## ✨ Planned Improvements
-   **Parameterize Scripts**: Modify all shell scripts to accept command-line arguments for paths, threads, and other settings, eliminating the need for manual editing.
-   **Develop a Workflow Manager Version**: Create a Nextflow or Snakemake version of the pipeline for improved portability, scalability, and workflow management.
-   **Add QC and Filtering**: Incorporate additional steps for quality control and filtering of fusion calls post-conversion.

---
## 📧 Contact
Developed in the Fungel Lab by Rajarshi.

For questions or issues, please open an issue on this GitHub repository.
