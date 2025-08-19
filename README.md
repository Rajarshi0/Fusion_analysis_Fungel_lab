# 🧬 High-Throughput Gene Fusion Analysis Pipeline – Fungel Lab

This repository contains a complete suite of scripts for performing high-throughput gene fusion analysis on RNA-seq data. It includes parallel execution wrappers for multiple fusion callers (**STAR-Fusion**, **Arriba**, **FusionCatcher**, and **Cicero/RNApeg**), Conda environment files for reproducible setup, and Python scripts for post-processing results into a standardized format.

---
## 📌 Table of Contents
- [Features](#-features)
- [System Requirements](#-system-requirements)
- [Setup & Installation](#-setup--installation)
- [Analysis Workflow](#-analysis-workflow)
- [Usage Instructions](#-usage-instructions)
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
## 🚀 Usage Instructions

**Note**: All `run_parallel*.sh` scripts are designed to be run from the directory containing your FASTQ files.

### 1️⃣ Running the STAR-Fusion / Arriba Pipeline
1.  Navigate to your directory of paired-end FASTQ files.
2.  Edit the paths inside `STAR_f_arriba_combined.sh`.
3.  Edit the `SCRIPT` variable inside `run_parallel.sh` to point to `STAR_f_arriba_combined.sh`.
4.  Execute the parallel wrapper:
    ```bash
    bash run_parallel.sh
    ```

### 2️⃣ Running the FusionCatcher Pipeline
1.  Navigate to your directory of paired-end FASTQ files.
2.  Edit the paths inside `fusioncatcher_parallel_2.sh`.
3.  Edit the `SCRIPT` variable inside `run_parallel_FC.sh` to point to `fusioncatcher_parallel_2.sh`.
4.  Execute the parallel wrapper:
    ```bash
    bash run_parallel_FC.sh
    ```

### 3️⃣ Running the Cicero / RNApeg Pipeline
1.  This pipeline uses BAM files generated from a previous alignment step (like the one in `STAR_f_arriba_combined.sh`).
2.  Edit all paths inside `run_cicero_parallel_docker.sh`, especially the `HOST_BASE_PATH` and reference locations.
3.  Execute the script:
    ```bash
    bash run_cicero_parallel_docker.sh
    ```

### 4️⃣ Converting Outputs to Cicero Format
1.  Activate the `fusion` Conda environment.
2.  Edit the `SEARCH_DIRECTORY` and `OUTPUT_DIRECTORY` paths in the Python scripts.
3.  Run the desired conversion script:
    ```bash
    # For STAR-Fusion results
    python starfusion2cicero.py

    # For Arriba results
    python automated_convert_arriba2cicero.py
    ```

---
## 📂 Script Descriptions

-   **`run_parallel.sh` / `run_parallel_FC.sh`**: Master scripts for launching batch analyses. They manage parallel job slots and monitor system memory to prevent overloading.
-   **`STAR_f_arriba_combined.sh`**: Per-sample script that runs `fastp` trimming, STAR alignment, and then STAR-Fusion and Arriba on the alignment outputs.
-   **`fusioncatcher_parallel_2.sh`**: Per-sample script that runs `fastp` trimming and then FusionCatcher.
-   **`run_cicero_parallel_docker.sh`**: Per-sample script that orchestrates the RNApeg and Cicero containers to perform fusion detection from existing BAM files.
-   **`starfusion2cicero.py` / `automated_convert_arriba2cicero.py`**: Post-processing scripts to standardize the outputs from STAR-Fusion and Arriba into the uniform Cicero format.
-   **`*.yml` / `*.txt`**: Conda environment files for ensuring a reproducible software stack.

---
## ✨ Planned Improvements
-   Consolidate parallel runners into a single, more flexible script.
-   Develop a Nextflow or Snakemake version of the pipeline for improved portability and workflow management.
-   Incorporate additional QC and fusion filtering steps post-conversion.

---
## 📧 Contact
Developed in the Fungel Lab by Rajarshi.

For questions or issues, please open an issue on this GitHub repository.
