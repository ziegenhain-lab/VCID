  # VCID

**V**ariant-based **C**ell-cell **I**nteraction **D**econvoluting(VCID)

## Abstract

This pipeline is called VCID, which stands for cell cell interaction deconvoluting. It is a pipeline for analyzing cell-cell interactions with single-cell transcriptome data.

Note: This pipeline is still under development 2025.6.21



## Usage Instructions for the VCID Pipeline

### 1.Prerequisites

Before using this pipeline, please ensure the following:

- Conda is installed on your system.

- The zUMIs software package is available (it will be invoked automatically if needed).

### 2.Overview of Pipeline Function

This pipeline is designed to deconvolute the transcriptome of doublets composed of two cells from distinct genetic backgrounds, based on SNP differences, and reconstruct the transcriptomes of individual cells.

**Use Cases**
The pipeline can be applied in several different input scenarios:

- Only RNA-seq FASTQ files are provided:
  - The pipeline will run zUMIs to generate BAM files.
  - SNP calling will be performed on the BAM files to produce VCF files.
  - Read separation will be performed based on the generated VCF.

- RNA-seq FASTQ files + a provided VCF file:
  - The pipeline will run zUMIs to generate BAM files.
  - Read separation will be performed using the provided VCF.

- BAM file + a provided VCF file:
  - The pipeline will skip zUMIs and SNP calling steps.
  - It will directly perform read separation using the provided files.

### 3.Workflow

#### 1) Edit the configuration file

Update config/sum_config.yaml with paths and parameters specific to your data.

#### 2) Create and activate a conda environment with Snakemake

```
conda create -n VCID_snakemake -c conda-forge -c bioconda snakemake
conda activate VCID_snakemake
```

#### 3) Run the pipeline within the VCID_snakemake environment

- All non-zUMIs rules are configured to use conda.

- Snakemake will automatically resolve and build the required environments.

#### 4) Navigate to the pipeline directory and run Snakemake

```
cd VCID4.0/algorithem/snakemake_pipeline
snakemake CCID_all --use-conda --cores 20 --jobs 1
```

⚠️ Note: Please set --jobs to 1. Running multiple zUMIs jobs in parallel may lead to conflicts.   
⚠️ Note: Do not set the number of cores higher than the num_threads specified in the YAML file.

***

### 4.Code and Test Data Used During Development
This repository contains scripts to generate test datasets during the development and testing.

#### 1) FASTQ Downsampling Using BBTools
Downsampled FASTQ files are stored in: `/mnt/data/user/zhouhui/downsampled_fastq`

These files were generated using the BBTools suite.
The corresponding script can be found at: `algorithem/other/bbtools_fastq_downsampling.sh`

These downsampled FASTQ files can be used for full-pipeline testing, including zUMIs-based workflows.

#### 2) Extracted BAM for Direct Testing
Extracted BAM files for controlled testing are located at: 
`/mnt/data/user/zhouhui/T_cell_dynamics_20kcells_results_extracted_bam`

These files were generated using:
`algorithem/other/extract_random_bam_by_condition_for_testing.sh`

These BAM files are specifically prepared for testing read separation pipelines that rely on external BAM and VCF files.

### Tips and tricks in use

When you modify the YAML file or some piece of code but believe that the change does not affect previous results — and you want to prevent Snakemake from automatically rerunning all previously executed rules — you can use the following command to "fake" the output file of a specific rule as up-to-date, thereby avoiding the issue:

```
snakemake -R your_rule_name --touch
```

## Flow chart of VCID
![alt text](algorithem/workflow_demo/mermaid-diagram-2025-06-18-130830.png)