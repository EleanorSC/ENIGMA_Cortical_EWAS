# LBC1936 DNAm QC and Cortical EWAS Pipeline

This script prepares LBC1936 Illumina 450K DNA methylation data for cortical EWAS analysis and runs CpG-level association models with cortical thickness and surface area phenotypes.

## Overview

The workflow:

1. Reads raw IDAT files and sample sheets using `minfi`
2. Creates and saves an `RGset`
3. Generates methylation QC reports
4. Performs Illumina and quantile normalisation
5. Removes low-quality, invariant, SNP-affected and sex chromosome probes
6. Generates methylation PCs and cell-count PCs
7. Merges DNAm, covariate and cortical imaging data
8. Runs EWAS models for cortical measures

## Main Inputs

```
Data/
├── idats/                         # Raw IDAT files and sample sheets
├── Input/
│   ├── Covariates.csv
│   └── CorticalMeasure_ENIGMA.csv
├── Imaging/
│   └── CorticalMeasuresENIGMA_ThickAvg.csv
└── Output/
```
##Main Outputs
```
Output/
├── RGset.rda
├── qcReport.pdf
├── Quan-norm.rda
├── Beta_Quantile.rda
├── fast_svd.rda
├── cellcount.rda
├── Methy.RData
├── PC_Beta.csv
├── pc_cell.csv
├── CorticalEWAS_covariates.csv
└── execution_log.txt
```

### R packages 
```
install.packages(c("tidyverse", "corpcor"))
install.packages("BiocManager")

BiocManager::install(c(
  "minfi",
  "minfiData",
  "FlowSorted.Blood.450k",
  "IlluminaHumanMethylation450kmanifest"
))

```

## Workflow
### 1. Read IDAT files

Raw methylation intensity data are loaded from IDAT files using cohort-specific sample sheets.
```
targets <- read.metharray.sheet(idat_dirs, sample_sheet, recursive = TRUE)

RGset <- read.metharray.exp(
  base = idat_dirs,
  targets = targets,
  verbose = TRUE
)
```
The resulting ```RGset``` is saved for reuse.

### 2. Generate QC reports

A standard ```minfi``` QC report is produced to identify problematic samples, batch structure and signal intensity issues
```
qcReport(
  RGset,
  sampNames = pd$Sentrix,
  sampGroups = pd$set,
  pdf = "Output/qcReport.pdf"
)
```
