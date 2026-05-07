# ENIGMA_Cortical_EWAS

Pipeline for conducting epigenome-wide association studies (EWAS) of cortical brain structure. This repository provides a harmonised workflow for analysing associations between DNA methylation (DNAm) and cortical morphology derived from MRI.

---

## Repository Structure
```text
ENIGMA_Cortical_EWAS/
│
├── Code/
│   └── Cortical_ROI/
│       └── Create ENIGMA2023_CortexROI_plots.R
│
├── DNAm_QC/
│   └── qcReport.pdf
│
├── Prep/
│   └── LBC1936_ENIGMA_DNAm_QC.R
│
├── Scripts/
│   ├── Association_analysis.R
│   ├── Debuggin_cortical_measures.R
│   ├── Debugging_script.R
│   ├── QC_and_reformat_methylation_data.R
│   ├── Remove_inflation_and_bias.R
│   ├── get_CT_and_SA_from_FreeSurfer.sh
│   ├── perform_EWAS.R
│   └── prep_cortical_measures.R
│
├── LICENSE
└── README.md
```
---

## Overview
```
DNAm → QC → PCs
           ↓
MRI → ROI extraction → harmonisation
           ↓
        EWAS
           ↓
   Meta-analysis / downstream
```

### 1. DNAm Quality Control

Script:  
```Scripts/QC_and_reformat_methylation_data.R```

Steps:
- Remove probes with SNPs (CpG, SBE, probe-level)
- Filter probes using detection p-values (p > 0.01)
- Remove probes failing in >20% of samples
- Remove invariant probes (β ≤ 0.2 or ≥ 0.8 across all samples)
- Remove sex chromosome probes (chrX, chrY)
- Generate cleaned beta matrix (Methy)
- Compute DNAm PCs and cell-type PCs

Outputs:
- QCed methylation matrix
- Invariant probe list
- Principal components

---

### 2. Cortical Measure Extraction

Script:  
```Scripts/get_CT_and_SA_from_FreeSurfer.sh```

- Extract cortical thickness and surface area from FreeSurfer outputs
- Generate ROI-level cortical measures

---

### 3. Cortical Data Preparation

Script:  
```Scripts/prep_cortical_measures.R```

- Format cortical measures
- Harmonise ROI naming (ENIGMA standards)
- Merge with subject IDs

---

### 4. Covariate Preparation

- Age and Age²  
- Sex  
- Intracranial volume (ICV)  
- DNAm PCs (1–4)  
- Cell-type PCs (1–2)  
- Batch / technical variables  

---

### 5. EWAS Analysis

Primary script:  
```Scripts/perform_EWAS.R```

Core model:
DNAm ~ covariates + cortical phenotype

Analyses performed:
- Full sample (with ICV)
- Full sample (without ICV)
- Sex-stratified analyses

Outputs:
- Effect sizes (β)
- Standard errors (SE)
- p-values
- Sample size (N)

---

### 6. Association & QC Checks

Scripts:
- ```Scripts/Association_analysis.R```
- ```Scripts/Remove_inflation_and_bias.R```

Includes:
- Post-EWAS quality control
- Inflation correction
- Model diagnostics

---

### 7. Debugging Utilities

Scripts:
- ```Scripts/Debugging_script.R```
- ```Scripts/Debuggin_cortical_measures.R```

Used for:
- Data alignment issues
- Missingness checks
- Pipeline validation

---

### 8. Visualisation

Script:
```Code/Cortical_ROI/Create ENIGMA2023_CortexROI_plots.R```

- Generates cortical maps
- Visualises EWAS results across ROIs

---

## Input Requirements
Each cohort must provide:

### DNAm Data
- Beta values (normalised)
- Sample IDs matching phenotype data

### Neuroimaging Data
- FreeSurfer-derived cortical measures (CT, SA)

### Covariates
- Age, Sex
- ICV 
- Technical covariates
- DNAm and cell PCs

---

## Outputs

Each EWAS run produces .RData files containing:
- ```Origin_Beta```
- ```Origin_SE```
- ```Origin_P```
- ```Origin_N```
- ```Model specification```
- ```Age summary statistics```
- ```Invariant probe list```

---

## Analytical Notes

- DNAm is treated as the dependent variable (standard EWAS convention)
- Linear regression is performed at each CpG × phenotype pair
- Multiple testing correction is applied downstream (meta-analysis stage)
- Stratified analyses are exploratory

---

## Dependencies

- R (≥ 4.0)
- Key packages:
  - minfi
  - stats
  - data.table

---

