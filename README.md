# ENIGMA_Cortical_EWAS

Pipeline for conducting epigenome-wide association studies (EWAS) of cortical brain structure. This repository provides a harmonised workflow for analysing associations between DNA methylation (DNAm) and cortical morphology derived from MRI.

---

## Overview

```
┌──────────────────────────────────────────────────────────────┐
│                    DNA METHYLATION (DNAm)                    │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
            IDAT import 
                │
                ▼
      ┌──────────────────────────────┐
      │        SAMPLE-LEVEL QC       │
      │  • qcReport (intensity)      │
      │  • sex check (predicted vs)  │
      │  • MDS / PCA outliers        │
      │  • remove controls/outliers  │
      └──────────────────────────────┘
                │
                ▼
      ┌──────────────────────────────┐
      │        PROBE-LEVEL QC        │
      │  • detection P filtering     │
      │  • SNP removal (CpG/SBE)     │
      │  • invariant probes          │
      │  • remove chrX / chrY        │
      └──────────────────────────────┘
                │
                ▼
      Quantile normalisation (preprocessQuantile)
                │
                ▼
      ┌──────────────────────────────┐
      │    FEATURE EXTRACTION        │
      │  • Beta matrix (CpGs × N)    │
      │  • DNAm PCs (SVD)            │
      │  • Cell-type proportions     │
      └──────────────────────────────┘
                │
                ▼
      OUTPUT: QCed DNAm + covariates



┌──────────────────────────────────────────────────────────────┐
│                        MRI (STRUCTURAL)                      │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
      FreeSurfer preprocessing
                │
                ▼
      ┌──────────────────────────────┐
      │      ROI EXTRACTION          │
      │  • Cortical thickness (CT)   │
      │  • Surface area (SA)         │
      └──────────────────────────────┘
                │
                ▼
      ┌──────────────────────────────┐
      │      MRI QC & CLEANING       │
      │  • visual QC                 │
      │  • remove segmentation fails │
      └──────────────────────────────┘
                │
                ▼
      ┌──────────────────────────────┐
      │      HARMONISATION           │
      │  • site/scanner (ComBat)     │
      │  • adjust: age, sex, ICV     │
      └──────────────────────────────┘
                │
                ▼
      OUTPUT: Harmonised CT / SA ROIs



┌──────────────────────────────────────────────────────────────┐
│                    DATA INTEGRATION                          │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
      ┌──────────────────────────────┐
      │    ID MATCHING / MERGING     │
      │  • DNAm ↔ MRI ↔ covariates   │
      │  • retain complete cases     │
      └──────────────────────────────┘
                │
                ▼
      FINAL ANALYSIS MATRICES:
      • DNAm (CpGs)
      • Cortical phenotypes (CT/SA)
      • Covariates (age, sex, ICV, PCs, cell counts)



┌──────────────────────────────────────────────────────────────┐
│                         EWAS                                 │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
      Model: DNAm ~ Cortical ROI + covariates
                │
                ▼
      Loop over:
      • CpGs (~500k–800k)
      • ROIs (global + regional)
                │
                ▼
      OUTPUT (per cohort):
      • β (effect sizes)
      • SE
      • P-values
      • N



┌──────────────────────────────────────────────────────────────┐
│                    META-ANALYSIS (ENIGMA)                    │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
      ┌──────────────────────────────┐
      │        QCEWAS                │
      │  • inflation correction      │
      │  • harmonise CpGs/models     │
      └──────────────────────────────┘
                │
                ▼
      Fixed-effects meta-analysis (METAL)
                │
                ▼
      • Discovery cohorts
      • Independent replication
      • Combined analysis



┌──────────────────────────────────────────────────────────────┐
│                    DOWNSTREAM ANALYSES                       │
└──────────────────────────────────────────────────────────────┘
                │
                ▼
      • Functional annotation (GREAT, eFORGE)
      • Genomic enrichment (pathways, chromatin)
      • PheWAS (cross-trait associations)
      • Mendelian randomisation (meQTL → DNAm → brain)
      • Sensitivity analyses (cell counts, lifestyle)



───────────────────────────────────────────────────────────────
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
### Example cohort DNAm QC
#### LBC1936 DNAm QC and Cortical EWAS 

Illumina 450K DNA methylation data for cortical EWAS analysis and runs CpG-level association models with cortical thickness and surface area phenotypes.

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
## 3. Normalise methylation data

The pipeline applies quantile normalisation using `minfi::preprocessQuantile`. This step performs stratified quantile normalisation and can automatically remove low-quality samples using the `removeBadSamples` argument.

```r
object <- preprocessQuantile(
  RGset,
  fixOutliers = TRUE,
  removeBadSamples = TRUE,
  badSampleCutoff = 10.5,
  quantileNormalize = TRUE,
  stratified = TRUE,
  mergeManifest = FALSE,
  sex = NULL,
  verbose = TRUE
)

save(object, file = "EleanorC/CorticalEWAS/Data/Output/Quan-norm.rda")
```

The normalised beta values are then extracted and sex chromosome probes are removed.

```r
beta <- getBeta(object)

keepIndex <- which(!seqnames(object) %in% c("chrX", "chrY"))
beta <- beta[keepIndex, ]

save(beta, file = "EleanorC/CorticalEWAS/Data/Output/Beta_Quantile.rda")
```

---

## 4. Generate methylation principal components

Principal components are generated from the normalised methylation beta matrix to capture unknown technical and biological structure.

```r
library(corpcor)

ss <- fast.svd(beta, tol = 0)

save(ss, file = "EleanorC/CorticalEWAS/Data/Output/fast_svd.rda")
```

The first four methylation PCs are used as covariates in downstream EWAS models.

```r
PC_Beta <- as.data.frame(ss$v[, 1:4])
PC_Beta$Subject <- Methy$Subject

write.csv(
  PC_Beta,
  "EleanorC/CorticalEWAS/Data/Output/PC_Beta.csv"
)
```

---

## 5. Estimate or load cell-type proportions

Estimated blood cell proportions are used to control for cellular heterogeneity in peripheral blood DNAm data.

```r
load("DNAm_3041_ppts_with_wave_and_cohort_13May2016.RData")

rownames(data) <- data$SampleID.1

cellcount <- data[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]

save(
  cellcount,
  file = "EleanorC/CorticalEWAS/Data/Output/cellcount.rda"
)
```

The first two cell-count PCs are generated and saved.

```r
tmp <- prcomp(cellcount)

pc_cell <- as.data.frame(tmp$x[, 1:2])
pc_cell$Subject <- rownames(pc_cell)

write.csv(
  pc_cell,
  "EleanorC/CorticalEWAS/Data/Output/pc_cell.csv"
)
```

---

## 6. Post-normalisation QC plots

MDS plots are generated after quantile normalisation to assess remaining batch effects.

```r
Cohort <- "LBC1936"

mydist <- dist(t(beta[sample(1:nrow(beta), 10000), ]))
mds <- cmdscale(mydist)

pd <- pData(dat)

pdf(
  paste(
    "EleanorC/CorticalEWAS/Data/Output/Batch_effect_Quantile_Normalised_mds_",
    Cohort,
    ".pdf",
    sep = ""
  ),
  width = 11,
  height = 4
)

plot(
  mds[, 1],
  col = as.numeric(as.factor(pd[, "array"])) + 1,
  xlab = "",
  ylab = "First Principal component",
  xaxt = "n"
)

dev.off()
```

---

## 7. Probe-level QC

If methylation data have not already been QCed using the ENIGMA-Epigenetics pipeline, probe-level filtering is performed.

```r
library(minfi)

load("EleanorC/CorticalEWAS/Data/Output/RGset.rda")
load("EleanorC/CorticalEWAS/Data/Output/Quan-norm.rda")
```

SNP information is added and probes overlapping SNPs are removed.

```r
objectWithSNPinfo <- addSnpInfo(object)

objectSNPQCed <- dropLociWithSnps(
  objectWithSNPinfo,
  snps = c("SBE", "CpG", "Probe"),
  maf = 0.05
)

rm(objectWithSNPinfo)
```

Detection p-values are calculated and probes failing detection thresholds are marked.

```r
detP <- detectionP(RGset)

Match1 <- match(colnames(objectSNPQCed), colnames(detP))
Match2 <- match(rownames(objectSNPQCed), rownames(detP))

detPSNPQCed <- detP[
  Match2[!is.na(Match2)],
  Match1[!is.na(Match1)]
]

failed <- detPSNPQCed > 0.01
```

Probes failing in more than 20% of samples are removed.

```r
beta <- getBeta(objectSNPQCed)

failedCG02 <- rowMeans(failed) > 0.2
```

Invariant probes are identified and saved.

```r
ProbeInvar <- (rowSums(beta <= 0.2) == ncol(beta)) |
              (rowSums(beta >= 0.8) == ncol(beta))

ListInvarProbe <- rownames(beta)[which(ProbeInvar)]
```

Sex chromosome probes and failed probes are removed.

```r
keepIndex <- !seqnames(objectSNPQCed) %in% c("chrX", "chrY")
keepIndex <- keepIndex & (!failedCG02)

beta[failed] <- NA

betaQC <- beta[which(keepIndex), ]
```

---

## 8. Reformat methylation data

The QCed beta matrix is transposed so that rows correspond to samples and columns correspond to CpG probes.

```r
Methy <- as.data.frame(t(betaQC))

Methy$Subject <- colnames(betaQC)

Methy <- Methy[, c(ncol(Methy), 1:(ncol(Methy) - 1))]

MethyName <- colnames(Methy)[-1]

save(
  Methy,
  file = "EleanorC/CorticalEWAS/Data/Output/Methy.RData"
)
```

---

## 9. Prepare cortical EWAS covariates

Imaging-derived intracranial volume and demographic covariates are merged with methylation PCs and cell-count PCs.

```r
covariates_examine2 <- read.csv(
  "EleanorC/CorticalEWAS/Data/Input/Covariates.csv"
)

thickness_examine <- read.csv(
  "EleanorC/CorticalEWAS/Data/Imaging/CorticalMeasuresENIGMA_ThickAvg.csv"
)

ICV <- thickness_examine %>% select("SubjID", "ICV")

RawCov <- merge(ICV, covariates_examine2, by = "SubjID")

names(RawCov)[names(RawCov) == "sex"] <- "Sex"
names(RawCov)[names(RawCov) == "ageyrs"] <- "Age"
names(RawCov)[names(RawCov) == "SubjID"] <- "Subject"

RawCov$Age_Square <- RawCov$Age^2
```

Subject IDs are matched across methylation and imaging datasets.

```r
ID_linkage <- data.frame(
  Subject_meth = dat$Sentrix,
  Subject = dat$ID
)

RawCov <- merge(RawCov, ID_linkage, by = "Subject")
```

Methylation and cell-count PCs are merged into the covariate file.

```r
names(PC_Beta)[names(PC_Beta) == "Subject"] <- "Subject_meth"
names(pc_cell)[names(pc_cell) == "Subject"] <- "Subject_meth"

Cov <- merge(RawCov, PC_Beta, by = "Subject_meth", all = FALSE)
Cov <- merge(Cov, pc_cell, by = "Subject_meth", all = FALSE)

Cov$Sex <- as.numeric(Cov$Sex)

write.csv(
  Cov,
  "EleanorC/CorticalEWAS/Data/Output/CorticalEWAS_covariates.csv"
)
```

---

## 10. Merge methylation, covariates and cortical phenotypes

The final analysis dataset is generated by matching methylation, covariate and cortical phenotype files.

```r
Cov <- read.csv(
  "EleanorC/CorticalEWAS/Data/Output/CorticalEWAS_covariates.csv"
)

Cov <- Cov %>% select(-c("X"))

Cov$ICV <- as.numeric(Cov$ICV)
Cov$Sex <- as.numeric(Cov$Sex)

CorticalMeasure <- read.csv(
  "EleanorC/CorticalEWAS/Data/Input/CorticalMeasure_ENIGMA.csv",
  header = TRUE
)

Data <- merge(Cov, CorticalMeasure, by = "Subject", all = FALSE)

Match <- match(Methy$Subject, Data$Subject_meth)
```

Analysis matrices are then created.

```r
Cov <- Cov %>% select(-c("Subject", "Subject_meth"))

CorticalMeasure <- Data %>%
  select(-c(
    "Subject",
    "Subject_meth",
    "ICV",
    "Sex",
    "Age",
    "Age_Square",
    "V1",
    "V2",
    "V3",
    "V4",
    "PC1",
    "PC2"
  ))

Methy <- Methy[!is.na(Match), -1]
```

Age summaries are saved for downstream reporting.

```r
AgeInfo <- data.frame(
  min(Cov$Age),
  max(Cov$Age),
  mean(Cov$Age)
)

names(AgeInfo) <- c("MinAge", "MaxAge", "MeanAge")
```

---

## 11. Run cortical EWAS

The EWAS tests each CpG against each cortical phenotype using linear regression.

Model:

```r
DNAm ~ covariates + cortical_measure
```

The output matrices are initialised as follows:

```r
Num_Methy <- ncol(Methy)
Num_Cov <- ncol(Cov)
Num_Cortical <- ncol(CorticalMeasure)

Origin_Beta <- matrix(NA, nrow = Num_Methy, ncol = Num_Cortical)
Origin_SE   <- matrix(NA, nrow = Num_Methy, ncol = Num_Cortical)
Origin_P    <- matrix(NA, nrow = Num_Methy, ncol = Num_Cortical)
Origin_N    <- matrix(NA, nrow = Num_Methy, ncol = Num_Cortical)

colnames(Origin_Beta) <- colnames(CorticalMeasure)
colnames(Origin_SE)   <- colnames(CorticalMeasure)
colnames(Origin_P)    <- colnames(CorticalMeasure)
colnames(Origin_N)    <- colnames(CorticalMeasure)

rownames(Origin_Beta) <- colnames(Methy)
rownames(Origin_SE)   <- colnames(Methy)
rownames(Origin_P)    <- colnames(Methy)
rownames(Origin_N)    <- colnames(Methy)
```

The CpG-by-phenotype regression loop is then run.

```r
start_time <- Sys.time()

Covar <- matrix(unlist(Cov), ncol = ncol(Cov), byrow = FALSE)

for (i in 1:Num_Methy) {
  for (j in 1:Num_Cortical) {
    Out <- summary(lm(Methy[, i] ~ Covar + CorticalMeasure[, j]))

    Origin_Beta[i, j] <- Out$coefficients[nrow(Out$coefficients), 1]
    Origin_SE[i, j]   <- Out$coefficients[nrow(Out$coefficients), 2]
    Origin_P[i, j]    <- Out$coefficients[nrow(Out$coefficients), 4]
    Origin_N[i, j]    <- Out$df[1] + Out$df[2]
  }
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time
```

Execution time is saved to a log file.

```r
log_file <- "execution_log.txt"

cat(
  "Linear regression execution time:",
  elapsed_time,
  "seconds",
  file = log_file
)
```

---

## Main Outputs

```text
Output/
├── Quan-norm.rda
├── Beta_Quantile.rda
├── fast_svd.rda
├── cellcount.rda
├── Methy.RData
├── PC_Beta.csv
├── pc_cell.csv
├── CorticalEWAS_covariates.csv
├── Batch_effect_Quantile_Normalised_mds_LBC1936.pdf
└── execution_log.txt
```

EWAS result objects:

```text
Origin_Beta   # Effect estimates
Origin_SE     # Standard errors
Origin_P      # P-values
Origin_N      # Sample sizes
```

---


## Final Outputs

Ultimately, per cohort, each EWAS run produces .RData files containing:
- ```Origin_Beta```
- ```Origin_SE```
- ```Origin_P```
- ```Origin_N```
- ```Model specification```
- ```Age summary statistics```
- ```Invariant probe list```

---

---
