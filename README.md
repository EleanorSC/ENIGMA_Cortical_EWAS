# ENIGMA_Cortical_EWAS
Brain imaging protocol for cortical EWAS and epigenetic age calculation 

This repository contains code for conducting a meta-analysis epigenome-wide association study (EWAS) investigating DNA methylation (DNAm) signatures associated with cortical morphology, focusing on two key measures:

Cortical thickness (CT)
Cortical surface area (SA)

The workflow is developed as part of the ENIGMA-Epigenetics Working Group, enabling large-scale, harmonized analyses across international cohorts.


## Repositry structure

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

The scripts are designed to:

Perform cohort-level quality control (QC) on DNAm and cortical measures, including aggregation of global and regional CT and SA metrics
Conduct cohort-level EWAS on cortical measures
Harmonise and quality-check EWAS results across multiple cohorts (QCEWAS)
Run fixed-effects meta-analyses using METAL on discovery cohorts, with replication testing in independent cohorts and combined datasets
Perform functional and genomic enrichment analyses, as well as phenome-wide association studies (PheWAS)
Conduct Mendelian randomisation analyses to investigate potential causal pathways from meQTLs → CT- or SA-associated DNAm → neuropsychiatric disorders or educational attainment

