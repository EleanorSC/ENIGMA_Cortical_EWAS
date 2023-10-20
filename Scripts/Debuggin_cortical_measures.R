# /CCACE_Shared/EleanorC/CorticalEWAS/Data/Input

setwd("/CCACE_Shared/EleanorC/CorticalEWAS/Data")

#BETAS <- readRDS("LBC_betas_3489_bloodonly.rds")


install.packages("tidyverse")
library(tidyverse)

cohort = "LBC1936"

## ----------------------------# 

Imaging_examine <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Imaging/ENIGMA_CNV_ImagingFile.csv")
thickness_examine <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Imaging/CorticalMeasuresENIGMA_ThickAvg.csv") # contains ICV
covariates_examine <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Imaging/LBC36_Covar_ENIGMACNV.csv")


ICV <- thickness_examine %>% select("SubjID", "ICV")

## ----------------------------#

# It looks like the .aparc files are just pre-versions of my files, separated by hemisphere?
# I will reformat our files to match the ENIGMA code

## NOTE: Rename Imaging_examine to 'LeftThickness' or similiar


LeftThickness <- Imaging_examine %>% select(c("ID", starts_with("lh_")))
RightThickness <- Imaging_examine %>% select(c("ID", starts_with("rh_")))
LeftArea <- Imaging_examine %>% select(c("ID", starts_with("lh_")))
RightArea <- Imaging_examine %>% select(c("ID", starts_with("rh_")))

## NOTE: First calculate LEFT ----------------------------#

# Calculate the sum of left frontal cortical thickness
LeftThickness$LeftFrontal <- rowMeans(LeftThickness[,c("lh_superiorfrontal_thickness","lh_rostralmiddlefrontal_thickness","lh_caudalmiddlefrontal_thickness","lh_parsopercularis_thickness","lh_parstriangularis_thickness","lh_parsorbitalis_thickness","lh_lateralorbitofrontal_thickness","lh_medialorbitofrontal_thickness","lh_precentral_thickness","lh_paracentral_thickness","lh_frontalpole_thickness")])

# Calculate the sum of left temporal cortical thickness
LeftThickness$LeftTemporal <- rowMeans(LeftThickness[,c("lh_superiortemporal_thickness","lh_bankssts_thickness","lh_fusiform_thickness","lh_transversetemporal_thickness","lh_entorhinal_thickness","lh_temporalpole_thickness","lh_parahippocampal_thickness")])

# Calculate the sum of left occipital cortical thickness
LeftThickness$LeftOccipital <- rowMeans(LeftThickness[,c("lh_lateraloccipital_thickness","lh_lingual_thickness","lh_cuneus_thickness","lh_pericalcarine_thickness")])

# Calculate the sum of left parietal cortical thickness
LeftThickness$LeftParietal <- rowMeans(LeftThickness[,c("lh_superiorparietal_thickness","lh_inferiorparietal_thickness","lh_supramarginal_thickness","lh_postcentral_thickness","lh_precuneus_thickness")])

# Calculate the sum of left cingulate cortical thickness
LeftThickness$LeftCingulate <- rowMeans(LeftThickness[,c("lh_rostralanteriorcingulate_thickness","lh_caudalanteriorcingulate_thickness","lh_posteriorcingulate_thickness","lh_isthmuscingulate_thickness")])

# To include the insula, as it doesn't belong to the lobes mentioned above
LeftThickness$LeftInsula <- LeftThickness$lh_insula_thickness

## NOTE: exactly the same but for RIGHT ----------------------------#
# Calculate the sum of right frontal cortical thickness
RightThickness$RightFrontal <- rowMeans(RightThickness[,c("rh_superiorfrontal_thickness","rh_rostralmiddlefrontal_thickness","rh_caudalmiddlefrontal_thickness","rh_parsopercularis_thickness","rh_parstriangularis_thickness","rh_parsorbitalis_thickness","rh_lateralorbitofrontal_thickness","rh_medialorbitofrontal_thickness","rh_precentral_thickness","rh_paracentral_thickness","rh_frontalpole_thickness")])
# Calculate the sum of right temporal cortical thickness
RightThickness$RightTemporal <- rowMeans(RightThickness[,c("rh_superiortemporal_thickness","rh_bankssts_thickness","rh_fusiform_thickness","rh_transversetemporal_thickness","rh_entorhinal_thickness","rh_temporalpole_thickness","rh_parahippocampal_thickness")])
# Calculate the sum of right occipital cortical thickness
RightThickness$RightOccipital <- rowMeans(RightThickness[,c("rh_lateraloccipital_thickness","rh_lingual_thickness","rh_cuneus_thickness","rh_pericalcarine_thickness")])
# Calculate the sum of right parietal cortical thickness
RightThickness$RightParietal <- rowMeans(RightThickness[,c("rh_superiorparietal_thickness","rh_inferiorparietal_thickness","rh_supramarginal_thickness","rh_postcentral_thickness","rh_precuneus_thickness")])
# Calculate the sum of right cingulate cortical thickness
RightThickness$RightCingulate <- rowMeans(RightThickness[,c("rh_rostralanteriorcingulate_thickness","rh_caudalanteriorcingulate_thickness","rh_posteriorcingulate_thickness","rh_isthmuscingulate_thickness")])
# To include the insula, as it doesn't belong to the lobes mentioned above
RightThickness$RightInsula <- RightThickness$rh_insula_thickness

# Create a dataset to merge cortical thickness of left and right hemispheres
MergedThickness <- merge(LeftThickness,
                         RightThickness,
                         by="ID")

#MergedThickness <- MergedThickness %>% rename(Subject = ID)
names(MergedThickness)[names(MergedThickness) == "ID"] <- "Subject"

### Save the mean cortical thickness ----
MeanThickness <- data.frame(Subject=MergedThickness$Subject)

# Calculate the mean cortical thickness across the left and right hemispheres
MeanThickness$MeanFrontalThickness <- rowMeans(MergedThickness[,c("LeftFrontal","RightFrontal")])
MeanThickness$MeanTemporalThickness <- rowMeans(MergedThickness[,c("LeftTemporal","RightTemporal")])
MeanThickness$MeanOccipitalThickness <- rowMeans(MergedThickness[,c("LeftOccipital","RightOccipital")])
MeanThickness$MeanParietalThickness <- rowMeans(MergedThickness[,c("LeftParietal","RightParietal")])
MeanThickness$MeanCingulateThickness <- rowMeans(MergedThickness[,c("LeftCingulate","RightCingulate")])
MeanThickness$MeanInsulaThickness <- rowMeans(MergedThickness[,c("LeftInsula","RightInsula")])
MeanThickness$LeftHemisphereThickness <- rowMeans(MergedThickness[,c("LeftFrontal","LeftTemporal","LeftOccipital","LeftParietal","LeftCingulate","LeftInsula")])
MeanThickness$RightHemisphereThickness <- rowMeans(MergedThickness[,c("RightFrontal","RightTemporal","RightOccipital","RightParietal","RightCingulate","RightInsula")])
MeanThickness$MeanThickness <- rowMeans(MeanThickness[,c("LeftHemisphereThickness","RightHemisphereThickness")])

###Read and format cortical surface area ####

# Calculate the sum of left frontal surface area
LeftArea$LeftFrontal <- rowSums(LeftArea[,c("lh_superiorfrontal_area","lh_rostralmiddlefrontal_area","lh_caudalmiddlefrontal_area","lh_parsopercularis_area","lh_parstriangularis_area","lh_parsorbitalis_area","lh_lateralorbitofrontal_area","lh_medialorbitofrontal_area","lh_precentral_area","lh_paracentral_area","lh_frontalpole_area")])
# Calculate the sum of left temporal surface area
LeftArea$LeftTemporal <- rowSums(LeftArea[,c("lh_superiortemporal_area","lh_bankssts_area","lh_fusiform_area","lh_transversetemporal_area","lh_entorhinal_area","lh_temporalpole_area","lh_parahippocampal_area")])
# Calculate the sum of left occipital surface area
LeftArea$LeftOccipital <- rowSums(LeftArea[,c("lh_lateraloccipital_area","lh_lingual_area","lh_cuneus_area","lh_pericalcarine_area")])
# Calculate the sum of left parietal surface area
LeftArea$LeftParietal <- rowSums(LeftArea[,c("lh_superiorparietal_area","lh_inferiorparietal_area","lh_supramarginal_area","lh_postcentral_area","lh_precuneus_area")])
# Calculate the sum of left cingulate surface area
LeftArea$LeftCingulate <- rowSums(LeftArea[,c("lh_rostralanteriorcingulate_area","lh_caudalanteriorcingulate_area","lh_posteriorcingulate_area","lh_isthmuscingulate_area")])
# To include the insula, as it doesn't belong to the lobes mentioned above
LeftArea$LeftInsula <- LeftArea$lh_insula_area


# Calculate the sum of right frontal surface area
RightArea$RightFrontal <- rowSums(RightArea[,c("rh_superiorfrontal_area","rh_rostralmiddlefrontal_area","rh_caudalmiddlefrontal_area","rh_parsopercularis_area","rh_parstriangularis_area","rh_parsorbitalis_area","rh_lateralorbitofrontal_area","rh_medialorbitofrontal_area","rh_precentral_area","rh_paracentral_area","rh_frontalpole_area")])
# Calculate the sum of right temporal surface area
RightArea$RightTemporal <- rowSums(RightArea[,c("rh_superiortemporal_area","rh_bankssts_area","rh_fusiform_area","rh_transversetemporal_area","rh_entorhinal_area","rh_temporalpole_area","rh_parahippocampal_area")])
# Calculate the sum of right occipital surface area
RightArea$RightOccipital <- rowSums(RightArea[,c("rh_lateraloccipital_area","rh_lingual_area","rh_cuneus_area","rh_pericalcarine_area")])
# Calculate the sum of right parietal surface area
RightArea$RightParietal <- rowSums(RightArea[,c("rh_superiorparietal_area","rh_inferiorparietal_area","rh_supramarginal_area","rh_postcentral_area","rh_precuneus_area")])
# Calculate the sum of right cingulate surface area
RightArea$RightCingulate <- rowSums(RightArea[,c("rh_rostralanteriorcingulate_area","rh_caudalanteriorcingulate_area","rh_posteriorcingulate_area","rh_isthmuscingulate_area")])
# To include the insula, as it doesn't belong to the lobes mentioned above
RightArea$RightInsula <- RightArea$rh_insula_area

# Create a dataset to merge surface area of left and right hemispheres
MergedArea <- merge(LeftArea, RightArea, by = "ID")

names(MergedArea)[names(MergedArea) == "ID"] <- "Subject"

### Save the mean surface area ----
TotalArea <- data.frame(Subject=MergedArea$Subject)

# Calculate the mean surface area across the left and right hemispheres
TotalArea$TotalFrontalArea <- rowSums(MergedArea[,c("LeftFrontal","RightFrontal")])
TotalArea$TotalTemporalArea <- rowSums(MergedArea[,c("LeftTemporal","RightTemporal")])
TotalArea$TotalOccipitalArea <- rowSums(MergedArea[,c("LeftOccipital","RightOccipital")])
TotalArea$TotalParietalArea <- rowSums(MergedArea[,c("LeftParietal","RightParietal")])
TotalArea$TotalCingulateArea <- rowSums(MergedArea[,c("LeftCingulate","RightCingulate")])
TotalArea$TotalInsulaArea <- rowSums(MergedArea[,c("LeftInsula","RightInsula")])
TotalArea$LeftHemisphereArea <- rowSums(MergedArea[,c("LeftFrontal","LeftTemporal","LeftOccipital","LeftParietal","LeftCingulate","LeftInsula")])
TotalArea$RightHemisphereArea <- rowSums(MergedArea[,c("RightFrontal","RightTemporal","RightOccipital","RightParietal","RightCingulate","RightInsula")])
TotalArea$TotalArea <- rowSums(TotalArea[,c("LeftHemisphereArea","RightHemisphereArea")])


### Merge and save all cortical surface measures in one dataset ----
CorticalMeasure_ENIGMA <- merge(MeanThickness,TotalArea,by="Subject")

# Save the merged cortical thickness and surface area measures for later use
write.csv(CorticalMeasure_ENIGMA, file="CorticalMeasure_ENIGMA.csv", na="", row.names = F)


## NOTE: now for covariates ----------------------------#

# Work out which covariates file is correct and contains greatest N
# covariates_examine  (n = 591), Subject ID, Sex = male / female
# CorticalCovariates.csv (n = 545)
# We will go with larger n

names(covariates_examine)[names(covariates_examine) == "SubjID"] <- "Subject"

RawCov <- merge(covariates_examine, CorticalMeasure_ENIGMA, by = "Subject")

 # add in ICV too
 names(ICV)[names(ICV) == "SubjID"] <- "Subject"

RawCov <- merge(RawCov, ICV, by = "Subject")

CorticalCovariates_ENIGMA <- merge(covariates_examine, ICV, by = "Subject")

write.csv(CorticalCovariates_ENIGMA, "CorticalCovariates_ENIGMA.csv")
######


#CorticalMeasure <- read.csv("CorticalMeasure.csv")
#CorticalMeasure <- read.csv("CorticalMeasure_ENIGMA.csv")
CorticalMeasure <- CorticalMeasure_ENIGMA

### Save column names of cortical surface measures, including thickness and area ----

CorticalMeasureName <- colnames(CorticalMeasure %>% select(-c("Subject")))

### CHECK THIS LIST ----
CorticalMeasureName

## Read and merge covariates ----
RawCov <- read.csv("./CorticalCovariates_ENIGMA.csv",header=T)

str(RawCov)

### Assume the IID that is equal to the subject ID from cortical and methylation data ----
RawCov <- RawCov[,-c(1)]

### Calculate age^2 and add to the last column of Raw_Cov
RawCov$Age_Square <- (RawCov$Age)^2 


## Merge Methylation covariates in here
#cellcount <- read.csv("Target_3045_with_age_cellcount_setid.csv")

#load(/CCACE_Shared/EleanorC/CorticalEWAS/Data/Input/cell-count.rda)
# PCA

## NOTE: You need to find cell counts, Gran missng
## NOTE: Need to set rownames to Subject / ID to run PCA
#cellcount <- Methcov %>% select("Bcell","CD4T","CD8T","Mono","NK","id")

#tmp <- prcomp(cell-count)

## CHECK STRUCTURE ----
str(RawCov)

############# need meth

### Re-name the ID column for easier merge ----
#names(RawCov)[1] <- c("Subject")
#names(CorticalMeasure)[1] <- c("Subject")

## CHECK STRUCTURE ----
str(RawCov)
str(CorticalMeasure)


load("DNAm_3041_ppts_with_wave_and_cohort_13May2016.RData")
data$id = data$SampleID.1
data = data[,c(14,15,16,17,18,59)]

#save(data, file = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Input/cell-count.rda")

d <- read.csv("Target_3045_with_age_cellcount_setid.csv", header=T)
d1 = d[,c(1,9,11,12,32,36,34,33,35,31,37)]
d1$id = paste(d1$Sentrix.Barcode, d1$Sample.Section, sep = "_") 
d1w2<-d1[which(d1$wave=='WAVE2'),]
colnames(d1w2)[1] <- "Subject"

Methcov = merge(data,d1w2,by=c("id"))
Methcov = Methcov[,c(1,7,2,3,4,5,6,9,12,10,15)]

## CHECK STRUCTURE ----
str(Methcov)
#Methcov2 <- Methcov %>% rename(Subject = id)

#names(Struc)[1] <- c("Subject")

#Cov <- merge(Raw_Cov, ICV, by="Subject", all=F) #combine covariates and ICV		   
Cov <- RawCov
Cov <- merge(Cov, Methcov, by="Subject", all=F) #combine covariates and methylation data covs

## NOTE: id is a huge factor variable - change now? ----
#CovName <- colnames(Cov)[-c(1)] #names of covariates

 ### NOTE: IDEA 1 change factored variables into character variables? ----

 Cov$id <- as.character(Cov$id)
 Cov$ICV <- as.numeric(Cov$ICV)
 Cov$Sex <- as.numeric(Cov$Sex)

 str(Cov)

 ### NOTE: IDEA 2 REMOVE ID and any other misc variables from Cov that could be contributing to the factor issue ----
Cov <- Cov %>% select(-c("GeneticID", "DiseaseType","PID",  "MID"    ,        "DiseaseType"   ,
 "AffectionStatus"    ,   "Affectionstatus2"  ,             "Sentrix.Barcode"   ,    "Plate"      ,           "Sample.Section"    ,
"Infinium.Date"      , "ScannerSite"))

 str(Cov)

CovName <- colnames(Cov) #names of covariates
CovName

### and now to remove those probes/samples from the methylation data frame (here, it is called x) ###

samp <- read.table("Samples_to_remove_mort.txt")
probe <- read.table("Probes_to_remove_mort.txt")


### NOTE: Loading this methylation dataframe will take some time, approx 3 minutes ###
load("Beta_3045_norm_bgcorrect_0.001BetaThreshold.RObject")


### NOTE: At the moment, Methy is our problem; 0 rows ###

# Load quantile normalised meth data - approx 1-2 minutes
Methy = as.data.frame(t(bnorm))

### NOTE: At the moment, Methy rownames = id (101130760149_R04C02 etc) ###
### NOTE: At the moment, Methy colnames = CpG (cg00011200 etc) ###

# add  ID to the meth data ...
#Methy$id = rownames(Methy)

#mk1 <- which(rownames(Methy) %in% samp[,1])
#mk2 <- which(colnames(Methy) %in% probe[,1]) # no probes that fail

#Methy <- Methy[-mk1,mk2]  #,-mk2

#ProbeInvar<-(rowSums(bnorm<=0.2)==ncol(bnorm))|(rowSums(bnorm>=0.8)==ncol(bnorm))
#ListInvarProbe<-rownames(bnorm)[which(ProbeInvar)]


### NOTE: At the moment, Methy rownames = id (101130760149_R04C02 etc) ###
# add  ID to the meth data ...
# Methy$id <- colnames(bnorm) #EITHER, this also works
Methy$id = rownames(Methy)

mk1 <- which(rownames(Methy) %in% samp[,1])
mk2 <- which(colnames(Methy) %in% probe[,1]) # no probes that fail

#Num_Samp <- nrow(mk1)

## NOTE: Methy =  'data.frame': 3045 obs. of  485513 variables;
## NOTE: Methy2 = 'data.frame':	3045 obs. of  485121 variables;
## NOTE: Methy3 = 'data.frame':	2653 obs. of  3143 variables (nb. some NaNs spotted)
## NOTE: Methy4 = 'data.frame':	3045 obs. of  482370 variables:


Methy2 <- Methy[-mk1]  #,-mk2
Methy3 <- Methy[-mk1,mk2]  #,-mk2
Methy4 <- Methy[-mk2]

## CHECK STRUCTURE ----
str(Methy)
str(Methy2)
str(Methy3)
str(Methy4)

## CHECK NaNs ----
library(dplyr)
nan_check <- Methy3 %>% summarise_all(~ sum(is.na(.)))

## IMPUTE NaNs or remove NaNs?
## Option 1: REMOVE
## NOTE: Methy_NaN_remove = 'data.frame':	140 obs. of  3143 variables:
## NOTE: Methy_NaN_remove = 'data.frame': rownames = id (3998499072_R04C02), colnames = CpG (cg00011200)

Methy_NaN_remove <- Methy3[complete.cases(Methy3), ]

str(Methy_NaN_remove)

## FOR IMMEDIATE CONSISTENCY WITH MICHELLLE CODE ----
#Methy <- Methy2

## FOR IMMEDIATE ATTEMPT TO RUN aka SMALLER no regressions ----
#Methy <- Methy3
Methy <- Methy_NaN_remove

## NOTE: Takes ~ 1-2 minutes to run
ProbeInvar<-(rowSums(bnorm<=0.2)==ncol(bnorm))|(rowSums(bnorm>=0.8)==ncol(bnorm))
ListInvarProbe<-rownames(bnorm)[which(ProbeInvar)]

ListInvarProbe


###Generating Files with Complete Data across All Data#####
## CHECK STRUCTURE ----
str(Methy) #'data.frame':	140 obs. of  3143 variables:
str(CorticalMeasure) #'data.frame':	652 obs. of  71 variables: (Subject included)
str(Cov) #'data.frame':	504 obs. of  11 variables: (Subject and id still included)

Data <- merge(Cov, CorticalMeasure, by="Subject", all=F) #combine covariates and structure data
# NOTE: there needs to be an id column in the methy data
Methy$id = rownames(Methy)
Match <- match(Methy$id, Data$id) #match methylation sample with combined covariates and structure sample

## CHECK STRUCTURE ----
str(Data) #'data.frame':	499 obs. of  81 variables:
#str(Match)

Cov <- Data[Match[!is.na(Match)],2:(length(CovName)+1)]  #generate covariates with matched individual
CorticalMeasure <- Data[Match[!is.na(Match)],-c(1:(length(CovName)+1))] #generate structure data with matched individual
Methy_Matched <- Methy[!is.na(Match),-c(1)] #generate methylation data with matched individual

## CHECK STRUCTURE ----
str(Cov) # NOTE: 'data.frame':	499 obs. of  11 variables: (17 obs)
str(CorticalMeasure) #NOTE: 'data.frame':	499 obs. of  69 variables: (17 obs)
str(Methy_Matched) # NOTE: 'data.frame':	499 obs. of  485120 variables: (17 obs. of  3143 )

Methy <- Methy_Matched

#### Generate a dataframe to save age information for further analysis
AgeInfo <- data.frame(min(Cov$Age),max(Cov$Age),mean(Cov$Age))
names(AgeInfo) <- c("MinAge","MaxAge","MeanAge")

## CHECK colnames ----
colnames(Cov) # consider changing GeneticID to Subject
colnames(CorticalMeasure)
#colnames(Methy)

## NOTE: with cortical measures, age and sex are still included, consider removing ----
Cov <- Cov %>% select(-c("id", "Mean_bankssts_surfavg"))

### REMOVE AGE AND SEX AS WE ADD THIS IN WITH THE COVARIATES LATER ----
#CorticalMeasure <- CorticalMeasure %>% select(-c("sex", "ageyrs"))

## CHECK colnames ----

colnames(Cov)

## cell counts PCA  ----




## EWAS with cortical measures ----

# NOTE: problem Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) :
 # contrasts can be applied only to factors with 2 or more levels

 # NOTE: initial troublshoot:
 # Covar = contains list of covarites, but we still have ID etc in here
 # ID is a factor variable with many levels, so two ideas
 # Idea 1: turn ID into a character variable
 # Idea 2: Remove ID from Covar / Cov entirely

 ### NOTE: IDEA 2 REMOVE ID and any other misc variables from Cov that could be contributing to the factor issue ----
#Cov <- Cov %>% select(-c("id", "GeneticID", "DiseaseType","PID",  "MID"    ,        "DiseaseType"   ,
 #"AffectionStatus"    ,   "Affectionstatus2"  ,             "Sentrix.Barcode"   ,    "Plate"      ,           "Sample.Section"    ,
#"Infinium.Date"      ,   "Mean_bankssts_surfavg", "ScannerSite"))

## EWAS with cortical measures ----
# Get the numbers of methylation probes, cortical measures and covariates
Num_Methy <- ncol(Methy)
Num_Cov <- ncol(Cov)
Num_Cortical <- ncol(CorticalMeasure)

Num_Methy
Num_Cov
Num_Cortical

## CHECK colnames ----

colnames(Cov)

##################################################################
##Association analysis: the following section should be run by ALL cohorts
##################################################################

######## EWAS with cortical measures: whole sample, controlling for ICV ##########

#CorticalMeasure <- Struc

## Prepare the output files for the cortical EWAS----
## NOTE: Here we are setting up space for our findings, with NA values

Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)

Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)

Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)

Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

#save(Origin_Beta, file = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Origin_Beta.csv")


## linear regression
## NOTE: The code here is part of a for loop that's used for performing linear regression analysis between DNA methylation data (Methy) and subcortical brain volume data (CorticalMeasure) for each combination of cortical brain volume and DNA methylation variable. This loop is structured as follows:

## Here's what the loop is doing:

##    (1). The outer loop iterates over each subcortical brain volume (indexed by i), and the inner loop iterates over each DNA methylation variable (indexed by j).

##    (2) For each combination of subcortical brain volume (i) and DNA methylation variable (j), it performs a linear regression analysis using the lm function.
##          The model is specified ## as Methy[, i] ~ Covar + CorticalMeasure[, j],
##          (a) where CorticalMeasure[, j] represents the cortical brain volume,
##          (b) Covar represents covariates,
##          (c) and Methy[, i] represents the DNA methylation variable.

##    The results of the linear regression analysis (including beta coefficients, standard errors, and p-values) are stored in the Out object.

 ##   The loop updates the Origin_Beta, Origin_SE, and Origin_P matrices with the
 ##     (a) beta coefficient,
 ##     (b) standard error,
 ##     (c) and p-value from the regression analysis, respectively.
 ## These matrices store the results for each combination of cortical brain volume and DNA methylation variable.

 ## NOTE:   THIS TAKES APPROX 10 MINUTES TO RUN

Covar <- matrix(unlist(Cov), ncol=ncol(Cov), byrow=F)
system.time({
for (i in 1:Num_Methy) {
  for (j in 1:Num_Cortical) {
    Out <- summary(lm(Methy[,i] ~ Covar + CorticalMeasure[,j]))
    Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
    Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
    Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
    Origin_N[i,j] <- Out$df[1]+Out$df[2]
  }
}
})
###########################################################################
###########################################################################

## NOTE: error message = Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) :
##  NA/NaN/Inf in 'y'
##In addition: Warning message:
##In storage.mode(v) <- "double" : NAs introduced by coercion


## Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov),collapse = "+"),"+CorticalMeasure[,j]")

## Save the output file of EWAS with cortical measures, in the whole sample ----
#save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures.RData",sep=""),compress=T)

save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Output_of_",cohort,"_EWAS_with_CorticalMeasures.RData",sep=""),compress=T)


######### EWAS with cortical measures: whole sample, not controlling for ICV ########
## Prepare the output files for the cortical EWAS----
Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)
Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)
Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)
Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

## linear regression
Cov_noICV <- Cov[, !(names(Cov) %in% c("ICV"))]
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)

for (i in 1:Num_Methy){
  for (j in 1:Num_Cortical){
    Out <- summary(lm(Methy[,i]~Covar + CorticalMeasure[,j]))
    Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
    Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
    Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
    Origin_N[i,j] <- Out$df[1]+Out$df[2]
  }
}

## Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov_noICV),collapse = "+"),"+CorticalMeasure[,j]")


## Save the output file of EWAS with cortical measures, in the whole sample ----
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Output_of_",cohort,"_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)




################ EWAS with cortical measures: stratified by sex ################
### Please make sure Sex was specified as follows: Males=1, Females=2
Num_Methy <- ncol(Methy) 
Num_Cov <- ncol(Cov_noICV)
Num_Cortical <- ncol(CorticalMeasure)

Gender <- colnames(Cov_noICV)=="Sex"
Category <- c("Male","Female")

### Generate a dataframe to save age information for further analysis
# to get mean ages for male and female samples separately
MaleCov <- subset(Cov_noICV, Sex==1)
MaleAgeInfo <- data.frame(min(MaleCov$Age),max(MaleCov$Age),mean(MaleCov$Age))
names(MaleAgeInfo) <- c("MaleMinAge","MaleMaxAge","MaleMeanAge")
FemaleCov <- subset(Cov_noICV, Sex==2)
FemaleAgeInfo <- data.frame(min(FemaleCov$Age),max(FemaleCov$Age),mean(FemaleCov$Age))
names(FemaleAgeInfo) <- c("FemaleMinAge","FemaleMaxAge","FemaleMeanAge")

# to run the stratified EWAS analysis
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)

### linear regression
for (k in 1:2) {
  Ind <- Covar[,Gender]==(k) 
  Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_Beta) <- colnames(CorticalMeasure)
  rownames(Origin_Beta) <- colnames(Methy)
  Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_SE) <- colnames(CorticalMeasure)
  rownames(Origin_SE) <- colnames(Methy)
  Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_P) <- colnames(CorticalMeasure)
  rownames(Origin_P) <- colnames(Methy)
  Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_N) <- colnames(CorticalMeasure)
  rownames(Origin_N) <- colnames(Methy) 
  for (i in 1:Num_Methy) {
    for (j in 1:Num_Cortical) {
      ### Remove gender from covariates
      Out <- summary(lm(Methy[Ind,i]~Covar[Ind,!Gender]+CorticalMeasure[Ind,j]))         
      Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
      Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
      Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
      Origin_N[i,j] <- Out$df[1]+Out$df[2]
    }
  }
  
  
  ###########################################################################
  ###########################################################################
  
  
  ### Save model description
  model <- paste0("Methy[Ind,i] ~ ",paste0(names(Cov_noICV)[!Gender],collapse = "+"),"+CorticalMeasure[Ind,j]")
  
  ### Save output file of EWAS with cortical measures, when the EWAS was conducted in each gender separately
  save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,MaleAgeInfo,FemaleAgeInfo,file=paste("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Output_of_",cohort,"_",Category[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
}



