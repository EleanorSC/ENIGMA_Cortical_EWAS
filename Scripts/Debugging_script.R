# /CCACE_Shared/EleanorC/CorticalEWAS/Data/Input

setwd("/CCACE_Shared/EleanorC/CorticalEWAS/Data")


install.packages("tidyverse")
library(tidyverse)

cohort = "LBC1936"

## ----------------------------# 

CorticalMeasure <- read.csv("CorticalMeasure.csv")

### Save column names of cortical surface measures, including thickness and area ----
#CorticalMeasureName <- colnames(CorticalMeasure)[-1]

### Save column names of cortical surface measures, including thickness and area ----
### ENSURE TO REMOVE SUBJECT, AGE AND SEX ----
#CorticalMeasureName <- colnames(CorticalMeasure)[,-c(1,72,73)]
#CorticalMeasureName <- colnames(CorticalMeasure)
#CorticalMeasureName <- colnames(CorticalMeasure)[-c("SubjID", "sex", "ageyrs")]
#CorticalMeasureName <- colnames(CorticalMeasure) %>% select(-c("SubjID", "sex", "ageyrs"))
#CorticalMeasureName <- colnames(CorticalMeasure)[,-which(names(CorticalMeasure) %in% c("SubjID", "sex", "ageyrs"))]

CorticalMeasureName <- colnames(CorticalMeasure %>% select(-c("SubjID", "sex", "ageyrs")))

### CHECK THIS LIST ----
CorticalMeasureName

### REMOVE AGE AND SEX AS WE ADD THIS IN WITH THE COVARIATES LATER ----
CorticalMeasure <- CorticalMeasure %>% select(-c("sex", "ageyrs"))

## Read and merge covariates ----
RawCov <- read.csv("./CorticalCovariates.csv",header=T)

### Assume the IID that is equal to the subject ID from cortical and methylation data ----
RawCov <- RawCov[,-c(1)]

### Calculate age^2 and add to the last column of Raw_Cov
RawCov$Age_Square <- (RawCov$Age)^2 


## CHECK STRUCTURE ----
str(RawCov)

############# need meth

### Re-name the ID column for easier merge ----
names(RawCov)[1] <- c("Subject")
names(CorticalMeasure)[1] <- c("Subject")

## CHECK STRUCTURE ----
str(RawCov)
str(CorticalMeasure)

#Struc <- CorticalMeasure


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

