## ---------------------------

setwd("/Users/eleanorc_worklaptop/desktop/ENIGMA/DATA")


## ----------------------------# 
# install relevant packages

cohort = "LBC1936"

## ----------------------------# 

CorticalMeasure <- read.csv("CorticalMeasure.csv")

### Save column names of cortical surface measures, including thickness and area ----
CorticalMeasureName <- colnames(CorticalMeasure)[-1]

## Read and merge covariates ----
RawCov <- read.csv("./CorticalCovariates.csv",header=T)

### Assume the IID that is equal to the subject ID from cortical and methylation data ----
RawCov <- RawCov[,-c(1)]

### Calculate age^2 and add to the last column of Raw_Cov
RawCov$Age_Square <- (RawCov$Age)^2 

############# need meth

### Re-name the ID column for easier merge ----
names(RawCov)[1] <- c("Subject")


Struc <- CorticalMeasure

### Re-name the ID column for easier merge ----
names(Struc)[1] <- c("Subject")

load("DNAm_3041_ppts_with_wave_and_cohort_13May2016.RData")
data$id = data$SampleID.1
data = data[,c(14,15,16,17,18,59)]


d <- read.csv("Target_3045_with_age_cellcount_setid.csv", header=T)
d1 = d[,c(1,9,11,12,32,36,34,33,35,31,37)]
d1$id = paste(d1$Sentrix.Barcode, d1$Sample.Section, sep = "_") 
d1w2<-d1[which(d1$wave=='WAVE2'),]
colnames(d1w2)[1] <- "Subject"

Methcov = merge(data,d1w2,by=c("id"))
Methcov = Methcov[,c(1,7,2,3,4,5,6,9,12,10,15)]


#Cov <- merge(Raw_Cov, ICV, by="Subject", all=F) #combine covariates and ICV		   
Cov <- RawCov
Cov <- merge(Cov, Methcov, by="Subject", all=F) #combine covariates and methylation data covs

CovName <- colnames(Cov)[-c(1)] #names of covariates



### and now to remove those probes/samples from the methylation data frame (here, it is called x) ###

samp <- read.table("Samples_to_remove_mort.txt")
probe <- read.table("Probes_to_remove_mort.txt")

load("Beta_3045_norm_bgcorrect_0.001BetaThreshold.RObject")
x = as.data.frame(t(bnorm))
x$id = rownames(x)

mk1 <- which(rownames(x) %in% samp[,1])
# mk2 <- which(colnames(x) %in% probe[,1]) # no probes that fail

x1 <- x[-mk1]  #,-mk2

ProbeInvar<-(rowSums(bnorm<=0.2)==ncol(bnorm))|(rowSums(bnorm>=0.8)==ncol(bnorm))
ListInvarProbe<-rownames(bnorm)[which(ProbeInvar)]

###Generating Files with Complete Data across All Data#####
Data <- merge(Cov, Struc, by="Subject", all=F) #combine covariates and structure data
Match <- match(x1$id, Data$id) #match methylation sample with combined covariates and structure sample
Cov <- Data[Match[!is.na(Match)],2:(length(CovName)+1)]  #generate covariates with matched individual
Struc <- Data[Match[!is.na(Match)],-c(1:(length(CovName)+1))] #generate structure data with matched individual
Methy <- x1[!is.na(Match),-c(1)] #generate methylation data with matched individual

##########
write.csv(Data, "Data.csv")
write.csv(Match, "Match.csv")
write.csv(Cov, "Cov.csv")
write.csv(Struc, "Struc.csv")
write.csv(Methy, "Methy.csv")


## ----------------------------# 
# attempt from where code failed on server

## ----------------------------# 



### Generate a dataframe to save age information for further analysis
AgeInfo <- data.frame(min(Cov$Age),max(Cov$Age),mean(Cov$Age))
names(AgeInfo) <- c("MinAge","MaxAge","MeanAge")

## EWAS with cortical measures ----
# Get the numbers of methylation probes, cortical measures and covariates
Num_Methy <- ncol(Methy) 
Num_Cov <- ncol(Cov)
Num_Cortical <- ncol(CorticalMeasure)


##################################################################
##Association analysis: the following section should be run by ALL cohorts 
##################################################################

######## EWAS with cortical measures: whole sample, controlling for ICV ##########

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
Covar <- matrix(unlist(Cov), ncol=ncol(Cov), byrow=F)
for (i in 1:Num_Methy) {
  for (j in 1:Num_Cortical) {
    Out <- summary(lm(Methy[,i] ~ Covar + CorticalMeasure[,j]))
    Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
    Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
    Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
    Origin_N[i,j] <- Out$df[1]+Out$df[2]
  }
}

###########################################################################
###########################################################################


## Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov),collapse = "+"),"+CorticalMeasure[,j]")

## Save the output file of EWAS with cortical measures, in the whole sample ----
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures.RData",sep=""),compress=T)



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
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)






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
  save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,MaleAgeInfo,FemaleAgeInfo,file=paste("./Output_of_",cohort,"_",Category[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
}

