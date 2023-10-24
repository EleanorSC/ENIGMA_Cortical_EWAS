
setwd("/CCACE_Shared/EleanorC/CorticalEWAS/Data")

#############################
### DNAm idats from LBC36 ###
#################''##########

install.packages("tidyverse")
library(tidyverse)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")


#To install the FlowSorted.Blood.450k package, which is necessary when using the function of “estimateCellCounts”, enter:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FlowSorted.Blood.450k")

#Highlighted portions of the instructions require you to make changes so that the commands work on your system and data.

##open R and copy the lines below
require(minfi)
require(minfiData)

library(RColorBrewer)
library(minfi)
library(limma)
library(FlowSorted.Blood.450k)

setwd("/CCACE_Shared/EleanorC/CorticalEWAS/Data")  # replace with your local working directory

#################################
### Reading the methylation data ####
#################################

## Creating the initial object of the minfi analysis that contains the raw intensities in the green and red channels
# Set your data directory containing IDAT files
idat_dirs <- c("/CCACE_Shared/EleanorC/CorticalEWAS/Data/idats")

#specify your local directory, which contains both raw IDAT files and the sample sheets (i.e., a .csv file, which includes the path for each sample’s IDAT file). The following scripts expect the sample sheet to include column names as shown in the figure below, where the header information (all lines up to [Data]) could be omitted (Variables such as sex, age, site and disease status can also be included as extra columns):

#####################################
### DNAm idats from LBC36 CREATED ###
#####################################

## Try creating RGset
#library(minfi)

targets=read.metharray.sheet(idat_dirs, "Deary Meth450K_plates 9_21_7-SampleSheet1220413.csv", recursive=T) #loads the corresponding .csv sample sheet file


targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates6-18-10_230413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates15_14_19_SampleSheet220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates16-08-12_220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates20_5_SampleSheet_220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates23_1_4_Angie_220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates25-26-13_220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates_2-11-22_220413.csv", recursive=T))
targets=rbind(targets, read.metharray.sheet(idat_dirs, "E11970_Meth450K_plates_17-3-24__220413.csv", recursive=T))


## NOTE: 1406 obs. of  8 variables ... not right?
## $ Sample_Name : chr  "LBC360070_36_W2" "LBC360055_36_W2" "LBC360056_36_W2" "LBC0242M_21_W4" ...
# $ Sample_Well : chr  NA NA NA NA ...
# $ Sample_Plate: chr  NA NA NA NA ...
# $ Sample_Group: logi  NA NA NA NA NA NA ...
# $ Pool_ID     : chr  NA NA NA NA ...
# $ Array       : chr  "R02C01" "R03C01" "R04C01" "R05C01" ...
# $ Slide       : chr  "8963303040" "8963303040" "8963303040" "8963303040" ...
# $ Basename    : chr  "/CCACE_Shared/EleanorC/CorticalEWAS/Data/idats/8963303040_R02C01" "/CCACE_Shared/EleanorC/CorticalEWAS/Data/idats/8963303040_R03C01" "/CCACE_Shared/# EleanorC/CorticalEWAS/Data/idats/8963303040_R04C01" "/CCACE_Shared/EleanorC/CorticalEWAS/Data/idats/8963303040_R05C01" ...


#targets=read.metharray.sheet(idat_dirs, "Sample_Sheets.csv", recursive=T) #loads the corresponding .csv sample sheet file
## NOTE incomplete final line found on '/CCACE_Shared/EleanorC/CorticalEWAS/Data/idats/Sample_Sheets_tidy.csv'

## NOTE: basename contains the file paths
# /CCACE_Shared/EleanorC/CorticalEWAS/Data/idats/8963303040_R03C01"

## 1. check for any duplicates in basenames

#targets <- targets[!duplicated(targets$Basename)]
targets <- targets %>% filter(Basename != "character(0)")

targets <- targets %>% filter(grepl("_W2", Sample_Name))



#### targets
## NOTE: later need sex variable included, so potentially use the targets_sex which contains cell types + sex info
#target_sex <- read.csv("targets_W2_sex.csv")

targets <- read.csv("new_targets_W2.csv")

RGset <- read.metharray.exp(
  base = idat_dirs,
  targets = targets,
  verbose = T
)

##NOTE this takes 20 minutes to run
save(RGset,file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/RGset.rda") #saves the object

### If the initial object with raw methylation data had already been created as above, it can be directly loaded with the command below (if the file hasn’t been generated yet or if you don’t know what it means, please follow the procedure above):
# If you accidentially kill the session
load("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/RGset.rda")

### Producing Quality Control plots
pd <- pData(RGset) #extracting the sample information (phenotype data) from the sample sheet

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kmanifest")

## NOTE: This took 10 minutes to run

qcReport(RGset, sampNames = pd$Sentrix, sampGroups=pd$set, pdf = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/qcReport.pdf") # produces a PDF QC report of common plots, colored by groups of samples. Check that the highlighted variables correspond to column names in your sample sheet, as illustrated above. In this example, the column named ‘Experimenter’ corresponds to the different waves by which our DNA samples were processed and hybridized.

#qcReport(RGset, sampNames = NULL, sampGroups=pd$set, pdf = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/qcReport_abrev.pdf", maxSamplesPerPage = 24)

#These plots are useful for identifying samples with data quality. They display summaries of signals from the array (e.g. density plots) as well as the values of several types of control probes included on the array. A good rule of thumb is to be wary of samples whose behaviour deviates from that of others in the same or similar experiments. The following figures display different types of plots, with no evidence for outliers in 6 samples.
# In case outliers are detected, their individuals’ Sample.ID should registered into the file ‘Outliers’, which will be used to remove outliers in a later step.

##################################
### Quality assessment of methylation data ####
##################################

object=preprocessIllumina(RGset) #We prefer preprocessIllumina because it considers the background correction as well as normalization to internal controls as described in the Illumina documentation. This will minimize the amount of variation between arrays.

object<-mapToGenome(object) #assign probes with its physical location on the genome
object=ratioConvert(object, type="Illumina") #convert raw methylation data
beta <- getBeta(object) #get the beta value for each probe
dat <- object #rename ‘object’ to ‘dat’
pd=pData(dat) #get phenotypes of methylation data


####Removing data from the X and Y chromosomes#######
keepIndex=which(!seqnames(dat)%in%c("chrX","chrY"))
beta <- beta[keepIndex,]

### Use MDS to check bad samples ####
mydist=dist(t(beta[sample(1:nrow(beta),10000),]))
mds=cmdscale(mydist)


##NOTE: Snetrix ID = Slide, Sentrix Position = Array; in our target data, 'array' is equivalent to 'Slide' and 'pos' = array
pdf("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Batch_effect_mds.pdf",width=11,height=4)
plot(mds[,1],col=as.numeric(as.factor(pd[,"array"]))+1,xlab="",ylab="First Principal component",xaxt="n")  #can replace with your own potential confounder variables.
dev.off()

### Use PCA to check batch effects ####
b <- beta - rowMeans(beta)

#Performing singular value decomposition (SVD), equivalent to principal component analysis with covariance matrix.
install.packages("corpcor")
library(corpcor)

##NOTE: takes ~15 minutes
ss <- fast.svd(b,tol=0) # this may take a while

#Looking at the percent of variance explained by each principal component
percvar <- ss$d^2/sum(ss$d^2) #calculates the variance explained by each component
pdf(file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/PCA_distribution.pdf")
plot(percvar,xlab="Principal Components",ylab="Variance explained")
dev.off()
#example of output



## Plotting Batch Effect over Methylation Components
for (i in c("pos","array")){ #can add your own potential confounder variables.
	pdf(file=paste("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/batch_",i,".pdf",sep=""))  #output filename
	variabletocolor=pd[,i] #assign different colours to batches
	bob=levels(factor(variabletocolor))
	colors=match(variabletocolor,bob)
	pairs(ss$v[,1:4], col=colors, labels=c("PC1", "PC2", "PC3", "PC4"))
	par(mfrow=c(2,2))
	boxplot(ss$v[,1]~pd[,i],ylab="PC1",xlab=i)
	boxplot(ss$v[,2]~pd[,i],ylab="PC2",xlab=i)
	boxplot(ss$v[,3]~pd[,i],ylab="PC3",xlab=i)
	boxplot(ss$v[,4]~pd[,i],ylab="PC4",xlab=i)
	dev.off()
}

#An example of the outputs:

#this figure shows no clear subpopulation based on the first 4 principal components

#this figure shows the box plot of slides effect against the first 4 principal components. It appears that the slide has some, but no gigantic, impact on PC1, but a huge effect on PC3, i.e. individuals on different slides show significantly different distributions on PC3.


### Mark individuals out of normal range (i.e. median+3SD or median-3SD) based on first 4 components ###
RM <- rep(FALSE, nrow(ss$v))
for (i in 1:4) {
Median <- median(ss$v[,i]) #calculate median of each component
SD <- sd(ss$v[,i]) #calculate the standard deviation of each component
RM <- RM|(abs(ss$v[,i]-Median)>3*SD) #mark individuals outside 3SD range as TRUE
}
# RM will be used to remove outliers at a later stage

##NOTE this doesn't work - no sex data encoded?

target_sex <- read.csv("targets_W2_sex.csv")
target_sex <- target_sex %>% select(Basename, age,   neut,   lymph,  mono,  eosin,  baso,  sex)

pd <- merge(pd, target_sex, by = "Basename")

###Predicting sex using methylation data###
predictedSex <- getSex(dat, cutoff = -2)$predictedSex # N.B., this function does not handle datasets with only females or only males)
pdf(file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Sex_Plot.pdf")
plotSex(getSex(dat, cutoff = -2))  #Plot of predicted gender information as follows
dev.off()

pdf(file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/predicted_gender.pdf")

Jitter1 <- jitter(as.numeric(as.factor(pd$sex))) #Here we assume the self-reported gender information is integrated in the RGset and therefore is available as pd$Sex. Otherwise, please replace pd$Sex with a vector containing the self-reported gender information.  Ideally, both self-reported and predictedSex should take the same format, i.e. “M” for male and “F” for female.

Jitter2 <- jitter(as.numeric(as.factor(predictedSex)))
plot(Jitter1, Jitter2,xlab="Sex (self reported)",ylab="Sex(predicted)",xaxt="n",yaxt="n")
axis(2,c(1,2),c("Female", "Male"))
axis(1,c(1,2),c("Female", "Male"))
#idx <- as.factor(pd[,"Sex"])==as.factor(predictedSex) #Change to your own vector of self-reported sex if pd$Sex is not available.
idx <- as.factor(pd[,"sex"])==as.factor(predictedSex) #Change to your own vector of self-reported sex if pd$Sex is not available.
Outliers <- pd[!idx, “Sentrix”] #mark IDs with mis-matched gender information
# text(Jitter1[!idx],Jitter2[!idx]-0.05,Outliers,cex=0.5) #If included, this line will add IDs to the mis-matched individuals. However, if there are too many outliers, the figure might look a bit messy.
dev.off()

# This figure shows one individual with wrong sex information to be removed from downstream analyses. The IDs of individuals with wrongly assigned sex information were registered into the variable ‘Outliers’.

#####Removing outlier individuals####
##We have listed outlier individuals based on gender discrepancy (graph above) and qcReport in variable ‘Outliers’, and we have also marked outlier individuals outside of 3 standard deviation (3SD) in variable ‘RM’.

##We also need to remove methylation control samples. If a ‘Note’ column in the Sample sheet provides such information, we can use the following code to mark these control individuals:

MC <- pd[,"Notes"] == "Methylation control"  #mark control individuals as TRUE
#Change the highlighted portion to the corresponding column name and identifier in your local Sample sheet .csv file that suits the type of samples you want to remove.

##However it could be the case that the ‘Note’ column does not exist, and therefore one needs to remove Methylation controls based on the ‘Sample.ID’ column, for example, if we know all control individuals share the same string ‘Meth’, we could use the following code:

Index <- grep("Meth", pd[,"Sample.ID"]) #find individuals as methylation control
MC <- rep(FALSE, length(pd[,"Sample.ID"])) #generate an all FALSE variable with length of sample size
MC[Index] <- TRUE #mark control individuals as TRUE
#Change the highlighted portion to the corresponding identifier in your local Sample sheet .csv file.

##Mark and remove all outliers:
RM <- RM|(pd[,"Sample.ID"]%in%Outliers) #combine outliers marked in RM and Outlier
RM <- RM|MC #merge outliers info from MC
pd = pd[-which(RM),] #remove outliers
RGset = read.metharray.exp(base = datadir, targets=pd, verbose=TRUE) #reload data without outliers
save(RGset, file="./RGset.rda") #save raw data with outliers removed

################## Assuming outlier removal

###if samples have been processed in batches, each batch can be assessed separately as above, but the data should be merged before further preprocessing. Example for merging 2 waves of data, and having been saved as ‘RGset1’ and ‘RGset2’ through save(RGset1, file="./RGset1.rda") and save(RGset2, file="./RGset2.rda"), respectively.

p1 <- pData(RGset1)
p2 <- pData(RGset2)
Match <- match(colnames(p2), colnames(p1)) #match columns between p2 and p1
p1  <- p1[,Match[!is.na(Match)]] #in p1, only keep columns where both p1 and p2 have values.
Match <- match(colnames(p1), colnames(p2)) #match columns between p1 and p2
p2  <- p2[,Match[!is.na(Match)]] # in p2, only keep columns where both p1 and p2 have values.
stopifnot(identical(colnames(p1),colnames(p2))) #compare if p1 and p2 have the same column names.
pd <- rbind(p1,p2) #merge p1 and p2 to form the combined phenotype file.

RGset = read.metharray.exp(targets=pd, verbose=TRUE) #reload databased on combined Sample sheet information.
save(RGset, file="./RGset.rda") #save combined RGsets after QC.


################################################
### Preprocessing (Quantile normalization) after initial QC ####
################################################
#We have used the preprocessQuantile function that implements stratified quantile normalization preprocessing for Illumina methylation microarrays (Epigenomics 4, 325-341 (2012). The algorithm used for such a normalization method relies on the assumptions necessary for quantile normalization to be applicable and thus is not recommended for cases where global changes are expected, such as in cancer–normal comparisons. When large-scale differences in methylation are expected (cancer/normal studies or between tissue studies) the functional normalization algorithm may be preferable. Our experience shows that using quantile normalization procedure for blood samples without inflammation disorders works well. Subset-quantile within array normalization (SWAN) is another method that worked quite well and could be recommended, but quantile normalization is more popular.

# Low-quality samples that normalization cannot correct are automatically filtered out using the removeBadSamples argument in the preprocessQuantile function.

##Running quantile normalization
##NOTE: This step takes around ~ 15-20 minutes

object=preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = TRUE, badSampleCutoff = 10.5, quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL, verbose = TRUE)  #Note this step removes samples with low intensity

save(object, file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Quan-norm.rda") #Saving the normalized data to a new file

## Get beta values of normalized data
beta <- getBeta(object) #Gets the normalized Beta Value from normalized data
keepIndex=which(!seqnames(object)%in%c("chrX","chrY")) #mark probles on X and Y chromosome
beta <- beta[keepIndex,] #only keep probes on autosomes
save(beta, file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Beta_Quantile.rda") #saves the normalized beta values to new file

## Generate PCA to control for unknown structure ####

#library(corpcor)
ss <- fast.svd(beta,tol=0) # this may take a while
save(ss,file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/fast_svd.rda")  #saves the output as a new file to be used as covariate in downstream analyses


## load cellcount


load("DNAm_3041_ppts_with_wave_and_cohort_13May2016.RData")
colnames(data)
rownames(data) <- data$SampleID.1
data = data[,c("CD8T", "CD4T", "NK" ,"Bcell", "Mono", "Gran")]

save(data, file = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/cellcount.rda")


### Estimate Cell-type Proportions ####
#fileName <- "./cellcount.rda"  #Output file name

#if (file.exists(fileName)) {
 #   load(fileName) #load existing output. Otherwise generate a new output
#} else {
#	cellcount <- estimateCellCounts(RGset) #estimate cell count
#	save(cellcount, file=fileName) #save output
#}

################################
#######QC plots after preprocessing######
################################

Cohort <- "LBC1936" #Please change to your own cohort name.
mydist=dist(t(beta[sample(1:nrow(beta),10000),]))
mds=cmdscale(mydist)
pd=pData(dat)
pdf(paste("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Batch_effect_Quantile_Normalised_mds_",Cohort, ".pdf",sep=""),width=11,height=4)
plot(mds[,1],col=as.numeric(as.factor(pd[,"array"]))+1,xlab="",ylab="First Principal component",xaxt="n") #highlighted can be changed to other possible batch effects
dev.off()


###########################################################

##Section 1: methylation data QC (skip this step if you have already QC’ed you methylation data following our ENIGMA-Epigenetics QC pipeline)

# Section 1: Methylation data quality check ----
# Prep1: Probe Quality Check ----

library(minfi)

## Load files ---- NOTE: takes approx 2 minutes to run
### Load RGset file ----
load('/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/RGset.rda')
### Load Quan-norm file ----
load('/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Quan-norm.rda')

## Add SNP information to the data NOTE takes approx 10 minutes to run ----
objectWithSNPinfo <- addSnpInfo(object)

## Drop probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension ----
objectSNPQCed <- dropLociWithSnps(objectWithSNPinfo, snps = c("SBE", "CpG", "Probe"), maf = 0.05)
rm(objectWithSNPinfo);

## A detection p-value is returned for every genomic position in every sample ----
detP <- detectionP(RGset)
Match1 <- match(colnames(objectSNPQCed),colnames(detP))
Match2 <- match(rownames(objectSNPQCed),rownames(detP))
detPSNPQCed <- detP[Match2[!is.na(Match2)],Match1[!is.na(Match1)]]
rm(Match1, Match2, detP);

## Positions with non-significant p-values (typically >0.01) should not be trusted ----
failed <- detPSNPQCed >0.01
rm(detPSNPQCed);

## Get beta values ---- NOTE: takes approx 2 minutes
beta <- getBeta(objectSNPQCed)

## Drop probes that failed quality control via the detection p-value in greater than 20% of samples ----
failedCG02 <- rowMeans(failed) >0.2

## Get the list of non-variable CpG sites i.e. those where beta values for all samples are ≤20% or ≥80% ----
#data                  0             1 0.129 FAL: 377417, TRU: 55958
ProbeInvar <- (rowSums(beta<=0.2)==ncol(beta))|(rowSums(beta>=0.8)==ncol(beta))

## Mark probes with either all beta value <=0.2 or all beta value>=0.80 and generate a list of probes marked as invariant ----
#This list will be included in the output file that you will send us
ListInvarProbe <- rownames(beta)[which(ProbeInvar)]
rm(ProbeInvar);

length(ListInvarProbe)
#[1] 55958

################### IF THIS DOESN'T WORK, CONSIDER FOLLOWING ----

samp <- read.table("Samples_to_remove_mort.txt")
probe <- read.table("Probes_to_remove_mort.txt")

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

################### IF THIS DOESN'T WORK, CONSIDER FOLLOWING ----

## Remove sex chromosome probes ----
### Mark probes on X and Y ----
keepIndex=!seqnames(objectSNPQCed)%in%c("chrX","chrY")
rm(objectSNPQCed);

## Remove failed probes ----
### Combine with failed probes in >20% of samples ----
keepIndex <- keepIndex&(!failedCG02)
rm(failedCG02);

### Remove all probes with detected P-value >0.01 ----
#str(beta) num [1:433375, 1:743]
beta[failed] <- NA
rm(failed);

### Remove marked failed probes ----
betaQC <- beta[which(keepIndex),]
rm(beta, keepIndex);

## OSTENSIBLY 11,818 REMOVED, BUT N REMAINS AT 421557... HUGE

# Prep2: Reformat beta value to match the format suggested above ----
## Load quantile normalized methylation data ---- NOTE: takes approx 1 minute to run
Methy <- as.data.frame(t(betaQC))
## Add subject ID to the methylation data as the last column ----
Methy$Subject <- colnames(betaQC)
rm(betaQC); gc()
## Re-order the column and making the subject ID as the first column ----
Methy <- Methy[,c(ncol(Methy),1:(ncol(Methy)-1))]
## Save probe names ----
MethyName <- colnames(Methy)[-c(1)]

## Save QCed methylation data for later use ----
save(Methy, file = "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Methy.RData")

# Prep3: format methylation covariates, including the first four PCs of beta values and the first 2 PCs of cell type ----
## Load principal components of beta value ----
load("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/fast_svd.rda")

## Generate a variable with first four principal components of beta value ----
PC_Beta <- as.data.frame(ss$v[,1:4])


# NOTE: Error in `$<-.data.frame`(`*tmp*`, Subject, value = c("201004900088_R02C02",  :
 # replacement has 743 rows, data has 792
 # PC_Beta (n =792)
 # colnames(b) = contains all of the 792 Subject IDs
 # Methy$Subject = 743 Subject IDs

 # We have to create a PC_Beta$Subject column based off colnames(b)
 # Remember that ss, and by association, PC_Beta, is...
 # save(beta, file="/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Beta_Quantile.rda") #saves the normalized beta values to new file

    # Generate PCA to control for unknown structure #
    #library(corpcor)
    #ss <- fast.svd(beta,tol=0) # this may take a while
    # Potentially need to take the Methy$Subject and only select those from Beta_Quantile

    load("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Beta_Quantile.rda")
    str(data)

 relevant_subjects <- colnames(b)

 ## To get n = 743....
 ## Note this takes approx... 20 minutes to run
 #Error in fast.svd(beta, tol = 0) : could not find function "fast.svd"
 library(corpcor)
 ss <- fast.svd(beta,tol=0)


 PC_Beta <- as.data.frame(ss$v[,1:4])

 ## Load principal components of beta value ----
load("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/Methy.RData")

## Add subject ID to the component variable ----
PC_Beta$Subject <- Methy$Subject

 write.csv(PC_Beta, "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/PC_Beta.csv")
## SUCCESS

## Load estimated cell type proportion ----
load("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/cellcount.rda")

cellcount <- data

## Generate principal components of cell count ----
tmp <- prcomp(cellcount)

## Generate a variable with subject ID and first 2 PCs of cell count ----
pc_cell <- as.data.frame(tmp$x[,1:2])

## Add row names ----
pc_cell$Subject <- rownames(pc_cell)


write.csv(pc_cell, "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/pc_cell.csv")

## save these

### n.b 3041 obs

## ----------------------------# use imaging script

#### N.b only n = 591 NOTE: check to see if higher n can be achieved

ID_linkage <- data.frame(Subject_meth = dat$Sentrix, Subject = dat$ID)

#Cortex <- read.csv("./CorticalMeasure_ENIGMA.csv",header=T)

covariates_examine2 <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Input/Covariates.csv")
thickness_examine <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Imaging/CorticalMeasuresENIGMA_ThickAvg.csv") # contains ICV
ICV <- thickness_examine %>% select("SubjID", "ICV")


#################################### RUN FROM HERE

# n = 636
RawCov <- merge(ICV, covariates_examine2, by = "SubjID")

names(RawCov)[names(RawCov) == "sex"] <- "Sex"
names(RawCov)[names(RawCov) == "ageyrs"] <- "Age"
names(RawCov)[names(RawCov) == "SubjID"] <- "Subject"

RawCov$Sex <- ifelse(RawCov$sex == 1,1,2)

### Calculate age^2 and add to the last column of Raw_Cov
RawCov$Age_Square <- (RawCov$Age)^2

### Re-name the ID column for easier merge ----
RawCov <- merge(RawCov, ID_linkage, by = "Subject")

## need RO1CO1
#names(RawCov)[1] <- c("Subject")

## Merge Methylation covariates in here
names(PC_Beta)[names(PC_Beta) == "Subject"] <- "Subject_meth"
names(pc_cell)[names(pc_cell) == "Subject"] <- "Subject_meth"


### Merge all covariates into one dataset ----
Cov <- merge(RawCov, PC_Beta, by ="Subject_meth",all=F)
Cov <- merge(Cov,pc_cell,by ="Subject_meth",all=F)

Cov$Sex <- as.numeric(Cov$Sex)

str(Cov)

write.csv(Cov, "/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/CorticalEWAS_covariates.csv")

# n = 541
####################################

## NOTE: THIS APPROACH LEAVES US WITH AN N OF 502 VS 541
## Read and merge covariates ---- # n = 591
#RawCov <- read.csv("./CorticalCovariates_ENIGMA.csv",header=T)
#str(RawCov)
### Assume the IID that is equal to the subject ID from cortical and methylation data ----
#RawCov <- RawCov[,-c(1)]
### Calculate age^2 and add to the last column of Raw_Cov
#RawCov$Age_Square <- (RawCov$Age)^2
### Re-name the ID column for easier merge ----
#RawCov <- merge(RawCov, ID_linkage, by = "Subject")
## need RO1CO1
#names(RawCov)[1] <- c("Subject")
## Merge Methylation covariates in here
#names(PC_Beta)[names(PC_Beta) == "Subject"] <- "Subject_meth"
#names(pc_cell)[names(pc_cell) == "Subject"] <- "Subject_meth"
#PC_Beta2 <- PC_Beta %>% rename(Subject_meth = Subject)
#pc_cell2 <- pc_cell %>% rename(Subject_meth = Subject)
### Merge all covariates into one dataset ----
#Cov <- merge(RawCov, PC_Beta, by ="Subject_meth",all=F)
#Cov <- merge(Cov,pc_cell,by ="Subject_meth",all=F)
## Change Sex to appropriate dummy
#Cov$Sex <- ifelse(Cov$Sex == 'male',1,2)
#str(Cov)
 #Cov$ICV <- as.numeric(Cov$ICV)
 #Cov$Sex <- as.numeric(Cov$Sex)
 #str(Cov)
 ### NOTE: IDEA 2 REMOVE ID and any other misc variables from Cov that could be contributing to the factor issue ----
#Cov <- Cov %>% select(-c("GeneticID", "DiseaseType","PID",  "MID"    ,        "DiseaseType"   ,
 #"AffectionStatus"    ,   "Affectionstatus2" , "ScannerSite" ))
####################################


Cov <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Output/CorticalEWAS_covariates.csv")

### Save column names of covariates for EWAS with cortical measures
Cov <- Cov %>% select(-c("X"))

 Cov$ICV <- as.numeric(Cov$ICV)
 Cov$Sex <- as.numeric(Cov$Sex)
 str(Cov)

Cov2 <- Cov %>% select(-c("Subject", "Subject_meth"))
CovName <- colnames(Cov2)


## Generate files with complete data across all methylation, cortical measures, and covariates with matched individuals ----
# Combine covariates and surface data

CorticalMeasure <- read.csv("/CCACE_Shared/EleanorC/CorticalEWAS/Data/Input/CorticalMeasure_ENIGMA.csv",header=T)
str(CorticalMeasure)

# Data (n = 541)
Data <- merge(Cov,CorticalMeasure, by="Subject", all=F)

# Match methylation data with combined covariates and brain data
Match <- match(Methy$Subject,Data$Subject_meth)

# Generate covariates with matched individual
#Cov <- Data[Match[!is.na(Match)],2:(length(CovName)+1)]
Cov <- Cov %>% select(-c("Subject", "Subject_meth"))

# Generate brain data with matched individual
#CorticalMeasure2 <- Data[Match[!is.na(Match)],-c(1:(length(CovName)+1))]
#CorticalMeasure2 <- Data[Match[!is.na(Match)]]
CorticalMeasure <- Data %>% select(-c("Subject", "Subject_meth","ICV","Sex","Age","Age_Square","V1","V2","V3","V4","PC1","PC2"))

#Generate methylation data with matched individual
Methy <- Methy[!is.na(Match),-c(1)]
str(Methy)
#'data.frame':	541 obs. of  421557 variables:

### Generate a dataframe to save age information for further analysis
AgeInfo <- data.frame(min(Cov$Age),max(Cov$Age),mean(Cov$Age))
names(AgeInfo) <- c("MinAge","MaxAge","MeanAge")

## EWAS with cortical measures ----
# Get the numbers of methylation probes, cortical measures and covariates
Num_Methy <- ncol(Methy)
Num_Cov <- ncol(Cov)
Num_Cortical <- ncol(CorticalMeasure)
####################################


## CHECK STRUCTURE ----
str(Cov) # NOTE: 'data.frame':	499 obs. of  11 variables: (17 obs) vs 541 obs of 10 variables
str(CorticalMeasure) #NOTE: 'data.frame':	499 obs. of  69 variables: (17 obs) vs 541 obs of 18 variables
str(Methy) # NOTE: 'data.frame':	499 obs. of  485120 variables: (17 obs. of  3143 ) vs 541 obs. of  421557 variables:

##################################################################
##Association analysis: the following section should be run by ALL cohorts
##################################################################

######## EWAS with cortical measures: whole sample, controlling for ICV ##########


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
## NOTE: The code here is part of a for loop that's used for performing linear regression analysis between DNA methylation data (Methy) and cortical brain volume data (CorticalMeasure) for each combination of cortical brain volume and DNA methylation variable. This loop is structured as follows:

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
#warning(paste("Warning: This code may take approx 10 minutes to run"))

 # Start measuring time
start_time <- Sys.time()

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

# Stop measuring time
end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Save the execution time to a log file
log_file <- "execution_log.txt"
cat("Linear regression execution time:", elapsed_time, "seconds", file = log_file)

# Print a message to the console
cat("Linear regression completed in", elapsed_time, "seconds. Log saved to", log_file, "\n")


###########################################################################
###########################################################################


