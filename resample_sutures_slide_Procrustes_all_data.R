# Title: Resampling suture LMs with the skull LMs

# Name: Heather White

# Date created: 16/05/21

# Last modified: 30/11/21

# License: MIT license


# Resampling the suture LMs whilst keeping them together with the skull LM data
# Code adapted from Ryan Felice

#######################################################################################################

rm(list=ls())

library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
#library("devtools") # Don't need to run this line
# devtools::install_github("rnfelice/SURGE")   # Don't need to run this line
#install.packages("remotes")  # Don't need to run this line
library(SURGE)
library(RColorBrewer) # for color palettes
library(magick)

#################################################################################################

# 1) Mirror the LM only data - see mirroring_LMs_using_3_LMs.R code
# This is mirrored only data
# Missing LMs have already been dealt with in the aformentioned code
# This dataset has been mirrored but NOT Procrusted afterwards

load("Data/Mirrored_LMs.Rdata")
# Variable for this is called comb.dataset

# Read in the species names data
species_names <- read.csv('Data/Specimen_names.csv', header = F) 


#################################################################################################

# 2) Resample the suture landmarks

# Import a table defining the suture curves
curve_table <- read_csv('Raw_Data/resampled_curves.csv')

# Locate the pts and ply folders
ptsfolder <- "Raw_Data/pts_LMs_and_curves" 
plyfolder <- "Raw_Data/plys_ascii.nosync/" # Need to rename the plys_ascii folder to end in .nosync so all the plys are downloaded

# Import the pts data (filenames and raw data)
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)
my_curves <- create_curve_info(curve_table, n_fixed = 69)
setwd(ptsfolder)

# Resample the curves based on the read in .csv for the sutures
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE)
subsampled.lm[subsampled.lm == 9999] <- NA

# Save the resampled curves
#save(subsampled.lm, file = "Data/LMs_and_resampled_curves.Rdata")
load("Data/LMs_and_resampled_curves.Rdata")


# Plot the resampled suture LMs onto a specimen to check they look OK
atarfa=ply2mesh(file="Raw_Data/plys_ascii/HW_Manis_tricuspis_91.363_skull_decimated.ply")
shade3d(atarfa,col='white')
spheres3d(subsampled.lm[c(1:69),,65], radius = 0.4, col=4)
spheres3d(subsampled.lm[c(70:569),,65], radius = 0.2, col=3)
spheres3d(subsampled.lm[c(570:1069),,65], radius = 0.2, col=6)
spheres3d(subsampled.lm[c(1070:1569),,65], radius = 0.2, col = 5)

atarfa=ply2mesh(file="Raw_Data/plys_ascii.nosync/Bettongia_penicillata_SAM-24228_80days.ply")
shade3d(atarfa,col='white')
spheres3d(subsampled.lm[c(1:69),,1], radius = 0.2, col=4)
spheres3d(subsampled.lm[c(70:569),,1], radius = 0.1, col=3)
spheres3d(subsampled.lm[c(570:1069),,1], radius = 0.1, col=6)
spheres3d(subsampled.lm[c(1070:1569),,1], radius = 0.1, col = 5)

# Photograph the snapshot if necessary
rgl.snapshot("Figures/Resampled_suture_LMs_Bettongia_46days.png")

#########

# Checking the resampled LMs

# Check to make sure curves look okay on each specimen
checkLM(subsampled.lm,path=plyfolder, pt.size = 0.2,suffix=".ply",render="s", begin = 1)

# Create a list that contains all the missing LMs (9999s)
newpts <- subsampled.lm
# Create missing list 
misslist<-createMissingList(dim(newpts)[3])
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 
# Fix the missing LMs 
# If there are no missing LMs there will be nothing to fix
newpts2<-fixLMtps(newpts)

# Check that pts and ply files match (i.e the names) 
# Should come back with '0' if the names match 
ptslist2<-gsub(pattern="\\.pts$","",ptslist)
plylist <-  dir(plyfolder, pattern='.ply', recursive=F)
plylist2<-gsub(pattern="\\.ply$","",plylist)
setdiff(plylist2,ptslist2)
names(misslist) <- ptslist2

setwd("~/Documents/PhD/PhD/R_Projects/mammal_suture_development")


####################################################################################################

# 3) Sliding the curve semiLMs

####################################

# Run the Slider3d_2 functions first

####################################


# Set working directory to the folder with the plys in
# Need to rename this folder from plys_ascii to plys_ascii.nosync to download ply files
setwd("Raw_Data/plys_ascii.nosync")

# Read the ply files into a list of mesh3d files and name this list with specimen names
temp = list.files(pattern="*.ply")
mymeshes = lapply(temp, ply2mesh)
names(mymeshes) <- temp

# Create a vector with the fixed LM points (1-69)
fix <- c(1:69)
# Create a vector with the curve points (70-1569)
surp <- c(1:nrow(subsampled.lm[,,1]))[-fix]

# Slide the semiLMs
# Make sure the working directory is set to the folder containing the ply files
slided4.all_3iterations <- slider3d_2(newpts2$out, SMvector= fix,
                                      surp=surp, sur.path = "./Raw_Data/plys_ascii.nosync", sur.name = temp, # "." here searches the current working directory
                                      meshlist = mymeshes, ignore = NULL,
                                      sur.type = "ply", tol = 1e-10, deselect = TRUE, inc.check = FALSE,
                                      recursive = TRUE, iterations = 3, initproc = TRUE,
                                      pairedLM = 0, mc.cores = 1, bending=TRUE,
                                      fixRepro = FALSE,stepsize=0.2,
                                      missingList=misslist)
dimnames(slided4.all_3iterations[["dataslide"]])[3]<-dimnames(newpts2$out)[3]

# Save the slid semiLMs to an array
slidedLMs <- slided4.all_3iterations$dataslide

#save(slidedLMs, file = "Data/slided_LMs_3iterations.Rdata")

setwd("~/Documents/PhD/PhD/R_Projects/mammal_suture_development")

# Read in the slid semiLMs
load("Data/slided_LMs_3iterations.Rdata")


####################################################################################################

# 4) Join the arrays

# Joining the two arrays from specimen-specific angle of rotation and 45 degree angle of rotation

# Delete the original landmarks from the subsampled.lm dataset - so these can then be replaced with the mirrored LMs
curves_only <- slidedLMs[-c(1:69),,]  # if sliding isn't done slidedLMs should be replaced by subsampled.lm

# Joining the two arrays
joint<-array(dim=c(1633,3,165))
for(i in 1:165)
{
  joint[,,i] <- rbind(comb.dataset$original[,,i], curves_only[,,i])
}
# Associate the array with species names
dimnames(joint)[3] = species_names


############################################################################################################

# 5) Implement Procrustes

# Create a vector containing the numbers of just the LMs only
LMs <- c(1:133)

# Perform Procrustes, using only the LMs in the vector above
Proc <- procSym(joint, use.lm = LMs, reflect = F)
# Subset out the Procrustes coordinates
Procrusted_LMs <- Proc$rotated
# Associate the array with species names
dimnames(Procrusted_LMs)[3] = species_names

#save(Procrusted_LMs, file = "Data/Procrusted_skull_and_suture_LMs.Rdata")


############################################################################################################

# 6) Deal with ABSENT bones - after Procrustes 

# DO NOT NEED TO DO THIS STAGE HERE!!!!
# As all LMs are being deleted
# Need to do this stage if I am using both LM and curve data though

# Dealing with ABSENT bones after Procrustes to slide LMs back on top of each other

# Create a new array that I can use to manipulate the absent bone LMs
# Want to input the Procrusted LMs here
missing_Proc <- Procrusted_LMs[-c(70:130),,]

# Read in a .csv with the absent bone information - only need to include the first LM for each absent bone
absent<-read.csv("Raw_Data/absent_LMs.csv")
# Read in a .csv with the LM information
LM_table <- read_csv('Raw_Data/LM_list.csv')

# Assign the LMs for the missing bones to variables
lm_jugal <- LM_table$LM[which(LM_table$Bone%in%c("jugal"))]
# or could do this: lm_jugal <- which(LM_table$Bone%in%c("jugal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_ventral_premax <- LM_table$LM[which(LM_table$Bone%in%c("ventral_premax"))]
lm_IP <- LM_table$LM[which(LM_table$Bone%in%c("interparietal"))]

# Run all these variably present lines together
# In the variably present jugal - make LMs 16-19 all the same for specimens missing the jugal
for (i in 1:nrow(absent)){
  if( !is.na(absent$Jugal[i]))
    missing_Proc[lm_jugal,c(1:3),i] <- matrix(Procrusted_LMs[16,c(1:3),i], nrow = length(lm_jugal), ncol=3,byrow=TRUE)
}
# In the variably present ventral premax - make LMs 67-69 all the same for specimens missing the ventral premax
for (i in 1:nrow(absent)){
  if( !is.na(absent$Ventral_premax[i]))
    missing_Proc[lm_ventral_premax,c(1:3),i] <- matrix(Procrusted_LMs[67,c(1:3),i], nrow = length(lm_ventral_premax), ncol=3,byrow=TRUE)
}
# In the variably present interparietal- make LMs 38-40 all the same for specimens missing the interparietal
for (i in 1:nrow(absent)){
  if( !is.na(absent$Interparietal[i]))
    missing_Proc[lm_IP,c(1:3),i] <- matrix(Procrusted_LMs[38,c(1:3),i], nrow = length(lm_IP), ncol=3,byrow=TRUE)
}

# Check all the missing LMs are now the same for the jugal of Manis specimens
missing_Proc[,,7]
missing_Proc[,,10]
missing_Proc[lm_jugal,,7]
# Check for the ventral_premax
missing_Proc[,,1]
missing_Proc[lm_ventral_premax,,1]
# Check for the interparietal
missing_Proc[,,15]
missing_Proc[lm_IP,,15]


# Can't check on a mesh after Procrustes, so do the mesh check before Procrustes

# Save the LMs only
missing_Proc_LMs <- missing_Proc[(1:69),,]
#save(missing_Proc_LMs, file = "New/absent_Procrusted_LMs_slid_skull_only_LMs.Rdata")


#################################################################################################

# 7) Subset out the data to extract the sagittal and coronal suture arrays

# Remove all the fixed LM points
Proc_curves_only <- Procrusted_LMs[-c(1:133),,] 

# Subset out the sagittal suture landmarks
sagittal <- Proc_curves_only[(1:500),,]
# Subset out the coronal suture landmarks
coronal <- Proc_curves_only[(501:1000),,]
# Subset out the interfrontal suture landmarks
IF <- Proc_curves_only[(1001:1500),,]


save(sagittal, file = "Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
save(coronal, file = "Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
save(IF, file = "Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")



