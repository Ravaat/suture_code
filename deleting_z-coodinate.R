# Title: Deleting the z-coordinate

# Name: Heather White

# Date created: 17/06/21

# Last modified: 30/11/21

# License: MIT license


# Converting 3D landmarks to 2D by deleting the z-coordinate after conducting the necessary rotation
# For the coronal suture this is done after the necessary rotation step
# For any sutures in the dorsal plane this is done without any rotation steps

#######################################################################################################

rm(list = ls())

library(tidyverse)

# Load the resampled semilandmarks
load(file = "Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
species_names <- read.csv('Data/Specimen_names.csv', header = F) 

folder <- 'Raw_Data/pts_LMs_and_curves'
# Number of specimens I have
nspecimens<-165
# Read in the file names to a list variable:
ptslist_cor <- list.files(path = folder, pattern = "*.pts", full.names = F)


# Convert the rotated and resampled landmarks to 2D
sagittal_2D <- array(dim = c(500,2,165))
#dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
# Add the rotated landmarks for each specimen to the array:
for(i in 1:length(ptslist_cor))
{
  sagittal_2D[,,i] <- sagittal[,-3,i]
}

# Associate the 2D rotated array with specimen names
dimnames(sagittal_2D)[3]<-species_names
str(sagittal_2D)
glimpse(sagittal_2D)

save(sagittal_2D, file = "Data/2D_sagittal_suture/Sagittal_2D_converted_LMs.Rdata")

write.csv(sagittal_2D, file = "Data/2D_sagittal_suture/Sagittal_2D_converted_LMs.csv")




