# Title: Rotating Coronal Suture

# Name: Heather White

# Date created: 15/06/21

# Last modified: 30/11/21

# License: MIT license



# Rotating the full dataset of coronal sutures using the rotation_matrix_test code - both approaches
# Approach 1 = specimen specific rotation angle to calculate the rotation matrix
# Approach 2 = rotate 45 degrees clockwise using a 45 degree rotation matrix

#######################################################################################################

rm(list = ls())


library(Rvcg)
library(rgl)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(tidyverse)

##############################################################################################

# To rotate the coronal curve by 45 degrees


# Read in the specimen dataset - LMs
folder <- 'Raw_Data/pts_LMs_and_curves'
# Read in the specimen names
species_names <- read.csv("Data/specimen_names.csv", header = F)


# Number of specimens I have
nspecimens<-165
# Read in the file names to a list variable:
ptslist_cor <- list.files(path = folder, pattern = "*.pts", full.names = T)
# Load the resampled semilandmarks
load(file = "Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")


# Creating a 45 degree rotation matrix that rotates clockwise around the x-axis
t2 <- (7*pi)/4
rotatmat <- cbind(c(1, 0, 0), c(0, cos(t2), sin(t2)), c(0, -sin(t2), cos(t2)))

# Test the matrix multiplication on one specimen
coronal[,,13] %*% rotatmat

# Multiply the resampled semilandmarks (n=500) for each specimen by the rotation matrix to rotate coronal suture by 45 degrees clockwise
# Create a blank array to fill
rotated_cor <- array(dim = c(500,3,165))
#dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
# Add the rotated landmarks for each specimen to the array:
for(i in 1:length(ptslist_cor))
{
  rotated_cor[,,i] <- coronal[,,i] %*% rotatmat
}

# Associating specimen names to the array
dimnames(rotated_cor)[3] <- species_names
str(rotated_cor)
glimpse(rotated_cor)


#save(rotated_cor, file = "Data/3D_Resampled_Procrusted_suture_LMs/45deg_rotated_coronal_suture_semilandmarks.Rdata")

#write.csv(rotated_cor, file = "Data/3D_Resampled_Procrusted_suture_LMs/45deg_rotated_coronal_suture_semilandmarks.csv")


#################################################################################################################

# Convert the 45 degree rotated landmarks to 2D - by deleting the z-coordinate
cor_2D <- array(dim = c(500,2,165))
#dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
# Add the rotated landmarks for each specimen to the array:
for(i in 1:length(ptslist_cor))
{
  cor_2D[,,i] <- rotated_cor[,-3,i]
}

# Associate the 2D rotated array with specimen names
dimnames(cor_2D)[3] = species_names
str(cor_2D)
glimpse(cor_2D)

#save(cor_2D, file = "Data/2D_coronal_suture/Coronal_2D_converted_LMs.Rdata")

#write.csv(cor_2D, file = "Data/2D_coronal_suture/Coronal_2D_converted_LMs.csv")



