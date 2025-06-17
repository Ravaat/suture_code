# Title: Calculating STFT and PSD

# Name: Heather White

# Date created: 05/06/19

# Last modified: 21/06/21

# License: MIT license



# Script to calculate the short-time Fourier transform and the subsequent power spectrum density

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(e1071)
library(psd)
library(geomorph)
library(Momocs)



#######################################################################################################

######
# Load the data - using the 3D converted to 2D semilandmarks 
######

# Load in the 2D suture semilandmark data
load(file = "Data/2D_sagittal_suture/Sagittal_2D_converted_LMs.Rdata")

# Reading in the species names .csv - this is in the same order as the specimens
species_names <- read.csv('Data/Specimen_names.csv', header = F)

folder <- 'Raw_Data/pts_LMs_and_curves'
# Number of specimens I have
nspecimens<-165
# Read in the file names to a list variable:
ptslist <- list.files(path = folder, pattern = "*.pts", full.names = F)



#######################################################################################################

# Scaling the data (to centroid size and line length)
# Procrustes analysis - remove non-shape aspects

# Performing Procrustes'
procrustes <- gpagen(sagittal_2D)
# Saving coordinate data to a new variable once normalised
procrustes_LMs <- procrustes$coords
# Access the raw size of the specimens - centroid size
size <- procrustes$Csize

plotAllSpecimens(procrustes_LMs) # Plot the Procrustes LMs

setwd("/Users/heatherwhite/Documents/PhD/PhD/R_Projects/mammal_suture_development")
#save(procrustes_LMs, file="Analysis/Complexity_results/Sagittal/Procrusted_LMs_2D_sagittal.Rdata")



#######################################################################################################

# STFT Analysis


STFT_array <- array(dim=c(39,64,165)) 
# 64 is the number of columns, determined by the number of coefficients and 39 is the number of rows for each coefficient, final number is the  number of specimens
#FD_summary <- array(79,9)
for(i in 1:165)
{
  run<-stft(procrustes_LMs[,,i], win=min(80,floor(length(procrustes_LMs[,,i])/10)), 
            inc=min(24, floor(length(procrustes_LMs[,,i])/30)), coef=64, wtype="hanning.window")
  STFT_array[,,i]<-run$values
}
# STFT results:
# rows = windows used in STFT; columns = Fourier coefficients (64) for each window

# Checking STFT_array has worked for a random specimen, number 78 in this case
STFT_array[,,16]

# Associating the species name with the correct matrix
# The number 3 is used to say associated the species names to the 3rd part of the array - in this case this is each matrix
dimnames(STFT_array)[3]=species_names

# Saving STFT_array in R
#save(STFT_array, file='Analysis/Complexity_results/Sagittal/STFT_results_sagittal.Rdata')

# Writing STFT_array to a csv
#write.csv(STFT_array, 'Analysis/Complexity_results/Sagittal/STFT_results_sagittal.csv')

##### Rerun above for STFT sagittal 3D


# Converting and saving the array to a 2D matrix
STFT_2D_array <- two.d.array(STFT_array, sep='.')
#save(STFT_2D_array, file='Analysis/Complexity_results/Sagittal/STFT_sagittal_2D_array.Rdata')



#############################################################################################################

# Power Spectrum Density (PSD)


# Calculating the average of the squared STFT coefficients over each frequency/local transforms
# This calculates the Power value for each harmonic 
STFT_average <- array(dim=c(39,1,165)) 
for(i in 1:165)
{
  STFT_average[,,i] <- rowMeans(STFT_array[,,i]^2, na.rm = FALSE, dims = 1)
}

# Viewing the results
STFT_average



# Summing the Power values at each harmonic to get a power value for each suture
PSD_array <- array(dim = c(165,1))
for(i in 1:165)
{
  PSD_array[i,] <- sum(STFT_average[,,i])
}

# Viewing the PSD values for each suture
PSD_array

# Associating the species name with the correct matrix
# The number 3 is used to say associated the species names to the 3rd part of the array - in this case this is each matrix
dimnames(PSD_array)[1]=species_names


# Saving PSD results to a csv
write.csv(PSD_array, "Analysis/Complexity_results/Sagittal/PSD_results_sagittal.csv")




