# Title: Calculating fractal dimension

# Name: Heather White

# Date created: 29/05/19

# Last modified: 30/11/21

# License: MIT license



# Script to calculate the fractal dimension of each specimen

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(fractaldim)
library(data.table)
library(pcaPP)
library(wavelets)
library(geomorph)



#######################################################################################################


######
# Load the data - using the 3D converted to 2D semilandmarks 
######

# Load in the 2D suture semilandmark data
load(file = "Data/2D_interfrontal_suture/Interfrontal_2D_converted_LMs.Rdata")

# Reading in the species names .csv - this is in the same order as the specimens
species_names <- read.csv('Data/Specimen_names.csv', header = F)

folder <- 'Raw_Data/pts_LMs_and_curves'
# Number of specimens I have
nspecimens<-165
# Read in the file names to a list variable:
ptslist <- list.files(path = folder, pattern = "*.pts", full.names = F)



##########################################################################################################

# Scaling using Procrustes'


# Performing Procrustes'
Procrustes<-gpagen(IF_2D)
# Saving coordinate data to a new variable once normalised
Procrusted_LMs<-Procrustes$coords



##########################################################################################################

# Calculating fractal dimension for all specimens - boxcounting method

par(mfrow=c(2,4))
FD_array_box <- array(dim=c(165,1))
for(i in 1:165)
{
  run_box<-fd.estim.boxcount(Procrusted_LMs[,,i], nlags = "all", plot.loglog = TRUE,
                         plot.allpoints = TRUE, main = i)
  FD_array_box[i]<-run_box$fd
}



# To get a summary of the FD_array variable
FD_array_box
# To save the FD_array as an R object
#save(FD_array_box, file="Analysis/FD_box_array.R")


# In FD_array get rownames to save as species names
# Have to use the function dimnames not rownames because it is an array
# The [1] accesses the first position of the array, when you do dim(FD_array) this gives you 79, 1. So using [1] renames the species names
dimnames(FD_array_box)[1]=species_names


# Writing FD summary results for each specimen to a csv
write.csv(FD_array_box, "Analysis/Complexity_results/Interfrontal/FD_box_results_interfrontal.csv")







