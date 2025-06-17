# Title: MANOVA

# Name: Heather White

# Date created: 01/09/21

# Last modified: 17/12/21

# License: MIT license



# MANOVAs and pMANOVAs to see if developmental age category, diet, dev strategy influence skull shape

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(ape)

#######################################################################################################

# Load the data

load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")

info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
Age <- info$Discrete_age
Species <- info$Species
Dev_strategy <- info$Dev_strategy
Diet <- info$Diet

# Get the raw Procrustes data into a geomorph dataframe
gdf <- geomorph.data.frame(sagittal, Age = info$Discrete_age, Species = info$Species, Dev_strat = info$Dev_strategy, Diet = info$Diet)
# name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","Age","Species", "Development_strategy", "Diet")


#######################################################################################################

# MANOVAs on the full dataset - species, age, dev, diet


# SPECIES

# MANOVA for differences in skull shape based on age category
speciesMVA <-procD.lm(Proc_coords ~ Species, f2 = NULL, f3 = NULL, logsz = TRUE, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(speciesMVA)
# This tells me that there is a significant effect of species on shape, the Rsq value is high 

##########

# AGE

# MANOVA for differences in skull shape based on age category
ageMVA <-procD.lm(Proc_coords ~ Age, f2 = NULL, f3 = NULL, logsz = TRUE, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(ageMVA)
# This tells me that there is a significant effect of age on shape, but the Rsq value is low in this categorical grouping

##########

# DEVELOPMENT STRATEGY

# MANOVA for differences in skull shape based on developmental strategy
devMVA <-procD.lm(Proc_coords ~ Development_strategy, f2 = NULL, f3 = NULL, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(devMVA)

##########

# DIET

# MANOVA for differences in skull shape based on diet
# Whole dataset
dietMVA <-procD.lm(Proc_coords ~ Diet, f2 = NULL, f3 = NULL, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(dietMVA)


sink("Analysis/Suture_morphology/MANOVAs/MANOVAs_sagittal_full_dataset.txt")
print("MANOVA for shape vs species - whole dataset")
print(summary(speciesMVA))
print("MANOVA for shape vs age - whole dataset")
print(summary(ageMVA))
print("MANOVA for shape vs developmental strategy - whole dataset")
print(summary(devMVA))
print("MANOVA for shape vs diet - whole dataset")
print(summary(dietMVA))
sink()

#######################################################################################################

# Phylogenetic MANOVA to assess the influence of developmental strategy and diet on shape 
# In adults only

# Read in phylogeny
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
# Load the adult only data
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs_adults.Rdata")
# Load adults only specimen details
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")
# Read in the specimen names to match the tree
tree_names <- read.csv("Data/tree_taxa_names.csv")
tree_names <- tree_names$Taxa_names

adults <- unname(adults_sagittal)
# Check the tree tip names match the shape data
shape.names <-dimnames(adults_sagittal)[[3]]
tree.names <-my_tree$tip.label
setdiff(shape.names,tree.names)
# Set the shape data to have the same names as the phylogeny
dimnames(adults_sagittal)[[3]] <- tree_names

# Add phylogeny to the gdf
gdf <- geomorph.data.frame(adults_sagittal, phy = my_tree, Species = adult_specimen_details$Species,
                           Dev_strategy = adult_specimen_details$Precocial_altricial_spectrum,
                           Diet = adult_specimen_details$Diet)
# Name the data frame 
names(gdf) <-c("Phy","Proc_coords","Species", "Dev_strategy", "Diet")

##############

# Phylogenetic MANOVA for diet
# Also no significance when collapsing the dietary categories into 3 categories - omnivore, herbivore, faunivore
dietPMVA <- procD.pgls(Proc_coords ~ Diet, phy = Phy, data = gdf, iter = 999)                          
summarydiet <- summary(dietPMVA)

# Phylogenetic MANOVA for developmental strategy
devPMVA <- procD.pgls(Proc_coords ~ Dev_strategy, phy = Phy, data = gdf, iter = 999)                          
summarydev <- summary(devPMVA)

# Save results
sink("Analysis/Suture_morphology/MANOVAs/pMANOVAs_sagittal_adults.txt")
print("pMANOVA for shape vs diet - adults")
print(summarydiet)
print("pMANOVA for shape vs developmental strategy - adults")
print(summarydev)
sink()



