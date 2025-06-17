# Title: Allometry

# Name: Heather White

# Date created: 25/08/21

# Last modified: 17/12/21

# License: MIT license



# Linear models looking at the relationship between Procrustes shape variables and CS - allometry

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)


#######################################################################################################


# Load the suture data - mirrored, resampled, slid, Procrusted
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")

info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
CS <- info$CS
Species <- info$Species

# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
sagittal <- unname(sagittal)

# Fit LM with CS - shape is correlated by size (to test allometry)
fit1 <- procD.lm(sagittal~log(CS)) # 1000 permutations automatic; CS is logged
summary(fit1)
# This tells me that shape is significantly correlated with size, but the Rsq value is low (0.16797)
# A low Rsq value tells me that only 16% of shape variation is associated with size
# This is because I have lots of different species

write.csv(fit1$aov.table, file = "Analysis/Suture_morphology/Allometry/allometry_results_sagittal_full_dataset.csv")



#######

# Divide out the data
Bettongia <- sagittal[,,1:10]
Bettongia_CS <- info[1:10,7]
Bradypus <- sagittal[,,11:17]
Bradypus_CS <- info[11:17,7]
Cebus <- sagittal[,,18:23]
Cebus_CS <- info[18:23,7]
Cyclopes <- sagittal[,,24:27]
Cyclopes_CS <- info[24:27,7]
Dasyprocta <- sagittal[,,28:35]
Dasyprocta_CS <- info[28:35,7]
Dasypus <- sagittal[,,36:41]
Dasypus_CS <- info[36:41,7]
Epomops <- sagittal[,,42:46]
Epomops_CS <- info[42:46,7]
Felis <- sagittal[,,46:53]
Felis_CS <- info[46:53,7]
Macroscelides <- sagittal[,,54:59]
Macroscelides_CS <- info[54:59,7]
Manis <- sagittal[,,60:66]
Manis_CS <- info[60:66,7]
Microcebus <- sagittal[,,67:70]
Microcebus_CS <- info[67:70,7]
Monodelphis <- sagittal[,,71:83]
Monodelphis_CS <- info[71:83,7]
Mus <- sagittal[,,84:89]
Mus_CS <- info[84:89,7]
Ornithorhynchus <- sagittal[,,90:93]
Ornithorhynchus_CS <- info[90:93,7]
Phacochoerus <- sagittal[,,94:104]
Phacochoerus_CS <- info[94:104,7]
Phascolarctos <- sagittal[,,105:115]
Phascolarctos_CS <- info[105:115,7]
Rattus <- sagittal[,,116:120]
Rattus_CS <- info[116:120,7]
Setifer <- sagittal[,,121:128]
Setifer_CS <- info[121:128,7]
Setonix <- sagittal[,,129:141]
Setonix_CS <- info[129:141,7]
Sminthopsis <- sagittal[,,142:148]
Sminthopsis_CS <- info[142:148,7]
Talpa <- sagittal[,,149:156]
Talpa_CS <- info[149:156,7]
Trichosaurus <- sagittal[,,157:165]
Trichosaurus_CS <- info[157:165,7]

# Fit the LM for shape by size (to test allometry)
fit_Phacochoerus <- procD.lm(Phacochoerus~log(Phacochoerus_CS)) # 1000 permutations automatic; CS is logged
summary(fit_Phacochoerus)
# Here we see that the allometry is significant 
# We also see that 73% of shape variation is associated with size
# Suggests that across species - taxonomic diversity contributes to skull diversity than size variation

fit_Bettongia <- procD.lm(Bettongia~log(Bettongia_CS)) 
summary(fit_Bettongia)
fit_Bradypus <- procD.lm(Bradypus~log(Bradypus_CS)) 
summary(fit_Bradypus)
fit_Cebus <- procD.lm(Cebus~log(Cebus_CS)) 
summary(fit_Cebus)
fit_Cyclopes <- procD.lm(Cyclopes~log(Cyclopes_CS)) 
summary(fit_Cyclopes)
fit_Dasyprocta <- procD.lm(Dasyprocta~log(Dasyprocta_CS)) 
summary(fit_Dasyprocta)
fit_Dasypus <- procD.lm(Dasypus~log(Dasypus_CS)) 
summary(fit_Dasypus)
fit_Epomops <- procD.lm(Epomops~log(Epomops_CS)) 
summary(fit_Epomops)
fit_Felis <- procD.lm(Felis~log(Felis_CS))
summary(fit_Felis)
fit_Macroscelides <- procD.lm(Macroscelides~log(Macroscelides_CS)) 
summary(fit_Macroscelides)
fit_Manis <- procD.lm(Manis~log(Manis_CS))
summary(fit_Manis)
fit_Microcebus <- procD.lm(Microcebus~log(Microcebus_CS)) 
summary(fit_Microcebus)
fit_Monodelphis <- procD.lm(Monodelphis~log(Monodelphis_CS)) 
summary(fit_Monodelphis)
fit_Mus <- procD.lm(Mus~log(Mus_CS)) 
summary(fit_Mus)
fit_Ornithorhynchus <- procD.lm(Ornithorhynchus~log(Ornithorhynchus_CS)) 
summary(fit_Ornithorhynchus)
fit_Phascolarctos <- procD.lm(Phascolarctos~log(Phascolarctos_CS)) 
summary(fit_Phascolarctos)
fit_Rattus <- procD.lm(Rattus~log(Rattus_CS))
summary(fit_Rattus)
fit_Setifer <- procD.lm(Setifer~log(Setifer_CS)) 
summary(fit_Setifer)
fit_Setonix <- procD.lm(Setonix~log(Setonix_CS)) 
summary(fit_Setonix)
fit_Sminthopsis <- procD.lm(Sminthopsis~log(Sminthopsis_CS))
summary(fit_Sminthopsis)
fit_Talpa <- procD.lm(Talpa~log(Talpa_CS)) 
summary(fit_Talpa)
fit_Trichosaurus <- procD.lm(Trichosaurus~log(Trichosaurus_CS)) 
summary(fit_Trichosaurus)

# Combine the above results into one dataframe
species_results <- rbind(fit_Bettongia$aov.table,fit_Bradypus$aov.table,fit_Cebus$aov.table,fit_Cyclopes$aov.table,
                         fit_Dasyprocta$aov.table, fit_Dasypus$aov.table,fit_Epomops$aov.table,fit_Felis$aov.table,
                         fit_Macroscelides$aov.table,fit_Manis$aov.table,fit_Microcebus$aov.table,fit_Monodelphis$aov.table,
                         fit_Mus$aov.table,fit_Ornithorhynchus$aov.table,fit_Phacochoerus$aov.table,
                         fit_Phascolarctos$aov.table, fit_Rattus$aov.table, fit_Setifer$aov.table, fit_Setonix$aov.table,
                         fit_Sminthopsis$aov.table, fit_Talpa$aov.table,fit_Trichosaurus$aov.table)

write.csv(species_results, file ="Analysis/Suture_morphology/Allometry/allometry_results_each_species_seperately_sagittal.csv")


######

# Fit the model with CS and species
fit2 <- procD.lm(sagittal~log(CS)*Species) # 1000 permutations automatic; CS is logged
summary(fit2)
# This suggests that 15% of variation is due to size and 65% of variation is due to species
write.csv(fit2$aov.table, file = "Analysis/Suture_morphology/Allometry/allometry_results_CS_and_species_sagittal.csv")


##########################################################################################################

# Plotting the allometry results

# Get the raw Procrustes data into a geomorph dataframe
gdf <- geomorph.data.frame(sagittal, CSize = info$CS, Species = info$Species, Specimens = info$Specimen_name)
# name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","CSize","Species", "Specimens")

species <- as.factor(gdf$Species)

# Assigning a colour to each species
col.species = c("paleturquoise1", "pink1", "palevioletred",
                "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
                 "green3", "turquoise1", "chartreuse4", "gold1", 
                "orangered1", "cyan3", "darkorange2", "darkgreen", "darkolivegreen1", "darkorchid4",
                "dodgerblue2", "blue2", "red3", "midnightblue")[gdf$Species]


# Plotting the allometry plot - coloured by species
plotAllometry(fit1, size = gdf$CSize, method = "CAC", pch = 16, col = col.species, cex = 1.5)
# Only care about the plot on the left (log(Size) and CAC) this is a regression between size and shape (CAC = common allometric componet)
# This can go in supplementary
# It tells me that each species has its own allometric trajectory within this plot

# Check assumptions
plot(fit1, type = "regression", predictor = CS)
par(mfrow=c(2, 2), mar=c(5,4,4,2))
plot(fit1, type = "diagnostics")




