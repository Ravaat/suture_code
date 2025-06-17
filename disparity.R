# Title: Disparity

# Name: Heather White

# Date created: 09/09/21

# Last modified: 20/12/21

# License: MIT license



# Disparity analysis for each discrete age category - for each suture

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(RRPP)
library(RColorBrewer)
library(gridExtra)
library(reshape2)

#######################################################################################################

# Load and setup the data

# Load the suture data - mirrored, resampled, slid, Procrusted
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")

specimen_info <- read.csv(file = "Data/Specimen_info_Ornithorhynchus_removed.csv", header =T, sep = ",")
specimen_info_Ornithorhynchus_included <- read.csv("Data/Specimen_info.csv", header = T, sep = ",")

# Load the suture data - mirrored, resampled, slid, Procrusted -adults only
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs_adults.Rdata")
adult_info <- read.csv(file = "Data/Specimen_info_adults.csv", header = T, sep =",")

take2 <- coronal[,,-90:-93]

# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
take2 <- unname(take2)
# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
coronal <- unname(coronal)

# Run allometry so that I can correct shape data for allometry
# Because 18% of my variance is correlated with size, meaning that a lot of the shape variance will be because of size
CS <- specimen_info$CS
fit1 <- procD.lm(take2~log(CS)) # 1000 permutations automatic; CS is logged
summary(fit1)
CS2 <- specimen_info_Ornithorhynchus_included$CS
fit2 <- procD.lm(coronal~log(CS2))
summary(fit2)
CS_adults <- adult_info$CS
fit_adults <- procD.lm(adults_coronal ~ log(CS_adults))
summary(fit_adults)


# Creating a datafrome with all the information
gdf <- geomorph.data.frame(Proc_coords = take2, Age_cat = specimen_info$CS_binned_group, Clade = specimen_info$Major_clades, 
                           Species = specimen_info$Species, Age_percent_CS = specimen_info$CS_percent_adult, CS = specimen_info$CS)

gdf_Ornithorhynchus_included <- geomorph.data.frame(Proc_coords = coronal, Age_cat = specimen_info_Ornithorhynchus_included$Discrete_age_group, Clade = specimen_info_Ornithorhynchus_included$Major_clades, 
                                                    Species = specimen_info_Ornithorhynchus_included$Species, Age_percent_CS = specimen_info_Ornithorhynchus_included$CS_percent_adult, 
                                                    CS = specimen_info_Ornithorhynchus_included$CS)

gdf_adults <- geomorph.data.frame(Proc_coords = adults_coronal, Dev_strategy = adult_info$Precocial_altricial_spectrum, CS = adult_info$CS)

#######################################################################################################

# Morphological disparity

# Calculate Procrustes variances and distances between groups, with p-value for each pair of groups
# How different are each group shapes?

# RUN 1: Disparity for the entire dataset 
# ~ 1 tells it to use the overall mean
disparity_full_dataset <- morphol.disparity(Proc_coords ~ 1, groups = NULL, iter = 999, data = gdf)

#Results and significance
summary(disparity_full_dataset)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_full_dataset.txt")
print(summary(disparity_full_dataset))
sink() 

########

# RUN 2: Disparity for the enture dataset 
# (same as above), accounting for allometry - shape data corrected for allometry
disparity_full_dataset_allometry <- morphol.disparity(f1 = fit1, groups = NULL, iter = 999, data = gdf)

#Results and significance
summary(disparity_full_dataset_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_full_dataset_allometry.txt")
print(summary(disparity_full_dataset_allometry))
sink() 

########

# RUN 3: Disparity between the different age categories, but not considering the species/clade
disparity_age_cat <- morphol.disparity(Proc_coords ~ 1, groups = ~ Age_cat, iter = 999, data = gdf)

#Results and significance
summary(disparity_age_cat)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_age_cat.txt")
print(summary(disparity_age_cat))
sink() 


########

# RUN 4: Disparity between the different age categories for the entire dataset 
# (same as above), accounting for allometry - shape data is corrected for allometry
disparity_age_cat_allometry <- morphol.disparity(f1 = fit1, groups = ~ Age_cat, iter = 999, data = gdf)

#Results and significance
summary(disparity_age_cat_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_age_cat_allometry.txt")
print(summary(disparity_age_cat_allometry))
sink() 

########

# RUN 5: Disparity between clades, considering entire dataset but not the different age categories
disparity_clade <- morphol.disparity(Proc_coords ~ 1, groups = ~ Clade, iter = 999, data = gdf)

#Results and significance
summary(disparity_clade)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade.txt")
print(summary(disparity_clade))
sink() 


########

# RUN 6: Disparity between clades, considering entire dataset but not the different age categories 
# (same as above), accounting for allometry - shape data is corrected for allometry
disparity_clade_allometry <- morphol.disparity(f1 = fit1, groups = ~ Clade, iter = 999, data = gdf)

#Results and significance
summary(disparity_clade_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade_allometry.txt")
print(summary(disparity_clade_allometry))
sink() 

########

# RUN 7: Disparity between clades, considering entire dataset but not the different age categories
disparity_species <- morphol.disparity(Proc_coords ~ 1, groups = ~ Species, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_species)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_species.txt")
print(summary(disparity_species))
sink() 


########

# RUN 8: Disparity between clades, considering entire dataset but not the different age categories
# (same as above), accounting for allometry - shape data corrected for allometry
disparity_species_allometry <- morphol.disparity(f1 = fit2, groups = ~ Species, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_species_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_species_allometry.txt")
print(summary(disparity_species_allometry))
sink() 


########

# RUN 9: Disparity between age categories, including entire dataset AND clade
disparity_clade_age_cat <- morphol.disparity(Proc_coords ~ 1, groups = ~ Clade*Age_cat, iter = 999, data = gdf)

#Results and significance
summary(disparity_clade_age_cat)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade_age_cat.txt")
print(summary(disparity_clade_age_cat))
sink() 


########

# RUN 10: Disparity between age categories, including entire dataset AND clade
# (same as above), accounting for allometry - shape data corrected for allometry
disparity_clade_age_cat_allometry <- morphol.disparity(f1 = fit1, groups = ~ Clade*Age_cat, iter = 999, data = gdf)

#Results and significance
summary(disparity_clade_age_cat_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade_age_cat_allometry.txt")
print(summary(disparity_clade_age_cat_allometry))
sink() 


########

# RUN 11: Disparity between age categories, including entire dataset AND species
disparity_species_age_cat <- morphol.disparity(Proc_coords ~ 1, groups = ~ Species*Age_cat, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_species_age_cat)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_species_age_cat.txt")
print(summary(disparity_species_age_cat))
sink() 


########

# RUN 12: Disparity between age categories, including entire dataset AND species
# (same as above), accounting for allometry - shape data is corrected for allometry
disparity_species_age_cat_allometry <- morphol.disparity(f1 = fit2, groups = ~ Species*Age_cat, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_species_age_cat_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_species_age_cat_allometry.txt")
print(summary(disparity_species_age_cat_allometry))
sink() 


########

# RUN 13: Disparity between clades, including entire dataset AND species
disparity_clade_species <- morphol.disparity(Proc_coords ~ 1, groups = ~ Clade*Species, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_clade_species)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade_species.txt")
print(summary(disparity_clade_species))
sink() 


########

# RUN 14: Disparity between clades, including entire dataset AND species
# Shape data is corrected for allometry
disparity_clade_species_allometry <- morphol.disparity(f1 = fit2, groups = ~ Clade*Species, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_clade_species_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_clade_species_allometry.txt")
print(summary(disparity_clade_species_allometry))
sink() 


########

# RUN 15: Disparity between the different age categories for the entire dataset 
# (same as above), accounting for allometry
# Ornithorhynchus included
disparity_age_cat_Ornithorhynchus_inc_allometry <- morphol.disparity(f1 = fit2, groups = ~ Age_cat, iter = 999, data = gdf_Ornithorhynchus_included)

#Results and significance
summary(disparity_age_cat_Ornithorhynchus_inc_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_age_cat_Ornithorhynchus_inc_allometry.txt")
print(summary(disparity_age_cat_Ornithorhynchus_inc_allometry))
sink() 


########

# RUN 16: Disparity between the different developmental strategies for the adults only
# (same as above), accounting for allometry
# Ornithorhynchus included
disparity_dev_strategy_allometry <- morphol.disparity(f1 = fit_adults, groups = ~ Dev_strategy, iter = 999, data = gdf_adults)

#Results and significance
summary(disparity_dev_strategy_allometry)

sink("Analysis/Suture_morphology/Disparity/Sagittal/disparity_dev_strategy_allometry.txt")
print(summary(disparity_dev_strategy_allometry))
sink() 




################################################################################################################

# Histogram plots for Procustes variance of genera at each stage

# Setup the data to plot histogram - clade vs age cat - accounting for allometry

# Save variances as object
PV_clade_age_cat_allometry <- data.frame(PV = disparity_clade_age_cat_allometry[["Procrustes.var"]])

PV_clade_age_cat_allometry # Look at order of genera and stages to create vectors
#List categories and genera
age_cat_list <- levels(gdf$Age_cat)
clade_list <- levels(gdf$Clade)

# Make columns with correct order and numuber of age categories and clades
PV_age_cat_allometry <- rep(age_cat_list, times = length(clade_list)) # number of times = number of clade - want to make a list with the correct number of comparisons
PV_clade_allometry <- rep(clade_list, each = length(age_cat_list)) # number of repeats = number of age categories

#Add labels and other attributes to tibble as columns
PV_clade_age_cat_allometry <- PV_clade_age_cat_allometry %>% 
  mutate(Age_cat = PV_age_cat_allometry, Clade = PV_clade_allometry)
glimpse(PV_clade_age_cat_allometry)


#########

# Setup the data to plot histogram - species vs age cat

# Save variances as object
PV_species_age_cat_allometry <- data.frame(PV = disparity_species_age_cat_allometry[["Procrustes.var"]])
#write.csv(PV_species_age_cat_allometry, file = "Data/Disparity_PV_species_age_cat_sagittal.csv")
PV_species_age_cat_allometry <- read.csv("Data/Disparity_PV_species_age_cat_coronal.csv")

# Order the data by clade and convert to factor for plotting and plot legend
PV_species_age_cat_allometry$Clade <- as.factor(PV_species_age_cat_allometry$Clade)
PV_species_age_cat_allometry <- PV_species_age_cat_allometry %>% arrange(Clade)
PV_species_age_cat_allometry


########

# Setup the data to plot histogram - species only - accounting for allometry

# Save variances as object
PV_species_allometry <- data.frame(PV = disparity_species_allometry[["Procrustes.var"]])
PV_species_allometry

species <- levels(gdf_Ornithorhynchus_included$Species)
PV_species_allometry <- PV_species_allometry %>% 
  mutate(Species = species)

# Setup the data to plot histogram - developmental strategy only - accounting for allometry
PV_dev_strategy_allometry <- data.frame(PV = disparity_dev_strategy_allometry[["Procrustes.var"]])
PV_dev_strategy_allometry
dev <- levels(gdf_adults$Dev_strategy)
PV_dev_strategy_allometry <- PV_dev_strategy_allometry %>% mutate(Dev_strategy = dev)


# Setup the data to plot histogram for age category only
PV_age_cat_allometry2 <- data.frame(PV = disparity_age_cat_allometry[["Procrustes.var"]])
PV_age_cat_allometry2
age <- levels(gdf_Ornithorhynchus_included$Age_cat)
PV_age_cat_allometry2 <- PV_age_cat_allometry2 %>% mutate(Age = age)


##########################################################################################################

# Plotting the histogram

# Create a colour palette for later use for clades
my_palette_clade = c("mediumpurple3", # Afrotheria
                     "#A6D854", # Euarchontoglires
                     "sandybrown", # Laurasiatheria
                     "cornflowerblue", # Marsupialia
                     "palevioletred") # Xenarthra
# View the colours of the colour palette
image(1:5,1, as.matrix(1:5), col = my_palette_clade,xlab="Oranges (sequential)", ylab = "", yaxt = "n")

# Create a colour palette for later use for ages
my_palette_age = c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF")
# View the colours of the colour palette
image(1:5,1, as.matrix(1:4), col = my_palette_age,xlab="Oranges (sequential)", ylab = "", yaxt = "n")

# Create a colour palete for later use for species
# Create a colour palette
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:28, 1, as.matrix(1:28), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")
# Order the species to match the colours above
PV_species_age_cat_allometry$Species <- factor(PV_species_age_cat_allometry$Species,
                              levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                         "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                         "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                         "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                         "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                         "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                         "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                         "Dasypus novemcinctus"),
                              ordered = TRUE)


########

# Histogram plot by age category and clade
PV_clade_age_cat_ggplot <- ggplot(PV_clade_age_cat_allometry, aes(Age_cat, PV, colour = Clade, fill = Clade)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_color_manual(name = "Clade", labels = c("Afrotheria", "Euarchontoglires", "Laurasiatheria", "Marsupialia", "Xenarthra"),
                     values = my_palette_clade, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  xlab("Growth stage")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by clade and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PV_clade_age_cat_ggplot

# Histogram plot by clade and age category
PV_age_cat_clade_ggplot <- ggplot(PV_clade_age_cat_allometry, aes(Clade, PV, colour = Age_cat, fill = Age_cat)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels = c("Afrotheria", "Euarchontoglires", "Laurasiatheria", "Marsupialia", "Xenarthra"))+
  scale_color_manual(name = "Growth Stage", labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                                                     "Group3" = "Sub-adult", "Group4" = "Adult"),
                     values = my_palette_clade, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  xlab("Clade")+
  ylab("PV (Procrustes Variance)")+
  ggtitle ("Disparity by age category and clade")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PV_age_cat_clade_ggplot

#########

# Histogram plot by species and age category
PV_species_age_cat_ggplot <- ggplot(PV_species_age_cat_allometry, aes(Age_cat, PV, colour = Species, fill = Species)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_color_manual(values = mypalette_species, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  xlab("Growth stage")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PV_species_age_cat_ggplot

# Histogram plot by species and age category
PV_age_cat_species_ggplot <- ggplot(PV_species_age_cat_allometry, aes(Species, PV, colour = Age_cat, fill = Age_cat)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                         "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                         "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                         "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                         "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                         "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                         "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                         "Dasypus novemcinctus"))+
  scale_color_manual(name = "Species",labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                                               "Group3" = "Sub-adult", "Group4" = "Adult"),
                     values = my_palette_age, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  xlab("Species")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PV_age_cat_species_ggplot


#########

# Histogram plotted for species only - accounting for allometry
PV_species_only_ggplot <- ggplot(PV_species_allometry, aes(Species, PV)) +
  geom_col(position = position_dodge2(padding = 0.2), fill = "grey")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  theme_classic(base_size = 12)+
  xlab("Species")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species only")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12))
PV_species_only_ggplot


#########

# Histogram plotted for age only - accounting for allometry
PV_age_only_ggplot <- ggplot(PV_age_cat_allometry2, aes(Age, PV)) +
  geom_col(position = position_dodge2(padding = 0.2), fill = "grey")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  theme_classic(base_size = 12)+
  xlab("Age")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by age category only")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12))
PV_age_only_ggplot


#########

# Histogram plotted for developmental strategy only - accounting for allometry
PV_dev_strategy_ggplot <- ggplot(PV_dev_strategy_allometry, aes(Dev_strategy, PV)) +
  geom_col(position = position_dodge2(padding = 0.2), fill = "grey")+
  scale_y_continuous(expand = c(0,0), limits= c(0,0.5))+
  theme_classic(base_size = 12)+
  xlab("Developmental strategy")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by developmental strategy only")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12))
PV_dev_strategy_ggplot


##########################################################################################################

# Heatmaps for disparity significant differences

# FUNCTION: Get lower triangle of the correlation matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}
# FUNCRION: Get upper triangle of the correlation matrix
get_upper_tri <- function(x){
  x[lower.tri(x)]<- NA
  return(x)
}
# FUNCTION: Reorder table
reorder_corr_table <- function(x){
  # Use correlation between variables as distance
  dd <- as.dist((1-x)/2)
  hc <- hclust(dd)
  x <-x[hc$order, hc$order]
}

# Create palette for heatmap plot
mypalette_seq <- brewer.pal(9,"Oranges")
image(1:9,1, as.matrix(1:9), col = mypalette_seq,xlab="Oranges (sequential)",
      ylab = "", yaxt = "n")

# Raw data
# Save p-values as object
disp_corr <- disparity_species$PV.dist
disp_pvals <- disparity_species$PV.dist.Pval

# Save row and col names as variables to change string - colnames = rownames for both
vars <- rownames(disp_corr)

# Set correct row and col names for both
rownames(disp_corr) <- vars
rownames(disp_pvals) <- vars
colnames(disp_corr) <- vars
colnames(disp_pvals) <- vars

# Get upper triangles only - half matrix, eliminates redundant info
disp_corr_upper_tri <- get_upper_tri(disp_corr)
disp_pvals_upper_tri <- get_upper_tri(disp_pvals)

# Melt to make table in the format needed for heatmap
disp_corr_melt <- melt(disp_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_pvals_melt <- melt(disp_pvals_upper_tri, value.name = "p", na.rm = TRUE)

# Add column to main table
disp_pvals_melt$corr <- disp_corr_melt$corr

# Create columns where only significant values are shown
disp_pvals_melt <- disp_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                              p_if_sig = ifelse(sig_p, p, NA),
                                              corr_if_sig = ifelse(sig_p, corr, NA))
glimpse(disp_pvals_melt)

# Plot the heatmap
disparity_species_heatmap <- ggplot(data = disp_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq[9], high = mypalette_seq[2], mid = mypalette_seq[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Disparity - heatmap for species")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 90, size = 14, vjust = 0, hjust = 1),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5))
disparity_species_heatmap


