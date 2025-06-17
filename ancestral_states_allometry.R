# Title: Ancestral states allometry

# Name: Heather White

# Date created: 01/10/21

# Last modified: 20/12/21

# License: MIT license



# Ancestral state reconstruction from species allometry - for each suture seperately
# Adapted from Morris et al. 2019

#######################################################################################################

rm(list = ls())

library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(RColorBrewer) 
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(png)
library(gridExtra)
library(phytools)
library(evomap)
library(abind)
library(dplyr)

#require(devtools)
#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

####################################################################################################

# 1) Read in the data

# Import trees in Nexus format - branch lengths needed
tree <- "Data/my_mammal_tree.nexus"  
# Read the tree for analysis
tree_species <- read.nexus(tree) 
plot(tree_species)
# Check names 
summary(tree_species)


# Load the PCA results data - PC scores
load(file ="Data/PC_scores_sagittal_all_specimens.Rdata")
load(file ="Data/PC_scores_coronal_all_specimens.Rdata")
load(file ="Data/PC_scores_interfrontal_all_specimens.Rdata")

# Load the suture data - mirrored, resampled, slid, Procrusted
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")

# Load the specimen info
info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
CS <- info$CS
logCS <- log(CS)
Phylogeny_species_list <- read.csv("Data/tree_taxa_names.csv")


#######################################################################################################

# 2) Manipulate the data into the correct format for analysis

# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
coronal <- unname(coronal)

# Combine data into dataframe
gdf <- geomorph.data.frame(coronal, CSize = CS, logCSize = logCS, Species = info$Species, Specimens = info$Specimen_name)
# Name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","CSize", "logCSize", "Species", "Specimens")

# Pull out the species names info
species_list <- factor(gdf$Species,
       levels = c("Bettongia penicillata","Bradypus tridactylus","Sapajus apella", 
                  "Cyclopes didactylus", "Dasyprocta leporina", "Dasypus novemcinctus",
                  "Epomops franqueti", "Felis catus", "Macroscelides proboscideus",
                  "Phataginus tricuspis", "Microcebus murinus", "Monodelphis domestica", 
                  "Mus musculus", "Ornithorhynchus anatinus", "Phacochoerus aethiopicus",
                  "Phascolarctos cincereus", "Rattus rattus", "Setifer setosus",
                  "Setonix brachyurus", "Sminthopsis macroura", "Talpa europaea",  
                  "Trichosaurus vulpecula"),
       ordered = TRUE)
species_list <- levels(species_list)


# Seperate the dataset by species
# Create a list containing a list for each species in the dataset
species_rows <- list()
for (i in 1:length(species_list)){
  species_rows[[i]] <- which(PCAresults_all_sagittal$Species == species_list[i])
}

# Seperate out the shape (coords) and size data by species
# A list for each species containing this information is produced
coords_list <- list()
logCsize_list <- list()
for (i in 1:length(species_list)){
  coords_list[[i]] <- gdf$Proc_coords[,,species_rows[[i]]]
  logCsize_list[[i]] <- gdf$logCSize[species_rows[[i]]]
}

# Assocaite species names with the subsetted dataset above
# This has to be the same names as the phylogeny - must match
names(coords_list) = Phylogeny_species_list$Taxa_names
names(logCsize_list) = Phylogeny_species_list$Taxa_names


#####################################################################################################

# 3) Analysis - allometry/regression

# Calculate allometric regression for each species seperately 
allometry_species <- list()
for (i in 1:22){
  allometry_species[[i]] <- procD.lm(coords_list[[i]] ~ logCsize_list[[i]], iter = 999)
}
allometry_species


# View the allometry results - seperate into a list for each species
summary_allometry_species <- list()
for (i in 1:22){
  summary_allometry_species[[i]] <- summary(allometry_species[[i]])
 }

# Assign species names to the allometry results
# These have to match the phylogeny
names(summary_allometry_species) <- Phylogeny_species_list$Taxa_names
summary_allometry_species

# Save results to file
sink("Analysis/Suture_morphology/Ancestral_trajectories/Sagittal/allometry_results_each_species_seperately_sagittal.txt", append = F)
summary_allometry_species
sink() 

######
# These results above are the same as the results from the allometry.R code
######


# Plot to obtain regscores for allometry
allometry_species_plot <- list()
for (i in 1:22){
  allometry_species_plot[[i]] <- plot(allometry_species[[i]], type = "regression", predictor = logCsize_list[[i]], reg.type = "RegScore",
                                      main = "Shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
}


#########################################################################################################

# 4) Pool the allometry results to be in a single table for plotting - obtain regression score information

# Pooling the regression scores (y) and logCS (x) from the allometry plots above into a dataframe 
# Calculated for each species individually
# Allometry done for each species seperately - this matches the regression lines plotted
x <- list()
for (i in 1:22){
x[[i]] <- c(allometry_species_plot[[i]][["plot.args"]][["y"]])
allometry_regscores <- unlist(x)
allometry_plot_data <- cbind(logCsize = logCS, RegScores = allometry_regscores)
}

# Convert data frame to tibble - to add information
allometry_plot_data <- as_tibble(allometry_plot_data)
allometry_plot_data <- allometry_plot_data %>% mutate(Specimen_names = info$Specimen_name, Species = info$Species, 
                            Clade = info$Major_clades)
glimpse(allometry_plot_data)


###################################################################################################################

# 5) Analysis - linear model to get slope/intercept info from allometry analysis

# Linear model for regression score (y) against logCS (x) - to get the line of best fit - slope and intercept
allometry_species_regline <- list()
for (i in 1:22){
  allometry_species_regline[[i]] <- lm(allometry_species_plot[[i]][["plot.args"]][["y"]] ~ allometry_species_plot[[i]][["plot.args"]][["x"]])
}


# Save the slope and intercept values from the linear model
SlopesList <- matrix()
InterceptsList <- matrix()
for (i in 1:22){
  SlopesList[[i]] <- allometry_species_regline[[i]][["coefficients"]][2]
  InterceptsList[[i]] <- allometry_species_regline[[i]][["coefficients"]][1]
}

# Assign species names to the allometry results
# These have to match the phylogeny
names(SlopesList) <- Phylogeny_species_list$Taxa_names
names(InterceptsList) <- Phylogeny_species_list$Taxa_names


########################################################################################################

# 6) Analysis - calculate slopes for ancestral states

# Calculate ancestral states with phytools - slope and intercept for each component
slope_anc.ML <- anc.ML(tree_species, SlopesList, CI = F)
int_anc.ML <- anc.ML(tree_species, InterceptsList, CI = F)

# Plot the phylogeny with node labels to check where each is
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5) # add ancestral node information to the tree
# Save the ancestral node numbers
nodes <- length(slope_anc.ML$ace)
# Create col names for ancestral data frame
anc_col_names <- c("Slope","Intercept")

# Create dataframe with ancestral state values, one for each component
anc_values_traj <- matrix(nrow = nodes, ncol = 2, dimnames = list(names(slope_anc.ML$ace),anc_col_names))
anc_values_traj[,"Slope"] <- slope_anc.ML$ace
anc_values_traj[,"Intercept"] <- int_anc.ML$ace
anc_values_traj


#######################################################################################################

# 7) Prepare the data for plotting

# Create dataframes for slope and intercept values
# To include both extant species and calculated ancestral states
slope <- data.frame(Slope = c(SlopesList,anc_values_traj[,1]))
intercept <- data.frame(Intercept = c(InterceptsList,anc_values_traj[,2]))
# Combine the slope and intercept data into one dataframe for plotting
values_traj <- cbind(intercept,slope)
values_traj

# Name the ancestral clade nodes
rownames(values_traj)[rownames(values_traj) == "23"] <- "Ancestral mammal"
rownames(values_traj)[rownames(values_traj) == "24"] <- "Ancestral therian mammal"
rownames(values_traj)[rownames(values_traj) == "25"] <- "Ancestral placental mammal"
rownames(values_traj)[rownames(values_traj) == "27"] <- "Ancestral Laurasiatheria"
rownames(values_traj)[rownames(values_traj) == "31"] <- "Ancestral Euarchontoglires"
rownames(values_traj)[rownames(values_traj) == "36"] <- "Ancestral Afrotheria"
rownames(values_traj)[rownames(values_traj) == "37"] <- "Ancestral Xenarthra"
rownames(values_traj)[rownames(values_traj) == "39"] <- "Ancestral Marsupialia"


# Make tibble for species and ancestral nodes wanting to plot
values_traj <- values_traj[c(1:25,27,31,36,37,39),]

# Save results to file
sink("Analysis/Suture_morphology/Ancestral_trajectories/Sagittal/ancestral_states_slopes_sagittal.txt", append = F)
values_traj
sink() 

write.csv(values_traj, file = "Analysis/Suture_morphology/Ancestral_trajectories/Sagittal/ancestral_states_slopes_sagittal.csv")

# Collate all the clade info so this can be appended to the tibble
clade_info <- c("Marsupialia", "Xenarthra", "Euarchontoglires", "Xenarthra", "Euarchontoglires",
                "Xenarthra", "Laurasiatheria", "Laurasiatheria", "Afrotheria", "Laurasiatheria",
                "Euarchontoglires", "Marsupialia", "Euarchontoglires", "Monotremata",
                "Laurasiatheria", "Marsupialia", "Euarchontoglires", "Afrotheria", "Marsupialia",
                "Marsupialia", "Laurasiatheria", "Marsupialia", "Ancestral", 
                "Ancestral", "Ancestral", "Ancestral", "Ancestral",
                "Ancestral", "Ancestral", "Ancestral")

# Set the different lines to have different alpha values
alpha_vector = rep(0.5, nrow(values_traj))
alpha_vector[c(3:10)] = 1
values_traj$alpha = alpha_vector

# Collate all the species info so this can be appended to the tibble
species_info <- c(species_list, "Ancestral mammal", "Ancestral therian mammal", "Ancestral placental mammal", 
                "Ancestral Laurasiatheria", "Ancestral Euarchontoglires", 
                "Ancestral Afrotheria", "Ancestral Xenarthra", "Ancestral Marsupialia")

# Collate all the extant/ancestral info so this can be appended to the tibble
extant_info <- c("Extant","Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant",
                 "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant",
                 "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Ancestral", "Ancestral", 
                 "Ancestral", "Ancestral", "Ancestral", "Ancestral", "Ancestral", "Ancestral")

# Append the specimen info to the tibble
values_traj <- values_traj %>% as_tibble() %>% mutate(Species_Phylogeny = rownames(values_traj), Species = species_info, Clade = clade_info, Extant = extant_info)
values_traj

# Order the data by clade and convert to factor for plotting and plot legend
values_traj$Clade <- as.factor(values_traj$Clade)
values_traj <- values_traj %>% arrange(Clade)
values_traj

###### 
# Data for the plot:
# values_traj = slope and intercept info for plotting the species lines
# allometry_plot_data = regression and logCS data for setting up the plot
######

# Check groups and variables
glimpse(allometry_plot_data)

# Create a colour palette
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "gray10", "gray35", "gray55", "darkorange2", "green3", "mediumpurple", "palevioletred", "dodgerblue2", # Ancestral
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:30, 1, as.matrix(1:30), col = mypalette_species, xlab = "Macroscelides proboscideus, Setifer setosus, Ancestral mammal, Ancestral Laurasiatheria, Ancestral Euarchontoglires, Ancestral Afrotheria, Ancestral Xenarthra, Ancestral Marsupialia, Cebus apella, Dasyprocta leporina, Microcebus murinus, Mus musculus, Rattus rattus, Epomops franqueti, Felis catus, Manis tricuspis, Phacochoerus aethiopicus, Talpa europaea, Bettongia penicillata, Monodelphis domestica, Phascolarctos cincereus, Setonix brachyurus, Sminthopsis macroura, Trichosaurus vulpecula, Ornithorhynchus anatinus, Bradypus tridactylus, Cyclopes didactylus, Dasypus novemcinctus",
      ylab = "", yaxt = "n")

mypalette_species_clade <- c("mediumpurple", "mediumpurple", # Afrotheria
                             "gray10", "gray35", "gray55", "darkorange2", "green3", "mediumpurple", "palevioletred", "dodgerblue2", # Ancestral
                             "green3", "green3", "green3", "green3", "green3", # Euarchontoglires
                             "darkorange2", "darkorange2", "darkorange2", "darkorange2", "darkorange2", # Laurasiatheria
                             "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", # Marsupialia
                             "gold1", # Monotremata
                             "palevioletred", "palevioletred", "palevioletred") # Xenarthra
image(1:30, 1, as.matrix(1:30), col = mypalette_species_clade, xlab = "Macroscelides proboscideus, Setifer setosus, Ancestral mammal, Ancestral Laurasiatheria, Ancestral Euarchontoglires, Ancestral Afrotheria, Ancestral Xenarthra, Ancestral Marsupialia, Cebus apella, Dasyprocta leporina, Microcebus murinus, Mus musculus, Rattus rattus, Epomops franqueti, Felis catus, Manis tricuspis, Phacochoerus aethiopicus, Talpa europaea, Bettongia penicillata, Monodelphis domestica, Phascolarctos cincereus, Setonix brachyurus, Sminthopsis macroura, Trichosaurus vulpecula, Ornithorhynchus anatinus, Bradypus tridactylus, Cyclopes didactylus, Dasypus novemcinctus",
      ylab = "", yaxt = "n")


values_traj$Species <- factor(values_traj$Species,
                             levels = c("Macroscelides proboscideus", "Setifer setosus", "Ancestral mammal",
                              "Ancestral therian mammal", "Ancestral placental mammal",
                              "Ancestral Laurasiatheria", "Ancestral Euarchontoglires", "Ancestral Afrotheria",      
                              "Ancestral Xenarthra", "Ancestral Marsupialia", "Sapajus apella",            
                              "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                              "Rattus rattus", "Epomops franqueti", "Felis catus",       
                              "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                              "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                              "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                              "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                              "Dasypus novemcinctus"),
                               ordered = TRUE)



############################################################################################################

# 8) Plot the species and ancestral allometry trajectories

# Plot allometric trajectories with abline for each group and ancestral node
allometry_anc_nodes_ggplot <- ggplot(allometry_plot_data, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant, alpha = Extant), size = 1)+
  scale_colour_manual(values = mypalette_species_clade)+    
  scale_alpha_discrete(range = c(1, 0.2))+
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Species allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-3.5,1.5)+
  xlim(2,8)
allometry_anc_nodes_ggplot


##########################################################################################################

# 9) Plot for each clade seperately

# Seperate out the various clades in the two datasets (values_traj and allometry_plot_data) 
# to make a clade-specific and corresponding ancestral state plot

values_traj_Afrotheria <- values_traj[c(1:3,8),]
values_traj_Xenarthra <- values_traj[c(3,9,28:30),] 
values_traj_Euarchontoglires <- values_traj[c(3,7,11:15),]
values_traj_Laurasiatheria <- values_traj[c(3,6,16:20),]
values_traj_Marsupialia <- values_traj[c(3,10,21:27),] # This includes the monotreme
values_traj_anc <- values_traj[3:10,]

allometry_plot_data_Afrotheria <- allometry_plot_data %>% filter(Clade == "Afrotheria")
allometry_plot_data_Xenarthra <- allometry_plot_data %>% filter(Clade == "Xenarthra")
allometry_plot_data_Laurasiatheria <- allometry_plot_data %>% filter(Clade == "Laurasiatheria")
allometry_plot_data_Euarchontoglires <- allometry_plot_data %>% filter(Clade == "Euarchontoglires")
allometry_plot_data_Marsupialia <- allometry_plot_data %>% filter(Clade == "Marsupialia" | Clade == "Monotremata") # This includes the monotreme

# Make species name a factor to ensure it is associated with the correct info on the plot
values_traj_Afrotheria$Species <- factor(values_traj_Afrotheria$Species,
                              levels = c("Macroscelides proboscideus", "Setifer setosus", "Ancestral mammal",          
                                          "Ancestral Afrotheria"),
                              ordered = TRUE)
values_traj_Xenarthra$Species <- factor(values_traj_Xenarthra$Species,
                              levels = c("Ancestral mammal",
                                         "Ancestral Xenarthra", "Bradypus tridactylus", "Cyclopes didactylus",      
                                         "Dasypus novemcinctus"),
                              ordered = TRUE)
values_traj_Laurasiatheria$Species <- factor(values_traj_Laurasiatheria$Species,
                              levels = c("Ancestral mammal",          
                                         "Ancestral Laurasiatheria", "Epomops franqueti", "Felis catus",       
                                         "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea"),
                              ordered = TRUE)
values_traj_Euarchontoglires$Species <- factor(values_traj_Euarchontoglires$Species,
                              levels = c("Ancestral mammal", "Ancestral Euarchontoglires", "Sapajus apella",            
                                         "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                         "Rattus rattus"),
                              ordered = TRUE)
values_traj_Marsupialia$Species <- factor(values_traj_Marsupialia$Species,
                              levels = c("Ancestral mammal", "Ancestral Marsupialia",           
                                         "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                         "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                         "Ornithorhynchus anatinus"),
                              ordered = TRUE)
values_traj_anc$Species <- factor(values_traj_anc$Species, levels = c("Ancestral mammal", "Ancestral therian mammal",
                                  "Ancestral placental mammal", "Ancestral Laurasiatheria", "Ancestral Euarchontoglires",
                                  "Ancestral Afrotheria", "Ancestral Xenarthra", "Ancestral Marsupialia"))

##########

# Colour palettes for each clade

mypalette_xenarthra <- c("gray10", "violetred", #  Ancestral
                       "palevioletred3", "violetred", "hotpink4") # Xenarthra
image(1:5, 1, as.matrix(1:5), col = mypalette_xenarthra, xlab = "Species",
      ylab = "", yaxt = "n")

mypalette_afrotheria <- c("mediumpurple", "mediumpurple4",  # Afrotheria
                          "gray10", "mediumpurple") # Ancestral
image(1:4, 1, as.matrix(1:4), col = mypalette_afrotheria, xlab = "Species",
      ylab = "", yaxt = "n")

mypalette_laurasiatheria <- c("gray10", "coral",  # Ancestral
                       "sandybrown", "coral", "darkorange2", "orangered", "red3") # Laurasiatheria
image(1:7, 1, as.matrix(1:7), col = mypalette_laurasiatheria, xlab = "Species",
      ylab = "", yaxt = "n")

mypalette_euarchontoglires <- c("gray10", "limegreen",  # Ancestral
                       "darkolivegreen1", "olivedrab3", "limegreen", "forestgreen", "darkgreen") # Euarchontoglires
image(1:7, 1, as.matrix(1:7), col = mypalette_euarchontoglires, xlab = "Species",
      ylab = "", yaxt = "n")

mypalette_marsupialia <- c("gray10", "steelblue2", # Ancestral
                       "paleturquoise2", "skyblue", "steelblue2", "royalblue1", "dodgerblue4", "navyblue", # Marsupialia
                       "gold1") # Monotremata
image(1:9, 1, as.matrix(1:9), col = mypalette_marsupialia, xlab = "Species",
      ylab = "", yaxt = "n")

mypalette_anc <- c("gray10", "gray35", "gray55", "coral", "limegreen", "mediumpurple", "violetred", "steelblue2")
image(1:8, 1, as.matrix(1:8), col = mypalette_anc, xlab = "Species",
      ylab = "", yaxt = "n")

##########

# Plotting ancestral allometry with each clade

# Ancestral
allometry_ancestral_ggplot <- ggplot(allometry_plot_data, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_anc, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_anc)+  
  scale_linetype_manual(name = NULL, labels = c("Ancestral"), values = c(1))+
  theme_classic(base_size = 12)+
  ggtitle("Ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2,8)
allometry_ancestral_ggplot


# Afrotheria
allometry_anc_afrotheria <- ggplot(allometry_plot_data_Afrotheria, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_Afrotheria, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_afrotheria)+           
  scale_linetype_manual(name = NULL, labels = c("Extant", "Ancestral"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Afrotheria allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2,8)
allometry_anc_afrotheria


# Xenarthra
allometry_anc_xenarthra <- ggplot(allometry_plot_data_Xenarthra, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_Xenarthra, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_xenarthra)+           
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Xenarthra allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2,8)
allometry_anc_xenarthra

# Laurasiatheria
allometry_anc_laurasiatheria <- ggplot(allometry_plot_data_Laurasiatheria, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_Laurasiatheria, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_laurasiatheria)+           
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Laurasiatheria allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2,8)
allometry_anc_laurasiatheria

# Euarchontoglires
allometry_anc_euarchontoglires <- ggplot(allometry_plot_data_Euarchontoglires, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_Euarchontoglires, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_euarchontoglires)+           
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Euarchontoglires allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2, 8)
allometry_anc_euarchontoglires

# Marsupialia
allometry_anc_marsupilia <- ggplot(allometry_plot_data_Marsupialia, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj_Marsupialia, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant), size = 1)+
  scale_colour_manual(values = mypalette_marsupialia)+           
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Marsupialia and Monotremata allometry and ancestral states allometery - sagittal suture")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-2,1.5)+
  xlim(2,8)
allometry_anc_marsupilia


##########################################################################################################

# Comparing the ancestral state to the data which was used to produce the ancestral state
# is circular and therefore logically flawed and should not be done.

# First I take the coordinate and CS data and perform a regression (allometry) for the species
# Then I plot the line for this regression using slope and intercept
# Then I infer the slope and intercept for the ancestral states from the species lines above
# To do a pairwise comparison of slopes I need to have regression data for each line
# Therefore I need to infer the points on the ancestral lines to regress against the species data
# I can infer this extra information, 
# but can't compare them as the ancestral data is based on the species data its inferred from
# Not two datasets involved, so can't really statistically compare - bias as the underlying data is the same




# Agnese code for pairwise comparison

# Make x values (size) based on extant taxa
# Generating a sequence of logCS values between the min and max for the extant species
anc_X <- expand.grid(logCS =  seq(min(logCS), max(logCS), length.out = 165)) #use min and max of x values as limits and use number of specimens as length of sequence        

# Calculate Y values (pc scores) for ancestral nodes using anc.ML intercept and slope values - using the values_traj object
values_traj

anc_mammal_Y <- -1.0438379 + (anc_X*0.20702975)
anc_laurasiatheria_Y <- -1.0775600 + (anc_X*0.21592353)
anc_euarchontoglires_Y <- -1.1308218 + (anc_X*0.22873477)
anc_afrotheria_Y <- -1.1083497 + (anc_X*0.22633799)
anc_xenarthra_Y <- -1.0273031 + (anc_X*0.20481911)
anc_marsupialia_Y <- -0.7201908 + (anc_X*0.16069839)


# Rename the column names to be the same as the allometry object (allometry_plot_data)
colnames(anc_X) <- "logCsize"
colnames(anc_mammal_Y) <- "RegScores"
colnames(anc_laurasiatheria_Y) <- "RegScores"
colnames(anc_euarchontoglires_Y) <- "RegScores"
colnames(anc_afrotheria_Y) <- "RegScores"
colnames(anc_xenarthra_Y) <- "RegScores"
colnames(anc_marsupialia_Y) <- "RegScores"

# Make data frames for Y values
anc_Y <- rbind(anc_mammal_Y,anc_laurasiatheria_Y,anc_euarchontoglires_Y, anc_afrotheria_Y, anc_xenarthra_Y,anc_marsupialia_Y)  

# Create groups and order variables for ancestral states data
anc_nodes <- cbind(c(rep("ancestral mammal", 165),rep("ancestral laurasiatheria", 165),rep("ancestral euarchontoglires", 165), 
                     rep("ancestral afrotheria", 165),rep("ancestral xenarthra", 165),rep("ancestral marsupialia", 165)))
anc_groups <- rep("ancestral nodes", 990)

# Create data frame with PC scores, logCS estimated and groups and orders for anc to match allometry tibble
allometry_nodes <- data.frame(logCS = anc_X, anc_Y, species = anc_nodes, group = anc_groups)
allometry_nodes <- as_tibble(allometry_nodes)
glimpse(allometry_nodes)

# Make smaller data frame with only relevant data
# Data for the tips of the phylogeny i.e. species in dataset
allometry_tips <- allometry_plot_data
allometry_tips <- allometry_tips %>% select(c(1:2))

# Add labels and other attributes to tibble as columns
allometry_tips <- allometry_tips %>% 
  mutate(group = info$Major_clades, species = gdf$Species)
glimpse(allometry_tips)

# Make one dataset combining tips and anc
allometry_anc_all <- bind_rows(allometry_tips, allometry_nodes)
glimpse(allometry_anc_all)

##Order values by genus, useful for plot legend
#Make factor for variable
allometry_anc_all$species <- factor(allometry_anc_all$species, 
                                         levels = c("Bettongia penicillata", "Bradypus tridactylus", "Cebus apella",           
                                                    "Cyclopes didactylus", "Dasyprocta leporina", "Dasypus novemcinctus",      
                                                    "Epomops franqueti", "Felis catus", "Macroscelides proboscideus",
                                                    "Manis tricuspis", "Microcebus murinus", "Monodelphis domestica", 
                                                    "Mus musculus", "Ornithorhynchus anatinus", "Phacochoerus aethiopicus",
                                                    "Phascolarctos cincereus", "Rattus rattus", "Setifer setosus",        
                                                    "Setonix brachyurus", "Sminthopsis macroura", "Talpa europaea",      
                                                    "Trichosaurus vulpecula", "ancestral afrotheria", "ancestral euarchontoglires",
                                                    "ancestral laurasiatheria", "ancestral mammal", "ancestral marsupialia",  
                                                    "ancestral xenarthra")) #copy from string printed with the code above
#Order
allometry_anc_all <- allometry_anc_all[order(allometry_anc_all$species),]
#Check
glimpse(allometry_anc_all)


# Iteratively create two species pairs and perform Procrustes ANOVA and save output data from comparison
PairwiseComparisons=list()
SlopeDifferences=list()
InterceptDifferences=list()
g=1
n=2
for (g in 1:length(allometry_anc_all$species)){
  for (h in n:length(allometry_anc_all$species)){
    if (g==length(allometry_anc_all$species)){
      break
    }else
      
      SpeciesA <- allometry_anc_all$species[g]
    SpeciesB <- allometry_anc_all$species[h]
    toMatch <- c(SpeciesA, SpeciesB)
    filename <- paste(SpeciesA, SpeciesB, sep=" vs. ")
    
    two.species <- abind(allometry_anc_all$RegScores[[SpeciesA]],allometry_anc_all$RegScores[[SpeciesB]])
    two.species.csize <- c(allometry_anc_all$logCsize[[SpeciesA]],allometry_anc_all$logCsize[[SpeciesB]])
    two.species.species <- as.factor(c(paste(allometry_anc_all$species[[SpeciesA]]),paste(allometry_anc_all$species[[SpeciesB]])))
    two.species.slope.diff<-SlopesList[[g]]-SlopesList[[h]]
    names(two.species.slope.diff) <- paste("PC", c(1:28), " slope diff.", sep="")
    two.species.intercept.diff<-InterceptsList[[g]]-InterceptsList[[h]]
    names(two.species.intercept.diff) <- paste("PC", c(1:28), " int. diff.", sep="")
    
    Output <- procD.lm(two.species~two.species.csize, ~two.species.species, logsz = TRUE, iter = 9999)##Remember to add more iterations##
    
    PairwiseComparisons[[filename]]<-Output
    SlopeDifferences[[filename]]<-two.species.slope.diff
    InterceptDifferences[[filename]]<-two.species.intercept.diff
    
  }
  n=n+1
}


##Pairwise comparison of regression model between groups
#Create models, with different slopes and int or just int
allometry_anc_all_null <- lm.rrpp(RegScores ~ logCsize,
                                  data = allometry_anc_all, print.progress = FALSE, iter = 999) 
allometry_anc_all_comb <- lm.rrpp(RegScores ~ logCsize + species,
                         data = allometry_anc_all, print.progress = FALSE, iter = 999) 
allometry_anc_all_int <- lm.rrpp(RegScores ~ logCsize * species,
                        data = allometry_anc_all, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_anc_all_null)
summary(allometry_anc_all_comb)
summary(allometry_anc_all_int)

#Anova for difference between models
anova(allometry_anc_all_null, allometry_anc_all_comb, allometry_anc_all_int)

#Pairwise comparison of slopes for the 2 models PC1
pairwise_allometry_anc <- pairwise(allometry_anc_all_int, fit.null = allometry_anc_all_comb, groups = allometry_anc_all$species, 
                                   covariate =  allometry_anc_all$logCsize, print.progress = FALSE) 
pairwise_allometry_anc

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "dist") 

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "DL") 

#Compare the dispersion around group slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "var")

#Univariate p-values - = dist
slope.diff <- sapply(pairwise_allometry_anc$slopes, function(x) as.vector(dist(x)))
rownames(slope.diff) <- rownames(summary(pairwise_allometry_anc)$summary.table)
apply(slope.diff, 1, RRPP:::pval) # P-values
apply(slope.diff, 1, RRPP:::effect.size) # Z-scores

#Save results to file
sink("Output/pairwise_allometry_ancestral_nodes_genera.txt")
print("ANOVA models")
anova(allometry_anc_all_null, allometry_anc_all_comb, allometry_anc_all_int)

print("1-Pairwise absolute distances slopes")
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "dist") 

print("1a-Univariate p-values and z-scores - = dist")
apply(slope.diff, 1, RRPP:::pval) # P-values
apply(slope.diff, 1, RRPP:::effect.size) # Z-scores

print("2-Distance between angles (slope directions)")
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "VC", angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "DL") 

print("4-Difference in dispersion around mean slope")
summary(pairwise_allometry_anc, confidence = 0.95, test.type = "var")
sink()

#Improve plot using new dataset
allometry_anc_all_ggplot <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores, colour = genus, fill = genus))+
  geom_point(size = 0, colour = "white", fill  = "white")+
  geom_smooth(method = "lm", aes(fill = genus, colour = genus, linetype = group, size = group),           #confidence intervals and reg line, before points
             alpha = 0.4)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = NULL, labels = c("ancestral node 5", "Balaenoptera", "ancestral node 6", "Stenella", "ancestral node 7", 
                                              "Delphinapterus", "Phocoena" ), #to be ordered as they appear in tibble
                      values = mypalette_nodes, aesthetics = c("color","fill"))+           
  scale_linetype_manual(values = c(1,3,2))+
  scale_size_manual(values = c(1.3,1,1))+
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("Regression Score")+
  ggtitle ("Allometry by genus with ancestral nodes - p-value = 0.001**")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.text = element_text(size = 10), 
        legend.position = "bottom", legend.direction = "horizontal")+
  guides(linetype = guide_legend(label = F, title = NULL, override.aes = list(shape = NA, linetype = NA, fill = NA)), 
         size = guide_legend(label = F, title = NULL))
allometry_anc_all_ggplot
#Add phylopics
allometry_anc_all_ggplot  <- 
  allometry_anc_all_ggplot  + 
  add_phylopic(Balaenoptera, alpha = 1, x = 4, y = 0.12, ysize = 0.11, color = mypalette_taxa[1])+
  add_phylopic(Delphinapterus, alpha = 1, x = 2.9, y = -0.27, ysize = 0.15, color = mypalette_taxa[2])+
  add_phylopic(Phocoena, alpha = 1, x = 3.25, y = -0.01, ysize = 0.135, color = mypalette_taxa[3])+
  add_phylopic(Stenella, alpha = 1, x = 2.5, y = -0.05, ysize = 0.11, color = mypalette_taxa[4])
allometry_anc_all_ggplot


