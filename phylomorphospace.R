# Title: Phylomorphospace

# Name: Heather White

# Date created: 07/10/21

# Last modified: 07/10/21

# License: MIT license



# Phylomorphospace - plotting adult mammal specimens only

#######################################################################################################

rm(list = ls())

library(ggplot2)
library(ggConvexHull)
library(ggpubr)
library(RColorBrewer)
library(ape)
library(geomorph)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(Morpho)
library(rgl)
library(ape)
library(paleomorph)
library(RRPP)
library(arrayhelpers)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(devtools)
library(ggConvexHull)
library(ggthemes)
library(ggfortify)
library(ggphylomorpho)
library(gginnards)


#######################################################################################################

# 1) Read in the data

# Read in phylogeny
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
# Read in taxa names that match the phylogeny
tree_names <- read.csv("Data/tree_taxa_names.csv")
tree_names <- tree_names$Taxa_names

# Check the phylogeny
plot(my_tree, type = "fan") # plotted in a circular fan shape

# Check the phylogeny with timescale 
plot(my_tree, cex = 0.5) # plotted in normal linear tree fashion
axisPhylo()

# Load the data
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs_adults.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs_adults.Rdata")

# Read in the specimen details that match the dataset
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")


#######################################################################################################

# 2) Manipulate the data to the correct format

# Check the tree tip names match the shape data
shape.names <-dimnames(adults_sagittal)[[3]]
tree.names <-my_tree$tip.label
setdiff(shape.names,tree.names)
# Set the shape data to have the same names as the phylogeny
dimnames(adults_sagittal)[[3]] <- tree_names


#######################################################################################################

# 3) Perform PCA

# Perform PCA
PCA<-gm.prcomp(adults_sagittal)
# Obtain the proportion of variance for each PC
PCA_summary <- summary(PCA) 
PCA_summary <- PCA_summary$PC.summary
# Obtain the PC scores for each species for each PC
PCscores <- PCA$x

# Convert PCA results to a tibble
# Associate specimen info with the PCA
PCAresults<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults<-PCAresults %>% mutate(Species = adult_specimen_details$Species, Phylogeny_names = tree_names, Subclass = adult_specimen_details$Subclass, Major_clade = adult_specimen_details$Major_clades) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Add tree name information to the tibble
pcscores_adults <- PCAresults %>%  mutate(.,ID = tree_names)
# Sort the table by custom Clade order 
pcscores_adults$Major_clade <- factor(pcscores_adults$Major_clade, levels = c("Afrotheria","Euarchontoglires", "Laurasiatheria", "Marsupialia", "Monotremata","Xenarthra"))


#######################################################################################################

# 4) Plot the phylomorphospace

# Create a palette for the clades
my_palette_clade <- c("mediumpurple3", # Afrotheria
                      "#A6D854", # Euarchontoglires
                      "sandybrown", # Laurasiatheria
                      "cornflowerblue", # Marsupialia
                      "#FFD92F", # Monotremata
                      "palevioletred") # Xenarthra


# Plot the tree and tip data
g <- ggphylomorpho(tree = my_tree, tipinfo = pcscores_adults, xvar=Comp1, yvar = Comp2, factorvar = Major_clade,labelvar = ID, tree.alpha = 0.7, edge.width = 0.5)
# Remove specimen IDs and point data - to make pretty in next step
g <- g %>%
  delete_layers("GeomPoint") %>%
  delete_layers("GeomTextRepel")
# Add information to make the plot pretty
g + geom_convexhull(data = pcscores_adults, aes(x=Comp1, y = Comp2, colour = Major_clade, fill = Major_clade), alpha = 0.4, show.legend = F) +
  geom_point(data = pcscores_adults, aes(x=Comp1, y = Comp2,colour = Major_clade, fill = Major_clade, shape = Subclass),size=3)+#,show.legend = FALSE)
  geom_text_repel(data = pcscores_adults, aes(x=Comp1, y=Comp2, label=Species, hjust=0,vjust=0, fontface ="italic"), size = 5, show.legend = FALSE) +
  scale_colour_manual(values = my_palette_clade) +
  scale_fill_manual(values = my_palette_clade) +
  ggtitle("Phylomorphospace - PC1 and PC2") +
  guides(color=guide_legend("Clade"), fill = FALSE) +
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  theme_classic(base_size = 18)

