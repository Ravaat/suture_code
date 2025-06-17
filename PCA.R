# Title: PCA

# Name: Heather White

# Date created: 20/08/21

# Last modified: 17/12/21

# License: MIT license



# PCA of all three sutures

#######################################################################################################


rm(list = ls())

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

# Load the data

# These are the suture LMs that have been mirrored, resampled, slid, and Procrusted based on LMs only
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")


# Reading in the specimen details
specimen_details <- read.csv("Data/Specimen_info.csv", header = T, sep=",")
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")


# Getting the adults only
adults_sagittal <- sagittal[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]
adults_coronal <- coronal[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]
adults_IF <- IF[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]

#save(adults_sagittal, file = "Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs_adults.Rdata")
#save(adults_coronal, file = "Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs_adults.Rdata")
#save(adults_IF, file = "Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs_adults.Rdata")



####################################################################################################

# PCA - for adults only


# To plot PCA quickly
PCA<-gm.prcomp(adults_coronal)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/Suture_morphology/PCA/PCA_sagittal_adults_results.csv")
# Write the PC loadings for every specimen to a .csv
write.csv(PCA$x, file = "Analysis/Suture_morphology/PCA/PC_loadings_sagittal_adults.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Getting the shape variation of the PC axes
shape <- mshape(adults_sagittal)
plotRefToTarget(shape, PCA$shapes$shapes.comp1$min) # PC1 min shape

plotRefToTarget(shape, PCA$shapes$shapes.comp1$max) # PC1 max shape


#######

# Prep colours for the PCA

# Create a colour palette
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:28, 1, as.matrix(1:28), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")

# Assigning a colour to each major clade
col.clade = c("mediumpurple3", # Afrotheria
              "#A6D854", # Euarchontoglires
              "sandybrown", # Laurasiatheria
              "cornflowerblue", # Marsupialia
              "#FFD92F", # Monotremata
              "palevioletred") # Xenarthra

# Age colours
col.age <- c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF")

col.dev <- c("sandybrown", "dodgerblue4", "tan4", "lightskyblue3", "bisque")

col.diet <- wes_palette("Cavalcanti1", n = 3)
image(1:3, 1, as.matrix(1:3), col = col.diet, xlab = "Diet",
      ylab = "", yaxt = "n")


#######

# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_adults<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_adults<-PCAresults_adults %>% mutate(Species = adult_specimen_details$Species, Subclass = adult_specimen_details$Subclass, Major_clade = adult_specimen_details$Major_clades, Dev_strategy = adult_specimen_details$Precocial_altricial_spectrum, Diet = adult_specimen_details$Diet, Diet3 = adult_specimen_details$Diet_3) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for adults only and labelling the specimens
PCA_adults_labelled <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, label=Species))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=2)+
  geom_text_repel(aes(fontface="italic"), size = 3)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  ggtitle("PCA adults labelled - sagittal suture")
PCA_adults_labelled

# Plot PCA for adults only and grouping specimens by subclass
PCA_adults_subclass <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Subclass, shape = Subclass))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=2)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("PC1 (29.37%)")+ # Found by using summary(PCA) above
  ylab("PC2 (18.37%)") +
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "Set2")+
  labs(color = "Subclass")+
  ggtitle("PCA adults divided by subclass")
PCA_adults_subclass

# Use this to view colour palette options:
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/


# Plot PCA for adults only and grouping specimens by major clade
PCA_adults_major_clade_subclass <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Major_clade, shape = Subclass))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=2)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "Set2")+
  labs(color = "Clade", shape = "Subclass")+
  ggtitle("PCA adults divided by subclass and major clade")
PCA_adults_major_clade_subclass

# Plot PCA for adults only and grouping specimens by major clade
PCA_adults_clade <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Major_clade, fill=Major_clade))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=2)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_fill_manual(values = col.clade)+
  scale_colour_manual(values = col.clade)+
  labs(color = "Clade")+
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  ggtitle("PCA adults divided by major clade")
PCA_adults_clade

# Plot PCA for adults only and grouping specimens by developmental strategy (altricial-precocial spectrum)
PCA_adults_dev_strategy <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Dev_strategy, fill = Dev_strategy))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=3)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = col.dev)+
  scale_fill_manual(values = col.dev)+
  labs(color = "Developmental Strategy")+
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  ggtitle("PCA adults divided by developmental strategy")
PCA_adults_dev_strategy

# Plot PCA for adults only and grouping specimens by diet type
PCA_adults_diet <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Diet, fill = Diet))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=2)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(color = "Diet")+
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  ggtitle("PCA adults divided by diet")
PCA_adults_diet

# Plot PCA for adults only and grouping specimens by diet type
PCA_adults_diet <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Diet3, fill = Diet3))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=3)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = col.diet)+
  scale_fill_manual(values = col.diet)+
  labs(color = "Diet")+
  ylim(-0.75,0.8)+
  xlim(-1.25,1)+
  ggtitle("PCA adults divided by diet")
PCA_adults_diet


########

# Histogram/screeplot of PCs for adults only

adult_histogram <- read.csv("Data/adult_histogram_data.csv")

# Method 1
PCA_adults_histogram <- ggplot(adult_histogram)+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_bar(stat = "identity", aes(x=Principal.Component, y = Proportion.of.Variance))+
  geom_line(aes(x=Principal.Component, y = Cumulative.Proportion, col ="red"))+
  geom_point(aes(x=Principal.Component, y = Cumulative.Proportion, col = "red"))+
  theme_classic()+
  ylim(0,1)+
  xlab("Principal Component")+ 
  ylab("Variances")+
  ggtitle("Screeplot")

# Method 2
barplot(adult_histogram[, 2], names.arg=1:nrow(adult_histogram), 
        main = "Screeplot",
        xlab = "Principal Components",
        ylab = "Variances",
        ylim = c(0,1.1),
        col ="grey")


####################################################################################################

# PCA - for all developmental specimens


# To plot PCA quickly
PCA<-gm.prcomp(IF)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/Suture_morphology/PCA/PCA_sagittal_all_specimen_results.csv")
# Write the PC loadings for every specimen to a .csv
write.csv(PCA$x, file = "Analysis/Suture_morphology/PCA/PC_loadings_sagittal_all_specimens.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Getting the shape variation of the PC axes
shape <- mshape(adults)
plotRefToTarget(PCA$shapes$shapes.comp3$min, shape) # PC1 min shape

plotRefToTarget(shape, PCA$shapes$shapes.comp3$max) # PC1 max shape

########

# Histogram/screeplot of PCs for adults only

full_histogram <- read.csv("Data/all_histogram_data.csv")

# Method 1
PCA_full_histogram <- ggplot(full_histogram)+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_bar(stat = "identity", aes(x=Principal.Component, y = Proportion.of.Variance))+
  geom_line(aes(x=Principal.Component, y = Cumulative.Variance, col ="red"))+
  theme_bw()+
  ylim(0,1)+
  xlab("Principal Component")+ 
  ylab("Variances")+
  ggtitle("Screeplot")

# Method 2
barplot(full_histogram[, 2], names.arg=1:nrow(full_histogram), 
        main = "Screeplot",
        xlab = "Principal Components",
        ylab = "Variances",
        ylim = c(0,1.1),
        col ="grey")


#######

# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_all<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_all<-PCAresults_all %>% mutate(Subclass = specimen_details$Subclass, Major_clade = specimen_details$Major_clades, Age = specimen_details$Age, Species = specimen_details$Species, Dev_strategy = specimen_details$Dev_strategy, New_Age = specimen_details$Discrete_age_group, Diet3 = specimen_details$Diet_3_cats) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Save PC scores with specimen info as an R object
save(PCAresults_all, file = "Data/PC_scores_sagittal_all_specimens.Rdata")

# Order the data by clade and convert to factor for plotting and plot legend
PCAresults_all$Major_clade <- as.factor(PCAresults_all$Major_clade)
PCAresults_all <- PCAresults_all %>% arrange(Major_clade)
PCAresults_all

# Order the species to match the colours above
PCAresults_all$Species <- factor(PCAresults_all$Species,
                                               levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                                          "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                                          "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                                          "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                                          "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                                          "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                                          "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                                          "Dasypus novemcinctus"),
                                               ordered = TRUE)



# Plot the PCA using ggplot all specimens - age not identifiable
PCA_all <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = Species))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = mypalette_species)+
  labs(color = "Clade")+
  xlim(-2,1.5)+
  ylim(-1,1.5)+
  ggtitle("PCA all species coloured")
PCA_all


# Plot the PCA using ggplot for all specimens and showing age
PCA_all_age1 <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "Set2")+
  labs(color = "Clade")+
  ggtitle("PCA all - divided by age only")
PCA_all_age1

# Plot PCA for all specimens showing age using groups E, I, SA, A
PCA_new_age <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour=New_Age, fill = New_Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=2)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = col.age, labels=c("Group1" = "Fetal", "Group2" = "Infant", 
                                                "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_fill_manual(values = col.age, labels=c("Group1" = "Fetal", "Group2" = "Infant", 
                                               "Group3" = "Sub-adult", "Group4" = "Adult"))+
  labs(color = "Age")+
  xlim(-2,1.5)+
  ylim(-1,1.5)+
  ggtitle("PCA all - divided by age category")
PCA_new_age

# Plot PCA for all specimens showing developmental strategy (altricial-precocial spectrum)
PCA_all_dev_strategy <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour=Dev_strategy, fill = Dev_strategy))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(size=2)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(color = "Developmental Strategy")+
  xlim(-2,1.5)+
  ylim(-1,1.5)+
  ggtitle("PCA all - divided by developmental strategy")
PCA_all_dev_strategy


# Plot the PCA using ggplot for all specimens and showing age and species
PCA_all_age2 <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_manual(values = mypalette_species)+
  labs(color = "Species")+
  xlim(-2,1.5)+
  ylim(-1,1.5)+
  ggtitle("PCA all - divided by age and coloured by species")
PCA_all_age2


# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_all_age3 <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = Major_clade, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Clade")+
  ggtitle("PCA all - divided by age and coloured by major clade")
PCA_all_age3



####################################################################################################

# PCA - including centroid size as the proxy for age

#######

# To plot PCA quickly
PCA<-gm.prcomp(coronal)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Getting the shape variation of the PC axes
shape <- mshape(adults)
plotRefToTarget(PCA$shapes$shapes.comp3$min, shape) # PC1 min shape

plotRefToTarget(shape, PCA$shapes$shapes.comp3$max) # PC1 max shape


#######

# PCA here is plotted with logged CS info (binned with the CS logged values)
# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_all<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_all<-PCAresults_all %>% mutate(Subclass = specimen_details$Subclass, Major_clade = specimen_details$Major_clades, Age = specimen_details$Age, Species = specimen_details$Species, CS = specimen_details$CS, CS_logged = specimen_details$CS_logged) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA
# Create 10 bins of the logged CS data
bins=cut(PCAresults_all$CS_logged, breaks=10, ordered_result = TRUE) 
# Add this binned info to the results spreadsheet
PCAresults_all<-PCAresults_all %>% mutate(CS_binned = bins)


# Plot the PCA using ggplot for all specimens and associating colour scale with logged centroid size
PCA_all_CS <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = CS_binned))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=2)+
  theme_classic(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_brewer(palette = "PuOr")+
  #labs(color = "Clade")+
  ggtitle("PCA all - CS across all specimens")
PCA_all_CS

# Centroid size logged for each species seperately produces the same values as when CS is logged with everything all together
# Nevertheless can see a seperation on the PCA for size of skulls (clearer when CS is logged than without logging)
# Otherwise the Phacochoerus skulls are the only big ones

#######

# PCA here is instead binned using the percent of adult (calculated using original CS)

# Associate specimen info with the PCA
PCAresults_all<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_all<-PCAresults_all %>% mutate(Subclass = specimen_details$Subclass, Major_clade = specimen_details$Major_clades, Age = specimen_details$Age, Species = specimen_details$Species, CS = specimen_details$CS, CS_percent_adult = specimen_details$CS_percent_adult) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA
# Create 10 bins of the logged CS data
bins=cut(PCAresults_all$CS_percent_adult, breaks=10, ordered_result = TRUE) 
# Add this binned info to the results spreadsheet
PCAresults_all<-PCAresults_all %>% mutate(CS_percent_binned = bins)


# Plot the PCA using ggplot for all specimens and associating colour scale with centroid size
PCA_all_CS_percent_adult <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = CS_percent_adult))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_viridis()+
  xlim(-2,1.5)+
  ylim(-1,1.5)+
  #labs(color = "Clade")+
  ggtitle("PCA all - CS (percent adult) across all specimens")
PCA_all_CS_percent_adult



################################################################################################################

# IF WANTED

# Plotting PCAs for the major clades seperately and including developmental information


Marsupialia_sagittal <- sagittal[,,c(1:10,71:83,90:93,105:115,129:148,157:165)]
Marsupialia_coronal <- coronal[,,c(1:10,71:83,90:93,105:115,129:148,157:165)]
Marsupialia_IF <- IF[,,c(1:10,71:83,90:93,105:115,129:148,157:165)]
Marsupialia_info <- specimen_details[c(1:10,71:83,90:93,105:115,129:148,157:165),]

Xenarthra_sagittal <- sagittal[,,c(11:17,24:27,36:41)]
Xenarthra_coronal <- coronal[,,c(11:17,24:27,36:41)]
Xenarthra_IF <- IF[,,c(11:17,24:27,36:41)]
Xenarthra_info <- specimen_details[c(11:17,24:27,36:41),]

Afrotheria_sagittal <- sagittal[,,c(54:59,121:128)]
Afrotheria_coronal <- coronal[,,c(54:59,121:128)]
Afrotheria_IF <- IF[,,c(54:59,121:128)]
Afrotheria_info <- specimen_details[c(54:59,121:128),]

Laurasiatheria_sagittal <- sagittal[,,c(42:53,60:66,94:104,149:156)]
Laurasiatheria_coronal <- coronal[,,c(42:53,60:66,94:104,149:156)]
Laurasiatheria_IF <- IF[,,c(42:53,60:66,94:104,149:156)]
Laurasiatheria_info <- specimen_details[c(42:53,60:66,94:104,149:156),]

Euarchontoglires_sagittal <- sagittal[,,c(18:23,28:35,67:70,84:89,116:120)]
Euarchontoglires_coronal <- coronal[,,c(18:23,28:35,67:70,84:89,116:120)]
Euarchontoglires_IF <- IF[,,c(18:23,28:35,67:70,84:89,116:120)]
Euarchontoglires_info <- specimen_details[c(18:23,28:35,67:70,84:89,116:120),]

Eutheria_sagittal <- sagittal[,,c(11:70,84:89,94:104,116:128,149:156)]
Eutheria_coronal <- coronal[,,c(11:70,84:89,94:104,116:128,149:156)]
Eutheria_IF <- IF[,,c(11:70,84:89,94:104,116:128,149:156)]
Eutheria_info <- specimen_details[c(11:70,84:89,94:104,116:128,149:156),]


######

# 1) Marsupialia - this is monotremes too - therefore Prototheria and Metatheria

# To plot PCA quickly
PCA<-gm.prcomp(Marsupialia_sagittal)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_marsupials_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_marsupials<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_marsupials<-PCAresults_marsupials %>% mutate(Species = Marsupialia_info$Species, Subclass = Marsupialia_info$Subclass, Major_clade = Marsupialia_info$Major_clades, Age = Marsupialia_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_marsupialia <- ggplot(PCAresults_marsupials, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Species")+
  ggtitle("PCA marsupials")

######

# 2) Xenarthra

# To plot PCA quickly
PCA<-gm.prcomp(Xenarthra)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_xenarthra_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_xenarthra<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_xenarthra<-PCAresults_xenarthra %>% mutate(Species = Xenarthra_info$Species, Subclass = Xenarthra_info$Subclass, Major_clade = Xenarthra_info$Major_clades, Age = Xenarthra_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_xenarthra <- ggplot(PCAresults_xenarthra, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Species")+
  ggtitle("PCA Xenarthra")

######

# 3) Afrotheria

# To plot PCA quickly
PCA<-gm.prcomp(Afrotheria)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_afrotheria_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_afrotheria<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_afrotheria<-PCAresults_afrotheria %>% mutate(Species = Afrotheria_info$Species, Subclass = Afrotheria_info$Subclass, Major_clade = Afrotheria_info$Major_clades, Age = Afrotheria_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_afrotheria <- ggplot(PCAresults_afrotheria, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Species")+
  ggtitle("PCA Afrotheria")

######

# 4) Laurasiatheria

# To plot PCA quickly
PCA<-gm.prcomp(Laurasiatheria)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_laurasiatheria_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_laurasiatheria<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_laurasiatheria<-PCAresults_laurasiatheria %>% mutate(Species = Laurasiatheria_info$Species, Subclass = Laurasiatheria_info$Subclass, Major_clade = Laurasiatheria_info$Major_clades, Age = Laurasiatheria_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_laurasiatheria <- ggplot(PCAresults_laurasiatheria, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Species")+
  ggtitle("PCA Laurasiatheria")

######

# 5) Euarchontoglires

# To plot PCA quickly
PCA<-gm.prcomp(Euarchontoglires)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_euarchontoglires_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_euarchontoglires<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_euarchontoglires<-PCAresults_euarchontoglires %>% mutate(Species = Euarchontoglires_info$Species, Subclass = Euarchontoglires_info$Subclass, Major_clade = Euarchontoglires_info$Major_clades, Age = Euarchontoglires_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_euarchontoglires <- ggplot(PCAresults_euarchontoglires, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+
  scale_color_brewer(palette = "Set2")+
  labs(color = "Species")+
  ggtitle("PCA Euarchontoglires")

######

# 6) Eutheria

# To plot PCA quickly
PCA<-gm.prcomp(Eutheria)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary
write.csv(PCA_summary, file ="Analysis/PCA_eutheria_results.csv")
plot(PCA, main = "PCA")
# Another way to plot PCA quickly
PCA.plot <-plot(PCA, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA$x[,1], y=PCA$x[,2], rownames(PCA$x))


# Plot PCA neatly using ggplot

# Associate specimen info with the PCA
PCAresults_eutheria<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults_eutheria<-PCAresults_eutheria %>% mutate(Species = Eutheria_info$Species, Subclass = Eutheria_info$Subclass, Major_clade = Eutheria_info$Major_clades, Age = Eutheria_info$Age) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Plot the PCA using ggplot for all specimens and showing age and Clade
PCA_eutheria <- ggplot(PCAresults_eutheria, aes(x=Comp1, y=Comp2, colour = Species, shape=Age))+ # To add in other groupings here, use shape =. This needs to have been added into a column in the above line first though
  geom_point(size=3)+
  theme_bw(base_size = 12)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  scale_shape_manual(values=c(16,4,1))+ 
  scale_color_manual(values = c("sandybrown", "darkolivegreen1","sienna2", "olivedrab3", 
                                "darkorange2", "lightsteelblue2", "skyblue", "plum3", "steelblue3", "darkolivegreen4", 
                                "forestgreen", "dodgerblue4", "darkgreen", "darkorchid3", "navyblue"))+
  labs(color = "Species")+
  ggtitle("PCA Eutheria")









