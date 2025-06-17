# Title: Suture complexity analysis

# Name: Heather White

# Date created: 17/12/19

# Last modified: 01/12/21

# License: MIT license


# Statistics on suture complexity, graphs to show suture complexity change with age etc.

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ape)
library(geomorph)
library(phytools)
library(ggpubr)
library(wesanderson)


#######################################################################################################

# Read in the .csv file of complexity results
complexity_results <- read.csv("Analysis/Complexity_results/combined_complexity_results_slided.csv")
complexity_results_adults <- read.csv("Analysis/Complexity_results/combined_complexity_results_adults_slided.csv")

# Read in the results where the complexity results ordered with the suture name listed - for model 13
complexity_results_suture <- read.csv("Analysis/Complexity_results/combined_complexity_results_by_suture_slided.csv")

# Read in the specimen info
specimen_info <- read.csv("Data/Specimen_info.csv")
adult_info <- read.csv("Data/Specimen_info_adults.csv")

# Combine the complexity results and specimen info into one dataframe
full_dataset <- cbind(complexity_results, specimen_info)


#######################################################################################################

# Summarise the data

# Calculate the mean, variance, and range for the sutures across the whole dataset
summary <- full_dataset %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal),
            varPSDsagittal = var(PSD_sagittal),
            varFDsagittal = var(FD_sagittal),
            varPSDcoronal = var(PSD_coronal),
            varFDcoronal = var(FD_coronal),
            varPSDinterfrontal = var(PSD_interfrontal),
            varFDinterfrontal = var(FD_interfrontal),
            rangePSDsagittal = range(PSD_sagittal),
            rangeFDsagittal = range(FD_sagittal),
            rangePSDcoronal = range(PSD_coronal),
            rangeFDcoronal = range(FD_coronal),
            rangePSDinterfrontal = range(PSD_interfrontal),
            rangeFDinterfrontal = range(FD_interfrontal))

write.csv(summary, file = "Analysis/Complexity_results/Statistics/complexity_summary_statistics_full_dataset.csv")


# Calculate the summary statistics for each species
summary_species <- full_dataset %>%
  group_by(Species) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal),
            varPSDsagittal = var(PSD_sagittal),
            varFDsagittal = var(FD_sagittal),
            varPSDcoronal = var(PSD_coronal),
            varFDcoronal = var(FD_coronal),
            varPSDinterfrontal = var(PSD_interfrontal),
            varFDinterfrontal = var(FD_interfrontal),
            rangePSDsagittal = range(PSD_sagittal),
            rangeFDsagittal = range(FD_sagittal),
            rangePSDcoronal = range(PSD_coronal),
            rangeFDcoronal = range(FD_coronal),
            rangePSDinterfrontal = range(PSD_interfrontal),
            rangeFDinterfrontal = range(FD_interfrontal))

write.csv(summary_species, file = "Analysis/Complexity_results/Statistics/complexity_summary_statistics_by_species.csv")

# Calculate the summary statistics for each species
summary_clade <- full_dataset %>%
  group_by(Major_clades) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal),
            varPSDsagittal = var(PSD_sagittal),
            varFDsagittal = var(FD_sagittal),
            varPSDcoronal = var(PSD_coronal),
            varFDcoronal = var(FD_coronal),
            varPSDinterfrontal = var(PSD_interfrontal),
            varFDinterfrontal = var(FD_interfrontal),
            rangePSDsagittal = range(PSD_sagittal),
            rangeFDsagittal = range(FD_sagittal),
            rangePSDcoronal = range(PSD_coronal),
            rangeFDcoronal = range(FD_coronal),
            rangePSDinterfrontal = range(PSD_interfrontal),
            rangeFDinterfrontal = range(FD_interfrontal))

write.csv(summary_clade, file = "Analysis/Complexity_results/Statistics/complexity_summary_statistics_by_clade.csv")


# Calculate the summary statistics for each species
summary_mars_plac <- full_dataset %>%
  group_by(Mars.plac) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal),
            varPSDsagittal = var(PSD_sagittal),
            varFDsagittal = var(FD_sagittal),
            varPSDcoronal = var(PSD_coronal),
            varFDcoronal = var(FD_coronal),
            varPSDinterfrontal = var(PSD_interfrontal),
            varFDinterfrontal = var(FD_interfrontal),
            rangePSDsagittal = range(PSD_sagittal),
            rangeFDsagittal = range(FD_sagittal),
            rangePSDcoronal = range(PSD_coronal),
            rangeFDcoronal = range(FD_coronal),
            rangePSDinterfrontal = range(PSD_interfrontal),
            rangeFDinterfrontal = range(FD_interfrontal))

write.csv(summary_mars_plac, file = "Analysis/Complexity_results/Statistics/complexity_summary_statistics_by_clade.csv")


# Calculate the summary statistics for discrete age categories - E, I, SA, A
summary_age <- full_dataset %>%
  group_by(Discrete_age) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal),
            varPSDsagittal = var(PSD_sagittal),
            varFDsagittal = var(FD_sagittal),
            varPSDcoronal = var(PSD_coronal),
            varFDcoronal = var(FD_coronal),
            varPSDinterfrontal = var(PSD_interfrontal),
            varFDinterfrontal = var(FD_interfrontal),
            rangePSDsagittal = range(PSD_sagittal),
            rangeFDsagittal = range(FD_sagittal),
            rangePSDcoronal = range(PSD_coronal),
            rangeFDcoronal = range(FD_coronal),
            rangePSDinterfrontal = range(PSD_interfrontal),
            rangeFDinterfrontal = range(FD_interfrontal))

write.csv(summary_age, file = "Analysis/Complexity_results/Statistics/complexity_summary_statistics_by_age.csv")

#######################################################################################################

# Compare the FD and PSD results 
# Are they correlated?

# Each suture seperately
sagittal_method_cor <- cor.test(complexity_results$PSD_sagittal, complexity_results$FD_sagittal, method = "spearman", conf.level = .95)
coronal_method_cor <- cor.test(complexity_results$PSD_coronal, complexity_results$FD_coronal, method = "spearman", conf.level = .95)
interfrontal_method_cor <- cor.test(complexity_results$PSD_interfrontal, complexity_results$FD_interfrontal, method = "spearman", conf.level = .95)

# Combining the sutures and just comparing PSD and FD
PSD_FD_cor <- cor.test(complexity_results_suture$PSD, complexity_results_suture$FD, method = "spearman", conf.level = .95)

sink("Analysis/Complexity_results/Statistics/correlation_between_PSD_and_FD_Spearmans_rank.txt")
print(sagittal_method_cor)
print(coronal_method_cor)
print(interfrontal_method_cor)
print(PSD_FD_cor)
sink() 


#######################################################################################################

# Statistics


##### Sagittal suture vs species

# Fit a linear model
model1 <- lm(PSD_sagittal ~ Species, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model1)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model1)
summary(model1) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Coronal suture vs species

# Fit a linear model
model2 <- lm(PSD_coronal ~ Species, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model2)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model2)
summary(model2) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Interfrontal suture vs species

# Fit a linear model
model3 <- lm(PSD_interfrontal ~ Species, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model3)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model3)
summary(model3) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Sagittal suture vs age

# Fit a linear model
model4 <- lm(PSD_sagittal ~ Discrete_age, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model4)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model4)
summary(model4) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Coronal suture vs age

# Fit a linear model
model5 <- lm(PSD_coronal ~ Discrete_age, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model5)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model5)
summary(model5) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Interfrontal suture vs age

# Fit a linear model
model6 <- lm(PSD_interfrontal ~ Discrete_age, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model6)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model6)
summary(model6) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Sagittal suture vs developmental strategy

# Fit a linear model
model7 <- lm(PSD_sagittal ~ Dev_strategy, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model7)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model7)
summary(model7) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Coronal suture vs developmental strategy

# Fit a linear model
model8 <- lm(PSD_coronal ~ Dev_strategy, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model8)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model8)
summary(model8) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Interfrontal suture vs developmental strategy

# Fit a linear model
model9 <- lm(PSD_interfrontal ~ Dev_strategy, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model9)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model9)
summary(model9) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Sagittal suture vs diet category

# Fit a linear model
model10 <- lm(PSD_sagittal ~ Diet_3_cats, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model10)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model10)
summary(model10) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Coronal suture vs diet category

# Fit a linear model
model11 <- lm(PSD_coronal ~ Diet_3_cats, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model11)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model11)
summary(model11) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Interfrontal suture vs diet category

# Fit a linear model
model12 <- lm(PSD_interfrontal ~ Diet_3_cats, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model12)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model12)
summary(model12) # This is saying which species lines are significantly different from 0 - a flat line

######


###### Sagittal suture vs coronoid

# Fit a linear model
model13 <- lm(PSD_sagittal ~ Coronoid, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model13)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model13)
summary(model13) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Coronal suture vs coronoid size

# Fit a linear model
model14 <- lm(PSD_coronal ~ Coronoid, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model14)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model14)
summary(model14) # This is saying which species lines are significantly different from 0 - a flat line

######

###### Interfrontal suture vs coronoid size

# Fit a linear model
model15 <- lm(PSD_interfrontal ~ Coronoid, data = full_dataset)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model15)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model15)
summary(model15) # This is saying which species lines are significantly different from 0 - a flat line


###### Difference across the 3 sutures

# Fit a linear model
model13 <- lm(PSD ~ Suture, data = complexity_results_suture)
# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(model13)
par(mfrow = c(1,1))
# Apply anova to model - anova table to get the answer to the question (sigificance)
anova(model13)
summary(model13) # This is saying which species lines are significantly different from 0 - a flat line



#######################################################################################################

# pANOVA

# Read in phylogeny
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
tree_names <- read.csv("Data/tree_taxa_names.csv")

# Save the factors for pANOVA as vectors and names them with the tree tip.labels
diet <- as.vector(complexity_results_adults$Diet)
names(diet) <- tree_names$Taxa_names
diet_3cats <- as.vector(complexity_results_adults$Diet_3)
names(diet_3cats) <- tree_names$Taxa_names
coronoid <- as.vector(complexity_results_adults$Coronoid)
names(coronoid) <- tree_names$Taxa_names
dev <- as.vector(complexity_results_adults$Precocial_altricial_spectrum)
names(dev) <- tree_names$Taxa_names
sag_PSD <- as.vector(complexity_results_adults$PSD_sagittal)
names(sag_PSD) <- tree_names$Taxa_names
cor_PSD <- as.vector(complexity_results_adults$PSD_coronal)
names(cor_PSD) <- tree_names$Taxa_names
IF_PSD <- as.vector(complexity_results_adults$PSD_interfrontal)
names(IF_PSD) <- tree_names$Taxa_names

# pANOVA for sagittal suture vs diet
pAN_sagittal_diet <- phylANOVA(my_tree, diet, sag_PSD, nsim=1000, posthoc=F)
pAN_sagittal_diet

# pANOVA for sagittal suture vs diet (divided into 3 categories)
pAN_sagittal_diet_3cats <- phylANOVA(my_tree, diet_3cats, sag_PSD, nsim=1000, posthoc=F)
pAN_sagittal_diet_3cats

# pANOVA for sagittal suture vs dev strategy
pAN_sagittal_dev <- phylANOVA(my_tree, dev, sag_PSD, nsim=1000, posthoc=F)
pAN_sagittal_dev 

# pANOVA for sagittal suture vs coronoid
pAN_sagittal_coronoid <- phylANOVA(my_tree, coronoid, sag_PSD, nsim=1000, posthoc=F)
pAN_sagittal_coronoid

#####

# pANOVA for coronal suture vs diet
pAN_coronal_diet <- phylANOVA(my_tree, diet, cor_PSD, nsim=1000, posthoc=F)
pAN_coronal_diet

# pANOVA for coronal suture vs diet (divided into 3 categories)
pAN_coronal_diet_3cats <- phylANOVA(my_tree, diet_3cats, cor_PSD, nsim=1000, posthoc=F)
pAN_coronal_diet_3cats

# pANOVA for coronal suture vs dev strategy
pAN_coronal_dev <- phylANOVA(my_tree, dev, cor_PSD, nsim=1000, posthoc=F)
pAN_coronal_dev 

# pANOVA for sagittal suture vs coronoid
pAN_coronal_coronoid <- phylANOVA(my_tree, coronoid, cor_PSD, nsim=1000, posthoc=F)
pAN_coronal_coronoid

#####

# pANOVA for interfrontal suture vs diet
pAN_IF_diet <- phylANOVA(my_tree, diet, IF_PSD, nsim=1000, posthoc=F)
pAN_IF_diet

# pANOVA for interfrontal suture vs diet
pAN_IF_diet_3cats <- phylANOVA(my_tree, diet_3cats, IF_PSD, nsim=1000, posthoc=F)
pAN_IF_diet_3cats

# pANOVA for coronal suture vs dev strategy
pAN_IF_dev <- phylANOVA(my_tree, dev, IF_PSD, nsim=1000, posthoc=F)
pAN_IF_dev 

# pANOVA for sagittal suture vs coronoid
pAN_IF_coronoid <- phylANOVA(my_tree, coronoid, IF_PSD, nsim=1000, posthoc=F)
pAN_IF_coronoid


sink("Analysis/Complexity_results/Statistics/pANOVAs.txt")
print("pANOVA for sagittal suture complexity vs diet category")
print(pAN_sagittal_diet)
print("pANOVA for sagittal suture complexity vs diet category (3 categories - herbivore, faunivore, omnivore)")
print(pAN_sagittal_diet_3cats)
print("pANOVA for sagittal suture complexity vs developmental strategy")
print(pAN_sagittal_dev)
print("pANOVA for coronal suture complexity vs diet category")
print(pAN_coronal_diet)
print("pANOVA for coronal suture complexity vs diet category (3 categories - herbivore, faunivore, omnivore)")
print(pAN_coronal_diet_3cats)
print("pANOVA for coronal suture complexity vs developmental strategy")
print(pAN_coronal_dev)
print("pANOVA for interfrontal suture complexity vs diet category")
print(pAN_IF_diet)
print("pANOVA for interfrontal suture complexity vs diet category (3 categories - herbivore, faunivore, omnivore)")
print(pAN_IF_diet_3cats)
print("pANOVA for interfrontal suture complexity vs developmental strategy")
print(pAN_IF_dev)
sink() 



# No difference in significance results if dividing the diet categories into 3 (herbivore, faunivore, carnivore)

#######################################################################################################

# Plot complexity results - scatter graphs

##########

# Setup the data for the plots

# Assigning a colour to each major clade - NEW PALETTE
col.clade = c("mediumpurple3", # Afrotheria
              "#A6D854", # Euarchontoglires
              "sandybrown", # Laurasiatheria
              "cornflowerblue", # Marsupialia
              "#FFD92F", # Monotremata
              "palevioletred") # Xenarthra
image(1:6, 1, as.matrix(1:6), col = col.clade, xlab = "Clade",
      ylab = "", yaxt = "n")

# Age colour
col.age <- c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF")

col.dev <- c("sandybrown", "dodgerblue4", "tan4", "lightskyblue3", "bisque")

col.diet <- wes_palette("Cavalcanti1", n = 3)
image(1:3, 1, as.matrix(1:3), col = col.diet, xlab = "Diet",
      ylab = "", yaxt = "n")

# Create a colour palette for species
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                        "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                        "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                        "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                        "gold1", # Monotremata
                        "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")
# Create a colour palette for species - embryos
mypalette_species_embryos <- c("mediumpurple1", # Afrotheria
                               "darkolivegreen1", "lawngreen", # Euarchontoglires
                               "burlywood1", "sandybrown", "darkorange2", "orangered1",  # Laurasiatheria
                               "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:16, 1, as.matrix(1:16), col = mypalette_species_embryos, xlab = "Species",
      ylab = "", yaxt = "n")
# Create a colour palette for species - infants
mypalette_species_infants <- c("mediumpurple1", "darkorchid4", # Afrotheria
                               "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                               "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                               "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                               "gold1", # Monotremata
                               "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species_infants, xlab = "Species",
      ylab = "", yaxt = "n")
# Create a colour palette for species - subadults
mypalette_species_SA <- c("mediumpurple1", "darkorchid4", # Afrotheria
                          "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                          "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                          "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred") # Xenarthra
image(1:20, 1, as.matrix(1:20), col = mypalette_species_SA, xlab = "Species",
      ylab = "", yaxt = "n")
# Create a colour palette for species - adults
mypalette_species_adults <- c("mediumpurple1", "darkorchid4", # Afrotheria
                              "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                              "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                              "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                              "gold1", # Monotremata
                              "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species_adults, xlab = "Species",
      ylab = "", yaxt = "n")

# Convert the dataset dataframe to a tibble so I can reorder the species
full_dataset_tbl <- as_tibble(full_dataset)
complexity_results_suture_tbl <- as_tibble(complexity_results_suture)

# Divide the data into the four age categories for plotting
PSD_embryos <- complexity_results_suture_tbl %>%
  filter(Age == "E")
PSD_embryos <- as_tibble(PSD_embryos)
PSD_infants <- complexity_results_suture_tbl %>%
  filter(Age == "I")
PSD_infants <- as_tibble(PSD_infants)
PSD_subadults <- complexity_results_suture_tbl %>%
  filter(Age == "SA")
PSD_subadults <- as_tibble(PSD_subadults)
PSD_adults <- complexity_results_suture_tbl %>%
  filter(Age == "A")
PSD_adults <- as_tibble(PSD_adults)

# Order the species to match the colours above
full_dataset_tbl$Species <- factor(full_dataset_tbl$Species,
                                 levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                            "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                            "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                            "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                            "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                            "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                            "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                            "Dasypus novemcinctus"),
                                 ordered = TRUE)
# Order the species to match the colours above
complexity_results_suture_tbl$Species <- factor(complexity_results_suture_tbl$Species,
                                   levels = c("Macroscelides proboscideus", "Setifer setosus", 
                                              "Sapajus apella", "Dasyprocta leporina", "Microcebus murinus", "Mus musculus", "Rattus rattus", 
                                              "Epomops franqueti", "Felis catus", "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                              "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                              "Ornithorhynchus anatinus", 
                                              "Bradypus tridactylus", "Cyclopes didactylus", "Dasypus novemcinctus"),
                                   ordered = TRUE)
# Order the species to match the colours above
PSD_embryos$Species <- factor(PSD_embryos$Species,
                                                levels = c("Macroscelides proboscideus", 
                                                           "Sapajus apella", "Dasyprocta leporina", 
                                                           "Epomops franqueti", "Felis catus", "Phataginus tricuspis", "Phacochoerus aethiopicus",          
                                                           "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                                           "Bradypus tridactylus", "Cyclopes didactylus", "Dasypus novemcinctus"),
                                                ordered = TRUE)
# Order the species to match the colours above
PSD_infants$Species <- factor(PSD_infants$Species,
                                                levels = c("Macroscelides proboscideus", "Setifer setosus", 
                                                           "Sapajus apella", "Dasyprocta leporina", "Microcebus murinus", "Mus musculus", "Rattus rattus", 
                                                           "Epomops franqueti", "Felis catus", "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                                           "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                                           "Ornithorhynchus anatinus", 
                                                           "Bradypus tridactylus", "Cyclopes didactylus", "Dasypus novemcinctus"),
                                                ordered = TRUE)
# Order the species to match the colours above
PSD_subadults$Species <- factor(PSD_subadults$Species,
                              levels = c("Macroscelides proboscideus", "Setifer setosus", 
                                         "Sapajus apella", "Dasyprocta leporina", "Microcebus murinus", "Mus musculus", "Rattus rattus", 
                                         "Felis catus", "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                         "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                         "Ornithorhynchus anatinus", 
                                         "Bradypus tridactylus", "Cyclopes didactylus"),
                              ordered = TRUE)
# Order the species to match the colours above
PSD_adults$Species <- factor(PSD_adults$Species,
                              levels = c("Macroscelides proboscideus", "Setifer setosus", 
                                         "Sapajus apella", "Dasyprocta leporina", "Microcebus murinus", "Mus musculus", "Rattus rattus", 
                                         "Epomops franqueti", "Felis catus", "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                         "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                         "Ornithorhynchus anatinus", 
                                         "Bradypus tridactylus", "Cyclopes didactylus", "Dasypus novemcinctus"),
                              ordered = TRUE)

#####################

# Scatter plots with continuous age on x-axis and PSD on y-axis, coloured by species
# One plot for each suture

# Scatter graph for sagittal suture PSD vs age (percent adult)
PSD_sagittal_species <- ggplot(full_dataset_tbl, aes(x = CS_percent_adult, y = PSD_sagittal, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 14)+
  scale_colour_manual(values = mypalette_species)+
  ylim(1.3, 2.3)+
  xlab("Age (% of adult CS)")+
  ylab("PSD score")+
  ggtitle("PSD complexity for sagittal suture")
PSD_sagittal_species

# Scatter graph for coronal suture PSD vs age (percent adult)
PSD_coronal_species <- ggplot(full_dataset_tbl, aes(x = CS_percent_adult, y = PSD_coronal, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 14)+
  scale_colour_manual(values = mypalette_species)+
  ylim(1.3, 2.3)+
  xlab("Age (% of adult CS)")+
  ylab("PSD score")+
  ggtitle("PSD complexity for coronal suture")
PSD_coronal_species

# Scatter graph for interfrontal suture PSD vs age (percent adult)
PSD_interfrontal_species <- ggplot(full_dataset_tbl, aes(x = CS_percent_adult, y = PSD_interfrontal, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 14)+
  scale_colour_manual(values = mypalette_species)+
  ylim(1.3, 2.3)+
  xlab("Age (% of adult CS)")+
  ylab("PSD score")+
  ggtitle("PSD complexity for interfrontal suture")
PSD_interfrontal_species


######################

# Scatter lots with three sutures on x-axis, PSD on y-axis, coloured by variables
# One plot for each variable (species, clade, age, diet, dev strategy)

# Scatter graph for the suture on x axis and species in colour
suture_comparison_species_scatter <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  ggtitle("PSD score for all three sutures, coloured by species")
suture_comparison_species_scatter

# Scatter graph for the suture on x axis and clade in colour
suture_comparison_clade_scatter <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, colour = Clade))+ 
  geom_point(size = 2)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_color_brewer(palette = "Set2")+
  ggtitle("PSD score for all three sutures, coloured by clade")
suture_comparison_clade_scatter

# Scatter graph for the suture on x axis and age in colour
suture_comparison_age_scatter <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, colour = Age_group))+ 
  geom_point(size = 2)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_color_brewer(palette = "Set2", labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                                                "Group3" = "Sub-adult", "Group4" = "Adult"))+
  ggtitle("PSD score for all three sutures, coloured by age")
suture_comparison_age_scatter

# Scatter graph for the suture on x axis and diet in colour
suture_comparison_diet_scatter <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, colour = Diet))+ 
  geom_point(size = 2)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_color_brewer(palette = "Set2")+
  ggtitle("PSD score for all three sutures, coloured by diet")
suture_comparison_diet_scatter

# Scatter graph for the suture on x axis and dev in colour
suture_comparison_dev_scatter <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, colour = Dev_strategy))+ 
  geom_point(size = 2)+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_color_brewer(palette = "Set2")+
  ggtitle("PSD score for all three sutures, coloured by developmental strategy")
suture_comparison_dev_scatter


###############

# Box plots showing the same as above - three sutures on x-axis, PSD on y-axis, coloured by variables
# One plot for each variable (species, clade, age, diet, dev strategy)

# Box plot for the suture on x axis and species in colour
suture_comparison_species_box <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, fill = Species))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_fill_manual(values = mypalette_species)+
  ggtitle("PSD score for all three sutures, coloured by species")
suture_comparison_species_box

# Box plot for the suture on x axis and clade in colour
suture_comparison_clade_box <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, fill = Clade))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_fill_manual(values = col.clade)+
  ggtitle("PSD score for all three sutures, coloured by clade")
suture_comparison_clade_box

# Box plot for the suture on x axis and age in colour
suture_comparison_age_box <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, fill = Age_group))+ 
  geom_boxplot()+
  theme_classic(base_size = 16)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_fill_manual(values = col.age, labels=c("Group1" = "Fetal", "Group2" = "Infant", 
                                               "Group3" = "Sub-adult", "Group4" = "Adult"))+
  ggtitle("PSD score for all three sutures, coloured by age")
suture_comparison_age_box

# Box plot for the suture on x axis and diet in colour
suture_comparison_diet_box <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, fill = Diet_3cats))+ 
  geom_boxplot()+
  theme_classic(base_size = 16)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_fill_manual(values = col.diet)+
  ggtitle("PSD score for all three sutures, coloured by diet")
suture_comparison_diet_box

# Box plot for the suture on x axis and dev in colour
suture_comparison_dev_box <- ggplot(complexity_results_suture_tbl, aes(x=Suture, y=PSD, fill = Dev_strategy))+ 
  geom_boxplot()+
  theme_classic(base_size = 16)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  scale_fill_manual(values = col.dev)+
  ggtitle("PSD score for all three sutures, coloured by developmental strategy")
suture_comparison_dev_box


#########

# Scatter plots dividing plot 'suture_comparison_species_box' into the 4 discrete age categories


# Scatter plot for the suture on x axis and species in colour, embyros only
suture_comparison_species_embryos <- ggplot(PSD_embryos, aes(x=Suture, y=PSD, fill = Species, colour = Species))+ 
  geom_dotplot(binwidth=0.012,binaxis = 'y', stackdir = 'center', position = position_dodge())+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  ylim(1.3, 2.3)+
  scale_fill_manual(values = mypalette_species_embryos)+
  scale_colour_manual(values = mypalette_species_embryos)+
  ggtitle("PSD score for all three sutures, coloured by species, embryos only")
suture_comparison_species_embryos

# Scatter plot for the suture on x axis and species in colour, infants only
suture_comparison_species_infants <- ggplot(PSD_infants, aes(x=Suture, y=PSD, fill = Species, colour = Species))+ 
  geom_dotplot(binwidth=0.012,binaxis = 'y', stackdir = 'center', position = position_dodge())+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  ylim(1.3, 2.3)+
  scale_fill_manual(values = mypalette_species_infants)+
  scale_colour_manual(values = mypalette_species_infants)+
  ggtitle("PSD score for all three sutures, coloured by species, infants only")
suture_comparison_species_infants

# Scatter plot for the suture on x axis and species in colour, subadults only
suture_comparison_species_subadults <- ggplot(PSD_subadults, aes(x=Suture, y=PSD, fill = Species, colour = Species))+ 
  geom_dotplot(binwidth=0.012,binaxis = 'y', stackdir = 'center', position = position_dodge())+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  ylim(1.3, 2.3)+
  scale_fill_manual(values = mypalette_species_SA)+
  scale_colour_manual(values = mypalette_species_SA)+
  ggtitle("PSD score for all three sutures, coloured by species, subadults only")
suture_comparison_species_subadults

# Scatter plot for the suture on x axis and species in colour, infants only
suture_comparison_species_adults <- ggplot(PSD_adults, aes(x=Suture, y=PSD, fill = Species, colour = Species))+ 
  geom_dotplot(binwidth=0.012,binaxis = 'y', stackdir = 'center', position = position_dodge())+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("PSD score")+
  ylim(1.3, 2.3)+
  scale_fill_manual(values = mypalette_species_adults)+
  scale_colour_manual(values = mypalette_species_adults)+
  ggtitle("PSD score for all three sutures, coloured by species, adults only")
suture_comparison_species_adults


#######################################################################################################

# Ancestral reconstruction

# Import trees in Nexus format - branch lengths needed
tree <- "Data/my_mammal_tree.nexus"  
# Read the tree for analysis
tree_species <- read.nexus(tree) 
Phylogeny_species_list <- read.csv("Data/tree_taxa_names.csv")

# Plot tree
plotTree(tree_species,ftype="i")

###########

# Ancestral reconstruction of suture complexity for adults only

sagittal <- complexity_results_adults$PSD_sagittal
coronal <- complexity_results_adults$PSD_coronal
interfrontal <- complexity_results_adults$PSD_interfrontal
# Give the data names that match the tree
names(sagittal) <- Phylogeny_species_list$Taxa_names
names(coronal) <- Phylogeny_species_list$Taxa_names
names(interfrontal) <- Phylogeny_species_list$Taxa_names

###########

# Sagittal suture

# Calculate ancestral states with phytools - slope and intercept for each component
# functions anc.ML and fastAnc give the same results here
anc_adult_sagittal_complexity <- anc.ML(tree_species, sagittal, CI = F)

# Plot the phylogeny with node labels to check where each is
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5) # add ancestral node information to the tree

# Combine the ancestral states with the original MSC_species results
combined <- data.frame(adult_sagittal_complexity = c(sagittal, anc_adult_sagittal_complexity$ace))

# Name the ancestral clade nodes
rownames(combined)[rownames(combined) == "23"] <- "Ancestral mammal"
rownames(combined)[rownames(combined) == "24"] <- "Ancestral therian mammal"
rownames(combined)[rownames(combined) == "25"] <- "Ancestral placental mammal"
rownames(combined)[rownames(combined) == "27"] <- "Ancestral Laurasiatheria"
rownames(combined)[rownames(combined) == "31"] <- "Ancestral Euarchontoglires"
rownames(combined)[rownames(combined) == "36"] <- "Ancestral Afrotheria"
rownames(combined)[rownames(combined) == "37"] <- "Ancestral Xenarthra"
rownames(combined)[rownames(combined) == "39"] <- "Ancestral Marsupialia"

# Save only the ancestral states I am interested in (the ones named above)
combined <- combined %>%
  slice(c(1:25,27,31,36,37,39))

write.csv(combined, file = "Analysis/Complexity_results/Ancestral_states/reconstructed_ancestral_states_sagittal_complexity_adults.csv")

# Map the traits onto the tree and plot
# contMap function maps the ancestral results using the fastAnc function
# The line above just gives the literal results
obj <- contMap(tree_species, sagittal, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))

##############

# Coronal suture

# Calculate ancestral states with phytools - slope and intercept for each component
# functions anc.ML and fastAnc give the same results here
anc_adult_coronal_complexity <- anc.ML(tree_species, coronal, CI = F)

# Plot the phylogeny with node labels to check where each is
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5) # add ancestral node information to the tree

# Combine the ancestral states with the original MSC_species results
combined <- data.frame(adult_coronal_complexity = c(coronal, anc_adult_coronal_complexity$ace))

# Name the ancestral clade nodes
rownames(combined)[rownames(combined) == "23"] <- "Ancestral mammal"
rownames(combined)[rownames(combined) == "24"] <- "Ancestral therian mammal"
rownames(combined)[rownames(combined) == "25"] <- "Ancestral placental mammal"
rownames(combined)[rownames(combined) == "27"] <- "Ancestral Laurasiatheria"
rownames(combined)[rownames(combined) == "31"] <- "Ancestral Euarchontoglires"
rownames(combined)[rownames(combined) == "36"] <- "Ancestral Afrotheria"
rownames(combined)[rownames(combined) == "37"] <- "Ancestral Xenarthra"
rownames(combined)[rownames(combined) == "39"] <- "Ancestral Marsupialia"

# Save only the ancestral states I am interested in (the ones named above)
combined <- combined %>%
  slice(c(1:25,27,31,36,37,39))

write.csv(combined, file = "Analysis/Complexity_results/Ancestral_states/reconstructed_ancestral_states_coronal_complexity_adults.csv")

# Map the traits onto the tree and plot
# contMap function maps the ancestral results using the fastAnc function
# The line above just gives the literal results
obj <- contMap(tree_species, coronal, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))

##############

# Interfrontal suture

# Calculate ancestral states with phytools - slope and intercept for each component
# functions anc.ML and fastAnc give the same results here
anc_adult_interfrontal_complexity <- anc.ML(tree_species, interfrontal, CI = F)

# Plot the phylogeny with node labels to check where each is
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5) # add ancestral node information to the tree

# Combine the ancestral states with the original MSC_species results
combined <- data.frame(adult_interfrontal_complexity = c(coronal, anc_adult_interfrontal_complexity$ace))

# Name the ancestral clade nodes
rownames(combined)[rownames(combined) == "23"] <- "Ancestral mammal"
rownames(combined)[rownames(combined) == "24"] <- "Ancestral therian mammal"
rownames(combined)[rownames(combined) == "25"] <- "Ancestral placental mammal"
rownames(combined)[rownames(combined) == "27"] <- "Ancestral Laurasiatheria"
rownames(combined)[rownames(combined) == "31"] <- "Ancestral Euarchontoglires"
rownames(combined)[rownames(combined) == "36"] <- "Ancestral Afrotheria"
rownames(combined)[rownames(combined) == "37"] <- "Ancestral Xenarthra"
rownames(combined)[rownames(combined) == "39"] <- "Ancestral Marsupialia"

# Save only the ancestral states I am interested in (the ones named above)
combined <- combined %>%
  slice(c(1:25,27,31,36,37,39))

write.csv(combined, file = "Analysis/Complexity_results/Ancestral_states/reconstructed_ancestral_states_interfrontal_complexity_adults.csv")

# Map the traits onto the tree and plot
# contMap function maps the ancestral results using the fastAnc function
# The line above just gives the literal results
obj <- contMap(tree_species, interfrontal, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))


#######################################################################################################

# HISTOGRAMS
# PROBABLY NOT THAT USEFUL!!!


# Plot complexity results - histograms

# Summarise the complexity data to get a mean value for each species discrete age category
summary_histogram <- full_dataset %>%
  group_by(Species, Discrete_age) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanFDsagittal = mean(FD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanFDcoronal = mean(FD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal),
            meanFDinterfrontal = mean(FD_interfrontal))

# Summary data for species only
summary_histogram_species <- full_dataset %>%
  group_by(Species) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal))

# Summary data for clade only
summary_histogram_clade <- full_dataset %>%
  group_by(Major_clades) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal))

# Summary data for age only
summary_histogram_age <- full_dataset %>%
  group_by(Discrete_age) %>%
  summarise(meanPSDsagittal = mean(PSD_sagittal),
            meanPSDcoronal = mean(PSD_coronal),
            meanPSDinterfrontal = mean(PSD_interfrontal))

# Histogram plot by species, coronal suture
PSD_species_histogram <- ggplot(summary_histogram, aes(Species, meanPSDcoronal)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_color_manual(values = mypalette_species, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  ylim(0,2.1)+
  xlab("Growth stage")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PSD_species_histogram

# Histogram plot by clade, coronal suture
PSD_clade_histogram <- ggplot(summary_histogram_clade, aes(Major_clades, meanPSDcoronal)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_color_manual(values = mypalette_species, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  ylim(0,2.1)+
  xlab("Growth stage")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PSD_clade_histogram

# Histogram plot by age, coronal suture
PSD_age_histogram <- ggplot(summary_histogram_age, aes(Discrete_age, meanPSDcoronal)) +
  geom_col(position = position_dodge2(padding = 0.2))+
  scale_x_discrete(labels=c("Group1" = "Embryo", "Group2" = "Infant", 
                            "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_color_manual(values = mypalette_species, aesthetics = c("color","fill"))+ 
  theme_classic(base_size = 12)+
  ylim(0,2.1)+
  xlab("Growth stage")+
  ylab("PV (Procrustes variances)")+
  ggtitle ("Disparity by species and age category")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
PSD_age_histogram


#########################################################################################################

# Correlation between suture complexity and skull size

#########

# ADULTS ONLY

# Plot the results - adults only
ggscatter(complexity_results_adults, x = "CS_logged", y = "PSD_sagittal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "PSD complexity score",
          title = "Correlation between suture sagittal suture complexity and skull size - adults only")

# Plot the results - adults only
ggscatter(complexity_results_adults, x = "CS_logged", y = "PSD_interfrontal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "PSD complexity score",
          title = "Correlation between suture interfrontal suture complexity and skull size - adults only")

# Plot the results - adults only
ggscatter(complexity_results_adults, x = "CS_logged", y = "PSD_coronal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PSD complexity score", 
          ylab = "Total suture closure score",
          title = "Correlation between suture coronal suture complexity and skull size - adults only")



# Spearmans rank correlation for adults only
cor_adults_sagittal <-cor.test(complexity_results_adults$PSD_sagittal, complexity_results_adults$CS_logged, method = "spearman", conf.level = .95, exact = F)
cor_adults_IF <-cor.test(complexity_results_adults$PSD_interfrontal, complexity_results_adults$CS_logged, method = "spearman", conf.level = .95, exact = F)
cor_adults_coronal <-cor.test(complexity_results_adults$PSD_coronal, complexity_results_adults$CS_logged, method = "spearman", conf.level = .95, exact = F)

sink("Analysis/Complexity_results/Skull_size/spearmans_rank_logCS_vs_suture_complexity_adults.txt")
print("Correlation between sagittal suture complexity and skull size")
print(cor_adults_sagittal)
print("Correlation between interfrontal suture complexity and skull size")
print(cor_adults_IF)
print("Correlation between coronal suture complexity and skull size")
print(cor_adults_coronal)
sink() 


#########

# ALL SPECIMENS

# Plot the results - all specimens
ggscatter(full_dataset, x = "CS_logged", y = "PSD_sagittal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "PSD complexity score",
          title = "Correlation between suture complexity and skull size - all specimens - sagittal suture")

# Plot the results - all specimens
ggscatter(full_dataset, x = "CS_logged", y = "PSD_interfrontal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "PSD complexity score",
          title = "Correlation between suture complexity and skull size - all specimens - interfrontal suture")

# Plot the results - all specimens
ggscatter(full_dataset, x = "CS_logged", y = "PSD_coronal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "PSD complexity score",
          title = "Correlation between suture complexity and skull size - all specimens - coronal suture")


# Spearmans rank correlation for all specimens
cor_all_sagittal <-cor.test(full_dataset$PSD_sagittal, full_dataset$CS_logged, method = "spearman", conf.level = .95)
cor_all_IF <-cor.test(full_dataset$PSD_interfrontal, full_dataset$CS_logged, method = "spearman", conf.level = .95)
cor_all_coronal <-cor.test(full_dataset$PSD_coronal, full_dataset$CS_logged, method = "spearman", conf.level = .95)

sink("Analysis/Complexity_results/Skull_size/spearmans_rank_logCS_vs_suture_complexity_all_specimens.txt")
print("Correlation between sagittal suture complexity and skull size")
print(cor_all_sagittal)
print("Correlation between interfrontal suture complexity and skull size")
print(cor_all_IF)
print("Correlation between coronal suture complexity and skull size")
print(cor_all_coronal)
sink() 








