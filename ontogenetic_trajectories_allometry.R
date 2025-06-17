# Title: Ontogenetic trajectories - allometry

# Name: Heather White

# Date created: 02/09/21

# Last modified: 20/12/21

# License: MIT license



# Ontogenetic trajectories plotted for each species - using allometry regressions

#######################################################################################################


rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(abind)

#######################################################################################################

# Load the data 

# Load the PCA results data - PC scores
load(file ="Data/PC_scores_sagittal_all_specimens.Rdata")
load(file ="Data/PC_scores_coronal_all_specimens.Rdata")
load(file ="Data/PC_scores_interfrontal_all_specimens.Rdata")

# Load the suture data - mirrored, resampled, slid, Procrusted
load("Data/3D_Resampled_Procrusted_suture_LMs/sagittal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/coronal_Procrusted_LMs.Rdata")
load("Data/3D_Resampled_Procrusted_suture_LMs/interfrontal_Procrusted_LMs.Rdata")


info <- read.csv(file = "Data/Specimen_info_trajectories.csv", header =T, sep = ",")
CS <- info$CS
logCS <- log(info$CS)


#######################################################################################################

# Format the data 

# Set colours
# Run PlottingValues_Function first to do this
mammal.PlottingValues <- PlottingValues(X=info,ColorGroup = "Species",ShapeGroup = "Major_clades")


# Perform PCA to get PCA results and output
PCA<-gm.prcomp(coronal)
PCA_summary <- summary(PCA) # This gives all the info for the PCA
PCA_summary <- PCA_summary$PC.summary

# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
coronal <- unname(coronal)

# Get the raw Procrustes data into a geomorph dataframe
gdf <- geomorph.data.frame(coronal, CSize = CS, logCsize = logCS, Species = info$Species, Specimens = info$Specimen_name)
# name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","CSize", "logCsize", "Species", "Specimens")

# Assigning a colour to each species
col.species = c("paleturquoise1", "pink1", "darkolivegreen1", "palevioletred",
                "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
                "darkorange2", "green3", "turquoise1", "chartreuse4", "gold1",
                "orangered1", "cyan3", "darkgreen", "darkorchid4",
                "dodgerblue2", "blue2", "red3", "midnightblue")[gdf$Species]

# Convert the mammal.PlottingValues to the colours I want - matching above
col.gp.1 <- c("paleturquoise1", "pink1", "darkolivegreen1", "palevioletred",
              "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
              "darkorange2", "green3", "turquoise1", "chartreuse4", "gold1",
              "orangered1", "cyan3", "darkgreen", "darkorchid4",
              "dodgerblue2", "blue2", "red3", "midnightblue") 

names( col.gp.1 )   <- levels(info$Species)
mammal.PlottingValues$legendcolor <- col.gp.1


# Subset dataset by species
GPAList=list()
CovariatesList=list()
SpeciesList=list()
CladeList=list()
PCList=list()
Shapes=list()
CSList=list()
LogCsize=list()
Colors=list()
Taxa=list()
for (g in 1:length(levels(info[["Species"]]))){
  x<-NULL
  x<-as.array(gdf$Proc_coords[,,grep(levels(info[["Species"]])[g],info[["Species"]])])
  y<-NULL
  y<-info[grep(levels(info[["Species"]])[g],info[["Species"]]),]
  z<-NULL
  z<-as.matrix(PCA$x[grep(levels(info[["Species"]])[g],info[["Species"]]),])
  {
    GPAList[[print(levels(info[["Species"]])[g])]]<-as.array(x)
    CSList[[levels(info[["Species"]])[g]]]<-info$CS[grep(levels(info[["Species"]])[g],info[["Species"]])]
    LogCsize[[levels(info[["Species"]])[g]]] <- info$CS_logged[grep(levels(info[["Species"]])[g],info[["Species"]])]
    CovariatesList[[levels(info[["Species"]])[g]]]<-y
    SpeciesList[[levels(info[["Species"]])[g]]]<-y$Species
    CladeList[[levels(info[["Species"]])[g]]]<-y$Major_clade
    Shapes[[levels(info[["Species"]])[g]]]<-mammal.PlottingValues$shape[grep(levels(info[["Species"]])[g],info[["Species"]])]
    Colors[[levels(info[["Species"]])[g]]]<-mammal.PlottingValues$color[grep(levels(info[["Species"]])[g],info[["Species"]])]
    Taxa<-names(GPAList)
  }
  if (ncol(z)==1) {
    z <- t(z)
    PCList[[levels(info[["Species"]])[g]]]<-z
  }else{
    PCList[[levels(info[["Species"]])[g]]]<-z
  }
}

# Remove rownames of matrix - otherwise this confuses later analysis
for (i in 1:22){
rownames(PCList[[i]]) <- NULL
}

appended_info_species <- list(name = paste(info$Species[-23]), taxa = Taxa[-23], coords = GPAList[-23], PCvalues = PCList[-23], CSize = CSList[-23], logCsize = LogCsize[-23], covariates = CovariatesList[-23], species = SpeciesList[-23], clades = CladeList[-23], color = Colors[-23], shape = Shapes[-23])

A <- list()
for (i in 1:length(appended_info_species$taxa)){
  if (!is.na(dim(appended_info_species$coords[[i]])[3]) & dim(appended_info_species$coords[[i]])[3]>2) {
    A$PCvalues[[appended_info_species$taxa[i]]] <- appended_info_species$PCvalues[[i]]
    A$CSize[[appended_info_species$taxa[i]]] <- appended_info_species$CSize[[i]]
    A$logCsize[[appended_info_species$taxa[i]]] <- appended_info_species$logCsize[[i]]
    A$coords[[appended_info_species$taxa[i]]] <- as.array(appended_info_species$coords[[i]])
    A$species[[appended_info_species$taxa[i]]] <- appended_info_species$species[[i]]
  }
}
A$taxa <- names(A$CSize)

#######################################################################################################################################

# Test whether species differ in ontogenetic trajectory: Pairwise Comparisons

# Iteratively generating Procrustes linear model for each species and save slope/intercept coefficients
TrajectoryList=list()
SlopesList=list()
InterceptsList=list()
ElevationsList=list()
for (j in 1:length(A$taxa)){
  x<-procD.lm(A$PCvalues[[j]]~A$logCsize[[j]],iter=999) # Using A$coords or A$PCvalues gives the same output for the pairwise comparisons, need to use A$PCvalues here to plot the trajectories
  TrajectoryList[[A$taxa[[j]]]]<-x
  SlopesList[[A$taxa[[j]]]]<-x$coefficients[2,]
  InterceptsList[[A$taxa[[j]]]]<-x$coefficients[1,]
  ElevationsList[[A$taxa[[j]]]]<- (x$coefficients[2,] * mean(log(A$CSize[[j]])) ) + x$coefficients[1,]
}


# Iteratively create two species pairs and perform Procrustes ANOVA and save output data from comparison
# Pairwise Procrustes ANOVAs of ontogenetic trajectories
PairwiseComparisons=list()
SlopeDifferences=list()
InterceptDifferences=list()
g=1
n=2
for (g in 1:length(A$taxa)){
  for (h in n:length(A$taxa)){
    if (g==length(A$taxa)){
      break
    }else
      
    SpeciesA <- A$taxa[g]
    SpeciesB <- A$taxa[h]
    toMatch <- c(SpeciesA, SpeciesB)
    filename <- paste(SpeciesA, SpeciesB, sep=" vs. ")
    
    two.species <- abind(A$coords[[SpeciesA]],A$coords[[SpeciesB]])
    two.species.csize <- c(A$CSize[[SpeciesA]],A$CSize[[SpeciesB]])
    two.species.species <- as.factor(c(paste(A$species[[SpeciesA]]),paste(A$species[[SpeciesB]])))
    two.species.slope.diff<-SlopesList[[g]]-SlopesList[[h]]
    names(two.species.slope.diff) <- paste("PC", c(1:164), " slope diff.", sep="")
    two.species.intercept.diff<-InterceptsList[[g]]-InterceptsList[[h]]
    names(two.species.intercept.diff) <- paste("PC", c(1:164), " int. diff.", sep="")
    
    Output <- procD.lm(two.species~two.species.csize, ~two.species.species, logsz = TRUE, iter = 9999)##Remember to add more iterations##
    # logsz above = TRUE which means this is done for logged CS data
    PairwiseComparisons[[filename]]<-Output
    SlopeDifferences[[filename]]<-two.species.slope.diff
    InterceptDifferences[[filename]]<-two.species.intercept.diff
    
  }
  n=n+1
}


# Create matrix and fill with values from pairwise comparisons
# Fills the matrix with all pairwise comparisons from PairwiseComparisons object using $aov.table - pairwise Procrustes ANOVA of ontogenetic trajectories
Pvalues = matrix (nrow = length(names(PairwiseComparisons)),
                  ncol = length(c(paste("AOV", names(PairwiseComparisons[[g]]$aov.table[1,]), sep = "_"), names(SlopeDifferences[[g]][1]), names(SlopeDifferences[[g]][2]), names(InterceptDifferences[[g]][1]), names(InterceptDifferences[[g]][2]))),
                  dimnames = list(c(names(PairwiseComparisons)), c(paste("AOV", names(PairwiseComparisons[[g]]$aov.table[2,]), sep = "_"), names(SlopeDifferences[[g]][1]), names(SlopeDifferences[[g]][2]), names(InterceptDifferences[[g]][1]), names(InterceptDifferences[[g]][2])))
)
for (i in 1:length(names(PairwiseComparisons))){
  Pvalues[i,] <- rbind(as.numeric(paste(c(PairwiseComparisons[[i]]$aov.table[1,], SlopeDifferences[[i]][1], SlopeDifferences[[i]][2], InterceptDifferences[[i]][1], InterceptDifferences[[i]][2]))))
}

write.csv(Pvalues, "Analysis/Suture_morphology/Ontogenetic_allometric_trajectories/Sagittal/Pairwise_comparison_results_sagittal.csv")

# Bonferroni corrected p-values to account for pairwise comparisons
# Create the empty matrix
CorrectedPvalues <- matrix(nrow = length(names(PairwiseComparisons)),
                           ncol = length(c("corrected AOV p.value")),
                           dimnames = list(c(names(PairwiseComparisons)), c("corrected_AOV_p.value"))
)
# Fill the matrix with the corrected p-values
CorrectedPvalues[,1] <- p.adjust(Pvalues[,"AOV_Pr(>F)"], method = "bonferroni", n = length(PairwiseComparisons))

write.csv(CorrectedPvalues, "Analysis/Suture_morphology/Ontogenetic_allometric_trajectories/Sagittal/Pairwise_comparisons_corrected_pvalues_sagittal.csv")

# Convert to a dataframe to add columns
Pvalues <- as.data.frame(Pvalues)
CorrectedPvalues <- as.data.frame(CorrectedPvalues) %>% mutate(AOV_Rsq = Pvalues$AOV_Rsq)

CorrectedPvalues_sig <- CorrectedPvalues %>% mutate(sig_p = ifelse(corrected_AOV_p.value < .05, T, F),
                                                  p_if_sig = ifelse(sig_p, corrected_AOV_p.value, NA),
                                                  Rsq_if_sig = ifelse(sig_p, AOV_Rsq, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

# Add the pairwise comparison names to seperate columns to be able to plot the heatmap
rownames <- as.character(rownames(CorrectedPvalues_sig))
rn <- do.call(rbind, strsplit(rownames, 'vs.'))
CorrectedPvalues_sig <- CorrectedPvalues_sig %>% mutate(Var1 = rn[,1], Var2 = rn[,2])
CorrectedPvalues_sig$Rsq_if_sig <- round(CorrectedPvalues_sig$Rsq_if_sig, digit =2)


# Create palette for heatmap plot
mypalette_seq <- brewer.pal(9,"Oranges")
image(1:9,1, as.matrix(1:9), col = mypalette_seq,xlab="Oranges (sequential)",
      ylab = "", yaxt = "n")

# Plot a heatmap for the pairwise comparison
Allometry_pairwise_heatmap <- ggplot(data = CorrectedPvalues_sig, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  #geom_text(aes(Var2, Var1, label = Rsq_if_sig), color = "white", size = 3) +
  scale_fill_gradient2(low = mypalette_seq[9], high = mypalette_seq[2], mid = mypalette_seq[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Ontogenetic allometry pairwise differences - sagittal suture")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 90, size = 14, vjust = 0, hjust = 1),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5))
Allometry_pairwise_heatmap


#######################################################################################################################################

# Plot the trajectory output
# log(CS) vs PC1 with regression line plotted for each species

# highlight all the lines below and run together
PC=1
#Plot PCs vs log(CS)
Xlim<-1.1*c(min(log(gdf$CSize)),max(log(gdf$CSize)))
Ylim<-1.1*c(min(PCA$x[,PC]),max(PCA$x[,PC]))

plot(0, 0, type = "n",
     xlim = c(2,max(Xlim)),
     ylim = c(-1.8,max(Ylim)),
     xlab = "log(Centroid Size)",
     ylab = xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")),
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

axis(1, c(2,3,4,5,6,7,8), pos=-1.8)
axis(2, c(-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4), pos=2)
clip(-1,max(Xlim),-1.8,max(Ylim))
abline(h=0, lty=3)

points(log(gdf$CSize), PCA$x[,PC], pch=21, bg=alpha(col.species, 0.75), asp=F)
#legend("right", legend = A$taxa, col = col.gp.2)

for (i in 1:length(TrajectoryList)){
  Line_coefficients <- TrajectoryList[[i]]$coefficients[,PC]
  Line_color <- mammal.PlottingValues$legendcolor[names(TrajectoryList)[i]]
  CS <- log(appended_info_species$CSize[[match(names(TrajectoryList)[i], appended_info_species$taxa)]])
  clip(min(CS),max(CS),-1.8,1)
  abline(Line_coefficients, col=Line_color, lwd=2)
}




