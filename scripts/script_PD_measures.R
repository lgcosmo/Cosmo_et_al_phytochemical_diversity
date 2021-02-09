#Performs PCA analysis of HPLC-MS processed data
#Extracts scores of PC1 and PC2 for comp-PC1 and comp-PC2 measures
#Calculates structural phytochemical diversity (PD) from NMR processed data
#Peforms within-subject centering to decompose comp-PC1, comp-PC2 and structural PD into within and between subjects components

#---------------------#
#-PCA of HPLC-MS data-#
#---------------------#

library(compositions)
library(speaq)
library(here)

ms<-read.table(here("data", "dataset_hplc_ms.txt"), header=TRUE) #Importing dataset

ms_names<-ms[, 1:2] #Separating sample and plant names
ms_data<-ms[, 3:18] #Separating matrix of compounds/intensity values

ms_data<-sweep(ms_data, MARGIN = 1, FUN = "/", STATS = rowSums(ms_data)) #Normalizing each row into proportions

ms_data<-as.data.frame(compositions::clr(ms_data+(min(ms_data[ms_data != 0])/10))) #Performing center log ratio transformation of proportions following 

ms_data<-speaq::SCANT(data.matrix = ms_data, type = c("pareto", "center")) #Scaling and centering before PCA

pca<-prcomp(ms_data) #Performing PCA analysis

scores<-as.data.frame(pca$x) #Retrieving PCA scores
comp_PC<-cbind(ms_names, scores[,1:2]) #FUll data frame with samples and plants with comp-PC1 and comp-PC2
names(comp_PC)[3:4]<-c("comp_PC1", "comp_PC2")

#--------------------------------------------------------------#
#-Calculating structural phytochemical diversity from NMR data-#
#--------------------------------------------------------------#

library(vegan)

nmr<-read.table("dataset_nmr.txt", header=TRUE) #Importing dataset
nmr_names<-nmr[,1:2] #Separating sample and plant names
nmr_data<-nmr[,3:403] #Separating matrix of chemical shifts/intensity values

structural_pd<-diversity(nmr_data, index="shannon") #Calculating SW-index for structural PD
structural_pd<-exp(structural_pd) #Calculating Hill's numbers of SW-index and obtaining final structural PD measure
structural_pd<-cbind(nmr_names, structural_pd)

#-------------------------------------------------------------------#
#-Performing within-subject centering of comp-PCs and structural PD-#
#-------------------------------------------------------------------#

#Compositional dimension of PD

among_comp_PCs <- data.frame(lapply(comp_PC[,3:4], function(x) 
  ave(x, comp_PC[,c("plant")], FUN = mean))) #Calculating subject (plants) means for comp-PCs

within_comp_PCs <- comp_PC[,3:4] - among_comp_PCs #Subtracting subject(plants) means from within-subject samples for comp-PCs

names(among_comp_PCs)<-c("among_comp_PC1", "among_comp_PC2") #Renaming dataframe columns
names(within_comp_PCs)<-c("within_comp_PC1", "within_comp_PC2") #Renaming dataframe columns

comp_PC<-cbind(comp_PC, within_comp_PCs, among_comp_PCs) #Merging all measures and obtaining final data frame for compositional dimension variables

#Structural dimension of PD

among_structural_PD<-ave(structural_pd[, 3], structural_pd[,c("plant")], FUN = mean) #Calculating subject (plants) means for structural PD
within_structural_PD<-structural_pd[, 3] - among_structural_PD #Subtracting subject(plants) means from within-subject samples for structural PD

structural_pd$among_structural_pd<-among_structural_PD #Putting among-plant measures into structural pd dataframe
structural_pd$within_structural_pd<-within_structural_PD #Putting within-plant measures into structural pd dataframe

PD_measures<-cbind(structural_pd, comp_PC[,3:8]) #Full dataframe with all PD measures
