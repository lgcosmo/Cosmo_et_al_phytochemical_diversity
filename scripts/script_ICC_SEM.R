#Calculates ICC
#Performs SEM analysis

library(rptR)
library(piecewiseSEM)
library(lme4)
library(nlme)
library(here)

dados<-read.table(here("data", "Cosmo_et_al_fulldataset.txt"), header=TRUE) #Importing dataset
dados$logitherb<-dados$herbivory/(1-dados$herbivory) #Creating logit variable of herbivory

#--------------------------------------------------------------------#
#-Using the package rptR to calculate ICC for all response variables-#
#--------------------------------------------------------------------#

icc_abundance<-rptPoisson(formula=catabund~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_richness<-rptPoisson(formula=richness~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_swindex<-rptGaussian(formula=sw_index~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_herbivory<-rptGaussian(formula=logitherb~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_structural_pd<-rptGaussian(formula=structural_pd~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_compPC1<-rptGaussian(formula=comp_PC1~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)
icc_compPC2<-rptGaussian(formula=comp_PC2~(1|plant), grname=c("plant"), data=dados, nboot=10000, parallel=TRUE)

#-------------------------------------------------#
#-Using the package piecewiseSEM to fit SEM model-#
#-------------------------------------------------#

sem<-psem(
  
  lme(scale(sw_index)~catabund+richness+scale(within_structural_pd)+scale(among_structural_pd)+scale(within_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC1)+scale(among_comp_PC2)+stratum+scale(nleaves)+scale(plant_height_cm)+season, random=~1|plant, data=dados),
  glmer(richness~catabund+scale(within_structural_pd)+scale(among_structural_pd)+scale(within_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC1)+scale(among_comp_PC2)+scale(plant_height_cm)+stratum+season+scale(nleaves)+(1|plant), family="poisson", data=dados),
  glmer(catabund~scale(within_structural_pd)+scale(among_structural_pd)+scale(within_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC1)+scale(among_comp_PC2)+scale(plant_height_cm)+stratum+season+scale(nleaves)+(1|plant), family="poisson", data=dados),
  lme(logitherb~scale(sw_index)+scale(within_structural_pd)+scale(among_structural_pd)+scale(within_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC1)+scale(among_comp_PC2)+stratum+season+scale(nleaves)+scale(plant_height_cm), random=~1|plant, data=dados),
  lme(scale(within_comp_PC1)~stratum+season, random=~1|plant, data=dados),
  lme(scale(within_comp_PC2)~stratum+season, random=~1|plant, data=dados),
  lm(scale(among_comp_PC1)~season+scale(plant_height_cm), data=dados),
  lm(scale(among_comp_PC2)~season+scale(plant_height_cm), data=dados),
  lme(scale(within_structural_pd)~stratum+season, random=~1|plant, data=dados),
  lm(scale(among_structural_pd)~season+scale(plant_height_cm), data=dados),
  scale(within_structural_pd) %~~% scale(within_comp_PC1),
  scale(within_structural_pd) %~~% scale(within_comp_PC2),
  scale(among_structural_pd) %~~% scale(among_comp_PC1),
  scale(among_structural_pd) %~~% scale(among_comp_PC2),
  
  data=dados
  
  
)

summary(sem, standardize="scale", direction=c(("scale(sw_index)<-catabund"), ("scale(sw_index)<-richness"), ("richness<-catabund"))) ##Coefficients and model fit
