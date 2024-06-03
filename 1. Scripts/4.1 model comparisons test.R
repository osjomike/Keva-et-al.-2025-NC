rm(list=ls())

## Author: Ossi Keva
## Date 03.06.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description: with this script we are comparing MixSIAR models with different covariates
## Suppelementary table 8

# Some packages to download
library(MixSIAR)
library(ggplot2)

## Selecting working directory
setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()

folders<-c("very long_Mod with lake+species_rand+PC1 w_0.23_randSlope",
           "very long_Mod with lake+species_rand+For_area w_0.23_randSlope",
           "very long_Mod with lake+species_rand+POM_NC_ratio w_0.23_randSlope",
           "long_Mod with lake+species w_0.23",
           "long_Mod with lake w_0.23",
           "long_empty")

jags.objects<-c("mixSIAR_model_Mod with lake+species_rand+PC1 w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species_rand+For_area w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species_rand+POM_NC_ratio w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species w_0.23.rds",
                "mixSIAR_model_Mod with lake w_0.23.rds",
                "mixSIAR_model_empty.rds")

model_list<-list()

for(i in 1:length(jags.objects)){
  setwd(file.path(orgfolder, "3. MixSIAR Models", folders[i]))
  model_list[[i]]<-readRDS(file=jags.objects[i])
  setwd(orgfolder)
}

names(model_list)<-c("PC1+(1+PC1|species)+(1|species)+(1|Lake)",
                     "For%+(1+For%|species)+(1|species)+(1|Lake)",
                     "(POM_N:C)+ (1+POM_N:C|species)+(1|species)+(1|Lake)",
                     "(1|species)+(1|Lake)",
                     "(1|species)",
                     "1")
compare_models(model_list, loo=TRUE)
#warnings()
#getAnywhere(compare_models)


### OK, now lets observe the covariate slope distribution in ILR space and how it affects the mixtures 

# making toy data first column is aq and second ter, This is proportional data 
# mixSIAR load_source_data() -function organises the sources alphabetically, 
# so this toy data emulates the situation we are presenting in the paper
test_data<-matrix(c(c(0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999),
                    c(0.999, 0.99, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0.001)), ncol = 2)
## Making ILR conversion of the the toy data to check how does the transformation impact the proportions
ILR_2<-c()
for ( i in 1:nrow(test_data)){
  ILR_2[i]<-sqrt(1/(1+1))*log(test_data[i,1]/test_data[i,2])
}

## We can observe that higher ILR scores indicate higher contribution of test_data[1] values
## In our case (manuscript) these higher ILR scores indicate higher aquatic dietary proportion in consumers
## and lower terrestrial contribution. In the other hand low ILR scores indicate high terrestrial contribution. 

## Then lets check what is the slope distribution in our model
PC1<-model_list[[1]]
CIs<-quantile(PC1[["BUGSoutput"]][["sims.list"]][["ilr.cont1"]], probs = c(0.025, 0.5, 0.975))
scaleFUN <- function(x) sprintf("%.3f", x)

library(ggplot2)
dat<-data.frame(slope=PC1[["BUGSoutput"]][["sims.list"]][["ilr.cont1"]])
distribution<-ggplot(dat, aes(x=slope))+ 
  geom_histogram(aes(y=..density..), binwidth = 0.01, colour="black", fill="white")+
  geom_density(fill="#FF6666", alpha=0.5)+ 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey")+
  coord_cartesian(xlim = c(-1,1))+
  annotate(geom="text", label=paste0("median (95% CI) = ", scaleFUN(CIs[2]), 
                                      " (", scaleFUN(CIs[1]), "--", scaleFUN(CIs[3]), ")"),
           x=0.5, y=4.5, size=4)+
  labs(y="Scaled density", x="Linear coefficient with PC1 in ILR space")+
  theme_classic()

png(file = "4. Figures/Linear coefficient with PC1 in ILR space.png",
    width = 7, height= 3, res=300, units="in")
distribution
dev.off()

### Now when comparing the ILR conversion and the slope distribution 
### we can conclude that increasing PCA increases aquatic source contribution in the consumer dietary 
### as higher ILR scores indicate higher aquatic contribution 
### and because of the coeffiecient have a positive slope distribution

#########################
### END OF THE SCRIPT ###
#########################