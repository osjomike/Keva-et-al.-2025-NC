rm(list=ls())

## Author: Ossi Keva
## Date 14.2.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description: With this script we are plotting covariate effect on consumer allochthony in different omega scenarios
## Article Extended data Figure 5

# Some packages to download#
library(dplyr) # data arrangement
library(ggplot2) # drawing of plots
library(cowplot) # get legend and plot arranging
library(gridExtra) # arranging of plots
library(MixSIAR) #

## Selecting working directory
setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()
list.files()
f<-"long_Mod with lake+species_rand+PC1 w_0.23_randSlope"
file.path(orgfolder,"3. MixSIAR Models", f, "MixSIAR_modelMod with lake+species_rand+PC1 w_0.23")
### uploading mixtures, sources, discrimination and jags object
mix<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, "mix_object.rds"))
source<-readRDS(file.path(orgfolder,"3. MixSIAR Models", f, "source_object.rds"))
discr<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, paste0("discr_object",".rds")))
jags<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, "mixSIAR_model_Mod with lake+species_rand+PC1 w_0.23_randSlope.rds"))

posterior_array<-jags[["BUGSoutput"]][["sims.list"]][["p.ind"]]
posterior_df<-reshape2::melt(posterior_array)
colnames(posterior_df)<-c("draw", "individual", "Source", "distribution")
posterior_df$Source2<-factor(posterior_df$Source, labels = source[["source_names"]])
posterior_df$Source2<-factor(posterior_df$Source2, levels=source[["source_names"]])

## Extracting the simulation individual names and lakes, multiplying them with n.sims
number.sim<-jags[["BUGSoutput"]][["n.sims"]]
number.sources<-length(source[["source_names"]])
gg<-cbind(rep(rep(mix[["data"]][["Lake"]], each=number.sim), times=number.sources),  # times=2 indicates number of sources
          rep(rep(mix[["data"]][["species"]], each=number.sim), times=number.sources), # times=2 indicates number of sources
          rep(rep(mix[["data"]][["PC1"]], each=number.sim), times=number.sources)) # times=2 indicates number of sources
colnames(gg)<-c("Lake", "species", "PC1")
posterior<-cbind(posterior_df, gg)

######
kldiv_matrix<-matrix(data=NA, nrow = length(unique(posterior$Lake)),ncol=length(unique(posterior$species)),
              dimnames = list(unique(posterior$Lake), unique(posterior$species)))
kldiv_up_matrix<-matrix(data=NA, nrow = length(unique(posterior$Lake)),ncol=length(unique(posterior$species)),
                        dimnames = list(unique(posterior$Lake), unique(posterior$species)))
kldiv_rel_matrix<-matrix(data=NA, nrow = length(unique(posterior$Lake)),ncol=length(unique(posterior$species)),
         dimnames = list(unique(posterior$Lake), unique(posterior$species)))

alpha <- rep(1, source$n.sources) #default prior values

Kullbeck<-function (x1, x2, nbreaks = 100, minx = min(c(x1, x2)), maxx = max(c(x1, x2)), small = 0.01) 
{
  dy <- (maxx - minx)/nbreaks
  breaks <- seq(minx, maxx, by = dy)
  x1d <- hist(x1, plot = F, breaks = breaks)
  x2d <- hist(x2, plot = F, breaks = breaks)
  x1p <- (x1d$counts + small)/sum(x1d$counts + small)
  x2p <- (x2d$counts + small)/sum(x2d$counts + small)
  vals <- x1p * log2(x1p/x2p)
  #up_vals<- (-x1p)*log2(x1p)
  kd <- sum(na.omit(vals))
  #up_bound_kd<-sum(na.omit(up_vals))
  #hout <- list(kd = kd, up_bound_kd=up_bound_kd, histograms = list(h1 = x1d, h2 = x2d))
  hout <- list(kd = kd,  histograms = list(h1 = x1d, h2 = x2d))
  class(hout) <- "kldiverg"
  return(hout)
}

Kullbeck_max<-function (x1, x2, nbreaks = 100, minx = min(c(x1, x2)), maxx = max(c(x1, x2)), small = 0.01) {
  dy <- (maxx - minx)/nbreaks
  breaks <- seq(minx, maxx, by = dy)
  x1d <- hist(x1, plot = F, breaks = breaks)
  x1p <- (x1d$counts + small)/sum(x1d$counts + small)
  up_vals<- (-x1p)*log2(x1p)
  up_bound_kd<-sum(na.omit(up_vals))
  hout <- list(up_bound_kd = up_bound_kd,  histograms = list(h1 = x1d))
  class(hout) <- "kldiverg"
  return(hout)
}

### R keeps on crashing, I ques to following for loop is not very efficient. 
## Lets just try to release memory from globenvironment, by removing couple of large items
jags<-NULL
mix<-NULL
source<-NULL
posterior_df<-NULL
p_post<-NULL
posterior_array<-NULL
gg<-NULL

## One could do this probably more efficiently with calculating the Kullbeck with arrays
for(i in rownames(kldiv_matrix)){
  for(j in colnames(kldiv_matrix)){

    p_post2 <- posterior[posterior$Lake==i &
                           posterior$species == j & 
                           posterior$Source2 =="Ter",]
    p_post_m<-as.matrix(p_post2)

    if(nrow(p_post_m)==0){
      kldiv_matrix[i,j]<-NA
    } ## End of if 
    if((nrow(p_post_m)!=0)){
    p_prior <- MCMCpack::rdirichlet(nrow(p_post_m), alpha) #draw prior samples
    kldiv_matrix[i,j]<-as.numeric(Kullbeck(as.numeric(p_post_m[,4]), p_prior[,1])$kd) ## x1=posterior, x2=prior
    }  ## End of if 
  } ## End of j loop
} ## End of i loop

for(i in rownames(kldiv_up_matrix)){
  for(j in colnames(kldiv_up_matrix)){
    
    p_post2 <- posterior[posterior$Lake==i &
                           posterior$species == j & 
                           posterior$Source2 =="Ter",]
    p_post_m<-as.matrix(p_post2)
    
    if(nrow(p_post_m)==0){
      kldiv_matrix[i,j]<-NA
    } ## End of if 
    if((nrow(p_post_m)!=0)){
      p_prior <- MCMCpack::rdirichlet(nrow(p_post_m), alpha) #draw prior samples
      kldiv_up_matrix[i,j]<-as.numeric(Kullbeck_max(as.numeric(p_post_m[,4]), p_prior[,1])$up_bound_kd)
    }  ## End of if 
  } ## End of j loop
} ## End of i loop


## transforming kldiv_matrix to dataframe
kldiv_df<-as.data.frame(kldiv_matrix)

colMeans(kldiv_df, na.rm=TRUE)
library(reshape2)
box<-melt(kldiv_df)
box$Lake<-rownames(kldiv_df)
colnames(box)[1]<-"species"


ZPL<-c("Bulk_ZPL", "Chaoborus", "Cladocera", "Copepods")
BMI<-c( "Bulk_BMI_littoral",
              "Asellus_littoral", "Chironomid_littoral",
               "Bulk_BMI_profundal", "Chironomid_Profundal")
Fish_plankt<-c("Perch_small", "Smelt", "Vendace", "Bleak")
Fish_mixed<-c("Perch_medium", "Roach_small", "Roach_large")
Fish_benthivor<-c("Ruffe")
Fish_piscivorous<-c("Perch_large", "Pike")

box$group<-NA ## making a new variable to box dataframe
box[box$species %in% ZPL,]$group<-"ZPL"
box[box$species %in% BMI,]$group<-"BMI"
box[box$species %in% Fish_plankt,]$group<-"Fish_plankt"
box[box$species %in% Fish_mixed,]$group<-"Fish_mixed"
box[box$species %in% Fish_benthivor,]$group<-"Fish_benthivor"
box[box$species %in% Fish_piscivorous,]$group<-"Fish_piscivorous"

box$species<-factor(box$species, levels = c(ZPL, BMI, Fish_benthivor, Fish_mixed, Fish_plankt, Fish_piscivorous))

ggplot(box)+
  geom_boxplot(data=box, aes(x=species, y=value, fill=group))+
  labs(y="Marginal Kullnack-Leibler \n divergence (gained information in bits)", x=NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
write.csv(kldiv_matrix, file="PC1_long_model_Kldivmatrix.csv")
