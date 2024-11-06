rm(list=ls())
## Author: Matthew Cobain & Ossi Keva
## Date 6.11.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description:
# Code to calculate a pseudo- conditional r2 for the models
# The relative proportion in isotope variation explained by the model
# Note though that mixing models are not constructed for optimal isotope reconstruction but to quantify the latent variable of diet proportions

### Preamble
library(MixSIAR)
library(coda)

WD <- "D:\\Keva et al. NC 2024/3. MixSIAR Models/"
folders<-c("very long_Mod with lake+species_rand+PC1 w_0.23_randSlope/",
           "very long_Mod with lake+species_rand+For_area w_0.23_randSlope/",
           "very long_Mod with lake+species_rand+POM_NC_ratio w_0.23_randSlope/",
           "long_Mod with lake+species w_0.23/",
           "long_Mod with lake w_0.23/",
           "long_empty/")

jags.objects<-c("mixSIAR_model_Mod with lake+species_rand+PC1 w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species_rand+For_area w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species_rand+POM_NC_ratio w_0.23_randSlope.rds",
                "mixSIAR_model_Mod with lake+species w_0.23.rds",
                "mixSIAR_model_Mod with lake w_0.23.rds",
                "mixSIAR_model_empty.rds")

Model_r2s <- data.frame(matrix(NA, nrow = length(jags.objects), ncol = 5))
colnames(Model_r2s) <- c("Model", "Median", "Mean", "95CI_L", "95CI_U")

### PC1 Model ###

#####
# Set working directory
setwd(paste0(WD, folders[1]))

# Load model, mix and source MixSIAR objects 
model <- readRDS(jags.objects[1])
mix <- readRDS("mix_object.rds")
source <- readRDS("source_object.rds")

# Following Stock et al. 2016, we predicted consumer i to have isotope composition X_j:
# X_ij = N(Sum_k{P_k.(Mu_jk + Lamda_jk)}, Sum_k{P^2_k.(W^2_jk + Tau^2_jk).Epsilon_j})
# Therefore predict based on means only (do not need to incorporate variances)

# Extract individual proportions from the model
P.ind <- model$BUGSoutput$sims.list$p.ind

# Extract original observations
mix.data <- mix$data
mix.data$Lake <- factor(mix.data$Lake, levels = mix[["FAC"]][[2]][["labels"]])
mix.data$species <- factor(mix.data$species, levels = mix[["FAC"]][[1]][["labels"]])

# Extract source means (by lake)
Source_means <- source[["MU_array"]]

# Construct source means array with source values for each individual consumer
Lakes.ind <- as.numeric(mix.data$Lake)
Source_means.ind <- t(Source_means[,,Lakes.ind])

# Multiply source means by proportions and sum for each individual applied across the 3000 draws
Pred.ind.Iso <- apply(P.ind, 1, FUN = function(x) {rowSums(x*Source_means.ind)})

### Estimating a pseudo-r2

# r2 = 1 - SS_res / SS_tot

# Estimate total sum of squares (observed - mean)
SS_tot = sum((mix$data_iso - mean(mix$data_iso))^2)

# Estimate residual sum of squares (observed - predicted)
SS_res = apply(Pred.ind.Iso, 2, FUN = function(x) {sum((mix$data_iso - x)^2)})

# r2
model.r2 <- 1 - SS_res/SS_tot 

# Store outputs
Model_r2s[1,1] <- jags.objects[1]
Model_r2s[1,2] <- signif(median(model.r2),3)
Model_r2s[1,3] <- signif(mean(model.r2),3)
Model_r2s[1,4] <- signif(HPDinterval(mcmc(model.r2))[1],3)
Model_r2s[1,5] <- signif(HPDinterval(mcmc(model.r2))[2],3)

# Clear environment
rm(list=setdiff(ls(), c("Model_r2s", "WD", "folders", "jags.objects")))

#####

### Forest Area Model ###

#####
# Set working directory
setwd(paste0(WD, folders[2]))

# Load model, mix and source MixSIAR objects 
model <- readRDS(jags.objects[2])
mix <- readRDS("mix_object.rds")
source <- readRDS("source_object.rds")

# Following Stock et al. 2016, we predicted consumer i to have isotope composition X_j:
# X_ij = N(Sum_k{P_k.(Mu_jk + Lamda_jk)}, Sum_k{P^2_k.(W^2_jk + Tau^2_jk).Epsilon_j})
# Therefore predict based on means only (do not need to incorporate variances)

# Extract individual proportions from the model
P.ind <- model$BUGSoutput$sims.list$p.ind

# Extract original observations
mix.data <- mix$data
mix.data$Lake <- factor(mix.data$Lake, levels = mix[["FAC"]][[2]][["labels"]])
mix.data$species <- factor(mix.data$species, levels = mix[["FAC"]][[1]][["labels"]])

# Extract source means (by lake)
Source_means <- source[["MU_array"]]

# Construct source means array with source values for each individual consumer
Lakes.ind <- as.numeric(mix.data$Lake)
Source_means.ind <- t(Source_means[,,Lakes.ind])

# Multiply source means by proportions and sum for each individual applied across the 3000 draws
Pred.ind.Iso <- apply(P.ind, 1, FUN = function(x) {rowSums(x*Source_means.ind)})

### Estimating a pseudo-r2

# r2 = 1 - SS_res / SS_tot

# Estimate total sum of squares (observed - mean)
SS_tot = sum((mix$data_iso - mean(mix$data_iso))^2)

# Estimate residual sum of squares (observed - predicted)
SS_res = apply(Pred.ind.Iso, 2, FUN = function(x) {sum((mix$data_iso - x)^2)})

# r2
model.r2 <- 1 - SS_res/SS_tot 

# Store outputs
Model_r2s[2,1] <- jags.objects[2]
Model_r2s[2,2] <- signif(median(model.r2),3)
Model_r2s[2,3] <- signif(mean(model.r2),3)
Model_r2s[2,4] <- signif(HPDinterval(mcmc(model.r2))[1],3)
Model_r2s[2,5] <- signif(HPDinterval(mcmc(model.r2))[2],3)

# Clear environment
rm(list=setdiff(ls(), c("Model_r2s", "WD", "folders", "jags.objects")))

#####

### POM NC Model ###

#####
# Set working directory
setwd(paste0(WD, folders[3]))

# Load model, mix and source MixSIAR objects 
model <- readRDS(jags.objects[3])
mix <- readRDS("mix_object.rds")
source <- readRDS("source_object.rds")

# Following Stock et al. 2016, we predicted consumer i to have isotope composition X_j:
# X_ij = N(Sum_k{P_k.(Mu_jk + Lamda_jk)}, Sum_k{P^2_k.(W^2_jk + Tau^2_jk).Epsilon_j})
# Therefore predict based on means only (do not need to incorporate variances)

# Extract individual proportions from the model
P.ind <- model$BUGSoutput$sims.list$p.ind

# Extract original observations
mix.data <- mix$data
mix.data$Lake <- factor(mix.data$Lake, levels = mix[["FAC"]][[2]][["labels"]])
mix.data$species <- factor(mix.data$species, levels = mix[["FAC"]][[1]][["labels"]])

# Extract source means (by lake)
Source_means <- source[["MU_array"]]

# Construct source means array with source values for each individual consumer
Lakes.ind <- as.numeric(mix.data$Lake)
Source_means.ind <- t(Source_means[,,Lakes.ind])

# Multiply source means by proportions and sum for each individual applied across the 3000 draws
Pred.ind.Iso <- apply(P.ind, 1, FUN = function(x) {rowSums(x*Source_means.ind)})

### Estimating a pseudo-r2

# r2 = 1 - SS_res / SS_tot

# Estimate total sum of squares (observed - mean)
SS_tot = sum((mix$data_iso - mean(mix$data_iso))^2)

# Estimate residual sum of squares (observed - predicted)
SS_res = apply(Pred.ind.Iso, 2, FUN = function(x) {sum((mix$data_iso - x)^2)})

# r2
model.r2 <- 1 - SS_res/SS_tot 

# Store outputs
Model_r2s[3,1] <- jags.objects[3]
Model_r2s[3,2] <- signif(median(model.r2),3)
Model_r2s[3,3] <- signif(mean(model.r2),3)
Model_r2s[3,4] <- signif(HPDinterval(mcmc(model.r2))[1],3)
Model_r2s[3,5] <- signif(HPDinterval(mcmc(model.r2))[2],3)

# Clear environment
rm(list=setdiff(ls(), c("Model_r2s", "WD", "folders", "jags.objects")))

#####

### Random Lake / Spp. Model ###

#####
# Set working directory
setwd(paste0(WD, folders[4]))

# Load model, mix and source MixSIAR objects 
model <- readRDS(jags.objects[4])
mix <- readRDS("mix_object.rds")
source <- readRDS("source_object.rds")

# Extract original observations
mix.data <- mix$data
mix.data$Lake <- factor(mix.data$Lake, levels = mix[["FAC"]][[2]][["labels"]])
mix.data$species <- factor(mix.data$species, levels = mix[["FAC"]][[1]][["labels"]])

# Extract source means (by lake)
Source_means <- source[["MU_array"]]

# Construct source means array with source values for each individual consumer
Lakes.ind <- as.numeric(mix.data$Lake)
Source_means.ind <- t(Source_means[,,Lakes.ind])


# Extract individual proportions from the model
# This needs to be done manually for this model
chain.len <- model$BUGSoutput$n.sims
n.sources <- source$n.sources

# Extract the isometric log ratios
# NB: for some reason, MixSIAR does not store ilr.global for random effects only models! So need to back calculaute from p.global
ilr.globalX <- sqrt(1/(1+1))*log(model$BUGSoutput$sims.list$p.global[,1]/model$BUGSoutput$sims.list$p.global[,2])
ilr.fac1X <- model$BUGSoutput$sims.list$ilr.fac1
ilr.fac2X <- model$BUGSoutput$sims.list$ilr.fac2 

ilr.ind <- ilr.globalX + ilr.fac1X[,as.numeric(mix.data$species),] + ilr.fac2X[,as.numeric(mix.data$Lake),]

# Transform ilr.fac2_tot from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
for(i in 1:(n.sources-1)){
  e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
  e[,i] <- e[,i]/sum(e[,i])
}
dim(ilr.ind)[2]
# dummy variables for inverse ILR calculation
cross <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2], n.sources))  
tmp <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2] ,n.sources))  
P.ind <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2], n.sources))  

for(d in 1:chain.len){
  for(l in 1:dim(ilr.ind)[2]){
    for(j in 1:(n.sources-1)){
      cross[d,l,] <- (e[,j]^ilr.ind[d,l])/sum(e[,j]^ilr.ind[d,l])
    }
    for(src in 1:n.sources){
      tmp[d,l,src] <- prod(cross[d,l,src])
    }
    for(src in 1:n.sources){
      P.ind[d,l,src] <- tmp[d,l,src]/sum(tmp[d,l,])
    }
  }
}

# Multiply source means by proportions and sum for each individual applied across the 3000 draws
Pred.ind.Iso <- apply(P.ind, 1, FUN = function(x) {rowSums(x*Source_means.ind)})

### Estimating a pseudo-r2

# r2 = 1 - SS_res / SS_tot

# Estimate total sum of squares (observed - mean)
SS_tot = sum((mix$data_iso - mean(mix$data_iso))^2)

# Estimate residual sum of squares (observed - predicted)
SS_res = apply(Pred.ind.Iso, 2, FUN = function(x) {sum((mix$data_iso - x)^2)})

# r2
model.r2 <- 1 - SS_res/SS_tot 

# Store outputs
Model_r2s[4,1] <- jags.objects[4]
Model_r2s[4,2] <- signif(median(model.r2),3)
Model_r2s[4,3] <- signif(mean(model.r2),3)
Model_r2s[4,4] <- signif(HPDinterval(mcmc(model.r2))[1],3)
Model_r2s[4,5] <- signif(HPDinterval(mcmc(model.r2))[2],3)

# Clear environment
rm(list=setdiff(ls(), c("Model_r2s", "WD", "folders", "jags.objects")))

#####

### Random Lake only Model ###

#####
# Set working directory
setwd(paste0(WD, folders[5]))

# Load model, mix and source MixSIAR objects 
model <- readRDS(jags.objects[5])
mix <- readRDS("mix_object.rds")
source <- readRDS("source_object.rds")

# Extract original observations
mix.data <- mix$data
mix.data$Lake <- factor(mix.data$Lake, levels = mix[["FAC"]][[1]][["labels"]]) # Note change of factor number becasue no species

# Extract source means (by lake)
Source_means <- source[["MU_array"]]

# Construct source means array with source values for each individual consumer
Source_means.ind <- t(Source_means[,,as.numeric(mix.data$Lake)])

# Extract individual proportions from the model
# This needs to be done manually for this model
chain.len <- model$BUGSoutput$n.sims
n.sources <- source$n.sources

# Extract the isometric log ratios
# NB: for some reason, MixSIAR does not store ilr.global for random effects only models! So need to back calculaute from p.global
ilr.globalX <- sqrt(1/(1+1))*log(model$BUGSoutput$sims.list$p.global[,1]/model$BUGSoutput$sims.list$p.global[,2])
ilr.fac1X <- model$BUGSoutput$sims.list$ilr.fac1
ilr.ind <- ilr.globalX + ilr.fac1X[,as.numeric(mix.data$Lake),]

# Transform ilr.fac2_tot from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
for(i in 1:(n.sources-1)){
  e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
  e[,i] <- e[,i]/sum(e[,i])
}

# dummy variables for inverse ILR calculation
cross <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2], n.sources))  
tmp <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2] ,n.sources))  
P.ind <- array(data=NA, dim=c(chain.len, dim(ilr.ind)[2], n.sources))  

for(d in 1:chain.len){
  for(l in 1:dim(ilr.ind)[2]){
    for(j in 1:(n.sources-1)){
      cross[d,l,] <- (e[,j]^ilr.ind[d,l])/sum(e[,j]^ilr.ind[d,l])
    }
    for(src in 1:n.sources){
      tmp[d,l,src] <- prod(cross[d,l,src])
    }
    for(src in 1:n.sources){
      P.ind[d,l,src] <- tmp[d,l,src]/sum(tmp[d,l,])
    }
  }
}

# Multiply source means by proportions and sum for each individual applied across the 3000 draws
Pred.ind.Iso <- apply(P.ind, 1, FUN = function(x) {rowSums(x*Source_means.ind)})

### Estimating a pseudo-r2

# r2 = 1 - SS_res / SS_tot

# Estimate total sum of squares (observed - mean)
SS_tot = sum((mix$data_iso - mean(mix$data_iso))^2)

# Estimate residual sum of squares (observed - predicted)
SS_res = apply(Pred.ind.Iso, 2, FUN = function(x) {sum((mix$data_iso - x)^2)})

# r2
model.r2 <- 1 - SS_res/SS_tot 

# Store outputs
Model_r2s[5,1] <- jags.objects[5]
Model_r2s[5,2] <- signif(median(model.r2),3)
Model_r2s[5,3] <- signif(mean(model.r2),3)
Model_r2s[5,4] <- signif(HPDinterval(mcmc(model.r2))[1],3)
Model_r2s[5,5] <- signif(HPDinterval(mcmc(model.r2))[2],3)

# Clear environment
rm(list=setdiff(ls(), c("Model_r2s", "WD", "folders", "jags.objects")))

#####

### Save Model R2 outputs
setwd(WD)
write.table(Model_r2s, file = "Model_R2_Comparisons.txt", row.names = FALSE)

