rm(list=ls())

## Author: Ossi Keva
## Date 03.06.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description:

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

### And upload some data
f<-c("very long_Mod with lake+species_rand+PC1 w_0.23_randSlope")
file.path(orgfolder,"3. MixSIAR Models", f, "MixSIAR_modelMod with lake+species_rand+PC1 w_0.23")
### uploading mixtures, sources, discrimination and jags object
mix<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, "mix_object.rds"))
source<-readRDS(file.path(orgfolder,"3. MixSIAR Models", f, "source_object.rds"))
discr<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, paste0("discr_object",".rds")))
jags<-readRDS(file = file.path(orgfolder,"3. MixSIAR Models", f, "mixSIAR_model_Mod with lake+species_rand+PC1 w_0.23_randSlope.rds"))


############# Then lets do prectic the consumer tracer values ########
mixtures <- jags[["BUGSoutput"]][["sims.list"]][["p.ind"]]
source_means <- source[["MU_array"]]
source_variances <- source[["SIG2_array"]]  # Variances provided as SIG2_array

# Check dimensions to ensure alignment
print(dim(source_means))    # Expected: [sources, tracers, lakes]
print(dim(source_variances)) # Expected: [sources, tracers, lakes]
print(dim(mixtures))        # Expected: [iterations, consumers, sources]

# Extract lake-specific data
lakes_p <- mix[["FAC"]][[2]][["values"]]

# Initialize lists to store predictions
pred_tracer <- list()
pred_tracer_mean <- list()

# Loop through consumer to calculate predicted tracer values
for (i in 1:length(lakes_p)) {
  lake_index <- lakes_p[i]
  
  # Initialize matrix to store predicted values for each iteration
  pred_values <- matrix(0, nrow = dim(mixtures)[1], ncol = dim(source_means)[2])
  
  for (tracer in 1:dim(source_means)[2]) {
    for (source in 1:dim(source_means)[1]) {
      # Simulate source values using the posterior means and standard deviations
      simulated_source_values <- rnorm(n = dim(mixtures)[1],
                                       mean = source_means[source, tracer, lake_index],
                                       sd = sqrt(source_variances[source, tracer, lake_index]))
      
      # Calculate the predicted values for the current tracer by weighting the simulated source values
      pred_values[, tracer] <- pred_values[, tracer] + mixtures[, i, source] * simulated_source_values
      ### Question for MAT!  Is this enough or have I missed something?
    }
  }
  
  # Store the predicted values for the current consumer and store the individual index as well
  pred_tracer[[i]] <- cbind(pred_values, i)
  
  ## This is probably not needed
  # Calculate the mean of the predicted values for each tracer across iterations
  pred_tracer_mean[[i]] <- cbind(colMeans(pred_values), i)
}


### First idea is just look how the mean predicted values correlate with the observed isotope values
predicted_d2h_mean<-do.call("rbind", pred_tracer_mean)
data_xx<-mix[["data"]]
data_xx$predicted_d2H<-predicted_d2h_mean[,1]

colnames(data_xx)
fig_list2<-list()
for(i in unique(data_xx$species)){
  data_sub<-data_xx[data_xx$species==i,]
  linear_model<-lm(d2H~predicted_d2H, data=data_sub)
  mod_sum<-summary(linear_model)

  fig_list2[[i]]<-ggplot(data=data_sub)+
    geom_point(aes(x=predicted_d2H, y=d2H))+
    geom_smooth(method = lm, aes(x=predicted_d2H, y=d2H))+
    annotate("text", label=paste0("adj.R2 = ", round(mod_sum$adj.r.squared, digits=3)), 
             x=mean(data_sub$predicted_d2H), y=mean(data_sub$d2H))+
    labs(title=i)
  
}
## THIS DOES NOT LOOK PROMISING


### Second idea is to plot the prediction distributions and plot observed d2H values on top of them
predicted_d2h<-do.call("rbind", pred_tracer)
#### Then lets create equally long columns including species and lake columns and cbind them with predicted_d2h
number.sim=3000
gg<-cbind(rep(mix[["data"]][["Lake"]], each=number.sim),  
          rep(mix[["data"]][["species"]], each=number.sim)) 
colnames(gg)<-c("Lake", "species")

predicted_d2h<-cbind(predicted_d2h, gg)

### Chance the column names and making it a dataframe
colnames(predicted_d2h)<-c("pred_d2H", "individual","Lake", "species")
predicted_d2h_df<-as.data.frame(predicted_d2h)
## Just ensuring that the predicted values are numeric
predicted_d2h_df$pred_d2H<-as.numeric(predicted_d2h_df$pred_d2H)

## Summary of the predicted tracer values
predicted_d2h_sum<-predicted_d2h_df %>%
  group_by(Lake, species) %>%
  summarise(y0.025 = quantile(pred_d2H, 0.025, na.rm = TRUE),
            y0.1= quantile(pred_d2H, 0.1, na.rm = TRUE),
            y0.25= quantile(pred_d2H, 0.25, na.rm = TRUE),
            y0.375= quantile(pred_d2H, 0.375, na.rm = TRUE),
            y0.4= quantile(pred_d2H, 0.4, na.rm = TRUE),
            y0.5= quantile(pred_d2H, 0.5, na.rm = TRUE),
            y0.625= quantile(pred_d2H, 0.625, na.rm = TRUE),
            y0.75= quantile(pred_d2H, 0.75, na.rm = TRUE),
            y0.8= quantile(pred_d2H, 0.8, na.rm = TRUE),
            y0.9= quantile(pred_d2H, 0.9, na.rm = TRUE),
            y0.975= quantile(pred_d2H,0.975, na.rm = TRUE))
predicted_d2h_sum<-as.data.frame(predicted_d2h_sum)

## Lets create a list for the figures
fig_list3<-list()
## And make a for loop to store the figures
for (i in unique(predicted_d2h_sum$species)){
  data_sub<-predicted_d2h_sum[predicted_d2h_sum$species==i,]
  fig_list3[[i]]<-ggplot(data=data_sub)+
    geom_linerange(aes(x=Lake, ymin=y0.025, ymax=y0.975))+
    geom_point(aes(x=Lake, y=y0.5))+
    geom_point(inherit.aes = FALSE, data=mix[["data"]][mix[["data"]]$species==i,], aes(x=Lake, y=d2H), shape=8, col="blue")+
    labs(title=i, y="d2H")+
    theme(axis.text.x = element_text(angle = 90))
  
}

## This looks quite ok, but I think there could be errors in predicting the consumer tracer values
## For example not sure if just resampling of the source values is enough?? Also I think we should ad the TDF SD to the sources
## TDF SD is 13


### Then we continue with the Kullbeck-leibeir statistics

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

Kullbeck<-function (x1, x2, nbreaks = 100, minx = min(c(x1, x2)), maxx = max(c(x1, x2)), small = 0.01) {
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


############ THEN SOME ALTERNATIVE STUFF: MC FADDEN PSEUDO R2 #### NOT SURE IF THIS IS VALID AT ALL

setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()
list.files()
full<-"test_Mod with lake+species_rand+PC1 w_0.23_randSlope"
full<-"long_Mod with lake+species_rand+PC1 w_0.23_randSlope"
empty<-"test_empty"

### uploading mixtures, sources, discrimination and jags object
jags_empty<-readRDS(file = file.path(orgfolder,empty, "mixSIAR_model_empty.rds"))
jags_full<-readRDS(file = file.path(orgfolder,full, "mixSIAR_model_Mod with lake+species_rand+PC1 w_0.23_randSlope.rds"))
mix<-readRDS(file=file.path(orgfolder, empty, "mix_object.rds"))

#### McFadden Pseudo R2 value for each individual
loglik_EMPTY<-jags_empty[["BUGSoutput"]][["sims.list"]][["loglik"]]
loglik_FULL<-jags_full[["BUGSoutput"]][["sims.list"]][["loglik"]]
dim(loglik_EMPTY)
dim(loglik_FULL)

### McFadden, D. (1987). Regression-based specification tests for the multinomial logit model. Journal of econometrics, 34(1-2), 63-82.
PseudoR<-list()
for (i in 1:1737){
  PseudoR[[i]]<-1-(loglik_FULL[1:1500,i]/loglik_EMPTY[,i])
}
PseudoR_df<-do.call("rbind", PseudoR)
melted_df<-reshape2::melt(PseudoR_df)
colnames(melted_df)<-c("Ind", "Draw", "PseudoR")


## Extracting the simulation individual names and lakes, multiplying them with n.sims
number.sim<-jags_full[["BUGSoutput"]][["n.sims"]]
number.sources<-2
gg<-cbind(rep(mix[["data"]][["Lake"]], times=number.sim),  # times=2 indicates number of sources
          rep(mix[["data"]][["species"]], times=number.sim)) # times=2 indicates number of sources
colnames(gg)<-c("Lake", "species")
PseudoR_df_env<-cbind(melted_df, gg)
library(ggplot2)
ggplot(data=PseudoR_df_env[PseudoR_df_env$PseudoR>-1 & PseudoR_df_env$PseudoR<2,], aes(x=PseudoR))+ 
  facet_wrap(~species, scales = "free_y")+
  geom_histogram(aes(y=..density..), binwidth = 0.1, colour="black", fill="white")+
  geom_density(fill="#FF6666", alpha=0.5)+ 
  labs(y="Scaled density", x="McFadde pseudo R2")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        plot.title=element_text(size=10))

