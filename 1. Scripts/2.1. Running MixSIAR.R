rm(list = ls()) # clear R's clobal envirnoment memory
set.seed(180959779) ## I have sampled 1:10^9 once and use that as a seed. #sample(1:10^9,1)

### Author: Ossi Keva
### Date 03.06.2024
### Project: Roger Jones - ALLOCARB , Antti Eloranta - FreshRestore

### Short description: Here with this R script we are 
### At Sections 1.1-1.7. uploading consumer, source and dscrimination files
### At Sections 2.1 Running the MixSIAR models

### 1.1. Some packages to download ####
## These include some unnecessary packages, but I did not have time to recheck which are not needed
library(plyr) 
library("dplyr")
library("nlstools")
library("reshape2") 
library(ggplot2)
library("ggrepel")
library(R2WinBUGS)
library(splancs)
library(rjags)
library(MixSIAR)

### 1.2. Setting working directory ####
setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()
list.files()

### 1.3. Loading the data sets #####
# Load mix or consumer data
mix<-list()

## loading mix data for lake
mix[["Mod with lake w_0.23"]] <- load_mix_data(filename= paste0("2. Data/Consumerdata_w_0.23newTP_2.1.csv"), 
                                                iso_names    = c("d2H"), 
                                                factors      = c("Lake"),
                                                fac_random   = c(TRUE),   
                                                fac_nested   = NULL,  
                                                cont_effects = NULL) 
## loading mix data for lake+species
mix[["Mod with lake+species w_0.23"]] <- load_mix_data(filename= paste0("2. Data/Consumerdata_w_0.23newTP_2.1.csv"), 
                                          iso_names    = c("d2H"), 
                                          factors      = c("species", "Lake"),
                                          fac_random   = c(TRUE, TRUE),   
                                          fac_nested   = c(FALSE, FALSE),   #  Only applies in cases with 2 factors. NULL otherwise
                                          cont_effects = NULL) 

## loading mix data for lake+species+continous variables (PC1, FOR%,  POM N:C)
loop_env_vars<-loop_env_vars<-c("PC1", "For_area", "POM_NC_ratio")
for ( i in loop_env_vars){
    mix[[paste("Mod with lake+species_rand+", i, " w_0.23_randSlope", sep="")]] <- load_mix_data(filename= paste0("2. Data/Consumerdata_w_0.23newTP_2.1.csv"), 
                                                                                  iso_names    = c("d2H"),
                                                                                  factors      = c("species", "Lake"), 
                                                                                  fac_random   = c(TRUE, TRUE),   
                                                                                  fac_nested   = c(FALSE, FALSE), 
                                                                                  cont_effects = i) 
}

### Loading mix data for lake+species+PC1 with differing omega values
wtot <- c(0.14, 0.32)
for ( w in unique(wtot)){
  mix[[paste("Mod with lake+species_rand+PC1 w_", w, "_randSlope", sep="")]] <- load_mix_data(filename= paste0("2. Data/Consumerdata_w_",w,"newTP_2.1.csv"), 
                                                                     iso_names    = c("d2H"), 
                                                                     factors      = c("species", "Lake"),
                                                                     fac_random   = c(TRUE, TRUE),   
                                                                     fac_nested   = c(FALSE, FALSE),
                                                                     cont_effects = "PC1") 
}

### 1.4. Loading sources and discrimination files ####
source   <- list()
discr    <- list()
n.mod<-names(mix)
for(mod in unique(n.mod)){
    source[[mod]] <- load_source_data(filename="2. Data/Source data.csv", 
                                      source_factors=c("Lake"),
                                      conc_dep=FALSE,
                                      data_type="means",
                                      mix[[mod]])
  
  discr[[mod]] <- load_discr_data(filename="2. Data/tef.csv", mix[[mod]]) # tef5 d2H discrimination 0, tef6 d2H discrimination sd is =-13 based on d15N TDF 3.4=-1 and clado base 2.0=-0.1 impacts on environmental water corrected d2H values
}


### adding purely empty model to the run list. The source data is modelled based on the average and stds.
mix[["empty"]]<- load_mix_data(filename= paste0("2. Data/Consumerdata_w_0.23newTP_2.1.csv"), 
                               iso_names    = c("d2H"), 
                               factors      = NULL,
                               fac_random   = NULL,   
                               fac_nested   = NULL,  
                               cont_effects = NULL) 
source[["empty"]] <- load_source_data(filename="2. Data/Source data2.csv", 
                                      source_factors=NULL,
                                      conc_dep=FALSE,
                                      data_type="raw",
                                      mix[[mod]])
discr[["empty"]] <- load_discr_data(filename="2. Data/tef.csv", mix[[mod]])


### 1.5. setting output options ######
output_options <- list(summary_save = TRUE,                 
                       summary_name = "summary_statistics", 
                       sup_post = TRUE,       ## suppress the posterior density plots           
                       plot_post_save_pdf = FALSE,           
                       plot_post_name = "posterior_density",
                       sup_pairs = TRUE,         ## Suppress pairs plot output    
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,   ## suppress xy/trace plots        
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE, ## Save diagnostics as .txt
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE, ## Save diagnostis as pdf
                       return_obj = TRUE)

### 1.6. Uploading some help functions to export data functions ####
setwd(file.path(orgfolder, c("1. Scripts"), c("Help functions")))
#source("output_diagnostics.R")
#source("output_posteriors.R")
source("run_model2.R") ## This is modified run_model() function. This includes species specific random slope and sig in the the output parameters.
setwd(orgfolder)


### 2.1. MixSIAR runs ####
### Setting the length of the run, try first with test
n.mod<-names(mix)
run<-c("very long") # , "test", "very short", "short", normal ,"very long" 

# for loop for each model we want to run, we have defined these previously to mix, source and dicr lists
for (mod in c(n.mod[3:8],n.mod[1:2])){  
  setwd(orgfolder) ## This is also here just in case
  jags.mod<-NULL ## removing jags.mod before new run
  
  ## Lets make separate folders for each of the models and save our results there
  mainDir <- getwd()
  subDir <- paste0(run,"_",  mod)
  dir.create(file.path(mainDir,"3. MixSIAR Models", subDir), showWarnings = FALSE)
  setwd(file.path(mainDir,"3. MixSIAR Models", subDir))
  
  ### Test if the model have no continous effect
  if(is.null(mix[[mod]][["cont_effects"]])){ ## If model does not have cont, then do this
    model_filename <- paste0("MixSIAR_model", mod, ".txt") 
    resid_err <- TRUE
    process_err <- TRUE
    write_JAGS_model(model_filename, process_err, resid_err, mix[[mod]], source[[mod]]) ## write jags model and print it to the working directory
    jags.mod <- run_model(run=run, mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1,
                          process_err = TRUE, resid_err = TRUE)
    
  }else { ## If the model have continous effect, then do this, 
    ## We have manually modified "MixSIAR_random_continuous_model_v.1.txt" to include random slope effect, 
    ## We have also manually modified run_model2()-function to be able to print the added random slope effects to jags.mod
    jags.mod <- run_model2(run=run, mix[[mod]], source[[mod]], discr[[mod]], model_filename = 
                             file.path(orgfolder, "3. MixSIAR Models","MixSIAR_random_continuous_model_v.1.txt"), alpha.prior=1,
                           process_err = TRUE, resid_err = TRUE)}
  
  ## Lets save the jags object for further purposes
  saveRDS(jags.mod, file = paste0("mixSIAR_model_", mod, ".rds"))
  
  ### lets export the output diagnostics as well, check the output options, 
  ### at least summary statistics and diagnostics should be exported. 
  options(max.print = 99999)
  output_JAGS(jags.mod, mix[[mod]], source[[mod]], output_options = output_options)
  
  ### Exporting the mix, source and discr objects to the same folder just in case
  saveRDS(mix[[mod]], file = paste0("mix_object",".rds"))
  saveRDS(source[[mod]], file = paste0("source_object",".rds"))
  saveRDS(discr[[mod]], file = paste0("discr_object",".rds"))
    
  setwd(orgfolder)
} # END of mod loop

#####################
### END OF SCRIPT ###
#####################
