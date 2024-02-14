rm(list = ls()) # clear R's clobal envirnoment memory

## Author: Ossi Keva
## Date 14.2.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

### Short description: Here with this R script we are 
### At Sections 1.1-1.7. running the mixSIAR with d2H values consumer values

### 1.1. Some packages to download ####
library(plyr) 
library(dplyr)
library(reshape2) 
library(ggplot2)
library(ggrepel)
library(R2WinBUGS)
library(splancs)
library(rjags)
library(MixSIAR)
#remotes::install_github("brianstock/MixSIAR", dependencies=T) ## This should include output_posteriors, but I was unable to download this package from github

### 1.2. setting working directory ####
setwd("D:\\Ossi Keva")
orgfolder<-getwd()
list.files()

### 1.3. Loading the data sets #####
# Load mix or consumer data
mix<-list()

## loading mix data for lake+species+continous variables (PC1, FOR%,  POM N:C)
loop_env_vars<-loop_env_vars<-c("PC1", "For_area", "POM_NC_ratio")
for ( i in loop_env_vars){
    mix[[paste("Mod with lake+species+", i, " w_0.23", sep="")]] <- load_mix_data(filename= paste0("2. Data/Consumer data_w_0.23.csv"), 
                                                                                  iso_names    = c("d2H"), ## This time only hydrogen
                                                                                  factors      = c("Lake", "species"), #factors=c("Lake", "taxa"),
                                                                                  fac_random   = c(TRUE, FALSE),   #fac_random = c(FALSE,TRUE), ## TRUE FLASE ; False-> fixed effect in the model
                                                                                  fac_nested   = c(FALSE, FALSE),   # fac_nested   = c(TRUE,FALSE), Only applies in cases with 2 factors. NULL otherwise
                                                                                  cont_effects = i) # i

}
### Loading mix data for lake+species+PC1 with differing omega value range 0.14-0.32
wtot <- c(0.14, 0.32)
for ( w in unique(wtot)){
  mix[[paste("Mod with lake+species+PC1 w_", w, sep="")]] <- load_mix_data(filename= paste0("2. Data/Consumer data_w_",w,".csv"), 
                                                                     iso_names    = c("d2H"), ## This time only hydrogen
                                                                     factors      = c("Lake", "species"), #factors=c("Lake", "taxa"),
                                                                     fac_random   = c(TRUE, FALSE),   #fac_random = c(FALSE,TRUE), ## TRUE FLASE ; False-> fixed effect in the model
                                                                     fac_nested   = c(FALSE, FALSE),   # fac_nested   = c(TRUE,FALSE), Only applies in cases with 2 factors. NULL otherwise
                                                                     cont_effects = "PC1") # i
  
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
  
    discr[[mod]] <- load_discr_data(filename="2. Data/tef.csv", mix[[mod]]) 
    # tef.csv, d2H discrimination sd is =-13 based on d15N TDF 3.4+-1 and clado base 2.0+-0.1 
    # impacts on environmental water corrected d2H values
}

### 1.5. Uploading some help functions to modify MixSIAR outputs ####
## just in case if your MixSIAR version do not have these functions
setwd(file.path(orgfolder, c("1. Scripts"), c("Help functions")))
source("output_diagnostics.R")
source("output_posteriors.R")
setwd(orgfolder)

### setting output options #####
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

### 1.6. Then running the models ####
run<-c("long") # , "test", "very short", "short", normal ,"very long"
for (j in unique(run)){
  for (mod in unique(n.mod)){
    jags.mod<-list() ## I placed this here as my intent is to save this object after each iteration, to prevent the R memory to get too full
    mainDir <- getwd()
    subDir <- paste0(j,"_",  mod)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))
    
    model_filename <- paste0("MixSIAR_model", mod, ".txt")
    resid_err <- TRUE
    process_err <- TRUE
    write_JAGS_model(model_filename, process_err, resid_err, mix[[mod]], source[[mod]])
    jags.mod[[mod]] <- run_model(run=j, mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1,
                                 process_err = TRUE, resid_err = TRUE)
    saveRDS(jags.mod[[mod]], file = paste0("mixSIAR_model_", mod, ".rds"))

    ### lets export the output_diagnostics as well, check the output options, 
    ### at least summary statistics and diagnostics should be exported. 
    options(max.print = 99999)
    output_JAGS(jags.mod[[mod]], mix[[mod]], source[[mod]], output_options = output_options)
    
    ### Exporting the mix, source and discr objects to the same folder just in case
    saveRDS(mix[[mod]], file = paste0("mix_object",".rds"))
    saveRDS(source[[mod]], file = paste0("source_object",".rds"))
    saveRDS(discr[[mod]], file = paste0("discr_object",".rds"))
    
    setwd(orgfolder)
  }
}

#########################
### END OF THE SCRIPT ###
#########################
