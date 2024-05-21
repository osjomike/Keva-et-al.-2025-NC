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

### Typing down The folder names we process
list.files()
Folders<-c("long_Mod with lake+species_rand+PC1 w_0.23_randSlope",
           "long_Mod with lake+species_rand+PC1 w_0.14_randSlope",
           "long_Mod with lake+species_rand+PC1 w_0.32_randSlope")

#Creating a list for data frames
df_list<-list()
f=Folders[1]
## Running a loop that extract linear coefficient impact on dietary proportions to df_list
for (f in unique(Folders[1])){
  #### 0 Reading some parameters from folder names and reading JAGS objects ####
  env<-sub(".*species_rand\\+(.*)( w_).*", "\\1", f) ## getting the environmental variable name from the model file name
  wtot<-sub(".*w_(.*)(_randSlope).*", "\\1", f) ## getting the used omega value from the model file name
  setwd(file.path(orgfolder,c("3. MixSIAR Models"), f))
  list.files()
  
  ### uploading mixtures, sources, discrimination and jags object
  mix<-readRDS(file = paste0("mix_object",".rds"))
  source<-readRDS(file = paste0("source_object",".rds"))
  discr<-readRDS(file = paste0("discr_object",".rds"))
  jags<-readRDS(file = paste0("mixSIAR_model_Mod with lake+species_rand+", env, " w_",wtot,"_randSlope.rds"))

  
  #### The intention here is to show lines for one source (Ter) + modifying the colours according to omega values
  
  #### 1.1 Continous variable plot for species #########
  R2jags::attach.jags(jags)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  alphaCI=0.05
  exclude_sources_below=0.001
  df_list<-list()
  ce=1
  ### FAC1 is species and FAC2 is lakes, both are random
  ## Here in this loop we are summing up global+rand_slope+species in ILR space to retrieve the continous variable impact
  for(f1 in 1:mix$FAC[[1]]$levels){ ### 
    fac.lab <- mix$FAC[[1]]$labels[f1]
    label <- mix$cont_effects[ce]
    cont <- mix$CE[[ce]]
    s_cont<-mix$CE[[ce]][mix$data$species==fac.lab,]
    ilr.cont <- get(paste("ilr.cont",ce,sep=""))
    
    n.plot = 200
    chain.len = dim(p.global)[1]
    Cont1.plot <- seq(from=round(min(s_cont),5), to=round(max(s_cont),5), length.out=n.plot)
    ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
    for(src in 1:n.sources-1){
      for(i in 1:n.plot){
        ## Note that in the MixSIAR model we used customed jags code and run_model2() function to retrieve ilr.cont1.r object
        ## In the custom model the slope can vary between species!
        ilr.plot[i,src,] <- ilr.global[,src] + (ilr.cont[,src]+ ilr.cont1.r[,f1,src])*Cont1.plot[i] + ilr.fac1[,f1,src]
      }
    }
    
    # Transform every draw from ILR-space to p-space
    e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
    for(i in 1:(n.sources-1)){
      e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
      e[,i] <- e[,i]/sum(e[,i])
    }
    # dummy variables for inverse ILR calculation
    cross <- array(data=NA,dim=c(n.plot, chain.len, n.sources, n.sources-1))  
    tmp <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
    p.plot <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
    for(i in 1:n.plot){
      for(d in 1:chain.len){
        for(j in 1:(n.sources-1)){
          cross[i,d,,j] <- (e[,j]^ilr.plot[i,j,d])/sum(e[,j]^ilr.plot[i,j,d])
        }
        for(src in 1:n.sources){
          tmp[i,d,src] <- prod(cross[i,d,src,])
        } # end of src loop
        for(src in 1:n.sources){
          p.plot[i,d,src] <- tmp[i,d,src]/sum(tmp[i,d,])
        }# end of src loop
      } #end of d loop
    } # end of i loop
    
    # now take quantiles, after ILR transform of every draw
    get_high <- function(x){return(quantile(x, 1-alphaCI/2))}
    get_low <- function(x){return(quantile(x, alphaCI/2))}
    p.low <- apply(p.plot, c(1,3), get_low)
    p.high <- apply(p.plot, c(1,3), get_high)
    p.median <- apply(p.plot, c(1,3), median)
    colnames(p.median) <- source_names
    
    Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
    df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
    colnames(df) <- c("source","median","x","low","high")
    df$source <- factor(df$source, levels=source_names)
    
    # remove sources from plot with very low proportions
    rm.srcs <- apply(p.median, 2, function(x) all(x < exclude_sources_below))
    df <- subset(df, source %in% source_names[!rm.srcs])
    
    ### Adding species and omega values for the dataframes
    df$species<-fac.lab
    df$omega<-wtot
    ## saving the dafs to the df_list
    df_list[[fac.lab]][[paste0("omega=", wtot)]]<-df
  } # End of f loop (species)
  
  
  #### 2.1 Continous variable plot for Global ####
  
  ### First we derive continous variable impact on allochthony, this will be exported to a dataframe 
  R2jags::attach.jags(jags)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  alphaCI=0.05 ## Set the Credibility interval alpha
  exclude_sources_below=0.001
  ce=1
  label <- mix$cont_effects[ce]
  cont <- mix$CE[[ce]]
  ilr.cont <- get(paste("ilr.cont",ce,sep=""))
  
  n.plot = 200
  chain.len = dim(p.global)[1]
  Cont1.plot <- seq(from=round(min(cont),5), to=round(max(cont),5), length.out=n.plot)
  ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
  for(src in 1:n.sources-1){
    for(i in 1:n.plot){
      ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i]
    }
  }
  
  # Transform regression lines from ILR-space to p-space
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }
  # dummy variables for inverse ILR calculation
  cross <- array(data=NA,dim=c(n.plot, chain.len, n.sources, n.sources-1))  
  tmp <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
  p.plot <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
  for(i in 1:n.plot){
    for(d in 1:chain.len){
      for(j in 1:(n.sources-1)){
        cross[i,d,,j] <- (e[,j]^ilr.plot[i,j,d])/sum(e[,j]^ilr.plot[i,j,d])
      }
      for(src in 1:n.sources){
        tmp[i,d,src] <- prod(cross[i,d,src,])
      }
      for(src in 1:n.sources){
        p.plot[i,d,src] <- tmp[i,d,src]/sum(tmp[i,d,])
      }
    }
  }
  # now take quantiles, after ILR transform of every draw
  get_high <- function(x){return(quantile(x, 1-alphaCI/2))}
  get_low <- function(x){return(quantile(x, alphaCI/2))}
  p.low <- apply(p.plot, c(1,3), get_low)
  p.high <- apply(p.plot, c(1,3), get_high)
  p.median <- apply(p.plot, c(1,3), median)
  colnames(p.median) <- source_names
  
  Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
  df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
  colnames(df) <- c("source","median","x","low","high")
  df$source <- factor(df$source, levels=source_names)
  
  # remove sources from plot with very low proportions
  rm.srcs <- apply(p.median, 2, function(x) all(x < exclude_sources_below))
  df <- subset(df, source %in% source_names[!rm.srcs])
  ## Adding omega and species names to the dataframe
  df$omega<-wtot
  df$species<-"Global"
  df_list[["Global"]][[paste0("omega=", wtot)]]<-df
  
  setwd(orgfolder)
} # End of folder loop

#### 3.1 Making the omega sensitivity plots ####
## Then we extract the cont effect dataframes from the df_list, and make the plots
plot_list<-list()

## Define omega colours and order
omega_order<-c("0.14","0.32", "0.23")
my_colours<-c("#EF7F2A", "#A79ECD", "brown")

for(sp in names(df_list)){
  p_df<-do.call("rbind", list(df_list[[sp]][["omega=0.14"]], 
                              df_list[[sp]][["omega=0.23"]], 
                              df_list[[sp]][["omega=0.32"]])) 

  p_df$omega<-factor(p_df$omega, levels=omega_order)
  plot_list[[sp]] <- ggplot(data=p_df[which(p_df$source=="Ter"),], aes(x=x,y=median)) +
    geom_line(aes(x=x, y=median, colour=omega), size=1.5) +
    geom_ribbon(data= p_df[which(p_df$source=="Ter" & p_df$omega!="0.23"),], ## Drawing 95% CI only to upper and lower omega
                aes(ymin=low, ymax=high, fill=omega), alpha=0.35)+
    labs(title = sp) +
    ylab("Proportion") +
    xlab(label) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
    scale_colour_manual(labels=omega_order,values = my_colours)+
    scale_fill_manual(labels=omega_order, values=my_colours)+
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"), axis.title=element_text(size=16), 
                   axis.text=element_text(size=14), legend.text=element_text(size=14), legend.position=c(.02,1), 
                   legend.justification=c(0,1), legend.title = element_blank())
  
}       
## Getting the legend
legend<-get_legend(plot_list[[1]]+labs(colour=c("Omega values"), fill=c("Omega values"))+
                     scale_fill_manual(labels=omega_order, values=my_colours,
                                       aesthetics = c("colour", "fill"),
                                       limits = omega_order)+ 
                     #labs(colour="Test", fill="test")+
                     theme(legend.position="right",
                           legend.box = "horizontal", 
                           legend.background = element_rect(linetype = 2, size = 0.5, colour = 1),
                           legend.title=element_text(),
                           legend.justification = NULL,
                           legend.margin=margin(c(5,5,5,5)),
                           plot.margin=unit(c(0,0,0,0), units="line")))
#plot_grid(legend) # just cheking if the extraction worked

fig_order<-c("Global","Bulk_BMI_littoral", "Perch_large", "Pike", 
             "Bulk_ZPL", "Asellus_littoral", "Perch_small", "Smelt",
             "Chaoborus", "Chironomid_littoral", "Vendace", "Bleak", 
             "Cladocera", "Bulk_BMI_profundal", "Perch_medium", "Roach_small",
             "Copepods", "Chironomid_Profundal", "Ruffe", "Roach_large")

comb_fig<-cowplot::plot_grid(plotlist = plot_list[fig_order], ncol=4)

#### 4.1 Exporting the created figure ######

### Nature quidelines 183 mm (double column) and the full depth of the page is 170 mm. max width=7.2in, height = 6.7in
#this was the furthest I could get with R at this time, modify axis titles, add y- and x-axis tick labels, subplot group names etc small stuff in vector graphics program
pdf(file=file.path(#orgfolder, "4. Figures", 
  paste0(env,"_Continous variable plots with posterior densities_test.pdf")), 
  height = 7, width = 6.5) 
if(grepl("PC1", env)){ ## for some reason grid.arrage did not like the ifs to be inside the function. It messed the figure ratios. Therefore I have now placed it like this. 
  grid.arrange(comb_fig, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
               top= textGrob(f ,gp = gpar(fontsize = 11)),
               bottom = textGrob(paste0("Forestry, DOC <--  ", env, " axis  --> Agriculture, nutrients"),gp = gpar(fontsize = 11)))
} else {
  grid.arrange(comb_fig, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
               top= textGrob(f ,gp = gpar(fontsize = 11)),
               bottom = textGrob(env,gp = gpar(fontsize = 11)))
}
dev.off()

#########################
### END OF THE SCRIPT ###
#########################