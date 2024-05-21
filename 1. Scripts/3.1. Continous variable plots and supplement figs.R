rm(list=ls())

## Author: Ossi Keva
## Date 14.2.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description: With this script we are plotting the covariate effects on consumer allochthony,
## including the lake and species specific posterior distribution 95CI for Terrestrial source
## Article Figures 3, Extended data Fig 3-4

# Some packages to download
library(plyr) # data arrangement
library(dplyr) # data arrangement
library(ggplot2) # drawing of plots
library(ggpubr)  
library(cowplot) # get legend and plot arranging
library(gridExtra) # arranging of plots
library(MixSIAR) #
library(grid)

## Selecting working directory
setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()

## Setting the colours to be used in the figures
my_colours1<-c("lightblue", "brown") ## Same order that occurs in the source files 1=aq=lightblue, 2=Ter=brown

### Typing down The folder names we process
list.files(file.path(orgfolder, "3. MixSIAR Models"))
Folders<-c("long_Mod with lake+species_rand+PC1 w_0.23_randSlope",
           "long_Mod with lake+species_rand+For_area w_0.23_randSlope",
           "long_Mod with lake+species_rand+POM_NC_ratio w_0.23_randSlope")

### With this loop, we are creating Article Fig 3 for different models (environmental factors)
### This produces Continous variable effect on consumer mixtures using plot_cont_var2() function from MixSIAR package (name in MixSIAR: plot_continous_var)
### And we supplement this figure with Lake and species specific posterior distributions
### In addition we make posterior density plots for for each lake and species

for (f in Folders){
  #### 0 Reading some parameters from folder names and reading JAGS objects ####
  env<-sub(".*species_rand\\+(.*)( w_).*", "\\1", f) ## getting the environmental variable name from the model file name
  #wtot<-sub(".*w_ *(.*?)", "\\1", f) ## getting the used omega value from the model file name 
  wtot<-0.23
  setwd(file.path(orgfolder, "3. MixSIAR Models", f))#file.path(orgfolder,c("3. MixSIAR Models"), f))
  list.files()
  
  ### uploading mixtures, sources, discrimination and jags object
  mix<-readRDS(file = paste0("mix_object",".rds"))
  source<-readRDS(file = paste0("source_object",".rds"))
  discr<-readRDS(file = paste0("discr_object",".rds"))
  jags<-readRDS(file = paste0("mixSIAR_model_Mod with lake+species_rand+", env, " w_",wtot,"_randSlope.rds"))
  #list.files()
  
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
    df$species<-fac.lab
    df_list[[fac.lab]]<-df ## saving the data to a list
  }

  #### 1.2 Extracting posterior densities to a dataframe ####
  ## Look the JAGS code, p.ind object is the one where all factors (including continous variable) is taken into account
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
            rep(rep(mix[["data"]][[label]], each=number.sim), times=number.sources)) # times=2 indicates number of sources
  colnames(gg)<-c("Lake", "species", label)
  posterior<-cbind(posterior_df, gg)
  
  ## making sure the levels are in correct order
  posterior$Lake<-factor(posterior$Lake, levels = unique(mix[["data"]][["Lake"]]))
  posterior$species<-factor(posterior$species, levels = unique(mix[["data"]][["species"]])[order(unique(mix[["data"]][["species"]]))])
  posterior$Source2<-factor(posterior$Source2, levels=source[["source_names"]])
  
  #colnames(posterior)
  #library(dplyr)
  ### Calculating summary statistics from posterior results
  posterior_summary<-posterior %>%
    group_by(Lake, species, Source2) %>%
    summarise(y0.025 = quantile(distribution, 0.025, na.rm = TRUE),
              y0.1= quantile(distribution, 0.1, na.rm = TRUE),
              y0.25= quantile(distribution, 0.25, na.rm = TRUE),
              y0.375= quantile(distribution, 0.375, na.rm = TRUE),
              y0.4= quantile(distribution, 0.4, na.rm = TRUE),
              y0.5= quantile(distribution, 0.5, na.rm = TRUE),
              y0.625= quantile(distribution, 0.625, na.rm = TRUE),
              y0.75= quantile(distribution, 0.75, na.rm = TRUE),
              y0.8= quantile(distribution, 0.8, na.rm = TRUE),
              y0.9= quantile(distribution, 0.9, na.rm = TRUE),
              y0.975= quantile(distribution,0.975, na.rm = TRUE))
  
  ### Adding the environmental variable used in the mixSIAR model and turning the blody tidy-objects to dfs
  posterior_summary$ENV<-factor(posterior_summary$Lake, labels =unique(mix[["data"]][[label]]))
  posterior_summary$ENV<-as.numeric(as.character(posterior_summary$ENV))
  colnames(posterior_summary)[colnames(posterior_summary)=="ENV"]<-env
  posterior_summary_env<-as.data.frame(posterior_summary)
  
  ## Unlisting df_list (this includes the proportion trends) to a data frame
  df_list_df<-do.call("rbind", df_list)
  
  ### Plottin all the different plots separately
  fig_list<-list()
  for (i in unique(df_list_df$species)){
    sub_trend<-df_list_df[df_list_df$species==i,]
    sub_points<-posterior_summary_env[posterior_summary_env$species==i,]
    
    fig_list[[i]]<-ggplot() +
      geom_line(data=sub_trend, aes(x=x, y=median,group=source,colour=source),size=1.5) +
      geom_ribbon(data=sub_trend, aes(x=x, ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
      scale_fill_manual(values=my_colours1)+
      scale_colour_manual(values=my_colours1)+
      labs(x=label, y="Proportion", title=i)+
      scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
      lims(x=c(min(min(mix$data[env]), min(df_list_df$x)), 
               max(max(mix$data[env]), max(df_list_df$x))))+
      geom_point(data=sub_points[sub_points$Source2=="Ter",], aes_string(x=label, y="y0.5"), colour="#7c2020")+
      geom_linerange(inherit.aes = FALSE, data=sub_points[sub_points$Source=="Ter" ,], 
                     aes_string(x=label, ymin="y0.025", ymax="y0.975"), show.legend = FALSE, colour="#7c2020", alpha=0.35)+
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black"), 
            axis.text=element_text(size=14),legend.position="none",
            axis.text.x=element_blank(), #remove x axis labels
            #axis.ticks.x=element_blank(), #remove x axis ticks
            axis.text.y=element_blank(),  #remove y axis labels
            #axis.ticks.y=element_blank(),  #remove y axis ticks
            axis.title = element_blank(),
            plot.title=element_text(size=10))
  }
  legend<-get_legend(ggplot() +
                     geom_line(data=sub_trend, aes(x=x, y=median,group=source,colour=source),size=1.5) +
                       geom_ribbon(data=sub_trend, aes(x=x, ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
                       scale_fill_manual(values = my_colours1, labels=c('Aquatic', 'Terrestrial'), name = "Dietary source")+ ## These are alphabetically orderder by their original labels which are Aq and Ter
                       scale_colour_manual(values = my_colours1, labels=c('Aquatic', 'Terrestrial'), name = "Dietary source")+ 
                       theme(legend.position="right", legend.box = "horizontal", 
                             legend.background = element_rect(linetype = 2, size = 0.5, colour = 1),
                             legend.title=element_text(),legend.justification = NULL,
                             legend.margin=margin(c(5,5,5,5)),plot.margin=unit(c(0,0,0,0), units="line")))
  legend1<-plot_grid(legend)

  ## Setting figure order
  fig_order<-c("Bulk_BMI_littoral", "Perch_large", "Pike", 
           "Bulk_ZPL", "Asellus_littoral", "Perch_small", "Smelt",
           "Chaoborus", "Chironomid_littoral", "Vendace", "Bleak", 
           "Cladocera", "Bulk_BMI_profundal", "Perch_medium", "Roach_small",
           "Copepods", "Chironomid_Profundal", "Ruffe", "Roach_large")
  
  new_figlist1<-list()
  new_figlist1<-fig_list[fig_order]
  new_figlist1<- append(new_figlist1,list(legend1), after = 0)
  comb_fig<-cowplot::plot_grid(plotlist = new_figlist1, ncol=4)

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
  
  
  #### 1.3 Continous variable plot for speciesthe global plot #######
  #### First we derive continous variable impact on allochthony, this will be exported to a dataframe  
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
  
  #### 2.1 Then we need to derive lake global mean values ####
  ### For Global lake estimate one needs to global+lake+continous effect
  ilr.globalX<-jags[["BUGSoutput"]][["sims.list"]][["ilr.global"]]
  ilr.fac2X<-jags[["BUGSoutput"]][["sims.list"]][["ilr.fac2"]]
  ilr.cont1X<-jags[["BUGSoutput"]][["sims.list"]][["ilr.cont1"]]
  ilr.cont1.rX<-jags[["BUGSoutput"]][["sims.list"]][["ilr.cont1.r"]]
  dim(ilr.cont1.rX)
  dim(ilr.fac2X)
  ## Dummy variable for p-conversion
  ilr.fac2_tot<-array(dim = c(chain.len,35)) ## 3000 is draws and 35 is lake
  for(i in 1:chain.len){
    ilr.fac2_tot[i,] <-ilr.globalX[i]+ilr.fac2X[i,,]+ (ilr.cont1X[i])*unique(mix$CE[[ce]]) ## Lets add the effects together for each draw
  }

  # Transform ilr.fac2_tot from ILR-space to p-space
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }
  
  # dummy variables for inverse ILR calculation
  cross <- array(data=NA,dim=c(chain.len,35, n.sources))  
  tmp <- array(data=NA,dim=c(chain.len,35 ,n.sources))  
  p.lake <- array(data=NA,dim=c(chain.len,35, n.sources))  
    for(d in 1:chain.len){
      for(l in 1:35){
        for(j in 1:(n.sources-1)){
          cross[d,l,] <- (e[,j]^ilr.fac2_tot[d,l])/sum(e[,j]^ilr.fac2_tot[d,l])
        }
        for(src in 1:n.sources){
          tmp[d,l,src] <- prod(cross[d,l,src])
        }
        for(src in 1:n.sources){
          p.lake[d,l,src] <- tmp[d,l,src]/sum(tmp[d,l,])
        }
    }
  }
 
  # Lake posterior means are now in p.lake
  ## Lets make a dataframe of it
  p.lake_df<-reshape2::melt(p.lake)
  colnames(p.lake_df)<-c("Draw", "Lake", "Source", "distribution")
  
  ## Here lets rename the source names based on the source object levels
  p.lake_df$Source2<-factor(p.lake_df$Source, labels = source[["source_names"]])
  p.lake_df$Source2<-factor(p.lake_df$Source2, levels=source[["source_names"]])
  
  ## Making another lake variable, names based on mix object
  p.lake_df$Lake2<-factor(p.lake_df$Lake, labels = unique(mix[["data"]][["Lake"]]))

  ### Calculating summary statistics from posterior results
  lake_post_sum<-p.lake_df %>%
    group_by(Lake2, Source2) %>%
    summarise(y0.025 = quantile(distribution, 0.025, na.rm = TRUE),
              y0.5= quantile(distribution, 0.5, na.rm = TRUE),
              y0.975= quantile(distribution,0.975, na.rm = TRUE))
  ## Adding env variable to the plot
  lake_post_sum$ENV<-factor(lake_post_sum$Lake2, labels =unique(mix[["data"]][[env]]))
  lake_post_sum$ENV<-as.numeric(as.character(lake_post_sum$ENV))
  colnames(lake_post_sum)[colnames(lake_post_sum)=="ENV"]<-env
  lake_post_sum_df<-as.data.frame(lake_post_sum)
  
  #### 2.2 Making a ggplot object of global #####
  global_fig<-ggplot()+
    geom_line(data=df, aes(x=x, y=median,group=source,colour=source),size=1.5) +
    geom_ribbon(data=df, aes(x=x, ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
    
    geom_point(inherit.aes = FALSE, data=lake_post_sum_df[lake_post_sum_df$Source2=="Ter",], aes_string(x=env, y="y0.5"), colour="#7c2020")+
    geom_linerange(inherit.aes = FALSE, data=lake_post_sum_df[lake_post_sum_df$Source2=="Ter",], 
                   aes_string(x=env, ymin="y0.025", ymax="y0.975"), show.legend = FALSE, colour="#7c2020", alpha=0.35)+
    
    scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
    lims(x=c(min(min(mix$data[env]), min(df$x)), 
             max(max(mix$data[env]), max(df$x))))+
    scale_fill_manual(values=my_colours1)+
    scale_colour_manual(values=my_colours1)+
    ylab("Proportion") +
    xlab(env)+
    labs(title="Global")+
    theme_bw() +
    theme(legend.position = c(0.9,0.9))
  
  pdf(file=file.path(orgfolder,"3. MixSIAR Models",  f, "Global continous variable.pdf"),
      height = 3.5, width = 3.5)
  print(global_fig)
  dev.off()

  #### 2.3 Ah we can combine the global and the species level plots!! #### 
  ##The original species stuff was in fig_list 
  ## And the order is fig_order
  
  new_figlist2<-fig_list[fig_order]
  new_figlist2<-append(new_figlist2, list(global_fig+   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                                              panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                              axis.line = element_line(colour = "black"), 
                                                              axis.text=element_text(size=14),legend.position="none",
                                                              axis.text.x=element_blank(), #remove x axis labels
                                                              #axis.ticks.x=element_blank(), #remove x axis ticks
                                                              axis.text.y=element_blank(),  #remove y axis labels
                                                              #axis.ticks.y=element_blank(),  #remove y axis ticks
                                                              axis.title = element_blank(),
                                                              plot.title=element_text(size=10))), 
                       after =0)
  comb_fig2<-cowplot::plot_grid(plotlist = new_figlist2, ncol=4)
  
  #### 3.1 Exporting the figure Article figure 3, ###### 
  ### Nature quidelines 183 mm (double column) and the full depth of the page is 170 mm. max width=7.2in, height = 6.7in
  #this was the furthest I could get with R at this time, modify axis titles, add y- and x-axis tick labels, subplot group names etc small stuff in vector graphics program
  pdf(file=file.path(orgfolder, "4. Figures", 
                     paste0(env,"_Continous variable plots with posterior densities.pdf")), 
      height = 7, width = 6.5) 
  if(grepl("PC1", env)){ ## for some reason grid.arrage did not like the ifs to be inside the function. It messed the figure ratios. Therefore I have now placed it like this. 
    grid.arrange(comb_fig2, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
                 top= textGrob(f ,gp = gpar(fontsize = 11)),
                 bottom = textGrob(paste0("Forestry, DOC <--  ", env, " axis  --> Agriculture, nutrients"),gp = gpar(fontsize = 11)))
  } else {
    grid.arrange(comb_fig2, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
                 top= textGrob(f ,gp = gpar(fontsize = 11)),
                 bottom = textGrob(env,gp = gpar(fontsize = 11)))
  }
  dev.off()
  
  
  
  
  #### 4.1 Lets plot then the species specific slope distributions #########
  #### This is Article supplementary figure 
  
  con<-jags[["BUGSoutput"]][["sims.list"]][["ilr.cont1"]]
  con.r<-jags[["BUGSoutput"]][["sims.list"]][["ilr.cont1.r"]]
  
  test_array<-array(NA, dim=c(3000, 19, 1))
  for(i in 1:3000){
    test_array[i,,]<-con.r[i,,]+con[i,]
  }  
  slopes<-reshape2::melt(test_array)
  colnames(slopes)<-c("draw", "species", "var", "Slope_ILR")
  slopes$species2<-factor(slopes$species, labels = unique(mix[["data"]][["species"]])[order(unique(mix[["data"]][["species"]]))])
  slopes$species2<-factor(slopes$species2, levels  = unique(mix[["data"]][["species"]])[order(unique(mix[["data"]][["species"]]))])
  
  slopes_summary<-slopes %>%
    group_by(species2) %>%
    summarise(y0.025 = quantile(Slope_ILR, 0.025, na.rm = TRUE),
              y0.5= quantile(Slope_ILR, 0.5, na.rm = TRUE),
              y0.975= quantile(Slope_ILR,0.975, na.rm = TRUE))
  slopes_summary<-as.data.frame(slopes_summary)
  
  scaleFUN <- function(x) sprintf("%.3f", x)
  slopes_summary$label<-paste0("median (95% CI) = \n", scaleFUN(slopes_summary[,3]), 
                               " (", scaleFUN(slopes_summary[,2]), " – ", scaleFUN(slopes_summary[,4]), ")")
  
  global_slope<-data.frame(Slope_ILR=jags[["BUGSoutput"]][["sims.list"]][["ilr.cont1"]])
  global_sum<-global_slope %>%
    summarise(y0.025 = quantile(Slope_ILR, 0.025, na.rm = TRUE),
              y0.5= quantile(Slope_ILR, 0.5, na.rm = TRUE),
              y0.975= quantile(Slope_ILR,0.975, na.rm = TRUE))
  global_sum<-as.data.frame(global_sum)
  global_sum$label<-paste0("median (95% CI) = \n", scaleFUN(global_sum[,2]), 
                               " (", scaleFUN(global_sum[,1]), " – " , scaleFUN(global_sum[,3]), ")")
  max(global_slope)
  global_slopefig<-ggplot(global_slope, aes(x=Slope_ILR))+
    geom_histogram(aes(y=..density..), binwidth = 0.1, colour="black", fill="white")+
    geom_density(fill="#FF6666", alpha=0.5)+ 
    geom_vline(xintercept = 0, linetype="dashed", colour="grey")+
    coord_cartesian(xlim = c(quantile(slopes$Slope_ILR, probs = 0.0001),quantile(slopes$Slope_ILR, probs = 0.9999)), ylim = c(0,4))+
    geom_text(data=global_sum, aes(label=label, x=-Inf, y=Inf), size=1.8, vjust=1, hjust=-0.05)+
    #labs(y="Scaled density", x="Linear coefficient with PC1 in ILR space")+
    labs(y=NULL, x=NULL, title="Global")+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          plot.title=element_text(size=10))
  
  slope_figs<-list()
  for (i in unique(slopes$species2)){
    slopes_sub<-slopes[slopes$species2==i,]
    slopes_summary_sub<-slopes_summary[slopes_summary$species2==i,]
    
    slope_figs[[i]]<-ggplot(slopes_sub, aes(x=Slope_ILR))+ 
      geom_histogram(aes(y=..density..), binwidth = 0.1, colour="black", fill="white")+
      geom_density(fill="#FF6666", alpha=0.5)+ 
      geom_vline(xintercept = 0, linetype="dashed", colour="grey")+
      coord_cartesian(xlim = c(quantile(slopes$Slope_ILR, probs = 0.0001),quantile(slopes$Slope_ILR, probs = 0.9999)), ylim = c(0,4))+
      geom_text(data=slopes_summary_sub, aes(label=label, x=-Inf, y=Inf), size=1.8, vjust=1, hjust=-0.05)+
      #labs(y="Scaled density", x="Linear coefficient with PC1 in ILR space")+
      labs(y=NULL, x=NULL, title=i)+
      theme_classic()+
      theme(axis.text.y = element_blank(),
            plot.title=element_text(size=10))
  }
  new_figlist<-list()
  new_figlist<-slope_figs[fig_order] ## Fig_order is saved earlier
  new_figlist<- append(new_figlist,list(global_slopefig), after = 0)
  comb_slope_fig<-cowplot::plot_grid(plotlist = new_figlist, ncol=4)
  
  pdf(file="slope distribution in ILR_space.pdf",
      height = 7, width = 6.5)
  grid.arrange(comb_slope_fig,
               left = textGrob("Scaled density",rot=90, gp = gpar(fontsize = 11)), 
               top= textGrob(f ,gp = gpar(fontsize = 11)),
               bottom = textGrob("Linear coefficient with PC1 in ILR space",gp = gpar(fontsize = 11)))  
  dev.off()

  #### 5.1 Plotting supplementary plots i.e. Lake and species specific posterior densities ####
  ### Lets first plot dietary distribution of species to different plots
  plot_list_temp<-list()
  ## With this function, we can control the number of digits shown in the ggplot figure axes
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  for(i in levels(posterior$species)){
    plot_list_temp[[i]]<-ggplot(posterior[which(posterior$species==i),]) +
      facet_wrap(~Lake) +
      geom_density(alpha=0.3 ,aes(x=distribution*100,  y=after_stat(scaled), colour=Source2, fill=Source2))+ # to get geom_text to work mapping=aes(), needs to be done here!
      scale_fill_manual(values=my_colours1)+
      scale_colour_manual(values=my_colours1)+
      scale_y_continuous(labels=scaleFUN)+
      labs(x="Modelled Dietary Contribution (%)", y="Scaled posterior density",
           colour="Source2", fill="Source2",
           title = paste(i,f, "Diet sources \n MixSIAR", sep="_"))+
      theme_classic()+
      theme(panel.border = element_rect(colour = "Black", fill = NA, linetype = 1),
            legend.direction = "horizontal",
            legend.position = "bottom",
            legend.background  = element_rect(colour = "black", fill=NA, linetype = 2),
            strip.placement = "outside",
            plot.title=element_text(size=12))
  } # end of i loop
  
  ### Then lets make a pdf that includes all the species
  pdf(file=file.path(orgfolder, "3. MixSIAR Models", f, "Dietary contribution_by lake and species.pdf"))
  print(plot_list_temp) ## ordering the list alphabetically
  dev.off()
  
  
  

  #### 6.1 Then lets produce Table S6 type for each of the runs #####

} # End of f loop  

#########################
### END OF THE SCRIPT ###
#########################