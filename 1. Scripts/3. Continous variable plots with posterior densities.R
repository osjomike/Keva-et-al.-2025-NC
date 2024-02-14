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

## Selecting working directory
setwd("D:\\Keva et al. NC 2024")
orgfolder<-getwd()

## Some help functions to use
setwd(file.path(orgfolder, "1. Scripts/Help functions"))
source("Theme_functions.r")
source("plot_continous_var2.R")
source("output_diagnostics.R")
source("output_posteriors.R")
setwd(orgfolder)

#list.files()

## Setting the colours to be used
my_colours1<-c("lightblue", "brown") ## Same order that occurs in the source files 1=aq=lightblue, 2=Ter=brown

### setting output options, this is for plot_continous_var2() function. 
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

### Typing down The folder names we process
list.files()
Folders<-c("long_Mod with lake+species+PC1 w_0.23",
           "long_Mod with lake+species+For_area w_0.23",
           "long_Mod with lake+species+POM_NC_ratio w_0.23")

### With this loop, we are creating Article Fig 3 for different models (environmental factors)
### This produces Continous variable effect on consumer mixtures using plot_cont_var2() function from MixSIAR package (name in MixSIAR: plot_continous_var)
### And we supplement this figure with Lake and species specific posterior distributions
### In addition we make posterior density plots for for each lake and species
for (f in unique(Folders)){
    env<-sub(".*species\\+(.*)( w_).*", "\\1", f) ## getting the environmental variable name from the model file name
    wtot<-sub(".*w_ *(.*?)", "\\1", f) ## getting the used omega value from the model file name 
    setwd(file.path(orgfolder,c("3. MixSIAR Models"), f))
    list.files()

    ### uploading mixtures, sources, discrimination and jags object
    mix<-readRDS(file = paste0("mix_object",".rds"))
    source<-readRDS(file = paste0("source_object",".rds"))
    discr<-readRDS(file = paste0("discr_object",".rds"))
    jags<-readRDS(file = paste0("mixSIAR_model_Mod with lake+species+", env, " w_",wtot,".rds"))
    list.files()
    
    ### Running the plot cont var function
    plot_list2_temp<-plot_continuous_var2(jags.1=jags, mix=mix, source=source, output_options=output_options,
                                                 alphaCI = 0.05, exclude_sources_below = 0.001)
    
    ## lets extract only the median continous variable figures pos1 [[1]] and then modify the colours and add posterior figures
    plot_list2<-list()
    for (k in 1:length(plot_list2_temp)){
      title<-plot_list2_temp[[k]][[1]][["labels"]][["title"]]
      plot_list2[[title]]<-plot_list2_temp[[k]][[1]]+
        scale_colour_manual(values = my_colours1)+
        scale_fill_manual(values=my_colours1)
    }
    
    ### Plotting initial plot, just for error checking
    png(file = paste0(env, "_Continous variable impact.png"), res=300, units="in", 
        width = 4*3, 3*length(plot_list2)/4)
    print(cowplot::plot_grid(plotlist=plot_list2, ncol=4))
    dev.off()
    
    #### Modifying Article Fig 3 ######
    
    ## Extracting posterior densities to a dataframe
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
              rep(rep(mix[["data"]][[env]], each=number.sim), times=number.sources)) # times=2 indicates number of sources
    colnames(gg)<-c("Lake", "species", env)
    posterior<-cbind(posterior_df, gg)
    
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
    
    ### Extracting the environmental variable used in the mixSIAR model
    gg<-distinct(data.frame(Lake=(mix[["data"]][["Lake"]]), env= (mix[["data"]][[env]])))
    colnames(gg)<-c("Lake", env) ## Just making sure the column names are correct
    ## combining env_dataset to posterior summary table
    posterior_summary_env<-merge(posterior_summary, gg, by.x="Lake", by.y="Lake") 
    
    #### Extracting a legend from the plots
    legend<-get_legend(plot_list2[[1]]+labs(colour=c("Diet sources"), fill=c("Diet sources"))+ # Renaming the Key title
                         scale_fill_manual(values = my_colours1, labels=c('Aquatic', 'Terrestrial'), name = "Dietary source")+ ## These are alphabetically orderder by their original labels which are Aq and Ter
                         scale_colour_manual(values = my_colours1, labels=c('Aquatic', 'Terrestrial'), name = "Dietary source")+ 
                         labs(colour="Test", fill="test")+
                         theme(legend.position="right",
                               legend.box = "horizontal", 
                               legend.background = element_rect(linetype = 2, size = 0.5, colour = 1),
                               legend.title=element_text(),
                               legend.justification = NULL,
                               legend.margin=margin(c(5,5,5,5)),
                               plot.margin=unit(c(0,0,0,0), units="line")))
    # plot_grid(legend) # just cheking if the extraction worked
    
    ## Making another figure list, where we export updated figures
    plot_list3<-list()
    ## Modifying the plots, removing axis names, adding lake specific posterior densities etc ## readjusting panel backgrounds in each subplotgroups 
    for(i in names(plot_list2)){
      plot_list3[[plot_list2[[i]][["labels"]][["title"]]]]<-plot_list2[[i]]+
        coord_cartesian(xlim=c(min(gg[env]), max(gg[env])))+ 
        scale_fill_manual(values=my_colours1)+     scale_colour_manual(values=my_colours1)+
        theme(legend.position="none",
              axis.text.x=element_blank(), #remove x axis labels
              #axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              #axis.ticks.y=element_blank(),  #remove y axis ticks
              axis.title = element_blank(),
              plot.title=element_text(size=10))+
        geom_linerange(inherit.aes = FALSE, data=posterior_summary_env[which(posterior_summary_env$Source=="Ter" & posterior_summary_env$species==plot_list2[[i]][["labels"]][["title"]]) ,], 
                       aes_string(x=env, ymin="y0.025", ymax="y0.975"), show.legend = FALSE, colour=my_colours1[2], alpha=0.35)+
        geom_point(inherit.aes = FALSE, data=posterior_summary_env[which(posterior_summary_env$Source=="Ter" & posterior_summary_env$species==plot_list2[[i]][["labels"]][["title"]]) ,], 
                   aes_string(x=env, y="y0.5"), show.legend = FALSE, colour=my_colours1[2])+
        theme(plot.background = element_rect(fill="grey85", colour="grey85"), panel.background = element_rect(fill="grey100", colour = "grey100"))
      
    }
    ##names(plot_list3) ##  just cheking if the i loop worked
    
    ### Extracting figures to different groups
    ZPL_list<-plot_list3[c("Bulk_ZPL", "Chaoborus", "Cladocera", "Copepods")]
    BMI_list<-plot_list3[c( "Bulk_BMI_littoral",
                            "Asellus_littoral", "Chironomid_littoral",
                            "Bulk_BMI_profundal", "Chironomid_Profundal")]
    Fish_plankt<-plot_list3[c("Perch_small", "Smelt", "Vendace", "Bleak")]
    Fish_mixed<-plot_list3[c("Perch_medium", "Roach_small", "Roach_large")]
    Fish_benthivor<-plot_list3[c("Ruffe")]
    Fish_piscivorous<-plot_list3[c("Perch_large", "Pike")]

    ### Drawing lines between the plot groups,
    ### I have drawn this figure by hand to paper first, so I know where to add the lines
    
    ZPL_fig2<-cowplot::plot_grid(plotlist = ZPL_list, ncol = 1) + theme(panel.border = theme_border(type=c("right"), size=2, linetype = 1, colour = "white"))
    BMI_fig2<-cowplot::plot_grid(plotlist = BMI_list, ncol = 1) + theme(panel.border = theme_border(type=c("right", "left"), size=2, linetype = 1, colour = "white"))
    Fish_plankt_fig2<-cowplot::plot_grid(plotlist = Fish_plankt, ncol = 2)+ theme(panel.border = theme_border(type=c("left", "bottom", "top"), size=2, linetype = 1, colour = "white"))
    Fish_benthivor_fig2<-cowplot::plot_grid(plotlist = Fish_benthivor, ncol=1, labels=c("f)"),label_x = 0.75)+ theme(panel.border = theme_border(type=c("left", "top"), size=2, linetype = 1, colour = "white"))
    ### awkward, just arranging of the plots
    Fish_benthivor2<-list()
    Fish_benthivor2[["Ruffe"]]<-Fish_benthivor_fig2
    Fish_mixed<-append(Fish_mixed, Fish_benthivor2, after=2)
    ## Fish mixed and piscivorous
    Fish_mixed_fig2<-cowplot::plot_grid(plotlist = Fish_mixed, ncol=2)+ theme(panel.border = theme_border(type=c("left", "top"), size=2, linetype = 1, colour = "white")) 
    Fish_piscivorous_fig2<-cowplot::plot_grid(plotlist=Fish_piscivorous,ncol=2)+theme(panel.border = theme_border(type=c("left", "bottom"), size=2, linetype = 1, colour = "white"))
    
    ### Combining the different plots And producing the final combination
    Comb_fig2<-plot_grid(plot_grid(plotlist = list(legend, ZPL_fig2) , ncol=1, rel_heights = c(1, 4), labels=c("","a)"),label_x = 0.75),
                         plot_grid(plotlist = list(BMI_fig2), ncol=1, labels=c("b)"),label_x = 0.75),
                         plot_grid(plotlist = list(Fish_piscivorous_fig2, Fish_plankt_fig2, Fish_mixed_fig2), rel_heights = c(1,2,2), ncol=1,  labels=c("c)", "d)", "e)"),label_x = 0.85), 
                         rel_widths = c(1,1,2), ncol=3)
    

    #### Plotting Article figure 3, ###### 
    ### Nature quidelines 183 mm (double column) and the full depth of the page is 170 mm. max width=7.2in, height = 6.7in
    #this was the furthest I could get with R at this time, modify axis titles, add y- and x-axis tick labels, subplot group names etc small stuff in vector graphics program
    pdf(file=paste0("4.Figures/",env,"Continous variable plots with posterior densities.pdf"), 
        height = 7, width = 6.5) 
    if(grepl("PC1", env)){ ## for some reason grid.arrage did not like the ifs to be inside the function. It messed the figure ratios. Therefore I have now placed it like this. 
      grid.arrange(Comb_fig2, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
                  top= textGrob(f ,gp = gpar(fontsize = 11)),
                  bottom = textGrob(paste0("Forestry, DOC <--  ", env, " axis  --> Agriculture, nutrients"),gp = gpar(fontsize = 11)))
    } else {
      grid.arrange(Comb_fig2, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
                   top= textGrob(f ,gp = gpar(fontsize = 11)),
                   bottom = textGrob(env,gp = gpar(fontsize = 11)))
    }
    dev.off()
    
    
    #### Plotting supplementary plots i.e. Lake and species specific posterior densities ####
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
    pdf(file=paste0("Dietary contribution_by lake and species_", f,".pdf"))
    print(plot_list_temp) ## ordering the list alphabetically
    dev.off()
    
    ## Setting the working directory back to original
    setwd(orgfolder)
} # End of f loop  

#########################
### END OF THE SCRIPT ###
#########################