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
Folders<-c("long_Mod with lake+species+PC1 w_0.23",
           "long_Mod with lake+species+PC1 w_0.14",
           "long_Mod with lake+species+PC1 w_0.32")

#Creating a list for data frames
df_list<-list()

## Running a loop that extract linear coefficient impact on dietary proportions to df_list
for (f in unique(Folders)){
  env<-sub(".*species\\+(.*)( w_).*", "\\1", f) ## getting the environmental variable name from the model file name
  wtot<-sub(".*w_ *(.*?)", "\\1", f) ## getting the used omega value from the model file name 
  setwd(file.path(orgfolder,c("3. MixSIAR Models"), f))
  list.files()
  
  ### uploading mixtures, sources, discrimination and jags object
  mix<-readRDS(file = paste0("mix_object",".rds"))
  source<-readRDS(file = paste0("source_object",".rds"))
  discr<-readRDS(file = paste0("discr_object",".rds"))
  jags<-readRDS(file = paste0("mixSIAR_model_Mod with lake+species+", env, " w_",wtot, ".rds"))
  #list.files()
  
  ### This is copied from plot_cont_var() function from MixSIAR package
  #### I had to do this as my intention is to only show lines for one source + modifying the colours according to omega values
  R2jags::attach.jags(jags)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  alphaCI=0.05
  exclude_sources_below=0.001
  
  for(ce in 1:mix$n.ce){
    if(mix$n.fe == 1){ 
      g <- vector("list", length = mix$FAC[[1]]$levels)
      for(f1 in 1:mix$FAC[[1]]$levels){
        g[[f1]] <- vector("list", 4)
        fac.lab <- mix$FAC[[1]]$labels[f1]
        label <- mix$cont_effects[ce]
        cont <- mix$CE[[ce]]
        ilr.cont <- get(paste("ilr.cont",ce,sep=""))
        
        n.plot = 200
        chain.len = dim(p.global)[1]
        Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
        ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
        for(src in 1:n.sources-1){
          for(i in 1:n.plot){
            ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i] + ilr.fac1[,f1,src]
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
              cross[i,d,,j] <- (e[,j]^ilr.plot[i,j,d])/sum(e[,j]^ilr.plot[i,j,d]);
            }
            for(src in 1:n.sources){
              tmp[i,d,src] <- prod(cross[i,d,src,]);
            }
            for(src in 1:n.sources){
              p.plot[i,d,src] <- tmp[i,d,src]/sum(tmp[i,d,]);
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
        df$omega<-wtot
        df_list[[fac.lab]][[paste0("omega=", wtot)]]<-df
      
      } # End of f1 loop
    } ## End of if statement
  } # End of ce loop
  setwd(orgfolder)
} # End of f loop 

## Then we extract the cont effect dataframes from the df_list, and make the plots
plot_list<-list()
for(sp in names(df_list)){
  p_df<-do.call("rbind", list(df_list[[sp]][["omega=0.14"]], 
                              df_list[[sp]][["omega=0.23"]], 
                              df_list[[sp]][["omega=0.32"]])) 
  omega_order<-c("0.14","0.32", "0.23")
  my_colours<-c("#EF7F2A", "#A79ECD", "brown")
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

## Making another figure list, where we export updated figures
plot_list3<-list()
## Modifying the plots, removing axis names, adding lake specific posterior densities etc ## readjusting panel backgrounds in each subplotgroups 
for(i in names(plot_list)){
  plot_list3[[i]]<-plot_list[[i]]+
    theme(legend.position="none",
          axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          #axis.ticks.y=element_blank(),  #remove y axis ticks
          axis.title = element_blank(),
          plot.title=element_text(size=10))+
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

### Drawing lines between the plotgroups,
### I have drawn this figure by hand to paper first, so I know where to add the lines
setwd(file.path(orgfolder, "1. Scripts/Help functions"))
source("Theme_functions.R")
setwd(orgfolder)

ZPL_fig2<-cowplot::plot_grid(plotlist = ZPL_list, ncol = 1) + theme(panel.border = theme_border(type=c("right"), size=2, linetype = 1, colour = "white"))
BMI_fig2<-cowplot::plot_grid(plotlist = BMI_list, ncol = 1) + theme(panel.border = theme_border(type=c("right", "left"), size=2, linetype = 1, colour = "white"))
Fish_plankt_fig2<-cowplot::plot_grid(plotlist = Fish_plankt, ncol = 2)+ theme(panel.border = theme_border(type=c("left", "bottom", "top"), size=2, linetype = 1, colour = "white"))
Fish_benthivor_fig2<-cowplot::plot_grid(plotlist = Fish_benthivor, ncol=1, labels=c("f)"),label_x = 0.75)+ theme(panel.border = theme_border(type=c("left",  "top"), size=2, linetype = 1, colour = "white"))
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

#### Plotting the omega sensitivy analysis, ###### 
### Nature quidelines 183 mm (double column) and the full depth of the page is 170 mm. max width=7.2in, height = 6.7in
#this was the furthest I could get with R at this time, modify axis titles, add y- and x-axis tick labels, subplot group names etc small stuff in vectorgraphics program
pdf(file="4. Figure/Omega Sensitivity analysis.pdf", 
    height = 7, width = 6.5) 
grid.arrange(Comb_fig2, left = textGrob("Dietary contribution",rot=90, gp = gpar(fontsize = 11)), 
            top= textGrob("Omega Sensitivity analysis" ,gp = gpar(fontsize = 11)),
            bottom = textGrob(paste0("Forestry, DOC <--  PC1 axis  --> Agriculture, nutrients"),gp = gpar(fontsize = 11)))
dev.off()

#########################
### END OF THE SCRIPT ###
#########################