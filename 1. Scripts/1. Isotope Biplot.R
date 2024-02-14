rm(list = ls()) # clear R's clobal envirnoment memory

## Author: Ossi Keva
## Date 14.2.2024
## Project: Roger I Jones -- Allocarb, Antti Eloranta -- FreshRestore

## Short description: With this script we are plotting the "isotope biplots". Article Fig 2.

## Some packages to download
library(dplyr) ## pipe operatore
library(reshape2) ## melt and dcast
library(ggplot2) ## plotting functions
library(ggnewscale) ## some colour layer functions
library(cowplot) ## plotting functions
library(gridExtra) ## plotting functions

## Setting working directory
setwd("D:\\Keva et al. NC 2024") ## Set correct wd here
orgfolder<-getwd()
## Patch for data
Data_dir<-file.path(orgfolder,"2. Data")

## downloading some data sets
consumers<-read.csv2(file = paste0(Data_dir, "/Consumer data_w_0.23.csv"), sep = ",", dec=".")
sources<-read.csv2(file = paste0(Data_dir, "/Source data.csv"), sep = ",", dec=".")
Lake_PCA<-read.csv2(file = paste0(Data_dir, "/lake_PCA_scores.csv"), sep = ",", dec=".", fileEncoding = "latin1")
water<-read.csv2(file = paste0(Data_dir, "/Lake water d2H.csv"), sep = ";", dec=".")  

# Checking column names
#colnames(consumers)
#colnames(sources)
#colnames(Lake_PCA)
#colnames(water)
## Modifying column names to match
colnames(sources)[1]<-"Source"
colnames(Lake_PCA)[2]<-"Lake"
## Changing the Scandic letters from lake names to ä->a, Ä->A, ö->o
Lake_PCA$Lake<-chartr("äÄö", "aAo", Lake_PCA$Lake)
water$Lake<-chartr("äÄö", "aAo", water$Lake)

## Renaming some variables
water$Source<-"water"
water$Meand2H<-water$"d2H_water_center"
water$SDd2H<-0 ## Setting water d2H SD to 0

## Subsetting water dataframe the files has other lakes as well.
lakes<-unique(Lake_PCA$Lake)
water_sub<-water[which(water$Lake %in% lakes),]

#colnames(sources)
source_water<-dplyr::bind_rows(sources[c("Lake", "Source","Meand2H", "SDd2H")], 
                               water_sub[c("Lake", "Source","Meand2H", "SDd2H")])

source_water_PC<-merge(x=source_water, y=Lake_PCA, by="Lake")


### OK THEN WE NEED TO DO TL LABELS
labels_n<-consumers %>%
  group_by(species) %>%
  summarise(n=sum(!is.na(d2H)), n_lakes=length(unique(Lake)), TL_=mean(TL, na.rm=TRUE), TL_SD=sd(TL, na.rm=TRUE))

### Creating a figure list and doing the plots
fig_list<-list()
for (i in unique(consumers$species)){
  fig_list[[i]]<-ggplot(consumers[which(consumers$species ==i),], aes(x=PC1, y=d2H))+
    geom_point(data=source_water_PC, aes(x=PC1, y=Meand2H, colour=Source), size=1.2, shape=4, inherit.aes = FALSE)+
    geom_linerange(data=source_water_PC, aes(x=PC1, ymin=Meand2H-SDd2H, ymax=Meand2H+SDd2H, group=Source, colour=Source), inherit.aes =FALSE) +
    scale_colour_manual("Key for colours", values=c("purple","brown","lightblue","black"),
                        limits=c("water","Ter","Aq","consumer"), ## one can modify the order of the legend with this
                        labels=c("Water","Allochtonous","Autochtonous",  "Consumer"))+

    new_scale_color() + ## This is from ggnewscale, one can add several colour palettes to different datasets using this, but here we are using only one colour still
    geom_point(aes(y=d2H, x=PC1),shape=10, colour="black")+
    geom_smooth(aes(y=d2H, x=PC1), show.legend=FALSE, colour="black")+
    #geom_smooth(method="lm",colour="black",show.legend = FALSE)+
    #scale_colour_manual(values = my_colours)+
    #scale_shape_manual(values = my_shapes)+
    geom_text(data=labels_n[which(labels_n$species ==i),], x=4.5,y=-250, hjust=1,vjust=0, size=2,
              aes(label=paste0("n=",n,"; lakes=", n_lakes, 
                               "; TL=", sprintf(TL_, fmt = '%#.1f'), "±", sprintf(TL_SD, fmt = '%#.1f'))))+
    
    labs(y=expression(paste(delta^{2}, "H (\u2030)"["cor, w=0.23"])),
         x=c("Forestry, DOC <------ PC1 ------> Agriculture, nutrients"),
         title = i)
}
## fig_list[[1]] just checking random figure from list, that the loop worked

### Getting the legend
legend<-get_legend(fig_list[[1]]+theme(legend.position="right",
                                       legend.box = "horizontal", 
                                       legend.background = element_rect(linetype = 2, size = 0.5, colour = 1),
                                       legend.title=element_text(),
                                       legend.justification = NULL,
                                       legend.margin=margin(c(5,5,5,5)),
                                       plot.margin=unit(c(0,0,0,0), units="line")))
# plot_grid(legend) # Checking the legend

## Uploading some help function to modify ggplot item borders
setwd(file.path(orgfolder, c("1. Scripts"), c("Help functions")))
source("Theme_functions.r")
setwd(orgfolder)
## Modifying the plots, removing axis names etc
plot_list3<-list()
for(i in names(fig_list)){
  plot_list3[[i]]<-fig_list[[i]]+
    #coord_cartesian(xlim=c(6, -3))+
    coord_cartesian(ylim=c(-250,-50))+ #, xlim=c(6,-4))+
    theme_classic()+
    theme(legend.position="none",
          axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          #axis.ticks.y=element_blank(),  #remove y axis ticks
          axis.title = element_blank(),
          plot.title=element_text(size=10),
          plot.background = element_rect(fill="grey85", colour="grey85"), panel.background = element_rect(fill="grey100", colour = "grey100"))
}

### Extracting figures to different groups
ZPL_list<-plot_list3[c("Bulk_ZPL", "Chaoborus", "Cladocera", "Copepods")]
BMI_list<-plot_list3[c("Bulk_BMI_littoral",
                      "Asellus_littoral", "Chironomid_littoral",
                      "Bulk_BMI_profundal", "Chironomid_Profundal")]
Fish_plankt<-plot_list3[c("Perch_small", "Smelt", "Vendace", "Bleak")]
Fish_mixed<-plot_list3[c("Perch_medium", "Roach_small", "Roach_large")]
Fish_benthivor<-plot_list3[c("Ruffe")]
Fish_piscivorous<-plot_list3[c("Perch_large", "Pike")]

### Creating plot groups, and drawing lines between the plot groups
ZPL_fig2<-cowplot::plot_grid(plotlist = ZPL_list, ncol = 1) + theme(panel.border = theme_border(type=c("right"), size=1, linetype = 1, colour = "white"))
BMI_fig2<-cowplot::plot_grid(plotlist = BMI_list, ncol = 1) + theme(panel.border = theme_border(type=c("right", "left"), size=2, linetype = 1, colour = "white"))
Fish_plankt_fig2<-cowplot::plot_grid(plotlist = Fish_plankt, ncol = 2)+ theme(panel.border = theme_border(type=c("left", "bottom", "top"), size=2, linetype = 1, colour = "white"))
Fish_benthivor_fig2<-cowplot::plot_grid(plotlist = Fish_benthivor, ncol=1, labels=c("f)"),label_x = 0.75)+theme(panel.border = theme_border(type=c("left", "right", "top"), size=2, linetype = 1, colour = "white"))
### awkward arranging of the plots
Fish_benthivor2<-list()
Fish_benthivor2[["Ruffe"]]<-Fish_benthivor_fig2
Fish_mixed<-append(Fish_mixed, Fish_benthivor2, after=2)
## Fish mixed and piscivorous
Fish_mixed_fig2<-cowplot::plot_grid(plotlist = Fish_mixed, ncol=2)+ theme(panel.border = theme_border(type=c("left", "top"), size=2, linetype = 1, colour = "white")) 
Fish_piscivorous_fig2<-cowplot::plot_grid(plotlist=Fish_piscivorous,ncol=2)+           theme(panel.border = theme_border(type=c("left", "bottom"), size=2, linetype = 1, colour = "white"))

### Combining the different plots and adjusting some widths (rel_widths) to make the left most figures similar sized compared to the rest of the plots
Comb_fig2<-plot_grid(plot_grid(plotlist = list(legend, ZPL_fig2) , ncol=1, rel_heights = c(1, 4), labels=c("","a)"),label_x = 0.75),
                     plot_grid(plotlist = list(BMI_fig2), ncol=1, labels=c("b)"),label_x = 0.75),
                     plot_grid(plotlist = list(Fish_piscivorous_fig2, Fish_plankt_fig2, Fish_mixed_fig2), rel_heights = c(1,2,2), ncol=1,  labels=c("c)", "d)", "e)"),label_x = 0.85), 
                     rel_widths = c(1,1,2), ncol=3)

### Finally doing the plot and adding axis legends. 
### Add the axis labels in incscape or similar, check labels e.g. here fig_list[[1]]

pdf(file="4. Figures/isotope biplots.pdf", 
    height = 7, width = 6.5)
grid.arrange(Comb_fig2, left = textGrob(expression(paste(delta^{2}, "H (\u2030)"["cor, w=0.23"])),rot=90, gp = gpar(fontsize = 11)), 
             bottom = textGrob("Forestry, DOC<--  PC1 axis  --> Agriculture, nutrients",gp = gpar(fontsize = 11)))
dev.off()

### Checking the samples outside the mixing envelope ####
mean<-dcast(sources, Lake~Source, value.var=c("Meand2H"))
sd<-dcast(sources, Lake~Source, value.var="SDd2H")

colnames(mean)<-c("Lake", "Aq_mean", "Ter_mean")
colnames(sd)<-c("Lake", "Aq_sd", "Ter_sd")
ss<-merge(mean, sd, by="Lake")
consumers_sources<-merge(consumers, ss, by="Lake")

## Samples outside source envelope source+-2SD
consumers_sources$diff_aq<-consumers_sources$d2H-(consumers_sources$Aq_mean-2*consumers_sources$Aq_sd)
consumers_sources$diff_Ter<-consumers_sources$d2H-(consumers_sources$Ter_mean+2*consumers_sources$Ter_sd)

errors<-consumers_sources %>%
  group_by(Lake, species) %>%
  #summarise_at(vars(sd), funs(grand_SD=grand.mean(.,N=n)))
  summarise(Diff_TER=sum(diff_Ter>0), Diff_aq=sum(diff_aq<0))

#If you are interested which consumers in which lakes were outside the envelope
x<-dcast(errors, Lake~species, value.var = "Diff_aq")
y<-dcast(errors, Lake~species, value.var = "Diff_TER")

### Checking how many samples were outside the envelope their percentage from total number of consumers
sum(errors$Diff_TER)+sum(errors$Diff_aq)
(sum(errors$Diff_TER)+sum(errors$Diff_aq))/length(consumers$d2H)

#########################
### END OF THE SCRIPT ###
#########################