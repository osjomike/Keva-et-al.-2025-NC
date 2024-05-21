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
consumers<-read.csv2(file = paste0(Data_dir, "/Consumerdata_w_0.23newTP.csv"), sep = ",", dec=".")
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
         title = i)+
    coord_cartesian(ylim=c(-250,-50))+ #, xlim=c(6,-4))+
    theme_classic()+
    theme(legend.position="none",
          axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          #axis.ticks.y=element_blank(),  #remove y axis ticks
          axis.title = element_blank(),
          plot.title=element_text(size=10),
          plot.background = element_blank(), panel.background = element_rect(fill="grey100", colour = "grey100"))
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
legend1<- plot_grid(legend) # Checking the legend

fig_order<-c("Bulk_BMI_littoral", "Perch_large", "Pike", 
             "Bulk_ZPL", "Asellus_littoral", "Perch_small", "Smelt",
             "Chaoborus", "Chironomid_littoral", "Vendace", "Bleak", 
             "Cladocera", "Bulk_BMI_profundal", "Perch_medium", "Roach_small",
             "Copepods", "Chironomid_Profundal", "Ruffe", "Roach_large")

new_figlist1<-list()
new_figlist1<-fig_list[fig_order]
new_figlist1<- append(new_figlist1,list(legend1), after = 0)
comb_fig<-cowplot::plot_grid(plotlist = new_figlist1, ncol=4)

### Finally doing the plot and adding axis legends. 
### Add the axis labels in incscape or similar, check labels e.g. here fig_list[[1]]

pdf(file="4. Figures/isotope biplots.pdf", 
    height = 7, width = 6.5)
grid.arrange(comb_fig, left = textGrob(expression(paste(delta^{2}, "H (\u2030)"["cor, w=0.23"])),rot=90, gp = gpar(fontsize = 11)), 
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