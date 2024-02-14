# Keva-et-al.-2024---NC
This is a data and script repository for Keva et al. 2024 article (Insert DOI here) "Land use drives terrestrial support of boreal lake food webs"

This repository includes:
1. Hydrogen isotope values of consumers (omega corrected), sources, trophic discrimination file
2. R scripts used to run: MixSIAR, isotope biplots, omega sensitivity analysis figures

# Folder Structure
1. Scripts # includes scripts for figures: isotope biplots, MixSIAR runs, Continous variable plots, Model comparisons, Omega sensitivity analysis
2. Data # Includes consumer omega corrected d2H data, source d2H data, tdf, lake PCA values, lake water d2H values
3. MixSIAR models # Includes the output folders from mixSIAR analysis
4. Figures # Includes manuscript figures
   
# Packages needed
1. Data arrangement: plyr, dplyr, reshape2 
2. Graphics: ggplot2, ggnewscale, cowplot, gridExtra, splancs, ggrepel, ggpubr
3. MixSIAR stuff: MixSIAR, rjags, R2WinBUGS

Remember to install JAGS as well
