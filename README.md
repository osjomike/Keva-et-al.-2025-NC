# Keva-et-al.-2024---NC
This is a data and script repository for Keva et al. 2024 article (Insert DOI here) "Land use drives terrestrial support of boreal lake food webs"

# This repository includes (Folder structure):
1. Scripts # includes scripts for figures: isotope biplots, MixSIAR runs, Continous variable plots, Model comparisons, Omega sensitivity analysis
2. Data # Includes consumer omega corrected d2H data, source d2H data, tdf, lake PCA values, lake water d2H values
3. MixSIAR models # Includes the output folders from mixSIAR analysis
4. Figures # Includes manuscript figures produced with the given scripts
   
# Packages needed
1. Data arrangement: plyr, dplyr, reshape2 
2. Graphics: ggplot2, ggnewscale, cowplot, gridExtra, splancs, ggrepel, ggpubr
3. MixSIAR stuff: MixSIAR, rjags, R2WinBUGS

Remember to install JAGS as well

# Jags objects
Jags objects are not stored in this repository as they are too large files (200MB each). One can retreive them from releases linked to this repository once this repository is made public. Just copy the large .rds objects to the model folders. OR alternatively rerun the MixSIAR code. 

For peer reviewers: unfortunately the github releases are not shown in GitFront mirrored repository. If you really need the large jags models, just rerun the mixSIAR code. The original GitHub repository will be made public after the manuscript has been accepted.
