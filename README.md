# Keva-et-al.-2024---NC
This is a data and script repository for Keva et al. 2024 article (Insert DOI here) "Land use drives terrestrial support of boreal lake food webs"

# This repository includes (Folder structure):
1. Scripts # includes scripts for figures: isotope biplots, MixSIAR runs, Continous variable plots, Model comparisons, Omega sensitivity analysis
2. Data # Includes consumer omega corrected d2H data, source d2H data, tdf, lake PCA values, lake water d2H values
3. MixSIAR models # Includes the output folders from mixSIAR analysis
4. Figures # Includes manuscript figures produced with the given scripts

In each of the folders there are more detailed readme-files

# Packages needed
1. Data arrangement: plyr, dplyr (version: 1.1.0), reshape2 
2. Graphics: ggplot2, ggnewscale, cowplot, gridExtra, splancs, ggrepel, ggpubr
3. MixSIAR stuff: MixSIAR (version: ), rjags, R2WinBUGS

Remember to install JAGS as well (e.g. 4.3.0)

# Mixiar jags model objects
The Mixiar jags model objects (output of run_model()-function) are stored in this repository through "Releases" feature.

# Licence

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
