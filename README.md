# Keva-et-al.-2025-NC
This is a data and script repository for Keva et al. article (Insert DOI here) "The role of land use in terrestrial support of boreal lake food webs "

# This repository includes (Folder structure):
1. Scripts # includes scripts for figures: isotope biplots, MixSIAR runs, Continous variable plots, Model comparisons, Omega sensitivity analysis
2. Data # Includes consumer omega corrected d2H data, source d2H data, tdf, lake PCA values, lake water d2H values
3. MixSIAR models # Includes the output folders from mixSIAR analysis
4. Figures # Includes manuscript figures produced with the given scripts

In each of the folders there are more detailed readme-files

# Packages needed
1. Data arrangement: plyr (version: 1.1.4), dplyr (1.1.0), reshape2 (1.1.4) 
2. Graphics: ggplot2 (3.4.1), ggnewscale (0.5.0), cowplot (1.1.3), gridExtra (2.3), splancs(2.01-45), ggrepel (0.9.6), ggpubr (0.6.0)
3. MixSIAR stuff: MixSIAR (3.1.12), rjags (4-16), R2WinBUGS (2.1-22.1)

Remember to install JAGS as well (e.g. 4.3.0)

# Mixiar jags model objects
The MixSiar jags model objects (output of run_model()-function) are stored in this repository through "Releases" feature.

# Licence

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
