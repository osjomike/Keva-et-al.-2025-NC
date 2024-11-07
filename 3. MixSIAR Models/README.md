# THIS FOLDER INCLUDES # 
Subfolders for MixSIAR runs with different covariates (PC1, Forest % or POM NC ratio), mixSIAR runs with PC1 as covariate and omega values of 0.14 and 0.32, mixSIAR runs without covariate, and empty model. 

The subfolders and their contents are created with MixSIAR running script (in folder 1. Scripts). 

Each of the subfolders contains following objects:
1. diagnostics.txt # Run diagnostics including e.g. Gweke statistics
2. summary.txt # Summary of posterior results including e.g. epsilon values
3. discr.object.rds # Discrimination object used in mixSIAR
4. mix_object.rds # Consumer mixture object used in mixSIAR
5. source.object.rds # Consumer mixture object used in mixSIAR
6. MixSIAR_modelMod with lake+species w_(value).txt # JAGS Model description (ONLY FOR Species or Species+Lake Models)

FOR Models including a continous covariate the JAGS model .txt file is in the 3. MixSIAR Models folder roots.
THIS is a custom JAGS script that includes Random Slope for factor 1. in our case random slope for each species. NOTE That it works only with the given model structure and is not a dynamic object.
