# THIS FOLDER INCLUDES # 
Folders for MixSIAR runs with different covariates (PC1, Forest % or POM NC ratio) and mixSIAR runs with PC1 as covariate and omega values of 0.14 and 0.32. 
The folders and their contents are created with MixSIAR running script (in folder 1. Scripts). Each of the folders contains following objects:

1. diagnostics.txt # Run diagnostics including e.g. Gweke statistics
2. summary.txt # Summary of posterior results including e.g. epsilon values
3. discr.object.rds # Discrimination object used in mixSIAR
4. mix_object.rds # Consumer mixture object used in mixSIAR
5. source.object.rds # Consumer mixture object used in mixSIAR
6. mixSIAR_model_Mod with lake+species+(covariate) w_(value).rds # jags-object including posterior results
7. MixSIAR_modelMod with lake+species w_(value).txt # JAGS Model description (ONLY FOR Species or Species+Lake Models)

FOR Models including a continous covariate the JAGS model .txt file is in the 3. MoxSIAR Models folder roots.
THIS is a custom JAGS script that includes Random Slope for factor 1. in our case random slope for each species. NOTE That it works only with the given model structure and is not a dynamic object.  