# This Folder includes following files # 

# Consumer data_w_(value).csv) # values in c(0, 0.14, 0.23, 0.32). 
Consumer d2H isotope values used in MixSIAR

Column names and descriptions: 
Lake: Lake names for each sample
species: Species or group for each sample
d15N: consumer d15N values
TL: Moedelled trophic level of the consumers (base cladocera TP=2.0, TDF of d15N = 3.4)
For_area: catchment area agriculture area coverage 
Agr_area: catchment area agriculture area coverage 
pH: Lake water pH for each lake
POM_NC_ratio: Lake water POM N:C ratio for each lake
PC1: Principal component 1 values for each lake
PC2: Principal component 2 values for each lake
PC3: Principal component 3 values for each lake
d2H: omega corrected consumer d2H values

# Source data.csv #
  Source d2H value dataframe used in MixSIAR

Column names and descriptions:  
Row names: Source identifier Aq (benthic algae) and Ter (inlet DOM) 
Lake: Lake names for each sample
Meand2H: d2H values for each lake and source
SDd2H: standard deviation for each sample. This is derived from spatiotemporal sampling of four lakes and applied to all lakes.
    
# tef.csv ## 
 Trophic discrimination factors for d2H used in MixSIAR

 Column names and descriptions:
 Row names: Source identifier Aq (benthic algae) and Ter (inlet DOM) 
 meand2H: TDF of d2H, set to 0 
 SDd2H: TDF SD is set to 13‰: To account variability in trophic level base line (cladocera TL=2.0±0.1) (Tanentzap et al., 2017) and TDF of nitrogen (3.4±1.0‰) (Post 2002), 
 we added additional uncertainty to consumer d2H values (±13‰, based on mean difference of above-mentioned scenarios), following a recent study (Vane et al., 2023)

# Lake water d2H.csv # 
data file including the lake water d2H values

Column names and descriptions:  
Lake: Lake names for each sample
d2H_water_center: d2H values of the analysed lake water

# lake_PCA_scores.csv # 
data file including the lake water d2H values used in isotope Biplots script

Column names and descriptions:  
Lake.name: Lake names 
Abbreviation: Abbreviated lake names
PC1 = Principal component 1 values for each lake
PC2 = Principal component 2 values for each lake
PC3 = Principal component 3 values for each lake

