# This Folder includes following files # 

# Consumer data_w_(value).csv # values in c(0, 0.14, 0.23, 0.32). 
Consumer d2H isotope values used in MixSIAR and raw d2H values of the samples (w=0)

Column names and descriptions: 
1. Lake: Lake names for each sample
2. species: Species or group for each sample
3. d15N: consumer d15N values
4. TL: Moedelled trophic level of the consumers (base cladocera TP=2.1, TDF of d15N = 3.4 for ZPL, profundal benthos and fish; asellus TP 2.1 TDF 3.4 for Littoral Benthos)
5. For_area: catchment area agriculture area coverage 
6. Agr_area: catchment area agriculture area coverage 
7. pH: Lake water pH for each lake
8. POM_NC_ratio: Lake water POM N:C ratio for each lake
9. PC1: Principal component 1 values for each lake
10. PC2: Principal component 2 values for each lake
11. PC3: Principal component 3 values for each lake
12. d2H: omega corrected consumer d2H values

# Source data.csv, Source data_epi.csv, Source data_bacw17.csv #
Source d2H value dataframes used in MixSIAR

Column names and descriptions:  
1. Row names: Source identifier Aq and Ter
    1.1. Source data.csv : Aq = modelled phytoplankton d2H values, Ter=inlet DOM d2H values
    1.2. Source data_epi.csv : Aq = benthic algae d2H values, Ter=inlet DOM d2H values
    1.3. Source data_bacw17.csv : Aq = modelled phytoplankton d2H values, Ter= modelled bacterial d2H values
3. Lake: Lake names for each sample
4. Meand2H: d2H values for each lake and source
5. SDd2H: standard deviation for each sample. This is derived from spatiotemporal sampling of four lakes and applied to all lakes.

There are also Source data2.csv that includes modelled source values based on the Source data.csv. This Source data2.csv was used only when running the empty MixSIAR model (e.g. Table S6). 
    
# tef.csv ## 
Trophic discrimination factors for d2H used in MixSIAR

Column names and descriptions:
1. Row names: Source identifier Aq (benthic algae) and Ter (inlet DOM) 
2. meand2H: TDF of d2H, set to 0 
 SDd2H: TDF SD is set to 13‰: To account variability in trophic level base line (cladocera TL=2.1±0.1 & asellus TL=2.1±0.1) (Tanentzap et al., 2017) and TDF of nitrogen (3.4±1.0‰) (Post 2002), 
 we added additional uncertainty to consumer d2H values (±13‰, based on mean difference of above-mentioned scenarios), following a recent study (Vane et al., 2023)

# Lake water d2H and source values.csv # 
data file including the lake water d2H values, benthic algae d2H values, inlet DOM d2H values, modelled bacterial and phytoplankton d2H values

Column names and descriptions:  
1. Lake: Lake names for each sample
2. d2H_water_center: d2H values of the analysed lake water
3. d2H_benthic_algae: d2H values of the analysed benthic algae samples
4. d2H_DOM: d2H values of the analysed inlet DOM samples
5. d2H_bact_w017_modelled: the modelled bacterial d2H values with omega value of 0.17 (Based on Fogela et al., 2016 table 1 mean value)
6. d2H_phyto_modelled: the modelled algal d2H values with fractination value of 109.7 (based on the mean difference between water and benthic algae d2H values)

# lake_PCA_scores.csv # 
data file including the lake water d2H values used in isotope Biplots script

Column names and descriptions:  
1. Lake.name: Lake names 
2. Abbreviation: Abbreviated lake names
3. PC1 = Principal component 1 values for each lake
4. PC2 = Principal component 2 values for each lake
5. PC3 = Principal component 3 values for each lake
