########################################################################
#  Peru CHAMBASA TDT paper
#  Brian J. Enquist
#  Load data packages from Allie Shenkin 
# Abundance weighted AZ Chem, LMA, and Photosytnehsis data
# Used to create Peru_Gradient_NPP_Merged6.csv 
#  4/20/16
########################################################################## 

load("/Users/brianjenquist/GitHub/R/Peru_Analyses/lma.RData")
load("/Users/brianjenquist/GitHub/R/Peru_Analyses/az_leaf_chem_moments.RData")
load("/Users/brianjenquist/GitHub/R/Peru_Analyses/photosyn_trait_moments.RData")

## view data tables - look in the Environment in R to see the names of each data file generated
lma
az_chem_moments
photosyn_trait_moments_list

## export the summary data tables
write.csv(lma , "Chambasa_lma_data.csv")
write.csv(az_chem_moments , "Chambasa_az_chem_moments_data.csv")
write.csv(photosyn_trait_moments_list , "Chambasa_photosyn_trait_moments_list.csv")
