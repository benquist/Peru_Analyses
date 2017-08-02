
#==========================================================================
# Code set 2
#
# Assessing trait-based scaling theory in tropical forests spanning a broad temperature gradient 
# 
#Authors:
#Brian J. Enquist1,2*, Lisa Patrick Bentley3, Alexander Shenkin3, Brian Maitner1, Van Savage4, Sean Michaletz1,5, Benjamin Blonder3, Vanessa Buzzard1, Tatiana Erika Boza Espinoza6, William Farfan-Rios 7,8 Chris Doughty9, Gregory R. Goldsmith10, Roberta E. Martin11, Norma Salinas3,8, Miles Silman8, Sandra Díaz12, Gregory P. Asner11 & Yadvinder Malhi3
#
# Main analyses and graphs
# 8/2/17
#==========================================================================

library(stargazer)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged7.csv",header=T)

#==========================================================================
#Calculated Variables 
#==========================================================================

##Calculate Boltzmann 1/kT ##

Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.+273.15))

Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

##Calculate the leaf N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))

##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPP_newperNMeanMean <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$NMeanMean))

Peru_Plot_Master.data$GPP_new_estimateperNMeanMean <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$NMeanMean))

#* calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)/(Peru_Plot_Master.data$RLeaf))


##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$n_percent.sun.mu.abundance)/ (Peru_Plot_Master.data$p_corrected_percent.sun.mu.abundance))
Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))


#Calculate NPP_newLeaf/RLeaf - Production per carbon respired
Peru_Plot_Master.data$NPP_newLeafperRLeaf <- ((Peru_Plot_Master.data$NPP_newLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST prediction
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_GPP1Upper <- ((Peru_Plot_Master.data$GPP_Malhi_2017_UpperRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_GPP1Lower <- ((Peru_Plot_Master.data$GPP_Malhi_2017_LowerRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP1Upper <- ((Peru_Plot_Master.data$NPP_Malhi_2017_UpperRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP1Lower <- ((Peru_Plot_Master.data$NPP_Malhi_2017_LowerRange)/(Peru_Plot_Master.data$MST_AGB))

theme_set(bw(base_size = 10))
#change of units to match
TransformedNpercent <- (Peru_Plot_Master.data$NMeanMean * 100)
TransformedNpercentLower <- (Peru_Plot_Master.data$NMeanLower *100)
TransformedNpercentUpper <- (Peru_Plot_Master.data$NMeanUpper *100)
TransformedPpercent <- (Peru_Plot_Master.data$PMeanMean * 100)
TransformedPpercentLower <- (Peru_Plot_Master.data$PMeanLower *100)
TransformedPpercentUpper <- (Peru_Plot_Master.data$PMeanUpper *100)
TransformedCpercent <- (Peru_Plot_Master.data$CMeanMean * 100)
TransformedCpercentLower <- (Peru_Plot_Master.data$CMeanLower *100)
TransformedCpercentUpper <- (Peru_Plot_Master.data$CMeanUpper *100)
Transformedmean_sla_lamina_petiole <- (Peru_Plot_Master.data$mean_sla_lamina_petiole * 1000)
TransformedSLALower <- (Peru_Plot_Master.data$SLAMeanLower * 1000)
TransformedSLAUpper <- (Peru_Plot_Master.data$SLAMeanUpper * 1000)

#==========================================================================
# Manuscript Tables
#==========================================================================

####################################################################################
#### Data Summary - Site and Plot Summary - Table 
#### http://jakeruss.com/cheatsheets/stargazer.html
#### If the summary option is set to FALSE then stargazer will instead print the contents of the data frame.

library("dplyr")
library("stargazer")
Table1 <- subset(Peru_Plot_Master.data, select=c("Plot_Name", "Latitude", "Longitude", "Elevation.m.", "SolarRadiation.GJ.m.2.yr.1.", "MeanAnnualAirTemperature.degC.", "Precipitation.mm.yr.1.", "Aboveground_biomass2", "GPP_Malhi_2017", "NPP_Malhi_2017" ))
#Table1 <- Peru_Plot_Master.data %>% slice(1:6)

Table2 <- stargazer(Table1, type = "html", summary = FALSE, rownames = FALSE, title="Peru CHAMBASA Plots", column.labels = c("Plot Name", "Latitude", "Longitude", "Elevation (m)", "Solar Radiation(GJ m2 yr-1)", "Mean Annual Air Temperature (degC)","Precipitation mm yr-1", "Aboveground_biomass (kg)", "GPP", "NPP" ), digits = 2)

write.csv(Table2, "Table_Summary2.html")

#==========================================================================
#Shift in mean traits of the most abundant species

linear.1 <- lm(n_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(p_corrected_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
linear.3 <- lm(1/PlotNtoP ~ Elevation.m., data=Peru_Plot_Master.data)
linear.4 <- lm(c_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
linear.5 <- lm(lma.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
linear.6 <- lm(amax.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotosynthesisPerLeafN ~ Elevation.m., data=Peru_Plot_Master.data)


summary(linear.3)
####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table1 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Most Abunant Species)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "P:N", "%C", "LMA", "Photo", "PNUE" ), digits = 2)
#no.space=TRUE
write.csv(table1 , "Table_1.html")

summary(linear.3)



#==========================================================================
#Shift in subsampled whole-community trait distributions
#Shift in community sampled trait means


linear.1 <- lm(TransformedNpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(TransformedPpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.3 <- lm(1/(PlotNtoPMean) ~ Elevation.m., data=Peru_Plot_Master.data)
linear.4 <- lm(TransformedCpercent  ~ Elevation.m., data=Peru_Plot_Master.data)
linear.5 <- lm(1/(SLAMeanMean) ~ Elevation.m., data=Peru_Plot_Master.data)
linear.6 <- lm(PhotoMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotoPerLeafNMean ~ Elevation.m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table2 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Trait Distribution Subsampling)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "P:N", "%C", "LMA", "Photo", "PNUE" ), digits = 2)
#no.space=TRUE
write.csv(table2,"Table2.html")

summary(linear.5)

#==========================================================================
#Shift in trait moments with elevation
#

theme_set(bw(base_size = 10))

linear.1 <- lm(TransformedNpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(NVarianceMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.3 <- lm(NSkewnessMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.4 <- lm(NKurtosisMean ~ Elevation.m., data=Peru_Plot_Master.data)

table1N <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                     title="Trait Drivers Theory - Assessing shifts in trait moments", 
                     single.row=FALSE, ci=TRUE,ci.level=0.95, 
                     omit.stat=c("f", "ser"), ci.separator = " - ", 
                     dep.var.labels.include = FALSE,
                     covariate.labels = c("Elevation (m)"),
                     model.numbers= FALSE, 
                     column.labels = c("%N Mean", "%N Variance", "%N Skewness", "%N Kurtosis"), digits = 2)
write.csv(table1N,"Table1N.html")


linear.1 <- lm(TransformedPpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(PVarianceMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.3 <- lm(PSkewnessMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.4 <- lm(PKurtosisMean ~ Elevation.m., data=Peru_Plot_Master.data)

table1P <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                     title="Trait Drivers Theory - Assessing shifts in trait moments", 
                     single.row=FALSE, ci=TRUE,ci.level=0.95, 
                     omit.stat=c("f", "ser"), ci.separator = " - ", 
                     dep.var.labels.include = FALSE,
                     covariate.labels = c("Elevation (m)"),
                     model.numbers= FALSE, 
                     column.labels = c("%P Mean", "%P Variance", "%P Skewness", "%P Kurtosis"), digits = 2)

write.csv(table1P,"Table1P.html")


linear.1 <- lm(TransformedCpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(CVarianceMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.3 <- lm(CSkewnessMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.4 <- lm(CKurtosisMean ~ Elevation.m., data=Peru_Plot_Master.data)

table1C <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                     title="Trait Drivers Theory - Assessing shifts in trait moments", 
                     single.row=FALSE, ci=TRUE,ci.level=0.95, 
                     omit.stat=c("f", "ser"), ci.separator = " - ", 
                     dep.var.labels.include = FALSE,
                     covariate.labels = c("Elevation (m)"),
                     model.numbers= FALSE, 
                     column.labels = c("%C Mean", "%C Variance", "%C Skewness", "%C Kurtosis"), digits = 2)
write.csv(table1C,"Table1C.html")


linear.1 <- lm(PhotoMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(PhotoVarianceMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.3 <- lm(PhotoSkewnessMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.4 <- lm(PhotoKurtosisMean ~ Elevation.m., data=Peru_Plot_Master.data)

table1Photo <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                          title="Trait Drivers Theory - Assessing shifts in trait moments", 
                          single.row=FALSE, ci=TRUE,ci.level=0.95, 
                          omit.stat=c("f", "ser"), ci.separator = " - ", 
                          dep.var.labels.include = FALSE,
                          covariate.labels = c("Elevation (m)"),
                          model.numbers= FALSE, 
                          column.labels = c("Photo Mean", "Photo Variance", "Photo Skewness", "Photo Kurtosis"), digits = 2)
write.csv(table1Photo,"Table1Photo.html")

linear.1 <- lm(1/(SLAMeanLower) ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(SLAVarianceMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.3 <- lm(SLASkewnessMean ~ Elevation.m. , data=Peru_Plot_Master.data)
linear.4 <- lm(SLAKurtosisMean ~ Elevation.m., data=Peru_Plot_Master.data)

table1SLA <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                       title="Trait Drivers Theory - Assessing shifts in trait moments", 
                       single.row=FALSE, ci=TRUE,ci.level=0.95, 
                       omit.stat=c("f", "ser"), ci.separator = " - ", 
                       dep.var.labels.include = FALSE,
                       covariate.labels = c("Elevation (m)"),
                       model.numbers= FALSE, 
                       column.labels = c("LMA Mean", "LMA Variance", "LMA Skewness", "LMA Kurtosis"), 
                       digits = 2)

write.csv(table1SLA,"Table1SLA.html")


#==========================================================================
#Table 2 MST Boltzmann fits to traits of sampled individuals of abundant species
#

#Calculate mass corrected NPP and GPP
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))

#Linear.1 <- lm((PCA1ScoresPlotTraits_a) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.2 <- lm(log(1/PlotNtoP) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.3 <- lm(log(PhotosynthesisPerLeafN) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(amax.sun.mu.abundance) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.5 <- lm(log(n_percent.sun.mu.abundance) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.6 <- lm(log(lma.sun.mu.abundance) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.7 <- lm(log(MST_NPP_new1) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.8 <- lm(log(MST_GPP1) ~ MAinvBT, data=Peru_Plot_Master.data)  ## note! use of GPP_new_estimate rather than GPP_new_estimate estimated here. 
#Linear.6 <- lm(log(Aboveground_biomass) ~ MAinvBT, data=Peru_Plot_Master.data)

table2MST <- stargazer(Linear.2, Linear.3, Linear.4, Linear.5, Linear.6,Linear.7,Linear.8,
                    type="html", 
                    title="Metabolic Scaling Theory - Boltzmann Temperature Results (Most abundant species)", 
                    single.row=FALSE, 
                    ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), 
                    ci.separator = " - ", 
                    dep.var.labels.include = FALSE, 
                    model.numbers= FALSE, 
                    column.labels = c("ln(P:N)", "ln(N effciency)", "ln(photosynthesis)", "ln(N)", "ln(LMA)", "ln(Mass Corrected NPP)", "ln(Mass Corrected GPP"), digits = 2)
write.csv(table2MST , "Table2_MST_.html")

summary(Linear.1)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table4 <- stargazer(Linear.1, Linear.2, Linear.3, Linear.4, Linear.5, Linear.6, Linear.7, Linear.8, type="html", 
                    title="Metabolic Scaling Theory - Boltzmann Temperature Results", 
                    single.row=FALSE, 
                    ci=TRUE,
                    ci.level=0.95, 
                    omit.stat=c("f", "ser"), 
                    ci.separator = " - ", 
                    dep.var.labels.include = FALSE, 
                    model.numbers= FALSE, 
                    column.labels = c("ln(P:N)", "ln(N effciency)", "ln(photosynthesis)", "ln(%N)", "ln(LMA)", "ln(Biomass(kg))", "ln(NPP_new)", "ln(GPP_new_estimate)"), 
                    digits = 2)

write.csv(table4, "Table_MST_bivariate_abundant.html")



#==========================================================================
#MST Boltzmann fits to subsampled trait distributions
#

stargazer(Peru_Plot_Master.data)
#linear.1 <- lm(log(PCA1ScoresTraitSample) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.2 <- lm(log(1/PlotNtoPMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.3 <- lm(log(PhotoPerLeafNMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(PhotoMeanMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log(NMeanMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.6 <- lm(log(1/SLAMeanMean) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.7 <- lm(log(MST_NPP_new1) ~ MAinvBT, data=Peru_Plot_Master.data)
Linear.8 <- lm(log(MST_GPP1) ~ MAinvBT, data=Peru_Plot_Master.data)  ## note! use of GPP_new_estimate rather than GPP_new_estimate estimated here. 

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:

## note to add linear.8 put it in for one of the other linear.X inputs to the table. 
table4 <- stargazer(linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, linear.8,
                    type="html",
                    title="Metabolic Scaling Theory - Boltzmann Temperature Results (Trait Distribution Subsampling)", 
                    single.row=FALSE, 
                    ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), 
                    ci.separator = " - ", 
                    dep.var.labels.include = FALSE, 
                    model.numbers= FALSE, 
                    column.labels = c("ln(P:N)", "ln(N effciency)", "ln(photosynthesis)", "ln(N)", "ln(LMA)", "ln(Mass Corrected NPP)", "ln(Mass Corrected GPP"), digits = 2)
write.csv(table4 , "Table_MST_bivariate_community.html")

#==========================================================================
#Table_S5a   MST GPP_new_estimate and NPP_new scaling fits to subsampled trait distributions
#

linear.1 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.2 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.3 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)

table5a <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", 
                    title="Metabolic Scaling Theory - Ecosystem NPP and GPP scaling", 
                    single.row=FALSE, ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), ci.separator = " - ", 
                    dep.var.labels.include = FALSE,
                    covariate.labels = c("Aboveground Biomass (kg)", "Temperature 1/kT"),
                    model.numbers= FALSE, 
                    column.labels = c("GPP", "NPP", "NPP", "GPP", "NPP"),
                    digits = 2)
write.csv(table5a , "Table_5a.html")
#write.csv(table5a , "Table5a.doc")


#==========================================================================
#Table_S5b   MST GPP_new_estimate and NPP_new scaling fits to subsampled trait distributions
#but with PNUE
linear.1 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.2 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.3 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.6 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.7 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + + MAinvBT + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.8 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)

table5b <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, linear.8, type="html", 
                    title="Metabolic Scaling Theory - Ecosystem NPP and GPP scaling", 
                    single.row=FALSE, ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), ci.separator = " - ", 
                    dep.var.labels.include = FALSE,
                    covariate.labels = c("Aboveground Biomass (kg)", "Temperature 1/kT", "Plot Mean PNUE"),
                    model.numbers= FALSE, 
                    column.labels = c("NPP", "GPP", "NPP", "GPP", "NPP", "GPP", "NPP", "GPP"),
                    digits = 2)
write.csv(table5b , "TableS5_b.html")


