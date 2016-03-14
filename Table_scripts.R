#############################################
#####   CHAMBASA Tables. 
#####   B. J. Enquist
#############################################

library(stargazer)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged5.csv",header=T)

######### Calculated Variables ###############
##Calculate Boltzmann 1/kT
Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.+273.15))

##Calculate the leaf N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$mean_n_percent))

Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))

##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPP_newperNMeanMean <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$NMeanMean))

Peru_Plot_Master.data$GPP_new_estimateperNMeanMean <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$NMeanMean))


#* calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))


##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$mean_n_percent)/ (Peru_Plot_Master.data$mean_p_percent))
Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

#Calculate NPP_newLeaf/RLeaf - Production per carbon respired
Peru_Plot_Master.data$NPP_newLeafperRLeaf <- ((Peru_Plot_Master.data$NPP_newLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST prediction
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)

Peru_Plot_Master.data$MST_GPP_new_estimate1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))

##################################
## assessing goodness of trait sampling vs. mean of abundant species 1:1?

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

####################################################################################
####Table_1  Shift in mean traits of the most abundant species
######

linear.1 <- lm(mean_n_percent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(mean_p_percent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.3 <- lm(PlotNtoP ~ Elevation.m., data=Peru_Plot_Master.data)
linear.4 <- lm(mean_c_percent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.5 <- lm(mean_sla_lamina_petiole ~ Elevation.m., data=Peru_Plot_Master.data)
linear.6 <- lm(mean_photosynthesis ~ Elevation.m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotosynthesisPerLeafN ~ Elevation.m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table1 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Most Abunant Species)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "N:P", "%C", "SLA", "Photo", "PNUE" ), digits = 2)
#no.space=TRUE
write.csv(table1 , "Table1.html")

#########################################################################
## Shift in subsampled whole-community trait distributions
####################################################################################
# Table 2 - shift in community sampled trait means
#change of units to match

linear.1 <- lm(TransformedNpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(TransformedPpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.3 <- lm(PlotNtoPMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.4 <- lm(TransformedCpercent  ~ Elevation.m., data=Peru_Plot_Master.data)
linear.5 <- lm(SLAMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.6 <- lm(PhotoMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotoPerLeafNMean ~ Elevation.m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table2 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Trait Distribution Subsampling)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "N:P", "%C", "SLA", "Photo", "PNUE" ), digits = 2)
#no.space=TRUE
write.csv(table2,"Table2.html")



####################################################################################
####Table_S1   TDT shift in trait moments with elevation
######

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



linear.1 <- lm(SLAMeanLower ~ Elevation.m., data=Peru_Plot_Master.data)
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
                       column.labels = c("SLA Mean", "SLA Variance", "SLA Skewness", "SLA Kurtosis"), digits = 2)

write.csv(table1SLA,"Table1SLA.html")






####################################################################################
####Table_S4 MST Boltzmann fits to traits of sampled individuals of abundant species
######

linear.1 <- lm(log(PlotNtoP) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.2 <- lm(log(PhotosynthesisPerLeafN) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.3 <- lm(log(mean_photosynthesis) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(TransformedNpercent) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log(Aboveground_biomass) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.6 <- lm(log(NPP_new) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.7 <- lm(log(GPP_new_estimate) ~ MAinvBT, data=Peru_Plot_Master.data)  ## note! use of GPP_new_estimate rather than GPP_new_estimate estimated here. 


####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table2 <- stargazer(linear.1, linear.2, linear.3, linear.4,linear.5, linear.6, linear.7, type="html", title="Metabolic Scaling Theory - Boltzmann Temperature Results", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("ln(N:P)", "ln(N effciency)", "ln(photosynthesis)", "ln(%N)", "ln(Biomass(kg))", "ln(NPP_new)", "ln(GPP_new_estimate)"), digits = 2)

write.csv(table2 , "Table2xx.html")

####################################################################################
####Table
######

stargazer(Peru_Plot_Master.data$elevation_m)
linear.1 <- lm(TransformedNpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.2 <- lm(TransformedPpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.3 <- lm(PlotNtoPMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.4 <- lm(TransformedCpercent ~ Elevation.m., data=Peru_Plot_Master.data)
linear.5 <- lm(Transformedmean_sla_lamina_petiole ~ Elevation.m., data=Peru_Plot_Master.data)
linear.6 <- lm(PhotoMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotoPerLeafNMean ~ Elevation.m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table3 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Trait Distribution Subsampling)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "N:P", "%C", "SLA", "Photo", "N efficiency" ), digits = 2)
#no.space=TRUE
write.csv(table3 , "Table3.html")


####################################################################################
####Table_S4   MST Boltzmann fits to subsampled trait distributions
######

stargazer(Peru_Plot_Master.data)
linear.1 <- lm(log(PlotNtoPMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.2 <- lm(log(PhotoPerLeafNMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.3 <- lm(log(PhotoMeanMean) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(TransformedNpercent) ~ MAinvBT, data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table4 <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", title="Metabolic Scaling Theory - Boltzmann Temperature Results (Trait Distribution Subsampling)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("ln(N:P)", "ln(N effciency)", "ln(photosynthesis)", "ln(%N)"), digits = 2)
write.csv(table4 , "Table4.html")

####################################################################################
####Table_S5a   MST GPP_new_estimate and NPP_new scaling fits to subsampled trait distributions
######

linear.2 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.3 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(PlotNtoP), data=Peru_Plot_Master.data)
linear.6 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(PlotNtoP), data=Peru_Plot_Master.data)

table5a <- stargazer(linear.2, linear.3, linear.4, linear.5, linear.6, type="html", 
                    title="Metabolic Scaling Theory - Ecosystem NPP and GPP scaling", 
                    single.row=FALSE, ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), ci.separator = " - ", 
                    dep.var.labels.include = FALSE,
                    covariate.labels = c("Aboveground Biomass (kg)", "Temperature 1/kT", "Plot Mean N:P"),
                    model.numbers= FALSE, 
                    column.labels = c("NPP", "GPP", "NPP", "GPP", "NPP", "GPP"),
                    digits = 2)
write.csv(table5a , "Table5a.html")


###################################################################################
####Table_S5b   MST GPP_new_estimate and NPP_new scaling fits to subsampled trait distributions
######  but with PNUE
linear.1 <- lm(log(NPP_new) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.2 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.3 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.6 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.7 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + + MAinvBT + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
linear.8 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)

table5b <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, linear.8, type="html", 
                    title="Metabolic Scaling Theory - Ecosystem NPP and GPP scaling", 
                    single.row=FALSE, ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), ci.separator = " - ", 
                    dep.var.labels.include = FALSE,
                    covariate.labels = c("Aboveground Biomass (kg)", "Temperature 1/kT", "Plot Mean PNUE"),
                    model.numbers= FALSE, 
                    column.labels = c("NPP", "GPP", "NPP", "GPP", "NPP", "GPP", "NPP", "GPP"),
                    digits = 2)
write.csv(table5b , "TableS5b.html")


