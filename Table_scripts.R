###############
##### summary stats
#####   B. J. Enquist
#####
install.packages("stargazer")
install.packages("dcolumn")
library(stargazer)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged4.csv",header=T)

######### Calculated Variables ###############
##Calculate Boltzmann 1/kT
Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$Mean.annual.air.temperature..degC.+273.15))

##Calculate the leaf N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$mean_n_percent))

Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))

##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPPperNMeanMean <- ((Peru_Plot_Master.data$NPP)/(Peru_Plot_Master.data$NMeanMean))

Peru_Plot_Master.data$GPPperNMeanMean <- ((Peru_Plot_Master.data$GPP)/(Peru_Plot_Master.data$NMeanMean))


#* calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))


##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$mean_n_percent)/ (Peru_Plot_Master.data$mean_p_percent))
Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

#Calculate NPPLeaf/RLeaf - Production per carbon respired
Peru_Plot_Master.data$NPPLeafperRLeaf <- ((Peru_Plot_Master.data$NPPLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST prediction
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP1 <- ((Peru_Plot_Master.data$NPP)/(Peru_Plot_Master.data$MST_AGB))

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
####Table_1 
######
stargazer(Peru_Plot_Master.data$elevation_m)
linear.1 <- lm(mean_n_percent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.2 <- lm(mean_p_percent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.3 <- lm(PlotNtoP ~ Elevation..m., data=Peru_Plot_Master.data)
linear.4 <- lm(mean_c_percent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.5 <- lm(mean_sla_lamina_petiole ~ Elevation..m., data=Peru_Plot_Master.data)
linear.6 <- lm(mean_photosynthesis ~ Elevation..m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotosynthesisPerLeafN ~ Elevation..m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table1 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Most Abunant Species)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "N:P", "%C", "SLA", "Photo", "N efficiency" ), digits = 2)
#no.space=TRUE
#write.csv(table1 , "Table1.html")

####################################################################################
####Table_2  MST Boltzmann fits to traits of sampled individuals of abundant species
######

stargazer(Peru_Plot_Master.data$elevation_m)
linear.1 <- lm(log(PlotNtoP) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.2 <- lm(log(PhotosynthesisPerLeafN) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.3 <- lm(log(mean_photosynthesis) ~ MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log(TransformedNpercent) ~ MAinvBT, data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table2 <- stargazer(linear.1, linear.2, linear.3, linear.4, type="html", title="Metabolic Scaling Theory - Boltzmann Temperature Results", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("ln(N:P)", "ln(N effciency)", "ln(photosynthesis)", "ln(%N)"), digits = 2)

####################################################################################
####Table_3 
######

stargazer(Peru_Plot_Master.data$elevation_m)
linear.1 <- lm(TransformedNpercent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.2 <- lm(TransformedPpercent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.3 <- lm(PlotNtoPMean ~ Elevation..m., data=Peru_Plot_Master.data)
linear.4 <- lm(TransformedCpercent ~ Elevation..m., data=Peru_Plot_Master.data)
linear.5 <- lm(Transformedmean_sla_lamina_petiole ~ Elevation..m., data=Peru_Plot_Master.data)
linear.6 <- lm(PhotoMeanMean ~ Elevation..m., data=Peru_Plot_Master.data)
linear.7 <- lm(PhotoPerLeafNMean ~ Elevation..m., data=Peru_Plot_Master.data)

####### Next, use R knitr to generate a pretty table. 
####### Next click File -> New File -> R HTML 
####### Next click Knit HTML and choose where to save the file. 
####### You'll get a HTML file,
####### Remove everything except the <htlm> tags. 
#######vCopy the output from stargazer function:
table3 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5, linear.6, linear.7, type="html", title="Change in Plot Traits with Elevation (Trait Distribution Subsampling)", single.row=FALSE, ci=TRUE,ci.level=0.95, omit.stat=c("f", "ser"), ci.separator = " - ", dep.var.labels.include = FALSE, model.numbers= FALSE, column.labels = c("%N", "%P", "N:P", "%C", "SLA", "Photo", "N efficiency" ), digits = 2)
#no.space=TRUE
#write.csv(table1 , "Table1.html")


####################################################################################
####Table_4   MST Boltzmann fits to subsampled trait distributions
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

####################################################################################
####Table_5   MST GPP and NPP scaling fits to subsampled trait distributions
######

linear.1 <- lm(log10(NPP) ~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.2 <- lm(log10(GPP_estimated) ~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
linear.3 <- lm(log10(NPP) ~ log10(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.4 <- lm(log10(GPP_estimated) ~ log10(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
linear.5 <- lm(log10(NPP) ~ log10(Aboveground_biomass2) + PlotNtoP, data=Peru_Plot_Master.data)
linear.6 <- lm(log10(GPP_estimated) ~ log10(Aboveground_biomass2) + PlotNtoP, data=Peru_Plot_Master.data)


table5 <- stargazer(linear.1, linear.2, linear.3, linear.4, linear.5,linear.6, type="html", 
                    title="Metabolic Scaling Theory - Ecosystem NPP and GPP scaling", 
                    single.row=FALSE, ci=TRUE,ci.level=0.95, 
                    omit.stat=c("f", "ser"), ci.separator = " - ", 
                    dep.var.labels.include = FALSE,
                    covariate.labels = c("Aboveground Biomass (kg)", "Temperature 1/kT", "Plot Mean N:P"),
                    model.numbers= FALSE, 
                    column.labels = c("NPP", "GPP", "NPP", "GPP", "NPP", "GPP"), digits = 2)

