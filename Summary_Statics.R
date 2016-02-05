###############
summary stats

install.packages("stargazer")
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



stargazer(Peru_Plot_Master.data$elevation_m)
linear.1 <- lm(TransformedNpercent ~mean_n_percent, data=Peru_Plot_Master.data)
linear.2 <- lm(TransformedPpercent ~ mean_p_percent, data=Peru_Plot_Master.data)

table <- stargazer(linear.1, linear.2, type="html",
          title="Regression Results", single.row=TRUE,
          ci=TRUE, ci.level=0.9, omit.stat=c("f", "ser"))
