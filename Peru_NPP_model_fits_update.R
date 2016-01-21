

##############################
#  Peru CHAMBASA plot analyses
#   Brian J. Enquist
#  12/20/15
#  
##############################

# summary notes - overall simple pairwise correlation of NPP with
# above ground biomass appears to explain most of the variation
# in NPP. Traits expalin little. Climate including temp.
# and solar radiation appear to also explain a lot of the 
# variation but not as much as biomass. Biomass exponents include
# MST predicted value. 

#Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged.csv",header=T)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged2.csv",header=T)

#####
install.packages ("ggplot2", dependencies = TRUE)
install.packages("plyr")
install.packages("ggthemes")
install.packages("reshape2")

library(smatr)
library(car)
library(lmSupport)
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape)
library(gtable)
library(AICcmodavg)


library(ggplot2)
library(reshape2)
library(ggthemes)
library(plyr)

names(Peru_Plot_Master.data)
str(Peru_Plot_Master.data)

######### Calculated Variables ###############
##Calculate 1/kT

Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$Mean.annual.air.temperature..degC.+273.15))

##Calculate the leaf carbon efficiency photosynthesis umol/m^2/s divided by RLeaf
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$mean_n_percent))

        #* calculate leaf carbon efficiency per leaf first before calculating the 
        # plot average?

##Calculate site leaf nitrogen use efficiency
Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))

##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$mean_n_percent)/ (Peru_Plot_Master.data$mean_p_percent))

#################
#####  Plots 
################

#NPP v Biomass
# http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations#log-and-sqrt-transformations

library(scales) # to access break formatting functions, x and y axis are transformed and formatted
myplot_NPP <- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass, NPP)) + geom_point(size = 3) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()

# log-log plot without log tick marks
myplot_NPP

# Show log tick marks
myplot_NPP + annotation_logticks() 



########
#GPP v Biomass

myplot_GPP <- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass, GPP)) + geom_point(size = 3) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()

# log-log plot without log tick marks
myplot_GPP

# Show log tick marks
myplot_GPP + annotation_logticks() 



########
#ResidenceTime v Biomass

myplot_Residence_time <- ggplot(data=Peru_Plot_Master.data, aes(x = Aboveground_biomass, y = Residence_time))
summary(myplot_Residence_time)

myplot_Residence_time_nice <- myplot_Residence_time + geom_point(size = 3)
myplot_Residence_time_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

########
#NPP v Solar Radiation
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Solar.radiation..GJ.m.2.yr.1., y = GPP))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 6)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


########
#elevation v leaf temperature

myplot_leaftemp <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_air_temp)) + geom_point(size = 3) 

  theme_bw()

# log-log plot without log tick marks
  myplot_leaftemp

# Show log tick marks
  myplot_leaftemp + annotation_logticks() 
#** leaf temps decrease with increasing elevation BUT  . . .

######## 
### leaf temperature v photosynthesis
  
myplot_leaftempphoto <- ggplot(Peru_Plot_Master.data, aes(mean_air_temp, mean_photosynthesis)) + geom_point(size = 3) 
  
  theme_bw()
  
  # log-log plot without log tick marks
  myplot_leaftempphoto
  
# No strong relationship between mean air temp of leaf chamber and mean photosyn . . in addtion . . .  
  
  ######## 
  ### site temperature v photosynthesis
  
myplot_plottempphoto <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_photosynthesis)) + geom_point(size = 3) 
  
myplot_plottempphoto <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_photosynthesis)) + geom_point(size = 3) 
  
  
theme_bw()
  
# log-log plot without log tick marks
myplot_plottempphoto
  
# ** there is an interesting positive correlation with site temperature and the plot mean photosynthesis indicating that colder plots have similar if not higher levels of mean photosynthesis. Hypothesis - shift in traits help enable similar rates of metabolism (photosynthesis)
  
  ######## 
  ### site temperature v SLA
  
  myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_sla_lamina_petiole)) + geom_point(size = 3) 
  
  theme_bw()
  
  # log-log plot without log tick marks
  myplot_sitetempsla
  #*weak positive correlation

###### site temperature v SLA mean_sla_lamina
  
myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_sla_lamina)) + geom_point(size = 3) 
  
#myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_sla_lamina)) + geom_point(size = 3) 

theme_bw()
  
# log-log plot without log tick marks
myplot_sitetempsla
  
          #*weak positive correlation
  
###### site temperature v mean_n_percent

myplot_sitetempN <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_n_percent)) + geom_point(size = 3) 

# myplot_sitetempN  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_n_percent)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempN

      #*very weak positive correlation


###### site temperature v mean_p_percent

myplot_sitetempP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_p_percent)) + geom_point(size = 3) 

myplot_sitetempP  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_p_percent)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempP

        #* no strong relationship. looks flat. What about respiration rates? and N:P ratio?

###### site temperature v mean plot N / mean plot P or plot N:P

myplot_sitetempPlotNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoP)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempPlotNtoP 



###### site temperature v RLeaf
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., RLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    #Leaf respiration increases with temperature. Looks like photosynthesis is modified by not leaf respiration.  If correct then the net carbon gain per leaf per carbon respiration changes with temperature. 


###### site temperature v PhotosynthesisPerRLeaf
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    # cool - mean plot Photosynthesis per RLeaf decreases with increasing site temperature. But why should this be?  Change in the nitrogen use efficiency?


###### site temperature v PhotosynthesisPerLeafN
myplot_sitetempLeafNEffic <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerLeafN)) + geom_point(size = 3) 

#myplot_sitetempLeafNEffic  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotosynthesisPerLeafN)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempLeafNEffic 

    ## ** impressive negative correlation - warmer temps have lower mean plot photosynthesis per unit leaf nitrogen. 


###### photosynthesis v PhotosynthesisPerLeafN
myplot_sitetempLeafNEffic <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerLeafN)) + geom_point(size = 3) 

#myplot_sitetempLeafNEffic  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotosynthesisPerLeafN)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempLeafNEffic 

    ## * cool - strong decrease in photosynthesis per leaf N with mean annual air temp. This indicates that colder sites have higher nitrogen use efficiency

###### plot temperature v leafN:P

myplot_sitetempNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoP)) + geom_point(size = 3) 

        #Boltzmann plot
myplot_siteBiomassNtoP  <- ggplot(Peru_Plot_Master.data, aes(x = Aboveground_biomass, y = PlotNtoP)) + geom_point(size = 3) 
theme_bw()
myplot_siteBiomassNtoP

    #* N:P is not related to above ground biomass. Constant?

### Plot mean P vs. plot mean NUE - testing Kerkhoff et al. 2005
myplot_sitePandNUE  <- ggplot(Peru_Plot_Master.data, aes(x = mean_p_percent, y = PhotosynthesisPerLeafN)) + geom_point(size = 3) 
theme_bw()
myplot_sitePandNUE

## exciting  - this seems to support findings from Kerkhoff et al. 2005 showing a positive correlation ebtween foliar P and NUE but here we show this for leaf traits. 

######################################################## 
######################################################## 
## Traits and ecosystem
#NPP v mean plot SLA
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = log10(mean_sla_lamina_petiole), y = NPP))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


#GPP v mean plot SLA
myplot_GPP <- ggplot(data=Peru_Plot_Master.data, aes(x = log10(mean_sla_lamina_petiole), y = GPP))
summary(myplot_GPP)

myplot_GPP_nice <- myplot_GPP + geom_point(size = 3)
myplot_GPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


#mean plot SLA v elevation
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = mean_sla_lamina_petiole))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#varplot SLA v elevation
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = var_sla_lamina_petiole))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#NPP v mean plot SLA lamina
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina, y = NPP))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean plot SLA lamina v elevation
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = mean_sla_lamina))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean plot SLA lamina v mean plot SLA w petiole
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina_petiole, y = mean_sla_lamina))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#NPP v mean plot N
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_n_percent, y = NPP))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean plot leaf N  v elevation
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = mean_n_percent))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean RLeaf  v MAinvBT
myplot_RLeaf <- ggplot(data=Peru_Plot_Master.data, aes(x = MAinvBT, y = RLeaf))
summary(myplot_NPP)

myplot_RLeaf_nice <- myplot_RLeaf + geom_point(size = 3)
myplot_RLeaf_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


### dual plots

myplot_NPP_nice + geom_point(aes(color = Elevation..m.), size =6) # explore effect of elevation
myplot_NPP_nice + geom_point(aes(color = Solar.radiation..GJ.m.2.yr.1.), size =6) # explore effect of solar radiation
myplot_NPP_nice + geom_point(aes(color = Soil.type), size =6) # explore effect of solar radiation
myplot_NPP_nice + geom_point(aes(color = Aboveground_biomass), size =6) # explore effect of above ground biomass
myplot_NPP_nice + geom_point(aes(color = Slope..deg.), size =6) # explore effect of above ground biomass

myplot_NPP_nice + geom_point(aes(color = Ptotal..mg.kg.1.), size =6) # explore effect of above ground biomass

# GPP
myplot_NPP_nice + geom_point(aes(color = mean_sla_lamina_petiole), size =6) # explore effect of mean SLA
myplot_GPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina_petiole, y = GPP))
summary(myplot_GPP)
myplot_GPP_nice <- myplot_GPP + geom_point(size = 4)

myplot_GPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


###################################################
#####  Fitting simple linear models 
##################################################


# NPP and aboveground biomass
m1NPP <- lm(log10(NPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data) 
summary(m1NPP) 
AIC(m1NPP)
modelEffectSizes(m1NPP)  #lm.sunSquares is depreciated
confint(m1NPP)

# GPP and aboveground biomass
m1GPP <- lm(log10(GPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data) 
summary(m1GPP) 
AIC(m1GPP)
modelEffectSizes(m1GPP)  #lm.sunSquares is depreciated
confint(m1GPP)

# NPP and solar radiation
m1NPP <- lm(log10(GPP)~ log10(Solar.radiation..GJ.m.2.yr.1.), data=Peru_Plot_Master.data) 
summary(m1NPP) 
AIC(m1NPP)
modelEffectSizes(m1NPP)  #lm.sunSquares is depreciated
confint(m1NPP)

# NPP and mean plot SLA
m1NPP_SLA <- lm(log10(NPP)~ log10(var_sla_lamina_petiole), data=Peru_Plot_Master.data)
summary(m1NPP_SLA) 
AIC(m1NPP_SLA)
modelEffectSizes(m1NPP_SLA)  #lm.sunSquares is depreciated
confint(m1NPP_SLA)

# NPP and mean plot N
m1NPP_N <- lm(log10(NPP)~ log10(mean_n_percent), data=Peru_Plot_Master.data)
summary(m1NPP_N) 
AIC(m1NPP_N)
modelEffectSizes(m1NPP_N)  #lm.sunSquares is depreciated
confint(m1NPP_N)

## NPP climate via Boltzman temp
m1NPP_temp <- lm(log10(NPP)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1NPP_temp) 
AIC(m1NPP_temp)
modelEffectSizes(m1NPP_temp)  #lm.sunSquares is depreciated
confint(m1NPP_temp)

## NPP climate via preciptiation
m1NPP_precip <- lm(log10(NPP)~ Precipitation..mm.yr.1., data=Peru_Plot_Master.data)
summary(m1NPP_precip) 
AIC(m1NPP_precip)
modelEffectSizes(m1NPP_precip)  #lm.sunSquares is depreciated
confint(m1NPP_precip)

### GPP
m1GPP <- lm(log10(GPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m1GPP) 
AIC(m1GPP)
modelEffectSizes(m1GPP)  #lm.sunSquares is depreciated
confint(m1GPP)

## temperature effect - Boltzmann temperature. Note, NPP is ln(NPP)
m1Temp <- lm(log(NPP)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1Temp) 
AIC(m1Temp)
modelEffectSizes(m1Temp)  #lm.sunSquares is depreciated
confint(m1Temp )


m2 <- lm(log10(GPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m2) 
AIC(m2)
modelEffectSizes(m2)  #lm.sunSquares is depreciated
avPlots(m2)
crPlots(m2)
confint(m2)
vif(m2) 

######################################### 
## Multiple regression 
######################################### 

#NPP and soil moisture and temperature
m3 <- lm(log10(NPP)~ Soil.moisture....+ MAinvBT + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m3) 
AIC(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

## multiple regression with log10 biomass and Boltzmann temperature. 
m4 <- lm(log10(GPP)~ MAinvBT + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 
  ## the fitted biomass scaling exponent is 0.59. Temperature is not important. Could argue that this is the general model to fit so as to extract out the allometric exponent

## multiple regression with log10 biomass and Boltzmann temperature. 
m4 <- lm(log10(GPP)~ MAinvBT + PhotosynthesisPerLeafN + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 

## multiple regression with kerkhoff et al. 2005 model log10 biomass and Boltzmann temperature and plot N:P. 
m4 <- lm(log10(GPP)~ PlotNtoP + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 

## multiple regression with log10 biomass and Solar Radiation. 
m4 <- lm(log10(GPP)~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

## multiple regression with log10 height and Boltzmann temperature. 
m5 <- lm(log10(GPP)~ MAinvBT + log10(Vegetation.height..m.), data=Peru_Plot_Master.data)
summary(m5) 
AIC(m5)
modelEffectSizes(m5)  #lm.sunSquares is depreciated
avPlots(m5)
crPlots(m5)
confint(m5)
vif(m5) 

## multiple regression with log10 biomass and mean plot SLA. 
m6 <- lm(log10(GPP)~ mean_sla_lamina_petiole + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m6) 
AIC(m6)
modelEffectSizes(m6)  #lm.sunSquares is depreciated
avPlots(m6)
crPlots(m6)
confint(m6)
vif(m6) 

m7 <- lm(log10(GPP)~ 1 + mean_sla_lamina_petiole + MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m6) 


m7 <- lm(log10(GPP)~ MAinvBT + log10(Aboveground_biomass) + log10(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
m7 <- lm(log10(GPP)~ log10(Aboveground_biomass) + log10(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)

summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

## multiple regression with log10 biomass and mean plot log10SLA. 
m7 <- lm(log10(NPP)~ log10(mean_sla_lamina_petiole) + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

## linear model on MST prediction with MAinvBT and PlotNtoP
m7 <- lm(log(PlotNtoP) ~ MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 



  #* if PlotNtoP covaries with temperature to compensate for kinetic effects of temp then we would expect that plot N:P scales as ln(N:P)~ 1/kT^0.6 or 0.33. This appears to be the case but confidence intervals are wide. Cold plots have low N:P (more P relative to N)

## linear model on MST prediction with MAinvBT and PlotNtoP
m7 <- lm(log(PhotosynthesisPerLeafN) ~ MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

    ### a similar result for N use efficiency of photosynthesis . .. But NUE increaes positively with MAinvBT. Cold plots have higher NUE

## is N:P ~ NUE ? kerkhoff et al. 2005 argues that N:P is NUE
m7 <- lm(PhotosynthesisPerLeafN ~ PlotNtoP, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

    #* ah, plot N:P is not related to NUE

######################################### 
##### pairs plots
######################################### 

pairs(~ Elevation..m. + Solar.radiation..GJ.m.2.yr.1.+ Precipitation..mm.yr.1.+ + Vegetation.height..m.+ MAinvBT+ log10(Aboveground_biomass) + log10(NPP) + log10(GPP), data=Peru_Plot_Master.data)

pairs(~ Elevation..m. + Solar.radiation..GJ.m.2.yr.1.+ Precipitation..mm.yr.1.+ mean_sla_lamina + mean_sla_lamina_petiole+ mean_n_percent+ var_sla_lamina_petiole + MAinvBT+ log10(Aboveground_biomass) + log10(NPP) + log10(GPP), data=Peru_Plot_Master.data)


options(digits=3)
cor(Peru_Plot_Master.data)


#############################
#############################
#############################
##  automated model fitting
#############################
library(lattice)
library(Hmisc)
library(MASS)
#fit <- lm(GPP ~ Solar.radiation..GJ.m.2.yr.1. + MAinvBT + log10(Aboveground_biomass) + mean_sla_lamina_petiole + var_sla_lamina_petiole, data=Peru_Plot_Master.data)

fit <- lm(NPP ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + mean_sla_lamina_petiole + var_sla_lamina_petiole, data=Peru_Plot_Master.data)
step <- stepAIC(fit, direction="both")
step$anova # display results

## glmulti finds what are the n best models (the confidence set of models) among all possible models (the candidate set, as specified by the user). Models are fitted with the specified fitting function (default is glm) and are ranked with the specified Information Criterion (default is aicc). The best models are found either through exhaustive screening of the candidates, or using a genetic algorithm, which allows very large candidate sets to be adressed. The output can be used for model selection, variable selection, and multimodel inference.

install.packages("glmulti")
#http://core.ac.uk/download/files/153/6303128.pdf
library(glmulti)
library(leaps)
library(MASS)
library(lme4)

## note, method is aicc for small sample corrected 
fit2 <- glmulti(CUE ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + Vegetation.height..m., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

### all variables - takes a while as there are 16900 models
#fit2 <- glmulti(log10(NPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + var_n_percent+ Vegetation.height..m. + Slope..deg.+ Mean.annual.air.temperature..degC. + Soil.total.N.... + Soil.total.C.... + Ptotal..mg.kg.1. + AG.vegetation.biomass + Soil.moisture...., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit2)
#That showed us the best model, now lets look at some of the others
weightable(fit2)
# while log10 biomass is coming out by far the best the model with mean_sla_lamina_petiole + MAinvBT is close at rank 5
tmp <- weightable(fit2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 10,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit2@objects[[1]])
plot(fit2)

#Variable Importance
plot(fit2, type="s")


# now let's do a random mixed effect model 
# define a function that takes a model formula and a dataset as input and then fits a random/mixed-effects meta-regression model to the given data using maximum likelihood estimation

#rma.glmulti <- function(formula, data, ...) {
#  rma(as.formula(paste(deparse(formula))), vi, data=data, method="ML", ...)
#}

### glmulti wrapper functions for mixed effects modelling
## lmer.glmulti from ?glmulti examples
#lmer.glmulti <- function (formula, data, random = "", ...) {
#  lmer(paste(deparse(formula), random), data = data, REML=FALSE, ...)
#} 

#fit2.rma <- glmulti(log10(NPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + Vegetation.height..m., data = Peru_Plot_Master.data, level=1, fitfunc=lmer.glmulti, crit="aicc", confsetsize=128)

fit3 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + log10(mean_n_percent) + log10(Vegetation.height..m.), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit3)
tmp <- weightable(fit3)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit3@objects[[1]])
plot(fit3)

#Variable Importance
plot(fit3, type="s")

#### Now fit many more leaf traits including mean plot photosynthesis
fit4 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + mean_n_percent + log10(Vegetation.height..m.) + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafN + PlotNtoP, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit4)
tmp <- weightable(fit4)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit4@objects[[1]])
plot(fit4)

#Variable Importance
plot(fit4, type="s")


#### Predicting total biomass
fit5 <- glmulti(log10(Aboveground_biomass) ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1. + Elevation..m. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + mean_n_percent + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafN + PlotNtoP, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit5)
tmp <- weightable(fit5)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit5@objects[[1]])
plot(fit5)

#Variable Importance
plot(fit5, type="s")

    #* interesting temperature and elevation comes out on top but the variable model average importance of terms is less than 0.8. Does this suggest that predicting total biomass is difficult? Importance of site history (land slides etc.)?- just the chance that there is a big tree in the sampled plot?




#### Predicting environmental temperature from plot traits
fit6 <- glmulti(MAinvBT ~ Aboveground_biomass + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + mean_n_percent + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafN + PlotNtoP + PhotosynthesisPerLeafN, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit6)
tmp <- weightable(fit6)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit6@objects[[1]])
plot(fit6)

#Variable Importance
plot(fit6, type="s")

    #* results suggest that Plot leafN:P ratio is the best predictor of site temperature. 

#evaluate variable importance, etc.
install.packages("MuMIn")
library(MuMIn)


##################
# Regression Tree
#################
install.packages("tree")
library(tree)
NPP.tree <- tree(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + Vegetation.height..m., data = Peru_Plot_Master.data, mindev=0)

NPP.tree
plot(residuals(NPP.tree) ~ predict(NPP.tree))
plot(NPP.tree, type = "uniform")
text(NPP.tree, cex = 0.5, all = T)


### Pruning tree
plot(prune.tree(NPP.tree))
NPP.tree.prune <- prune.tree(NPP.tree, best = 4)
plot(NPP.tree.prune, type = "uniform")
text(NPP.tree.prune, cex = 0.5, all = T)

