
#==========================================================================
# Code set 1
#
# Assessing trait-based scaling theory in tropical forests spanning a broad temperature gradient 
# 
#Authors:
#Brian J. Enquist1,2*, Lisa Patrick Bentley3, Alexander Shenkin3, Brian Maitner1, Van Savage4, Sean Michaletz1,5, Benjamin Blonder3, Vanessa Buzzard1, Tatiana Erika Boza Espinoza6, William Farfan-Rios 7,8 Chris Doughty9, Gregory R. Goldsmith10, Roberta E. Martin11, Norma Salinas3,8, Miles Silman8, Sandra Díaz12, Gregory P. Asner11 & Yadvinder Malhi3
#
# Main analyses and graphs
# 8/2/17
#==========================================================================

# load chambasa data
Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged7.csv",header=T)


#==========================================================================
#### load libraries

install.packages ("ggplot2", dependencies = TRUE)
install.packages("plyr")
install.packages("ggthemes")
install.packages("reshape2")
install.packages("glmulti")
install.packages("grid")
install.packages("gtable")

library(car)
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape)
library(gtable)
library(AICcmodavg)
library(gtable)
library(gridExtra)
library(glmulti)
library(corrplot)
library(grid)

library(reshape2)
library(ggthemes)
library(plyr)

names(Peru_Plot_Master.data)
str(Peru_Plot_Master.data)

#==========================================================================
#Calculated Variables 
#==========================================================================

#Calculate Boltzmann 1/kT
Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.+273.15))

Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$AboveGroundBiomass_Malhi_2017)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

##Calculate the leaf PNUE or N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean*100))

# calculation amax per unit mass
Peru_Plot_Master.data$amax.sun.massSpecific <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)*(1/(Peru_Plot_Master.data$lma.sun.mu.abundance)))

#calculate PNUE on a mass basis 
Peru_Plot_Master.data$Photosynthesis_massPerLeafN <- ((Peru_Plot_Master.data$amax.sun.massSpecific)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPP_newperNMeanMean <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

Peru_Plot_Master.data$GPPperNMeanMean <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

#Calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
#Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))

##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$n_percent.sun.mu.abundance)/ (Peru_Plot_Master.data$p_corrected_percent.sun.mu.abundance))

Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

#Calculate NPP_newLeaf/RLeaf - Production per carbon respired
#Peru_Plot_Master.data$NPP_LeafperRLeaf <- ((Peru_Plot_Master.data$NPP_newLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST predictions
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$AboveGroundBiomass_Malhi_2017)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_GPP1Upper <- ((Peru_Plot_Master.data$GPP_Malhi_2017_UpperRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_GPP1Lower <- ((Peru_Plot_Master.data$GPP_Malhi_2017_LowerRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP1Upper <- ((Peru_Plot_Master.data$NPP_Malhi_2017_UpperRange)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP1Lower <- ((Peru_Plot_Master.data$NPP_Malhi_2017_LowerRange)/(Peru_Plot_Master.data$MST_AGB))



#==========================================================================
#### CHAMBASA - Plotting Functions
### Panel Plot - define multiplot function first
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#==========================================================================
# Community trait PCA analyses
#==========================================================================

# Traits from most abundant individuals 
Photo_a <- Peru_Plot_Master.data$mean_photosynthesis
PN_a <- 1/Peru_Plot_Master.data$PlotNtoP
P_a <- Peru_Plot_Master.data$mean_p_percent
N_a <- Peru_Plot_Master.data$mean_n_percent
C_a <- Peru_Plot_Master.data$mean_c_percent
LMA_a <- 1/Peru_Plot_Master.data$mean_sla_lamina_petiole
#SLA <- Peru_Plot_Master.data$mean_sla_lamina_petiole
PNUE_a <- Peru_Plot_Master.data$PhotosynthesisPerLeafN

## Subsampling total plot distribution method
Photo <- Peru_Plot_Master.data$PhotoMeanMean
PN <- 1/Peru_Plot_Master.data$PlotNtoPMean
P <- Peru_Plot_Master.data$PMeanMean
N <- Peru_Plot_Master.data$NMeanMean
C <- Peru_Plot_Master.data$CMeanMean
LMA <- 1/Peru_Plot_Master.data$SLAMeanMean
#SLA <- Peru_Plot_Master.data$SLAMeanMean
PNUE <- Peru_Plot_Master.data$PhotoPerLeafNMean

chambasapca_a <- princomp(~ LMA_a + N_a + C_a + P_a + Photo_a + PN_a + PNUE_a, cor=TRUE)
summary(chambasapca_a)
loadings(chambasapca_a)
loadings_a <- loadings(chambasapca_a)
write.csv(loadings_a, "PCA_a_CommunityTrait.csv")

chambasapca_c <- princomp(~ LMA + N + C + P + Photo + PN + PNUE, cor=TRUE)
summary(chambasapca_c)
loadings(chambasapca_c)
loadings_c <- loadings(chambasapca_c)
write.csv(loadings_c, "PCA_c_CommunityTrait.csv")


png("Figure_PCA_a_gradient_plot.png", units="in", width=5, height=4, pointsize=9, res=900)
PCA_Plot_a <- biplot(chambasapca_a, col=c("gray","black"),cex=c(0.8,0.5))
PCA_Plot_a
dev.off()

png("Figure_PCA_c_gradient_plot.png", units="in", width=5, height=4, pointsize=9, res=900)
PCA_Plot_c <- biplot(chambasapca, col=c("gray","black"),cex=c(0.8,0.5))
PCA_Plot_c
dev.off()

screeplot(chambasapca)
chambasapca$scores
chambasapca$loadings
#* first two principle components explain about 72% of the variation. 

first.component.scores <- chambasapca$scores[,1]
summary(first.component.scores)
length(first.component.scores)
#write.csv(first.component.scores, "first.component.scores.abund.csv")
#write.csv(first.component.scores, "first.component.scores.community.csv")

second.component.scores <- chambasapca$scores[,2]
summary(second.component.scores)
length(second.component.scores)
#write.csv(second.component.scores, "second.component.scores.abund.csv")
#write.csv(first.component.scores, "second.component.scores.community.csv")

nrow(Peru_Plot_Master.data)
length(Peru_Plot_Master.data$Elevation.m.)
elev <- (Peru_Plot_Master.data$Elevation.m.)
length(elev)


########### Plots of PCA abundant axes vs. environmental variation
PCA1a <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PCA1ScoresPlotTraits_a)) +
  geom_point(size = 4, color="darkgrey") +
  xlab(expression(paste("Mean Annual Temperature", " ", C^{o}))) + ylab("PCA1 Abundant Species") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth(method = "lm", se=FALSE, color="black")+
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
PCA1a

ModelPCA1_temp <- lm(PCA1ScoresPlotTraits_a ~ MeanAnnualAirTemperature.degC., data=Peru_Plot_Master.data)
summary(ModelPCA1_temp)  # y ~ x

ModelPCA1_solar <- lm(PCA1ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1., data=Peru_Plot_Master.data)
summary(ModelPCA1_solar)


PCA2a <- ggplot(Peru_Plot_Master.data, aes(Precipitation.mm.yr.1., PCA2ScoresPlotTraits_a)) +
  xlab(expression(paste("Mean Annual Precipitation", " ", "mm", " ", yr^{-1}))) + ylab("PCA2 Abundant Species") +
  geom_point(size = 4, color="darkgrey") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth(method = "lm", se=FALSE, color="black")+
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
PCA2a

Model_PCA2a <- lm(PCA2ScoresPlotTraits_a ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
summary(Model_PCA2)

png("Figure_PCA_a_multi.png", units = "px", width=1500, height=900, res=200)

multiplot(PCA1a, PCA2a, cols=2)
dev.off()

#########################################################
########### Plots of PCA community axes vs. environmental variation
PCA1c <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PCA1ScoresTraitSample)) +
  xlab(expression(paste("Mean Annual Temperature", " ", C^{o}))) + ylab("PCA1 Community") +
  geom_point(size = 4, color="darkgrey") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth(method = "lm", se=FALSE, color="black")+
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
PCA1c

ModelPCA1 <- lm(PCA1ScoresPlotTraits_a ~ MeanAnnualAirTemperature.degC., data=Peru_Plot_Master.data)
summary(ModelPCA1)  # y ~ x


PCA2c <- ggplot(Peru_Plot_Master.data, aes(Soil.moisture...., PCA2ScoresTraitSample)) +
  xlab(expression(paste("Soil Moisture", " ", "%"))) + ylab("PCA2 Community") +
  geom_point(size = 4, color="darkgrey") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth(method = "lm", se=FALSE, color="black")+
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
PCA2c 

PCA2_2 <- ggplot(Peru_Plot_Master.data, aes(Precipitation.mm.yr.1., PCA2ScoresTraitSample)) +
  xlab(expression(paste("Precipitation", " ", "mm", " ", yr^{-1} ))) + ylab("PCA2 Community") +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
PCA2_2 

Model_PCA2 <- lm(PCA2ScoresPlotTraits_a ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
summary(Model_PCA2)

png("Figure_PCA_c_multi.png", units = "px", width=1500, height=900, res=200)

multiplot(PCA1c, PCA2c, cols=2)
dev.off()



png("Figure_PCA_a_multi.png", units = "px", width=1500, height=900, res=200)

multiplot(PCA_a, PCA2_2, cols=2)
dev.off()

#==========================================================================
# Best predictors of PCA1 and PCA2 variation, climate?
#==========================================================================

PCAfita <- glmulti(PCA1ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + Soil.moisture.... + Slope..deg. + Aspect..deg., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfita)
tmpPCA1a <- weightable(PCAfit)
tmpPCA1a <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmpPCA1a
summary(PCAfita@objects[[1]]) 
summary(PCAfita@objects[[2]])
summary(PCAfita@objects[[3]])
# a best model includes both air temp and solar radiation although inclusion of solar radiation is marginally significant. Finding does somewhat support recent Fyllas et al. 2017 paper but tempreature is main driver of trait shifts with secondary signal of solar radiation 
summary(PCAfita@objects[[4]])
plot(PCAfita)
#Variable Importance
plot(PCAfita, type="s")
#PCA1 is mainly temperature
Model_1 <- lm(PCA1ScoresPlotTraits_a ~ MeanAnnualAirTemperature.degC. + SolarRadiation.GJ.m.2.yr.1., data=Peru_Plot_Master.data)
summary(Model_1)

PCA2afit <- glmulti(PCA2ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + Soil.moisture.... + Slope..deg. + Aspect..deg., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCA2afit)
#tmp <- weightable(PCA2fit)
#tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
#tmp
summary(PCA2afit@objects[[1]])
summary(PCA2afit@objects[[2]])
summary(PCA2afit@objects[[3]])
plot(PCA2afit)
#Variable Importance
plot(PCA2afit, type="s")


########  Whole Community PCA scores

PCAfitSample <- glmulti(PCA1ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + Soil.moisture.... + Slope..deg. + Aspect..deg., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfitSample)
tmp <- weightable(PCAfitSample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfitSample@objects[[1]])
summary(PCAfitSample@objects[[2]])
plot(PCAfitSample)
#Variable Importance
plot(PCAfitSample, type="s")
#PCA1 is mainly mean annual tempreature

PCAfit2Sample <- glmulti(PCA2ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + Soil.moisture.... + Slope..deg. + Aspect..deg., data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit2Sample)
tmp <- weightable(PCAfit2Sample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit2Sample@objects[[1]])
summary(PCAfit2Sample@objects[[2]])
plot(PCAfit2Sample)
#Variable Importance
plot(PCAfit2Sample, type="s")

## pairwise correlations with abiotic environment and PCA scores
#install.packages("leaps")
#library(leaps)
#leaps=regsubsets(PCA1ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, nbest=10)
#plot(leaps, scale="adjr2")
#plot(leaps, scale="bic")
#summary(leaps)
#library(car)
#subsets(leaps, statistic="bic")
#subsets(leaps, statistic="adjr2")


#require(lattice)
#require(ggplot2)
#PCA1 <- Peru_Plot_Master.data$PCA1ScoresPlotTraits_a
#SolarRad. <- Peru_Plot_Master.data$SolarRadiation.GJ.m.2.yr.1.
#Precip. <- Peru_Plot_Master.data$Precipitation.mm.yr.1.
#Elev. <- Peru_Plot_Master.data$Elevation.m.
#Temp. <- Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.
#pairs(~ PCA1 + SolarRad. + Precip. + Temp. + Elev., pch = 21)

## correlation matrix of abiotic correlations with PCA1 showing that temperature is best pairwise predictor


#==========================================================================
# Corrplots of PCA abiotic correlations abundant taxa
#==========================================================================

library(corrplot)
#ftp://cran.r-project.org/pub/R/web/packages/corrplot/vignettes/corrplot-intro.html

## add p-values on no significant coefficient
# corrplot(M, p.mat = res1[[1]], insig = "p-value")

lpp_axis_a <- data.frame("PCA1 Dominant" = Peru_Plot_Master.data$PCA1ScoresPlotTraits_a,
                                 "PCA2 Dominant" = Peru_Plot_Master.data$PCA2ScoresPlotTraits_a,
                                 "Elevation" = Peru_Plot_Master.data$Elevation.m.,
                                 "Solar Radiation" = Peru_Plot_Master.data$SolarRadiation.GJ.m.2.yr.1., 
                                 "MAT"= Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.,
                                 "Precipitation" = Peru_Plot_Master.data$Precipitation.mm.yr.1.,
                                 "Soil Moisture" = Peru_Plot_Master.data$Soil.moisture....)
#my_data_a <- Peru_Plot_Master.data[, c(39,40,6,18,19,20,21)]
res_a <- cor(na.omit(lpp_axis_a))
#corrplot(res_a, method="circle")
#corrplot(res_a, method="ellipse")
#corrplot(res_a, method="number", type="upper", insig = "blank")
#corrplot(res_a, insig = "blank")
#corrplot.mixed(res_a, lower="ellipse", upper="number")

png("CorrTable_PCA_a_climate.png", units = "px", width=1000, height=1500, res=150)
corrplot(res_a, method="number", type="upper", insig = "blank")
dev.off()

#### Corrplots of PCA abiotic correlations whole-community  ##############

lpp_axis <- data.frame("PCA1 Community" = Peru_Plot_Master.data$PCA1ScoresTraitSample,
                                 "PCA2 Community" = Peru_Plot_Master.data$PCA2ScoresTraitSample,
                                 "PCA1 Dominant" = Peru_Plot_Master.data$PCA1ScoresPlotTraits_a,
                                 "PCA2 Dominant" = Peru_Plot_Master.data$PCA2ScoresPlotTraits_a,
                                 "Elevation" = Peru_Plot_Master.data$Elevation.m.,
                                 "Solar Radiation" = Peru_Plot_Master.data$SolarRadiation.GJ.m.2.yr.1., 
                                 "MAT"= Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.,
                                 "Precipitation" = Peru_Plot_Master.data$Precipitation.mm.yr.1.,
                                 "Soil Moisture" = Peru_Plot_Master.data$Soil.moisture....)

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(lpp_axis,0.95)
res2 <- cor.mtest(lpp_axis,0.99)
## specialized the insignificant value according to the significant level

#my_data <- Peru_Plot_Master.data[, c(41,42,6,18,19,20,21)]
#my_data <- Peru_Plot_Master.data[, c(PCA1_Community,PCA2_Community,Elevation,Solar_Radiation,MAT,Precipitation,Soil_Moisture)]
res_Community <- cor(na.omit(lpp_axis))
#corrplot(res_Community, method="circle")
#corrplot(res_Community, method="ellipse")
#corrplot(res_Community, p.mat = res1[[1]], insig = "p-value")
#corrplot(res_Community, method="number", type="upper", insig = "blank", p.mat = res_Community[[1]], insig = "p-value")
#corrplot(res_Community, method="number", type="upper", insig = "p-value")
#corrplot(res, insig = "blank")
#corrplot.mixed(res, lower="ellipse", upper="number")

png("CorrTable_PCA_Community_climate.png", units = "px", width=1000, height=1500, res=150)
corrplot(res_Community, method="number", type="upper", p.mat = res1[[1]], sig.level=0.05)
dev.off()


#==========================================================================
###  Trait distributions across elevation
###   Relationship w dominant (abundant) community species traits and elevation
#==========================================================================

png("Figure_Plot_CMW_a_Traits.png", units = "px", width=1800, height=1800, res=200)
#dev.off()
T1 <- myplot_TDTNAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., n_percent.sun.mu.abundance)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) + 
  geom_errorbar(aes(ymin=n_percent.sun.mu.abundance-n_percent.sun.se.abundance, ymax=n_percent.sun.mu.abundance+n_percent.sun.se.abundance), width=.01,
                position=position_dodge(.05)) +
  geom_smooth(method=lm)
T1

Model_T1 <- lm(n_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T1)


T2 <- myplot_TDTPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., p_corrected_percent.sun.mu.abundance)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf P") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=p_corrected_percent.sun.mu.abundance-p_corrected_percent.sun.se.abundance, ymax=p_corrected_percent.sun.mu.abundance+p_corrected_percent.sun.se.abundance), width=.01,
                position=position_dodge(.05))
T2

Model_T2 <- lm(p_corrected_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T2)

P_N <- (1/Peru_Plot_Master.data$PlotNtoP)

T3 <- myplot_TDTNtoPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., P_N)) +
  xlab("Elevation (m)") + ylab("Plot Leaf P:N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTNtoPAve
T3

#PlotNtoP ,P_N
Model_T3 <- lm(PlotNtoP ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T3)


T4 <- myplot_TDTCAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., c_percent.sun.mu.abundance)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf C ") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=c_percent.sun.mu.abundance-c_percent.sun.se.abundance, ymax=c_percent.sun.mu.abundance+c_percent.sun.se.abundance), width=.01,
                position=position_dodge(.05))
T4

Model_T4 <- lm(c_percent.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T4)


T5 <- myplot_TDTSLAAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., lma.sun.mu.abundance)) +
  xlab("Elevation (m)") + 
  ylab(bquote('Plot LMA ('~ g~m^-2*')'))+
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=lma.sun.mu.abundance-lma.sun.se.abundance, ymax=lma.sun.mu.abundance+lma.sun.se.abundance), width=.01,
                position=position_dodge(.05)) +
  geom_smooth(method=lm)
#myplot_TDTSLAAve
T5

Model_T5 <- lm(lma.sun.se.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T5)


T6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., amax.sun.mu.abundance)) +
  
  #T6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., Peru_Plot_Master.data$amax.sun.massSpecific)) +
  #xlab("Elevation (m)") + 
  #ylab(bquote('Plot Amass ('*mu~ 'mol' ~CO[2]~ g^-1~s^-1*')')) + 
  
  xlab("Elevation (m)") + 
  ylab(bquote('Leaf Photosynthesis ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=amax.sun.mu.abundance-amax.sun.se.abundance, ymax=amax.sun.mu.abundance+amax.sun.se.abundance), width=.01,
                position=position_dodge(.05))
geom_smooth(method=lm)
#myplot_TDTPhotoAve
T6

Model_T6 <- lm(amax.sun.mu.abundance ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T6)


T7 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotosynthesisPerLeafN)) +
  xlab("Elevation (m)") + 
  ylab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ m^-2~gN-1~s^-1*')')) + 
  
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)

T7

Model_T7 <- lm( PhotosynthesisPerLeafN ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T7)

#Model_T7 <- lm( (log(PhotosynthesisPerLeafN)) ~ MAinvBT, data=Peru_Plot_Master.data)
#summary(Model_T7)

#

T8 <- myplot_PhotovPNUE <- ggplot(Peru_Plot_Master.data, aes(PhotosynthesisPerLeafN, amax.sun.mu.abundance)) +
  xlab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ g^-1~s^-1*')')) + 
  ylab(bquote('Plot Leaf Photo. ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotovPNUE
#wow - highly positive . . variation in mean photosynthesis appears to be driven by PNUE
T8

Model_T8 <- lm(amax.sun.mu.abundance ~ PhotosynthesisPerLeafN, data=Peru_Plot_Master.data)
summary(Model_T8)

T9 <- myplot_NtoPvPNUE <- ggplot(Peru_Plot_Master.data, aes(PlotNtoP, PhotosynthesisPerLeafN)) +
  xlab("Plot Leaf N:P") + 
  ylab(bquote('Plot Leaf N efficiency ('*~CO[2]~ g^-1~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NtoPvPNUE
#huh . . . no correlation. Due to covariation?
T9

# use PhotosynthesisPerLeafN or PhotoPerLeafNMean 

T10 <- myplot_NtoPvPNUE <- ggplot(Peru_Plot_Master.data, aes(p_corrected_percent.sun.mu.abundance, PhotosynthesisPerLeafN)) +
  #T10 <- myplot_NtoPvPNUE <- ggplot(Peru_Plot_Master.data, aes(PMeanMean, PhotoPerLeafNMean)) +
  xlab("Plot Leaf P") + 
  ylab(bquote('Plot Leaf N efficiency ('*~CO[2]~ g^-1~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)

T10
#myplot_NtoPvPNUE

T11 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., Peru_Plot_Master.data$Photosynthesis_massPerLeafN)) +
  xlab("Elevation (m)") + 
  ylab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ g^-1~gN-1~s^-1*')')) + 
  
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
#geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
#position=position_dodge(0.05)) +
#geom_smooth(method=lm)  # not significant on a mass basis
#myplot_TDTPhotoAve
T11

Model_T11 <- lm(Photosynthesis_massPerLeafN ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T11)

multiplot(T1,T2,T3,T4,T5,T6,T7,T8, cols=3)
dev.off()

#==========================================================================
# Trait distributions across elevation
# Relationship between community species traits estiamated by subsampling routine and elevation
#==========================================================================


png("Figure_Plot_CMW_c_Traits.png", units = "px", width=1800, height=1800, res=200)
#dev.off()
C1 <- myplot_TDTNAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NMeanMean*100)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) + 
  geom_errorbar(aes(ymin=NMeanLower*100, ymax=NMeanUpper*100), width=.01,
                position=position_dodge(.05)) +
  geom_smooth(method=lm)
C1

Model_C1 <- lm(NMeanMean*100 ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C1)


C2 <- myplot_TDTPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PMeanMean*100)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf P") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=PMeanLower*100, ymax=PMeanUpper*100), width=.01,
                position=position_dodge(.05))
C2

Model_C2 <- lm(PMeanMean*100 ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C2)


C3 <- myplot_TDTNtoPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., 1/PlotNtoPMean)) +
  xlab("Elevation (m)") + ylab("Plot Leaf P:N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTNtoPAve
C3

Model_C3 <- lm(PlotNtoPMean ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C3)


C4 <- myplot_TDTCAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CMeanMean*100)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf C ") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=CMeanLower*100, ymax=CMeanUpper*100), width=.01,
                position=position_dodge(.05))
C4

Model_C4 <- lm(CMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C4)


C5 <- myplot_TDTSLAAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., 1/SLAMeanMean)) +
  xlab("Elevation (m)") + 
  ylab(bquote('Plot LMA ('~ kg~m^-2*')'))+
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=1/SLAMeanLower, ymax=1/SLAMeanUpper), width=.01,
                position=position_dodge(.05)) +
  geom_smooth(method=lm)
#myplot_TDTSLAAve
C5

Model_C5 <- lm(1/SLAMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C5)


C6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoMeanMean)) +
  
  xlab("Elevation (m)") + 
  ylab(bquote('Leaf Photosynthesis ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.01,
                position=position_dodge(.05))
geom_smooth(method=lm)

C6

Model_C6 <- lm(PhotoMeanMean ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C6)


C7 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoPerLeafNMean)) +
  xlab("Elevation (m)") + 
  ylab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ m^-2~gN-1~s^-1*')')) + 
  
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)

C7

Model_C7 <- lm( PhotoPerLeafNMean ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_C7)


C8 <- myplot_PhotovPNUE <- ggplot(Peru_Plot_Master.data, aes(PhotoPerLeafNMean, PhotoMeanMean)) +
  xlab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ g^-1~s^-1*')')) + 
  ylab(bquote('Plot Leaf Photo. ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotovPNUE
#wow - highly positive . . variation in mean photosynthesis appears to be driven by PNUE
C8

Model_C8 <- lm(PhotoPerLeafNMean ~ PhotoMeanMean, data=Peru_Plot_Master.data)
summary(Model_C8)

multiplot(C1,C2,C3,C4,C5,C6,C7,C8, cols=3)
dev.off()



#==========================================================================
# Trait Distributions from subsampled trait distribuitons
#==========================================================================

# Nitrogen 

#png("Figure_Sampled_Nitrogen2.png", units = "px", width=1800, height=1800, res=200)
#dev.off()
#theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
N1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NMeanMean)) +
  xlab("Elevation (m)") + ylab("Mean % N") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()
#myplot_TDTNMean
N1
#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean

ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
N2 <- myplot_TDTNVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NVarianceMean)) +
  xlab("Elevation (m)") + ylab("  Variance % N ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NVarianceLower, ymax=NVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()
#myplot_TDTNVar
N2
ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
N3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NSkewnessMean)) +
  xlab("Elevation (m)") + ylab("  Skewness % N ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NSkewnessLower, ymax=NSkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm) +
  #geom_smooth()
  expand_limits(y=c(-1.0,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTSkew 
N3
ModelTDTNSkew <- lm(Elevation.m. ~ NSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTNSkew)
confint(ModelTDTNSkew)
coef(ModelTDTNSkew)

#Plot Kurtosis
N4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., (NKurtosisMean-3))) +
  xlab("Elevation (m)") + ylab("Kurtosis % N ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=(NKurtosisLower-3), ymax=(NKurtosisUpper-3)), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  
  #geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,10)) +
  geom_hline(yintercept = 0)

N4

ModelTDTNVar <- lm(Elevation.m. ~ NKurtosisMean, Peru_Plot_Master.data)
summary(NKurtosisMean)
confint(NKurtosisMean)
coef(NKurtosisMean)

# Multipanel plot
#multiplot(N1, N2, N3, N4, cols=2)
#dev.off()


##################
#### Leaf Carbon
################

#png("Figure_Sampled_Carbon2.png", units = "px", width=1800, height=1800, res=200)

#theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))

C1 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CMeanMean)) +
  xlab("Elevation (m)") + ylab("  Mean % Leaf C") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()
#myplot_TDTCMean
C1         

#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean

ModelTDTCMean <- lm(Elevation.m. ~ CMeanMean, Peru_Plot_Master.data)
summary(ModelTDTCMean)
confint(ModelTDTCMean)
coef(ModelTDTCMean)

#Plot Variance
C2 <-  ggplot(Peru_Plot_Master.data, aes(Elevation.m., CVarianceMean)) +
  xlab("Elevation (m)") + ylab("  Variance % C ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=CVarianceLower, ymax=CVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0))
#geom_smooth(method=lm)
#geom_smooth()
#myplot_TDTCVar

ModelTDTCVar <- lm(Elevation.m. ~ CVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTCVar)
confint(ModelTDTCVar)
coef(ModelTDTCVar)
C2         


#Plot Skewness

C3 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CSkewnessMean)) +
  xlab("Elevation (m)") + ylab("  Skewness % C ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=CSkewnessLower, ymax=CSkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+
  expand_limits(y=c(-1.0,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTCSkew 
C3        
ModelTDTCSkew <- lm(Elevation.m. ~ CSkewnessMean, Peru_Plot_Master.data)
summary(ModelTDTCSkew)
confint(ModelTDTCSkew)
coef(ModelTDTCSkew)


#Plot Kurtosis
C4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., (CKurtosisMean-3))) +
  xlab("Elevation (m)") + ylab("Kurtosis % C ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=(CKurtosisLower-3), ymax=(CKurtosisUpper-3)), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,10)) +
  geom_hline(yintercept = 0)
#myplot_TDTCKurtosis 
C4         

ModelTDTCKurt <- lm(Elevation.m. ~ CKurtosisMean, Peru_Plot_Master.data)
summary(ModelTDTCKurt)
confint(ModelTDTCKurt)
coef(ModelTDTCKurt)

# Multipanel plot
#multiplot(C1, C2, C3, C4, cols=2)
#dev.off()


######################
#### Leaf Phosphorus
################

#png("Figure_Sampled_Phosphorus2.png", units="in", width=5, height=4, pointsize=12, res=900)
#theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))

P1 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PMeanMean)) +
  xlab("Elevation (m)") + ylab("  Mean % Leaf P") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PMeanLower, ymax=PMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0))
#geom_smooth(method=lm)+
#geom_smooth()
#myplot_TDTPMean

#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean

ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
P2 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PVarianceMean)) +
  xlab("Elevation (m)") + ylab("  Variance % P ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PVarianceLower, ymax=PVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0))
#geom_smooth(method=lm)+
#geom_smooth()
#myplot_TDTPVar

ModelTDTPVar <- lm(Elevation.m. ~ PVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTPVar)
confint(ModelTDTPVar)
coef(ModelTDTPVar)

#Plot Skewness

P3 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PSkewnessMean)) +
  xlab("Elevation (m)") + ylab("Skewness % P ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PSkewnessLower, ymax=PSkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+
  expand_limits(y=c(-1.0,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTPSkew 

ModelTDTNVar <- lm(Elevation.m. ~ PSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
P4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., (PKurtosisMean-3))) +
  xlab("Elevation (m)") + ylab("Kurtosis % P ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=(PKurtosisLower-3), ymax=(PKurtosisUpper-3)), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,10)) +
  geom_hline(yintercept = 0)
#myplot_TDTPKurtosis 
P4

ModelTDTNVar <- lm(Elevation.m. ~ PKurtosisMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
#multiplot(P1, P2, P3, P4, cols=2)
#dev.off()

############################
#### Leaf Photosynthesis
################

#png("Figure_Sampled_Photo.png", units="in", width=5, height=4, pointsize=12, res=900) 
#theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
Photo1 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoMeanMean)) +
  xlab("Elevation (m)") + ylab("Mean % Leaf Photo") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
#geom_smooth(method=lm)+
#geom_smooth()
#myplot_TDTCMean

#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean

ModelTDTNMean <- lm(Elevation.m. ~ PhotoMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
Photo2 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoVarianceMean)) +
  xlab("Elevation (m)") + ylab("Variance Leaf Photo") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoVarianceLower, ymax=PhotoVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0))
#geom_smooth(method=lm)+
#geom_smooth()
#myplot_TDTPhotoVar

ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTPhotoVar)
confint(ModelTDTPhotoVar)
coef(ModelTDTPhotoVar)

#Plot Skewness
Photo3 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoSkewnessMean)) +
  xlab("Elevation (m)") + ylab("Skewness Leaf Photo") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoSkewnessLower, ymax=PhotoSkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm) +
  #geom_smooth()+
  expand_limits(y=c(-1.0,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTCSkew 

ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTPhotoVar)
confint(ModelTDTPhotoVar)
coef(ModelTDTPhotoVar)

#Plot Kurtosis
Photo4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., (PhotoKurtosisMean-3))) +
  xlab("Elevation (m)") + ylab("Kurtosis Leaf Photo") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=(PhotoKurtosisLower-3), ymax=(PhotoKurtosisUpper-3)), width=.2, position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,10)) +
  geom_hline(yintercept = 0)
Photo4
#myplot_TDTCKurtosis 


ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoKurtosisMean, Peru_Plot_Master.data)
summary(ModelTDTPhotoVar)
confint(ModelTDTPhotoVar)
coef(ModelTDTPhotoVar)

# Multipanel plot
# multiplot(Photo1, Photo2, Photo3, Photo4, cols=2)
# dev.off()


################################################################################
#### SLA
################

#png("Figure_Sampled_SLA2.png", units="in", width=5, height=4, pointsize=12, res=900)

#theme_set(theme_gray(base_size = 15))

LMA1 <- myplot_SLANMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., 1/SLAMeanMean)) +
  xlab("Elevation (m)") +
  ylab(bquote('Subsampling Mean LMA ('~ g^-1~m^-2*')')) +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=1/SLAMeanLower, ymax=1/SLAMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()
LMA1
#myplot_SLANMean

LMA_data <- 1/(Peru_Plot_Master.data$SLAMeanMean)

ModelTDTLMAMean <- lm(Elevation.m. ~ LMA_data, Peru_Plot_Master.data)

summary(ModelTDTLMAMean)
confint(ModelTDTLMAMean)
coef(ModelTDTLMAMean)  

#Plot Variance
LMA2 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAVarianceMean)) +
  xlab("Elevation (m)") + ylab("Variance LMA") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=SLAVarianceLower, ymax=SLAVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
#geom_smooth(method=lm)+
#geom_smooth()
LMA2              

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
LMA3 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLASkewnessMean)) +
  xlab("Elevation (m)") + ylab("  Skewness LMA ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=SLASkewnessLower, ymax=SLASkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+
  expand_limits(y=c(-1.0,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTSkew

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis

LMA4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., (SLAKurtosisMean-3))) +
  xlab("Elevation (m)") + ylab("  Kurtosis LMA") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=(SLAKurtosisLower-3), ymax=(SLAKurtosisUpper-3)), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_smooth(method=lm)+
  #geom_smooth()+
  expand_limits(y=c(-10,10)) +
  geom_hline(yintercept = 0)
LMA4             

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
#multiplot(LMA1, LMA2, LMA3, LMA4, cols=2)
#dev.off()

#png("Figure_Sampled_Traits_Combined.png", units = "px", width=2400, height=1800, res=170)
#multiplot(N1, N2, N3, N4, C1, C2, C3, C4, P1, P2, P3, P4, Photo1, Photo2, Photo3, Photo4, LMA1, LMA2, LMA3, LMA4, cols=5)

png("Figure_Sampled_Traits_Combined.png", units = "px", width=2400, height=1800, res=170)
multiplot(N2, N3, N4, C2, C3, C4, P2, P3, P4, Photo2, Photo3, Photo4, LMA2, LMA3, LMA4, cols=5)
dev.off() 


####################  subsampled community mean P:N
          
PlotPtoNMean<- 1/Peru_Plot_Master.data$PlotNtoPMean
PN1 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PlotPtoNMean)) +
  xlab("Elevation (m)") +
  ylab("Mean % P:N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
PN1

PNUE_elev <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoPerLeafNMean)) +
  xlab("Elevation (m)") + 
  ylab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ m^-2~gN-1~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
PNUE_elev

PNUE_Photo <- myplot_PhotovPNUE <- ggplot(Peru_Plot_Master.data, aes(PhotoPerLeafNMean, PhotoMeanMean)) +
  xlab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ g^-1~s^-1*')')) + 
  ylab(bquote('Plot Leaf Photo. ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#myplot_PhotovPNUE
#wow - highly positive . . variation in mean photosynthesis appears to be driven by PNUE
PNUE_Photo
##
#Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean*100))

    
png("Figure_Sampled_Traits_Mean_Combined.png", units = "px", width=1800, height=1800, res=200)

multiplot(N1,P1,PN1,C1,LMA1,Photo1,PNUE_elev, PNUE_Photo, cols=3)
          
dev.off() 
          
 

#==========================================================================
##  MST predictions ln (NPP/M_Tot^0.6) ~ 1/kT ln(b)
##
#==========================================================================

Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_Malhi_2017)/(Peru_Plot_Master.data$MST_AGB1))

png("Figure_Plot_MST_Plots.png", units = "px", width=1900, height=600, res=300)

theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
MST <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=log(MST_GPP1))) +
  #xlab(italic("1/kT") + ylab(expression(ln(GPP/M[tot]^{0.6}))) +
  labs(x=expression(Temperature~1/kT),
       y=expression(ln(GPP/M[tot]^{0.6}))) +
  geom_point(size = 3, color="darkgrey") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=log(MST_GPP1Lower), ymax=log(MST_GPP1Upper)), width=0.03, position=position_dodge(0.05)) +
  geom_smooth(method = "lm", se=FALSE, color="black") +
  xlim(38.9, 41)+ylim(-0.5, 1.0) +
  geom_abline(intercept = 25.9, slope = -0.65, linetype=3)
MSTa <- MST + coord_fixed(ratio = 0.9)
MSTa
coef(lm(log(MST_GPP1) ~ MAinvBT, data = Peru_Plot_Master.data))
#(Intercept)     MAinvBT 
#6.9286608  -0.1576762 



#p <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, log(MST_GPP1))) + 
#geom_point() +
#p + geom_hline(yintercept = 0.8)
#p + geom_abline(intercept = 6.928, slope = -0.1)
#xlim(38.7, 40)+ylim(0.3, 1.3)
#p
#p + geom_abline(intercept = 26.5, slope = -0.65) +
#geom_smooth(method = "lm", se = FALSE)


library(grid)
grob <- grobTree(textGrob("A", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
MST1 <- MSTa + annotation_custom(grob)

mMST1 <- lm(log(MST_GPP1)~ MAinvBT, data=Peru_Plot_Master.data)
summary(mMST1) 
lm.sumSquares(mMST1)
AIC(mMST1)
AICc(mMST1)
modelEffectSizes(mMST1)  #lm.sunSquares is depreciated
varImp(mMST1, scale = FALSE)
#avPlots(mMST1)
#crPlots(mMST1)
confint(mMST1)
vif(mMST1) 


MST2 <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=log(MST_NPP_new1))) +
  # xlab(italic("1/kT") + ylab(expression(ln(NPP/M[tot]^{0.6}))) +
  labs(x=expression(Temperature~1/kT),
       y=expression(ln(NPP/M[tot]^{0.6}))) +
  geom_point(size = 3, color="darkgrey") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=log(MST_NPP1Lower), ymax=log(MST_NPP1Upper)), width=0.03, position=position_dodge(0.05)) +
  geom_smooth(method = "lm", se=FALSE, color="black") +
  xlim(38.9, 41)+ylim(-1.5, 0) +
  geom_abline(intercept = 24.95, slope = -0.65, linetype=3)
MST2a <- MST2 + coord_fixed(ratio = 0.9)
MST2a

library(grid)
grob <- grobTree(textGrob("B", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
MST3 <- MST2a + annotation_custom(grob)


mMST2 <- lm(log(MST_NPP_new1)~ MAinvBT, data=Peru_Plot_Master.data)
summary(mMST2) 
lm.sumSquares(mMST2)
AIC(mMST2)
AICc(mMST2)
modelEffectSizes(mMST2)  #lm.sunSquares is depreciated
varImp(mMST2, scale = FALSE)
#avPlots(mMST1)
#crPlots(mMST1)
confint(mMST2)
vif(mMST2) 


multiplot(MST1, MST3, cols=2)
dev.off()


png("Figure_Plot_MST_Plots_Celcius.png", units = "px", width=2400, height=1000, res=300)

theme_set(theme_gray(base_size = 15))

#--Figure 6a
panel_6a <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=MST_GPP1)) + 
  geom_abline(intercept = 25.9, slope = -0.65, linetype=3) + 
  geom_smooth(method = "lm", color="black") +
  geom_point(size=3, color="darkgrey") +
  geom_errorbar(aes(ymin=MST_GPP1Lower, ymax=MST_GPP1Upper), width=0.02, position=position_dodge(0.05)) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(paste("Inverse temperature, 1/kT (", eV^{-1}, ")"))) + 
  ylab(expression(paste("GPP / M " [tot]^{0.6}, "  (Mg C Mg ", M^{-0.6}, " ", ha^{-1}, yr^{-1}, ")"))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Temperature (", degree, C, ")")))) +
  theme_bw(base_size=12)
print(panel_6a)

library(grid)
grob <- grobTree(textGrob("A", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
panel_6_A <- panel_6a + annotation_custom(grob)


#--Figure 6b
panel_6b <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=MST_NPP_new1)) + 
  geom_abline(intercept = 24.95, slope = -0.65, linetype=3) +
  geom_smooth(method = "lm", color="black") +
  geom_point(size = 3, color="darkgrey") +
  geom_errorbar(aes(ymin=MST_NPP1Lower, ymax=MST_NPP1Upper), width=0.02, position=position_dodge(0.05)) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(paste("Inverse temperature, 1/kT (", eV^{-1}, ")"))) + 
  ylab(expression(paste("NPP / M " [tot]^{0.6}, "  (Mg C Mg ", M^{-0.6}, " ", ha^{-1}, yr^{-1}, ")"))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Temperature (", degree, C, ")")))) +
  #xlim(38.9, 40.8)+ylim(0.4, 3) +
  theme_bw(base_size=12)
print(panel_6b)

library(grid)
grob <- grobTree(textGrob("B", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
panel_6_B <- panel_6b + annotation_custom(grob)


grid.arrange(arrangeGrob(panel_6_A,panel_6_B,ncol=2))

dev.off()




#==========================================================================
# Assessing trait bimodality or unimodality 
# The bimodality coefficient (BC; SAS Institute, 1989)

# http://psych.nyu.edu/freemanlab/pubs/2012_BRM.pdf see also https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791391/
#The BC is based on an empirical relationship between bimodality and the third and fourth statistical moments of a distribution (skewness and kurtosis). It is proportional to the division of squared skewness with uncorrected kurtosis, BC ∝ (s2 + 1)/k, with the underlying logic that a bimodal distribution will have very low kurtosis, an asymmetric character, or both; all of these conditions increase BC. The values range from 0 and 1, with those exceeding .555 (the value representing a uniform distribution) suggesting bi- modality (SAS Institute, 1989).
# this is modified as BC ∝ (s2 + 1)/k+3

BC_LMA <- ((Peru_Plot_Master.data$SLASkewnessMean)^2 +1)/((Peru_Plot_Master.data$SLAKurtosisMean)+3)

BC_N <- ((Peru_Plot_Master.data$NSkewnessMean)^2 +1)/((Peru_Plot_Master.data$NKurtosisMean)+3)

BC_P <- ((Peru_Plot_Master.data$PSkewnessMean)^2 +1)/((Peru_Plot_Master.data$PKurtosisMean)+3)

BC_Photo <- ((Peru_Plot_Master.data$PhotoSkewnessMean)^2 +1)/((Peru_Plot_Master.data$PhotoKurtosisMean)+3)

BC_C <- ((Peru_Plot_Master.data$CSkewnessMean)^2 +1)/((Peru_Plot_Master.data$CKurtosisMean)+3)

png("Figure_Plot_BC_Plots.png", units = "px", width=1800, height=1800, res=200)

LMA_bimodal <- ggplot(Peru_Plot_Master.data, aes(x=Elevation.m., y=BC_LMA)) +
  xlab("Elevation m ") + ylab("LMA Bimodality coefficient, BC") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm) +
  geom_hline(yintercept = 0.555)+ 
  xlim(0, 3100)+ylim(0, 1) 
LMA_bimodal + coord_fixed(ratio = 0.9)
LMA_bimodal

N_bimodal <- ggplot(Peru_Plot_Master.data, aes(x=Elevation.m., y=BC_N)) +
  xlab("Elevation m ") + ylab("Nitrogen Bimodality coefficient, BC") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm) +
  geom_hline(yintercept = 0.555)+ 
  xlim(0, 3100)+ylim(0, 1) 
N_bimodal + coord_fixed(ratio = 0.9)
N_bimodal

P_bimodal <- ggplot(Peru_Plot_Master.data, aes(x=Elevation.m., y=BC_P)) +
  xlab("Elevation m ") + ylab("Phosphorus Bimodality coefficient, BC") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm) +
  geom_hline(yintercept = 0.555)+ 
  xlim(0, 3100)+ylim(0, 1) 
P_bimodal + coord_fixed(ratio = 0.9)
P_bimodal

C_bimodal <- ggplot(Peru_Plot_Master.data, aes(x=Elevation.m., y=BC_C)) +
  xlab("Elevation m ") + ylab("Carbon Bimodality coefficient, BC") +
  geom_point(size = 4, color="darkgrey") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method = "lm", se=FALSE, color="black") +
  geom_hline(yintercept = 0.555)+ 
  xlim(0, 3100)+ylim(0, 1) 
C_bimodal + coord_fixed(ratio = 0.9)
C_bimodal

Photo_bimodal <- ggplot(Peru_Plot_Master.data, aes(x=Elevation.m., y=BC_Photo)) +
  xlab("Elevation m ") + ylab("Photosynthesis Bimodality coefficient, BC") +
  geom_point(size = 4, color="darkgrey") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method = "lm", se=FALSE, color="black") +
  geom_hline(yintercept = 0.555)+ 
  xlim(0, 3100)+ylim(0, 1) 
Photo_bimodal + coord_fixed(ratio = 0.9)
Photo_bimodal

multiplot(LMA_bimodal, Photo_bimodal, C_bimodal, N_bimodal, P_bimodal, cols=2)
dev.off()


#==========================================================================
# Does the functinal response of b0 respond to tempeature? combine all traits in the model to extimate how b_0 changes 
## h*PNUE*SLA*P:N ~ b_0
# include Beta L? likely value is on the order of 0.00001

CUE <- Peru_Plot_Master.data$NPP_Malhi_2017/Peru_Plot_Master.data$GPP_Malhi_2017

png("Figure_b0_Temp.png", units = "px", width=900, height=600, res=170)

#rho_a <- (Peru_Plot_Master.data$PhotosynthesisPerLeafN * Peru_Plot_Master.data$lma.sun.mu.abundance * Peru_Plot_Master.data$PlotNtoP)

#rho_c <- (Peru_Plot_Master.data$PhotoPerLeafNMean * Peru_Plot_Master.data$SLAMeanMean * (1/Peru_Plot_Master.data$PlotNtoPMean))

#rho_c_cue <- (CUE * Peru_Plot_Master.data$PhotoPerLeafNMean * Peru_Plot_Master.data$SLAMeanMean * ((1/Peru_Plot_Master.data$PlotNtoPMean)))

# including an estimate of Beta_L in our estiamte of b_0. 
rho_c_cue_betaL <- (CUE * Peru_Plot_Master.data$PhotoPerLeafNMean * Peru_Plot_Master.data$SLAMeanMean * ((1/Peru_Plot_Master.data$PlotNtoPMean)) * 0.1)

#m_rho_c <- lm( log(rho_c) ~ MAinvBT, data=Peru_Plot_Master.data)
#summary(m_rho_c)

#m_rho_c_cue <- lm( log(rho_c_cue) ~ MAinvBT, data=Peru_Plot_Master.data)
#summary(m_rho_c_cue)

#m_rho_a <- lm( log(rho_a) ~ MAinvBT, data=Peru_Plot_Master.data)
#summary(m_rho_a)

m_rho_c_cue_betaL <- lm( log(rho_c_cue_betaL) ~ MAinvBT, data=Peru_Plot_Master.data)
summary(m_rho_c_cue_betaL)


rho_c_cue_plot  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, rho_c_cue_betaL)) +
  labs(x=expression(Temperature~1/kT),
       y=expression(paste("Allometric Tree Growth Normalization (", b[0],")"))) +
       #y=expression(Allometric Growth Normalization b[0]))+
  geom_point(size = 4, color="darkgrey") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  ylim(0, 0.15) 
#xlim(0.3, 1.2) +
#geom_smooth() 
#myplot_TDTPKurtosis 

rho_c_cue_plot

#multiplot(rho_c_plot, cols=1)
dev.off()


#==========================================================================
##  Additional analyses not reported in the paper
## assessing trait sampling vs. mean of abundant species to the 1:1 line
#==========================================================================

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


png("SampledvsMeasuredTraitMeans2.png", units="in", width=10,height=8, pointsize=12,res=500)

p1 <- ggplot(Peru_Plot_Master.data, aes(mean_n_percent, TransformedNpercent)) + 
  ylab("% Nitrogen Subsampling") + xlab("% Nitrogen Plot Abundant Species") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)
p1

Modelp1 <- lm(mean_n_percent ~ TransformedNpercent, Peru_Plot_Master.data)
summary(Modelp1)
confint(Modelp1)
coef(Modelp1)

p2 <- ggplot(Peru_Plot_Master.data, aes(mean_p_percent, TransformedPpercent)) +
  ylab("% Phosphorus Subsampling") + xlab("% Phosphorus Plot Abundant Species") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=TransformedPpercentLower, ymax=TransformedPpercentUpper), width=.001,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)
Modelp2 <- lm(mean_p_percent ~ TransformedPpercent, Peru_Plot_Master.data)
summary(Modelp2)
confint(Modelp2)
coef(Modelp2)

p3 <- ggplot(Peru_Plot_Master.data, aes(mean_c_percent, TransformedCpercent)) +
  ylab("% Carbon Subsampling") + xlab("% Carbon Plot Abundant Species") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=TransformedCpercentLower, ymax=TransformedCpercentUpper), width=.01,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)

Modelp3 <- lm(mean_c_percent ~ TransformedCpercent, Peru_Plot_Master.data)
summary(Modelp3)
confint(Modelp3)
coef(Modelp3)

p4 <- ggplot(Peru_Plot_Master.data, aes(mean_photosynthesis, PhotoMeanMean)) +
  
  ylab(bquote('Subsampling ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  xlab(bquote('Plot Photosynthesis Abundant Species('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.01,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)

Modelp4 <- lm(mean_photosynthesis ~ PhotoMeanMean, Peru_Plot_Master.data)
summary(Modelp4)
confint(Modelp4)
coef(Modelp4)


p5 <- ggplot(Peru_Plot_Master.data, aes(Transformedmean_sla_lamina_petiole, SLAMeanMean)) +
  ylab(bquote('Subsampling SLA ('~ m^2~kg^-1*')')) +
  xlab(bquote('Plot SLA Abundant Species ('~ m^2~kg^-1*')')) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=SLAMeanLower, ymax=SLAMeanUpper), width=.01,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)
Modelp5 <- lm(SLAMeanMean ~ Transformedmean_sla_lamina_petiole, Peru_Plot_Master.data)
summary(Modelp5)
confint(Modelp5)
coef(Modelp5)


p6 <- ggplot(Peru_Plot_Master.data, aes(PlotNtoP,PlotNtoPMean)) +
  ylab("% N:P Subsampling") + xlab("% N:P Plot Abundant Species") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedCpercentLower, ymax=TransformedCpercentUpper), width=.01,
  # position=position_dodge(0.05)) +
  geom_smooth(method=lm)

multiplot(p1, p2, p3, p4, p5, p6, cols=2)
dev.off()


