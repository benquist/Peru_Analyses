
##########################################
#### CHAMBASA - load data 

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged5.csv",header=T)

##########################################
#### CHAMBASA - load libraries

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

library(reshape2)
library(ggthemes)
library(plyr)

names(Peru_Plot_Master.data)
str(Peru_Plot_Master.data)


##########################################
#### CHAMBASA -  Calculated Variables ###############

##Calculate Boltzmann 1/kT
Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.+273.15))

##Calculate the leaf N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$mean_n_percent))

Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))

##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPP_newperNMeanMean <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$NMeanMean))

Peru_Plot_Master.data$GPPperNMeanMean <- ((Peru_Plot_Master.data$GPP)/(Peru_Plot_Master.data$NMeanMean))


#* calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))

##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$mean_n_percent)/ (Peru_Plot_Master.data$mean_p_percent))
Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

#Calculate NPP_newLeaf/RLeaf - Production per carbon respired
Peru_Plot_Master.data$NPP_LeafperRLeaf <- ((Peru_Plot_Master.data$NPP_newLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST prediction
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))


##########################################
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


############################################################################  
# PCA analyses
#######################################

#chambasapca <- princomp(~mean_sla_lamina_petiole + NMeanMean + PMeanMean + CMeanMean + mean_photosynthesis + PhotoPerLeafNMean + PlotNtoPMean, data=Peru_Plot_Master.data, cor=TRUE)

# Traits from most abundant individuals 
Photo <- Peru_Plot_Master.data$mean_photosynthesis
NP <- Peru_Plot_Master.data$PlotNtoP
P <- Peru_Plot_Master.data$mean_p_percent
N <- Peru_Plot_Master.data$mean_n_percent
C <- Peru_Plot_Master.data$mean_c_percent
SLA <- Peru_Plot_Master.data$mean_sla_lamina_petiole
PNUE <- Peru_Plot_Master.data$PhotosynthesisPerLeafN

## Subsampling total plot distribution method
Photo <- Peru_Plot_Master.data$PhotoMeanMean
NP <- Peru_Plot_Master.data$PlotNtoPMean
P <- Peru_Plot_Master.data$PMeanMean
N <- Peru_Plot_Master.data$NMeanMean
C <- Peru_Plot_Master.data$CMeanMean
SLA <- Peru_Plot_Master.data$SLAMeanMean
PNUE <- Peru_Plot_Master.data$PhotoPerLeafNMean


chambasapca <- princomp(~ SLA + N + C + P + Photo + NP + PNUE, data=Peru_Plot_Master.data, cor=TRUE)

summary(chambasapca)
loadings(chambasapca)
loadings <- loadings(chambasapca)
#write.csv(loadings, "PCA_CommunityTrait.csv")

png("Figure_PCA2_gradient_plot.png", units="in", width=5, height=4, pointsize=9, res=900)
PCA <- biplot(chambasapca, col=c("gray","red"),cex=c(0.8,0.5))
PCA
dev.off()

screeplot(chambasapca)
dev.off()
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

#merge(elev,first.component.scores) 

png("Figure_PCA_gradient_Abund.png", units = "px", width=1800, height=1800, res=300)

#########  Best predictors of PCA1 and PCA2 variation, 
#PCAfit <- glmulti(PCA1Scores ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + Elevation.m. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

PCAfit <- glmulti(PCA1ScoresPlotTraits ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit)
tmp <- weightable(PCAfit)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit@objects[[1]])
plot(PCAfit)
#Variable Importance
plot(PCAfit, type="s")
#PCA1 is mainly temperature

PCA2fit <- glmulti(PCA2ScoresPlotTraits ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCA2fit)
tmp <- weightable(PCA2fit)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCA2fit@objects[[1]])
plot(PCA2fit)
#Variable Importance
plot(PCA2fit, type="s")


## Whole community trait variation in PCA space
PCAfitSample <- glmulti(PCA1ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfitSample)
tmp <- weightable(PCAfitSample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfitSample@objects[[1]])
plot(PCAfitSample)
#Variable Importance
plot(PCAfitSample, type="s")
#PCA1 is mainly solar radiation

PCAfit2Sample <- glmulti(PCA2ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit2Sample)
tmp <- weightable(PCAfit2Sample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit2Sample@objects[[1]])
plot(PCAfit2Sample)
#Variable Importance
plot(PCAfit2Sample, type="s")
#PCA1 is mainly solar radiation

#########################################################
########### Plots of PCA axes vs. environmental variation
PCA1 <- myplot_PCA1vElev<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PCA1ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA1vElev

ModelPCA1 <- lm(PCA1ScoresPlotTraits ~ MeanAnnualAirTemperature.degC., data=Peru_Plot_Master.data)
summary(ModelPCA1)  # y ~ x


PCA2 <- myplot_PCA2vElev<- ggplot(Peru_Plot_Master.data, aes(Precipitation.mm.yr.1., PCA2ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA2vElev

Model_PCA2 <- lm(PCA2ScoresPlotTraits ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
summary(Model_PCA2)

png("Figure_PCA_multi.png", units = "px", width=1500, height=900, res=200)

multiplot(PCA1, PCA2, cols=2)
dev.off()

#############################################
#Best predictors of PCA1 and PCA2 variation, climate?

#PCAfit <- glmulti(PCA1Scores ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + Elevation.m. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

PCAfit <- glmulti(PCA1ScoresPlotTraits ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit)
#tmpPCA1 <- weightable(PCAfit)
#tmpPCA1 <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
#tmpPCA1
summary(PCAfit@objects[[1]])
summary(PCAfit@objects[[2]])
plot(PCAfit)
#Variable Importance
plot(PCAfit, type="s")
#PCA1 is mainly temperature


PCA2fit <- glmulti(PCA2ScoresPlotTraits ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCA2fit)
#tmp <- weightable(PCA2fit)
#tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
#tmp
summary(PCA2fit@objects[[1]])
summary(PCA2fit@objects[[2]])

plot(PCA2fit)
#Variable Importance
plot(PCA2fit, type="s")


PCAfitSample <- glmulti(PCA1ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfitSample)
tmp <- weightable(PCAfitSample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfitSample@objects[[1]])
plot(PCAfitSample)
#Variable Importance
plot(PCAfitSample, type="s")
#PCA1 is mainly solar radiation

PCAfit2Sample <- glmulti(PCA2ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit2Sample)
tmp <- weightable(PCAfit2Sample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit2Sample@objects[[1]])
plot(PCAfit2Sample)
#Variable Importance
plot(PCAfit2Sample, type="s")
#PCA1 is mainly solar radiation



########################################################################
########################################################################
###  Relationship between dominant community species traits and elevation
###

png("Figure_Plot_Traits.png", units = "px", width=1800, height=1800, res=300)
#dev.off()
T1 <- myplot_TDTNAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_n_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf N") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
T1
Model_T1 <- lm(mean_n_percent ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T1)


T2 <- myplot_TDTPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_p_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf P") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
#myplot_TDTPAve
T2
Model_T2 <- lm(mean_p_percent ~ Elevation.m., data=Peru_Plot_Master.data)
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
Model_T3 <- lm(PlotNtoP ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T3)


T4 <- myplot_TDTCAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_c_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf C ") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
#geom_smooth()
#myplot_TDTCAve
T4
Model_T4 <- lm(mean_c_percent ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T4)



T5 <- myplot_TDTSLAAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_sla_lamina_petiole)) +
  xlab("Elevation (m)") + 
  ylab(bquote('Plot SLA ('~ m^2~kg^-1*')'))+
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) 
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
#myplot_TDTSLAAve
T5
Model_T5 <- lm(mean_sla_lamina_petiole ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T5)


T6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_photosynthesis)) +
  xlab("Elevation (m)") + 
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
#myplot_TDTPhotoAve
T6
Model_T6 <- lm(mean_photosynthesis ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T6)

T7 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotosynthesisPerLeafN)) +
  xlab("Elevation (m)") + 
  ylab(bquote('PNUE ('*mu~ 'mol' ~CO[2]~ g^-1~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTPhotoAve
T7
Model_T7 <- lm( PhotosynthesisPerLeafN ~ Elevation.m., data=Peru_Plot_Master.data)
summary(Model_T7)

T8 <- myplot_PhotovPNUE <- ggplot(Peru_Plot_Master.data, aes(PhotosynthesisPerLeafN, mean_photosynthesis)) +
  
  xlab(bquote('Plot Leaf N efficiency ('*~CO[2]~ g^-1 m^2 ~s^-1*')')) +
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

Model_T8 <- lm(mean_photosynthesis ~ PhotosynthesisPerLeafN, data=Peru_Plot_Master.data)
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

T10 <- myplot_NtoPvPNUE <- ggplot(Peru_Plot_Master.data, aes(mean_p_percent, PhotosynthesisPerLeafN)) +
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

multiplot(T1,T2,T3,T4,T5,T6,T7,T8, cols=3)
dev.off()


##################################
###### Distributions from subsampled trait distribuitons
##################################


########################################################
########################################################
# Nitrogen 
#
######

png("Figure_Sampled_Nitrogen2.png", units="in", width=5, height=4, pointsize=12, res=900)
#dev.off()
theme_set(theme_gray(base_size = 15))
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
  #geom_smooth(method=lm)+
  geom_smooth()
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
  #geom_smooth(method=lm)+
  geom_smooth()
#myplot_TDTNVar

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
  #geom_smooth(method=lm)+
  geom_smooth()+
  expand_limits(y=c(-2.5,2.5)) +
  geom_hline(yintercept = 0)
#myplot_TDTSkew 

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
N4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NKurtosisMean)) +
  xlab("Elevation (m)") + ylab("  Kurtosis % N ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NKurtosisLower, ymax=NKurtosisUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)+
  geom_smooth()+ 
  expand_limits(y=c(-10,30)) +
  geom_hline(yintercept = 0)
#myplot_TDTKurtosis 

#myplot_TDTKurtosis

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
multiplot(N1, N2, N3, N4, cols=2)
dev.off()


#####################################
#####################################
### assessing TDT predictions  
# Nitrogen v Temperature

png("Nitrogen_MST2", units="in", width=5, height=4, pointsize=12, res=900)
theme_set(theme_gray(base_size = 15)
          #theme_set(theme_classic(base_size = 30))
          N1 <- myplot_MSTtempvN <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) +
            xlab(bquote('Average annual temp., <1/kT> ('*~ EV^-1*')')) +
            ylab("  % Leaf Nitrogen") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.01,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            geom_smooth(method=lm)
          #myplot_TDTNMean
          
          ModelTDTNMean <- lm(MAinvBT ~ NMeanMean, Peru_Plot_Master.data)
          summary(ModelTDTNMean)
          confint(ModelTDTNMean)
          coef(ModelTDTNMean)
          
          
          #Plot Variance
          N2 <- myplot_MSTtempvVarN  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NVarianceMean)) +
            xlab(bquote('Average annual temp., <1/kT> ('*~ EV^-1*')')) +
            ylab("  Variance % Leaf Nitrogen") +
            geom_point(size = 3, color="red") +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            geom_errorbar(aes(ymin=NVarianceLower, ymax=NVarianceUpper), width=.01,
                          position=position_dodge(0.05)) +
            geom_smooth(method=lm)
          #myplot_MSTtempvVarN 
          #myplot_TDTNVar
          
          ModelTDTNVar <- lm(NVarianceMean ~ log(MAinvBT), Peru_Plot_Master.data)
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Skewness
          N3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NSkewnessMean)) +
            xlab(bquote('Average annual temp., <1/kT> ('*~ EV^-1*')')) +
            ylab("  Skewness % Leaf Nitrogen") +
            geom_point(size = 3, color="red") +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            geom_errorbar(aes(ymin=NSkewnessLower, ymax=NSkewnessLower), width=.01,
                          position=position_dodge(0.05)) +
            expand_limits(y=c(-2.5,4)) +
            geom_hline(yintercept = 0) +
            geom_smooth()
          #myplot_TDTSkew
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Kurtosis
          N4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NKurtosisMean)) +
            geom_point(size = 3, color="red")+
            xlab(bquote('Average annual temp., <1/kT> ('*~ EV^-1*')')) +
            ylab("  Kurtosis % Leaf Nitrogen") +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            geom_errorbar(aes(ymin=NKurtosisLower, ymax=NKurtosisUpper), width=.05, position=position_dodge(0.01))+
            expand_limits(y=c(-10,30)) +
            geom_hline(yintercept = 0) +
            #geom_smooth(method=lm)
            geom_smooth() 
          #myplot_TDTKurtosis
          
          ModelTDTNVar <- lm(MAinvBT ~ NVarianceMean, Peru_Plot_Master.data)
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          # Multipanel plot
          multiplot(N1, N2, N3, N4, cols=2)
          dev.off()
          
          #* shifts in mean and variance with elevation . . . all skewness values are greater than zero indicating that either all plots are shifting or there is asymetric colonization into plots. Also, all plots have postive kurtosis - indicative of strong stabilizing filtering around an optimal value . .  TDT would predict distributions that have more positive kurtosis
          
          
################################################################
################################################################
#### Leaf Carbon
################
          
          png("Figure_Sampled_Carbon2.png", units="in", width=5, height=4, pointsize=12, res=900)
          
          theme_set(theme_gray(base_size = 15))
          #theme_set(theme_classic(base_size = 30))
          C1 <- myplot_TDTCMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CMeanMean)) +
            xlab("Elevation (m)") + ylab("  Mean % Leaf C") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTCMean
          
          #myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
          #myplot_TDTNMean
          
          ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNMean)
          confint(ModelTDTNMean)
          coef(ModelTDTNMean)
          
          
          #Plot Variance
          C2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CVarianceMean)) +
            xlab("Elevation (m)") + ylab("  Variance % C ") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=CVarianceLower, ymax=CVarianceUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTCVar
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Skewness
          C3 <- myplot_TDTCSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CSkewnessMean)) +
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
            geom_smooth()+
            expand_limits(y=c(-2.5,2.5)) +
            geom_hline(yintercept = 0)
          #myplot_TDTCSkew 
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Kurtosis
          C4 <- myplot_TDTCKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CKurtosisMean)) +
            xlab("Elevation (m)") + ylab("  Kurtosis % C ") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=CKurtosisLower, ymax=CKurtosisUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
            #position=position_dodge(0.05)) +
            #geom_smooth(method=lm)+
            geom_smooth()+ 
            expand_limits(y=c(-10,30)) +
            geom_hline(yintercept = 0)
          #myplot_TDTCKurtosis 
          
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          # Multipanel plot
          multiplot(C1, C2, C3, C4, cols=2)
          dev.off()
          
          ################################################################
          ################################################################
          #### Leaf Phosphorus
          ################
          
          png("Figure_Sampled_Phosphorus2.png", units="in", width=5, height=4, pointsize=12, res=900)
          
          theme_set(theme_gray(base_size = 15))
          #theme_set(theme_classic(base_size = 30))
          P1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PMeanMean)) +
            xlab("Elevation (m)") + ylab("  Mean % Leaf P") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=PMeanLower, ymax=PMeanUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTPMean
          
          #myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
          #myplot_TDTNMean
          
          ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNMean)
          confint(ModelTDTNMean)
          coef(ModelTDTNMean)
          
          
          #Plot Variance
          P2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PVarianceMean)) +
            xlab("Elevation (m)") + ylab("  Variance % P ") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=PVarianceLower, ymax=PVarianceUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTPVar
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Skewness
          P3 <- myplot_TDTPSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PSkewnessMean)) +
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
            geom_smooth()+
            expand_limits(y=c(-2.5,2.5)) +
            geom_hline(yintercept = 0)
          #myplot_TDTPSkew 
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Kurtosis
          P4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PKurtosisMean)) +
            xlab("Elevation (m)") + ylab("Kurtosis % P ") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=PKurtosisLower, ymax=PKurtosisUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
            #position=position_dodge(0.05)) +
            #geom_smooth(method=lm)+
            geom_smooth()+ 
            expand_limits(y=c(-10,30)) +
            geom_hline(yintercept = 0)
          #myplot_TDTPKurtosis 
          
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          # Multipanel plot
          multiplot(P1, P2, P3, P4, cols=2)
          dev.off()
          
################################################################
################################################################
#### Leaf Photosynthesis
################
          
png("Figure_Sampled_Photo.png", units="in", width=5, height=4, pointsize=12, res=900) 
          
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
Photo1 <- myplot_TDTPhotoMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoMeanMean)) +
            xlab("Elevation (m)") + ylab("Mean % Leaf Photo") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTCMean
          
          #myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
          #myplot_TDTNMean
          
          ModelTDTNMean <- lm(Elevation.m. ~ PhotoMeanMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNMean)
          confint(ModelTDTNMean)
          coef(ModelTDTNMean)
          
          
          #Plot Variance
          Photo2 <- myplot_TDTPhotoVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoVarianceMean)) +
            xlab("Elevation (m)") + ylab("Variance Leaf Photo") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=PhotoVarianceLower, ymax=PhotoVarianceUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTPhotoVar
          
          ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTPhotoVar)
          confint(ModelTDTPhotoVar)
          coef(ModelTDTPhotoVar)
          
          #Plot Skewness
          Photo3 <- myplot_TDTPhotoSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoSkewnessMean)) +
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
            #geom_smooth(method=lm)+
            geom_smooth()+
            expand_limits(y=c(-2.5,2.5)) +
            geom_hline(yintercept = 0)
          #myplot_TDTCSkew 
          
          ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTPhotoVar)
          confint(ModelTDTPhotoVar)
          coef(ModelTDTPhotoVar)
          
 #Plot Kurtosis
Photo4 <- myplot_TDTPhotoKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Kurtosis Leaf Photo") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoKurtosisLower, ymax=PhotoKurtosisUpper), width=.2, position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
            #position=position_dodge(0.05)) +
            #geom_smooth(method=lm)+
            geom_smooth()+ 
            expand_limits(y=c(-10,30)) +
            geom_hline(yintercept = 0)
          #myplot_TDTCKurtosis 
          
          
          ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoVarianceMean, Peru_Plot_Master.data)
          summary(ModelTDTPhotoVar)
          confint(ModelTDTPhotoVar)
          coef(ModelTDTPhotoVar)
          
          # Multipanel plot
          multiplot(Photo1, Photo2, Photo3, Photo4, cols=2)
          dev.off()
          
          
          
################
#### SLA
################
          
png("Figure_Sampled_SLA2.png", units="in", width=5, height=4, pointsize=12, res=900)
          
theme_set(theme_gray(base_size = 15))
          
SLA1 <- myplot_SLANMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAMeanMean)) +
            xlab("Elevation (m)") +
            ylab(bquote('Mean SLA ('~ m^2~kg^-1*')')) +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=SLAMeanLower, ymax=SLAMeanUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_SLANMean
          ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNMean)
          confint(ModelTDTNMean)
          coef(ModelTDTNMean)  
          
          #Plot Variance
          SLA2 <- myplot_TDTSLAVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAVarianceMean)) +
            xlab("Elevation (m)") + ylab("Variance SLA") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=SLAVarianceLower, ymax=SLAVarianceUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_smooth(method=lm)+
            geom_smooth()
          #myplot_TDTSLAVar
          
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Skewness
          SLA3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLASkewnessMean)) +
            xlab("Elevation (m)") + ylab("  Skewness SLA ") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=SLASkewnessLower, ymax=SLASkewnessUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
            #position=position_dodge(0.05)) +
            #geom_smooth(method=lm)+
            geom_smooth()+
            expand_limits(y=c(-2.5,2.5)) +
            geom_hline(yintercept = 0)
          #myplot_TDTSkew
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          #Plot Kurtosis
          SLA4 <- myplot_TDTSLAKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAKurtosisMean)) +
            xlab("Elevation (m)") + ylab("  Kurtosis SLA") +
            geom_point(size = 3, color="red") +
            geom_errorbar(aes(ymin=SLAKurtosisLower, ymax=SLAKurtosisUpper), width=.2,
                          position=position_dodge(0.05)) +
            theme_bw() +
            theme(axis.text = element_text(size = 9),
                  axis.text.y = element_text(size = rel(1.5), angle = 0)) +
            theme(axis.text = element_text(size = 9),
                  axis.text.x = element_text(size = rel(1.5), angle = 0)) +
            #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
            #position=position_dodge(0.05)) +
            #geom_smooth(method=lm)+
            geom_smooth()+
            expand_limits(y=c(-10,30)) +
            geom_hline(yintercept = 0)
          
          myplot_TDTSLAKurtosis
          
          ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
          summary(ModelTDTNVar)
          confint(ModelTDTNVar)
          coef(ModelTDTNVar)
          
          # Multipanel plot
          multiplot(SLA1, SLA2, SLA3, SLA4, cols=2)
          dev.off()
          
###############################################################################
#########
########
          ##########################
library(glmulti)
library(leaps)
library(MASS)
library(lme4)        
library(MuMIn)
          
# Here are the main models hypothesized . . .
          m1 <- lm(log(GPP_new_estimate) ~ MAinvBT, data=Peru_Plot_Master.data) #climate models
          m2 <- lm(log(GPP_new_estimate) ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
          m3 <- lm(log(GPP_new_estimate) ~ log(Precipitation.mm.yr.1.), data=Peru_Plot_Master.data)
          m4 <- lm(log(GPP_new_estimate) ~ log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m5 <- lm(log(GPP_new_estimate) ~ log(SolarRadiation.GJ.m.2.yr.1.), data=Peru_Plot_Master.data)
          m6 <- lm(log(GPP_new_estimate) ~ log(SolarRadiation.GJ.m.2.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m7 <- lm(log(GPP_new_estimate) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m8 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data) #Biomass model
          m9 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data) # Biomass and climate
          m10 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m11 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data) # TDT Biomass and Traits
          m12 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)
          m13 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
          m14 <- lm(log(GPP_new_estimate) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
          m15 <- lm(log(GPP_new_estimate) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
          m16 <- lm(log(GPP_new_estimate) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(NVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
          m17 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
          m18 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN) + log(PlotNtoP), data=Peru_Plot_Master.data)
          
          
          
          TableGPPAICc <-  model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = aicc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
          write.table(TableGPPAICc, "TableGPP_AICc.xls", sep = "\t")
          
          # worst models m7, m15, m18, m16, m2, m3, m6, . . .  # Best models m12, m9,  then m8, m19, m14, m11,  
          #to pull out other info see 
          # summary(fit)$coefficients[,4]   summary(fit)$r.squared  summary(m12)$adj.r.squared
          
          # using QAIC by changing rank statement to rank = QAIC, . . . best models m16, m12, m15, m7
          
          TableGPPQAICc <- model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = QAIC, rank.args = list(chat = deviance(m12) / df.residual(m12)))
          
          write.table(TableGPPQAICc, "TableGPPQAICc.xls", sep = "\t")
          
          
          ##### NPP ######
          # Here are the main models hypothesized . . .
          m1 <- lm(log(NPP_new) ~ MAinvBT, data=Peru_Plot_Master.data) #climate models
          m2 <- lm(log(NPP_new) ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
          m3 <- lm(log(NPP_new) ~ log(Precipitation.mm.yr.1.), data=Peru_Plot_Master.data)
          m4 <- lm(log(NPP_new) ~ log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m5 <- lm(log(NPP_new) ~ log(SolarRadiation.GJ.m.2.yr.1.), data=Peru_Plot_Master.data)
          m6 <- lm(log(NPP_new) ~ log(SolarRadiation.GJ.m.2.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m7 <- lm(log(NPP_new) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m8 <- lm(log(NPP_new) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data) #Biomass model
          m9 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data) # Biomass and climate
          m10 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m11 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data) # TDT Biomass and Traits
          m12 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)
          m13 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
          m14 <- lm(log(NPP_new) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
          m15 <- lm(log(NPP_new) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
          m16 <- lm(log(NPP_new) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(NVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
          m17 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
          m18 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN) + log(PlotNtoP), data=Peru_Plot_Master.data)
          
          
          TableAICcNPP <-model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = aicc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
          # worst models m15, m16, m7, m18, m12 . . .  # Best models m9, m14, m10  but note that m16 has the highest logLik value
          
          write.table(TableAICcNPP, "TableAICcNPPc.xls", sep = "\t")
          
          model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = logLik, rank.args = list(chat = deviance(m12) / df.residual(m12)))
          # highest logLik m16, m15, m15, m12  lowest logLik  m2, m3, m1, m8
          
          TableQAICcNPP <-model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = QAICc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
          
          write.table(TableQAICcNPP, "TableQAICcNPPc.xls", sep = "\t")
          
          # worst models m15, m16, m7, m18, m12 . . .  # Best models m9, m14, m10
          #summary(m12)$adj.r.squared
          
          m1 <- lm(log(NPP_new) ~  log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)
          m2 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
          m3 <- lm(log(NPP_new) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
          m4 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          m5 <- lm(log(NPP_new) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
          m6 <- lm(log(NPP_new) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
          m7 <- lm(log(NPP_new) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
          
          
          library(MuMIn)
          model.sel(m1,m2,m3,m4,m6,m7, rank = qaicc(), rank.args = list(chat = deviance(m2) / df.residual(m2)))
          
          
          
          
####################################################################
##  MST predictions ln (NPP/M_Tot^0.6) ~ 1/kT ln(b)
##
###################################################################

Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB1))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))

png("Figure_Plot_MST_Plots.png", units = "px", width=1900, height=600, res=300)
          
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
MST <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=log(MST_GPP1))) +
#xlab(italic("1/kT") + ylab(expression(ln(GPP/M[tot]^{0.6}))) +
xlab("Inverse Temperature (1/kT)") + ylab("ln(GPP/M[tot]^{0.6}") +
geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
MST + coord_fixed(ratio = 0.9)
MST

library(grid)
grob <- grobTree(textGrob("A", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
MST1 <- MST + annotation_custom(grob)



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
xlab("Inverse Temperature (1/kT)") + ylab("ln(NPP/M[tot]^{0.6}") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
MST2 + coord_fixed(ratio = 0.9)
library(grid)
grob <- grobTree(textGrob("B", x=0.9,  y=0.85, hjust=0,
                          gp=gpar(col="Black", fontsize=15, fontface="bold")))
# Plot
MST3 <- MST2 + annotation_custom(grob)

                        
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



