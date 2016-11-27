########################################################################
#  Peru CHAMBASA TDT paper
#   Brian J. Enquist
#  10/21/16
########################################################################
 

#load data 
Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged6.csv",header=T)

update.packages(checkBuilt = TRUE)

#load libraries
install.packages ("ggplot2", dependencies = TRUE)
install.packages("plyr")
install.packages("ggthemes")
install.packages("reshape2")
install.packages("glmulti")
install.packages("grid")
install.packages("gtable")

# << NOTE: Make sure to wait until ggplot installation says "DONE" before proceeding!!! >>
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_1.0.0.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
# verify installation of v1.0.0
#packageDescription("ggplot2")$Version

library(smatr)
library(car)
library(lmSupport)
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape)
library(gtable)
library(AICcmodavg)
library(gtable)
library(gridExtra)
library(glmulti)

library(reshape2)
library(ggthemes)
library(plyr)

names(Peru_Plot_Master.data)
str(Peru_Plot_Master.data)

##################### Calculated Variables #######################

#Calculate Boltzmann 1/kT
Peru_Plot_Master.data$MAinvBT <- 1/(0.00008617*(Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.+273.15))

Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB1))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))

##Calculate the leaf PNUE or N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$PhotosynthesisPerLeafN <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))
#Peru_Plot_Master.data$PhotoPerLeafNMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))
 
      # calculation amax per unit mass
Peru_Plot_Master.data$amax.sun.massSpecific <- ((Peru_Plot_Master.data$amax.sun.mu.abundance)*(1/(Peru_Plot_Master.data$lma.sun.mu.abundance)))

    #calculate PNUE on a mass basis 
Peru_Plot_Master.data$Photosynthesis_massPerLeafN <- ((Peru_Plot_Master.data$amax.sun.massSpecific)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))


##Calculate the plot N productivity umol/m^2/s divided by foliar N
Peru_Plot_Master.data$NPP_newperNMeanMean <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))

Peru_Plot_Master.data$GPPperNMeanMean <- ((Peru_Plot_Master.data$GPP)/(Peru_Plot_Master.data$n_percent.sun.mu.abundance))


#Calculate leaf carbon efficiency per leaf first before calculating the plot average?

##Calculate site leaf carbon production efficiency
#Peru_Plot_Master.data$PhotosynthesisPerRLeaf <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$RLeaf))

##Calculate site N:P 
Peru_Plot_Master.data$PlotNtoP <- ((Peru_Plot_Master.data$n_percent.sun.mu.abundance)/ (Peru_Plot_Master.data$p_corrected_percent.sun.mu.abundance))

Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

#Calculate NPP_newLeaf/RLeaf - Production per carbon respired
Peru_Plot_Master.data$NPP_LeafperRLeaf <- ((Peru_Plot_Master.data$NPP_newLeaf)/ (Peru_Plot_Master.data$RLeaf))

#MST prediction
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))


##########################################
#### Plotting Functions
### Panel Plot - define multiplot function first
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

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


### ggcorplot function
library(ggplot2)

#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
  for(i in rev(new_order)){
    x=relevel(x,ref=i)
  }
  return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
  # normalize data
  for(i in 1:length(data)){
    data[,i]=(data[,i]-mean(data[,i]))/sd(data[,i])
  }
  # obtain new data frame
  z=data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      temp=as.data.frame(cbind(x,y))
      temp=cbind(temp,names(data)[i],names(data)[j])
      z=rbind(z,temp)
      j=j+1
    }
  }
  names(z)=c('x','y','x_lab','y_lab')
  z$x_lab = ezLev(factor(z$x_lab),names(data))
  z$y_lab = ezLev(factor(z$y_lab),names(data))
  z=z[z$x_lab!=z$y_lab,]
  #obtain correlation values
  z_cor = data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      x_mid = min(x)+diff(range(x))/2
      y_mid = min(y)+diff(range(y))/2
      this_cor = cor(x,y)
      this_cor.test = cor.test(x,y)
      this_col = ifelse(this_cor.test$p.value<.05,'<.05','>.05')
      this_size = (this_cor)^2
      cor_text = ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b=as.data.frame(cor_text)
      b=cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor=rbind(z_cor,b)
      j=j+1
    }
  }
  names(z_cor)=c('cor','x_mid','y_mid','p','rsq','x_lab','y_lab')
  z_cor$x_lab = ezLev(factor(z_cor$x_lab),names(data))
  z_cor$y_lab = ezLev(factor(z_cor$y_lab),names(data))
  diag = z_cor[z_cor$x_lab==z_cor$y_lab,]
  z_cor=z_cor[z_cor$x_lab!=z_cor$y_lab,]
  #start creating layers
  points_layer = layer(
    geom = 'point'
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_line_layer = layer(
    geom = 'line'
    , geom_params = list(colour = 'red')
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_ribbon_layer = layer(
    geom = 'ribbon'
    , geom_params = list(fill = 'green', alpha = .5)
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  cor_text = layer(
    geom = 'text'
    , data = z_cor
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=cor
      , size = rsq
      , colour = p
    )
  )
  var_text = layer(
    geom = 'text'
    , geom_params = list(size=var_text_size)
    , data = diag
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=x_lab
    )
  )
  f = facet_grid(y_lab~x_lab,scales='free')
  o = opts(
    panel.grid.minor = theme_blank()
    ,panel.grid.major = theme_blank()
    ,axis.ticks = theme_blank()
    ,axis.text.y = theme_blank()
    ,axis.text.x = theme_blank()
    ,axis.title.y = theme_blank()
    ,axis.title.x = theme_blank()
    ,legend.position='none'
  )
  size_scale = scale_size(limits = c(0,1),to=cor_text_limits)
  return(
    ggplot()+
      points_layer+
      lm_ribbon_layer+
      lm_line_layer+
      var_text+
      cor_text+
      f+
      o+
      size_scale
  )
}

#set up some fake data
library(MASS)
N=100

#first pair of variables
variance1=1
variance2=2
mean1=10
mean2=20
rho = .8
Sigma=matrix(c(variance1,sqrt(variance1*variance2)*rho,sqrt(variance1*variance2)*rho,variance2),2,2)
pair1=mvrnorm(N,c(mean1,mean2),Sigma,empirical=T)

#second pair of variables
variance1=10
variance2=20
mean1=100
mean2=200
rho = -.4
Sigma=matrix(c(variance1,sqrt(variance1*variance2)*rho,sqrt(variance1*variance2)*rho,variance2),2,2)
pair2=mvrnorm(N,c(mean1,mean2),Sigma,empirical=T)

my_data=data.frame(cbind(pair1,pair2))

ggcorplot(
  data = my_data
  , var_text_size = 30
  , cor_text_limits = c(2,30)
)

############################################################################  
# PCA analyses
#######################################

#chambasapca <- princomp(~mean_sla_lamina_petiole + NMeanMean + PMeanMean + CMeanMean + mean_photosynthesis + PhotoPerLeafNMean + PlotNtoPMean, data=Peru_Plot_Master.data, cor=TRUE)

# Traits from most abundant individuals 
Photo_a <- Peru_Plot_Master.data$amax.sun.mu.abundance
PN_a <- 1/(Peru_Plot_Master.data$PlotNtoP)
P_a <- Peru_Plot_Master.data$p_corrected_percent.sun.mu.abundance
N_a <- Peru_Plot_Master.data$n_percent.sun.mu.abundance
C_a <- Peru_Plot_Master.data$c_percent.sun.mu.abundance
LMA_a <- Peru_Plot_Master.data$lma.sun.mu.abundance
PNUE_a <- Peru_Plot_Master.data$PhotosynthesisPerLeafN

## Subsampling total plot subsampling distribution method
Photo <- Peru_Plot_Master.data$PhotoMeanMean
PN <- 1/(Peru_Plot_Master.data$PlotNtoPMean)
P <- Peru_Plot_Master.data$PMeanMean
N <- Peru_Plot_Master.data$NMeanMean
C <- Peru_Plot_Master.data$CMeanMean
LMA <- 1/(Peru_Plot_Master.data$SLAMeanMean)
PNUE <- Peru_Plot_Master.data$PhotoPerLeafNMean


#chambasapca <- princomp(~ LMA + N + C + P + Photo + PN + PNUE, data=Peru_Plot_Master.data, cor=TRUE)

chambasapca <- princomp(~ LMA_a + N_a + C_a + P_a + Photo_a + PN_a + PNUE_a, data=Peru_Plot_Master.data, cor=TRUE)

summary(chambasapca)
loadings(chambasapca)
loadings <- loadings(chambasapca)
#write.csv(loadings, "PCA_CommunityTrait.csv")
#* first two principle components explain about 72% of the variation.

png("Figure_PCA_Gradient_plot.png", units="in", width=5, height=4, pointsize=9, res=900)
PCA <- biplot(chambasapca, col=c("gray","red"),cex=c(0.8,0.5))
PCA
dev.off()

screeplot(chambasapca)
#dev.off()
chambasapca$scores
chambasapca$loadings
 

first.component.scores <- chambasapca$scores[,1]
summary(first.component.scores)
length(first.component.scores)
write.csv(first.component.scores, "First.component.scores.abund.csv")


second.component.scores <- chambasapca$scores[,2]
summary(second.component.scores)
length(second.component.scores)
write.csv(second.component.scores, "Second.component.scores.abund.csv")


#################################
#########  Best predictors of PCA1 and PCA2 variation, 
#PCAfit <- glmulti(PCA1Scores ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + Elevation.m. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

PCA_fit <- glmulti(PCA1ScoresPlotTraits_a~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp + Aspect..deg.+ Slope..deg.+ Soil.moisture.... + Elevation.m., data = Peru_Plot_Master.data, crit=BIC, level=1, fitfunc=glm, method="h")
# prior had method ="h"
summary(PCA_fit)
tmp <- weightable(PCA_fit)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCA_fit@objects[[1]])
plot(PCA_fit)
print(PCA_fit)
#Variable Importance
plot(PCA_fit, type="s")
summary(PCA_fit@objects[[1]])

### using PCA1ScoresPlotTraits_a PCA 1 is maintly soil moisture

PCA2fit <- glmulti(PCA2ScoresPlotTraits ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + MeanAnnualAirTemperature.degC. + mean_air_temp + Aspect..deg.+ Slope..deg.+ Soil.moisture.... + Elevation.m.,, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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


PCAfit2Sample <- glmulti(PCA2ScoresTraitSample ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")


install.packages("leaps")
library(leaps)
leaps=regsubsets(PCA1ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, nbest=10)
plot(leaps, scale="adjr2")
plot(leaps, scale="bic")
summary(leaps)
library(car)
subsets(leaps, statistic="bic")


require(lattice)
require(ggplot2)
PCA1 <- Peru_Plot_Master.data$PCA1ScoresPlotTraits_a
SolarRad. <- Peru_Plot_Master.data$SolarRadiation.GJ.m.2.yr.1.
Precip. <- Peru_Plot_Master.data$Precipitation.mm.yr.1.
Elev. <- Peru_Plot_Master.data$Elevation.m.
Temp. <- Peru_Plot_Master.data$MeanAnnualAirTemperature.degC.
pairs(~ PCA1 + SolarRad. + Precip. + Temp. + Elev., pch = 21)

#df<-data.frame(PCA1,SolarRad.,Precip.,Temp.,Elev.)



# pairwise correlation PCA1a
# help("cor") shows that pearson is default
my_data <- Peru_Plot_Master.data[, c(41,6,18,19,20,21)]
#head(my_data, 6)
res <- cor(na.omit(my_data))
round(res, 2)

## correlation matrix of abiotic correlations with PCA1 showing that temperature is best pairwise predictor

library(corrplot)
res <- cor(na.omit(my_data))
corrplot(res, method="circle")
corrplot(res, method="ellipse")
corrplot(res, method="number", type="upper", insig = "blank")
#corrplot(res, insig = "blank")
#corrplot.mixed(res, lower="ellipse", upper="number")


# pairwise correlation PCA2a
my_data2 <- Peru_Plot_Master.data[, c(42,6,18,19,20,21)]
#head(my_data, 6)
res2 <- cor(na.omit(my_data2))
round(res2, 2)

## correlation matrix of abiotic correlations with PCA1 showing that temperature is best pairwise predictor

library(corrplot)
res2 <- cor(na.omit(my_data2))
corrplot(res2, method="circle")
corrplot(res2, method="ellipse")
corrplot(res2, method="number", type="upper", insig = "blank")
#corrplot(res2, insig = "blank")
#corrplot.mixed(res, lower="ellipse", upper="number")



######### PCA subsampling ##############
# pairwise correlation PCA1 subsample
my_data <- Peru_Plot_Master.data[, c(43,6,18,19,20,21)]
#head(my_data, 6)
res <- cor(na.omit(my_data))
round(res, 2)

## correlation matrix of abiotic correlations with PCA1 showing that temperature is best pairwise predictor

library(corrplot)
res <- cor(na.omit(my_data))
corrplot(res, method="circle")
corrplot(res, method="ellipse")
corrplot(res, method="number", type="upper", insig = "blank")
#corrplot(res, insig = "blank")
#corrplot.mixed(res, lower="ellipse", upper="number")

# pairwise correlation PCA2 subsample
my_data2 <- Peru_Plot_Master.data[, c(44,6,18,19,20,21)]
#head(my_data, 6)
res2 <- cor(na.omit(my_data2))
round(res2, 2)

## correlation matrix of abiotic correlations with PCA1 showing that temperature is best pairwise predictor

library(corrplot)
res2 <- cor(na.omit(my_data2))
corrplot(res2, method="circle")
corrplot(res2, method="ellipse")
corrplot(res2, method="number", type="upper", insig = "blank")
#corrplot(res2, insig = "blank")
#corrplot.mixed(res, lower="ellipse", upper="number")

## precip comes out very strong on PCA2 subsample




#########################################################
########### Plots of PCA axes vs. environmental variation
myplot_PCA1vElev<- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PCA1ScoresPlotTraits_a)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA1vElev

myplot_PCA1vTemp<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PCA1ScoresPlotTraits_a)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA1vTemp


#myplot_PCA1vSolar<- ggplot(Peru_Plot_Master.data, aes(SolarRadiation.GJ.m.2.yr.1., PCA1ScoresPlotTraits_a)) +
  #geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
#myplot_PCA1vSolar

#myplot_PCA1vSoilMoist<- ggplot(Peru_Plot_Master.data, aes(Soil.moisture...., PCA1ScoresPlotTraits_a)) +
  #geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
#myplot_PCA1vSoilMoist


ModelPCA1a <- lm(PCA1ScoresPlotTraits_a ~ MeanAnnualAirTemperature.degC., data=Peru_Plot_Master.data)
summary(ModelPCA1a)  # y ~ x

ModelPCA1 <- lm(PCA1ScoresTraitSample ~ MeanAnnualAirTemperature.degC., data=Peru_Plot_Master.data)
summary(ModelPCA1)  # y ~ x


PCA2 <- myplot_PCA2vElev<- ggplot(Peru_Plot_Master.data, aes(SolarRadiation.GJ.m.2.yr.1., PCA2ScoresPlotTraits_a)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA2vElev

Model_PCA2a <- lm(PCA2ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1., data=Peru_Plot_Master.data)
summary(Model_PCA2a)

Model_PCA2 <- lm(PCA2ScoresPlotTraits_a ~ Soil.moisture...., data=Peru_Plot_Master.data)
summary(Model_PCA2)

Model_PCA2 <- lm(PCA2ScoresTraitSample~ Soil.moisture...., data=Peru_Plot_Master.data)
summary(Model_PCA2)

PCA2subsample <- myplot_PCA2vElevSubsample<- ggplot(Peru_Plot_Master.data, aes(Soil.moisture...., PCA2ScoresPlotTraits_a)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
myplot_PCA2vElevSubsample

png("Figure_PCA_Multi.png", units = "px", width=1500, height=900, res=170)

multiplot(PCA1, PCA2, cols=2)
dev.off()

#############################################
#Best predictors of PCA1 and PCA2 variation, climate?

#PCAfit <- glmulti(PCA1Scores ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + Elevation.m. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

PCAfit <- glmulti(PCA1ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1.+ MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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


PCA2fit <- glmulti(PCA2ScoresPlotTraits_a ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + MeanAnnualAirTemperature.degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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
###  Relationship between dominant community species traits and elevation
###


png("Figure_Plot_CMW_Traits.png", units = "px", width=1800, height=1800, res=200)
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

Model_T3 <- lm(P_N ~ Elevation.m., data=Peru_Plot_Master.data)
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
  ylab(bquote('Plot LMA ('~ kg^-1~m^2*')'))+
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
  ylab(bquote('Plot Aarea ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
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


multiplot(T1,T2,T3,T4,T5,T6,T7,T8,T11, cols=3)
dev.off()



T12 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SolarRadiation.GJ.m.2.yr.1.)) +
  xlab("Elevation (m)") + 
  ylab(bquote("SolarRadiation.GJ.m.2.yr.1.")) + 
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
T12

####################################################################
##################################
##################################
###### Distributions from subsampled trait distribuitons
##################################


########################################################
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
  expand_limits(y=c(-1.0,10.0)) +
  geom_hline(yintercept = 0)
#myplot_TDTSkew 
N3
ModelTDTNSkew <- lm(Elevation.m. ~ NSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTNSkew)
confint(ModelTDTNSkew)
coef(ModelTDTNSkew)

#Plot Kurtosis
N4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Kurtosis % N ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NKurtosisLower, ymax=NKurtosisUpper), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  
  #geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,20)) +
  geom_hline(yintercept = 0)

N4

ModelTDTNVar <- lm(Elevation.m. ~ NKurtosisMean, Peru_Plot_Master.data)
summary(NKurtosisMean)
confint(NKurtosisMean)
coef(NKurtosisMean)

# Multipanel plot
#multiplot(N1, N2, N3, N4, cols=2)
#dev.off()


################################################################
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
  expand_limits(y=c(-1.0,10.0)) +
  geom_hline(yintercept = 0)
#myplot_TDTCSkew 
C3        
ModelTDTCSkew <- lm(Elevation.m. ~ CSkewnessMean, Peru_Plot_Master.data)
summary(ModelTDTCSkew)
confint(ModelTDTCSkew)
coef(ModelTDTCSkew)


#Plot Kurtosis
C4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CKurtosisMean)) +
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
  #geom_smooth()+ 
  expand_limits(y=c(-10,20)) +
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


################################################################
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
  expand_limits(y=c(-1.0,10.0)) +
  geom_hline(yintercept = 0)
#myplot_TDTPSkew 

ModelTDTNVar <- lm(Elevation.m. ~ PSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
P4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PKurtosisMean)) +
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
  #geom_smooth()+ 
  expand_limits(y=c(-10,20)) +
  geom_hline(yintercept = 0)
#myplot_TDTPKurtosis 


ModelTDTNVar <- lm(Elevation.m. ~ PKurtosisMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
#multiplot(P1, P2, P3, P4, cols=2)
#dev.off()

################################################################
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
  expand_limits(y=c(-1.0,10.0)) +
  geom_hline(yintercept = 0)
#myplot_TDTCSkew 

ModelTDTPhotoVar <- lm(Elevation.m. ~ PhotoSkewnessMean, Peru_Plot_Master.data)

summary(ModelTDTPhotoVar)
confint(ModelTDTPhotoVar)
coef(ModelTDTPhotoVar)

#Plot Kurtosis
Photo4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PhotoKurtosisMean)) +
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
  geom_smooth(method=lm)+
  #geom_smooth()+ 
  expand_limits(y=c(-10,20)) +
  geom_hline(yintercept = 0)
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
                 expand_limits(y=c(-1.0,10.0)) +
                 geom_hline(yintercept = 0)
               #myplot_TDTSkew
               
ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)
               
#Plot Kurtosis
               
 LMA4 <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAKurtosisMean)) +
   xlab("Elevation (m)") + ylab("  Kurtosis LMA") +
   geom_point(size = 3, color="red") +
   geom_errorbar(aes(ymin=SLAKurtosisLower, ymax=SLAKurtosisUpper), width=.2,
                 position=position_dodge(0.05)) +
                 theme_bw() +
                 theme(axis.text = element_text(size = 9),
                       axis.text.y = element_text(size = rel(1.5), angle = 0)) +
                 theme(axis.text = element_text(size = 9),
                       axis.text.x = element_text(size = rel(1.5), angle = 0)) +
                 #geom_smooth(method=lm)+
                 #geom_smooth()+
                 expand_limits(y=c(-10,20)) +
                 geom_hline(yintercept = 0)
LMA4             
               
               ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)
               summary(ModelTDTNVar)
               confint(ModelTDTNVar)
               coef(ModelTDTNVar)
               
               # Multipanel plot
               #multiplot(LMA1, LMA2, LMA3, LMA4, cols=2)
               #dev.off()
               
png("Figure_Sampled_Traits_Combined.png", units = "px", width=2400, height=1800, res=170)
multiplot(N1, N2, N3, N4, C1, C2, C3, C4, P1, P2, P3, P4, Photo1, Photo2, Photo3, Photo4, LMA1, LMA2, LMA3, LMA4, cols=5)
dev.off() 



##############################################################################
# Subsampling method - derived plot-level traits
###########

# P:N 
Peru_Plot_Master.data$PlotNtoPMean <- ((Peru_Plot_Master.data$NMeanMean)/ (Peru_Plot_Master.data$PMeanMean))

PlotPtoNMean  <- (1/Peru_Plot_Master.data$PlotNtoPMean)

PN <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PlotPtoNMean)) +
  xlab("Elevation (m)") + ylab("Plot P:N") +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=SLAKurtosisLower, ymax=SLAKurtosisUpper), width=.2,
                #position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
  #geom_smooth()+
  #expand_limits(y=c(-10,20)) +
  #geom_hline(yintercept = 0)
PN  

ModelPN <- lm(Elevation.m. ~ PlotPtoNMean, Peru_Plot_Master.data)
summary(ModelPN)
confint(ModelPN)
coef(ModelPN)

# PNUE

Peru_Plot_Master.data$PNUE_Mean_subsample <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$NMeanMean))

PNUE <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PNUE_Mean_subsample)) +
  xlab("Elevation (m)") + ylab("Plot PNUE g × g leaf N^−1 × time^−1") +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=SLAKurtosisLower, ymax=SLAKurtosisUpper), width=.2,
  #position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()+
#expand_limits(y=c(-10,20)) +
#geom_hline(yintercept = 0)
PNUE 

ModelPNUE <- lm(Elevation.m. ~ PNUE_Mean_subsample, Peru_Plot_Master.data)
summary(ModelPNUE)
confint(ModelPNUE)
coef(ModelPNUE)


PNUEPhoto <- ggplot(Peru_Plot_Master.data, aes(PlotPtoNMean, PNUE_Mean_subsample)) +
  xlab("Plot P:N") + ylab("Plot PNUE g × g leaf N^−1 × time^−1") +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=SLAKurtosisLower, ymax=SLAKurtosisUpper), width=.2,
  #position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_smooth(method=lm)
#geom_smooth()+
#expand_limits(y=c(-10,20)) +
#geom_hline(yintercept = 0)
PNUEPhoto 

ModelPNUEPhoto  <- lm(elevation.m. ~ PNUE_Mean_subsample, Peru_Plot_Master.data)
summary(ModelPNUEPhoto)
confint(ModelPNUEPhoto)
coef(ModelPNUEPhoto)

png("Figure_PN_PNUE.png", units = "px", width=1500, height=900, res=200)
multiplot(PN, PNUE, cols=2)
dev.off() 



#=========================================================================================
# Draw Figure 6 with Celcius and 1/kT temperature axes
#=========================================================================================

png("Figure_MST_Plots.png", units = "px", width=2500, height=1000, res=300)

library(grid)

panel_6a.1 <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=MST_GPP1)) + 
  #ggtitle("a") + 
  geom_smooth(method = "lm") +
  geom_point(size=3, color="red") +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(paste("Inverse temperature, 1/kT (", eV^{-1}, ")"))) + 
  ylab(expression(paste("GPP / M " [tot]^{0.6}, "  (Mg C Mg ", M^{-0.6}, " ", ha^{-1}, yr^{-1}, ")"))) +
  theme_bw(base_size=12)
  #panel_6a.1 + geom_vline(xintercept = 39)

panel_6a.1 <- panel_6a.1 + geom_abline (intercept = 37, slope = -0.6, color="black")
 
#panel_6a.1 
  
grob <- grobTree(textGrob("A", x=0.9,  y=0.85, hjust=0, gp=gpar(col="Black", fontsize=15, fontface="bold"))) # Plot
panel_6A.1  <- panel_6a.1  + annotation_custom(grob)

#print(panel_6A.1)


panel_6b.1 <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=MST_NPP_new1)) +
  geom_smooth(method = "lm") +
  geom_point(size = 3, color="red") +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(paste("Inverse temperature, 1/kT (", eV^{-1}, ")"))) + 
  ylab(expression(paste("NPP / M " [tot]^{0.6}, "  (Mg C Mg ", M^{-0.6}, " ", ha^{-1}, yr^{-1}, ")"))) +
  theme_bw(base_size=12)

grob <- grobTree(textGrob("B", x=0.9,  y=0.85, hjust=0, gp=gpar(col="Black", fontsize=15, fontface="bold"))) # Plot
panel_6B.1  <- panel_6b.1  + annotation_custom(grob)

#print(panel_6B.1)

grid.arrange(arrangeGrob(panel_6A.1,panel_6B.1,ncol=2))
dev.off()


#########################
#########################
#########################
# predicting GPP - best model appears to be gpp ~ mass^0.6 x PNUE

install.packages("caret")
install.packages("AICcmodavg")
library(caret)
library(AICcmodavg)


# PhotosynthesisPerLeafN

m_1 <- lm( log(GPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)
confint(m_1)
coef(m_1)

AIC(m_1)
AICc(m_1)
modelEffectSizes(m_1)  #lm.sunSquares is depreciated
varImp(m_1, scale = FALSE)
avPlots(m_1)
crPlots(m_1)
confint(m_1)
vif(m_1) 

# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean
#m_1 <- lm( log(NPP_new) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
m_1 <- lm( log(NPP_new) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
summary(m_1)
confint(m_1)
coef(m_1)

AIC(m_1)
AICc(m_1)
modelEffectSizes(m_1)  #lm.sunSquares is depreciated
varImp(m_1, scale = FALSE)
avPlots(m_1)
crPlots(m_1)
confint(m_1)
vif(m_1) 


# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean
#m_1 <- lm( log(PhotosynthesisPerLeafN) ~  log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
m_1 <- lm( log(PhotosynthesisPerLeafN) ~   MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)  #PhotoPerLeafNMean estimate of E ~ 0.1539  PhotosynthesisPerLeafN ~ 0.1915
confint(m_1)
coef(m_1)

AIC(m_1)
AICc(m_1)
modelEffectSizes(m_1)  #lm.sunSquares is depreciated
varImp(m_1, scale = FALSE)
avPlots(m_1)
crPlots(m_1)
confint(m_1)
vif(m_1) 

# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean PlotNtoP
#m_1 <- lm( log(1/PlotNtoPMean) ~  log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
m_1 <- lm( log(1/PlotNtoPMean) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)  #PPlotNtoPMean estimate of E ~ 0.256  PhotosynthesisPerLeafN ~ 0.155
confint(m_1)
coef(m_1)

# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean PlotNtoP
#m_1 <- lm( log(1/PlotNtoPMean) ~  log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
m_1 <- lm( log(CMeanMean) ~  Elevation.m., data=Peru_Plot_Master.data)
summary(m_1)  #PPlotNtoPMean estimate of E ~ 0.256  PhotosynthesisPerLeafN ~ 0.155
confint(m_1)
coef(m_1)

CUE_C <- ((Peru_Plot_Master.data$CUE/ Peru_Plot_Master.data$CMeanMean))
#m_1 <- lm( log(1/PlotNtoPMean) ~  log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
m_1 <- lm(log(CUE_C) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)


m_1 <- lm( log(SLAMeanMean) ~  log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)

# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean SLAMeanMean  lma.sun.mu.abundance

#m_1 <- lm( log(SLAMeanMean) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)  #SLAMeanMean estimate of E ~ -0.1914  1/lma.sun.mu.abundance ~ -0.2923
confint(m_1)
coef(m_1)

AIC(m_1)
AICc(m_1)
modelEffectSizes(m_1)  #lm.sunSquares is depreciated
varImp(m_1, scale = FALSE)
avPlots(m_1)
crPlots(m_1)
confint(m_1)
vif(m_1) 


#### combine all traits in the model
rho_a <- (Peru_Plot_Master.data$PhotoPerLeafNMean * Peru_Plot_Master.data$SLAMeanMean * (1/(Peru_Plot_Master.data$PlotNtoPMean))

m_1 <- lm( log(rho) ~ MAinvBT, data=Peru_Plot_Master.data)
# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean SLAMeanMean  lma.sun.mu.abundance
#m_1 <- lm( log(SLAMeanMean) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)  # invariant with temperature!
confint(m_1)
coef(m_1)

rho_m <- (Peru_Plot_Master.data$Photosynthesis_massPerLeafN * (1/(Peru_Plot_Master.data$PlotNtoPMean)))
m_2 <- lm( log(rho_m) ~ MAinvBT, data=Peru_Plot_Master.data)
# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean SLAMeanMean  lma.sun.mu.abundance
#m_1 <- lm( log(SLAMeanMean) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_2)  # not significant
confint(m_1)
coef(m_1)


m_1 <- lm( log(NPP_new) ~ log(Aboveground_biomass2) + log(rho_m) + MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)

rho_a <- (Peru_Plot_Master.data$PhotoPerLeafNMean * (Peru_Plot_Master.data$SLAMeanMean))
m_3 <- lm( log(rho_a) ~ MAinvBT, data=Peru_Plot_Master.data)
# PhotosynthesisPerLeafN MAinvBT PhotosynthesisPerLeafN lma.sun.mu.abundance  PlotNtoPMean MAinvBT  PhotoPerLeafNMean SLAMeanMean  lma.sun.mu.abundance
#m_1 <- lm( log(SLAMeanMean) ~  MAinvBT, data=Peru_Plot_Master.data)
summary(m_3)  
confint(m_1)
coef(m_1)

m_1 <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
summary(m_1)
