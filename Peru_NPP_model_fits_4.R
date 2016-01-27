

########################################################################
#  Peru CHAMBASA preliminary NPP, Trait, climate plot analyses
#   Brian J. Enquist
#  12/20/15
########################################################################

# summary notes - overall simple pairwise correlation of NPP with
# above ground biomass appears to explain most of the variation
# in NPP. Additional multiple variable model selection approaches points to biomass as key driver of NPP and that Traits and climate expalin little variation in NPP and GPP. Climate including temp and solar radiation appear to also explain a lot of the variation but not as much as biomass. Biomass exponents include MST predicted value. Interestingly, plot level shifts in foliar N:P appear to covary with plot temperature in a way that is consistent with adaptive shifts in N:P -> increase foliar N-productivity

#Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged.csv",header=T)

#Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged2.csv",header=T)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged3.csv",header=T)

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


### Plotting Functions
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


####################################################################
#####  Exploratory Plots 
# Bivariate approach first - multivariate model competition below
########################

#################################################
##### Traits of dominant speies sampled within each plot ########
### assessing TDT predictions N, P, C, SLA, Photo -> mean_p_percent, mean_c_percent, mean_n_percent, mean_sla_lamina_petiole, mean_photosynthesis
### *** need 95%CIs ****


png("Figure_Plot_Traits.png")
#dev.off()
T1 <- myplot_TDTNAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_n_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf Nitrogen Plot") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
                #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
  #geom_smooth()
#myplot_TDTNAve


T2 <- myplot_TDTPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_p_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf Phosphorus") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTPAve

T3 <- myplot_TDTNtoPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PlotNtoP)) +
  xlab("Elevation (m)") + ylab("Plot Leaf N:P") +
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

T4 <- myplot_TDTCAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_c_percent)) +
  xlab("Elevation (m)") + ylab("Plot % Leaf C ") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
  #geom_smooth()
#myplot_TDTCAve

T5 <- myplot_TDTSLAAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_sla_lamina_petiole)) +
  xlab("Elevation (m)") + 
  ylab(bquote('Plot SLA ('~ m^2~kg^-1*')'))+
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTSLAAve

T6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., mean_photosynthesis)) +
  xlab("Elevation (m)") + 
  ylab(bquote('Plot Leaf Photosyn ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
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

T7 <- myplot_NtoPvPhoto <- ggplot(Peru_Plot_Master.data, aes(PlotNtoP, mean_photosynthesis)) +
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
#myplot_NtoPvPhoto
  #*note mean photosynthesis is not related to Plot N:P

T8 <- myplot_PhotovNeff <- ggplot(Peru_Plot_Master.data, aes(PhotosynthesisPerLeafN, mean_photosynthesis)) +

  xlab(bquote('Plot Leaf N efficiency ('*~CO[2]~ g^-1~s^-1*')')) +
  ylab(bquote('Plot Leaf Photosynthesis ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotovNeff
  #wow - highly positive . . variation in mean photosynthesis appears to be driven by Nefficiency

T9 <- myplot_NtoPvNeff <- ggplot(Peru_Plot_Master.data, aes(PlotNtoP, PhotosynthesisPerLeafN)) +
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
#myplot_NtoPvNeff
#huh . . . no correlation. Due to covariation?

multiplot(T1,T2,T3,T4,T5,T6, cols=2)
dev.off()



##################################
###### Distributions from subsampled trait distribuitons
##################################


########################################################
########################################################
# Nitrogen 
#
######

png("Figure_Sampled_Nitrogen.png")
#dev.off()
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
N1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., NMeanMean)) +
  xlab("Elevation (m)") + ylab("Sampled Mean % Leaf Nitrogen") +
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

#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean

ModelTDTNMean <- lm(Elevation..m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
N2 <- myplot_TDTNVar <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., NVarianceMean)) +
  xlab("Elevation (m)") + ylab("Sampled Variance % N ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
N3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., NSkewnessMean)) +
  xlab("Elevation (m)") + ylab("Sampled Skewness % N ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
N4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., NKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Sampled Kurtosis % N ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
multiplot(N1, N2, N3, N4, cols=2)
dev.off()



#####################################
### assessing TDT predictions  
# Nitrogen v Temperature

png("Nitrogen_MST")
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
N1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth()
#myplot_TDTNMean

#myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NMeanMean)) + geom_point(size = 3, color="red") + geom_errorbar(aes(ymin=NMeanLower, ymax=NMeanUpper), width=.2, position=position_dodge(0.05))
#myplot_TDTNMean


ModelTDTNMean <- lm(MAinvBT ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
N2 <- myplot_TDTNVar <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NVarianceMean)) +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=NVarianceLower, ymax=NVarianceUpper), width=.2,
                position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth()

#myplot_TDTNVar

ModelTDTNVar <- lm(MAinvBT ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
N3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NSkewnessMean)) +
  geom_point(size = 3, color="red")+
  geom_errorbar(aes(ymin=NSkewnessLower, ymax=NSkewnessUpper), width=.2,
                position=position_dodge(0.05)) +
  #geom_smooth(method=lm)
  geom_smooth() + 
  expand_limits(y=c(-2.5,4)) +
  geom_hline(yintercept = 0)
#myplot_TDTSkew

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
N4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, NKurtosisMean)) +
  geom_point(size = 3, color="red")+
  geom_errorbar(aes(ymin=NKurtosisLower, ymax=NKurtosisUpper), width=.2,
                position=position_dodge(0.05))+
  #geom_smooth(method=lm)
  geom_smooth() + 
  expand_limits(y=c(-10,30)) +
  geom_hline(yintercept = 0)

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

png("Figure_Sampled_Carbon.png")
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
C1 <- myplot_TDTCMean <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., CMeanMean)) +
  xlab("Elevation (m)") + ylab("Sampled Mean % Leaf C") +
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

ModelTDTNMean <- lm(Elevation..m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
C2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., CVarianceMean)) +
  xlab("Elevation (m)") + ylab("Sampled Variance % C ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
C3 <- myplot_TDTCSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., CSkewnessMean)) +
  xlab("Elevation (m)") + ylab("Sampled Skewness % C ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
C4 <- myplot_TDTCKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., CKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Sampled Kurtosis % C ") +
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


ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)
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

png("Figure_Sampled_Phosphorus.png")
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
P1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PMeanMean)) +
  xlab("Elevation (m)") + ylab("Sampled Mean % Leaf P") +
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

ModelTDTNMean <- lm(Elevation..m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
P2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PVarianceMean)) +
  xlab("Elevation (m)") + ylab("Sampled Variance % P ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
P3 <- myplot_TDTPSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PSkewnessMean)) +
  xlab("Elevation (m)") + ylab("Sampled Skewness % P ") +
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

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
P4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Sampled Kurtosis % P ") +
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


ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
multiplot(P1, P2, P3, P4, cols=2)
dev.off()

################
#### SLA
################

png("Figure_Sampled_SLA.png")
theme_set(theme_gray(base_size = 15))

SLA1 <- myplot_SLANMean <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., SLAMeanMean)) +
  xlab("Elevation (m)") +
  ylab(bquote('Subsampling Mean SLA ('~ m^2~kg^-1*')')) +
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
myplot_SLANMean
ModelTDTNMean <- lm(Elevation..m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)  
  
#Plot Variance
SLA2 <- myplot_TDTSLAVar <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., SLAVarianceMean)) +
  xlab("Elevation (m)") + ylab("Subsampling Variance SLA") +
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
myplot_TDTSLAVar


ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
SLA3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., SLASkewnessMean)) +
  xlab("Elevation (m)") + ylab("Sampled Skewness SLA ") +
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
myplot_TDTSkew

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
SLA4 <- myplot_TDTSLAKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation..m., SLAKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Sampled Kurtosis SLA ") +
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
#myplot_TDTKurtosis 

myplot_TDTSLAKurtosis

ModelTDTNVar <- lm(Elevation..m. ~ NVarianceMean, Peru_Plot_Master.data)
summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

# Multipanel plot
multiplot(SLA1, SLA2, SLA3, SLA4, cols=2)
dev.off()

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


png("SampledvsMeasuredTraitMeans.png")
p1 <- ggplot(Peru_Plot_Master.data, aes(mean_n_percent, TransformedNpercent)) + 
  xlab("% Nitrogen Subsampling") + ylab("% Nitrogen Plot") +
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


p2 <- ggplot(Peru_Plot_Master.data, aes(mean_p_percent, TransformedPpercent)) +
  xlab("% Phosphorus Subsampling") + ylab("% Phosphorus Plot") +
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


p3 <- ggplot(Peru_Plot_Master.data, aes(mean_c_percent, TransformedCpercent)) +
  xlab("% Carbon Subsampling") + ylab("% Carbon Plot") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=TransformedCpercentLower, ymax=TransformedCpercentUpper), width=.005,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)

p4 <- ggplot(Peru_Plot_Master.data, aes(mean_photosynthesis, PhotoMeanMean)) +
  
  xlab(bquote('Subsampling ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) + 
  ylab(bquote('Plot ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.005,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)


p5 <- ggplot(Peru_Plot_Master.data, aes(Transformedmean_sla_lamina_petiole, SLAMeanMean)) +
  xlab(bquote('Subsampling SLA ('~ m^2~kg^-1*')')) +
  ylab(bquote('Plot SLA ('~ m^2~kg^-1*')')) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=SLAMeanLower, ymax=SLAMeanUpper), width=.2,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)


multiplot(p1, p2, p3, p4, p5, cols=2)
dev.off()


####################################################
#############  shift in derived traits N:P and N per Photo
#N:P and N productivity across elevation  using subsampled trait values
###################################################

png("Figure_Stoich_MST.png")
theme_set(theme_gray(base_size = 10))

a_1 <- myplot_TDTNP<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1)+
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTNP

ModelTDTNP <- lm(Elevation..m. ~ PlotNtoPMean, Peru_Plot_Master.data)
summary(ModelTDTNP)
confint(ModelTDTNP)
coef(ModelTDTNP)

a_2 <- myplot_PhotoPerN<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotoPerLeafNMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NperNPP

ModelPhotoPerN <- lm(Elevation..m. ~ PhotoPerLeafNMean, Peru_Plot_Master.data)
summary(ModelPhotoPerN)
confint(ModelPhotoPerN)

coef(ModelTDTNP)

multiplot(a1, a2, cols=2)


a_3 <- myplot_TvPhotoPerN<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotoPerLeafNMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TvPhotoPerN

a_4 <- myplot_TvPhotoMeanMean<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotoMeanMean)) +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoMeanLower, ymax=PhotoMeanUpper), width=.2,
  position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvNtoP


multiplot(a_1, a_2, a_3, a_4, cols=2)
dev.off()
  #* foliar photosynthesis is largely independent of elevation and temperature. if anything slightly higher at highest elevations


## N:P . . . plotting both values calculated for the abundant species sampled on plot and the species subsampled. So, PlotNtoP is the foliar N:P from plot samples and PlotNtoPMean comes from local and BIEN traits for the rare species. 

png("Figure_StoichMST_2.png")
b1 <- myplot_NPvElev<- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvElev

b_1 <- myplot_NPvElevAve<- ggplot(Peru_Plot_Master.data, aes(Elevation..m., PlotNtoP)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvElevAve

b2 <- myplot_NPvTemp<- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvTemp

b_2 <- myplot_NPvTempAve<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoP)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvTempAve

b3 <- myplot_NPvBTemp<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvBTemp

ModelNPvBTemp <- lm(log(PlotNtoPMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelNPvBTemp)
confint(ModelNPvBTemp)
coef(ModelNPvBTemp)

  #* Next, using the average values of the most abundant species
b_3 <- myplot_NPvBTempAve<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_photosynthesis)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvBTempAve

ModelNPvBTempAve <- lm(log(PlotNtoP) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelNPvBTempAve)
confint(ModelNPvBTempAve)
coef(ModelNPvBTempAve)
  #pattern is much more clear and in accordance with MST predictions suggesting that there are discrepencies between inferences using subsampled traits and triats measured on the ground

b4 <- myplot_PhotoPerLeafNvBTemp<- ggplot(Peru_Plot_Master.data, aes(MAinvBT,PhotoPerLeafNMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotoPerLeafNvBTemp

ModelPhotoPerLeafNMeanvBTemp <- lm(log(PhotoPerLeafNMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafNMeanvBTemp)
confint(ModelPhotoPerLeafNMeanvBTemp)
coef(ModelPhotoPerLeafNMeanvBTemp)


b_4 <- myplot_PhotoPerLeafNvBTempAve<- ggplot(Peru_Plot_Master.data, aes(MAinvBT,PhotosynthesisPerLeafN)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotoPerLeafNvBTempAve

ModelPhotoPerLeafNMeanvBTempAve <- lm(log(PhotoPerLeafNMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafNMeanvBTempAve)
confint(ModelPhotoPerLeafNMeanvBTempAve)
coef(ModelPhotoPerLeafNMeanvBTempAve)

## Calculate photosynthetic efficiency per P

Peru_Plot_Master.data$PhotoPerPMeanMean <- ((Peru_Plot_Master.data$PhotoMeanMean)/(Peru_Plot_Master.data$PMeanMean))

Peru_Plot_Master.data$PhotoPerPAve <- ((Peru_Plot_Master.data$mean_photosynthesis)/(Peru_Plot_Master.data$mean_p_percent))

b5 <- myplot_PhotoPerLeafPvBTemp<- ggplot(Peru_Plot_Master.data, aes(MAinvBT,PhotoPerPMeanMean )) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotoPerLeafPvBTemp

ModelPhotoPerLeafNMeanvBTemp <- lm(log(PhotoPerLeafPMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafPMeanvBTemp)
confint(ModelPhotoPerLeafPMeanvBTemp)
coef(ModelPhotoPerLeafPMeanvBTemp)

b_5 <- myplot_PhotoPerLeafPvBTempAve<- ggplot(Peru_Plot_Master.data, aes(MAinvBT,PhotoPerPAve )) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PhotoPerLeafPvBTempAve

ModelPhotoPerLeafNMeanvBTempAve <- lm(log(PhotoPerPAve) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafPMeanvBTempAve)
confint(ModelPhotoPerLeafPMeanvBTempAve)
coef(ModelPhotoPerLeafPMeanvBTempAve)

  #* we do see support that increases in P relative to N can increase the photosynthesis per foliar N so that higher elevation sites have 

b6 <- myplot_PlotSLAvBTemp<- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC.,SLAMeanMean )) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PlotSLAvBTemp

ModelPhotoPerLeafNMeanvBTemp <- lm(log(PhotoPerLeafPMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafPMeanvBTemp)
confint(ModelPhotoPerLeafPMeanvBTemp)
coef(ModelPhotoPerLeafPMeanvBTemp)

b_6 <- myplot_PlotSLAvBTempAve<- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC.,mean_sla_lamina )) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PlotSLAvBTempAve

ModelPhotoPerLeafNMeanvBTemp <- lm(log(PhotoPerLeafPMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafPMeanvBTemp)
confint(ModelPhotoPerLeafPMeanvBTemp)
coef(ModelPhotoPerLeafPMeanvBTemp)

multiplot(b_1,b_2,b_3,b_4,b_5,b_6, cols=2)
dev.off()

##############################
#  Ecosystem measures NPP and GPP
###

png("Figure_Ecosystem_MST.png")
theme_set(theme_gray(base_size = 10))
NPP1 <- myplot_NPP<- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass2,NPP)) +
  geom_point(size = 3, color="red") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm) +
  annotation_logticks() 
myplot_NPP

GPP1 <- myplot_GPP<- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass,GPP)) +
  geom_point(size = 3, color="red") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm) +
  annotation_logticks()
myplot_GPP


#NPPperN <- myplot_NPPperN <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., NPPperNMeanMean)) + 
  #geom_point(size = 3, color="red") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth(method=lm)
myplot_NPPperN


#GPPperN <- myplot_GPPperN <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., GPPperNMeanMean)) + 
  #geom_point(size = 3, color="red") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth(method=lm)
#myplot_GPPperN

  # No strong relationship with NPP per foliar N and temperature . . . remember there is covariation with elevation and total biomass so these results likely driven by a shift in total biomass which could overly influence this expected increase in N productivity per total biomass

MST_NPP1 <- myplot_MST_NPP1<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, MST_NPP1)) +
  geom_point(size = 3, color="red") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm)
#myplot_MST_NPP1
ModelMST_NPP1 <- lm(MAinvBT ~ log(MST_NPP1), Peru_Plot_Master.data)
summary(ModelMST_NPP1)
confint(ModelMST_NPP1)
coef(ModelMST_NPP1)


MST_GPP1 <- myplot_MST_GPP1 <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, MST_GPP1)) + 
  geom_point(size = 3, color="red") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm)
#myplot_MST_GPP1
ModelMST_GPP1 <- lm(MAinvBT ~ log(MST_GPP1), Peru_Plot_Master.data)
summary(ModelMST_GPP1)
confint(ModelMST_GPP1)
coef(ModelMST_GPP1)
#hist(Peru_Plot_Master.data$MSTNPP)

multiplot(NPP1, GPP1, cols=2)
dev.off()

png("Figure_Ecosystem_MST_2.png")
multiplot(MST_GPP1, MST_NPP1, cols=2)
dev.off()








################################################################
################################################################
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

# model fits - linear model
ModelGPP <- lm(log10(GPP) ~ log10(Aboveground_biomass), Peru_Plot_Master.data)

summary(ModelGPP)
confint(ModelGPP)
coef(ModelGPP)


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
#myplot_sitetempP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., mean_p_percent)) + geom_point(size = 3) 
myplot_sitetempP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PMeanMean)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempP

        #* no strong relationship. looks flat. What about respiration rates? and N:P ratio?


###### site temperature v mean plot N / mean plot P or plot N:P
myplot_sitetempPlotNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoPMean)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP 

myplot_sitetempPlotNtoP  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoPMean)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempPlotNtoP 
#MST model fit
ModelNtoPvMAinvBT <- lm(log(PlotNtoP) ~ MAinvBT, Peru_Plot_Master.data)
ModelNtoPvMAinvBT <- lm(log(PlotNtoPMean) ~ MAinvBT, Peru_Plot_Master.data)

summary(ModelNtoPvMAinvBT)
AICc(ModelNtoPvMAinvBT)
confint(ModelNtoPvMAinvBT)
coef(ModelNtoPvMAinvBT)

  #*Boltzman fit is around -0.3 indicating that N:P

###### site temperature v RLeaf
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., RLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    #Leaf respiration increases with temperature. Looks like photosynthesis is modified by not leaf respiration.  If correct then the net carbon gain per leaf per carbon respiration changes with temperature. Could this be because of more biomass in warmer plots? 


###### site temperature v PhotosynthesisPerRLeaf
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    # cool - mean plot Photosynthesis per RLeaf decreases with increasing site temperature. But why should this be?  Change in the nitrogen use efficiency?  change in biomass?

###### site temperature v NPPLeafperRLeaf
myplot_sitetempNPPLeafperRLeaf <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., NPPLeafperRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 
theme_bw()

# log-log plot without log tick marks
myplot_sitetempNPPLeafperRLeaf 

  #* no strong relationship between MAT and NPPLeaf per RLeaf suggesting that changes in RLeaf with elevation is due to just more Leaf biomass in warm areas?

###### plot biomass v NPPLeafperRLeaf
myplot_Aboveground_biomassRLeaf <- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass, NPPLeafperRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 
theme_bw()

# log-log plot without log tick marks
myplot_Aboveground_biomassRLeaf 

  #* No, aparantly . . . no relationship between Aboveground biomass and NPPLeafPerR Leaf



###### site temperature v PhotosynthesisPerLeafN
myplot_sitetempLeafNEffic <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerLeafN)) + geom_point(size = 3) 

myplot_sitetempLeafNEffic <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PhotosynthesisPerLeafNMean)) + geom_point(size = 3) 

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

    ## * cool - strong decrease in photosynthesis per leaf N with mean annual air temp. This indicates that colder sites have higher N-productivity

###### plot temperature v leafN:P

myplot_sitetempNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempNtoP <- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PlotNtoPMean)) + geom_point(size = 3) 

myplot_sitetempNtoP

###### plot biomass v leafN:P
        #Boltzmann plot
myplot_siteBiomassNtoP  <- ggplot(Peru_Plot_Master.data, aes(x = Aboveground_biomass, y = PlotNtoP)) + geom_point(size = 3) 
theme_bw()

myplot_siteBiomassNtoP

    #* N:P is not related to above ground biomass. Constant?

### Plot mean P vs. plot mean N productivity - testing Kerkhoff et al. 2005
myplot_sitePandNProductivity  <- ggplot(Peru_Plot_Master.data, aes(x = mean_p_percent, y = PhotosynthesisPerLeafN)) + geom_point(size = 3) 


myplot_sitePandNProductivity

myplot_sitePandNProductivity   <- ggplot(Peru_Plot_Master.data, aes(x = PMeanMean, y = PhotosynthesisPerLeafNMean)) + geom_point(size = 3) 

myplot_sitePandNProductivity

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
#myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_n_percent, y = NPP))
myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = NMeanMean, y = NPP))
summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

##### mean plot leaf N  v elevation
#myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_n_percent, y = NMeanMean))

#myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = mean_n_percent))

myplot_NPP <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation..m., y = NMeanMean))

summary(myplot_NPP)

myplot_NPP_nice <- myplot_NPP + geom_point(size = 3)
myplot_NPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean RLeaf  v MAinvBT
myplot_RLeaf <- ggplot(data=Peru_Plot_Master.data, aes(x = MAinvBT, y = RLeaf))
summary(myplot_NPP)

myplot_RLeaf_nice <- myplot_RLeaf + geom_point(size = 3)
myplot_RLeaf_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


  #unlike photosynthesis, foliar respiration does show a temperature dependency.

#### PhotosynthesisPerLeafN v log10(mean_sla_lamina_petiole)
myplot_NProductionvsSLA <- ggplot(Peru_Plot_Master.data, aes(PhotosynthesisPerLeafNMean, mean_sla_lamina_petiole)) + geom_point(size = 3)
  
  myplot_NProductionvsSLA

  # model fits - linear model
ModelNProductionvsSLA <- lm(mean_sla_lamina_petiole ~ PhotosynthesisPerLeafNMean, Peru_Plot_Master.data)

summary(ModelNProductionvsSLA)
confint(ModelNProductionvsSLA)
coef(NProductionvsSLA)



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
m1NPP_N <- lm(log10(NPP)~ log10(NMeanMean), data=Peru_Plot_Master.data)
summary(m1NPP_N) 
AIC(m1NPP_N)
modelEffectSizes(m1NPP_N)  #lm.sunSquares is depreciated
confint(m1NPP_N)

## NPP climate via Boltzman temp,
m1NPP_temp <- lm(log(NPP)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1NPP_temp) 
AIC(m1NPP_temp)
AICc(m1NPP_temp)
modelEffectSizes(m1NPP_temp)  #lm.sunSquares is depreciated
confint(m1NPP_temp)

## NPP climate via preciptiation
m1NPP_precip <- lm(log10(NPP)~ Precipitation..mm.yr.1., data=Peru_Plot_Master.data)
summary(m1NPP_precip) 
AIC(m1NPP_precip)
AICc(m1NPP_precip)
modelEffectSizes(m1NPP_precip)  #lm.sunSquares is depreciated
confint(m1NPP_precip)

### GPP
m1GPP <- lm(log10(GPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m1GPP) 
AIC(m1GPP)
AICc(m1GPP)
modelEffectSizes(m1GPP)  #lm.sunSquares is depreciated
confint(m1GPP)

## temperature effect - Boltzmann temperature. Note, NPP is ln(NPP)
m1Temp <- lm(log(NPP)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1Temp) 
AIC(m1Temp)
AICc(m1Temp)
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
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

## multiple regression with just log10 Biomass
m4 <- lm(log10(GPP)~ log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 



## multiple regression with log10 biomass and Boltzmann temperature. 
m4 <- lm(log10(GPP)~ MAinvBT + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 
  ## the fitted biomass scaling exponent is 0.59. Temperature is not important. Could argue that this is the general model to fit so as to extract out the allometric exponent

## multiple regression with log10 biomass, Boltzmann temperature, and N-productivity 
m4 <- lm(log10(GPP)~ MAinvBT + PhotosynthesisPerLeafNMean + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 
  #* Impressive model fits but vif factors are too high!


## multiple regression with log10 biomass and N-productivity 
m4 <- lm(log10(GPP)~ PhotosynthesisPerLeafNMean + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 
#* Impressive model fits but vif factors are too high!


## multiple regression with log10 biomass, Boltzmann temperature, and foliar N:P 
m4 <- lm(log10(GPP)~ MAinvBT + PlotNtoPMean + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)
  #* model not significant plus vif values are higher

## multiple regression with kerkhoff et al. 2005 model log10 biomass and Boltzmann temperature and plot N:P. 
m4 <- lm(log10(GPP)~ PlotNtoP + log10(Aboveground_biomass), data=Peru_Plot_Master.data)
summary(m4) 
AICc(m4)
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 

  #* results suggest that biomass has the most effect on GPP variation. Also covariation between temperature and foliar N:P may cancel out leaving just importance of biomass

#### 
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
m7 <- lm(log(PlotNtoPMean) ~ MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AICc(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 



  #* if PlotNtoP covaries with temperature to compensate for kinetic effects of temp then we would expect that plot N:P scales as ln(N:P)~ 1/kT^0.6 or 0.33. This appears to be the case but confidence intervals are wide. Cold plots have low N:P (more P relative to N)

## linear model on MST prediction with MAinvBT and PlotNtoP
m7 <- lm(log(PhotosynthesisPerLeafNMean) ~ MAinvBT, data=Peru_Plot_Master.data)
#m7 <- lm(log(PhotosynthesisPerLeafN) ~ MAinvBT, data=Peru_Plot_Master.data)

summary(m7) 
AICc(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

    ### a similar result for N use efficiency of photosynthesis . .. But NUE increaes positively with MAinvBT. Cold plots have higher NUE

## is N:P ~ NUE ? kerkhoff et al. 2005 argues that N:P is NUE
m7 <- lm(PhotosynthesisPerLeafNMean ~ PlotNtoPMean, data=Peru_Plot_Master.data)
summary(m7) 
AICc(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

    #* ah, plot N:P is not related to NUE BUT! Kerkhoff et al predicts that changes in P then drives changes in NUE. True for chambasa?

#m7 <- lm(log10(PhotosynthesisPerLeafN) ~ log10(mean_p_percent), data=Peru_Plot_Master.data)

m7 <- lm(log10(PhotosynthesisPerLeafNMean) ~ log10(PMeanMean), data=Peru_Plot_Master.data)

summary(m7) 
AICc(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7)

    ## the positive trend reported by Kerkhoff et al. 2005 is there but it is not significant. 


#### Leaf respiration - RLeaf  v MAinvBT

m7 <- lm(log(RLeaf) ~ MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
AICc(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7)



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

# Using traits fromt the most abundant indivduals sampled within plot - no trait averaging for unknown species

fit2 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + mean_c_percent + mean_p_percent + PhotosynthesisPerLeafN + mean_photosynthesis + var_photosynthesis, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit2)
tmp <- weightable(fit2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit2@objects[[1]])
plot(fit2)
#Variable Importance
plot(fit2, type="s")


### Just biomass and traits of most abundant species
fit2_2 <- glmulti(log10(GPP) ~ log10(Aboveground_biomass) + mean_sla_lamina_petiole + var_sla_lamina_petiole + mean_n_percent + mean_c_percent + mean_p_percent + PhotosynthesisPerLeafN + mean_photosynthesis + var_photosynthesis, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit2_2)
tmp <- weightable(fit2_2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit2_2@objects[[1]])
plot(fit2_2)
#Variable Importance
plot(fit2_2, type="s")



fit3 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + NMeanMean + log10(Vegetation.height..m.), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

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
#fit4 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + Precipitation..mm.yr.1. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + mean_n_percent + log10(Vegetation.height..m.) + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafN + PlotNtoP, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")


## Models predicting GPP from environemntal and plot trait
fit4 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + log10(Aboveground_biomass) + log10(Vegetation.height..m.) + Precipitation..mm.yr.1. + MAinvBT + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

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

# Predicting GPP but excluding biomass
## Models predicting GPP from environemntal and plot trait
fit4 <- glmulti(log10(GPP) ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1. + MAinvBT + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

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

# Predicting GPP - just biomass and traits from subsampling
fit5 <- glmulti(log10(GPP) ~ log10(Aboveground_biomass) + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

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


#### Predicting total biomass
fit5 <- glmulti(log10(Aboveground_biomass) ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1. + Elevation..m. + SLAMeanMean + SLAMeanVariance + MAinvBT + NMeanMean + PMeanMean + PhotoMeanMean + PhotoVariance + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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
fit6 <- glmulti(MAinvBT ~ Aboveground_biomass + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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


############################################################################  
# PCA analyses
#######################################

#chambasapca <- princomp(~mean_sla_lamina_petiole + NMeanMean + PMeanMean + CMeanMean + mean_photosynthesis + PhotoPerLeafNMean + PlotNtoPMean, data=Peru_Plot_Master.data, cor=TRUE)

# Traits from most abundant individuals PCA 1
Photo <- Peru_Plot_Master.data$mean_photosynthesis
NP <- Peru_Plot_Master.data$PlotNtoP
P <- Peru_Plot_Master.data$mean_p_percent
N <- Peru_Plot_Master.data$mean_n_percent
C <- Peru_Plot_Master.data$mean_c_percent
SLA <- Peru_Plot_Master.data$mean_sla_lamina_petiole
Neff <- Peru_Plot_Master.data$PhotosynthesisPerLeafN

## Subsampling total plot distribution method
Photo <- Peru_Plot_Master.data$PhotoMeanMean
NP <- Peru_Plot_Master.data$PlotNtoPMean
P <- Peru_Plot_Master.data$PMeanMean
N <- Peru_Plot_Master.data$NMeanMean
C <- Peru_Plot_Master.data$CMeanMean
SLA <- Peru_Plot_Master.data$SLAMeanMean
Neff <- Peru_Plot_Master.data$PhotoPerLeafNMean

png("Figure_PCA_Plot_gradient.png")
chambasapca <- princomp(~ SLA + N + C + P + Photo + NP + Neff, data=Peru_Plot_Master.data, cor=TRUE)

summary(chambasapca)
loadings(chambasapca)
biplot(chambasapca, col=c("gray","red"),cex=c(0.4,0.5))
screeplot(chambasapca)
chambasapca$scores
chambasapca$loadings
  #* first two principle components explain about 72% of the variation. 
dev.off()


first.component.scores <- chambasapca$scores[,1]
summary(first.component.scores)
length(first.component.scores)
write.csv(first.component.scores, "first.component.scores.csv")

second.component.scores <- chambasapca$scores[,2]
summary(second.component.scores)
length(second.component.scores)
write.csv(second.component.scores, "second.component.scores.csv")


nrow(Peru_Plot_Master.data)
length(Peru_Plot_Master.data$Elevation..m.)
elev <- (Peru_Plot_Master.data$Elevation..m.)
length(elev)

#merge(elev,first.component.scores) 

png("Figure_PCA_gradient_Plot.png")
PCA1 <- myplot_PCA1vElev<- ggplot(Peru_Plot_Master.data, aes(Mean.annual.air.temperature..degC., PCA1ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PCA1vElev


PCA2 <- myplot_PCA2vElev<- ggplot(Peru_Plot_Master.data, aes(Precipitation..mm.yr.1., PCA2ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PCA2vElev


multiplot(PCA1, PCA2, cols=2)
dev.off()
#####
#Best predictors of PCA1 and PCA2 variation, climate?

#PCAfit <- glmulti(PCA1Scores ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1. + Elevation..m. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

PCAfit <- glmulti(PCA1ScoresPlotTraits ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1.+ Mean.annual.air.temperature..degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit)
tmp <- weightable(PCAfit)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit@objects[[1]])
plot(PCAfit)
#Variable Importance
plot(PCAfit, type="s")
#PCA1 is mainly temperature


PCA2fit <- glmulti(PCA2ScoresPlotTraits ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1. + Mean.annual.air.temperature..degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCA2fit)
tmp <- weightable(PCA2fit)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCA2fit@objects[[1]])
plot(PCA2fit)
#Variable Importance
plot(PCA2fit, type="s")


PCAfitSample <- glmulti(PCA1ScoresTraitSample ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1.+ Mean.annual.air.temperature..degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfitSample)
tmp <- weightable(PCAfitSample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfitSample@objects[[1]])
plot(PCAfitSample)
#Variable Importance
plot(PCAfitSample, type="s")
#PCA1 is mainly solar radiation

PCAfit2Sample <- glmulti(PCA2ScoresTraitSample ~ Solar.radiation..GJ.m.2.yr.1. + Precipitation..mm.yr.1.+ Mean.annual.air.temperature..degC. + mean_air_temp, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(PCAfit2Sample)
tmp <- weightable(PCAfit2Sample)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(PCAfit2Sample@objects[[1]])
plot(PCAfit2Sample)
#Variable Importance
plot(PCAfit2Sample, type="s")
#PCA1 is mainly solar radiation



library(corrplot)
library(qgraph)
qgraph(cor(Peru_Plot_Master))


qg.pca <- qgraph.pca(Peru_Plot_Master.data, factors = 2, rotation = "varimax")
