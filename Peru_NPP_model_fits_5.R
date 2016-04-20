

########################################################################
#  Peru CHAMBASA preliminary NPP_new, Trait, climate plot analyses
#   Brian J. Enquist
#  12/20/15
########################################################################

# summary notes - overall simple pairwise correlation of NPP_new with
# above ground biomass appears to explain most of the variation
# in NPP_new. Additional multiple variable model selection approaches points to biomass as key driver of NPP_new and that Traits and climate expalin little variation in NPP_new and GPP. Climate including temp and solar radiation appear to also explain a lot of the variation but not as much as biomass. Biomass exponents include MST predicted value. Interestingly, plot level shifts in foliar N:P appear to covary with plot temperature in a way that is consistent with adaptive shifts in N:P -> increase foliar N-productivity

#Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_new_Merged.csv",header=T)

#Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_new_Merged2.csv",header=T)

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged5.csv",header=T)

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

Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP)/(Peru_Plot_Master.data$MST_AGB))

Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))


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
##  MST predictions ln (NPP/M_Tot^0.6) ~ 1/kT ln(b)
##
###################################################################
Peru_Plot_Master.data$MST_AGB1 <- ((Peru_Plot_Master.data$Aboveground_biomass2)^0.6)
Peru_Plot_Master.data$MST_GPP1 <- ((Peru_Plot_Master.data$GPP_new_estimate)/(Peru_Plot_Master.data$MST_AGB))
Peru_Plot_Master.data$MST_NPP_new1 <- ((Peru_Plot_Master.data$NPP_new)/(Peru_Plot_Master.data$MST_AGB))

png("Figure_Plot_MST_Plots.png", units = "px", width=1800, height=1800, res=300)

MST1 <- ggplot(Peru_Plot_Master.data, aes(x=MAinvBT, y=log(MST_GPP1))) +
  xlab(italic("1/kT") + ylab(expression(ln(NPP/M[tot]^{0.6}))) +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
geom_smooth(method=lm)
MST1a <- MST1 + coord_fixed(ratio = 1.3)
MST1a

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
  xlab(italic("1/kT") + ylab(expression(ln(NPP/M[tot]^{0.6}))) +
  geom_point(size = 3, color="red") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  #geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
MST2a <- MST2 + coord_fixed(ratio = 1.1)
MST2a

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

multiplot(MST1a, MST2a, cols=2)
dev.off()

####################################################################
#####  Exploratory Plots 
# Bivariate approach first - multivariate model competition below
########################

#################################################
##### Traits of dominant speies sampled within each plot ########
### assessing TDT predictions N, P, C, SLA, Photo -> mean_p_percent, mean_c_percent, mean_n_percent, mean_sla_lamina_petiole, mean_photosynthesis
### *** need 95%CIs ****


png("Figure_Plot_Traits.png", units = "px", width=1800, height=1800, res=300)
#dev.off()
T1 <- myplot_TDTNAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_n_percent)) +
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
T1

T2 <- myplot_TDTPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_p_percent)) +
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
T2

T3 <- myplot_TDTNtoPAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PlotNtoP)) +
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

T4 <- myplot_TDTCAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_c_percent)) +
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

T5 <- myplot_TDTSLAAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_sla_lamina_petiole)) +
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

T6 <- myplot_TDTPhotoAve <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_photosynthesis)) +
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

T8 <- myplot_PhotovPNUE <- ggplot(Peru_Plot_Master.data, aes(PhotosynthesisPerLeafN, mean_photosynthesis)) +

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
#myplot_PhotovPNUE
  #wow - highly positive . . variation in mean photosynthesis appears to be driven by PNUE

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
#huh . . . no correlation using the subsampling method. Due to covariation?  signs of weak positive forrelation between P and PNUE using the traits from the most abundant species

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

png("Figure_Sampled_Nitrogen2.png", units="in", width=5, height=4, pointsize=12, res=900)
#dev.off()
theme_set(theme_gray(base_size = 15))
#theme_set(theme_classic(base_size = 30))
N1 <- myplot_TDTNMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NMeanMean)) +
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

ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
N2 <- myplot_TDTNVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NVarianceMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
N3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NSkewnessMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
N4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., NKurtosisMean)) +
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
  ylab("Subsampled % Leaf Nitrogen") +
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
  ylab("Subsampled Variance % Leaf Nitrogen") +
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
  ylab("Subsampled Skewness % Leaf Nitrogen") +
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
  ylab("Subsampled Kurtosis % Leaf Nitrogen") +
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

ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
C2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CVarianceMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
C3 <- myplot_TDTCSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CSkewnessMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
C4 <- myplot_TDTCKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., CKurtosisMean)) +
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

ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)


#Plot Variance
P2 <- myplot_TDTPVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PVarianceMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
P3 <- myplot_TDTPSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PSkewnessMean)) +
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

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
P4 <- myplot_TDTKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PKurtosisMean)) +
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
  xlab("Elevation (m)") + ylab("Sampled Mean % Leaf C") +
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
  xlab("Elevation (m)") + ylab("Sampled Variance % C ") +
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
  xlab("Elevation (m)") + ylab("Sampled Skewness % C ") +
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
  xlab("Elevation (m)") + ylab("Sampled Kurtosis % C ") +
  geom_point(size = 3, color="red") +
  geom_errorbar(aes(ymin=PhotoKurtosisLower, ymax=PhotoKurtosisUpper), width=.2,
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

png("Figure_Sampled_SLA2.png", units="in", width=10,height=8, pointsize=12,res=500)

theme_set(theme_gray(base_size = 15))

SLA1 <- myplot_SLANMean <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAMeanMean)) +
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
#myplot_SLANMean
ModelTDTNMean <- lm(Elevation.m. ~ NMeanMean, Peru_Plot_Master.data)

summary(ModelTDTNMean)
confint(ModelTDTNMean)
coef(ModelTDTNMean)  
  
#Plot Variance
SLA2 <- myplot_TDTSLAVar <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAVarianceMean)) +
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
#myplot_TDTSLAVar


ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Skewness
SLA3 <- myplot_TDTSkew <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLASkewnessMean)) +
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
#myplot_TDTSkew

ModelTDTNVar <- lm(Elevation.m. ~ NVarianceMean, Peru_Plot_Master.data)

summary(ModelTDTNVar)
confint(ModelTDTNVar)
coef(ModelTDTNVar)

#Plot Kurtosis
SLA4 <- myplot_TDTSLAKurtosis <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., SLAKurtosisMean)) +
  xlab("Elevation (m)") + ylab("Sampled Kurtosis SLA") +
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


png("SampledvsMeasuredTraitMeans2.png", units="in", width=10,height=8, pointsize=12,res=500)

p1 <- ggplot(Peru_Plot_Master.data, aes(mean_n_percent, TransformedNpercent)) + 
  ylab("% Nitrogen Subsampling") + xlab("% Nitrogen Plot") +
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

Modelp1 <- lm(mean_n_percent ~ TransformedNpercent, Peru_Plot_Master.data)
summary(Modelp1)
confint(Modelp1)
coef(Modelp1)

p2 <- ggplot(Peru_Plot_Master.data, aes(mean_p_percent, TransformedPpercent)) +
  ylab("% Phosphorus Subsampling") + xlab("% Phosphorus Plot") +
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
  ylab("% Carbon Subsampling") + xlab("% Carbon Plot") +
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
  xlab(bquote('Plot ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
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
  xlab(bquote('Plot SLA ('~ m^2~kg^-1*')')) + 
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
  ylab("% N:P Subsampling") + xlab("% N:P Plot") +
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

# relationships between plot variance of abundant species traits with subsampled values. Really the only fair comparison is between the mean values
p6 <- ggplot(Peru_Plot_Master.data, aes(var_photosynthesis, PhotoVarianceMean)) + 
  xlab("Variance Photosynthesis Subsampling") + ylab("Variance Photosynthesis Plot") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=TransformedNpercentLower, ymax=TransformedNpercentUpper), width=.01,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)+
  expand_limits(y=c(0,35)) +
  geom_hline(yintercept = 0)
#p6

p7 <- ggplot(Peru_Plot_Master.data, aes(var_p_percent, PVarianceMean)) + 
  xlab("Variance % Phosphorus Subsampling") + ylab("Variance % Phosphorus Plot") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=PVarianceLower, ymax=PVarianceUpper), width=.0001,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)
  #expand_limits(y=c(0,35)) +
  #geom_hline(yintercept = 0)
#p7


p8 <- ggplot(Peru_Plot_Master.data, aes(var_n_percent, NVarianceMean)) + 
  xlab("Variance % Nitrogen Subsampling") + ylab("Variance % Nitrogen Plot") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=NVarianceLower, ymax=NVarianceUpper), width=.0001,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)
  #expand_limits(y=c(0,0.001)) +
  #geom_hline(yintercept = 0)
#p8


p9 <- ggplot(Peru_Plot_Master.data, aes(var_sla_lamina_petiole, SLAVarianceMean)) + 
  xlab("Variance SLA Subsampling") + ylab("Variance SLA Plot") +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = rel(1.5), angle = 0)) +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(size = rel(1.5), angle = 0)) +
  geom_errorbar(aes(ymin=SLAVarianceLower, ymax=SLAVarianceUpper), width=.0001,
                position=position_dodge(0.05)) +
  geom_smooth(method=lm)+ expand_limits(y=c(0,10)) +
  geom_hline(yintercept = 0)
#p9


####################################################
#############  shift in derived traits N:P and N per Photo
#N:P and N productivity across elevation  using subsampled trait values
###################################################

png("Figure_Stoich_MST2.png", units="in", width=5, height=4, pointsize=12, res=900) 

theme_set(theme_gray(base_size = 10))

a_1 <- myplot_TDTNP<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  geom_abline(intercept = 0, slope = 1)+
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_TDTNP

ModelTDTNP <- lm(Elevation.m. ~ PlotNtoPMean, Peru_Plot_Master.data)
summary(ModelTDTNP)
confint(ModelTDTNP)
coef(ModelTDTNP)

a_2 <- myplot_PhotoPerN<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotoPerLeafNMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NperNPP_new

ModelPhotoPerN <- lm(Elevation.m. ~ PhotoPerLeafNMean, Peru_Plot_Master.data)
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

png("Figure_StoichMST_3.png", units="in", width=5, height=4, pointsize=12, res=900))

b1 <- myplot_NPvElev<- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PlotNtoPMean)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvElev

b_1 <- myplot_NPvElevAve<- ggplot(Peru_Plot_Master.data, aes(Elevation.m., PlotNtoP)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_NPvElevAve

b2 <- myplot_NPvTemp<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PlotNtoPMean)) +
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

b6 <- myplot_PlotSLAvBTemp<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC.,SLAMeanMean )) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PlotSLAvBTemp

ModelPhotoPerLeafNMeanvBTemp <- lm(log(PhotoPerLeafPMean) ~ MAinvBT, Peru_Plot_Master.data)
summary(ModelPhotoPerLeafPMeanvBTemp)
confint(ModelPhotoPerLeafPMeanvBTemp)
coef(ModelPhotoPerLeafPMeanvBTemp)

b_6 <- myplot_PlotSLAvBTempAve<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC.,mean_sla_lamina )) +
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

##########################################################################################
#  Ecosystem measures NPP_new and GPP
###


MSTGPP <- ((Peru_Plot_Master.data$GPP_new_estimate)*(1/Peru_Plot_Master.data$Aboveground_biomass2))

MSTNPP <- ((Peru_Plot_Master.data$NPP_new)*(1/Peru_Plot_Master.data$Aboveground_biomass2))
  
png("Figure_Ecosystem_MST.png", units="in", width=5, height=4, pointsize=12, res=900)

theme_set(theme_gray(base_size = 10))
NPP_new1 <- myplot_NPP_new<- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass2,NPP_new)) +
  geom_point(size = 3, color="red") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm) +
  annotation_logticks() 
myplot_NPP_new

ModelNPP_new1 <- lm(log10(NPP_new) ~ log10(Aboveground_biomass2) , Peru_Plot_Master.data)
summary(ModelNPP_new1)
confint(ModelNPP_new1)
coef(ModelNPP_new1)


GPP1 <- myplot_GPP<- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass2, GPP_new_estimate)) +
  geom_point(size = 3, color="red") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm) +
  annotation_logticks()
myplot_GPP

ModelGPP_new1 <- lm(log10(GPP_new_estimate) ~ log10(Aboveground_biomass2) , Peru_Plot_Master.data)
summary(ModelGPP_new1)
confint(ModelGPP_new1)
coef(ModelGPP_new1)

#NPP_newperN <- myplot_NPP_newperN <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., NPP_newperNMeanMean)) + 
  #geom_point(size = 3, color="red") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth(method=lm)
#myplot_NPP_newperN


#GPPperN <- myplot_GPPperN <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., GPPperNMeanMean)) + 
  #geom_point(size = 3, color="red") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth(method=lm)
#myplot_GPPperN

  # No strong relationship with NPP_new per foliar N and temperature . . . remember there is covariation with elevation and total biomass so these results likely driven by a shift in total biomass which could overly influence this expected increase in N productivity per total biomass

MST_NPP_new1 <- myplot_MST_NPP_new1<- ggplot(Peru_Plot_Master.data, aes(MAinvBT, MSTNPP)) +
  geom_point(size = 3, color="red") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm)
#myplot_MST_NPP_new1

ModelMST_NPP_new1 <- lm(MAinvBT ~ log(MSTNPP), Peru_Plot_Master.data)
summary(ModelMST_NPP_new1)
confint(ModelMST_NPP_new1)
coef(ModelMST_NPP_new1)


MST_GPP1 <- myplot_MST_GPP1 <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, MSTGPP)) + 
  geom_point(size = 3, color="red") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method=lm)
#myplot_MST_GPP1

ModelMST_GPP1 <- lm(MAinvBT ~ log(MSTGPP), Peru_Plot_Master.data)
summary(ModelMST_GPP1)
confint(ModelMST_GPP1)
coef(ModelMST_GPP1)
#hist(Peru_Plot_Master.data$MSTNPP_new)

multiplot(NPP_new1, GPP1, cols=2)
dev.off()

png("Figure_Ecosystem_MST_2.png", units="in", width=5, height=4, pointsize=12, res=900)

multiplot(MST_GPP1, MST_NPP_new1, cols=2)
dev.off()








################################################################
################################################################
#GPP v Biomass

myplot_GPP <- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass2, GPP)) + geom_point(size = 3) +
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
ModelGPP <- lm(log10(GPP_new_estimate) ~ log10(Aboveground_biomass2), Peru_Plot_Master.data)

summary(ModelGPP)
confint(ModelGPP)
coef(ModelGPP)


########
#ResidenceTime v Biomass

myplot_Residence_time <- ggplot(data=Peru_Plot_Master.data, aes(x = Aboveground_biomass2, y = Residence_time))
summary(myplot_Residence_time)

myplot_Residence_time_nice <- myplot_Residence_time + geom_point(size = 3)
myplot_Residence_time_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

########
#NPP_new v Solar Radiation
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = SolarRadiation.GJ.m.2.yr.1., y = GPP_new_estimate))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 6)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


########
#elevation v leaf temperature

myplot_leaftemp <- ggplot(Peru_Plot_Master.data, aes(Elevation.m., mean_air_temp)) + geom_point(size = 3) 

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
  
myplot_plottempphoto <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., mean_photosynthesis)) + geom_point(size = 3) 
  
myplot_plottempphoto <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_photosynthesis)) + geom_point(size = 3) 
  
  
theme_bw()
  
# log-log plot without log tick marks
myplot_plottempphoto
  
  # ** there is an interesting positive correlation with site temperature and the plot mean photosynthesis indicating that colder plots have similar if not higher levels of mean photosynthesis. Hypothesis - shift in traits help enable similar rates of metabolism (photosynthesis)
  
######## 
### site temperature v SLA
  
  myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., mean_sla_lamina_petiole)) + geom_point(size = 3) 
  
  theme_bw()
  
  # log-log plot without log tick marks
  myplot_sitetempsla
  #*weak positive correlation

###### site temperature v SLA mean_sla_lamina
  
myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., mean_sla_lamina)) + geom_point(size = 3) 
  
#myplot_sitetempsla <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_sla_lamina)) + geom_point(size = 3) 

theme_bw()
  
# log-log plot without log tick marks
myplot_sitetempsla
  
          #*weak positive correlation
  
###### site temperature v mean_n_percent
myplot_sitetempN <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., mean_n_percent)) + geom_point(size = 3) 

# myplot_sitetempN  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, mean_n_percent)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempN

      #*very weak positive correlation


###### site temperature v mean_p_percent
#myplot_sitetempP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., mean_p_percent)) + geom_point(size = 3) 
myplot_sitetempP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PMeanMean)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempP

        #* no strong relationship. looks flat. What about respiration rates? and N:P ratio?


###### site temperature v mean plot N / mean plot P or plot N:P
myplot_sitetempPlotNtoP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempPlotNtoP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PlotNtoPMean)) + geom_point(size = 3) 

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
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., RLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    #Leaf respiration increases with temperature. Looks like photosynthesis is modified by not leaf respiration.  If correct then the net carbon gain per leaf per carbon respiration changes with temperature. Could this be because of more biomass in warmer plots? 


###### site temperature v PhotosynthesisPerRLeaf
myplot_sitetempRLeaf <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PhotosynthesisPerRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempRLeaf 

    # cool - mean plot Photosynthesis per RLeaf decreases with increasing site temperature. But why should this be?  Change in the nitrogen use efficiency?  change in biomass?

###### site temperature v NPP_LeafperRLeaf
myplot_sitetempNPP_LeafperRLeaf <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., NPP_LeafperRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 
theme_bw()

# log-log plot without log tick marks
myplot_sitetempNPP_LeafperRLeaf 

  #* no strong relationship between MAT and NPP_newLeaf per RLeaf suggesting that changes in RLeaf with elevation is due to just more Leaf biomass in warm areas?

###### plot biomass v NPP_LeafperRLeaf
myplot_Aboveground_biomass2RLeaf <- ggplot(Peru_Plot_Master.data, aes(Aboveground_biomass2, NPP_LeafperRLeaf)) + geom_point(size = 3) 

#myplot_sitetempRLeaf  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, RLeaf)) + geom_point(size = 3) 
theme_bw()

# log-log plot without log tick marks
myplot_Aboveground_biomass2RLeaf 

  #* No, aparantly . . . no relationship between Aboveground biomass and NPP_newLeafPerR Leaf



###### site temperature v PhotosynthesisPerLeafN
myplot_sitetempLeafPNUE <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PhotosynthesisPerLeafN)) + geom_point(size = 3) 

myplot_sitetempLeafPNUE <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PhotosynthesisPerLeafNMean)) + geom_point(size = 3) 

#myplot_sitetempLeafPNUEic  <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotosynthesisPerLeafN)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempLeafPNUE 

    ## ** impressive negative correlation - warmer temps have lower mean plot photosynthesis per unit leaf nitrogen. 


###### photosynthesis v PhotosynthesisPerLeafN
myplot_sitetempLeafPNUE <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PhotosynthesisPerLeafN)) + geom_point(size = 3) 

#myplot_sitetempLeafPNUE <- ggplot(Peru_Plot_Master.data, aes(MAinvBT, PhotosynthesisPerLeafN)) + geom_point(size = 3) 

theme_bw()

# log-log plot without log tick marks
myplot_sitetempLeafPNUE 

    ## * cool - strong decrease in photosynthesis per leaf N with mean annual air temp. This indicates that colder sites have higher N-productivity

###### plot temperature v leafN:P

myplot_sitetempNtoP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PlotNtoP)) + geom_point(size = 3) 

myplot_sitetempNtoP <- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PlotNtoPMean)) + geom_point(size = 3) 

myplot_sitetempNtoP

###### plot biomass v leafN:P
        #Boltzmann plot
myplot_siteBiomassNtoP  <- ggplot(Peru_Plot_Master.data, aes(x = Aboveground_biomass2, y = PlotNtoP)) + geom_point(size = 3) 
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
#NPP_new v mean plot SLA
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = log10(mean_sla_lamina_petiole), y = NPP_new))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


#GPP v mean plot SLA
myplot_GPP <- ggplot(data=Peru_Plot_Master.data, aes(x = log10(mean_sla_lamina_petiole), y = GPP_new_estimate))
summary(myplot_GPP)

myplot_GPP_nice <- myplot_GPP + geom_point(size = 3)
myplot_GPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


#mean plot SLA v elevation
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation.m., y = mean_sla_lamina_petiole))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#varplot SLA v elevation
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation.m., y = var_sla_lamina_petiole))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#NPP_new v mean plot SLA lamina
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina, y = NPP_new))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean plot SLA lamina v elevation
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation.m., y = mean_sla_lamina))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean plot SLA lamina v mean plot SLA w petiole
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina_petiole, y = mean_sla_lamina))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#NPP_new v mean plot N
#myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_n_percent, y = NPP_new))
myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = NMeanMean, y = NPP_new))
summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

##### mean plot leaf N  v elevation
#myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_n_percent, y = NMeanMean))

#myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation.m., y = mean_n_percent))

myplot_NPP_new <- ggplot(data=Peru_Plot_Master.data, aes(x = Elevation.m., y = NMeanMean))

summary(myplot_NPP_new)

myplot_NPP_new_nice <- myplot_NPP_new + geom_point(size = 3)
myplot_NPP_new_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#mean RLeaf  v MAinvBT
myplot_RLeaf <- ggplot(data=Peru_Plot_Master.data, aes(x = MAinvBT, y = RLeaf))
summary(myplot_NPP_new)

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

myplot_NPP_new_nice + geom_point(aes(color = Elevation.m.), size =6) # explore effect of elevation
myplot_NPP_new_nice + geom_point(aes(color = SolarRadiation.GJ.m.2.yr.1.), size =6) # explore effect of solar radiation
myplot_NPP_new_nice + geom_point(aes(color = Soil.type), size =6) # explore effect of solar radiation
myplot_NPP_new_nice + geom_point(aes(color = Aboveground_biomass2), size =6) # explore effect of above ground biomass
myplot_NPP_new_nice + geom_point(aes(color = Slope..deg.), size =6) # explore effect of above ground biomass

myplot_NPP_new_nice + geom_point(aes(color = Ptotal..mg.kg.1.), size =6) # explore effect of above ground biomass

# GPP
myplot_NPP_new_nice + geom_point(aes(color = mean_sla_lamina_petiole), size =6) # explore effect of mean SLA
myplot_GPP <- ggplot(data=Peru_Plot_Master.data, aes(x = mean_sla_lamina_petiole, y = GPP_new_estimate))
summary(myplot_GPP)
myplot_GPP_nice <- myplot_GPP + geom_point(size = 4)

myplot_GPP_nice + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


###################################################
#####  Fitting simple linear models 
##################################################

# NPP_new and aboveground biomass
m1NPP_new <- lm(log10(NPP_new)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data) 
summary(m1NPP_new) 
AIC(m1NPP_new)
modelEffectSizes(m1NPP_new)  #lm.sunSquares is depreciated
confint(m1NPP_new)


# NPP_new and aboveground biomass
m1NPP_new <- lm(log10(NPP_new)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data) 
summary(m1NPP_new) 
AIC(m1NPP_new)
modelEffectSizes(m1NPP_new)  #lm.sunSquares is depreciated
confint(m1NPP_new)

# GPP and aboveground biomass
m1GPP <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data) 
summary(m1GPP) 
AIC(m1GPP)
modelEffectSizes(m1GPP)  #lm.sunSquares is depreciated
confint(m1GPP)

# NPP and solar radiation
m1NPP <- lm(log10(GPP_new_estimate)~ log10(SolarRadiation.GJ.m.2.yr.1.), data=Peru_Plot_Master.data) 
summary(m1NPP_new) 
AIC(m1NPP_new)
modelEffectSizes(m1NPP_new)  #lm.sunSquares is depreciated
confint(m1NPP_new)

# NPP_new and mean plot SLA
m1NPP_new_SLA <- lm(log10(NPP_new)~ log10(var_sla_lamina_petiole), data=Peru_Plot_Master.data)
summary(m1NPP_new_SLA) 
AIC(m1NPP_new_SLA)
modelEffectSizes(m1NPP_new_SLA)  #lm.sunSquares is depreciated
confint(m1NPP_new_SLA)

# NPP_new and mean plot N
m1NPP_new_N <- lm(log10(NPP_new)~ log10(NMeanMean), data=Peru_Plot_Master.data)
summary(m1NPP_new_N) 
AIC(m1NPP_new_N)
modelEffectSizes(m1NPP_new_N)  #lm.sunSquares is depreciated
confint(m1NPP_new_N)

## NPP_new climate via Boltzman temp,
m1NPP_new_temp <- lm(log(NPP_new)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1NPP_new_temp) 
AIC(m1NPP_new_temp)
AICc(m1NPP_new_temp)
modelEffectSizes(m1NPP_new_temp)  #lm.sunSquares is depreciated
confint(m1NPP_new_temp)

## NPP_new climate via preciptiation
m1NPP_new_precip <- lm(log10(NPP_new)~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
summary(m1NPP_new_precip) 
AIC(m1NPP_new_precip)
AICc(m1NPP_new_precip)
modelEffectSizes(m1NPP_new_precip)  #lm.sunSquares is depreciated
confint(m1NPP_new_precip)

### GPP
m1GPP <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m1GPP) 
AIC(m1GPP)
AICc(m1GPP)
modelEffectSizes(m1GPP)  #lm.sunSquares is depreciated
confint(m1GPP)

## temperature effect - Boltzmann temperature. Note, NPP_new is ln(NPP_new)
m1Temp <- lm(log(NPP_new)~ MAinvBT, data=Peru_Plot_Master.data)
summary(m1Temp) 
AIC(m1Temp)
AICc(m1Temp)
modelEffectSizes(m1Temp)  #lm.sunSquares is depreciated
confint(m1Temp )


m2 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m2) 
AIC(m2)
modelEffectSizes(m2)  #lm.sunSquares is depreciated
avPlots(m2)
crPlots(m2)
confint(m2)
vif(m2) 




######################################### 
######################################### 
######################################### 
## Multiple regression 
######################################### 

#NPP_new and soil moisture and temperature
m3 <- lm(log10(NPP_new)~ Soil.moisture....+ MAinvBT + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

### NPP_new and biomass and boltzman temp  **use
m3 <- lm(log(NPP_new)~ MAinvBT + log(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m3) 
lm.sumSquares(m3)
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

install.packages("caret")
library(caret)
### GPP and biomass and boltzman temp  **use
m3 <- lm(log(GPP_new_estimate)~ MAinvBT + log(Aboveground_biomass2) +, data=Peru_Plot_Master.data)
summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

### GPP and biomass and boltzman temp  and mean photosynthesis  **use
m3 <- lm(log(GPP_new_estimate)~ log(Aboveground_biomass2) + MAinvBT+ log(mean_photosynthesis), data=Peru_Plot_Master.data)

#m3 <- lm(log(GPP_new_estimate)~ log(Aboveground_biomass2) + MAinvBT+ log(PhotoMeanMean), data=Peru_Plot_Master.data)


summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 





## GPP and biomass and boltzman temp  **use
m3 <- lm(log10(GPP_new_estimate)~1 + log10(Aboveground_biomass2) + Precipitation.mm.yr.1.+ MAinvBT, data=Peru_Plot_Master.data)
summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

#########  
## multiple regression GPP with just log10 Biomass  **use
m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 

## multiple regression GPP with just log10 Biomass  **use
m4 <- lm(log10(NPP_new)~ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

#########
## Using abundant species trait measures
##  GPP vs. biomass and mean photosynthesis

m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(mean_photosynthesis), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

###############################################
## Using abundant species trait measures   **use
##  GPP vs. biomass and mean photosynthesis and variance of photosynthesis from subsampled distribution

m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(mean_photosynthesis) +log10(PhotoVarianceMean), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

##  GPP vs. biomass and variance of photosynthesis from subsampled distribution **use
m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PhotoVarianceMean), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

##  GPP vs. biomass and variance of photosynthesis from subsampled distribution and plot N:P from abundant species

m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PhotoVarianceMean) + PlotNtoP, data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

##  GPP vs. biomass and plot N:P from abundant species **use

m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PlotNtoP), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

##  GPP vs. biomass and plot N:P from abundant species
m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PlotNtoP), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

##  GPP vs. biomass and plot N:P from subsampled distribution
m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PlotNtoPMean), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

m4 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

## multiple regression with log10 biomass and Boltzmann temperature. 
m4 <- lm(log10(GPP_new_estimate)~ MAinvBT + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
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
m4 <- lm(log10(GPP_new_estimate)~ MAinvBT + PhotosynthesisPerLeafNMean + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m4) 
AICc(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4) 
  #* Impressive model fits but vif factors are too high!


## multiple regression with log10 biomass and N-productivity 
m_4 <- lm(log10(GPP_new_estimate) ~ PhotosynthesisPerLeafNMean + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m_4) 
AICc(m_4)
modelEffectSizes(m_4)  #lm.sunSquares is depreciated
avPlots(m_4)
crPlots(m_4)
confint(m_4)
vif(m_4) 
#* Impressive model fits but vif factors are too high!


## multiple regression with log10 biomass, Boltzmann temperature, and foliar N:P 
m4 <- lm(log10(GPP_new_estimate)~ MAinvBT + PlotNtoPMean + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
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
m4 <- lm(log10(GPP_new_estimate)~ PlotNtoP + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
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
m4 <- lm(log10(GPP_new_estimate)~ SolarRadiation.GJ.m.2.yr.1. + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(m4) 
AIC(m4)
modelEffectSizes(m4)  #lm.sunSquares is depreciated
avPlots(m4)
crPlots(m4)
confint(m4)
vif(m4)

## multiple regression with log10 height and Boltzmann temperature. 
m5 <- lm(log10(GPP_new_estimate)~ MAinvBT + log10(Vegetation.height..m.), data=Peru_Plot_Master.data)
summary(m5) 
AIC(m5)
modelEffectSizes(m5)  #lm.sunSquares is depreciated
avPlots(m5)
crPlots(m5)
confint(m5)
vif(m5) 

## multiple regression with log10 biomass and mean plot SLA. 
m6 <- lm(log10(GPP_new_estimate)~ SLAVarianceMean+ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
#m6 <- lm(log10(GPP_new_estimate)~ mean_sla_lamina_petiole + SLAVarianceMean+ log10(Aboveground_biomass2), data=Peru_Plot_Master.data)

summary(m6) 
AIC(m6)
modelEffectSizes(m6)  #lm.sunSquares is depreciated
avPlots(m6)
crPlots(m6)
confint(m6)
vif(m6) 

m7 <- lm(log10(GPP_new_estimate)~ 1 + mean_sla_lamina_petiole + MAinvBT, data=Peru_Plot_Master.data)
summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m6) 


m7 <- lm(log10(GPP_new_estimate)~ MAinvBT + log10(Aboveground_biomass2) + log10(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
m7 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)

summary(m7) 
AIC(m7)
modelEffectSizes(m7)  #lm.sunSquares is depreciated
avPlots(m7)
crPlots(m7)
confint(m7)
vif(m7) 

## multiple regression with log10 biomass and mean plot log10SLA. 
m7 <- lm(log10(NPP_new)~ log10(mean_sla_lamina_petiole) + log10(Aboveground_biomass2), data=Peru_Plot_Master.data)
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

pairs(~ Elevation.m. + SolarRadiation.GJ.m.2.yr.1.+ Precipitation.mm.yr.1.+ + Vegetation.height..m.+ MAinvBT+ log10(Aboveground_biomass2) + log10(NPP_new) + log10(GPP_new_estimate), data=Peru_Plot_Master.data)

pairs(~ Elevation.m. + SolarRadiation.GJ.m.2.yr.1.+ Precipitation.mm.yr.1.+ mean_sla_lamina + mean_sla_lamina_petiole+ mean_n_percent+ var_sla_lamina_petiole + MAinvBT+ log10(Aboveground_biomass2) + log10(NPP_new) + log10(GPP_new_estimate), data=Peru_Plot_Master.data)


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
#fit <- lm(GPP ~ SolarRadiation.GJ.m.2.yr.1. + MAinvBT + log10(Aboveground_biomass2) + mean_sla_lamina_petiole + var_sla_lamina_petiole, data=Peru_Plot_Master.data)

fit <- lm(NPP_new ~ SolarRadiation.GJ.m.2.yr.1. + log10(Aboveground_biomass2) + mean_sla_lamina_petiole + var_sla_lamina_petiole, data=Peru_Plot_Master.data)
step <- stepAIC(fit, direction="both")
step$anova # display results

## glmulti finds what are the n best models (the confidence set of models) among all possible models (the candidate set, as specified by the user). Models are fitted with the specified fitting function (default is glm) and are ranked with the specified Information Criterion (default is aicc). The best models are found either through exhaustive screening of the candidates, or using a genetic algorithm, which allows very large candidate sets to be adressed. The output can be used for model selection, variable selection, and multimodel inference.

install.packages("glmulti")
#http://core.ac.uk/download/files/153/6303128.pdf
library(glmulti)
library(leaps)
library(MASS)
library(lme4)

############## GPP ###################
fit2 <- glmulti(log(GPP_new_estimate) ~ log(SolarRadiation.GJ.m.2.yr.1.) + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + log(mean_sla_lamina_petiole) + log(SLAVarianceMean) + MAinvBT + log(mean_n_percent) + log(mean_c_percent) + log(mean_p_percent) + log(PlotNtoP) + (PhotosynthesisPerLeafN) + log(mean_photosynthesis) + log(mean_p_percent) + log(PhotoVarianceMean), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit2)
tmp <- weightable(fit2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(fit2@objects[[1]])
plot(fit2)
plot(fit2, type="s")

# here are the best models ranked by aicc
#1  log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent)
#2  log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2) + MAinvBT
#3  log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2)
#4  log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT 
#5  log(GPP_new_estimate) ~ 1 + MAinvBT + log(PhotoVarianceMean)


Model1 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) +  log(mean_n_percent), data=Peru_Plot_Master.data)

summary(Model1)
AIC(Model1)
AICc(Model1)
modelEffectSizes(Model1)  #lm.sunSquares is depreciated
varImp(Model1, scale = FALSE)
confint(Model1)
avPlots(Model1)
crPlots(Model1)
confint(Model1)
vif(Model1)

Model2 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
summary(Model2)
AIC(Model2)
AICc(Model2)
modelEffectSizes(Model2)  #lm.sunSquares is depreciated
varImp(Model2, scale = FALSE)
confint(Model2)
avPlots(Model2)
crPlots(Model2)
confint(Model2)
vif(Model2)

Model3 <- lm(log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2), data=Peru_Plot_Master.data)
summary(Model3)
AIC(Model3)
AICc(Model3)
modelEffectSizes(Model3)  #lm.sunSquares is depreciated
varImp(Model3, scale = FALSE)
confint(Model3)

Model4 <- lm(log(GPP_new_estimate) ~ 1 + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
summary(Model4)
AIC(Model4)
AICc(Model4)
modelEffectSizes(Model4)  #lm.sunSquares is depreciated
varImp(Model4, scale = FALSE)
confint(Model4)

############## NPP ###################
fit3 <- glmulti(log(NPP_new) ~ log(SolarRadiation.GJ.m.2.yr.1.) + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + log(mean_sla_lamina_petiole) + log(SLAVarianceMean) + MAinvBT + log(mean_n_percent) + log(mean_c_percent) + log(mean_p_percent) + log(PlotNtoP) + log(PhotosynthesisPerLeafN) + log(mean_photosynthesis) + log(mean_p_percent) + log(PhotoVarianceMean), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit3)
tmp <- weightable(fit3)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 5,]
tmp
summary(fit3@objects[[1]])
plot(fit3)
plot(fit3, type="s")

# the best models are 
#1 log(NPP_new) ~ 1 + log(mean_c_percent) + log(PhotosynthesisPerLeafN) + log(PhotoVarianceMean)
#2 log(NPP_new) ~ 1 + log(Aboveground_biomass2) + MAinvBT
#3 log(NPP_new) ~ 1 + log(PhotosynthesisPerLeafN) + log(PhotoVarianceMean)
#4 log(NPP_new) ~ 1 + MAinvBT + log(PhotoVarianceMean)
#5 log(NPP_new) ~ 1 + log(mean_sla_lamina_petiole) + MAinvBT

Model1NPP <- lm(log(NPP_new) ~ log(mean_c_percent) + log(PhotosynthesisPerLeafN) + log(PhotoVarianceMean), data=Peru_Plot_Master.data)

summary(Model1NPP)
AIC(Model1NPP)
AICc(Model1NPP)
modelEffectSizes(Model1NPP)  #lm.sunSquares is depreciated
varImp(Model1NPP, scale = FALSE)
confint(Model1NPP)
avPlots(Model1NPP)
crPlots(Model1NPP)
confint(Model1NPP)
vif(Model1NPP)

Model2NPP <- lm(log(NPP_new) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)

summary(Model2NPP)
AIC(Model2NPP)
AICc(Model2NPP)
modelEffectSizes(Model2NPP)  #lm.sunSquares is depreciated
varImp(Model2NPP, scale = FALSE)
confint(Model2NPP)
avPlots(Model2NPP)
crPlots(Model2NPP)
confint(Model2NPP)
vif(Model2NPP)


Model3NPP <- lm(log(NPP_new) ~ 1 + log(PhotosynthesisPerLeafN) + log(PhotoVarianceMean), data=Peru_Plot_Master.data)

summary(Model3NPP)
AIC(Model3NPP)
AICc(Model3NPP)
modelEffectSizes(Model3NPP)  #lm.sunSquares is depreciated
varImp(Model3NPP, scale = FALSE)
confint(Model3NPP)
avPlots(Model3NPP)
crPlots(Model3NPP)
confint(Model3NPP)
vif(Model3NPP)

############################

fit2 <- glmulti(log(NPP_new) ~ SolarRadiation.GJ.m.2.yr.1. + log(Aboveground_biomass2) + Precipitation.mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + mean_c_percent + mean_p_percent + PhotosynthesisPerLeafN + mean_photosynthesis + var_photosynthesis, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit2)
tmp <- weightable(fit2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(fit2@objects[[1]])
plot(fit2)
plot(fit2, type="s")

####################################################
### Just biomass and traits of most abundant species
fit2_2 <- glmulti(log10(GPP_new_estimate) ~ log10(Aboveground_biomass2) + mean_sla_lamina_petiole + var_sla_lamina_petiole + mean_n_percent + mean_c_percent + mean_p_percent + PhotosynthesisPerLeafN + mean_photosynthesis + var_photosynthesis, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit2_2)
tmp <- weightable(fit2_2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit2_2@objects[[1]])
plot(fit2_2)
#Variable Importance
plot(fit2_2, type="s")


fit3 <- glmulti(log(GPP_new_estimate) ~ SolarRadiation.GJ.m.2.yr.1. + log(Aboveground_biomass2) + Precipitation.mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + NMeanMean + log(Vegetation.height..m.), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit3)
tmp <- weightable(fit3)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit3@objects[[1]])
plot(fit3)

#Variable Importance
plot(fit3, type="s")

#### Now fit many more leaf traits including mean plot photosynthesis
#fit4 <- glmulti(log10(GPP_new_estimate) ~ SolarRadiation.GJ.m.2.yr.1. + log10(Aboveground_biomass2) + Precipitation.mm.yr.1. + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + MAinvBT + mean_n_percent + log10(Vegetation.height..m.) + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafN + PlotNtoP, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")


## Models predicting GPP from environemntal and plot trait
fit4 <- glmulti(log10(GPP_new_estimate) ~ SolarRadiation.GJ.m.2.yr.1. + log10(Aboveground_biomass2) + MeanAnnualAirTemperature.degC. + Precipitation.mm.yr.1. + MAinvBT + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit4)
tmp <- weightable(fit4)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit4@objects[[1]])
plot(fit4)
#Variable Importance
plot(fit4, type="s")

# Predicting GPP but excluding biomass
## Models predicting GPP from environemntal and plot trait
fit4 <- glmulti(log10(GPP_new_estimate) ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + MAinvBT + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit4)
tmp <- weightable(fit4)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit4@objects[[1]])
plot(fit4)
#Variable Importance
plot(fit4, type="s")

##################################################################
# Predicting GPP - just biomass and traits from subsampling
fit5 <- glmulti(log10(GPP_new_estimate) ~ log10(Aboveground_biomass2) + SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit5)
tmp <- weightable(fit5)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit5@objects[[1]])
plot(fit5)
#Variable Importance
plot(fit5, type="s")

####### but w no biomass, just plot traits from most abundant species
fit5 <- glmulti(log10(GPP_new_estimate) ~  SLAMeanMean + SLAVarianceMean + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + PhotoMeanMean + PhotoVarianceMean + CMeanMean + CVarianceMean + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")

summary(fit5)
tmp <- weightable(fit5)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7
#So, we could now examine the "best" model in closer detail with:
summary(fit5@objects[[1]])
plot(fit5)
#Variable Importance
plot(fit5, type="s")




#### Predicting total biomass
fit5 <- glmulti(log10(Aboveground_biomass2) ~ SolarRadiation.GJ.m.2.yr.1. + Precipitation.mm.yr.1. + Elevation.m. + SLAMeanMean + SLAMeanVariance + MAinvBT + NMeanMean + PMeanMean + PhotoMeanMean + PhotoVariance + PhotoPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit5)
tmp <- weightable(fit5)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7

#So, we could now examine the "best" model in closer detail with:
summary(fit5@objects[[1]])
plot(fit5)

#Variable Importance
plot(fit5, type="s")

    #* interesting temperature and elevation comes out on top but the variable model average importance of terms is less than 0.8. Does this suggest that predicting total biomass is difficult? Importance of site history (land slides etc.)?- just the chance that there is a big tree in the sampled plot?



#### Predicting environmental temperature from plot traits
fit6 <- glmulti(MAinvBT ~ Aboveground_biomass2 + log10(mean_sla_lamina_petiole) + var_sla_lamina_petiole + NMeanMean + NVarianceMean + PMeanMean + PVarianceMean + mean_photosynthesis + var_photosynthesis + PhotosynthesisPerLeafNMean + PlotNtoPMean, data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit6)
tmp <- weightable(fit6)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
#  note log10(GPP_new_estimate) ~ 1 + var_sla_lamina_petiole comes out as #7

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
NPP_new.tree <- tree(log10(GPP_new_estimate) ~ SolarRadiation.GJ.m.2.yr.1. + log10(Aboveground_biomass2) + Precipitation.mm.yr.1. + mean_sla_lamina_petiole + var_sla_lamina_petiole + MAinvBT + mean_n_percent + Vegetation.height..m., data = Peru_Plot_Master.data, mindev=0)

NPP_new.tree
plot(residuals(NPP_new.tree) ~ predict(NPP_new.tree))
plot(NPP_new.tree, type = "uniform")
text(NPP_new.tree, cex = 0.5, all = T)


### Pruning tree
plot(prune.tree(NPP_new.tree))
NPP_new.tree.prune <- prune.tree(NPP_new.tree, best = 4)
plot(NPP_new.tree.prune, type = "uniform")
text(NPP_new.tree.prune, cex = 0.5, all = T)


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

png("Figure_PCA_gradient_plot.png", units="in", width=5, height=4, pointsize=12, res=900)

biplot(chambasapca, col=c("gray","red"),cex=c(0.8,0.5))
dev.off()
png("Figure_PCA_plot_gradien_plot.png", units="in", width=5, height=4, pointsize=12, res=900))
screeplot(chambasapca)
dev.off()
chambasapca$scores
chambasapca$loadings
  #* first two principle components explain about 72% of the variation. 

first.component.scores <- chambasapca$scores[,1]
summary(first.component.scores)
length(first.component.scores)
write.csv(first.component.scores, "first.component.scores.csv")

second.component.scores <- chambasapca$scores[,2]
summary(second.component.scores)
length(second.component.scores)
write.csv(second.component.scores, "second.component.scores.csv")


nrow(Peru_Plot_Master.data)
length(Peru_Plot_Master.data$Elevation.m.)
elev <- (Peru_Plot_Master.data$Elevation.m.)
length(elev)

#merge(elev,first.component.scores) 

#png("Figure_PCA_gradient_Plot.png")

PCA1 <- myplot_PCA1vElev<- ggplot(Peru_Plot_Master.data, aes(MeanAnnualAirTemperature.degC., PCA1ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PCA1vElev


PCA2 <- myplot_PCA2vElev<- ggplot(Peru_Plot_Master.data, aes(Precipitation.mm.yr.1., PCA2ScoresPlotTraits)) +
  geom_point(size = 3, color="red") +
  #geom_errorbar(aes(ymin=CMeanLower, ymax=CMeanUpper), width=.2,
  #position=position_dodge(0.05)) +
  geom_smooth(method=lm)
#myplot_PCA2vElev


multiplot(PCA1, PCA2, cols=2)
dev.off()
#####
#Best predictors of PCA1 and PCA2 variation, climate?

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



