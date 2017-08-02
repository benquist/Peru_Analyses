
#==========================================================================
# Code set 3
#
# Assessing trait-based scaling theory in tropical forests spanning a broad temperature gradient 
#
#Authors:
#Brian J. Enquist1,2*, Lisa Patrick Bentley3, Alexander Shenkin3, Brian Maitner1, Van Savage4, Sean Michaletz1,5, Benjamin Blonder3, Vanessa Buzzard1, Tatiana Erika Boza Espinoza6, William Farfan-Rios 7,8 Chris Doughty9, Gregory R. Goldsmith10, Roberta E. Martin11, Norma Salinas3,8, Miles Silman8, Sandra Díaz12, Gregory P. Asner11 & Yadvinder Malhi3
#
# Multiple regression analyses and plots, model competition analyse
# 8/2/17
#==========================================================================

Peru_Plot_Master.data <- read.csv(file="/Users/brianjenquist/GitHub/R/Peru_Analyses/Peru_Gradient_NPP_Merged7.csv",header=T)

#==========================================================================

install.packages("glmulti")
#http://core.ac.uk/download/files/153/6303128.pdf
library(glmulti)
library(leaps)
library(MASS)
library(lme4)

#==========================================================================
############## GPP ###################

fit2 <- glmulti(log(GPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + log(Malhi_2017SoilCstock.MgCha.1.from0to30cm.) + log(Soil.moisture....) + log(1/mean_sla_lamina_petiole) + log(SLAVarianceMean) + MAinvBT + log(mean_n_percent) + log(mean_c_percent) + log(mean_p_percent) + log(PlotNtoP) + log(PhotosynthesisPerLeafN) + log(mean_photosynthesis) + log(mean_p_percent) + log(PhotoVarianceMean), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit2)
tmp <- weightable(fit2)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 20,]
tmp
summary(fit2@objects[[1]])
plot(fit2)
plot(fit2, type="s")

# here are the best models ranked by aicc model
#1  log(GPP_Malhi_2017) ~ 1 + MAinvBT + log(mean_c_percent) + log(mean_p_percent) + log(mean_photosynthesis)
#2  log(GPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2) + log(SLAVarianceMean) + log(mean_n_percent)
#3  log(GPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2) + log(mean_n_percent)
#4  log(GPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2) + MAinvBT
#5  log(GPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2)
#6  log(GPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent)

Model1GPP <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)

summary(Model1GPP)
AIC(Model1GPP)
AICc(Model1GPP)
modelEffectSizes(Model1GPP)  #lm.sunSquares is depreciated
varImp(Model1GPP, scale = FALSE)
confint(Model1GPP)
avPlots(Model1GPP)
crPlots(Model1GPP)
confint(Model1GPP)
vif(Model1GPP)

## Use avPlots to calculate adjusted values and store as lists
panel_4a.data <- avPlots(Model1GPP, ~ MAinvBT, col="gray70", col.lines="black", main="", xlab=expression(atop(paste("Adjusted average"), paste("annual temperature, <1/kT>" [gs], " (", eV^{-1}, ")"))), ylab=expression(paste("Adjusted ln(GPP" [yr], ")")))

panel_4b.data <- avPlots(Model1GPP, ~log(Aboveground_biomass2), col="gray70", col.lines="black", main="", xlab=expression("Adjusted ln(stand biomass)"), ylab=expression(paste("Adjusted ln(GPP" [yr], ")")))

## Store adjusted values in a dataframe
DataFig4 <- data.frame(panel_4a.data, panel_4b.data)
colnames(DataFig4) <- c("AdjMAinvBT","AdjlnGPP_MAinvBT", "AdjlnBiomass", "AdjlnGPP_Biomass")

## Convert ln'ed data back to absolute numbers for plotting
DataFig4$AdjBiomass <- exp(DataFig4$AdjlnBiomass)
DataFig4$AdjGPP_MAinvBT <- exp(DataFig4$AdjlnGPP_MAinvBT)
DataFig4$AdjGPP_Biomass <- exp(DataFig4$AdjlnGPP_Biomass)


## Plot panels
panel_4a <-ggplot(DataFig4, aes(x=AdjMAinvBT, y=AdjGPP_MAinvBT)) + 
  ggtitle("a") + 
  geom_point(size = 4, color="darkgrey") + 
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3),
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(atop(paste("Adjusted average annual"), paste("temperature, <1/kT>" [gs], " (", eV^{-1}, ")")))) + 
  ylab(expression(paste("Adjusted GPP" [yr], " (kg ", ha^{-1}, " ", yr^{-1}, ")"))) +
  #geom_smooth(method=lm) +
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
print(panel_4a)

panel_4b <-ggplot(DataFig4, aes(x=AdjBiomass, y=AdjGPP_Biomass)) + 
  ggtitle("b") + 
  geom_point(size = 4, color="darkgrey") + 
  scale_x_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=4), 
                     labels = trans_format("log", math_format(e^.x))) + 
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  xlab(paste("Adjusted stand biomass (kg)")) + 
  ylab(expression(paste("Adjusted GPP" [yr], " (kg ", ha^{-1}, " ", yr^{-1}, ")"))) +
  #geom_smooth(method=lm) +
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
print(panel_4b)


#png("Figure_Adjusted_GPP.png", units = "px", width=900, height=600, res=150)
#multiplot(panel_4a, panel_4b, cols=2)
#dev.off()

#==========================================================================
############## NPP ###################
fit3 <- glmulti(log(NPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + log(Malhi_2017SoilCstock.MgCha.1.from0to30cm.) + log(Soil.moisture....) + log(1/mean_sla_lamina_petiole) + log(SLAVarianceMean) + MAinvBT + log(mean_n_percent) + log(mean_c_percent) + log(mean_p_percent) + log(PlotNtoP) + log(PhotosynthesisPerLeafN) + log(mean_photosynthesis) + log(mean_p_percent) + log(PhotoVarianceMean), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
summary(fit3)
tmp <- weightable(fit3)
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 5,]
tmp
summary(fit3@objects[[1]])
plot(fit3)
plot(fit3, type="s")

# the best models are 
#1 log(NPP_Malhi_2017) ~ 1 + log(Aboveground_biomass2) + MAinvBT
#2 log(NPP_Malhi_2017) ~ 1 + MAinvBT + log(PhotoVarianceMean)
#3 log(NPP_Malhi_2017) ~ 1 + log(mean_sla_lamina_petiole) + MAinvBT
#4 log(NPP_Malhi_2017) ~ 1 + MAinvBT + log(mean_c_percent)
#5 log(NPP_Malhi_2017) ~ 1 + log(SLAVarianceMean) + MAinvBT
#6 log(NPP_Malhi_2017) ~ 1 + log(mean_photosynthesis) + log(PhotoVarianceMean)
#7 log(NPP_Malhi_2017) ~ 1 + MAinvBT + log(mean_photosynthesis)

Model1NPP <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)

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

## Use avPlots to calculate adjusted values and store as lists
panel_5a.data <- avPlots(Model1NPP, ~ MAinvBT, col="gray70", col.lines="black", main="", xlab=expression(atop(paste("Adjusted average"), paste("annual temperature, <1/kT>" [gs], " (", eV^{-1}, ")"))), ylab=expression(paste("Adjusted ln(NPP" [yr], ")")))

panel_5b.data <- avPlots(Model1NPP, ~log(Aboveground_biomass2), col="gray70", col.lines="black", main="", xlab=expression("Adjusted ln(stand biomass)"), ylab=expression(paste("Adjusted ln(NPP" [yr], ")")))

## Store adjusted values in a dataframe
DataFig5 <- data.frame(panel_5a.data, panel_5b.data)
colnames(DataFig5) <- c("AdjMAinvBT","AdjlnNPP_MAinvBT", "AdjlnBiomass", "AdjlnNPP_Biomass")

## Convert ln'ed data back to absolute numbers for plotting
DataFig5$AdjBiomass <- exp(DataFig5$AdjlnBiomass)
DataFig5$AdjNPP_MAinvBT <- exp(DataFig5$AdjlnNPP_MAinvBT)
DataFig5$AdjNPP_Biomass <- exp(DataFig5$AdjlnNPP_Biomass)

## Plot panels
panel_5a <-ggplot(DataFig5, aes(x=AdjMAinvBT, y=AdjNPP_MAinvBT)) + 
  ggtitle("c") + 
  geom_point(size = 4, color="darkgrey") + 
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3),
                     labels = trans_format("log", math_format(e^.x))) + 
  xlab(expression(atop(paste("Adjusted average annual"), paste("temperature, <1/kT>" [gs], " (", eV^{-1}, ")")))) + 
  ylab(expression(paste("Adjusted NPP" [yr], " (kg ", ha^{-1}, " ", yr^{-1}, ")"))) +
  #geom_smooth(method=lm) +
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
print(panel_5a)

panel_5b <-ggplot(DataFig5, aes(x=AdjBiomass, y=AdjNPP_Biomass)) + 
  ggtitle("d") + 
  geom_point(size = 4, color="darkgrey") +
  scale_x_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=4), 
                     labels = trans_format("log", math_format(e^.x))) + 
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  xlab(paste("Adjusted stand biomass (kg)")) + 
  ylab(expression(paste("Adjusted NPP" [yr], " (kg ", ha^{-1}, " ", yr^{-1}, ")"))) +
  #geom_smooth(method=lm) +
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  theme_bw(base_size=12) + 
  theme(legend.position="none", plot.title=element_text(hjust=0.94, vjust=-1.8))
print(panel_5b)


png("Figure_Adjusted_GPP_NPP.png", units = "px", width=900, height=600, res=100)
multiplot(panel_4a, panel_4b, panel_5a, panel_5b, cols=2)
dev.off()


#==========================================================================
##########################
install.packages("MuMIn")
library(MuMIn)
#m1 <- lm(log(GPP_Malhi_2017) ~  log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)

#m2 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
#m3 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
#m4 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
#m5 <- lm(log(GPP_Malhi_2017) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
#m6 <- lm(log(GPP_Malhi_2017) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
#m7 <- lm(log(GPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)

#==========================================================================
# Here are the main models hypothesized . . .
m1 <- lm(log(GPP_Malhi_2017) ~ MAinvBT, data=Peru_Plot_Master.data) #climate models
m2 <- lm(log(GPP_Malhi_2017) ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
m3 <- lm(log(GPP_Malhi_2017) ~ log(Precipitation.mm.yr.1.), data=Peru_Plot_Master.data)
m4 <- lm(log(GPP_Malhi_2017) ~ log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m5 <- lm(log(GPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.), data=Peru_Plot_Master.data)
m6 <- lm(log(GPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m7 <- lm(log(GPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m8 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data) #Biomass model
m9 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data) # Biomass and climate
m10 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m11 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data) # TDT Biomass and Traits
m12 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)
m13 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
m14 <- lm(log(GPP_Malhi_2017) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
m15 <- lm(log(GPP_Malhi_2017) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
m16 <- lm(log(GPP_Malhi_2017) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(NVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
m17 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
m18 <- lm(log(GPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN) + log(PlotNtoP), data=Peru_Plot_Master.data)



TableAICc <-  model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = aicc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
write.table(TableAICc, "TableAICc.xls", sep = "\t")

summary(m15)$adj.r.squared

#to pull out other info see summary(fit)$coefficients[,4]   summary(fit)$r.squared  summary(m12)$adj.r.squared


#TableQAICc <- model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = QAIC, rank.args = list(chat = deviance(m12) / df.residual(m12)))

#write.table(TableAICc, "TableQAICc.xls", sep = "\t")

#==========================================================================
##### NPP ######
# Here are the main models hypothesized . . .
m1 <- lm(log(NPP_Malhi_2017) ~ MAinvBT, data=Peru_Plot_Master.data) #climate models
m2 <- lm(log(NPP_Malhi_2017) ~ Precipitation.mm.yr.1., data=Peru_Plot_Master.data)
m3 <- lm(log(NPP_Malhi_2017) ~ log(Precipitation.mm.yr.1.), data=Peru_Plot_Master.data)
m4 <- lm(log(NPP_Malhi_2017) ~ log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m5 <- lm(log(NPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.), data=Peru_Plot_Master.data)
m6 <- lm(log(NPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m7 <- lm(log(NPP_Malhi_2017) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m8 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data) #Biomass model
m9 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data) # Biomass and climate
m10 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m11 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data) # TDT Biomass and Traits
m12 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)
m13 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole), data=Peru_Plot_Master.data)
m14 <- lm(log(NPP_Malhi_2017) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
m15 <- lm(log(NPP_Malhi_2017) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
m16 <- lm(log(NPP_Malhi_2017) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(NVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
m17 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)
m18 <- lm(log(NPP_Malhi_2017) ~ log(Aboveground_biomass2) + log(PhotosynthesisPerLeafN) + log(PlotNtoP), data=Peru_Plot_Master.data)


TableAICcNPP <-model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = aicc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
# worst models m15, m16, m7, m18, m12 . . .  # Best models m9, m14, m10  but note that m16 has the highest logLik value

write.table(TableAICcNPP, "TableAICcNPPc.xls", sep = "\t")

#model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = logLik, rank.args = list(chat = deviance(m12) / df.residual(m12)))
# highest logLik m16, m15, m15, m12  lowest logLik  m2, m3, m1, m8

#TableQAICcNPP <-model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = QAICc, rank.args = list(chat = deviance(m12) / df.residual(m12)))

#write.table(TableAICcNPP, "TableQAICcNPPc.xls", sep = "\t")

summary(m2)$adj.r.squared

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


