
#################################
### Multiple Regression fits
### and tables of model output from
### competed models
#################################

install.packages("caret")
library(caret)

#GPP_new_estimate
# NPP_new and biomass and boltzman temp  **use
m3 <- lm(log(GPP_new_estimate)~ MAinvBT  + log(Aboveground_biomass2), data=Peru_Plot_Master.data)
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


m3 <- lm(log(GPP_new_estimate)~ log(mean_sla_lamina_petiole)  + log(mean_n_percent) + log(Aboveground_biomass2), data=Peru_Plot_Master.data)
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

# NPP_new and biomass and boltzman temp  **use
m3 <- lm(log(GPP_new_estimate)~ MAinvBT + log(Aboveground_biomass2), data=Peru_Plot_Master.data)
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


# NPP_new and biomass and boltzman temp  **use
m3 <- lm(log(NPP_new)~ MAinvBT + log(Precipitation.mm.yr.1.)+ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
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



### with traits
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


### GPP and biomass and boltzman temp  and plot N:P 
m3 <- lm(log(GPP_new_estimate)~ log(Aboveground_biomass2) + MAinvBT + log(PlotNtoP), data=Peru_Plot_Master.data)

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

# from glumati lm(formula = log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data = Peru_Plot_Master.data)

m3 <- lm(formula = log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data = Peru_Plot_Master.data)

summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 

# WOW!!! low vifs - low AICc, high AdjR2~0.95  too good to be true?


m3 <- lm(formula = log(NPP_new) ~ log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data = Peru_Plot_Master.data)

summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3) 



## no temperature - just biomass and traits . . . models do not do as well as models with tmperature

#photo
m3 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(mean_photosynthesis), data=Peru_Plot_Master.data)

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

#N:P
m3 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2)  + log10(PlotNtoP), data=Peru_Plot_Master.data)


summary(m3) 
AIC(m3)
AICc(m3)
modelEffectSizes(m3)  #lm.sunSquares is depreciated
varImp(m3, scale = FALSE)
avPlots(m3)
crPlots(m3)
confint(m3)
vif(m3)

### PhotoPerLeafN
m3 <- lm(log10(GPP_new_estimate)~ log10(Aboveground_biomass2) + log10(PhotosynthesisPerLeafN), data=Peru_Plot_Master.data)

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


##########################

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
fit3 <- glmulti(log(NPP_new) ~ SolarRadiation.GJ.m.2.yr.1. + log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + log(mean_sla_lamina_petiole) + log(SLAVarianceMean) + MAinvBT + log(mean_n_percent) + log(mean_c_percent) + log(mean_p_percent) + log(PlotNtoP) + log(PhotosynthesisPerLeafN) + log(mean_photosynthesis) + log(mean_p_percent) + log(PhotoVarianceMean), data = Peru_Plot_Master.data, crit=aicc, level=1, fitfunc=glm, method="h")
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


##########################

library(MuMIn)
m1 <- lm(log(GPP_new_estimate) ~  log(Aboveground_biomass2) + log(mean_sla_lamina_petiole) + log(mean_n_percent), data=Peru_Plot_Master.data)

m2 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + MAinvBT, data=Peru_Plot_Master.data)
m3 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2), data=Peru_Plot_Master.data)
m4 <- lm(log(GPP_new_estimate) ~ log(Aboveground_biomass2) + log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)
m5 <- lm(log(GPP_new_estimate) ~ MAinvBT + log(PhotoVarianceMean), data=Peru_Plot_Master.data)
m6 <- lm(log(GPP_new_estimate) ~ log(mean_sla_lamina_petiole) + log(mean_n_percent) + log(SLAVarianceMean)+ log(Aboveground_biomass2) , data=Peru_Plot_Master.data)
m7 <- lm(log(GPP_new_estimate) ~ log(SolarRadiation.GJ.m.2.yr.1.) +  log(Precipitation.mm.yr.1.) + MAinvBT, data=Peru_Plot_Master.data)


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



TableAICc <-  model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = aicc, rank.args = list(chat = deviance(m12) / df.residual(m12)))
write.table(TableAICc, "TableAICc.xls", sep = "\t")

summary(m12)$adj.r.squared

# worst models m7, m15, m18, m16, m2, m3, m6, . . .  # Best models m12, m9,  then m8, m19, m14, m11,  
#to pull out other info see 
# summary(fit)$coefficients[,4]   summary(fit)$r.squared  summary(m12)$adj.r.squared

# using QAIC by changing rank statement to rank = QAIC, . . . best models m16, m12, m15, m7

TableQAICc <- model.sel(m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18, rank = QAIC, rank.args = list(chat = deviance(m12) / df.residual(m12)))

write.table(TableAICc, "TableQAICc.xls", sep = "\t")


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

write.table(TableAICcNPP, "TableQAICcNPPc.xls", sep = "\t")

summary(m9)$adj.r.squared

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


