#require(private.BRI)
require(MASS)#Needed for distribution fitting
require(fitdistrplus)
#need to find a new distribution-or use "while" to re-draw from normal dist. until all values are positive
#sla dist needs to be 0 to pos. infin
#cnp are percent, so use dist bounded by 0 and 1: beta seems likely


#############################
#Due to data format differences in bien and peru, I'm going to try and do this one trait at a time instead:
bien_traits<-read.csv("bien_traits_12192015")
all_fp_trees<-read.csv("data_12192015/all_fp_trees.csv")
all_trees<-read.csv("data_12192015/all_trees.csv")
photosyn<-read.csv("data_12192015/photosyn.csv",colClasses = "character")
#convert photosyn to match BIEN format
#We can focus on the main  functional traits - leaf N, P, C and leaf SLA or LMA. 
#Remember LMA is just the inverse of SLA. Also, leaf photosynthesis per unit area and per unit mass. 
#We can also look at wood density but that  Is another dataset that needs to be found.

Leaf_Nmass<-(as.numeric(photosyn$n_percent))#npc values will be in percent, BIEN data converted to match
Leaf_Pmass<-(as.numeric(photosyn$p_corrected_percent))
Leaf_Cmass<-(as.numeric(photosyn$c_percent))
Specific_leaf_area_SLA<-as.numeric(photosyn$sla_lamina_petiole)
photosyn<-cbind(photosyn,Leaf_Nmass,Leaf_Pmass,Leaf_Cmass,Specific_leaf_area_SLA)
rm(Leaf_Cmass,Leaf_Pmass,Leaf_Nmass,Specific_leaf_area_SLA)
photosyn$fp_species_id<-as.character(photosyn$fp_species_id)
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet"),]#prune out the observations which are "indetermined"
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet indet"),]#prune out the observations which are "indetermined"
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet indet"),]#prune out the observations which are "indetermined"
#out<-photosyn[which(photosyn$fp_species_name=="Indet indet"),]#prune out the observations which are "indetermined"



#there are duplicate measurements for many leaf values, so duplicate values and unwanted columns will be removed

photosyn2<-as.data.frame(cbind(photosyn$leaf_id,photosyn$fp_species_name,photosyn$fp_genus_name,photosyn$fp_family_name,photosyn$plot_code,photosyn$Leaf_Nmass,photosyn$Leaf_Pmass,photosyn$Leaf_Cmass,photosyn$Specific_leaf_area_SLA,deparse.level = 2))
names(photosyn2)<-c("leaf_id","fp_species_name","fp_genus_name","fp_family_name","plot_code","Leaf Nmass","Leaf Pmass","Leaf Cmass","Specific leaf area (SLA)")
photosyn3<-unique(photosyn2)
photosyn<-photosyn3
rm(photosyn2,photosyn3)


#remove indets from trees data
all_fp_trees<-all_fp_trees[which(all_fp_trees$fp_species_name!="Indet indet"),]
###

#standardize BIEN cnp measurements to percent

bien_traits$trait_value<-as.numeric(as.character(bien_traits$trait_value))
for(i in 1:length(bien_traits[,1])){
  #print(i)
  if(is.na(bien_traits$trait_name[i])==FALSE){
  if(bien_traits$trait_name[i]=='Leaf Pmass'|bien_traits$trait_name[i]=='Leaf Nmass'){
  
    if(bien_traits$unit[i]=="%"){
      bien_traits$trait_value[i]<-bien_traits$trait_value[i]*0.01  
      }#if  %  
    
  if(bien_traits$unit[i]=="mg/g"){
    bien_traits$trait_value[i]<-bien_traits$trait_value[i]*0.001  
    bien_traits$unit[i]<-"%"
  }#if  mg/g  
  }#if n or p mass
  }#name is not NA
}#bien standardizing loop

#1) Plot information
plots<-unique(all_fp_trees$plot_code)
sp_by_plot<-unique(cbind(as.character(all_fp_trees$plot_code),as.character(all_fp_trees$fp_species_name),as.character(all_fp_trees$fp_family_name)))
##################################################

##NMass distribution fitting
output_Leaf_Nmass<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
    traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$`Leaf Nmass`))
  
  
    
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
    traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$`Leaf Nmass`))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
    traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$`Leaf Nmass`))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
    traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$`Leaf Nmass`))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
    traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$`Leaf Nmass`))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
    traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$`Leaf Nmass`))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Nmass"),]
    traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
 
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Nmass"),]
    traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
 
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Nmass"),]  
    traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))  
  
      
      
  traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
  
  if(length(traits_i)>1){
  traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  }
  
  if(length(traits_i)>1){
  #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
  print(length(traits_i))
  #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
  #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
  dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mme"))
  shape1_i<-dist_i$estimate[1]
  shape2_i<-dist_i$estimate[2]
  #mean_i<-dist_i$estimate[1]
  #sd_i<-dist_i$estimate[2]
  output_Leaf_Nmass<-rbind(output_Leaf_Nmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
  
  #calculate
  
  }else{#if>1
  
  if(length(traits_i)==1){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i[,1]))
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
    #dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mle"))
    shape1_i<-as.numeric(traits_i[1])
    shape2_i<-NA
    #mean_i<-dist_i$estimate[1]
    #sd_i<-dist_i$estimate[2]
    output_Leaf_Nmass<-rbind(output_Leaf_Nmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
    
    #calculate
    
  }#if=1
}#else
  
  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plots,shape1_i,shape2_i)

############################################

#Next,CMass

output_Leaf_Cmass<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$`Leaf Cmass`))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$`Leaf Cmass`))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$`Leaf Cmass`))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$`Leaf Cmass`))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$`Leaf Cmass`))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$`Leaf Cmass`))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Cmass"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Cmass"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Cmass"),]  
  traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))  
  
  
  
  traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
  
  if(length(traits_i)>1){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  }
  
  if(length(traits_i)>1){
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i))
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
    dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mme"))
    shape1_i<-dist_i$estimate[1]
    shape2_i<-dist_i$estimate[2]
    #mean_i<-dist_i$estimate[1]
    #sd_i<-dist_i$estimate[2]
    output_Leaf_Cmass<-rbind(output_Leaf_Cmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
    
    #calculate
    
  }else{#if>1
    
    if(length(traits_i)==1){
      traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
      #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
      print(length(traits_i[,1]))
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
      #dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mle"))
      shape1_i<-as.numeric(traits_i[1])
      shape2_i<-NA
      #mean_i<-dist_i$estimate[1]
      #sd_i<-dist_i$estimate[2]
      output_Leaf_Cmass<-rbind(output_Leaf_Cmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
      
      #calculate
      
    }#if=1
  }#else
  
  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plot_i,shape1_i,shape2_i)
##############################################
#PMass

output_Leaf_Pmass<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$`Leaf Pmass`))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$`Leaf Pmass`))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$`Leaf Pmass`))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$`Leaf Pmass`))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$`Leaf Pmass`))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$`Leaf Pmass`))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Pmass"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Pmass"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Pmass"),]  
  traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))  
  
  
  
  traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
  
  if(length(traits_i)>1){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  }
  
  if(length(traits_i)>1){
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i))
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
    dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mme"))
    shape1_i<-dist_i$estimate[1]
    shape2_i<-dist_i$estimate[2]
    #mean_i<-dist_i$estimate[1]
    #sd_i<-dist_i$estimate[2]
    output_Leaf_Pmass<-rbind(output_Leaf_Pmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
    
    #calculate
    
  }else{#if>1
    
    if(length(traits_i)==1){
      traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
      #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
      print(length(traits_i[,1]))
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
      #dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mle"))
      shape1_i<-as.numeric(traits_i[1])
      shape2_i<-NA
      #mean_i<-dist_i$estimate[1]
      #sd_i<-dist_i$estimate[2]
      output_Leaf_Pmass<-rbind(output_Leaf_Pmass,cbind(plot_i,sp_i,shape1_i,shape2_i))
      
      #calculate
      
    }#if=1
  }#else
  
  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plot_i,shape1_i,shape2_i)

#############################################
#photosyn$`Specific leaf area (SLA)`

output_Leaf_SLA<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$`Specific leaf area (SLA)`))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$`Specific leaf area (SLA)`))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$`Specific leaf area (SLA)`))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$`Specific leaf area (SLA)`))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$`Specific leaf area (SLA)`))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$`Specific leaf area (SLA)`))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]  
  traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))  
  
  
  
  traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
  
  if(length(traits_i)!=0){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i[,1]))
    dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    mean_i<-dist_i$estimate[1]
    sd_i<-dist_i$estimate[2]
    output_Leaf_SLA<-rbind(output_Leaf_SLA,cbind(plot_i,sp_i,mean_i,sd_i))
    
    #calculate
    
  }else{
    #traits_i<-NA  
    #do something here if trait is NA?
  }
  
  
  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,mean_i,sd_i,plot_i,plots)



######################################################

#Now, we need a list of species per plot with the occurrences numbers

occurrences<-NULL
for(i in 1:length(sp_by_plot[,1])){
sp_i<-sp_by_plot[,2][i]  
plot_i<-sp_by_plot[,1][i]
occurrences_i<-length(which(all_fp_trees$plot_code==plot_i & all_fp_trees$fp_species_name==sp_i))
occurrences<-c(occurrences,occurrences_i)
  
}

sp_by_plot<-cbind(sp_by_plot,occurrences)
rm(occurrences,i,occurrences_i,plot_i,sp_i)

#####################################################

#Now, we just need to draw values according to the trait distributions
plots<-unique(all_fp_trees$plot_code)
CMass<-list()

for(i in 1:length(plots)){
plot_i<-as.character(plots[i])
occurrences_i<-sp_by_plot[sp_by_plot[,1]==plot_i,]
trait_draws<-NULL
for(s in 1:length(occurrences_i[,1])){
sp_s<-occurrences_i[,2][s]
occ_s<-as.numeric(occurrences_i[,4][s])    
vals<-output_Leaf_Cmass[which(output_Leaf_Cmass[,1]==plot_i & output_Leaf_Cmass[,2]==sp_s),]
shape1_s<-as.numeric(vals[3])
shape2_s<-as.numeric(vals[4])

if(is.na(shape2_s)==FALSE){
#trait_vals_s<-rnorm(n=occ_s,mean=mean_s,sd = sd_s)
trait_vals_s<-rbeta(n=occ_s,shape1=shape1_s,shape2 = shape2_s)
trait_draws<-c(trait_draws,trait_vals_s)

}

if(is.na(shape2_s)==TRUE){


trait_vals_s<-matrix(shape1_s,occ_s) 
trait_draws<-c(trait_draws,trait_vals_s)

}



}#trait draw for plot i
CMass[[i]]<-trait_draws
  
#add code here to add trait draws to list  
}#cmass trait draw

rm(i,occ_s,plot_i,s,shape1_s,shape2_s,sp_s,trait_draws,trait_vals_s,vals)

#hist(CMass[[3]])

#########################################

#N Mass
#Now, we just need to draw values according to the trait distributions
Nmass<-list()

for(i in 1:length(plots)){
  plot_i<-as.character(plots[i])
  occurrences_i<-sp_by_plot[sp_by_plot[,1]==plot_i,]
  trait_draws<-NULL
  for(s in 1:length(occurrences_i[,1])){
    sp_s<-occurrences_i[,2][s]
    occ_s<-as.numeric(occurrences_i[,4][s])    
    vals<-output_Leaf_Nmass[which(output_Leaf_Nmass[,1]==plot_i & output_Leaf_Nmass[,2]==sp_s),]
    shape1_s<-as.numeric(vals[3])
    shape2_s<-as.numeric(vals[4])
    
    if(is.na(shape2_s)==FALSE){
      #trait_vals_s<-rnorm(n=occ_s,mean=mean_s,sd = sd_s)
      trait_vals_s<-rbeta(n=occ_s,shape1=shape1_s,shape2 = shape2_s)
      trait_draws<-c(trait_draws,trait_vals_s)
      
    }
    
    if(is.na(shape2_s)==TRUE){
      
      
      trait_vals_s<-matrix(shape1_s,occ_s) 
      trait_draws<-c(trait_draws,trait_vals_s)
      
    }
    
    
    
  }#trait draw for plot i
  Nmass[[i]]<-trait_draws
  
  #add code here to add trait draws to list  
}#Nmass trait draw

rm(i,occ_s,plot_i,s,shape1_s,shape2_s,sp_s,trait_draws,trait_vals_s,vals)

#hist(Nmass[[3]])


#################################


#P Mass
#Now, we just need to draw values according to the trait distributions
Pmass<-list()

for(i in 1:length(plots)){
  plot_i<-as.character(plots[i])
  occurrences_i<-sp_by_plot[sp_by_plot[,1]==plot_i,]
  trait_draws<-NULL
  for(s in 1:length(occurrences_i[,1])){
    sp_s<-occurrences_i[,2][s]
    occ_s<-as.numeric(occurrences_i[,4][s])    
    vals<-output_Leaf_Pmass[which(output_Leaf_Pmass[,1]==plot_i & output_Leaf_Pmass[,2]==sp_s),]
    shape1_s<-as.numeric(vals[3])
    shape2_s<-as.numeric(vals[4])
    
    if(is.na(shape2_s)==FALSE){
      #trait_vals_s<-rnorm(n=occ_s,mean=mean_s,sd = sd_s)
      trait_vals_s<-rbeta(n=occ_s,shape1=shape1_s,shape2 = shape2_s)
      trait_draws<-c(trait_draws,trait_vals_s)
      
    }
    
    if(is.na(shape2_s)==TRUE){
      
      
      trait_vals_s<-matrix(shape1_s,occ_s) 
      trait_draws<-c(trait_draws,trait_vals_s)
      
    }
    
    
    
  }#trait draw for plot i
  Pmass[[i]]<-trait_draws
  
  #add code here to add trait draws to list  
}#Pmass trait draw

rm(i,occ_s,plot_i,s,shape1_s,shape2_s,sp_s,trait_draws,trait_vals_s,vals)

#hist(Pmass[[5]])

