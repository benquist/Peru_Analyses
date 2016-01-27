
##########################
#Units stuff

#LMA =  g/m^2

#SLA =  m^2/g

#photosynthesis = umol/m^2/s

#lamina area = m^2#one of these is wrong, since it doesn't match the reported sla/lma

#lamina mass = g# one of these has to be wrong

#photosyn<-read.csv("data_12192015/photosyn.csv",colClasses = "character")

#plot(as.numeric(as.character(photosyn$laminapetiole_area))/as.numeric(as.character(photosyn$laminapetiole_drymass))~as.numeric(as.character(photosyn$sla_lamina_petiole)))
#as per the plot, the units differ by orders of magnitude, but are otherwise identical.  
#the sla and lma seem to be correct, since they agree with the bien output

#photosynthesis = umol/m^2/s

#The reported LMA and SLA values appear to come from a separate project associated with CHAMBASA.  In short, the lma values in the photosyn table were not automatically generated from the other columns directly.  These values of sla and lma come from a student who went through lots of images, judging how good they were, etc.and made corrections using the quality flags in the leaf area columns, and filtering based on those flags. For quality, there is a quality field in that table where "e” = exclude, "c" = caution, and "g" = good.  However, for the purpose of these analyses we may just want to report the LMA distribution using the ‘lma_lamina_petiole’ field as well as just recalculating lma which we can call ‘ma_lamina_petiole_calculated’ where you just divide the lamina dry mass by the lamina area as reported in those separate fields.  Make sense?

#Do we have enough info now to plot out the LMA and photosynthesis distributions across the gradient?

##########################

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
#We can focus on the main  functional traits - leaf N, P, C and leaf SLA or LMA. This also allows us to prune out duplicate measurements
#Remember LMA is just the inverse of SLA. Also, leaf photosynthesis per unit area and per unit mass. 
#We can also look at wood density but that  Is another dataset that needs to be found.

photosyn$fp_species_id<-as.character(photosyn$fp_species_id)
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet"),]#prune out the observations which are "indetermined"
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet indet"),]#prune out the observations which are "indetermined"
photosyn<-photosyn[which(photosyn$fp_species_name!="Indet indet"),]#prune out the observations which are "indetermined"
#out<-photosyn[which(photosyn$fp_species_name=="Indet indet"),]#prune out the observations which are "indetermined"

#convert units on cnp
photosyn$n_percent<-as.numeric(as.character(photosyn$n_percent))*.01
photosyn$p_corrected_percent<-as.numeric(as.character(photosyn$p_corrected_percent))*.01
photosyn$c_percent<-as.numeric(as.character(photosyn$c_percent))*.01


#there are duplicate measurements for many leaf values, so duplicate values and unwanted columns will be removed
photosyn2<-photosyn[c(-1,-26:-56)]# drop columns associated with photosynthesis measurements for now, since they lead to duplicates on other values

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
rm(i)
#Convert Peru SLA data to BIEN units

photosyn$sla_lamina_petiole<-(as.numeric(as.character(photosyn$sla_lamina_petiole))*1000)
photosyn$lma_lamina_petiole<-(as.numeric(as.character(photosyn$lma_lamina_petiole))/1000)
photosyn_cnpsla<-photosyn

########


#Calculate SLA and LMA from raw peru data: units still seem off by a factor of 10, but otherwise identical, no reason to use both.
################


#1) Plot information
plots<-unique(all_fp_trees$plot_code)
sp_by_plot<-unique(cbind(as.character(all_fp_trees$plot_code),as.character(all_fp_trees$fp_species_name),as.character(all_fp_trees$fp_family_name)))
##################################################

##NMass distribution fitting
output_Leaf_Nmass<-NULL
output_metadata<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
    traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$n_percent))
  
  
    
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
    traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$n_percent))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
    traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$n_percent))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
    traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$n_percent))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
    traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$n_percent))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
    traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$n_percent))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Nmass"),]
    traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
 
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Nmass"),]
    traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
 
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Nmass"),]  
     traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))        
     names(traits_i)<-c("a","b","c","d","e","f","g","h","i")  
  
      
      
  traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
  
  level_used<-names(traits_i[1])
  sample_size<-length(traits_i[[1]])
  output_metadata<-rbind(output_metadata,cbind(plot_i,sp_i,level_used,sample_size))
  
  if(length(traits_i)>0){
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
  
Nmass_fitting_metadata<-output_metadata  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i,level_used,output_metadata,plot_i,sample_size)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plots,shape1_i,shape2_i)

############################################

#Next,CMass

output_Leaf_Cmass<-NULL
output_metadata<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$c_percent))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$c_percent))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$c_percent))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$c_percent))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$c_percent))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$c_percent))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Cmass"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Cmass"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Cmass"),]  
  traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))  
  
  names(traits_i)<-c("a","b","c","d","e","f","g","h","i")
  
traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data      
level_used<-names(traits_i[1])   

if(length(traits_i)==0){
sample_size=0  
}

if(length(traits_i)>0){
sample_size<-length(traits_i[[1]])   
}

output_metadata<-rbind(output_metadata,cbind(plot_i,sp_i,level_used,sample_size))
  
  if(length(traits_i)>0){
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
Cmass_fitting_metadata<-output_metadata  

  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plot_i,shape1_i,shape2_i,level_used,sample_size)
##############################################
#PMass

output_Leaf_Pmass<-NULL
output_metadata<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$p_corrected_percent))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$p_corrected_percent))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$p_corrected_percent))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$p_corrected_percent))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$p_corrected_percent))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$p_corrected_percent))
  
  #BIEN bits need work
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Leaf Pmass"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Leaf Pmass"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Leaf Pmass"),]  
   traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))        
   names(traits_i)<-c("a","b","c","d","e","f","g","h","i")  
  
  
  
traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data      
level_used<-names(traits_i[1])   
if(length(traits_i)==0){ sample_size=0   }  

if(length(traits_i)>0){ sample_size<-length(traits_i[[1]])    }

output_metadata<-rbind(output_metadata,cbind(plot_i,sp_i,level_used,sample_size))
  
  if(length(traits_i)>0){
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
Pmass_fitting_metadata<-output_metadata  

  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,plot_i,shape1_i,shape2_i,level_used,sample_size)

#############################################
#photosyn$sla_lamina_petiole

output_Leaf_SLA<-NULL
output_metadata<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$sla_lamina_petiole))
  
  
  
  b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$sla_lamina_petiole))
  
  
  c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$sla_lamina_petiole))
  
  
  
  d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$sla_lamina_petiole))
  
  
  e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$sla_lamina_petiole))
  
  f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$sla_lamina_petiole))
  
  #BIEN bits
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Specific leaf area (SLA)"),]  
   traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))        
   names(traits_i)<-c("a","b","c","d","e","f","g","h","i")  
  
  
  
traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data      
level_used<-names(traits_i[1])   
if(length(traits_i)==0){ sample_size=0   }  
if(length(traits_i)>0){ sample_size<-length(traits_i[[1]])    }
output_metadata<-rbind(output_metadata,cbind(plot_i,sp_i,level_used,sample_size))
 
  ###
  if(length(traits_i)>0){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  }
  
  if(length(traits_i)>1){
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i))
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
    dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "gamma",method="BFGS")
    shape_i<-dist_i$estimate[1]
    rate_i<-dist_i$estimate[2]
    #mean_i<-dist_i$estimate[1]
    #sd_i<-dist_i$estimate[2]
    output_Leaf_SLA<-rbind(output_Leaf_SLA,cbind(plot_i,sp_i,shape_i,rate_i))
    
    #calculate
    
  }else{#if>1
    
    if(length(traits_i)==1){
      traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
      #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
      print(length(traits_i[,1]))
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
      #dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mle"))
      shape_i<-as.numeric(traits_i[1])
      rate_i<-NA
      #mean_i<-dist_i$estimate[1]
      #sd_i<-dist_i$estimate[2]
      output_Leaf_SLA<-rbind(output_Leaf_SLA,cbind(plot_i,sp_i,shape_i,rate_i))
      
      #calculate
      
    }#if=1
  }#else
  
  ### 
sla_fitting_metadata<-output_metadata  

  
  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,shape_i,rate_i,plot_i,level_used,sample_size)



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

trait_list<-CMass
plot(density(na.omit(trait_list[[1]])),col="red",ylim=c(0,25),main = "C Mass",xlab="Percent")  
lines(density(na.omit(trait_list[[2]])),col="orange")  
lines(density(na.omit(trait_list[[3]])),col="yellow")  
lines(density(na.omit(trait_list[[4]])),col="light green")  
lines.default(density(na.omit(trait_list[[5]])),col="dark green")  
lines(density(na.omit(trait_list[[6]])),col="light blue")  
lines(density(na.omit(trait_list[[7]])),col="dark blue")  
lines(density(na.omit(trait_list[[8]])),col="violet")  
lines(density(na.omit(trait_list[[9]])),col="purple")  
lines(density(na.omit(trait_list[[10]])),col="black")  

#red SPD-01 
#orange TAM-06
#yellow ACJ-01 
#light green ESP-01 
#dark green PAN-02 
#light blue PAN-03 
#dark blue SPD-02 
#violet TRU-04 
#purple WAY-01 
#black TAM-05

trait_list<-Nmass
plot(density(na.omit(trait_list[[1]])),col="red",ylim=c(0,120),main = "N Mass",xlab="Percent")  
lines(density(na.omit(trait_list[[2]])),col="orange")  
lines(density(na.omit(trait_list[[3]])),col="yellow")  
lines(density(na.omit(trait_list[[4]])),col="light green")  
lines.default(density(na.omit(trait_list[[5]])),col="dark green")  
lines(density(na.omit(trait_list[[6]])),col="light blue")  
lines(density(na.omit(trait_list[[7]])),col="dark blue")  
lines(density(na.omit(trait_list[[8]])),col="violet")  
lines(density(na.omit(trait_list[[9]])),col="purple")  
lines(density(na.omit(trait_list[[10]])),col="black")  

trait_list<-Pmass
plot(density(na.omit(trait_list[[1]])),col="red",ylim=c(0,1500),main = "P Mass",xlab = "Percent")  
lines(density(na.omit(trait_list[[2]])),col="orange")  
lines(density(na.omit(trait_list[[3]])),col="yellow")  
lines(density(na.omit(trait_list[[4]])),col="light green")  
lines.default(density(na.omit(trait_list[[5]])),col="dark green")  
lines(density(na.omit(trait_list[[6]])),col="light blue")  
lines(density(na.omit(trait_list[[7]])),col="dark blue")  
lines(density(na.omit(trait_list[[8]])),col="violet")  
lines(density(na.omit(trait_list[[9]])),col="purple")  
lines(density(na.omit(trait_list[[10]])),col="black")  
rm(trait_list)
###############################
peru_draws_beta_distribution<-function(output_file,nreps){
plots<-unique(all_fp_trees$plot_code)
draws<-list()
for(i in 1:length(plots)){
  plot_i<-as.character(plots[i])
  occurrences_i<-sp_by_plot[sp_by_plot[,1]==plot_i,]
  trait_draws<-NULL
  for(s in 1:length(occurrences_i[,1])){
    sp_s<-occurrences_i[,2][s]
    occ_s<-as.numeric(occurrences_i[,4][s])    
    #vals<-output_Leaf_Cmass[which(output_Leaf_Cmass[,1]==plot_i & output_Leaf_Cmass[,2]==sp_s),]
    vals<-output_file[which(output_file[,1]==plot_i & output_file[,2]==sp_s),]
    shape1_s<-as.numeric(vals[3])
    shape2_s<-as.numeric(vals[4])
     trait_vals_s<-NULL
  for(r in 1:nreps){    
    if(is.na(shape2_s)==FALSE){
      #trait_vals_s<-rnorm(n=occ_s,mean=mean_s,sd = sd_s)
      trait_vals_r<-rbeta(n=occ_s,shape1=shape1_s,shape2 = shape2_s)
      trait_vals_s<-rbind(trait_vals_s,trait_vals_r)
      
    }
    
    if(is.na(shape2_s)==TRUE){
      
      
      trait_vals_r<-matrix(shape1_s,occ_s,nrow=1) 
      trait_vals_s<-rbind(trait_vals_s,trait_vals_r)
      
    }
    
    
    
  }#trait draw for plot i
  trait_draws<-cbind(trait_draws,trait_vals_s) 
      
  }#s occurrences
  draws[[i]]<-trait_draws
  #add code here to add trait draws to list    
}#cmass trait draw
return(draws)
rm(i,occ_s,plot_i,s,shape1_s,shape2_s,sp_s,trait_draws,trait_vals_s,vals)
}
###########
peru_draws_gamma_distribution<-function(output_file,nreps){
  plots<-unique(all_fp_trees$plot_code)
  draws<-list()
  for(i in 1:length(plots)){
    plot_i<-as.character(plots[i])
    occurrences_i<-sp_by_plot[sp_by_plot[,1]==plot_i,]
    trait_draws<-NULL
    for(s in 1:length(occurrences_i[,1])){
      sp_s<-occurrences_i[,2][s]
      occ_s<-as.numeric(occurrences_i[,4][s])    
      #vals<-output_Leaf_Cmass[which(output_Leaf_Cmass[,1]==plot_i & output_Leaf_Cmass[,2]==sp_s),]
      vals<-output_file[which(output_file[,1]==plot_i & output_file[,2]==sp_s),]
      
      
      
      shape_s<-as.numeric(vals[3])
      rate_s<-as.numeric(vals[4])
      trait_vals_s<-NULL
      for(r in 1:nreps){    
        if(is.na(rate_s)==FALSE){
          #trait_vals_s<-rnorm(n=occ_s,mean=mean_s,sd = sd_s)
          trait_vals_r<-rgamma(n=occ_s,shape=shape_s,rate = rate_s)
          #trait_vals_r<-rgamma(n=occ_s,shape1=shape1_s,shape2 = shape2_s)
          trait_vals_s<-rbind(trait_vals_s,trait_vals_r)
          
        }
        
        if(is.na(rate_s)==TRUE){
          
          
          trait_vals_r<-matrix(shape_s,occ_s,nrow=1) 
          trait_vals_s<-rbind(trait_vals_s,trait_vals_r)
          
        }
        
        
        
      }#trait draw for plot i
      trait_draws<-cbind(trait_draws,trait_vals_s) 
      
    }#s occurrences
    draws[[i]]<-trait_draws
    #add code here to add trait draws to list    
  }#cmass trait draw
  return(draws)
  rm(i,occ_s,plot_i,s,shape1_s,shape2_s,sp_s,trait_draws,trait_vals_s,vals)
}

#########
cmass_draws<-peru_draws_beta_distribution(output_file = output_Leaf_Cmass,nreps=1000)
pmass_draws<-peru_draws_beta_distribution(output_file = output_Leaf_Pmass,nreps=1000)
nmass_draws<-peru_draws_beta_distribution(output_file = output_Leaf_Nmass,nreps=1000)
sla_draws<-peru_draws_gamma_distribution(output_file = output_Leaf_SLA,nreps = 1000)

require(moments)
plots<-unique(all_fp_trees$plot_code)

peru_draw_analysis<-function(draws_file){
  output<-NULL
  for(i in 1:length(draws_file)){
  #draws_i<-draws_file[[i]]
  #draws_naomit<-draws_file[[i]][ , ! apply( draws_file[[i]] , 2 , function(x) all(is.na(x)) ) ]#remove na columns
    
  draws_naomit<-draws_file[[i]][ , ! apply( draws_file[[i]] , 2 , function(x) all(is.na(x)) ) ]#remove na columns
  plot<-as.character(plots[i])
  mean<-apply(draws_naomit,1,mean)
  variance<-apply(draws_naomit,1,var)
  skewness<-apply(draws_naomit,1,skewness)
  kurtosis<-apply(draws_naomit,1,kurtosis)  
  mean95<-sort(mean,decreasing = FALSE)[(.025*length(mean)+1):(0.975*length(mean))]
  variance95<-sort(variance,decreasing = FALSE)[(.025*length(variance)+1):(0.975*length(variance))]
  skewness95<-sort(skewness,decreasing = FALSE)[(.025*length(skewness)+1):(0.975*length(skewness))]
  kurtosis95<-sort(kurtosis,decreasing = FALSE)[(.025*length(kurtosis)+1):(0.975*length(kurtosis))]
  mean<- cbind(min(mean95),mean(mean),max(mean95))
  variance<- cbind(min(variance95),mean(variance),max(variance95))
  skewness<- cbind(min(skewness95),mean(skewness),max(skewness95))
  kurtosis<- cbind(min(kurtosis95),mean(kurtosis),max(kurtosis95))
  
  #output_i<-rbind(mean,variance,skewness,kurtosis)
  output_i<-cbind(plot,mean,variance,skewness,kurtosis)
  
  
  colnames(output_i)<-c("Plot","Mean Lower","Mean Mean","Mean Upper",
                        
                        "Variance Lower","Variance Mean","Variance Upper",
                        "Skewness Lower","Skewness Mean","Skewness Upper",
                        "Kurtosis Lower","Kurtosis Mean","Kurtosis Upper")
  output<-rbind(output,output_i)
  
}#for i loop
    row.names(output)<-NULL
    output<-as.data.frame(output)
    return(output)
  
}

peru_moments_cmass<-peru_draw_analysis(draws_file = cmass_draws)
peru_moments_pmass<-peru_draw_analysis(draws_file = pmass_draws)
peru_moments_nmass<-peru_draw_analysis(draws_file = nmass_draws)
peru_moments_sla<-peru_draw_analysis(draws_file = sla_draws)

##############
#BIEN photosyn measurements
#Area-based photosynthesis (Aarea)
#Mass-based photosynthesis (Amass)

bien_area_photo<-BIEN.trait.trait("Area-based photosynthesis (Aarea)")
bien_mass_photo<-BIEN.trait.trait("Mass-based photosynthesis (Amass)")
unique(bien_area_photo$unit) #all in the units "Âµmol.m-2.s-1"
unique(bien_mass_photo$unit) #all in the units "Âµmol.g-1.s-1"

#Peru photosyn units umol/m^2/s
#BIEN histogram
hist(as.numeric(as.character(bien_area_photo$trait_value)))#mean about 10, range about 0-40
hist(as.numeric(as.character(bien_mass_photo$trait_value)))#mean about .2, range about 0-.7

#need to re-load photosyn since I dropped the relevant columns earlier
photosyn<-read.csv("data_12192015/photosyn.csv",colClasses = "character")
hist(as.numeric(as.character(photosyn$photosynthesis))) #peru looks like area-based photo

photosyn_amax<-photosyn[which(photosyn$pm_type=="AMAX"),]
photosyn_asat<-photosyn[which(photosyn$pm_type=="ASAT"),]
hist(as.numeric(as.character(photosyn_amax$photosynthesis)))
hist(as.numeric(as.character(photosyn_asat$photosynthesis)))


##########################
#photosynthesis

#bien_area)photo
#photosyn_amax
output_Leaf_photo<-NULL
output_metadata<-NULL
for( i in 1:length(sp_by_plot[,1])){
  sp_i<-sp_by_plot[,2][i]  
  plot_i<-sp_by_plot[,1][i]
  genus_i<-strsplit(sp_i,split = " ")[[1]][1]
  family_i<-sp_by_plot[,3][i]
  
  traits_i<-list()
  #look within plot
  
  a_plot_data_sp_i<-photosyn_amax[which(photosyn_amax$plot_code==plot_i & photosyn_amax$fp_species_name==sp_i) ,]
  traits_i[[1]]<-na.omit(as.matrix(a_plot_data_sp_i$photosynthesis))
  
  
  
  b_study_data_sp_i<-photosyn_amax[which(photosyn_amax$fp_species_name==sp_i) ,]#neat trick= using "which" prevents NA lines from showing up
  traits_i[[2]]<-na.omit(as.matrix(b_study_data_sp_i$photosynthesis))
  
  
  c_plot_data_genus_i<-photosyn_amax[which(photosyn_amax$plot_code==plot_i & photosyn_amax$fp_genus_name==genus_i) ,]
  traits_i[[3]]<-na.omit(as.matrix(c_plot_data_genus_i$photosynthesis))
  
  
  
  d_study_data_genus_i<-photosyn_amax[which(photosyn_amax$fp_genus_name==genus_i) ,]
  traits_i[[4]]<-na.omit(as.matrix(d_study_data_genus_i$photosynthesis))
  
  
  e_plot_data_family_i<-photosyn_amax[which(photosyn_amax$plot_code==plot_i & photosyn_amax$fp_family_name==family_i) ,]
  traits_i[[5]]<-na.omit(as.matrix(e_plot_data_family_i$photosynthesis))
  
  f_study_data_family_i<-photosyn_amax[which(photosyn_amax$fp_family_name==family_i) ,]
  traits_i[[6]]<-na.omit(as.matrix(f_study_data_family_i$photosynthesis))
  
  #BIEN bits
  
  g_bien_data_sp_i<-bien_traits[which(bien_traits$species==sp_i & bien_traits$trait_name=="Area-based photosynthesis (Aarea)"),]
  traits_i[[7]]<-na.omit(as.matrix(g_bien_data_sp_i$trait_value))
  
  h_bien_data_genus_i<-bien_traits[which(bien_traits$genus==genus_i & bien_traits$trait_name=="Area-based photosynthesis (Aarea)"),]
  traits_i[[8]]<-na.omit(as.matrix(h_bien_data_genus_i$trait_value))
  
  i_bien_data_family_i<-bien_traits[which(bien_traits$family==family_i & bien_traits$trait_name=="Area-based photosynthesis (Aarea)"),]  
   traits_i[[9]]<-na.omit(as.matrix(i_bien_data_family_i$trait_value))        
   names(traits_i)<-c("a","b","c","d","e","f","g","h","i")  
  
  
  
traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data      
level_used<-names(traits_i[1])   
if(length(traits_i)==0){ sample_size=0   }  
if(length(traits_i)>0){ sample_size<-length(traits_i[[1]])    }
output_metadata<-rbind(output_metadata,cbind(plot_i,sp_i,level_used,sample_size))
  
  ###
  if(length(traits_i)>0){
    traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  }
  
  if(length(traits_i)>1){
    #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
    print(length(traits_i))
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
    #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
    dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "gamma",method="BFGS")
    shape_i<-dist_i$estimate[1]
    rate_i<-dist_i$estimate[2]
    #mean_i<-dist_i$estimate[1]
    #sd_i<-dist_i$estimate[2]
    output_Leaf_photo<-rbind(output_Leaf_photo,cbind(plot_i,sp_i,shape_i,rate_i))
    
    #calculate
    
  }else{#if>1
    
    if(length(traits_i)==1){
      traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
      #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
      print(length(traits_i[,1]))
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
      #dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Beta",list(shape1=1,shape2=1))
      #dist_i<-fitdist(data=(as.numeric(as.vector(traits_i))),distr = "beta", method = ("mle"))
      shape_i<-as.numeric(traits_i[1])
      rate_i<-NA
      #mean_i<-dist_i$estimate[1]
      #sd_i<-dist_i$estimate[2]
      output_Leaf_photo<-rbind(output_Leaf_photo,cbind(plot_i,sp_i,shape_i,rate_i))
      
      #calculate
      
    }#if=1
  }#else
  
  ### 
  
photosynthesis_fitting_metadata<-output_metadata  

  
}#for loop for fitting distributions

rm(a_plot_data_sp_i,b_study_data_sp_i,c_plot_data_genus_i,d_study_data_genus_i,e_plot_data_family_i,f_study_data_family_i,
   g_bien_data_sp_i,h_bien_data_genus_i,i_bien_data_family_i)
rm(traits_i,sp_i,genus_i,family_i,i,dist_i,shape_i,rate_i,plot_i,level_used,sample_size)

#######
photosynthesis_draws<-peru_draws_gamma_distribution(output_file = output_Leaf_photo,nreps = 1000)

peru_moments_cmass<-peru_draw_analysis(draws_file = cmass_draws)
peru_moments_pmass<-peru_draw_analysis(draws_file = pmass_draws)
peru_moments_nmass<-peru_draw_analysis(draws_file = nmass_draws)
peru_moments_sla<-peru_draw_analysis(draws_file = sla_draws)
peru_moments_photosynthesis<-peru_draw_analysis(draws_file = photosynthesis_draws)

#######
#save(pmass_draws,file = "pmass_draws")
#save(cmass_draws,file = "cmass_draws")
#save(nmass_draws,file = "nmass_draws")
#save(photosynthesis_draws,file = "photosynthesis_draws")
#save(sla_draws,file = "sla_draws")

#write.csv(peru_moments_photosynthesis,"peru_moments_photosynthesis.csv")
#write.csv(peru_moments_sla,"peru_moments_sla.csv")
#write.csv(peru_moments_nmass,"peru_moments_nmass.csv")
#write.csv(peru_moments_pmass,"peru_moments_pmass.csv")
#write.csv(peru_moments_cmass,"peru_moments_cmass.csv")

#write.csv(Cmass_fitting_metadata,"Cmass_fitting_metadata.csv")
#write.csv(Nmass_fitting_metadata,"Nmass_fitting_metadata.csv")
#write.csv(Pmass_fitting_metadata,"Pmass_fitting_metadata.csv")
#write.csv(sla_fitting_metadata,"sla_fitting_metadata.csv")
#write.csv(photosynthesis_fitting_metadata,"photosynthesis_fitting_metadata.csv")
#a_plot_data_sp_i
#b_study_data_sp_i
#c_plot_data_genus_i
#d_study_data_genus_i
#e_plot_data_family_i
#f_study_data_family_i
#g_bien_data_sp_i
#h_bien_data_genus_i
#i_bien_data_family_i


load("pmass_draws")
load("cmass_draws")
load("nmass_draws")
load("sla_draws")
load("photosynthesis_draws")

peru_moments_cmass<-read.csv("peru_moments_cmass.csv")
peru_moments_nmass<-read.csv("peru_moments_nmass.csv")
peru_moments_pmass<-read.csv("peru_moments_pmass.csv")
peru_moments_photosynthesis<-read.csv("peru_moments_photosynthesis.csv")
peru_moments_sla<-read.csv("peru_moments_sla.csv")
#######
#stuff to add

#1) it would be great to have a multi panel plot showing the continuous distributions 
#for each of the traits across the gradient.

#2) for each plot we will need metrics including :  
    #total number of species, 
    #total number of species with locally measures within plot traits,  
    #number of species with substituted traits (coming from somewhere on the gradient, and 
    #number of species with a BIEN trait, and number of species with a genus and family trait mean).  In short some sort of measure of how much filling was done.
    #number of species without data or where dist. couldn't be fit


#####


#PLot order:
#plot_code	Elevation (m)
#TAM-06	215 1
#TAM-05	223 2
#PAN-02	595 3
#PAN-03	859 4
#SPD-02	1494  5
#SPD-01	1713  6
#TRU-04	2719  7
#ESP-01	2868  8
#WAY-01	3045  9
#ACJ-01	3537  10

#[Order in lists]:plot:order in increasing elevation
#[1] "SPD-01" 6
#[2] "TAM-06" 1
#[3] "ACJ-01" 10
#[4] "ESP-01" 8
#[5] "PAN-02" 3
#[6] "PAN-03" 4
#[7] "SPD-02" 5
#[8] "TRU-04" 7
#[9] "WAY-01" 9
#[10] "TAM-05"  2


#add

require(grDevices)
colfunc <- colorRampPalette(c("red", "blue"))
colors_for_plots<-colfunc(10)


trait_list<-cmass_draws
max_value<-max(na.omit(as.numeric(as.character(photosyn$c_percent))))
min_value<-min(na.omit(as.numeric(as.character(photosyn$c_percent))))
lwd=2
dev.new(width=5, height=4)
plot(density(na.omit(as.vector(trait_list[[2]][which(trait_list[[2]]>min_value &  trait_list[[2]]<min_value)]))),col=colors_for_plots[1],ylim=c(0,25),xlim=c(0.35,0.6),main = "Carbon",xlab="Percent",lwd=lwd,lty=1)  
lines(density(na.omit(as.vector(trait_list[[10]]))),col=colors_for_plots[2],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[5]]))),col=colors_for_plots[3],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[6]]))),col=colors_for_plots[4],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[7]]))),col=colors_for_plots[5],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[1]]))),col=colors_for_plots[6],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[8]]))),col=colors_for_plots[7],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[4]]))),col=colors_for_plots[8],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[9]]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[3]]))),col=colors_for_plots[10],lwd=lwd,lty=1)  

legend("topright",lty=c(1,2,3,4,5,5,4,3,2,1),col=colors_for_plots[10:1],legend=c(
  "ACJ-01:3537m","WAY-01:3045m","ESP-01:2868m","TRU-04:2719m","SPD-01:1713m",
 "SPD-02:1494m","PAN-03:859m","PAN-02:595m","TAM-05:223m","TAM-06:215m"
  ),lwd=lwd)

trait_list<-nmass_draws
lwd=2
dev.new(width=5, height=4)
plot(density(na.omit(as.vector(trait_list[[2]]))),col=colors_for_plots[1],ylim=c(0,180),xlim=c(0,.05),main = "Nitrogen",xlab="Percent",lwd=lwd,lty=1)  
lines(density(na.omit(as.vector(trait_list[[10]]))),col=colors_for_plots[2],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[5]]))),col=colors_for_plots[3],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[6]]))),col=colors_for_plots[4],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[7]]))),col=colors_for_plots[5],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[1]]))),col=colors_for_plots[6],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[8]]))),col=colors_for_plots[7],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[4]]))),col=colors_for_plots[8],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[9]]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[3]]))),col=colors_for_plots[10],lwd=lwd,lty=1)  

legend("topright",lty=c(1,2,3,4,5,5,4,3,2,1),col=colors_for_plots[10:1],legend=c(
  "ACJ-01:3537m","WAY-01:3045m","ESP-01:2868m","TRU-04:2719m","SPD-01:1713m",
  "SPD-02:1494m","PAN-03:859m","PAN-02:595m","TAM-05:223m","TAM-06:215m"
),lwd=lwd)


trait_list<-pmass_draws
lwd=2
dev.new(width=5, height=4)
plot(density(na.omit(as.vector(trait_list[[2]]))),col=colors_for_plots[1],ylim=c(0,1800),xlim=c(0,.004),main = "Phosphorus",xlab="Percent",lwd=lwd,lty=1)  
lines(density(na.omit(as.vector(trait_list[[10]]))),col=colors_for_plots[2],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[5]]))),col=colors_for_plots[3],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[6]]))),col=colors_for_plots[4],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[7]]))),col=colors_for_plots[5],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[1]]))),col=colors_for_plots[6],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[8]]))),col=colors_for_plots[7],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[4]]))),col=colors_for_plots[8],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[9]]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[3]]))),col=colors_for_plots[10],lwd=lwd,lty=1)  

legend("topright",lty=c(1,2,3,4,5,5,4,3,2,1),col=colors_for_plots[10:1],legend=c(
  "ACJ-01:3537m","WAY-01:3045m","ESP-01:2868m","TRU-04:2719m","SPD-01:1713m",
  "SPD-02:1494m","PAN-03:859m","PAN-02:595m","TAM-05:223m","TAM-06:215m"
),lwd=lwd)

trait_list<-sla_draws
lwd=2
dev.new(width=5, height=4)
plot(density(na.omit(as.vector(trait_list[[2]]))),col=colors_for_plots[1],xlim=c(0,30),ylim=c(0,.5),main = "SLA",xlab="m^2/g",lwd=lwd,lty=1)  
lines(density(na.omit(as.vector(trait_list[[10]]))),col=colors_for_plots[2],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[5]]))),col=colors_for_plots[3],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[6]]))),col=colors_for_plots[4],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[7]]))),col=colors_for_plots[5],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[1]]))),col=colors_for_plots[6],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[8]]))),col=colors_for_plots[7],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[4]]))),col=colors_for_plots[8],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[9]]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
#lines(density(na.omit(as.vector(trait_list[[9]][which(trait_list[[9]]<200)]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[3]]))),col=colors_for_plots[10],lwd=lwd,lty=1)  

legend("topright",lty=c(1,2,3,4,5,5,4,3,2,1),col=colors_for_plots[10:1],legend=c(
  "ACJ-01:3537m","WAY-01:3045m","ESP-01:2868m","TRU-04:2719m","SPD-01:1713m",
  "SPD-02:1494m","PAN-03:859m","PAN-02:595m","TAM-05:223m","TAM-06:215m"
),lwd=lwd)

####
#photosyn_cnpsla use for c,n,p,sla
#photosun_amax use for amax

outlier_trimmer<-function(draws_file,min,max){
for(i in 1:length(draws_file)){
draws_file[[i]][which(draws_file[[i]]>max )]<-max
draws_file[[i]][which(draws_file[[i]]<min )]<-min
}#i loop
rm(i)
return(draws_file)
}#edit fx

#cmass-remove outliers
cmass_max<-max(na.omit(photosyn_cnpsla$c_percent))  
cmass_min<-min(na.omit(photosyn_cnpsla$c_percent))
cmass_draws_edited<-outlier_trimmer(draws_file = cmass_draws,min = cmass_min,max = cmass_max)
peru_moments_cmass_edited<-peru_draw_analysis(draws_file = cmass_draws_edited)

#nmass-remove outliers
nmass_max<-max(na.omit(photosyn_cnpsla$n_percent))  
nmass_min<-min(na.omit(photosyn_cnpsla$n_percent))
nmass_draws_edited<-outlier_trimmer(draws_file = nmass_draws,min = nmass_min,max = nmass_max)
peru_moments_nmass_edited<-peru_draw_analysis(draws_file = nmass_draws_edited)

#pmass-remove outliers
pmass_max<-max(na.omit(photosyn_cnpsla$p_corrected_percent))  
pmass_min<-min(na.omit(photosyn_cnpsla$p_corrected_percent))
pmass_draws_edited<-outlier_trimmer(draws_file = pmass_draws,min = pmass_min,max = pmass_max)
peru_moments_pmass_edited<-peru_draw_analysis(draws_file = pmass_draws_edited)

#sla-remove outliers
sla_max<-max(na.omit(photosyn_cnpsla$sla_lamina_petiole))  
sla_min<-min(na.omit(photosyn_cnpsla$sla_lamina_petiole))
sla_draws_edited<-outlier_trimmer(draws_file = sla_draws,min = sla_min,max = sla_max)
peru_moments_sla_edited<-peru_draw_analysis(draws_file = sla_draws_edited)

#photosynthesis-remove outliers
photosynthesis_max<-as.numeric(max(na.omit(photosyn_amax$photosynthesis)))  
photosynthesis_min<-as.numeric(min(na.omit(photosyn_amax$photosynthesis)))
photosynthesis_draws_edited<-outlier_trimmer(draws_file = photosynthesis_draws,min = photosynthesis_min,max = photosynthesis_max)
peru_moments_photosynthesis_edited<-peru_draw_analysis(draws_file = photosynthesis_draws_edited)

#write.csv(peru_moments_photosynthesis_edited,"peru_moments_photosynthesis_edited.csv")
#write.csv(peru_moments_sla_edited,"peru_moments_sla_edited.csv")
#write.csv(peru_moments_nmass_edited,"peru_moments_nmass_edited.csv")
#write.csv(peru_moments_pmass_edited,"peru_moments_pmass_edited.csv")
#write.csv(peru_moments_cmass_edited,"peru_moments_cmass_edited.csv")




############
trait_list<-sla_draws_edited
lwd=2
dev.new(width=5, height=4)
plot(density(na.omit(as.vector(trait_list[[2]]))),col=colors_for_plots[1],xlim=c(0,30),ylim=c(0,.35),main = "SLA",xlab="m^2/g",lwd=lwd,lty=1)  
lines(density(na.omit(as.vector(trait_list[[10]]))),col=colors_for_plots[2],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[5]]))),col=colors_for_plots[3],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[6]]))),col=colors_for_plots[4],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[7]]))),col=colors_for_plots[5],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[1]]))),col=colors_for_plots[6],lwd=lwd,lty=5)  
lines(density(na.omit(as.vector(trait_list[[8]]))),col=colors_for_plots[7],lwd=lwd,lty=4)  
lines(density(na.omit(as.vector(trait_list[[4]]))),col=colors_for_plots[8],lwd=lwd,lty=3)  
lines(density(na.omit(as.vector(trait_list[[9]]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
#lines(density(na.omit(as.vector(trait_list[[9]][which(trait_list[[9]]<200)]))),col=colors_for_plots[9],lwd=lwd,lty=2)  
lines(density(na.omit(as.vector(trait_list[[3]]))),col=colors_for_plots[10],lwd=lwd,lty=1)  

legend("topright",lty=c(1,2,3,4,5,5,4,3,2,1),col=colors_for_plots[10:1],legend=c(
  "ACJ-01:3537m","WAY-01:3045m","ESP-01:2868m","TRU-04:2719m","SPD-01:1713m",
  "SPD-02:1494m","PAN-03:859m","PAN-02:595m","TAM-05:223m","TAM-06:215m"
),lwd=lwd)
