#require(private.BRI)
require(MASS)

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

Leaf_Nmass<-as.numeric(photosyn$n_percent)*as.numeric(photosyn$laminapetiole_drymass)
Leaf_Pmass<-as.numeric(photosyn$p_corrected_percent)*as.numeric(photosyn$laminapetiole_drymass)
Leaf_Cmass<-as.numeric(photosyn$c_percent)*as.numeric(photosyn$laminapetiole_drymass)
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

#First up, Leaf N

#remove indets from trees data
all_fp_trees<-all_fp_trees[which(all_fp_trees$fp_species_name!="Indet indet"),]



#1) Plot information
plots<-unique(all_fp_trees$plot_code)
sp_by_plot<-unique(cbind(as.character(all_fp_trees$plot_code),as.character(all_fp_trees$fp_species_name),as.character(all_fp_trees$fp_family_name)))
max_vals<-NULL
min_vals<-NULL
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
  
  if(length(traits_i)!=0){
  traits_i<-as.matrix(traits_i[[1]])#pulls out the highest ranking level of traits
  #traits_i<-as.data.frame(traits_i)#convert to data frame to allow easier indexing
  print(length(traits_i[,1]))
  dist_i<-fitdistr(x=as.numeric(as.vector(traits_i)),densfun = "Normal")
  mean_i<-dist_i$estimate[1]
  sd_i<-dist_i$estimate[2]
  output_Leaf_Nmass<-rbind(output_Leaf_Nmass,cbind(plot_i,sp_i,mean_i,sd_i))
  
  #calculate
  
  }else{
  #traits_i<-NA  
  #do something here if trait is NA?
  }

  
  
}#for loop for fitting distributions
