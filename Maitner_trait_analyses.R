#require(private.BRI)

bien_traits<-read.csv("bien_traits_12192015")
all_fp_trees<-read.csv("data_12192015/all_fp_trees.csv")
all_trees<-read.csv("data_12192015/all_trees.csv")
photosyn<-read.csv("data_12192015/photosyn.csv")

#Required data

#1) Plot information
plots<-unique(all_fp_trees$plot_code)
sp_by_plot<-unique(cbind(as.character(all_fp_trees$plot_code),as.character(all_fp_trees$fp_species_name),as.character(all_fp_trees$fp_family_name)))

for( i in 1:length(sp_by_plot[,1])){
sp_i<-sp_by_plot[,2][i]  
plot_i<-sp_by_plot[,1][i]
genus_i<-strsplit(sp_i,split = " ")[[1]][1]
family_i<-sp_by_plot[,3][i]

traits_i<-NULL
#look within plot
traits_i[[1]]<-as.matrix(a_plot_data_sp_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_species_name==sp_i) ,])
traits_i[[2]]<-as.matrix(b_study_data_sp_i<-photosyn[which(photosyn$fp_species_name==sp_i) ,])#neat trick= using "which" prevents NA lines from showing up
traits_i[[3]]<-as.matrix(c_plot_data_genus_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_genus_name==genus_i) ,])
traits_i[[4]]<-as.matrix(d_study_data_genus_i<-photosyn[which(photosyn$fp_genus_name==genus_i) ,])
traits_i[[5]]<-as.matrix(e_plot_data_family_i<-photosyn[which(photosyn$plot_code==plot_i & photosyn$fp_family_name==family_i) ,])
traits_i[[6]]<-as.matrix(f_study_data_family_i<-photosyn[which(photosyn$fp_family_name==family_i) ,])
#g_bien_data_sp_i<-
#h_bien_data_genus_i<-
#i_bien_data_family_i<-  

traits_i<-traits_i[which(lengths(traits_i)>0)]#prunes list down to hierarchical levels with data
traits_i_2<-traits_i[[1]]


}#for loop fitting distributions





#2) Species abundances in each plot
#3) Required/available trait data
#4) Distribution parms for each species/plot combination
#     -Preference of data: species data from plot>all plots>bien >genus data from plot>all plots>bien >family data from plot>all plots>bien
#5) Sample traits from distribution according to abundances


