#require(private.BRI)

bien_traits<-read.csv("bien_traits_12192015")
all_fp_trees<-read.csv("data_12192015/all_fp_trees.csv")
all_trees<-read.csv("data_12192015/all_trees.csv")
photosyn<-read.csv("data_12192015/photosyn.csv")

#Required data

#1) Plot information
#2) Species abundances in each plot
#3) Required/available trait data
#4) Distribution parms for each species/plot combination
#     -Preference of data: species data from plot>all plots>bien >genus data from plot>all plots>bien >family data from plot>all plots>bien
#5) Sample traits from distribution according to abundances


