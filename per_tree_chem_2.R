#######################################
##  Peru Chambasa data  12/7/15
## Brian J. Enquist
## Modified code from Allie Shenkin
##
######################################


library(gemtraits)

con = connect_gemtraits_db()
photosyn = get_photosyn(con)

library(plyr)

## original code
#
# per_tree_chem = ddply(photosyn, .(tree_id), summarize, 
#                   mean_c_percent = mean(c_percent, na.rm = T),
#                   mean_n_percent = mean(n_percent, na.rm = T),
#                   dbh = dbh[1])

##  Longer list of associated traits per tree

per_tree_chem_2 = ddply(photosyn, .(tree_id), summarize, mean_c_percent = mean(c_percent, na.rm = T),
                        mean_n_percent = mean(n_percent, na.rm = T),
                        mean_p_percent = mean(p_corrected_percent, na.rm = T),
                        mean_c13_delta = mean(c13_delta, na.rm = T),
                        dbh = dbh[1])

######  Average trait values per species
# note subsitute csp_full_name for fp_species_name to get full species name

per_species_chem_2 = ddply(photosyn, .(fp_species_name), summarize, mean_c_percent = mean(c_percent, na.rm = T),
                        mean_n_percent = mean(n_percent, na.rm = T),
                        mean_p_percent = mean(p_corrected_percent, na.rm = T),
                        mean_c13_delta = mean(c13_delta, na.rm = T),
                        dbh = dbh[1])

######  Average trait values per genus

per_genus_chem_2 = ddply(photosyn, .(fp_genus_name), summarize, 
                         mean_c_percent = mean(c_percent, na.rm = T),
                         mean_n_percent = mean(n_percent, na.rm = T),
                         mean_p_percent = mean(p_corrected_percent, na.rm = T),
                         mean_c13_delta = mean(c13_delta, na.rm = T),
                         dbh = dbh[1])

######  Average trait values per plot

per_plot_chem_2 = ddply(photosyn, .(plot_code), summarize, 
                          mean_c_percent = mean(c_percent, na.rm = T),
                          var_c_percent = var(c_percent, na.rm = T),
                          mean_n_percent = mean(n_percent, na.rm = T),
                          var_n_percent = var(n_percent, na.rm = T),
                          mean_p_percent = mean(p_corrected_percent, na.rm = T),
                          var_p_percent = var(p_corrected_percent, na.rm = T),
                          mean_c13_delta = mean(c13_delta, na.rm = T),
                          var_c13_percent = var(c13_delta, na.rm = T),
                          dbh = dbh[1])




### traits to add  
## photosynthesis c13_delta n15_delta lamina_area lamina_drymass lma_lamina lma_lamina_petiole  fp_species_name


############################
##
## add non-chambasa trees to chem dataframe
## Modified code from Allie Shenkin
############################

library(gemtraits)

con = connect_gemtraits_db()
photosyn = get_photosyn(con)

library(plyr)

per_tree_chem = ddply(photosyn, .(tree_id), summarize, 
                      mean_c_percent = mean(c_percent, na.rm = T),
                      mean_n_percent = mean(n_percent, na.rm = T),
                      dbh = dbh[1])

# Chambasa sometimes sampled trees outside of the plot, or sampled small lianas.  
# These "Individuals" are not sampled as part of the regular census, and hence do 
# not show up in the #forestplots database, and should not be used when calculating stem
# -weighted figures (though #their traits may be used to calculate means, etc).

So, all_fp_trees just includes trees in the forestplots census. all_trees includes censused trees as well as those just sampled in chambasa.

all_trees = get_trees(con)
all_fp_trees = subset(all_trees, ! grepl("^I", all_trees$tree_code))
all_trees_with_chem = join(all_fp_trees[,c("tree_id","plot_code","tree_code","fp_species_name","dbh")], per_tree_chem, by = "tree_id", type = "left")



