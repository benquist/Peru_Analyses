##############################
##  Peru Chambasa data  12/7/15
## Modified code from Allie Shenkin
##
#############################


library(gemtraits)

con = connect_gemtraits_db()
photosyn = get_photosyn(con)

library(plyr)

## original code
#
# per_tree_chem = ddply(photosyn, .(tree_id), summarize, mean_c_percent = mean(c_percent, na.rm = T),
#                                                       mean_n_percent = mean(n_percent, na.rm = T),
#                                                       dbh = dbh[1])

##  Longer list of associated traits per tree

per_tree_chem_2 = ddply(photosyn, .(tree_id), summarize, mean_c_percent = mean(c_percent, na.rm = T),
                        mean_n_percent = mean(n_percent, na.rm = T),
                        mean_p_percent = mean(p_corrected_percent, na.rm = T),
                        mean_c13_delta = mean(c13_delta, na.rm = T),
                        dbh = dbh[1])

######  Average trait values per species

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
                          mean_n_percent = mean(n_percent, na.rm = T),
                          mean_p_percent = mean(p_corrected_percent, na.rm = T),
                          mean_c13_delta = mean(c13_delta, na.rm = T),
                          dbh = dbh[1])




### traits to add  
## photosynthesis c13_delta n15_delta lamina_area lamina_drymass lma_lamina lma_lamina_petiole  fp_species_name

