library(gemtraits)

con = connect_gemtraits_db()
photosyn = get_photosyn(con)

library(plyr)

per_tree_chem = ddply(photosyn, .(tree_id), summarize, mean_c_percent = mean(c_percent, na.rm = T),
                                                       mean_n_percent = mean(n_percent, na.rm = T))
