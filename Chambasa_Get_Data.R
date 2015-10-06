# install postgreSQL, install.packages("foo", dependencies=...)
install.packages("RPostgreSQL")
install.packages("devtools")
library(devtools)

install.packages("/Users/brianjenquist/GitHub/gemtraits", type="source", repos=NULL)

library(gemtraits)

# See instructions in R package gemtraits.pdf

# first thing that you do is to create a connection with the database
# type con = connect_gemtraits_db()

con <- connect_gemtraits_db()

#queries run by commands in the R file db_queries.r
#trees_crowns = get_trees_with_crowns(con)
trees_dbh <- get_trees_with_dbh(con)
LeafChem = get_leaf_chemistry(con)
WoodDensity = get_wood_density(con)
LeafParts = get_leaf_parts(con)
