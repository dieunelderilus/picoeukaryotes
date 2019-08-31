#author : Dieunel  Derilus
#this script will generate the pairwise distance matrix  usung camberra index.
# be sure that the input contains the  count numbers of each edge in column and the spcies ID in row
#open the edge abundance table
a1 <- read.table("edge_abun_transpose.txt", header = TRUE)
#generate the pairwise distance matrix with the the method of interest
dm<-vegdist(a1,diag=T, method="canberra")
# save the distance matrix in a table, that will be used
