# This script  is used to plot a  newick tree with and aabrplot showing the proteome size of each species
library("ape")
library(phytools)
p1=read.table("proteome.txt",header =TRUE, col.names = 1)
tree=read.tree("mrbayes_2019.nwk")
cw<-reorder(tree)
bsize<-setNames(p1[,1],rownames(p1))

par(lab=c(5,10,7))
plotTree.barplot(cw,p1,list(fsize=1,lwd=2),
                 list(col="black",space=0.7,xlab="Proteome size(codons)"))
add.scale.bar(cex = 0.7, font = 2, col = "black")
