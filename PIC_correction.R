#this script run Phylogenetic independent contrast
# load datat
data <- read.csv("pC.csv", row.names = 1, header = T)
load library
library(ape)
#load the tree and th
tree=read.tree("mrabyes_enero17.nwk")

# Extract columns
p1 <- data[, "Genome_size"]
g1 <- data[, "p141"] # the pathway of interest
names(p1) <- names(g1) <- rownames(data)


#make a model
pPic <- pic(p1, tree)
gPic <- pic(g1, tree)

picModel <- lm(gPic ~ pPic)
# summarize model
summary(picModel)
