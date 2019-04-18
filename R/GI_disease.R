# GI_disease.R

setwd("C:/R-files/disease")   # set working directory, Yours will be different
source("GI_user_functions.R")  

library(linkcomm)
library(igraph)
#library(dplyr)

disease <- read.csv("C06data.csv")  # read CSV file of disease to proteins

# some basic frequency counts of what we have. Only C06 diseases with implicated genes are used
length(unique(disease$geneName))
length(unique(disease$diseaseName))


lc <- getLinkCommunities(disease, hcmethod="single")  # Figure out the link communities/clusters

lc2 <- newLinkCommsAt(lc, cutat=0.7)  # MAke the cut in dendrogram to create more clusters
  
# now stick in the commands for the plots etc.
getNodesIn(lc,clusterids = 1)
plot(lc,type="graph",clusterids=1)

# some nice plots
plot(lc,type="members")
plot(lc,type="summary")

# How related are our clusters? very similar ones could be noise
cr <- getClusterRelatedness(lc,hcmethod = "ward.D")

# phase two: biological pathway involvement
# https://bioconductor.org/packages/3.7/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
require(clusterProfiler)   # "require" is another way of installing/loading packages
require(ReactomePA)

proteinlist <- getNodesIn(lc,clusterids=1)  # select cluster to analyse, 
mypaths <- kegg_analysis(proteinlist)
dim(mypaths)
head(mypaths,10)
barplot(mypaths, showCategory=10)




# phase three: complex network analysis
# This converts dataframe to an igraph object, we can do a load of stuff with this later on.
C06_graph <- graph.data.frame(d = disease, directed = FALSE)
plot(C06_graph)  # a complete black out!


# Use this later
tkid <- tkplot(C06_graph) #tkid is the id of the tkplot that will open
l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
tk_close(tkid, window.close = T)
plot(C06_graph, layout=l)

# Another thing is to see how the drugs are involved, could they be repurposed from other diseases and 
# vice-versa? 

