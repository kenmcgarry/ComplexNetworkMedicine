# drugstructure_gi.R

library(ChemmineR)
library(ChemmineOB)
library(ape)
library(sparcl)
library(cluster) # used for kmeans and silhoutte plot
library(xtable)
library(gplots) 

setwd("C:/R-files/bigfiles")
sdfset <- read.SDFset("structures.sdf") # load in huge file of chemical structures (approx 7,000)
valid <- validSDF(sdfset) # ensure invalid data removed
sdfset <- sdfset[valid]

## Assign drugbank IDs from datablock to cid slot
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
cid(sdfset) <- as.character(blockmatrix[,"DRUGBANK_ID"])

## Generate APset and FPset (note FPset: has better search performance)
apset <- sdf2ap(sdfset)
fpset <- desc2fp(apset, descnames=512, type="FPset")

# Overwrite drug names with our candidate drugs using drugbank_id
drugs <- unique(drug_list$drugbank_id)

# we have problems with some drugs that are not in FP e.g. DB05436, DB01115, DB05291, SO remove them!
badlist<- c("DB05436","DB01115","DB05291","DB07886","DB01398")
drugs <- setdiff(drugs, badlist)
fpdrugs <- fpset[drugs]   # extract our drugs chemical signatures from the many. We now have 189 drugs left

fpdrugs <- sample(fpdrugs)#randomize order of drugs
params <- genParameters(fpdrugs)  # params is used to calculate similarity scores
fpSim(fpset[[1]], fpset, top=25, parameters=params) 
clusters <- cmp.cluster(db=fpdrugs, save.distances="distmat.rda",cutoff = c(0.7, 0.8, 0.9), quiet = TRUE)

#dummy <- cmp.cluster(db=fpdrugs, cutoff=0, save.distances="distmat.rda", quiet=TRUE) 
load("distmat.rda") 
hc <- hclust(as.dist(distmat), method="single") 
hc[["labels"]] <- cid(fpdrugs) # Assign correct item labels 
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE) 

simMA <- sapply(cid(fpdrugs), function(x) fpSim(fpdrugs[x], fpdrugs, sorted=FALSE))
hc <- hclust(as.dist(1-simMA), method="single") 
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=FALSE)

heatmap.2(1-distmat, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), 
          col=colorpanel(40, "darkblue", "yellow", "white"), 
          density.info="none", trace="both")

# NEEDS TO BE DEBUGGED FROM THIS POINT ON TO GET CHEM SIMILARITY
results1 <- fpSim(fpdrugs[[71]], fpdrugs, top=77, parameters=params,method="Tanimoto") 
results2 <- fpSim(fpdrugs[[73]], fpdrugs, top=77, parameters=params,method="Tanimoto") 
results3 <- fpSim(fpdrugs[[74]], fpdrugs, top=77, parameters=params,method="Tanimoto") 
results1 <- cbind(id2name(rownames(results1)),results1) # USE drug names and NOT drugbank ID's
results2 <- cbind(id2name(rownames(results2)),results2)
results3 <- cbind(id2name(rownames(results3)),results3)
colnames(results1)[1] <- "name"
colnames(results2)[1] <- "name"
colnames(results3)[1] <- "name"
colnames(results1)[2] <- "sim1"
colnames(results2)[2] <- "sim2"
colnames(results3)[2] <- "sim3"

results1 <- as_tibble(results1)
results2 <- as_tibble(results2)
results3 <- as_tibble(results3)

chemsim1 <- results1 %>%
  dplyr::select(name,sim1)
chemsim1

chemsim2 <- results2 %>%
  dplyr::select(name,sim2)
chemsim2

chemsim3 <- results3 %>%
  dplyr::select(name,sim3)
chemsim3

chemsim <- chemsim1  %>%  dplyr::select(name,sim1) %>%
  dplyr::full_join(chemsim2,by="name")
chemsim

chemsim <- chemsim3  %>%  dplyr::select(name,sim3) %>%
  dplyr::full_join(chemsim,by="name")
chemsim

#calculate joint similarity score for the candidate drugs
chemsinjoint <- ((chemsim$sim3 * 33.33) + (chemsim$sim1 * 33.33) + (chemsim$sim2* 33.33))/100
chemsim <- cbind(chemsim,chemsinjoint)
chemsim <- arrange(chemsim,desc(chemsinjoint))

print.xtable(xtable(chemsim)) # displays tables for paper. 
#print.xtable(xtable(results2))
#print.xtable(xtable(results3))

# Convert drugbank ids to drugnames.
drugnames <- ID2name(drugs)

# Creates the similarity score matrix and cluster them.  
simMA <- sapply(cid(fpdrugs), function(x) fpSim(x=fpdrugs[x], fpdrugs, sorted=TRUE)) 
colnames(simMA)<-drugnames
rownames(simMA)<-drugnames


cl <- kmeans(simMA,10,nstart=50) #cl <- kmeans(simMA,10,nstart=5)
sk <- silhouette(cl$cl,dist(simMA))
plot(sk)

#par(mar=c(3, 3, 3, 3))
#hc <- hclust(as.dist(1-simMA), method="complete")
#plot(as.phylo(hc), cex = 0.9, label.offset = 0.01)#

y <- cutree(hc,20) #10
ColorDendrogram(hc,y=y,labels=drugnames,branchlength = 0.7,cex = 2)





