# gi_functions.R
library(ChemmineR)
library(ChemmineOB)
library(ape)
library(sparcl)
library(cluster) # used for kmeans and silhoutte plot
library(xtable)
library(gplots) 
library(scatterplot3d) 
library(igraph)
library(ROCR)
library(VennDiagram)
library(ggplot2)
library(linkcomm)

## --------------------- FUNCTION DEFINITIONS -----------------------

# get_drug_names() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drug_names <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls) 
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  if(length(ilist) > 0){
    for (j in 1:length(ilist)){
      cat("\ndrug",j,"is", ilist[j])
    }
    return(ilist)
  }else{
    cat("\n","Sorry, no drugs found...check umls code is correct for your disease")
    return(NULL)}
}

# get_drugs_plus() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# The restricted list of drugs we cant use is passed to this function.
# It returns the drug names plus other information.
get_drugs_plus <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls) 
  mydrugs <- ilist # temp storage
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  
  if(length(ilist) == 0){  # stuck in a 2nd null test as zero length data creeping in and causing crashes.
    cat("\n","No drugs found...check umls code is correct for your disease or maybe no drugs for this disease")
    return(NULL)}
  
    options(warn=-1)
    ilist <- filter(mydrugs,drugbank_name == ilist) # causes warning message
    options(warn=0)
    ilist <- ilist[!duplicated(ilist$drugbank_name),]
    
  if(nrow(ilist) > 0){
    for (j in 1:nrow(ilist)){
      cat("\ndrug",j,"is", ilist[j,2])
    }
    ilist <- select(ilist,drugbank_id,drugbank_name,umls_cui_from_meddra,meddra_name)
    return(ilist)
  }else{
    cat("\n","No drugs found...check umls code is correct for your disease or maybe no drugs for this disease")
    #cat("\n i is...",j)
    return(NULL)}
}

# Convert the ID to umls code in order to access indications and obtain drug_id and drug_name, 
# digestive dataframe is passed as id; use ID to match with MeshID and get umls
id2umls <- function(id){
  x <- vector(mode="character",length=nrow(id))
  y <- mappings[1,] # instantiate a temporary vector, 
    
  for (i in 1:nrow(id)){
    x[i] <- id$ID[i] 
    tempy <- filter(mappings, meshId == x[i])  # if tempy has meshid then keep, else ignore and do save it
    if(nrow(tempy) > 0){
      y[i,] <- tempy[1,]
      y[i,]$umls <- gsub("umls:","",y[i,]$umls) # remove "umls:" from string
      #cat("\ni is ...",i)
    }
   
  }
  y <- na.omit(y)  # get rid of records with NA where no umls Id's exist for the Mesh id's
  return(y)
}


# Adds the MESH code into drugs_list dataframe, useful info when relating drugs to disease categories
# 
add_meshcode <- function(drugs){
  x <- vector(mode="character",length=nrow(drugs))
  y <- vector(mode="character",length=nrow(drugs))
  #drugplus <- 0
  
  for (i in 1:nrow(drugs)){
    tempx <- paste("umls:",drugs[i,3],sep="") # add the bloody "umls:" back in, and some magic numbers.
    tempy <- filter(mappings, umls == tempx) 
    tempz <- filter(digestive, ID == tempy$meshId)
    x[i] <- tempz$MeSH
    y[i] <- tempz$ID
  }
  
  drugsplus <- cbind(drugs,x)   # add the x or Mesh values
  drugsplus <- cbind(drugsplus,y)   # add the y or ID values
  colnames(drugsplus)[5] <- "MeSH"   # change from x and y to better names
  colnames(drugsplus)[6] <- "ID"
  
  return(drugsplus)
}


# R provides a tail and head command to view last six and first six elements, so why not the middle six?
middle <- function(mydata) {
  len <- nrow(mydata)
  startpoint <- round(len/2)
  endpoint <- startpoint+5
  mydata[startpoint:endpoint,]
  
}


# Supply get_disease_genes() with a "drug_list" which has "umls_cui_from_meddra" field to tie 
# with "diseaseId" in "disgene", For each disease see what genes are implicated and return dataframe
get_disease_genes <- function(mydrugs){
  implicated <-0 # instantiate, unfortunelay makes a zero entry, delete this at end of function.
  
  disgene$diseaseId <- gsub("umls:","",disgene$diseaseId) # get rid of the bloody "umls:" from disgene for good!
  mydrugs <- mydrugs[!duplicated(mydrugs[,c('umls_cui_from_meddra','meddra_name')]),]  # just keep unique diseases
  mydrugs <- select(mydrugs,umls_cui_from_meddra,meddra_name, MeSH,ID)  # keep name and ids etc
  
  for (i in 1:nrow(mydrugs)){
    tempx <- mydrugs[i,1]  # seek out the umls code
    tempy <- filter(disgene, diseaseId == tempx) 
    if(nrow(tempy) > 0){
      cat("\nFound gene(s)!")
      tempy <- select(tempy,diseaseId, geneName, diseaseName)
      implicated <- rbind(tempy,implicated)
    }
  }
  
  implicated <- implicated[-nrow(implicated),]            # last entry is zero so remove it
  return(implicated)
}


# Convert drugbank IDs (DB00035) to drugnames (Desmopressin)
ID2name <- function(DBid){
  thenames <- vector(mode="character",length=length(DBid))
  
  for (i in 1:length(DBid)){
    dname <- filter(indications, drugbank_id == DBid[i])
    dname <- dname[!duplicated(dname[,'drugbank_name']),]
    thenames[i] <- dname$drugbank_name
  }
  return(thenames)
}

# plot_chemsim() creates plots generated from the data computed by drugstructure_gi.R code .
plot_chemsim <- function(){
  plot.new()
  
  hc <- hclust(as.dist(1-simMA), method="complete") 
  #par(cex=0.1)
  heatmap.2((1-simMA), Rowv=as.dendrogram(hc), 
            Colv=as.dendrogram(hc), 
            #col=greenred(10),
            keysize = 2,
            key=TRUE,
            #col=bluered(256),
            col=colorpanel(40, "white","yellow", "darkblue"), 
            density.info="none", trace="none")
  
  # Creates the similarity score matrix and cluster them.  
  colnames(simMA)<-drugnames
  rownames(simMA)<-drugnames
  hc <- hclust(as.dist(1-simMA), method="complete") 
  plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=FALSE)
  
  cl <- kmeans(simMA,10,nstart=50) #cl <- kmeans(simMA,10,nstart=5)
  sk <- silhouette(cl$cl,dist(simMA))
  plot(sk)
  
  y <- cutree(hc,25) #10
  par(cex=0.8)
  ColorDendrogram(hc,y=y,labels=drugnames,branchlength = 0.7,cex = 0.7)  
  
}


# get diseases that are not C06 but connected to C06 by shared genes 
get_linked_diseases <- function(dgenes){

  linked_diseases <- disgene[1,] # instantiate before use
  
  for (i in 1:length(dgenes)){
    #cat("\nlength of dgenes is ",length(dgenes))
    gene <- dgenes[i]
    glist <- filter(disgene, geneName == gene)
    if(nrow(glist) > 0){
      tempdis <- setdiff(glist$diseaseName,disease_umls$meshId)
      for (j in 1:length(tempdis)){
        glist <- filter(glist, diseaseName == tempdis[j])
        linked_diseases <- rbind(glist,linked_diseases)
      }
      #tempdis <- select(tempdis,diseaseId, geneName, diseaseName)
      #linked_diseases <- rbind(tempdis,linked_diseases)
    }
  }
  
  linked_diseases <- linked_diseases[-nrow(linked_diseases),]     # last entry is zero so remove it
  linked_diseases <- linked_diseases[!duplicated(linked_diseases[,'diseaseName']),]   # get rid of the many duplicates
  linked_diseases <- select(linked_diseases,diseaseId,geneId,geneName,diseaseName)  # drop "score" variable
  linked_diseases <- arrange(linked_diseases,diseaseName)  # sort alphabetically
  return(linked_diseases)
}


# input  a gene or list of genes and get all diseases asscoiated with these genes 
# really best just using JUST one gene.
get_all_linked_diseases <- function(dgenes){
  all_diseases <- disgene[1,] # instantiate before use
  
  for (i in 1:length(dgenes)){
    gene <- dgenes[i]
    glist <- filter(disgene, geneName == gene)
    if(nrow(glist) > 0){
        all_diseases <- rbind(glist,all_diseases)
    }
  }
  
  all_diseases <- all_diseases[-nrow(all_diseases),]     # last entry is fake so remove it
  all_diseases <- all_diseases[!duplicated(all_diseases[,'diseaseName']),]   # get rid of the many duplicates
  all_diseases <- select(all_diseases,diseaseId,geneId,geneName,diseaseName)  # drop "score" variable
  all_diseases <- arrange(all_diseases,diseaseName)  # sort alphabetically
  return(all_diseases)
}




