# gi_functions.R
library(dplyr) 
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
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(GOplot)
library(scales)
library(ontologySimilarity)
library(ontologyIndex)
library(rentrez)
library(stringr)
# huge amount of data gets loaded in here onwards!
library(ontologySimilarity)
library(ontologyIndex)
library(infotheo)
data(go)
data(gene_GO_terms)
data(GO_IC)

## --------------------- FUNCTION DEFINITIONS -----------------------

# Makes first letter uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# the cosine formula gives a similarity. It is 1 if the vectors are pointing in the same direction. 
# Distance measures need the value to be 0 when vectors are the same so 1 - similarity = distance. 
# Many uses need distance rather than similarity 
cosineDist1 <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}


# This vesrion uses the dot product matrix, to compute the cosine similarity matrix: 
# Input S is the matrix of dot product, dt is dataset. Divide the dot product by the norms of vectors. 
cosineDist2 <- function(S){
  doc_norm <- apply(as.matrix(dt),1,function(x) norm(as.matrix(x),"f")) 
  divide_one_norm <- S/doc_norm 
  cosine <- t(divide_one_norm)/doc_norm
  return (cosine)
}


# get_drug_names() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drug_names <- function(umls,rlist) {
  umls <- gsub("umls:","",umls) # get rid of the bloody "umls:" in the string
  cat("\nCode is..",umls)
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


# get_linked_diseases() that are not C06 but connected to C06 by shared genes. Modified to include all
# diseases associated with a disorder, 30/11/17 Need to inclued drugs associated with disorders.
get_linked_diseases <- function(dgenes){
  linked_diseases <- disgene[1,] # instantiate before use
  
  for (i in 1:length(dgenes)){
    #cat("\nlength of dgenes is ",length(dgenes))
    gene <- dgenes[i] # get genes individually and see what disorders they linked with
    glist <- filter(disgene, geneName == gene)  # This bit is OK
    if(nrow(glist) > 0){
      linked_diseases <- rbind(linked_diseases,glist) }
  }
  
  linked_diseases <- arrange(linked_diseases,diseaseName)  # sort alphabetically
  linked_diseases <- linked_diseases[,c(1,2,4,6)]            # keep only key variables
  # The complicated line below removes duplicate entries, there are quite a few and Im not sure how they got in.
  # but if diseasename and gene name in the same row occur then keep only one copy,
  linked_diseases  <- linked_diseases[!(duplicated(linked_diseases[c("diseaseName","geneName")]) | duplicated(linked_diseases[c("diseaseName","geneName")], fromLast = TRUE)), ]
  linked_diseases <- linked_diseases[-nrow(linked_diseases),]     # last entry is zero so remove it
  linked_diseases <- anti_join(linked_diseases, disease_umls, by="diseaseName") # If C06 disorders appear , remove them.
  
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

# print tables in LaTex format for inclusion into paper
print_tables <- function(){
  temp_table <- disgene[,c(1,4:6)]
  
  tli.table <- xtable(head(temp_table))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable((digestive[86:101,]))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(head(digestive))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(middle(drug_list))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
 
  tli.table <- xtable(filter(drug_list, MeSH == "C06.405"))
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(filter(gene_list, diseaseId == "C0017178"))
  print(tli.table,floating=FALSE)
  
}

# goanalysis() will enrich a gene with GO terms
# depends on clusterprofiler library and several other things...
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis
go_analysis <- function(yourgenes,ontotype){
  cat("\n",yourgenes)
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ego <- enrichGO(gene          = eg[,2],
                  #universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = ontotype, # one of CC, BP or MF
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  return(ego)
}


# KEGG over-representation test
kegg_analysis <- function(yourgenes){
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  kk <- enrichKEGG(gene         = eg[,2],
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)

  return(kk)
}


# getDiseaseModules() pass it the linkcomm structure, the appropriate data will be
# converted into a dataframe - ready to be passed to GO enrichment functions before
# sending it to GObubble for analysis and display.
getDiseaseModules <- function(linkdata,batch){
  category <- c("MF","BP","CC")
  enrich <- c("GO:0017091", "AU-rich element binding","1/1", "23/16982", "0.001354375","RU12","FU",99)
  
  # remove modules with fewer than 20 genes - as per Menche 2015 paper
  linkdata$clusters <- Filter(function(x)length(x) > 20, linkdata$clusters)
  cat("\nFound ",length(linkdata$clusters), " usable modules.")

  if(is.null(batch)){
    cat("\nERROR: you need to enter ''all'' for all modules or a range in quotes e.g. ''77:89''")
    return(NULL)}
  
  if(batch =="all") {  # process all diseasemodules if batch ==all
    indexstart <- 1
    indexend <- length(linkdata$clusters)
  }else{
    tempy <- unlist(strsplit(batch,":"))  # its a range
    indexstart <- as.numeric(tempy[1])
    indexend <- as.numeric(tempy[2])
  }
  
  #for (i in 1:length(linkdata$clusters)){        # i=num of disease modules
  temp_i <- vector(mode="integer", length=indexend-indexstart);
  z <-1
  
  for (i in indexstart:indexend){        # i=num of disease modules
    items <- getNodesIn(linkdata, clusterids = i)
    temp_i[z] <- i
    z <- z +1
    for (k in 1:length(items)){              # k=num genes in disease module
      for (j in 1:length(category)){  # enrich from MF, BP and CC
        temp_enrich <- go_analysis(items[k],category[j])
        if(!is.null(temp_enrich) && nrow(temp_enrich)>0){
          temp_enrich <- temp_enrich[,c(1:5,8)]
          temp_enrich[,7] <- category[j]       # add the category e.g MF
          temp_enrich[,8] <- i              # add the disease module number (i.e. cluster number)
          enrich <- rbind(enrich,temp_enrich)}
      }
    }
  }
  
  # Fix the dataframe: add z-score, rename V7, add disease module number
  # To match what GOplot i.e. GObubble expects to find
  names(enrich)[names(enrich)=="Description"] <- "term"
  names(enrich)[names(enrich)=="pvalue"] <- "adj_pval"
  names(enrich)[names(enrich)=="V7"] <- "category"
  names(enrich)[names(enrich)=="V8"] <- "DiseaseModule"
  names(enrich)[names(enrich)=="geneID"] <- "genes"
  # convert p-value into "zscore", rename "pvalue" to "adj_pval" 
  namevector <- "zscore"
  enrich[ , namevector] <- qnorm(1 - as.numeric(enrich$adj_pval)/2)
  namevector <- "logFC"
  enrich[ , namevector] <- runif(nrow(enrich), -2, 3)#10#qnorm(1 - as.numeric(enrich$adj_pval)/2)
  namevector <- "count"
  enrich[ , namevector] <- 0  # Number of genes attached to this term.
  enrich$adj_pval <- as.numeric(enrich$adj_pval)
  enrich <- enrich[-1, ]     # 1st entry is rubbish so remove it
 
  # Set "count" for each term
  enrich <- setcount(enrich,temp_i)

  return(enrich)
}

# setcount() gets a count of the terms assigned to each disease module.
setcount <- function(dms,ind){
  #countn <- unique(dms$DiseaseModule) # How many disease modules are there?
  countn <- length(ind) # How many disease modules are there? based on index range?
  cat("\nThere are ",length(ind)," disease modules..numbered from",ind)
  for (j in 1:length(ind)){
    cat("\nJ is now",j)
    temp_dms <- filter(dms,DiseaseModule == (ind[j]))
    nterm <- unique(temp_dms$term) # How many unique terms do we have for this disease module?
    for (k in 1:length(nterm)){
      tcount <- nrow(filter(dms,term == nterm[k]))
      dms$count[dms$term == nterm[k] & dms$DiseaseModule == (ind[j])] <- tcount
      #cat("\nFor DM",ind[j]," tcount is ",tcount)
    }
  }
  
  return(dms)
} 

# This is a far quicker version of getDiseaseModules(). This version uses the ontologySimilarity packages
# by Daniel Green. NB for the moment cannot do bubbleplot as we cannot as yet generate p-values or zscores.
createDiseaseModules <- function(linkdata){
  tempgenes <- names(gene_GO_terms)
  enrich <- data.frame(ID="GO:0000666", genes="RU12",DiseaseModule=666,adj_pval=0.001,zscore=6.001,
                       category="FU",term="Satanic like behaviour",stringsAsFactors=FALSE) #instantiate.
  # remove modules with fewer than 20 genes - as per Menche 2015 paper
  newclusters <- Filter(function(x)length(x) > 20, linkdata$clusters)
  cat("\nFound ",length(newclusters), " usable modules.")
  j<- 0  # set counter for clusters bigger than 20 
  
  for (i in 1:length(linkdata$clusters)){
    tempnodes <- getNodesIn(linkdata, clusterids = i,type="names")
    if(length(tempnodes) >= 20){
      j <- j+1
      #cat("\nj =...",j)
      tempnodes <- tempnodes[tempnodes %in% tempgenes]  # remove genes that do not exist in GO data 
      tempgo <- gene_GO_terms[tempnodes]
      cc <- go$id[go$name == "cellular_component"]
      bp <- go$id[go$name == "biological_process"]
      mf <- go$id[go$name == "molecular_function"] 
      temp_cc <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=cc, x))
      temp_bp <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=bp, x))
      temp_mf <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=mf, x))
      tmp_enrich_c <- data.frame(unlist(temp_cc),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...CC")
      tmp_enrich_c <- cbind(tmp_enrich_c,rownames(tmp_enrich_c),stringsAsFactors=FALSE) # gene names
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(j,length(unlist(temp_cc))))               # diseasemodule number
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(0.01,length(unlist(temp_cc))))            # adj_pval
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(3.2,length(unlist(temp_cc))))             # zscore
      tmp_enrich_c <- cbind(tmp_enrich_c,rep("CC",length(unlist(temp_cc))),stringsAsFactors=FALSE) # category
      tmp_enrich_c <- cbind(tmp_enrich_c,unname(go$name[unlist(temp_cc)]),stringsAsFactors=FALSE)
      
      tmp_enrich_b <- data.frame(unlist(temp_bp),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...BP")
      tmp_enrich_b <- cbind(tmp_enrich_b,rownames(tmp_enrich_b),stringsAsFactors=FALSE) # gene names
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(j,length(unlist(temp_bp))))               # diseasemodule number
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(0.01,length(unlist(temp_bp))))            # adj_pval
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(3.2,length(unlist(temp_bp))))             # zscore
      tmp_enrich_b <- cbind(tmp_enrich_b,rep("BP",length(unlist(temp_bp))),stringsAsFactors=FALSE) # category
      tmp_enrich_b <- cbind(tmp_enrich_b,unname(go$name[unlist(temp_bp)]),stringsAsFactors=FALSE)
      
      tmp_enrich_m <- data.frame(unlist(temp_mf),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...MF")
      tmp_enrich_m <- cbind(tmp_enrich_m,rownames(tmp_enrich_m),stringsAsFactors=FALSE) # gene names
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(j,length(unlist(temp_mf))))               # diseasemodule number
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(0.01,length(unlist(temp_mf))))            # adj_pval
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(3.2,length(unlist(temp_mf))))             # zscore
      tmp_enrich_m <- cbind(tmp_enrich_m,rep("MF",length(unlist(temp_mf))),stringsAsFactors=FALSE) # category
      tmp_enrich_m <- cbind(tmp_enrich_m,unname(go$name[unlist(temp_mf)]),stringsAsFactors=FALSE)
      
      colnames(tmp_enrich_b)[1] <- "ID"; colnames(tmp_enrich_b)[2] <- "genes"; colnames(tmp_enrich_b)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_b)[4] <- "adj_pval"; colnames(tmp_enrich_b)[5] <- "zscore"; colnames(tmp_enrich_b)[6] <- "category" 
      colnames(tmp_enrich_b)[7] <- "term" 
      rownames(tmp_enrich_b) <- c()
      
      colnames(tmp_enrich_c)[1] <- "ID"; colnames(tmp_enrich_c)[2] <- "genes"; colnames(tmp_enrich_c)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_c)[4] <- "adj_pval"; colnames(tmp_enrich_c)[5] <- "zscore"; colnames(tmp_enrich_c)[6] <- "category" 
      colnames(tmp_enrich_c)[7] <- "term" 
      rownames(tmp_enrich_c) <- c()
      
      colnames(tmp_enrich_m)[1] <- "ID"; colnames(tmp_enrich_m)[2] <- "genes"; colnames(tmp_enrich_m)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_m)[4] <- "adj_pval"; colnames(tmp_enrich_m)[5] <- "zscore"; colnames(tmp_enrich_m)[6] <- "category" 
      colnames(tmp_enrich_m)[7] <- "term" 
      rownames(tmp_enrich_m) <- c()
      
      #cat("\nRBIND....CC , BF & MF to enrich")
      enrich <- rbind(enrich,tmp_enrich_c)
      enrich <- rbind(enrich,tmp_enrich_b)
    }
  }
  enrich <- enrich[-1, ]     # 1st entry is rubbish so remove it
  tmp_count <- table(enrich$ID)  # obtain a count of how many times each GO ID appears.
  idcount <- vector(mode="integer",length=nrow(enrich))
  
  idcount <- as.vector(tmp_count)[match(enrich$ID, names(tmp_count))] # almost too clever
  enrich <- cbind(enrich,idcount)
  colnames(enrich)[8] <- "count"   

  return(enrich)
}


# getdrugs() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drugs <- function(umls,rlist) {
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


# print_dm_table() will generate the latex stuff based on annoations and ranking 
# methods to create the disease module table for the paper, containing:
#   C06/DX0 numbers
#   GO enrichment counts
#   KEGG enrichment counts
#   Biological plausibility ranking
#   Current drugs / DX0 drug reposition candidates
#   components/complexity

print_dm_table <- function(){

  dm.table <- xtable(head(temp_table))
  #digits(tli.table)[c(2,6)] <- 0
  print(dm.table,floating=FALSE)
  
}











