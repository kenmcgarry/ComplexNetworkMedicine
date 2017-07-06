## gastrointestinal_hdn.R     6/12/16
## load in files for Gastrointestinal disease (hdn=human disease network) and gene lists
##                 source("http://bioconductor.org/biocLite.R")
##                biocLite("GOstats")

## NOTE: tree structure is generated from parstree.R file.

library(igraph)
library(linkcomm)
library(bipartite)

library(dplyr)
library(tidyr)
library(xtable)
library(ggplot2)

library(XLConnect) # for reading and writing excel spreadsheets
library(rentrez)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)

# from DisGeneNet site # http://www.disgenet.org/web/DisGeNET/menu/home
mappings <- file.path('C://R-files//disease//','map_umls_2_mesh.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
disgene <- file.path('C://R-files//disease//','curated_gene_disease_associations.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
disgene <- file.path('C://R-files//disease//','all_gene_disease_associations.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")

# https://www.datacamp.com/community/tutorials/r-tutorial-read-excel-into-r#gs.ax73e9Q
wb <- loadWorkbook("C://R-files//disease//stomachcancer_genealacart6-12-16.xlsx")
df <- readWorksheet(wb, sheet=3)

string <- file.path('C://R-files//disease//','string_interactions.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
string <- data.frame(lapply(string, as.character), stringsAsFactors=FALSE)

#----------------- use OMIM database (via rentrez) to generate list of implicated disease genes for each disease -----------
entrez_db_searchable("omim")
Dproteins <- c("TNF","HRH2","APC","MLN","CDKN2A")   # crud<-c("7124","3274","324","1029","9318","79026","1654","65003","6240","3476","6238","3836")

x <- entrez_search(db="omim", term="126850[GENE]",retmax=2000)

gene <- c("11171", "8243", "112464", "2194",
          "9318", "79026", "1654", "65003",
          "6240", "3476", "6238", "3836",
          "4176", "1017", "249",
          "7124","3274","324","1029")

genelist <- bitr(Dproteins,fromType = "SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
genelist <- genelist[2]
genelist <- sapply(genelist, as.character)
genelist <- as.character(unname(genelist, force = TRUE))
genelist<- as.character(genelist)
yy = enrichPathway(genelist, pvalueCutoff=0.1,organism="human",readable=TRUE)
head(yy)

# ============= FUNCTIONS HERE ==================================
# Calculate some statistics about the disease gene network
get_gstatistics <- function(gt) {
  gstats <- data.frame(
    modularity=modularity(gt, membership(cluster_walktrap(gt))),
    avepath=average.path.length(gt),
    nedges=ecount(gt),
    nverts=vcount(gt),
    transit=transitivity(gt),
    avedegree=mean(degree(gt)),
    diameter=diameter(gt,weights=NA),
    connect=is.connected(gt),
    closeness=closeness(gt),
    betweenness=betweenness(gt,directed=FALSE),
    density=graph.density(gt),
    hubness=hub_score(gt)$vector,
    authority=authority.score(gt)$vector)
  #power=bonpow(gt))
  return(gstats)
}

# Count how many interactions each protein has.
count_interactions <- function(CVDP) {
  
  for (i in 1:length(CVDP)){
    
    print(CVDP[i])
    pname <- paste(CVDP[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids[1])
    plist<-unique(plist[13])
    
    ptemp <- cbind(CVDP[i],nrow(plist))
    
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  
  return(ppi)
}

# For each disease protein get the proteins they interact with. Uses NCBI2R package.
get_interactions <- function(CVDP){
  
  for (i in 1:length(CVDP)){
    print(CVDP[i])
    pname <- paste(CVDP[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids)
    
    plist<-unique(plist[13])
    #print(plist)
    cvd<-rep(CVDP[i],nrow(plist))
    ptemp <- cbind(cvd,plist)
    
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  
  return(ppi)
}



# See how many research articles are written about our proteins. Uses rentrez package.
count_articles <- function (CVDP){
  for (i in 1:length(CVDP)){
    
    print(CVDP[i])
    pname <- paste(CVDP[i],'[GENE]) AND (Homo sapiens[ORGN])',sep="")
    ids<-entrez_search(db="pubmed", term=pname,retmax=40000)
    atemp <- cbind(CVDP[i],length(ids$ids))
    
    if(i!=1){
      articles <- rbind(articles,atemp)} 
    else{
      articles <- atemp}
  }
  return(articles)
}

# scrap the interaction table from the NCBI web page.
get.ppiNCBI <- function(g.n) {
  require(XML)
  ppi <- data.frame()
  for(i in 1:length(g.n)){
    o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
    # check if interaction table exists
    exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
    if(exist){
      p <- getNodeSet(o, "//table")
      ## need to know which table is the good one
      for(j in 1:length(p)){
        int <- readHTMLTable(p[[j]])
        if(colnames(int)[2]=="Interactant"){break}
      }
      ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
    }
    # play nice! and avoid being kicked out from NCBI servers
    Sys.sleep(1)
  }
  if(dim(ppi)[1]>0){
    ppi <- unique(ppi)
    print(paste(dim(ppi)[1], "interactions found"))
    return(ppi)
  } else{
    print("No interaction found")
  }
}


