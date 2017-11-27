# get_drugs.R
# Ken McGarry updated: 23/11/17

# packages are loaded in by gi_functions.R

setwd("C:/R-files/disease")    # point to where my code lives
load("C06disease-27thNov2017.RData") # load in required data - the contents will change regulary
source("gi_functions.R")  # load in the functions required for finding lists of drugs 
source("gi_run.R")   # some routine code to load in.

# Work with following dataframes: indications; mappings; digestive; disgene; 
# save(disease_umls,indications,restrictedlist,mappings,digestive,disgene,gene_list,drug_list,simMA,fpdrugs,drugnames,drugids, file = "C06disease-17thNov2017.RData")
# the key to linking these files is the "ID" in digestive also called "meshId" in mappings

# ----- Stages:  ----------
# 1. Use "digestive" dataframe to parse through list of digestive diseases using the ID.
# 2. We need to get the drugs associated with each disease.
# 3. Get genes associated (if known) with each disease from "disgene" dataframe.
# 4. Get chemical structures of drugs and create fingerprints.
# 5. Get 2nd shell protein & drug interactions
# 6. Use 2nd shell stuff to investigate linkages with other diseases

# Code below perfoms stages 1 and 2 ----------
disease_umls <- id2umls(digestive)  # get umls codes for each disease
drug_list <- data.frame(drugbank_id=character(),drugbank_name=character(),umls_cui_from_meddra=character(),meddra_name=character())

for (i in 1:nrow(disease_umls)){
  drug_temp <- get_drugs_plus(disease_umls[i,1],restrictedlist) # get drugs used for each disease
  if (!is.null(drug_temp)){
    drug_list <- rbind(drug_temp,drug_list)} # add next drug to end of "drug_list" dataframe
}

drug_list <- add_meshcode(drug_list) # adds the old mesh code and unique id to drug_list, useful to relate drugs to diseases.
unique(drug_list$drugbank_name) # how many unique drugs do we have?

# Code for stage 3.--------------------------
gene_list <- get_disease_genes(drug_list)  # not all diseases will have implicated genes, so list will shorter!
                                                                          # 55 in fact

# Code for stage 4. ------------------------
# See drugstructure_gi.R for details of processing........
# generated chemical FingerPrints of 1024 bits for 189 drugs using SDF data
# generated dendrogram plot and heatmap plot for the 189 drugs.
# dataframes are: simMA; drugids; drugnames; fpdrugs; 
# need the following packages: library(ChemmineR); library(ChemmineOB); library(ape); library(sparcl)
# library(cluster); library(xtable); library(gplots); library(scatterplot3d) 

# THE ABOVE CODE IS ONLY PART COMPLETED AS WE NEED TO EXAMINE 2ND SHELL DRUGS STRUCTURE AND COMPARE


# Code for stage 5. ------------------------
# Collect 2nd shell proteins known to interact with a. our drugs b. our proteins and c. targets of the drugs 
# This will provide the basis for linkages to other diseases

# Make a text file to be uploaded to STITCH webpage
# 892 unique genes, sorted alphbetically.Removed quotation marks manually as it screws up STITCH
write.table(unique(sort(drugnames)),"C:\\R-files\\disease\\C06drugs.txt",sep=",",row.names = FALSE,col.names = FALSE)
write.table(unique(sort(gene_list$geneName)),"C:\\R-files\\disease\\C06genes.txt",sep=",",row.names = FALSE,col.names = FALSE)

# STITCH can only find 820 out of 892 C06 genes connected to 907 2nd shell genes
# with 15,736 interactions.

# load in the mesh data, keep only first three variables.
temptree <- file.path('C://R-files//disease//','meshtreefull.csv') %>% read.delim(na.strings='',sep=',',header=TRUE,comment.char="#")
temptree <- temptree[,1:3]

# load in 2nd shell genes from STITCH search
#shell2_genes <- read.csv("C:\\R-files\\disease\\C06-shell2-low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
shell2_genes <- read.csv("C:\\R-files\\disease\\shell2_genes_low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false

# Keep interacting gene columns and confidence score, discard the other 13 variables.
shell2_genes <- shell2_genes[,c("X.node1","node2","combined_score")]
shell2_genes <- c(shell2_genes$X.node1,shell2_genes$node2)
shell2_genes <- unique(shell2_genes)  
shell2_genes <- setdiff(shell2_genes,gene_list$geneName) # we are left with 109 genes unique to shell2 and not related to C06
write.table(unique(sort(shell2_genes)),"C:\\R-files\\disease\\shell2genes.txt",sep=",",row.names = FALSE,col.names = FALSE)

# Now seek out the two shell levels of diseases
shell1 <- get_linked_diseases(gene_list$geneName)  # diseases directly linked to C06 disease genes

shell2 <- get_linked_diseases(shell2_genes)  # diseases indirectly linked through 2nd shell genes

         #get_linked_diseases("SMAD2")  # diseases indirectly linked through 2nd shell genes

# Load in drug interactions, majority are drug-2-drug interactions with a few genes thrown in.
drug_interactions <- read.csv("C:\\R-files\\disease\\drug_interactions.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
drug_interactions <- drug_interactions[,1:2]

# load 1st shell interactions
shell1_interactions <- read.csv("C:\\R-files\\disease\\C06-shell1-low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
shell1_interactions <- shell1_interactions[,1:2]

# load 1st shell interactions
shell2_interactions <- read.csv("C:\\R-files\\disease\\shell2_genes_low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
shell2_interactions <- shell2_interactions[,1:2]


# get_all_linked_diseases("BRCA1") # test out new function

##################################################################
# Using LINKCOMM to detect disease modules
# 1. drug modules
# 2. 1st shell gene modules
# 3. 2nd shell gene modules

d1 <- getLinkCommunities(drug_interactions, hcmethod = "single")
s1 <- getLinkCommunities(shell1_interactions, hcmethod = "single")
s2 <- getLinkCommunities(shell2_interactions, hcmethod = "single")

d1mods <- getDiseaseModules(d1)
s1mods <- getDiseaseModules(s1)
s2mods <- getDiseaseModules(s2)




print(s2)
head(s2$numclusters)
plot(s2, type = "graph", shownodesin = 2, node.pies = TRUE)
plot(s2, type = "summary")
oc <- getOCG.clusters(shell2_interactions) 
plot(oc, type = "graph", shownodesin = 7, scale.vertices = 0.1)
cent <- getCommunityCentrality(s2)
plot(cent)

plot_pdens(d1$pdens)
plot_pdens(s1$pdens)
plot_pdens(s2$pdens)

commod <- getCommunityConnectedness(d1, conn = "modularity")
comcon<- getCommunityConnectedness(d1, conn = "conn")
plot_com(commod,comcon)

commod <- getCommunityConnectedness(s1, conn = "modularity")
comcon<- getCommunityConnectedness(s1, conn = "conn")
plot_com(commod,comcon)

commod <- getCommunityConnectedness(s2, conn = "modularity")
comcon<- getCommunityConnectedness(s2, conn = "conn")
plot_com(commod,comcon)

cent <- getCommunityCentrality(d1)
plot_centrality(cent,"drug index")

cent <- getCommunityCentrality(s1)
plot_centrality(cent,"gene index")

cent <- getCommunityCentrality(s2)
plot_centrality(cent,"gene index")


# matrix plot
plotLinkCommMembers(s1, nodes = head(names(lc$numclusters), 20),
                    pal = brewer.pal(11, "Spectral"), shape = "rect", total = TRUE,
                    fontsize = 11, nspace = 3.5, maxclusters = 20)


lc <- getLinkCommunities(drug_interactions, hcmethod = "single")
plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", shownodesin = 2, node.pies = TRUE)
plot(lc, type = "summary")
plot(lc, type = "con")


cr <- getClusterRelatedness(lc, hcmethod = "ward.D")
mc <- meta.communities(lc, hcmethod = "ward.D", deepSplit = 0)
cent <- getCommunityCentrality(s2)
cm <- getCommunityConnectedness(lc, conn = "modularity")

head(sort(cc, decreasing = TRUE))

plot(lc, type = "commsumm", summary = "modularity")
plot(lc, type = "commsumm", summary = "con")


# matrix plot
plotLinkCommMembers(lc, nodes = head(names(lc$numclusters), 20),
                    pal = brewer.pal(11, "Spectral"), shape = "rect", total = TRUE,
                    fontsize = 11, nspace = 3.5, maxclusters = 20)



 

##################################################################
## GO and KEGG enrichment
goa <- go_analysis(shell2_genes)
barplot(goa, drop=TRUE, showCategory=20)

kega <- kegg_analysis(shell2_genes)
barplot(kega, drop=TRUE, showCategory=20)

# Tables for paper. 
tempgoa <- head(goa,row.names=FALSE)
tempgoa <- tempgoa[,1:5]
print(xtable(tempgoa, display=c("s","s","s","s","s","g")), math.style.exponents = TRUE,include.rownames = FALSE)

tempkega <- tail(kega,row.names=FALSE)
tempkega <- tempkega[,1:5]
print(xtable(tempkega, display=c("s","s","s","s","s","g")), math.style.exponents = TRUE,include.rownames = FALSE)

dismods <- getDiseaseModules(s2)
#dismods <- enrich
dismods <- dplyr::select(enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
GOBubble(dismods, labels = 3)
GOBubble(circ, labels = 3)
#rm(goa,kega,tli.table,EC)










