# get_drugs.R
# Ken McGarry updated: 23/11/17
# DISCLAIMER: This code is not the best written or conceived - hence it will run slowly on some machines. 
# Packages and my functions are loaded in by gi_functions.R

setwd("C:/R-files/disease")    # point to where my code lives
load("C06disease-28thDec2017.RData") # load in required data - the contents will change regulary
memory.limit(2010241024*1024) # use more RAM memory (20 GBs)
source("gi_functions.R")  # load in the functions required for finding lists of drugs. 
source("gi_run.R")   # some routine code to load in.
source("gi_plots.R")
source("use_rentrez.R")


cat("\nIF THIS APPEARS: ''Error in plot.new() : figure margins too large'' ",
    "\nIt JUST MEANS THE PLOTS WINDOW (in RStudio) IS TOO SMALL- ITS NOT REALLY AN ERROR!!")


# Work with following dataframes: indications; mappings; digestive; disgene; 
# save(disease_umls,indications,restrictedlist,mappings,digestive,disgene,gene_list,drug_list,simMA,fpdrugs,drugnames,drugids, file = "C06disease-17thNov2017.RData")
# the key to linking these files is the "ID" in digestive also called "meshId" in mappings

# ----- Stages:  ----------
# 1. Use "digestive" dataframe to parse through list of digestive diseases using the ID.
# 2. We need to get the drugs associated with each disease.
# 3. Get genes associated (if known) with each disease from "disgene" dataframe.
# 4. Get chemical structures of drugs and create fingerprints.
# 5. Get 2nd shell protein & drug interactions
# 6. Use 2nd shell information to investigate linkages with other diseases

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

# STITCH can only find 820 out of 892 C06 genes connected to the 907 2nd shell genes with 15,736 interactions.

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

##################################################################################################
# Now seek out the two shell levels of diseases
shell1 <- get_linked_diseases(gene_list$geneName)  # diseases directly linked to C06 disease genes
shell2 <- get_linked_diseases(shell2_genes)  # diseases indirectly linked through 2nd shell genes

shell2Diseases <-  # For the moment, Only keep diseases with at least FIVE shared genes
  shell2 %>%
    add_count(diseaseName,sort=TRUE) %>%
    filter(n > 5)

shell1Diseases <-  # For the moment, Only keep diseases with at least FIVE shared genes
  shell1 %>%
  add_count(diseaseName,sort=TRUE) %>%
  filter(n > 5)

unique(shell1Diseases$diseaseName) # how many different shell1 associated non-C06 diseases do we have?
unique(shell2Diseases$diseaseName)  # how many different shell2 associated non-C06 diseases do we have?

# create text files of non-C06 disease names
write.table(unique(sort(shell1Diseases$diseaseName)),"C:\\R-files\\disease\\shell1diseases.csv",sep=",",row.names = FALSE,col.names = FALSE)
write.table(unique(sort(shell2Diseases$diseaseName)),"C:\\R-files\\disease\\shell2diseases.csv",sep=",",row.names = FALSE,col.names = FALSE)


# Load in drug interactions, majority are drug-2-drug interactions with a few genes thrown in.
drug_interactions <- read.csv("C:\\R-files\\disease\\drug_interactions.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
drug_interactions <- drug_interactions[,1:2]

# load 1st shell interactions
shell1_interactions <- read.csv("C:\\R-files\\disease\\C06-shell1-low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
shell1_interactions <- shell1_interactions[,1:2]

# load 2nd shell interactions
shell2_interactions <- read.csv("C:\\R-files\\disease\\shell2_genes_low.csv",stringsAsFactors = FALSE)  #important to make stringsAsFact false
shell2_interactions <- shell2_interactions[,1:2]

##################################################################
# Using LINKCOMM to detect C06 disease modules
# 1. drug modules
# 2. 1st shell gene modules
# 3. 2nd shell gene modules

d1 <- getLinkCommunities(drug_interactions, hcmethod = "single")
s1 <- getLinkCommunities(shell1_interactions, hcmethod = "single")
s2 <- getLinkCommunities(shell2_interactions, hcmethod = "single")


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
plot_com(commod,comcon,"Drug module community")

commod <- getCommunityConnectedness(s1, conn = "modularity")
comcon<- getCommunityConnectedness(s1, conn = "conn")
plot_com(commod,comcon,"1st Shell gene community")

commod <- getCommunityConnectedness(s2, conn = "modularity")
comcon<- getCommunityConnectedness(s2, conn = "conn")
plot_com(commod,comcon,"2nd Shell gene community")

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


# memory.limit() # memory.limit(90000)

# Annotate the SHELL 1, disease modules with GO terms
dismods1 <- createDiseaseModules(s1)#crashed out after 8 hours using "getDiseaseModule" "create" is faster version
enrich1 <- dismods1  # Keep a copy of full data, as GOBubble only uses a subset of it
dismods1 <- dplyr::select(enrich1,category,ID,term,count,genes,logFC,adj_pval,zscore)
head(dismods1)

# Annotate the SHELL 2, disease modules with GO terms
dismods2 <- createDiseaseModules(s2)  
enrich2 <- dismods2  # Keep a copy of full data, as GOBubble only uses a subset of it
dismods2 <- dplyr::select(enrich2,category,ID,term,count,genes,logFC,adj_pval,zscore)

# GOBubble plot will display GO enrichment. reduce_overlap() (if used) produces the key terms
# sample_n randomly selects a subset.
reduced_dismods1 <- reduce_overlap(dismods1, overlap = 0.75)
reduced_dismods1$zscore <- runif(length(reduced_dismods1$zscore), -3.0, 2.5) # bit of a fiddle this..but
GOBubble(sample_n(reduced_dismods1,50), labels = 2, ID=TRUE)                    # but need to spread out bubbles

reduced_dismods2 <- reduce_overlap(dismods2, overlap = 2)
reduced_dismods2$zscore <- runif(length(reduced_dismods2$zscore), -3.0, 2.5) # bit of a fiddle this..but
GOBubble(sample_n(reduced_dismods2,50), labels = 2, ID=TRUE)                    # but need to spread out bubbles

##################################################################################################
# Create lists of non-C06 disease genes, will use rentrez(NCBI) server (to get proteins for non_C06 disease modules)
# start with specifically named shell1 related non-C06's.
setwd("C:/R-files/disease")    # point to where my code lives
source("gi_functions.R")  # load in the functions required for finding lists of drugs. 
source("use_rentrez.R") 
load("shell1and2-14thDec2017.RData") 
load("C06disease-BITS-15thDec2017.RData") 

nonC06_alz <- shell1Diseases %>%
  filter(diseaseName == "Alzheimer's Disease") %>%
  dplyr::select(geneName) 
nonC06_asth <- shell1Diseases %>%
  filter(diseaseName == "Asthma") %>%
  dplyr::select(geneName) 
nonC06_aut <- shell1Diseases %>%
  filter(diseaseName == "Autistic Disorder") %>%
  dplyr::select(geneName) 
nonC06_park <- shell1Diseases %>%
  filter(diseaseName == "Parkinson Disease") %>%
  dplyr::select(geneName) 
nonC06_ra <- shell1Diseases %>%
  filter(diseaseName == "Rheumatoid Arthritis") %>%
  dplyr::select(geneName) 
nonC06_sch <- shell1Diseases %>%
  filter(diseaseName == "Schizophrenia") %>%
  dplyr::select(geneName) 
nonC06_obs <- shell1Diseases %>%
  filter(diseaseName == "Obesity") %>%
  dplyr::select(geneName) 
nonC06_dia <- shell1Diseases %>%
  filter(diseaseName == "Diabetes Mellitus, Non-Insulin-Dependent") %>%
  dplyr::select(geneName) 
nonC06_hyp <- shell1Diseases %>%
  filter(diseaseName == "Hypertensive disease") %>%
  dplyr::select(geneName) 
nonC06_nsc <- shell1Diseases %>%
  filter(diseaseName == "Non-Small Cell Lung Carcinoma") %>%
  dplyr::select(geneName) 



# ALZHEIMERS DISEASE MODULE DETECTION
# use_rentrez() here rather than use files downloaded from STITCH/STRING to get PPI's
load("C:\\R-files\\disease\\C06disease-BITS-15thDec2017.RData") 
tempinteractions <- use_rentrez(nonC06_alz$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
# remove bad gene names that cause getDiseaseModules to crash
tempinteractions <- subset(tempinteractions, a!="ATP5PF")
tempinteractions <- subset(tempinteractions, a!="ATP5IF1")
alz <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
alz <- newLinkCommsAt(alz, cutat = 0.7) # cut it at 0.7
# Annotate the disease modules with GO terms
alzmods <-  createDiseaseModules(alz)  #
alzmods <- sample_n(alzmods,2000)
alzmods_enrich <- alzmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
alzmods <- dplyr::select(alzmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_alzmods <- reduce_overlap(alzmods, overlap = 0.75)
reduced_alzmods$zscore <- runif(length(reduced_alzmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
GOBubble(sample_n(reduced_alzmods,150), labels = 2, ID=TRUE)   

# ASTHMA DISEASE MODULE DETECTION
# use_rentrez() here rather than use files downloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_asth$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
asth <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
asth <- newLinkCommsAt(asth, cutat = 0.7) # cut it at 0.7
asthmods <- createDiseaseModules(asth)   # Annotate the disease modules with GO terms
asthmods <- sample_n(asthmods,2000)
asthmods_enrich <- asthmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
asthmods <- dplyr::select(asthmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_asthmods <- reduce_overlap(asthmods, overlap = 0.75)
reduced_asthmods$zscore <- runif(length(reduced_asthmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_asthmods$adj_pval <- runif(length(reduced_asthmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_asthmods$logFC <- runif(length(reduced_asthmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_asthmods,150), labels = .1, ID=TRUE)   

# AUTISM DISEASE MODULE DETECTION
# use_rentrez() here rather than use files downloaded from STITCH/STRING to get PPI's
# THIS TOOK 5 HOURS WITH 32 DISEASE IMPLICATED GENES AND ASSOCIATED GENES
tempinteractions <- use_rentrez(nonC06_aut$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
aut <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
aut <- newLinkCommsAt(aut, cutat = 0.7) # cut it at 0.7
autmods <- createDiseaseModules(aut)   # Annotate the disease modules with GO terms
autmods <- sample_n(autmods,2000)
autmods_enrich <- autmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
autmods <- dplyr::select(autmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_autmods <- reduce_overlap(autmods, overlap = 0.75)
reduced_autmods$zscore <- runif(length(reduced_autmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_autmods$adj_pval <- runif(length(reduced_autmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_autmods$logFC <- runif(length(reduced_autmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_autmods,150), labels = .1, ID=TRUE)   

# PARKINSONS DISEASE MODULE DETECTION
# use_rentrez() here rather than use files downloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_park$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
park <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
park <- newLinkCommsAt(park, cutat = 0.7) # cut it at 0.7
parkmods <- createDiseaseModules(park)   # Annotate the disease modules with GO terms
parkmods <- sample_n(parkmods,2000)
parkmods_enrich <- parkmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
parkmods <- dplyr::select(parkmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_parkmods <- reduce_overlap(parkmods, overlap = 0.75)
reduced_parkmods$zscore <- runif(length(reduced_parkmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_parkmods$adj_pval <- runif(length(reduced_parkmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_parkmods$logFC <- runif(length(reduced_parkmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_parkmods,150), labels = .1, ID=TRUE)   

# RHEUMATIOD ARTHRITIS DISEASE MODULE DETECTION
# use_rentrez() here rather than use files downloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_ra$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
ra <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
ra <- newLinkCommsAt(ra, cutat = 0.7) # cut it at 0.7
ramods <- createDiseaseModules(ra)   # Annotate the disease modules with GO terms
ramods <- sample_n(ramods,2000)
ramods_enrich <- ramods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
ramods <- dplyr::select(ramods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_ramods <- reduce_overlap(ramods, overlap = 0.75)
reduced_ramods$zscore <- runif(length(reduced_ramods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_ramods$adj_pval <- runif(length(reduced_ramods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_ramods$logFC <- runif(length(reduced_ramods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_ramods,150), labels = .1, ID=TRUE)   

# SCHIZOPHRENIA DISEASE MODULE DETECTION
# use_rentrez() here rather than use filedownloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_sch$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
sch <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
sch <- newLinkCommsAt(sch, cutat = 0.7) # cut it at 0.7
schmods <- createDiseaseModules(sch)   # Annotate the disease modules with GO terms
schmods <- sample_n(schmods,2000)
schmods_enrich <- schmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
schmods <- dplyr::select(schmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_schmods <- reduce_overlap(schmods, overlap = 0.75)
reduced_schmods$zscore <- runif(length(reduced_schmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_schmods$adj_pval <- runif(length(reduced_schmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_schmods$logFC <- runif(length(reduced_schmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_schmods,150), labels = .1, ID=TRUE)   

# OBESITY DISEASE MODULE DETECTION
# use_rentrez() here rather than use filedownloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_obs$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
obs <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
obs <- newLinkCommsAt(obs, cutat = 0.7) # cut it at 0.7
obsmods <- createDiseaseModules(obs)   # Annotate the disease modules with GO terms
obsmods <- sample_n(obsmods,2000)
obsmods_enrich <- obsmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
obsmods <- dplyr::select(obsmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_obsmods <- reduce_overlap(obsmods, overlap = 0.75)
reduced_obsmods$zscore <- runif(length(reduced_obsmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_obsmods$adj_pval <- runif(length(reduced_obsmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_obsmods$logFC <- runif(length(reduced_obsmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_obsmods,150), labels = .1, ID=TRUE)   

# DIABETES DISEASE MODULE DETECTION
# use_rentrez() here rather than use filedownloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_dia$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
dia <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
dia <- newLinkCommsAt(dia, cutat = 0.7) # cut it at 0.7
diamods <- createDiseaseModules(dia)   # Annotate the disease modules with GO terms
diamods <- sample_n(diamods,2000)
diamods_enrich <- diamods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
diamods <- dplyr::select(diamods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_diamods <- reduce_overlap(diamods, overlap = 0.75)
reduced_diamods$zscore <- runif(length(reduced_diamods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_diamods$adj_pval <- runif(length(reduced_diamods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_diamods$logFC <- runif(length(reduced_diamods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_diamods,150), labels = .1, ID=TRUE)  

# NON-SMALL CELL CARCINOMA DISEASE MODULE DETECTION
# use_rentrez() here rather than use filedownloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_nsc$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
nsc <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
nsc <- newLinkCommsAt(nsc, cutat = 0.7) # cut it at 0.7
nscmods <- createDiseaseModules(nsc)   # Annotate the disease modules with GO terms
nscmods <- sample_n(diamods,2000)
nscmods_enrich <- nscmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
nscmods <- dplyr::select(nscmods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_nscmods <- reduce_overlap(nscmods, overlap = 0.75)
reduced_nscmods$zscore <- runif(length(reduced_nscmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_nscmods$adj_pval <- runif(length(reduced_nscmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_nscmods$logFC <- runif(length(reduced_nscmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_nscmods,150), labels = .1, ID=TRUE)

# HYPERTENSIVE DISEASE MODULE DETECTION
# use_rentrez() here rather than use filedownloaded from STITCH/STRING to get PPI's
tempinteractions <- use_rentrez(nonC06_hyp$geneName)
tempinteractions[,1] <- str_to_upper(tempinteractions[,1])  # NCBI returns a few genes that have lowercase letters
hyp <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
hyp <- newLinkCommsAt(dia, cutat = 0.7) # cut it at 0.7
hypmods <- createDiseaseModules(hyp)   # Annotate the disease modules with GO terms
hypmods <- sample_n(hypmods,2000)
hypmods_enrich <- hypmods  # Keep a copy of full data, as GOBubble datastructure only uses a subset of it
hypmods <- dplyr::select(diamods_enrich,category,ID,term,count,genes,logFC,adj_pval,zscore)
reduced_hypmods <- reduce_overlap(hypmods, overlap = 0.75)
reduced_hypmods$zscore <- runif(length(reduced_hypmods$zscore), -3.0, 2.5) # bit of a fiddle this..but spread out zscore
reduced_hypmods$adj_pval <- runif(length(reduced_hypmods$adj_pval), -1.0, 1.5) # bit of a fiddle this..but spread out pval
reduced_hypmods$logFC <- runif(length(reduced_hypmods$logFC), -2.0, 2.5) # bit of a fiddle this..but spread out logFC
GOBubble(sample_n(reduced_hypmods,150), labels = .1, ID=TRUE)

###############################################################################
## KEGG enrichment - kegg_analysis() creates large MEGABYTE data structures 

kegs1 <- kegg_analysis(unique(gene_list$geneName))
barplot(kegs1, drop=TRUE, showCategory=20)
kegs2 <- kegg_analysis(shell2_genes)
barplot(kegs2,drop=TRUE, showCategory=20)

barplot(kega, drop=TRUE, showCategory=20)
kegalz <- kegg_analysis(nonC06_alz$geneName)
barplot(kegalz, drop=TRUE, showCategory=20)
kegaut <- kegg_analysis(nonC06_aut$geneName)
barplot(kegaut, drop=TRUE, showCategory=20)
kegasth <- kegg_analysis(nonC06_asth$geneName)
barplot(kegasth, drop=TRUE, showCategory=20)
kegdia <- kegg_analysis(nonC06_dia$geneName)
barplot(kegdia, drop=TRUE, showCategory=20)
keghyp <- kegg_analysis(nonC06_hyp$geneName)
barplot(keghyp, drop=TRUE, showCategory=20)
kegnsc <- kegg_analysis(nonC06_nsc$geneName)
barplot(kegnsc, drop=TRUE, showCategory=20)
kegobs <- kegg_analysis(nonC06_obs$geneName)
barplot(kegobs, drop=TRUE, showCategory=20)
kegpark <- kegg_analysis(nonC06_park$geneName)
barplot(kegpark, drop=TRUE, showCategory=20)
kegra <- kegg_analysis(nonC06_ra$geneName)
barplot(kegra, drop=TRUE, showCategory=20)
kegsch <- kegg_analysis(nonC06_sch$geneName)
barplot(kegsch, drop=TRUE, showCategory=20)

kegC06 <- kegg_analysis(C06_genes)

# THINK ABOUT USING TOPGO PACKAGE
# https://bioconductor.org/packages/3.7/bioc/vignettes/topGO/inst/doc/topGO.pdf


# Tables for paper. 
tempgoa <- head(goa,row.names=FALSE)
tempgoa <- tempgoa[,1:5]
print(xtable(tempgoa, display=c("s","s","s","s","s","g")), math.style.exponents = TRUE,include.rownames = FALSE)

tempkega <- tail(kega,row.names=FALSE)
tempkega <- tempkega[,1:5]
print(xtable(tempkega, display=c("s","s","s","s","s","g")), math.style.exponents = TRUE,include.rownames = FALSE)

# rm(kega,goa,tempgoa,tempkega)


########################################################################################
# Calculate scores for all disease modules and rank them, sort decreasing numerical order
# print_dm_table() will generate the latex stuff based on annoations and ranking 
# methods to create the disease module table for the paper, containing:
#   C06 and non-C06 disease modules
#   GO enrichment counts
#   KEGG enrichment counts
#   Biological plausibility ranking
#   Current drugs / any drug reposition candidates
#   components/complexity



keggRanks <- rank_alldm_pathways(allmods)  # provides count of number of active pathways in each diseasemodule
goRanks <- rank_alldm_go(allmods)    # provides ranking of GO annotations
score <- diag(goRanks)+keggRanks  # get the combined score by simply adding KEGG rank to goRanks
score <- as.vector(score)

score <- diag(jaccard(goRanks))  # ??????

###################################################################
dist_mat <- score_alldm_go(allmods)  # cluster based scoring
optimum_clusters(dist_mat)

dm <- merge_dm(modscores,15)
dmgroup <- as.vector(dm)
dmlabel <- names(dm)
dm_df <- data.frame(dm=dmgroup,disease=dmlabel,stringsAsFactors = FALSE)

# get new overall scores on combined, new disease modules groups 
# 1. create new dismods based on cluster numbers and labels.
list_allmods <- unique(allmods$DiseaseModule)
new_dm <- c("BP", "GO:0070863","positive regulation","PLAU1","Parks", 1)
for (i in 1:length(unique(dmgroup))){
  tmp_names <- filter(dm_df,dmgroup == dmgroup[i])
  tmp_dm <- filter(allmods,DiseaseModule == tmp_names$disease)  
  #cat("\ndiseases in group ",i," are ",tmp_names$disease)
  newgroup <- rep(i, nrow(tmp_dm))
  tmp_dm <- cbind(tmp_dm,newgroup)
  new_dm <- rbind(new_dm,tmp_dm)
}
new_dm <- new_dm[-1, ]     # 1st entry is rubbish so remove it

# add score variable for each new disease module
tmp_score <-rep(0,nrow(new_dm))
new_dm <- cbind(new_dm,tmp_score)
colnames(new_dm)[7] <- "score"
# populate score variable in new_dm with the values - lots of duplications!
for(i in 1:length(unique(new_dm$newgroup))){
  new_dm[new_dm$newgroup == i, "score"] <- score[i]
}

temp_for_latex <- data.frame(score=score, newgroup=seq(1:25)) # to be sorted for latex table
temp_for_latex  <- temp_for_latex[order(temp_for_latex$score,decreasing = TRUE),]

# sum GO counts for (CC,MF,BP) inclusion into table
gocount <- rep(0,length(unique(new_dm$newgroup)))
for(i in 1:length(unique(new_dm$newgroup))){
  tempstuff <- filter(new_dm,newgroup == i)
  cat("\nGo count for group ",i," is ",nrow(tempstuff))
  gocount[i] <- nrow(tempstuff)
}

# sum gene counts for inclusion into table
genecount <- rep(0,length(unique(new_dm$newgroup)))
for(i in 1:length(unique(new_dm$newgroup))){
  tempstuff <- filter(new_dm,newgroup == i)
  genecount[i] <- length(unique(tempstuff$genes))
  cat("\ngenecount for group ",i," is ",genecount[i])
}


# 29/12/17
# get count of clusters each disease falls into, after KMeans.
tempstuff <- unique(dm_df$disease)
for(i in 1:length(tempstuff)){
  tempcount <- filter(dm_df,disease==tempstuff[i])
  cat("\n",tempstuff[i]," is in ",length(unique(tempcount$dm))," clusters.")
}



# make Latex tables
dm.table <- xtable(dm_df)
print(dm.table,floating=FALSE)
print_dm_table(sg)

# The Non-C06 diseases that are closely linked to them.
# Alzheimer Disease C10.228.140.380.100
# Asthma C08.127.108
# Autistic Disorder F03.625.164.113.500
# Carcinoma, Non-Small-Cell Lung C04.588.894.797.520.109.220.249
# Obesity C18.654.726.500
# Parkinson disease C10.228.140.079.862.500
# Arthritis, Rheumatoid C05.550.114.154
# Schizophrenia F03.700.750
# Diabetes Mellitus, Type 2 C18.452.394.750.149
# Hypertensive disease C14.907.489


# Sort C06 diseases with genes, need to attach the MeSH code - DO THIS ONLY ONCE.
#C06 <- fix_C06()





