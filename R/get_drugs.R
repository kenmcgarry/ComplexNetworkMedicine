# get_drugs.R
# Ken McGarry 23/10/17

library(dplyr)

setwd("C:/R-files/disease")    # point to where my code lives
load("C06disease-15thNov2017.RData") # load in required data - the contents will change regulary
source("gi_functions.R")  # load in the functions required for finding lists of drugs 

# Work with following dataframes: indications; mappings; digestive; disgene; 
# save(indications,restrictedlist,mappings,digestive,disgene,gene_list,drug_list,simMA,fpdrugs,drugnames,drugids, file = "C06disease-15thNov2017.RData")
# the key to linking these files is the "ID" in digestive also called "meshId" in mappings

# ----- Stages:  ----------
# 1. Use "digestive" dataframe to parse through list of digestive diseases using the ID.
# 2. We need to get the drugs associated with each disease.
# 3. Get genes associated (if known) with each disease from "disgene" dataframe.
# 4. Get chemical structures of drugs and create fingerprints.
# 5. Get 2nd shell protein interactions

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


# Code for stage 4. ------------------------
# See drugstructure_gi.R for details
# generated FingerPrints of 1024 bits for 189 drugs with chemical data
# generated dendrogram plot and heatmap plot for the 189 drugs.
# dataframes are: simMA; drugids; drugnames; fpdrugs; 
# need the following packages: library(ChemmineR); library(ChemmineOB); library(ape); library(sparcl)
# library(cluster); library(xtable); library(gplots); library(scatterplot3d) 

plot_chemsim()


# THE ABOVE CODE IS ONLY PART COMPLETED AS WE NEED TO EXAMINE 2ND SHELL DRUGS AND COMPARE

# Code for stage 5. ------------------------
# Collect proteins known to interact with our long list of drugs and with the targeted proteins of those drugs




