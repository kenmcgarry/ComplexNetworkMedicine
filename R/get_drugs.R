# get_drugs.R
# Ken McGarry 23/10/17

library(dplyr)

setwd("C:/R-files/disease")    # point to where my code lives
load("C06disease.RData") # load in required data - the contents will change regulary
source("gi_functions.R")  # load in the functions required for finding lists of drugs 

# Create a structure to hold all known drugs treating your particular disease of interest
# Obviously instantiate these functions before calling them.

drug_list <- get_drug_names("C0002395",restrictedlist)  # umls code for Alzheimers
drug_list <- get_drug_names("C0020179",restrictedlist)  # umls code for Huntingdons
drug_list <- get_drug_names("C0149925",restrictedlist)  # umls code for small cell lung cancer
drug_list <- get_drug_names("C0007131",restrictedlist)  # umls code for non small cell lung cancer
drug_list <- get_drug_names("C0024141",restrictedlist)  # umls code for systemic lupus erythematosus
drug_list <- get_drug_names("C0029408",restrictedlist)  # umls code for Osteoarthritis
drug_list <- get_drug_names("C0003873",restrictedlist)  # umls code for Rheumatoid arthritis

drug_list <- get_drugs_plus("C0013295",restrictedlist)  # Duodenal Ulcer.

# work with following dataframes: indications; mappings; digestive; disgene; 
# save(indications, restrictedlist, mappings, digestive, disgene, file = "C06disease.RData")
# the key to linking these files is the "ID" in digestive also called "meshId" in mappings

# 1. use digestive dataframe to parse through list of digestive diseases using the ID,
# 2. we need to get the drugs associated with each disease.

disease_umls <- id2umls(digestive)  # get umls codes for each disease

for (i in 1:nrow(disease_umls)){
  drug_list <- get_drugs_plus(disease_umls[i,1],restrictedlist)  
  cat("\n i is now ...",i)
}


drug_list <- get_drugs_plus("C0333355",restrictedlist)







