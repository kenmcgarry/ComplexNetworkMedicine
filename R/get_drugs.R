# get_drugs.R
# Ken McGarry 23/10/17
#  

library(dplyr)

setwd("C:/R-files/disease")    # point to where my code lives
load("C06disease.RData") # load in required data - the contents will change regulary
source("gi_functions.R")  # load in the functions required for finding lists of drugs 

# Create a structure to hold all known drugs treating your particular disease of interest
# Obviously instantiate these functions before calling them.

drug_list <- get_drugs("C0002395",restrictedlist)  # umls code for Alzheimers
drug_list <- get_drugs("C0020179",restrictedlist)  # umls code for Huntingdons
drug_list <- get_drugs("C0149925",restrictedlist)  # umls code for small cell lung cancer
drug_list <- get_drugs("C0007131",restrictedlist)  # umls code for non small cell lung cancer
drug_list <- get_drugs("C0024141",restrictedlist)  # umls code for systemic lupus erythematosus
drug_list <- get_drugs("C0029408",restrictedlist)  # umls code for Osteoarthritis
drug_list <- get_drugs("C0003873",restrictedlist)  # umls code for Rheumatoid arthritis

drug_list <- get_drugs_plus("C0013295",restrictedlist)  # Duodenal Ulcer.

# work with following dataframes: indications; mappings; digestive; disgene; 
# save(indications, restrictedlist, mappings, digestive, disgene, file = "C06disease.RData")






