# gi_functions.R

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
    tempx <- paste("umls:",drugs[i,3],sep="") # add the bloody "umls:" back, and magic numbers: beware
    tempy <- filter(mappings, umls == tempx) 
    tempz <- filter(digestive, ID == tempy$meshId)
    x[i] <- tempz$MeSH
    y[i] <- tempz$ID
  }
  
  drugsplus <- cbind(drugs,x)   # add the x or Mesh values
  drugsplus <- cbind(drugsplus,y)   # add the y or ID values
  
  colnames(drugsplus)[5] <- "MeSH"
  colnames(drugsplus)[6] <- "ID"
  
  return(drugsplus)
}





