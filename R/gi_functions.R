# gi_functions.R

## --------------------- FUNCTION DEFINITIONS -----------------------

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

get_drugs_plus <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls) 
  mydrugs <- ilist # temp storage
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
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
    return(NULL)}
}



