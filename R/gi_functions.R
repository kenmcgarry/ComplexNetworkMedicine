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