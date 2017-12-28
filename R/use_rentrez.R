# use_rentrez.R
# use internet connection to NCBI to collect interactions rather than clumsy STITCH/STRING process, it
# assumes library(rentrez) and library(stringr) are loaded.
# https://github.com/ropensci/rentrez/wiki/Find-genes-known-to-interact-with-a-given-gene

use_rentrez <- function(mygenes){
  interactionList <- data.frame(a="RU12",b="RFU2")
  for (i in 1:length(mygenes)){
    onegene <- mygenes[i]                  # convert from name to entrez ID
    gene_search <- entrez_search(db="gene", term=str_c("(",onegene,"[GENE]) AND (Homo sapiens[ORGN])"))
    if(gene_search$count > 0){
      if(!is.null(gene_search$ids)){
        templist <- interactions_from_gene(gene_search$ids[1])  # get interaction partners for this single gene
        
        
          #cat("\nFound numbers only gene.....in module",i)
          #cat("\ntemplist=",templist)
          #readline(prompt="Press [enter] to continue")
      
          cat("\nonegene is ",onegene)
          n <- length(templist)
          genevec <- rep(onegene,n)
          tempvec <- cbind(templist,genevec)
          colnames(tempvec)<- c("a","b") 
          interactionList <- rbind(interactionList,tempvec)
          tempvec <- NULL}
      }else{
        cat("\nNo interaction partners for ",onegene)}
  }
  interactionList <- interactionList[-1,] # remove silly entry initializing interactionList
  return(interactionList)
}


interactions_from_gene <- function(gene_id){
  
  xmlrec <- entrez_fetch(db="gene", id=gene_id, rettype="xml", parsed=TRUE) 
  if(is.null(xmlrec)){return}else{
  XML::xpathSApply(xmlrec,"//Gene-commentary[Gene-commentary_heading[./text()='Interactions']]//Other-source[Other-source_src/Dbtag/Dbtag_db[./text()='GeneID']]//Other-source_anchor",
                   XML::xmlValue) }                                                                                                                                 
}
