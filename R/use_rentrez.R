# use_rentrez.R
# use internet connection to NCBI to collect interactios rather than clumsy STITCH/STRING process, it
# assumes library(rentrez) and library(stringr) are loaded.
# https://github.com/ropensci/rentrez/wiki/Find-genes-known-to-interact-with-a-given-gene

use_rentrez <- function(mygenes){
  interactionList <- data.frame(a="RU12",b="FU2")
  
  for (i in 1:length(mygenes)){
    onegene <- mygenes[i]
    #cat("\nlength mygenes=",length(mygenes))
    gene_search <- entrez_search(db="gene", term=str_c("(",onegene,"[GENE]) AND (Homo sapiens[ORGN])"))
  
    if(!is.null(gene_search$ids)){
      templist <- interactions_from_gene(gene_search$ids)
      cat("\nonegene is ",onegene)
      n <- length(templist)
      #cat("\nlength templist=",n)
      genevec <- rep(onegene,n)
      #cat("\nlength genevec=",length(genevec))
      tempvec <- cbind(templist,genevec)
      colnames(tempvec)<- c("a","b") 
      #cat("\nnames tempvec",names(tempvec))
      #cat("\nnames interactionList",names(interactionList))
      
      interactionList <- rbind(interactionList,tempvec)
      tempvec <- NULL
    }else{
      cat("\nNo interaction partners for ",onegene)
    }
    
  }
  interactionList <- interactionList[-1,] # remove silly entry initializing interactionList
  return(interactionList)
}


interactions_from_gene <- function(gene_id){                                                                                                                                                                         
  xmlrec <- entrez_fetch(db="gene", id=gene_id, rettype="xml", parsed=TRUE)                                                                                                                                          
  XML::xpathSApply(xmlrec,                                                                                                                                                                                           
                   "//Gene-commentary[Gene-commentary_heading[./text()='Interactions']]//Other-source[Other-source_src/Dbtag/Dbtag_db[./text()='GeneID']]//Other-source_anchor",
                   XML::xmlValue)                                                                                                                                  
}