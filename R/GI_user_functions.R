# GI_user_functions.R
# requires(reactomePA)
# requires(clusterprofiler)

kegg_analysis <- function(yourgenes){
  #cat("\nyourgenes are: ",yourgenes)
  sym2ent <- bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db");
  sym2ent <- sym2ent[,2]
  kk <- enrichKEGG(gene= sym2ent, organism= 'hsa', pvalueCutoff = 0.05)
  
  return(kk)
}



