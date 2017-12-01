# drug_targets.R
# Obtained known drug targets from 
# http://www.guidetopharmacology.org/lists.jsp

drug_targets <- drug_targets[,c(1,4,6)]
drug_targets <- file.path('C://R-files//disease','drug.target.interaction.tsv.gz') %>% read.delim(na.strings='',header =TRUE,stringsAsFactors = FALSE) 
names(drug_targets)[names(drug_targets)=="DRUG_NAME"] <- "DrugName"
names(drug_targets)[names(drug_targets)=="TARGET_CLASS"] <- "TargetClass"
names(drug_targets)[names(drug_targets)=="GENE"] <- "Gene"
drug_targets <- drug_targets[,c(1,4,6)]  # We only need three variables
drug_targets$DrugName <- firstup(drug_targets$DrugName)   # convert first letter to uppercase

# This batch of code is not hard coded - it is used to generate the details for table 2 in the paper!
# Enter drug name and obtain details for table.
tli.table <- xtable(filter(drug_targets,DrugName == "Erythromycin"))
print(tli.table,floating=FALSE)
filter(drug_list,drugbank_name == "Erythromycin")




