# diseasecomorbid.R
# C:\\R-files\\disease
# This is about disease comordbity, side-effects and Pharmacovigilance
# Pharmacovigilance: integrating association rules mined from disease comorbidity data with biological 
# information for automated knowledge generation
#
# Ken McGarry  3/7/15, 26/7/16
# Human disease network showing interconnectedness of many diseases
# data is downloaded from http://barabasilab.neu.edu/projects/hudine/resource/data/data.html#Data\
# I have used the bulk download ICD 3-digit data
# You will also need to download their paper:
# "Hidalgo, Blumm, Barabasi and Christakis, A Dynamic Network Approach for the study of Human Phenotypes, 
# Plos Computational Biology, April 2009, Vol5, issue 4, e1000353".

library(linkcomm)
library(icd) # package is used to convert ICD codes into text descriptions, e.g 2.0 -> "Typhoid and paratyphoid fevers"
library(igraph)

# loading in file takes about 5 minutes, 6 million obs of 10 variables.
disease <- read.table('C:\\R-files\\disease\\AllNet5.net',header=FALSE);

# ----Column headings----
# [1] ICD-9 code disease 1
# [2] ICD-9 code disease 2
# [3] Prevalence disease 1
# [4] Prevalence disease 2
# [5] Co-ocurrence between diseases 1 and 2
# [6] Relative Risk
# [7] Relative Risk 99% Conf. Interval (left)
# [8] Relative Risk 99% Conf. Interval (right)
# [9] Phi-correlation
# [10] t-test value

#icd9Explain(icd9DecimalToShort(2.0),isShort = TRUE)
#icd9Explain(icd9DecimalToShort(disease[5000,1]),isShort = TRUE) # patient 5000, 1st disease
#icd9Explain(icd9DecimalToShort(disease[5000,2]),isShort = TRUE) # patient 5000, 2nd disease

# Scott Wilkes: Ovarian cancer risk with ovulation induction with clomiphene citrate.
# In clinical practice because no direct causal link is established (i.e. only association) 
# we continue to use it (responsibly) a lot.

comorbidity <- disease[disease$V1 == 183,] 

# Top ten diseases associated with YOUR disease
# Here I've sorted the data using order() on the V4 variable i.e. the prevalance of disease 2, '-' sign means descending order
newcomorbidity <- comorbidity[order(-comorbidity$V4),]










