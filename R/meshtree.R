# meshtree.R
# formally parsetree.R
# 9/12/2016, restarted 6/7/2017. parses the Mesh (Medical Subject Headings) datafile, 
# Collects subheadings based on intial main heading string
# First we need to parse it into a sensible structure and create a tree: 
# This software needs a rewrite in functional form as its brute force repetitive code.

library(data.tree)
library(rentrez)
library(dplyr)
library(tidyr)

meshtree <- file.path('C://R-files//disease//','meshtreefull.csv') %>% read.delim(na.strings='',sep=',',header=TRUE,comment.char="#")
meshtree <- meshtree[,1:3]

level1<-"C06"
digestive <- subset(meshtree,grepl(level1,meshtree$MeSH)) # search based on level 1 e.g. "C06" is digestive system

rownames(digestive)=NULL
digestive <- data.frame(lapply(digestive, as.character), stringsAsFactors=FALSE)
digestive <- digestive[!duplicated(digestive[,2]),] # keep only unique entries of diseases, remove duplicates

# Obtain Level 2 names
level2<-" "
n=7 # because the 2nd level is seven characters e.g. "C06.130"
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n)
    level2[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n),seq(n,nchar(digestive[i,1])+n-1,n)) 
}

level2 <- unique(level2)
x <- nchar(level2)
level2 <- as.character(na.omit(level2[ x %in% n ]))
# ---------------------------------------------------------------------

# Obtain Level 3 names
n3=11 # because the 3rd level is 11 characters including decimal points e.g. "C06.130.120"
level3=" "
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n3)
    level3[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n3),seq(n3,nchar(digestive[i,1])+n3-1,n3)) 
}

level3 <- unique(level3)
x <- nchar(level3)
level3 <- as.character(na.omit(level3[ x %in% n3 ]))
# ---------------------------------------------------------------------

# Obtain Level 4 names
n4=15 # because the 4th level is 15 characters including decimal points e.g. "C06.130.120"
level4=" "
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n4)
    level4[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n4),seq(n4,nchar(digestive[i,1])+n4-1,n4)) 
}

level4 <- unique(level4)
x <- nchar(level4)
level4 <- as.character(na.omit(level4[ x %in% n4 ])) # na.omit(your.data.frame)
# ---------------------------------------------------------------------

# Obtain Level 5 names
n5=19 # because the 5th level is 19 characters including decimal points e.g. "C06.130.120"
level5=""
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n5)
    level5[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n5),seq(n5,nchar(digestive[i,1])+n5-1,n5)) 
}

level5 <- unique(level5)
x <- nchar(level5)
level5 <- as.character(na.omit(level5[ x %in% n5 ])) # na.omit(your.data.frame)

# ---------------------------------------------------------------------

# Obtain Level 6 names
n6=23 # because the 6th level is 23 characters including decimal points e.g. "C06.130.120"
level6=" "
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n6)
    level6[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n6),seq(n6,nchar(digestive[i,1])+n6-1,n6)) 
}

level6 <- unique(level6)
x <- nchar(level6)
level6 <- as.character(na.omit(level6[ x %in% n6 ])) # na.omit(your.data.frame)

# ---------------------------------------------------------------------

# Obtain Level 7 names
n7=27 # because the 7th level is 27 characters including decimal points e.g. "C06.130.120"
level7=" "
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n7)
    level7[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n7),seq(n7,nchar(digestive[i,1])+n7-1,n7)) 
}

level7 <- unique(level7)
x <- nchar(level7)
level7 <- as.character(na.omit(level7[ x %in% n7 ])) # na.omit(your.data.frame)

# ---------------------------------------------------------------------

# Obtain Level 8 names
n8=31 # because the 8th level is 31 characters including decimal points e.g. "C06.130.120"
level8=" "
for ( i in 1:nrow(digestive)){
  if(nchar(digestive[i,1])== n8)
    level8[i] <- substring(digestive[i,1],seq(1,nchar(digestive[i,1]),n8),seq(n8,nchar(digestive[i,1])+n8-1,n8)) 
}

level8 <- unique(level8)
x <- nchar(level8)
level8 <- as.character(na.omit(level8[ x %in% n8 ])) # na.omit(your.data.frame)

# ---------------------------------------------------------------------

# BUILD THE TREE
# create root node (level1) and add level2 nodes

dsd <- Node$new(digestive[1,1])
for(i in 1:length(level2)){
  dsd$AddChild(level2[i])
}

# now add level 3 nodes
for (i in 1:length(level2)){
  x <- FindNode(dsd,level2[i]) # get a level2 node
  nodestoadd <- subset(level3,grepl(level2[i],level3))
  for(j in 1:length(nodestoadd))
    x$AddChild(nodestoadd[j])
}

# now add level 4 nodes
for (i in 1:length(level3)){
  x <- FindNode(dsd,level3[i]) # get a level3 node
  nodestoadd <- subset(level4,grepl(level3[i],level4))
  if(length(nodestoadd>0))
    for(j in 1:length(nodestoadd))
      x$AddChild(nodestoadd[j])
}

# now add level 5 nodes
for (i in 1:length(level4)){
  x <- FindNode(dsd,level4[i]) # get a level4 node
  nodestoadd <- subset(level5,grepl(level4[i],level5))
  if(length(nodestoadd>0))
    for(j in 1:length(nodestoadd))
      x$AddChild(nodestoadd[j])
}

# now add level 6 nodes
for (i in 1:length(level5)){
  x <- FindNode(dsd,level5[i]) # get a level5 node
  nodestoadd <- subset(level6,grepl(level5[i],level6))
  if(length(nodestoadd>0))
    for(j in 1:length(nodestoadd))
      x$AddChild(nodestoadd[j])
}


# now add level 7 nodes
for (i in 1:length(level6)){
  x <- FindNode(dsd,level6[i]) # get a level5 node
  nodestoadd <- subset(level7,grepl(level6[i],level7))
  if(length(nodestoadd>0))
    for(j in 1:length(nodestoadd))
      x$AddChild(nodestoadd[j])
}

# now add level 8 nodes
for (i in 1:length(level7)){
  x <- FindNode(dsd,level7[i]) # get a level5 node
  nodestoadd <- subset(level8,grepl(level7[i],level8))
  if(length(nodestoadd>0))
    for(j in 1:length(nodestoadd))
      x$AddChild(nodestoadd[j])
}

dsd$diseasename <- digestive[1,3]
dsd$id <- digestive[1,2]

# ADD NAMES TO TREE
# LEVEL 2: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level2)){
  thename <- digestive[which(digestive[,1] == level2[i]),3]
  x <- FindNode(dsd,level2[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level2[i]),2] 
}

# LEVEL 3: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level3)){
  thename <- digestive[which(digestive[,1] == level3[i]),3]
  x <- FindNode(dsd,level3[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level3[i]),2]
}

# LEVEL 4: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level4)){
  thename <- digestive[which(digestive[,1] == level4[i]),3]
  x <- FindNode(dsd,level4[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level4[i]),2]
}

# LEVEL 5: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level5)){
  thename <- digestive[which(digestive[,1] == level5[i]),3]
  x <- FindNode(dsd,level5[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level5[i]),2]
}

# LEVEL 6: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level6)){
  thename <- digestive[which(digestive[,1] == level6[i]),3]
  x <- FindNode(dsd,level6[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level6[i]),2]
}

# LEVEL 7: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level7)){
  thename <- digestive[which(digestive[,1] == level7[i]),3]
  x <- FindNode(dsd,level7[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level7[i]),2]
}

# LEVEL 8: now add the actual disease names to the MeSH id's as a variable name
thename=NULL
for (i in 1:length(level8)){
  thename <- digestive[which(digestive[,1] == level8[i]),3]
  x <- FindNode(dsd,level8[i])
  x$diseasename <- thename
  x$id <- digestive[which(digestive[,1] == level8[i]),2]
}

print(dsd,"diseasename","id",limit=NULL)

# tidy up please!
rm(level1,level2,level3,level4,level5,level6,level7,level8,n3,n4,n5,n6,n7,n8,thename,nodestoadd)

# -------------------------------------------------------------------------------------
# BUILD LISTS OF DISEASE GENES FOR EACH ID IF THEY EXIST 
# e.g. "Liver Cirrhosis, Biliary" -> "D008105" -> "umls:C0238065" -> "MERTK" and four other genes

# However, two excel files need to be used:
# map-umls2mesh.csv, use the ID to lookup UMLS code.
# curated_disease_gene_associtions.csv is indexed by UMLS code for any disease genes

# from DisGeneNet site # http://www.disgenet.org/web/DisGeNET/menu/home
mappings <- file.path('C://R-files//disease//','map_umls_2_mesh.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
disgene <- file.path('C://R-files//disease//','curated_gene_disease_associations.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
mappings <- data.frame(lapply(mappings, as.character), stringsAsFactors=FALSE) # needs to be strings NOT factors
disgene <- data.frame(lapply(disgene, as.character), stringsAsFactors=FALSE) # needs to be strings NOT factors
disgene <- disgene[,1:6]

# Convert id's ("D008105") into umls codes (umls:C0238065) using mappings structure
umls=""
for(i in 1:nrow(digestive)){
  if(length(mappings[which(mappings[,3] == digestive[i,2]),1])>0) # check for zero length returns when umls maps dont exist
     umls[i] <- mappings[which(mappings[,3] == digestive[i,2]),1]
  else
    umls[i]="UNKNOWN" # 68 out of 299 are "unknown!!"
  
}

umlscode = data.frame(umls, digestive,stringsAsFactors=FALSE) # add the umls to mesh ids

dglist=""
# Find any disease genes using disgene structure and the umls id, e.g. "umls:C0000768" 
for(i in 1:nrow(umlscode)){
  if(length(disgene[which(disgene[,1] == umlscode[i,1]),4])>0){ # check for zero length returns when umls maps dont exist
    temp <- disgene[which(disgene[,1] == umlscode[i,1]),4]
    dglist[i] <- paste(temp, collapse = ";")
    x <- FindNode(dsd,umlscode[i,2])
    x$genes <- temp
    x$num <- length(temp)
    cat(length(temp), x$diseasename,"\n")
  }
  else
    dglist[i]="NO ASSOCIATED GENES" # 
  
}

umlscode = data.frame(umlscode, dglist,stringsAsFactors=FALSE) # add the umls to mesh ids

print(dsd,"diseasename","id","num",limit=NULL)

# ============================================================================================
# Now create protein to disease, complex networks
library(igraph)

# grab each protein  and make a list
thegenes<-array("")
for(i in 1:nrow(umlscode)){
  if(grepl("NO ASSOCIATED GENES",umlscode[i,5]) != TRUE){
    #cat(unlist(strsplit(umlscode[i,5], ";")),"\n")
    temp <- array(unlist(strsplit(umlscode[i,5], ";")))
    thegenes <- append(thegenes,temp)
    }
}

thegenes <- unique(thegenes) # remove duplicates
thegenes <- thegenes[thegenes !=""] # remove empty string, not sure how it got there in first place!

thediseases <- data.frame(Disease=character(),Protein=character(),stringsAsFactors=FALSE) 

k=1; #index for thediseases dataframe
for(i in 1:nrow(umlscode)){ # for every disease make entry for each gene its associated with.
  if(grepl("NO ASSOCIATED GENES",umlscode[i,5]) != TRUE){
    temp <- array(unlist(strsplit(umlscode[i,5], ";")))
      for(j in 1:length(temp)){
        thediseases[k,1] <- umlscode[i,4]
        thediseases[k,2] <- temp[j]
        k=k+1
      }
  }
} # we now have a list of the diseases and associated genes (where they exist)


# what we really need is a pairwise list of diseases linked by common proteins
df <- data.frame(DiseaseA=character(),DiseaseB=character(),Protein=character(),stringsAsFactors=FALSE) 

# Identify those genes that appear twice or more (we should have 316 out of 1,905). 
# We have to use some table functions from dplyr package to get counts.
gene_tbl <- tbl_df(thediseases)
mygenes <- tally(group_by(gene_tbl, Protein))
mydiseases <- tally(group_by(gene_tbl,Disease))
mydiseases <- mydiseases[mydiseases$n > 1.0,] # diseases must appear more than once to be in a network of diseases!
mygenes <- mygenes[mygenes$n > 1.0,] # genes must appear more than once to be in a network of diseases!

netdiseases <- data.frame(Disease=character(),Protein=character(),stringsAsFactors=FALSE) 

for(i in 1:length(mydiseases$Disease)){ # for every gene make an entry for each disease its associated with.
  netgenes <- thediseases[which(thediseases[,1] == mydiseases$Disease[i]),] #look for protein(s)
  netdiseases <- rbind(netdiseases,netgenes)
  #netdiseases[k,2] <- temp[j]
} # we now have a SMALLER list of the connected diseases and associated genes

# Save workspace at this point in R-files/disease/ and its called....
# GI-disease-Oct21st2017.RData

# Build a SMALL graph of the multiple (2+diseases) connected diseases
# use igraph for conventional gene and diseases linkages
gs <- graph_from_data_frame(netdiseases, directed=FALSE)

# use TKPLOT
gs <- graph.edgelist(as.matrix(netdiseases),directed=FALSE)
ad <- get.adjacency(gs)
nodesize=5 # change the number if you want bigger nodes

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
x <- 1:ncol(ad)

nodelabel<-V(gs)$name
d1<-as.character(netdiseases[,1])
d2<-as.character(netdiseases[,2])

for ( i in 1:ncol(ad)){
  z<-(sapply(nodelabel[i],grep,d2))
  if(is.integer(z)){
    nodecolor[i]<-"pink"}
  y<-(sapply(nodelabel[i],grep,d1)) 
  if(is.integer(y)){   
    nodecolor[i]<-"lightblue"}                              
}

#--- TKPLOT allows you drag nodes around and create a better graph --------
tkplot(gs,layout = layout.kamada.kawai,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,vertex.label.cex=0.85,
       edge.arrow.size=0, edge.curved=FALSE)



# Build a graph of the LARGER network with single unconnected diseases as well as connected
# use igraph for conventional gene and diseases linkages
g <- graph_from_data_frame(thediseases, directed=FALSE)

# use TKPLOT
g <- graph.edgelist(as.matrix(thediseases),directed=FALSE)
ad <- get.adjacency(g)
nodesize=5 # change the number if you want bigger nodes

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
x <- 1:ncol(ad)

nodelabel<-V(g)$name
d1<-as.character(thediseases[,1])
d2<-as.character(thediseases[,2])

for ( i in 1:ncol(ad)){
  z<-(sapply(nodelabel[i],grep,d2))
  if(is.integer(z)){
    nodecolor[i]<-"pink"}
  y<-(sapply(nodelabel[i],grep,d1)) 
  if(is.integer(y)){   
    nodecolor[i]<-"tomato"}                              
} 

#--- TKPLOT allows you drag nodes around and create a better graph --------
tkplot(g,layout = layout.kamada.kawai,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,vertex.label.cex=0.85,
       edge.arrow.size=0, edge.curved=FALSE)

#============== link disease proteins with rest of PPI network =================
proteins<-(unique(thediseases[1001:nrow(thediseases)]))
# Write CSV in R
write.csv(proteins, file = "C:\\R-files\\disease\\proteins\\upload.csv", row.names=FALSE,col.names = FALSE)

#--------------------------------------------
# seek shared proteins between selected diseases
disA <- thediseases[which(thediseases[,1] == "Celiac Disease"),2] #look for the protein(s)
disB <- thediseases[which(thediseases[,1] == "Gastroparesis"),2] #look for the protein(s)
disC <- thediseases[which(thediseases[,1] == "Mucositis"),2] #look for the protein(s)

universe <- unique(c(disA,disB,disC))

GroupA <-universe %in% disA
GroupB <-universe %in% disB
GroupC <-universe %in% disC

## proteins that are in GroupA and in GroupB and C
sharedp <- universe[GroupA & GroupB & GroupC]




