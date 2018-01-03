

disma <- hamming.distance(mat) 
hr <- hclust(as.dist(disma)) 
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=T)


listInput <- list(alz=c(nonC06_alz$geneName),asth=c(nonC06_asth$geneName),
                  diab=c(nonC06_dia$geneName),hype=c(nonC06_hyp$geneName),
                  park=c(nonC06_park$geneName),ra=c(nonC06_ra$geneName),
                  schizo=c(nonC06_sch$geneName),obes=c(nonC06_obs$geneName))

upset(fromList(listInput), order.by = "freq",nsets = 8, number.angles = 30)

# This stckoverflow post achieves the gene overlap problem.
# https://codereview.stackexchange.com/questions/17905/compute-intersections-of-all-combinations-of-vectors-in-a-list-of-vectors-in-r
overlap <- function(l) {
  results <- lapply(l, unique)
  # combinations of m elements of list l
  for (m in seq(along=l)[-1]) {
    # generate and iterate through combinations of length m
    for (indices in combn(seq(length(l)), m, simplify=FALSE)) {
      
      # make name by concatenating the names of the elements
      # of l that we're intersecting
      name_1 <- paste(names(l)[indices[-m]], collapse="_")
      name_2 <- names(l)[indices[m]]
      name <- paste(name_1, name_2, sep="_")
      
      results[[name]] <- intersect(results[[name_1]], results[[name_2]])
      
    }
  }
  return(results)
}


