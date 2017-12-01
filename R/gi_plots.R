# gi_plots.R
# some plots()

# This plot is for the pdens Versus clustering threshold
# assumes LinkComm has been used to create S2 structure
plot_pdens <- function(ddata){
     plot(ddata,
     xlab="Clustering threshold",
     ylab="Partition density",
     lty=1,
     col="blue",
     type="b",
     lwd=1,
     cex=1.5,
     cex.lab=1.7,
     cex.axis=1.7,
     pch=16,
     panel.first=grid())
}

# This plot is for the Community Connectedness Versus Community Modularity
# assumes LinkComm has been used to create S2 structure
plot_com <- function(cm,cc,ntitle){
  cm <- sort(cm)
  cc <- sort(cc)
  plot(cm,cc,
     xlab="Community modularity",
     ylab="Community connectedness",
     main=ntitle,
     cex.main=2,
     lty=1,
     col="red",
     type="b",
     lwd=1,
     cex=1.5,
     cex.lab=1.7,
     cex.axis=1.7,
     pch=16,
     panel.first=grid())
}

# plot community centrality for each member
plot_centrality <- function(ddata,xtext){
  plot(ddata,
       xlab=xtext,
       ylab="centrality",
       lty=1,
       col="green",
       type="b",
       lwd=1.7,
       cex=1.5,
       cex.lab=1.7,
       cex.axis=1.7,
       pch=16,
       panel.first=grid())
}





