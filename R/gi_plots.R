# gi_plots.R
# some plots()

# This plot is for the pdens Versus clustering threshold
# assumes LinkComm has been used to create S2 structure
plot(d1$pdens,
     xlab="Clustering threshold",
     ylab="Partition density",
     lty=1,
     col="blue",
     type="b",
     lwd=1,
     cex=1,
     pch=16,
     panel.first=grid())

# This plot is for the Community Connectedness Versus Community Modularity
# assumes LinkComm has been used to create S2 structure
plot(cmd1,ccd1,
     xlab="Community modularity",
     ylab="Community connectedness",
     lty=1,
     col="red",
     type="b",
     lwd=1,
     cex=1,
     pch=16,
     panel.first=grid())






