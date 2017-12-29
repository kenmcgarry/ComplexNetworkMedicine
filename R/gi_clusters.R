# gi_clusters.R
# 27/12/17  Using new ranking techniques
library(factoextra)
library(NbClust)

df <- dist_mat

fviz_nbclust(df, kmeans, method = "wss", k.max = 20) +
  geom_vline(xintercept = 10, linetype = 2) +
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette",k.max = 20) +
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
nb <- fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 500, k.max = 20) +
        labs(subtitle = "Gap statistic method")

nb <- NbClust(df, distance = "euclidean", min.nc = 2, max.nc = 20, method = "kmeans")
fviz_nbclust(nb)
