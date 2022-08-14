#### Input ####----

# All the variables are the ones after the general workflows

# VARIABLE(S): Denyer_CellDistribution, Denyer_ClusterLabel
## e.g. Denyer_CellDistribution <- Denyer_CS
## e.g. Denyer_ClusterLabel <- Denyer_NP
### Denyer_NP labeled with Denyer_CS
### Denyer_CS labeled with Denyer_NP
### Denyer_whole labeled with Denyer_CS
### Denyer_whole labeled with Denyer_NP



#### Applying the labels to another clustering plot ####----

Denyer_CellDistribution@active.ident <- Denyer_ClusterLabel@active.ident
Denyer_CellDistribution@meta.data$RNA_snn_res.0.5 <- Denyer_ClusterLabel@meta.data$RNA_snn_res.0.5
Denyer_CellDistribution@meta.data$seurat_clusters <- Denyer_ClusterLabel@meta.data$seurat_clusters
# The key factor is "active.ident" (rather than "RNA_snn_res.0.5" or "seurat_clusters")

DimPlot(Denyer_CellDistribution, dims = c(1, 2), reduction = "umap")


### EXCLUDING the label of cluster 0 :
Denyer_CellDistribution_NoClus0 <- subset(Denyer_CellDistribution, subset = (seurat_clusters != 0))

DimPlot(Denyer_CellDistribution_NoClus0, dims = c(1, 2), reduction = "umap")
# not so meaningful for "Denyer_NP labeled with Denyer_CS" & "Denyer_whole labeled with Denyer_CS" without cluster 0