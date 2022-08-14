#### Input ####----

# All the variables are the ones after the general workflows

# VARIABLE(S): Denyer, cluster_markers, target_cluster

## Denyer <- Denyer_CS_Res0.1

## cluster_markers <- cluster_markers_CS_Res0.1

## target_cluster: from 0 to 5, manually



#### Ordering the marker genes within clusters ####----

num_clusters <- length(levels(Denyer@meta.data$seurat_clusters))
cluster_markers_reo <- data.frame()

for (x in 0:(num_clusters-1)){
  each_cluster <- cluster_markers[cluster_markers$cluster == x, ]
  order_saved <- order(each_cluster$avg_log2FC, decreasing = TRUE)
  each_reo <- each_cluster[order_saved, ]
  cluster_markers_reo <- rbind(cluster_markers_reo, each_reo)
}

cluster_markers_reo



#### Calculation ####----

how_many_top_markers <- 10  # We choose 10

marker_chosen <- cluster_markers_reo[cluster_markers_reo$cluster == target_cluster, ]$gene[1:how_many_top_markers]
result_chosen <- data.frame(cluster_number = rep(NA, times = how_many_top_markers), marker_gene = rep(NA, times = how_many_top_markers), specificity_value = rep(NA, times = how_many_top_markers), p_val = rep(NA, times = how_many_top_markers), avg_log2FC = rep(NA, times = how_many_top_markers))


for (a in 1:how_many_top_markers) {
  target_marker <- marker_chosen[a]
  
  given_marker <- Denyer[Denyer@assays$RNA@counts@Dimnames[[1]] == target_marker]
  expression_given_MarkerAndCluster <- data.frame(cluster_number = rep(NA, times = num_clusters), averaged_gene_count = rep(NA, times = num_clusters))
  
  for (b in (0:(num_clusters-1))) {
    given_cluster <- subset(given_marker, subset = seurat_clusters == b)
    how_many_cells <- length(given_cluster@assays$RNA@counts@Dimnames[[2]])
    averaged_gene_count <- ((sum(given_cluster@assays$RNA@counts@x)) / how_many_cells)
    expression_given_MarkerAndCluster[b+1, ] <- c(b, averaged_gene_count)
  }
  
  expression_target_cluster <- expression_given_MarkerAndCluster[target_cluster+1, "averaged_gene_count"]
  expression_column <- expression_given_MarkerAndCluster[ ,"averaged_gene_count"]
  expression_highest <- max(expression_column)
  expression_2nd_highest <- max(expression_column[expression_column != expression_highest])
  p_value <- cluster_markers_reo[(cluster_markers_reo$gene == target_marker & cluster_markers_reo$cluster == target_cluster), ]$p_val
  averaged_log_fold_change <- cluster_markers_reo[(cluster_markers_reo$gene == target_marker & cluster_markers_reo$cluster == target_cluster), ]$avg_log2FC
  
  if(expression_target_cluster == expression_highest){
    speci_val <- (expression_target_cluster - expression_2nd_highest) / expression_target_cluster
    result_chosen[a, ] <- c(target_cluster, target_marker, speci_val, p_value, averaged_log_fold_change)
  }else{
    result_chosen[a, ] <- c(target_cluster, target_marker, "bizarre", p_value, averaged_log_fold_change) 
  }
}

order_saved <- order(result_chosen$specificity_value, decreasing = TRUE)
result_chosen <- result_chosen[order_saved, ]

result_chosen



#### Identifying marker genes (manually) ####----

# Then, map these marker genes with the atlas from Wendrich et al. (https://bioit3.irc.ugent.be/plant-sc-atlas/) to determine which are the most appropriate

# also check their expression profiles:
## VARIABLE(S): Name_of_marker_gene (e.g. AT3G23830 for community 0)
FeaturePlot(Denyer, features = "Name_of_marker_gene", reduction = "umap")



#### Labeling with cell types ####----

new.cluster.ids <- c("Non-vascular tissue", "Lateral root cap", "Stele", "Ground tissue", "Atrichoblast", "Trichoblast")
names(new.cluster.ids) <- levels(Denyer)  # levels: 0, 1, 2, 3, 4, 5
Denyer <- RenameIdents(Denyer, new.cluster.ids)

DimPlot(Denyer, reduction = "umap")