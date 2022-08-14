#### Input ####----

# All the variables are the ones after the general workflows

# VARIABLE(S): Denyer, cluster_markers

## Denyer <- Denyer_CS
## Denyer <- Denyer_NP

## cluster_markers <- cluster_markers_CS
## cluster_markers <- cluster_markers_NP



#### The limitation owing to the number of non-coding markers ####----

nrow(cluster_markers)
# CS: 14386
# NP: 42

how_many_markers <- 42


#### Ordering the marker genes altogether ####----

order_saved <- order(cluster_markers$avg_log2FC, decreasing = TRUE)
cluster_markers_reo <- cluster_markers[order_saved, ]



#### Calculation ####----

result_limited <- cluster_markers_reo[1:how_many_markers, ]
num_clusters <- length(levels(Denyer@meta.data$seurat_clusters))
specificity_list <- c()
cluster_most_specific_instead <- c()


for (a in 1:nrow(result_limited)) {
  target_marker <- result_limited[a, "gene"]
  target_cluster <- result_limited[a, "cluster"]
  
  given_marker <- Denyer[Denyer@assays$RNA@counts@Dimnames[[1]] == target_marker]
  expression_given_MarkerAndCluster <- data.frame(cluster_number = rep(NA, times = num_clusters), averaged_gene_count = rep(NA, times = num_clusters))
  
  for (b in (0:(num_clusters-1))) {
    given_cluster <- subset(given_marker, subset = seurat_clusters == b)
    how_many_cells <- length(given_cluster@assays$RNA@counts@Dimnames[[2]])
    averaged_gene_count <- ((sum(given_cluster@assays$RNA@counts@x)) / how_many_cells)
    expression_given_MarkerAndCluster[b+1, ] <- c(b, averaged_gene_count)
  }
  
  expression_target_cluster <- expression_given_MarkerAndCluster[expression_given_MarkerAndCluster$cluster_number == target_cluster, ]$averaged_gene_count
  expression_column <- expression_given_MarkerAndCluster[ ,"averaged_gene_count"]
  expression_highest <- max(expression_column)
  expression_2nd_highest <- max(expression_column[expression_column != expression_highest]) 
  
  speci_val <- (expression_highest - expression_2nd_highest) / expression_highest
  specificity_list[a] <- speci_val
  
  if(expression_target_cluster == expression_highest){
    cluster_most_specific_instead[a] <- "not bizarre"
  }else{
    cluster_most_specific_instead[a] <- expression_given_MarkerAndCluster[which.max(expression_column), "cluster_number"]
  }
}

result_limited$specificity_value <- specificity_list
result_limited$cluster_most_specific_instead <- cluster_most_specific_instead

result_limited



#### Statistics ####----

mean(result_limited$specificity_value)
# CS: 0.8993138
# NP: 0.314443

var(result_limited$specificity_value)
# CS: 0.02507433
# NP: 0.07261533

hist(result_limited$specificity_value, breaks = 10)


# PREVIOUS VARIABLE(S): result_limited

## result_limited -> result_limited_CS
## result_limited -> result_limited_NP

ks.test(result_limited_CS$specificity_value, result_limited_NP$specificity_value)
# p-value: 2.2e-16 (< 0.05)