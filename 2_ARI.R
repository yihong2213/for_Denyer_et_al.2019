# Originally written by Dr. Olivier MARTIN
# Simplified for the dataset from Denyer et al (https://doi.org/10.1016/j.devcel.2019.02.022) by me (Yi-Hong TU)


#### Input ####----

# All the variables are the ones after the general workflows

# VARIABLE(S): Denyer_1st_input, Denyer_2nd_input
## e.g. Denyer_1st_input <- Denyer_CS
## e.g. Denyer_2nd_input <- Denyer_NP
### Denyer_whole_Louvain vs Denyer_whole_Leiden
### Denyer_CS_PC100 vs Denyer_CS_PC150
### Denyer_NP_PC100 vs Denyer_NP_PC150
### Denyer_NP_PC150 vs Denyer_NP_PC200
### Denyer_CS vs Denyer_NP



#### Calculation ####----

clustering1 <- Denyer_1st_input@meta.data$seurat_clusters
clustering2 <- Denyer_2nd_input@meta.data$seurat_clusters

Ncells <- length(Denyer_1st_input@assays$RNA@counts@Dimnames[[2]])
# must be the same value as the "length" of the 2nd input

my_data_frame <- data.frame(clustering1, clustering2)
colnames(my_data_frame) <- c("cluster1_index","cluster2_index")

my_table_ns <- table(my_data_frame)

my_table_nnm1s <- my_table_ns*(my_table_ns-1)   # n*(n-1)
tmy_table_ns <- t(my_table_ns)   # work on transpose

rowtotals <- rowSums(my_table_ns)
coltotals <- rowSums(tmy_table_ns)

sum_in_row <- sum(rowtotals*(rowtotals-1))
sum_in_col <- sum(coltotals*(coltotals-1))

numerator <- sum(my_table_nnm1s) - sum_in_row*sum_in_col/(Ncells*(Ncells-1))
denominator <- 0.5*(sum_in_row + sum_in_col) - sum_in_row*sum_in_col/(Ncells*(Ncells-1))

ARI <- numerator/denominator
ARI
### Denyer_whole_Louvain vs Denyer_whole_Leiden : 0.9169847
### Denyer_CS_PC100 vs Denyer_CS_PC150 : 0.8703654
### Denyer_NP_PC100 vs Denyer_NP_PC150 : 0.6964339
### Denyer_NP_PC150 vs Denyer_NP_PC200 : 0.9178307
### Denyer_CS vs Denyer_NP : 0.02662705