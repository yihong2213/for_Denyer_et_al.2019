#### Input ####----

# All the variables are the ones after the general workflows

# VARIABLE(S): Denyer
## Denyer <- Denyer_CS
## Denyer <- Denyer_NP
## Denyer <- Denyer_CS_Res0.1



#### Calculation ####----

cluster_markers <- FindAllMarkers(Denyer, only.pos = TRUE)
# "only.pos = TRUE" : because we only want the gene expression in 1 given cluster higher than the others
# "logfc.threshold" : default is 0.25
# "min.pct" : default is 0.1


### ONLY for the NP subset:
cluster_markers <- FindAllMarkers(Denyer, only.pos = TRUE, logfc.threshold = 0.1)
# lower "logfc.threshold", otherwise there's no marker for cluster 2 in the NP subset



# VARIABLES FOR DOWNSTREAM ANALYSES:

## cluster_markers -> cluster_markers_CS_Res0.5
## cluster_markers -> cluster_markers_NP_Res0.5
## cluster_markers -> cluster_markers_CS_Res0.1
