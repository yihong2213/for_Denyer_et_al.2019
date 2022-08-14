#### Library ####----

library(dplyr)
library(Seurat)
library(patchwork)



#### Loading & QC ####----

setwd("")  # import the 1st replicate

Denyer.data <- ReadMtx(
  mtx = "matrix.mtx",
  cells = "barcodes.tsv",
  features = "merger.tsv",  # "merger.tsv" is the feature file added with Dr. Blein's gene annotations
  cell.sep = "\t",
  feature.sep = "\t",
  skip.feature = 1,
  feature.column = 1,  # Otherwise, R doesn't read the column we want since the default is 2
)
# https://satijalab.org/seurat/reference/readmtx

# cell ID (barcodes) : 6411 (1st replicate) // 5314 (2nd replicate)
# gene name (features) : 83581 (1st replicate) // 83581 (2nd replicate)
# columns in the matrix (mtx) files: gene number (gene ID), cell number (cell ID), RNA counts (RNA transcripts)


Denyer <- CreateSeuratObject(counts = Denyer.data, min.cells=0, min.features=0, project = "replicate_1")
# The recommended parameters in the Seurat tutorial of PBMC are "min.cells=3, min.features=200" (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)


Denyer[["percent.mt"]] <- PercentageFeatureSet(Denyer, pattern = "^ATM")
# ATM: realted to mitochondria
Denyer[["percent.ch"]] <- PercentageFeatureSet(Denyer, pattern = "^ATC")
# ATC: related to chloroplasts

# The steps below are the same for chloroplasts, mitochondria, nCount:
## "hist(Denyer@meta.data$nFeature_RNA)" to see the distribution
## "summary(Denyer@meta.data$nFeature_RNA)" to know the general info
## "quantile(Denyer@meta.data$nFeature_RNA, probs = c(0.05, 0.25, 0.75, 0.95))"
## probs = 0.05 & 0.95 : to filter out 5%
## probs = 0.25 & 0.75 : to see the 1st & 3rd quantile for checking, should be the same as what we see in "summary"

# Thus, for replicate 1 we choose :
## no filters for mt & ch since most of them are 0 (in agreement with the hypothetical result using nuclear scRNA-seq (no mitochondria & chloroplasts)) [i.e. no effect for lower quantiles]
Denyer <- subset(Denyer, subset = nFeature_RNA < 7503.0 & nCount_RNA < 54341.5)
# 6056 cells, 83581 genes


Denyer1 <- Denyer
rm(Denyer, Denyer.data)



setwd("")  # import the 2nd replicate

Denyer.data <- ReadMtx(
  mtx = "matrix.mtx",
  cells = "barcodes.tsv",
  features = "merger.tsv",  # "merger.tsv" is the feature file added with Dr. Blein's gene annotations
  cell.sep = "\t",
  feature.sep = "\t",
  skip.feature = 1,
  feature.column = 1,  # Otherwise, R doesn't read the column we want since the default is 2
)


Denyer <- CreateSeuratObject(counts = Denyer.data, min.cells=0, min.features=0, project = "replicate_2")


Denyer[["percent.mt"]] <- PercentageFeatureSet(Denyer, pattern = "^ATM")
Denyer[["percent.ch"]] <- PercentageFeatureSet(Denyer, pattern = "^ATC")

# Thus, for replicate 2 we choose :
## no filters for mt & ch since most of them are 0 (in agreement with the hypothetical result using nuclear scRNA-seq (no mitochondria & chloroplasts)) [i.e. no effect for lower quantiles]
Denyer <- subset(Denyer, subset = nFeature_RNA < 7865.4 & nCount_RNA < 63057.80)
# 5013 cells, 83581 genes


Denyer2 <- Denyer
rm(Denyer, Denyer.data)



#### Merging replicate 1 & 2 ####----

Denyer_whole = merge(Denyer1, Denyer2, merge.data = TRUE, project = "whole", add.cell.ids = c("R1", "R2"))
# merge: merge Seurat objects
# 11069 cells, 83581 genes

rm(Denyer1, Denyer2)



#### Separation into the CS & NP subsets ####----

setwd("")
tag <- read.table(file = "gene_type.tsv", header = TRUE)  # "gene_type.tsv" is the file of gene annotations from Dr. Blein

coding <- subset(tag, classification == "coding")
structural <- subset(tag, classification == "structural")
CS <- rbind(coding, structural)  # rbind: bind 2 dataframes

noncoding <- subset(tag, classification == "noncoding")
pseudogene <- subset(tag, classification == "pseudogene")
NP <- rbind(noncoding, pseudogene)

# the numbers of genes by viewing the file directly:
## coding: 29077
## structural: 1051
## noncoding: 18129
## pseudogene: 813
## We don't want TE (transposable element)

Denyer_CS <- subset(Denyer_whole, features = CS$gene_id)
# 11069 cells, 30055 genes
Denyer_NP <- subset(Denyer_whole, features = NP$gene_id)
# 11069 cells, 6792 genes



#### Input ####----

# VARIABLE(S): Denyer

## Denyer <- Denyer_whole
## Denyer <- Denyer_CS
## Denyer <- Denyer_NP



#### Normalization ####----

Denyer <- NormalizeData(Denyer)



#### Identification of highly variable features (feature selection) ####----

Denyer <- FindVariableFeatures(Denyer, selection.method = "vst", nfeatures = 2000)



#### Scaling ####----

Denyer <- ScaleData(Denyer)  # default: variable features



#### PCA ####----

Denyer <- RunPCA(Denyer, features = VariableFeatures(object = Denyer), npcs = 150)
# PC-relevant parameter "npc": 100, 150, 200
## -> Denyer_CS_PC100, Denyer_CS_PC150
## -> Denyer_NP_PC100, Denyer_NP_PC150, Denyer_NP_PC200

ElbowPlot(Denyer, ndims = 150)
# visualization: elbow plot
# PC-relevant parameter "ndims"

# decision: PC = 150



#### UMAP ####----

Denyer <- RunUMAP(Denyer, reduction="pca", dims = 1:150, n.components = 4L)
# PC-relevant parameter "dims"



#### Community detection (Clustering) ####----

Denyer <- FindNeighbors(Denyer, reduction = "pca", dims = 1:150, k.param = 20)  
# PC-relevant parameter "dims"

Denyer <- FindClusters(Denyer, resolution = 0.5, algorithm = 2)
# "algorithm" = 4 for the Leiden
## -> Denyer_whole_Leiden
# "resolution" = 0.1 for part 3 of [Results & Discussion]
## -> Denyer_CS_Res0.1



#### Visualization ####----

DimPlot(Denyer, dims = c(1, 2), reduction = "umap")
# visualization: clustering plot


### ONLY for replicate 1 vs replicate 2:
DimPlot(Denyer, dims = c(1, 2), reduction = "umap", group.by = "orig.ident")
# "orig.ident": the "project" names I gave to replicate 1 & 2



# VARIABLES FOR DOWNSTREAM ANALYSES:

## Denyer -> Denyer_whole
### (i.e. Denyer_whole_PC150_Res0.5_Louvain OR Denyer_whole_Louvain)
## Denyer -> Denyer_whole_Leiden
### (i.e. Denyer_whole_PC150_Res0.5_Leiden)

## Denyer -> Denyer_CS_PC100
### (i.e. Denyer_CS_PC100_Res0.5_Louvain)
## Denyer -> Denyer_CS_PC150
### (i.e. Denyer_CS_PC150_Res0.5_Louvain OR Denyer_CS)

## Denyer -> Denyer_NP_PC100
### (i.e. Denyer_NP_PC100_Res0.5_Louvain)
## Denyer -> Denyer_NP_PC150
### (i.e. Denyer_NP_PC150_Res0.5_Louvain OR Denyer_NP)
## Denyer -> Denyer_NP_PC200
### (i.e. Denyer_NP_PC200_Res0.5_Louvain)

## Denyer -> Denyer_CS_Res0.1
### (i.e. Denyer_CS_PC150_Res0.1_Louvain)