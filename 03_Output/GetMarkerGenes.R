# Sarah Goggin
# 2020 08 04
# Identify marker genes for clusters as starting point for cluster annotation

print("Start GetMarkerGenes.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(data.table)
library(Seurat)

#####################################################################FUNCTION
#' @param exprs_mat_filepath filepath to expression matrix
#' @param clusters_in_filepath filepath to cluster ids
#' @param test Test to make comparison and determine significance of marker genes, relevant options for us are "roc", "t"
findMarkerGenesCytof <- function(exprs_mat, clusters_in, test) {
  # Find marker genes via DE 
  # First read in cluster ids
  cluster.ids <- data.frame("cluster_id" = as.character(clusters_in))
  # Expression matrix currently has genes as columns and cells as rows
  # Transpose and create matrix of data (may already be matrix)
  exprs_mat <- as.matrix(t(exprs_mat))
  # Add "names" for cells, which are now cols...
  colnames(exprs_mat) <- c(1:dim(exprs_mat)[2])
  # Initialize the Seurat object
  exprs_seurat <- CreateSeuratObject(counts = exprs_mat, 
                                     meta.data = cluster.ids)
  # Set default identity to be cluster.id
  Idents(object = exprs_seurat) <- exprs_seurat@meta.data$cluster_id
  # Find markers for all clusters, set do.print=TRUE to output progress
  markers.all <- FindAllMarkers(exprs_seurat,test.use = test, do.print = TRUE)
  
  return(markers.all)
}
#####################################################################FUNCTION


## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
MARKER.SIG.TEST <- "t" #"roc"

MARKERS.OUT.FILENAME <- "Clusters_Extracted.csv"

print("Finished reading input parameters, reading files")

## Read needed files ==============================================================================
concat.transformed <- read.concat.transformed(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)

print("Finished reading needed files, identifying marker genes")

## Subset cells based on cluster id ===============================================================
markers.all <- findMarkerGenesCytof(exprs_mat = concat.transformed,
                                    clusters_in = clusters.in, 
                                    test = MARKER.SIG.TEST)

print("Finished identifying marker genes, writing output file")

## Write output file ==============================================================================
markers.all$gene <- gsub("-[^-]*$", "", markers.all$gene)
fwrite(markers.all,file = MARKERS.OUT.FILENAME,row.names = FALSE)


print("Completed GetMarkerGenes.R")