# Sarah Goggin
# 2020 08 04
# Identify marker genes for clusters as starting point for cluster annotation

print("Start GetMarkerGenes.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(data.table)
library(Seurat)

#Functions for Zunder lab pipeline
#' Metadata file input
#'
#' Reads metadata file and produces format friendly to scripts in Zunder lab pipeline. Metadata
#' strategy inspired by Nowicka et al., 2017, F1000 Research
#' @param input.folder directory containing metadata file
#' @param md.filename metadata filename, defaults to metadata.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.metadata <- function(input.folder,md.filename = "/metadata.csv"){
  md <- read.csv(paste0(input.folder,md.filename))
  #make filenames character vectors
  md$file_name <- as.character(md$file_name)
  return(md)
}

#' Panel file input
#'
#' Reads panel file and produces format friendly to scripts in Zunder lab pipeline
#' @param input.folder directory containing metadata file
#' @param panel.filename panel filename, defaults to panel.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.panel <- function(input.folder,panel.filename = "/panel.csv"){
  panel <- read.csv(paste0(input.folder,panel.filename))
  panel$Antigen <- gsub("-", "_", panel$Antigen)
  panel$Antigen <- gsub("\\.", "_", panel$Antigen)
  return(panel)
}

#' Concat_Transformed file input
#'
#' Reads Concat_Transformed file for use in Zunder lab pipeline
#' @param input.folder directory containing concat_transformed file
#' @param concat.transformed.filename concat_transformed filename, defaults to
#' Concat_Transformed.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.concat.transformed <- function(input.folder,
                                    concat.transformed.filename = "/Concat_Tranformed.csv"){
  concat.transformed <- read.csv(paste0(input.folder,concat.transformed.filename))
  colnames(concat.transformed) <- gsub("-", "_", colnames(concat.transformed))
  colnames(concat.transformed) <- gsub(".", "_", colnames(concat.transformed), fixed = TRUE)
  return(concat.transformed)
}

#' Layout file input
#'
#' Reads uMAP or other 2D layout for plotting
#' @param input.folder directory containing layout file
#' @param layout.filename layout filename, defaults to UMAP_layout.csv as outputted by generation
#'  script
#' @importFrom utils read.csv
#' @export
read.layout <- function(input.folder,layout.filename = "/UMAP_layout.csv"){
  layout <- read.csv(paste0(input.folder,layout.filename))
  return(layout)
}

#' Cluster file input
#'
#' Reads clusters from pipeline
#' @param input.folder directory containing layout file
#' @param clusters.filename clustering replicate filename, defaults to clusters.csv as outputted by
#'  generation script
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @export
read.clusters <- function(input.folder,clusters.filename = "/clusters.csv"){
  clusters.in <- fread(paste0(input.folder,clusters.filename))
  if (min(clusters.in)==0){
    #Add 1 since clustering done in Python begins indexing at zero
    clusters <- as.vector(t(clusters.in) + 1)
  } else{
    clusters <- as.vector(t(clusters.in))
  }
  
  return(clusters)
}

#' Clustering replicate file input
#'
#' Reads clustering replicates for consensus clustering and stability analysis
#' @param input.folder directory containing layout file
#' @param clustering.rep.filename clustering replicate filename, defaults to
#' ClusterStabilityCheck.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.clustering.rep <- function(input.folder,clustering.rep.filename = "/ClusterStabilityCheck.csv")
{
  #Add 1 since clustering done in Python begins indexing at zero
  clustering.rep <- read.csv(paste0(input.folder,clustering.rep.filename),header = FALSE) + 1
  return(clustering.rep)
}

#' Cluster stability file input
#'
#' Reads cluster stability values
#' @param input.folder directory containing layout file
#' @param cluster.stability.filename cluster stability filename, defaults to
#' Final_Cluster_stability.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.cluster.stability <- function(input.folder,
                                   cluster.stability.filename = "/Final_Cluster_Stability.csv"){
  cluster.stability <- as.vector(read.csv(paste0(input.folder,cluster.stability.filename),
                                          header = FALSE)[[1]])
  return(cluster.stability)
}
#' Pull out markers for transformation
#'
#' Finds markers in panel file marked for either clustering or plotting
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.markers <- function(panel){
  transform.markers <- as.character(panel$Metal[panel$Clustering == 1 | panel$Plotting == 1])
  return(transform.markers)
}

#' Get marker names in more legibile format for Concat_Transformed file
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.annotate <- function(panel){
  transform.markers.annotate <- as.character(paste0(panel$Antigen[panel$Clustering == 1 |
                                                                    panel$Plotting == 1],
                                                    "_",panel$Metal[panel$Clustering == 1 |
                                                                      panel$Plotting == 1]))
  return(transform.markers.annotate)
}

#' Get clustering variables in legible format
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.clustering.annotate <- function(panel){
  clustering.markers.annotate <- as.character(paste0(panel$Antigen[panel$Clustering == 1],"_",
                                                     panel$Metal[panel$Clustering == 1]))
  return(clustering.markers.annotate)
}

#' Get panel variables in legible format
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.plotting.annotate <- function(panel){
  plotting.markers.annotate <- as.character(paste0(panel$Antigen[panel$Plotting == 1],"_",
                                                   panel$Metal[panel$Plotting == 1]))
  return(plotting.markers.annotate)
}

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