#Corey Williams, University of Virginia
#14 May, 2020
#Downsample files to equally represent clusters

print("Start downsampleByCluster.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(data.table)
library(dplyr)

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

## Input parameters ===============================================================================
SUBSAMPLES <- 1000 #Number of cells per cluster. If a cluster contains less, all cells will be kept

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
DOWNSAMPLE.OUT.FILENAME <- "Concat_Transformed_.csv" #rename
SUBSAMPLE.ID.OUT.FILENAME <- "Subsample_ID_by_Cluster.csv" #Relative to concat_transformed file, rename
LAYOUT.OUT.FILENAME <- "UMAP_layout_Subsampled.csv" #rename
CLUSTERS.OUT.FILENAME <- "clusters_Subsampled.csv" #will be needed to cooperate with pipeline, rename

print("Finished reading input parameters, reading files")

## Read needed files ==============================================================================
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCATTRANSFORMED.FILENAME), stringsAsFactors = F)
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)

print("Finished reading needed files, downsampling")

## Perform cluster-based downsampling =============================================================
#Initialize lists for outputs
numcluster <- length(unique(clusters.in))
subsample.ids <- vector(mode = "list",length = numcluster)
expression.out <- subsample.ids
layout.out <- subsample.ids
clusters.out <- subsample.ids

#loop through each cluster
for (this_cluster in sort(unique(clusters.in))) {
  
  #Check if subsampling is needed or if all cells for cluster will be used
  if (sum(clusters.in == this_cluster) < SUBSAMPLES) {
    #get ids of cells in cluster file
    subsample.ids[[this_cluster]] <- which(clusters.in == this_cluster)
    #use subsample ids to subsample from expression info
    expression.out[[this_cluster]] <- concat.transformed[subsample.ids[[this_cluster]],]
    # store layout
    layout.out[[this_cluster]] <- layout.in[subsample.ids[[this_cluster]],]
    #store cluster id's
    clusters.out[[this_cluster]] <- clusters.in[subsample.ids[[this_cluster]]]
    
  } else {
    #Get vector of all ids for this cluster
    all.ids <- which(clusters.in == this_cluster)
    #randomly subsample from the selection of cells
    subsample.ids[[this_cluster]] <- sample(all.ids,size = SUBSAMPLES)
    #use subsample ids to subsample from expression info
    expression.out[[this_cluster]] <- concat.transformed[subsample.ids[[this_cluster]],]
    #store subsampled layout
    layout.out[[this_cluster]] <- layout.in[subsample.ids[[this_cluster]],]
    #store cluster id's
    clusters.out[[this_cluster]] <- clusters.in[subsample.ids[[this_cluster]]]
    
  }
}

#Move outputs from list format to vector or array to output
subsample.ids.out <- unlist(subsample.ids)
expression.out <- rbindlist(expression.out)
layout.out <- rbindlist(layout.out)
clusters.out <- unlist(clusters.out)

print("Finished downsampling, writing output file")

## Write output file ==============================================================================
fwrite(list(subsample.ids.out),file = SUBSAMPLE.ID.OUT.FILENAME,col.names = FALSE,sep = ",")
fwrite(expression.out,file = DOWNSAMPLE.OUT.FILENAME,row.names = FALSE)
fwrite(layout.out,file = LAYOUT.OUT.FILENAME)
fwrite(list(clusters.out),file = CLUSTERS.OUT.FILENAME,col.names = FALSE,sep = ",")

print(table(clusters.out))

print("Completed downsampleByCluster.R")