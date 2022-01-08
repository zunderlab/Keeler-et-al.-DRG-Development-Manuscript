#Austin Keeler, adapated from Corey Williams, University of Virginia
#01 August, 2020
#Plot cells over threshold MARKER EXPRESSION


print("Start PlotBySilhouetteScore_Threshold.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ggfortify)
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

print("libraries loaded")

## Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CLUSTERS.FILENAME <- "/clusters.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CONCAT.TRANSFORMED <- "/Concat_Transformed.csv"
CONCAT.OUT.FILENAME <- "Concat_Transformed_.csv"
LAYOUT.OUT.FILENAME <- "UMAP_layout_.csv"
CLUSTERS.OUT.FILENAME <- "clusters_.csv" #will be needed to cooperate with pipeline
EXPRESSION.VALUE <-  #value to extract, change line 57 to the correct marker
  
  print("input parameters loaded, reading needed files")

## Read needed files
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)
concat.transformed.in <- fread(paste0(INPUT.FOLDER,CONCAT.TRANSFORMED), stringsAsFactors = F)
colnames(concat.transformed.in) <- gsub("-", "_", colnames(concat.transformed.in))
colnames(concat.transformed.in) <- gsub(".", "_", colnames(concat.transformed.in), fixed = TRUE)
concat.transformed.in <- as.data.frame(concat.transformed.in)

#assign cell ID number to silhouettescore
cell_id <- 1:nrow(layout.in)
cell_id_df <- data.frame("cell_id" = cell_id)
layout.id <- cbind(layout.in,cell_id_df)
layout.clusters.id <- cbind(clusters.in,layout.id)

#assign cell ID number to concat_transform
cell_id <- 1:nrow(layout.in)
cell_id_df <- data.frame("cell_id" = cell_id)
layout.id <- cbind(layout.in,cell_id_df)
layout.concat.transform.id <- cbind(layout.id,concat.transformed.in)

print("needed files read, prepping data to threshold")

## Prep dataframe for thresholding
#thresholding.df <- merge(layout.clusters.id,layout.concat.transform.id)
#colnames(thresholding.df) <- c("cell_id","clusters","umap_x","umap_y","silhouettescore")
thresholding.ct.df <- merge(layout.clusters.id,layout.concat.transform.id)
thresholded.ct.df <- thresholding.ct.df[thresholding.ct.df$TuJ1_Y89Di <= EXPRESSION.VALUE, ]

print("data ready to write")

## Pull thresholded.ct.df apart into layout, cluster, and concat.transform

layout.out <- thresholded.ct.df[,c(1:2)]
clusters.out <- thresholded.ct.df[,c(4)]
concat.transformed.out <- thresholded.ct.df[,c(5:46)]

## Write output file ==============================================================================
fwrite(concat.transformed.out,file = CONCAT.OUT.FILENAME,row.names = FALSE)
fwrite(layout.out,file = LAYOUT.OUT.FILENAME)
fwrite(list(clusters.out),file = CLUSTERS.OUT.FILENAME,col.names = FALSE,sep = ",")