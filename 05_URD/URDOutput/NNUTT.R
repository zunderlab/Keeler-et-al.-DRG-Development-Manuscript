#Originally created by Emily Puleo. Modified on 10/19/21 by Austin Keeler


rm(list = ls())

# Load packages
#library(tidyverse)
library(ggfortify)
library(data.table)
library(URD)
library(ggplot2)
library(dplyr)
library(tidyr)

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


# Define constants
INPUT_FOLDER <- getwd()
WRITE_PATH<-getwd()
# these are the cells you want to map
RA_FILENAME1 <- "Concat_Transformed_1.csv" #ie, dataset1
#check segmentlist on line 88 and change as needed

#this is the URD you are using
IV_FILENAME <- "Concat_Transformed.csv" #ie, dataset2
URD_IN <- "/URD.RData"
#^open the URD workspace into the environment manually

# define distance function
euc_dist <- function(row1, row2) {
  row1 <- as.numeric(row1)
  row2 <- as.numeric(row2)
  # make sure lengths are equal
  if (length(row1) == length(row2)) { 
    num_params <- length(row1)
  } else {
    stop("Length of row1 does not equal length of row2.")
  }
  ind_sum <- vector(length = num_params)
  for (param_i in seq_len(num_params)) { 
    #iterate over parameters, store each component of dist in ind_sum
    ind_sum[param_i] <- (row1[param_i] - row2[param_i])^2
  }
  dist <- sqrt(sum(ind_sum))
  return(dist)
}
# Read in files =====================================================================================
load(paste0(INPUT_FOLDER,URD_IN))
dataset_RA <- fread(file.path(INPUT_FOLDER, RA_FILENAME1), 
                    stringsAsFactors = F,
                    data.table = F)

iv_dataset <- fread(file.path(INPUT_FOLDER, IV_FILENAME),
                    stringsAsFactors = F,
                    data.table = F)

##FILE 1 ======================================
#isolate file
expr_RA <- dataset_RA

# combine urd expr data w/ comp tSNE layout coords
iv_dataset <- cbind(iv_dataset, URD_Object@tsne.y$tSNE1,URD_Object@tsne.y$tSNE2)
iv_expr <- select(iv_dataset, -c("File","URD_Object@tsne.y$tSNE1","URD_Object@tsne.y$tSNE2"))

#make sure panel is the same
expr_RA<-expr_RA[,colnames(expr_RA)==colnames(iv_expr)]
iv_expr<-iv_expr[,colnames(expr_RA)==colnames(iv_expr)]
colnames(expr_RA)==colnames(iv_expr)
#^ check to make sure all TRUE

#make dists vector
dists <- vector(mode = "list", length = nrow(expr_RA))

#calc dist
for (i in 1:nrow(expr_RA)) {
  cell <- expr_RA[i,]
  dists[[i]] <- apply(iv_expr, 1, function(comp_cell) {
    euc_dist(comp_cell, cell)
  }
  )
}


tsne_ind<-rep(0,length(dists))
for (i in 1:length(dists)){
  tsne_ind[i]<-which(dists[[i]]==min(dists[[i]]))
}

segmentlist<-c() #list all segments in your data
new_index<-matrix(data=NA, nrow=length(tsne_ind),ncol=1)
for (k in 1:length(segmentlist)){
  for (i in 1:length(tsne_ind)){
    for (j in 1:length(URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]])){
      if (paste0("V",tsne_ind[i]) == URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]][j]){
        new_index[i]<-tsne_ind[i]
        #which(paste0("V",tsne_ind[i]) == URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]][j])]
      }
    }}}

new_index<-data.frame(new_index)
new_index %>% drop_na()
new_index<-as.list(new_index)
#new_index<-as.matrix(new_index)
new_index<-as.numeric(unlist(new_index))

#plotURD
similar.cells1<-URD_Object@group.ids$init[new_index] #most similar points
#URD_Object.tree.new<-urdSubset(URD_Object, cells.keep=whichCells(URD_Object,"init",similar.cells))
png("NNUTT.png")
plotTree(URD_Object.tree,"Stage",cells.highlight=whichCells(URD_Object,"init",similar.cells1), cell.alpha=0, cells.highlight.alpha = 1, cells.highlight.size=2)
dev.off()