#Corey Williams, University of Virginia
#16 Feb, 2020
#Make metadata csv for URD
#Sets up metadata to have columns for file, genotype, and stage of each cell

print("Start GetUrdMetadata.R")

rm(list = ls())
.libPaths( c( .libPaths(), "/project/zunderlab/R/3.6") )

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

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
STAGE.TAG <- "stage" #column name in metadata for developmental stage
METADATA.FILENAME <- "/metadata.csv" #Must include "stage" column for URD
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv" #.csv file for pipeline clustering
URD.META.FILENAME <- "metadata_URD.csv" #.txt file For URD - unused if no subsampling

print("got input parameters, loading input files")

## Read input files ===============================================================================
md <- read.metadata(INPUT.FOLDER,METADATA.FILENAME)
concat.transformed <- read.concat.transformed(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME)

print("finished reading input files, getting metadata")

## Get standard metadata ==========================================================================
#initialize array
urd.md <- matrix(0,nrow = nrow(concat.transformed),ncol = 3)
#get metadata for each cell
urd.md <- md[concat.transformed[,"File"],c("file_name",STAGE.TAG)]

## Get Root identities ============================================================================
if ("Root_clusters" %in% colnames(md)){
  
  #Get file ID(s) (as corresponds to metadata and concat transformed) that have roots
  root_files <- which(md$Root_clusters != "")
  
  #Get root cluster filenames if name independently. 
  #Default is name of fcs file (sans .fcs) + "_clusters.csv", so cells.fcs is cells_clusters.csv
  root_filenames <- vector(mode = "character",length = length(root_files))
  for (this_file in seq_along(root_files)) {
    #Blank filenames will be NA if all left blank, will be blank string if only one left blank
    if (md$Root_cluster_filename[root_files[this_file]] == "" | 
        is.na(md$Root_cluster_filename[root_files[this_file]]) == TRUE) {
      #Generate cluster filename
      this_filename <- md$file_name[root_files[this_file]]
      root_filenames[this_file] <- paste0(substr(this_filename,start=1,
                                                 stop=nchar(this_filename)-4),
                                          "_clusters.csv")
    } else {
      #Give inputted filename
      root_filenames[this_file] <- as.character(md$Root_cluster_filename[root_files[this_file]])
    }
  }
  
  #Read in clusters from root files
  root_clusters_in <- lapply(root_filenames, function(this_file){
    read.clusters(INPUT.FOLDER,paste0("/",this_file))})
  
  #Get root cluster ID's from metadata
  root_clusters <- lapply(root_files, function(this_file){
    md$Root_clusters[this_file] %>%
      as.character %>%
      strsplit(",") %>% 
      unlist() %>% 
      as.numeric
  })
  
  #Identify cells in each fcs file corresponding to the root cluster ID's
  root_cells <- lapply(seq_along(root_clusters_in),function(this_file){
    root_clusters_in[[this_file]] %in% root_clusters[[this_file]]
  })
  
  #Make full vector of logicals to include all cells
  root_cells_out <- vector(mode = "logical",length = nrow(urd.md))
  for (this_file in seq_along(root_files)) {
    root_cells_out[which(concat.transformed$File==root_files[this_file])]<-root_cells[[this_file]]
  }
  
  #append root info to metadata
  urd.md <- data.frame(urd.md,root_cells_out)
  colnames(urd.md)[ncol(urd.md)] <- "roots"
}

## Get tip identities =============================================================================
if ("Tip_clusters" %in% colnames(md)){
  
  #Get file ID(s) (as corresponds to metadata and concat transformed) that have tips
  tip_files <- which(md$Tip_clusters != "")
  
  #Get tip cluster filenames if name independently. 
  #Default is name of fcs file (sans .fcs) + "_clusters.csv", so cells.fcs is cells_clusters.csv
  tip_filenames <- vector(mode = "character",length = length(tip_files))
  for (this_file in seq_along(tip_files)) {
    #Blank filenames will be NA if all left blank, will be blank string if only one left blank
    if (md$Tip_cluster_filename[tip_files[this_file]] == "" | 
        is.na(md$Tip_cluster_filename[tip_files[this_file]]) == TRUE) {
      #Generate cluster filename
      this_filename <- md$file_name[tip_files[this_file]]
      tip_filenames[this_file] <- paste0(substr(this_filename,start=1,
                                                stop=nchar(this_filename)-4),
                                         "_clusters.csv")
    } else {
      #Give inputted filename
      tip_filenames[this_file] <- as.character(md$Tip_cluster_filename[tip_files[this_file]])
    }
  }
  
  #Read in clusters from tip files
  tip_clusters_in <- lapply(tip_filenames, function(this_file){
    read.clusters(INPUT.FOLDER,paste0("/",this_file))})
  
  #Get tip cluster ID's from metadata
  tip_clusters <- lapply(tip_files, function(this_file){
    md$Tip_clusters[this_file] %>%
      as.character %>%
      strsplit(",") %>% 
      unlist() %>% 
      as.numeric
  })
  
  #Identify cells in each fcs file corresponding to the tip cluster ID's
  cluster_count <- 0
  tip_cells <- list()
  for (this_file in seq_along(tip_clusters_in)){
    tip_cell_clusters <- vector(mode = "integer",length = length(tip_clusters_in[[this_file]]))
    
    #find which cells correspond to which cluster
    for (this_cluster in seq_along(tip_clusters[[this_file]])){
      tip_cell_clusters[which(tip_clusters_in[[this_file]] == tip_clusters[[this_file]][this_cluster])] <- this_cluster + cluster_count
    }
    tip_cells[[this_file]] <- tip_cell_clusters
    
    #keep track of clusters recorded so that numbering is coordinated across multiple files
    cluster_count <- cluster_count + this_cluster
  }
  
  #Make full vector of logicals to include all cells
  tip_cells_out <- vector(mode = "integer",length = nrow(urd.md))
  for (this_file in seq_along(tip_files)) {
    tip_cells_out[which(concat.transformed$File == tip_files[this_file])] <- tip_cells[[this_file]]
  }
  
  #append tip info to metadata
  urd.md <- data.frame(urd.md,tip_cells_out)
  colnames(urd.md)[ncol(urd.md)] <- "tips"
}

## Write urd metadata file ========================================================================
fwrite(urd.md,URD.META.FILENAME)