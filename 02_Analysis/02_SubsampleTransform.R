#Corey Williams, University of Virginia
#17 May, 2019
#Subsample and make csv file containing all points asinh transformed

print("Start subsampleTransform.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(flowCore)
library(data.table)

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
METADATA.FILENAME <- "/metadata.csv"
PANEL.FILENAME <- "/panel.csv"
SUBSAMPLE.ID.FILENAME <- "Subsample_Ids.csv" #unused if no subsampling
CONCATTRANSFORMED.FILENAME <- "Concat_Transformed.csv"
ASINH.FACTOR <- NULL #set to NULL if set in panel
SUBSAMPLES <- NULL #change this value if you want subsampling. Number of cells per file

print("got input parameters, loading fcs files")

## Read metadata ==================================================================================
md <- read.metadata(INPUT.FOLDER,METADATA.FILENAME)

## Read fcs files =================================================================================
fcs.in <- read.flowSet(md$file_name, path=INPUT.FOLDER, transformation = FALSE,
                       truncate_max_range = FALSE)

print("loaded fcs files")

## Read panel =====================================================================================
panel <- read.panel(INPUT.FOLDER, PANEL.FILENAME)

## Pull out markers for transformation ============================================================
#Transform all markers marked for clustering or plotting
transform.markers <- get.transform.markers(panel)
transform.markers.annotate <- get.transform.annotate(panel)

## Get asinh factors
if (is.null(ASINH.FACTOR) == TRUE) {
  asinh.factor <- panel$asinh.factor[panel$Metal %in% transform.markers]
} else {
  asinh.factor <- ASINH.FACTOR
}

##Subsample & make concatenated exprs arrays ======================================================
if (is.numeric(SUBSAMPLES) == TRUE){
  #Get subsample id's for each file
  subsample.ids <- matrix(0,nrow = length(fcs.in),ncol = SUBSAMPLES)
  subsample.ids <- fsApply(fcs.in,function(x) sample(1:nrow(x),SUBSAMPLES))
  #add file to left of subsample id's for clarity
  rownames(subsample.ids) <- md$file_name
  #Get concatenated exprs array
  concat.exprs.list <- lapply(seq_along(fcs.in),function(x) 
    cbind(exprs(fcs.in[[x]])[subsample.ids[x,],transform.markers],rep(x,SUBSAMPLES)))
  concat.exprs <- do.call(rbind,concat.exprs.list)
  #give column names
  colnames(concat.exprs) <- c(transform.markers.annotate,"File")
  #perform asinh transform
  concat.exprs[,transform.markers.annotate] <- asinh(t(sapply(
    seq_along(1:nrow(concat.exprs)), 
    function(x){concat.exprs[x,transform.markers.annotate]/asinh.factor})))
  #save concatenated exprs array
  fwrite(subsample.ids,file = SUBSAMPLE.ID.FILENAME,col.names = FALSE,sep = ",")
  fwrite(concat.exprs,file = CONCATTRANSFORMED.FILENAME,row.names = FALSE)
  
  print("saved subsample ID file & concatenated and transformed file")
} else {
  #Get concatenated exprs array
  concat.exprs.list <- lapply(seq_along(fcs.in),function(x) 
    cbind(exprs(fcs.in[[x]])[,transform.markers],rep(x,nrow(fcs.in[[x]]))))
  concat.exprs <- do.call(rbind,concat.exprs.list)
  #give column names
  colnames(concat.exprs) <- c(transform.markers.annotate,"File")
  #perform asinh transform
  concat.exprs[,transform.markers.annotate] <- asinh(t(sapply(
    seq_along(1:nrow(concat.exprs)), 
    function(x){concat.exprs[x,transform.markers.annotate]/asinh.factor})))
  #save concatenated exprs array
  fwrite(concat.exprs,file = CONCATTRANSFORMED.FILENAME,row.names = FALSE)
  
  print("saved concatenated and transformed file (no subsampling)")
}