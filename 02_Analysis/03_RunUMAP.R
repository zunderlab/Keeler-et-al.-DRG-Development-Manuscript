#Corey Williams, University of Virginia
#02 Apr, 2019
#Run UMAP dimensionality reduction

print("Start RunUMAP.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(umap)
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
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
PANEL.FILENAME <- "/panel.csv"
UMAP.LAYOUT.FILENAME <- "UMAP_layout.csv" 
UMAP.INDEXES.FILENAME <- "UMAP_knn_indexes.csv"
UMAP.DISTANCES.FILENAME <- "UMAP_knn_distances.csv"
UMAP_N_NEIGHBORS <- 15
UMAP_MIN_DIST <- 0.000001 #UMAP default is 0.1

#Advanced settings
UMAP_N_COMPONENTS <- 2
UMAP_METRIC <- "euclidean"
UMAP_N_EPOCHS <- 1000
UMAP_INPUT <- "data"
UMAP_INIT <- "spectral"
UMAP_SET_OP_MIX_RATIO <- 1
UMAP_LOCAL_CONNECTIVITY <- 1
UMAP_BANDWIDTH <- 1
UMAP_ALPHA <- 1
UMAP_GAMMA <- 1
UMAP_NEGATIVE_SAMPLE_RATE <- 5
UMAP_A <- NA
UMAP_B <- NA
UMAP_SPREAD <- 1
UMAP_RANDOM_STATE <- 1
UMAP_TRANSFORM_STATE <- NA
UMAP_KNN_REPEATS <- 1
UMAP_VERBOSE <- FALSE
UMAP_LEARN_ARGS <- NA

print("got input parameters, loading data files")

##Read Concat_Transformed =========================================================================
concat.exprs <- fread(paste0(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME), stringsAsFactors = F)
colnames(concat.exprs) <- gsub("-", "_", colnames(concat.exprs))
colnames(concat.exprs) <- gsub(".", "_", colnames(concat.exprs), fixed = TRUE)
concat.exprs <- as.data.frame(concat.exprs)

##Pull out clustering variables ===================================================================
panel <- read.panel(INPUT.FOLDER,PANEL.FILENAME)
clustering.vars <- get.clustering.annotate(panel)
umap.in <- as.matrix(concat.exprs[,clustering.vars])

print("files loaded, running UMAP")

##Run UMAP ========================================================================================
#Set up UMAP settings
umap.settings <- umap.defaults
umap.settings$n_neighbors <- UMAP_N_NEIGHBORS
umap.settings$n_components <- UMAP_N_COMPONENTS
umap.settings$metric <- UMAP_METRIC
umap.settings$n_epochs <- UMAP_N_EPOCHS
umap.settings$input <- UMAP_INPUT
umap.settings$init <- UMAP_INIT
umap.settings$min_dist <- UMAP_MIN_DIST
umap.settings$set_op_mix_ratio <- UMAP_SET_OP_MIX_RATIO
umap.settings$local_connectivity <- UMAP_LOCAL_CONNECTIVITY
umap.settings$bandwidth <- UMAP_BANDWIDTH
umap.settings$alpha <- UMAP_ALPHA
umap.settings$gamma <- UMAP_GAMMA
umap.settings$negative_sample_rate <- UMAP_NEGATIVE_SAMPLE_RATE
umap.settings$a <- UMAP_A
umap.settings$b <- UMAP_B
umap.settings$spread <- UMAP_SPREAD
umap.settings$random_state <- UMAP_RANDOM_STATE
umap.settings$transform_state <- UMAP_TRANSFORM_STATE
umap.settings$knn_repeats <- UMAP_KNN_REPEATS
umap.settings$verbose <- UMAP_VERBOSE
umap.settings$umap_learn_args <- UMAP_LEARN_ARGS
#Run UMAP
umap.out <- umap(umap.in,config = umap.settings)

print("ran UMAP, outputting files")

## Save UMAP layout and knn graph =================================================================
fwrite(umap.out$knn$indexes,file = UMAP.INDEXES.FILENAME,row.names = FALSE,col.names = FALSE,
       sep = ",")
fwrite(umap.out$knn$distances,file = UMAP.DISTANCES.FILENAME,row.names = FALSE,col.names = FALSE,
       sep = ",")
colnames(umap.out$layout) <- c("umap_x","umap_y")
fwrite(umap.out$layout,file = UMAP.LAYOUT.FILENAME,row.names = FALSE,col.names = TRUE, sep = ",")

print("UMAP files outputted")
print("Finish RunUMAP.R")