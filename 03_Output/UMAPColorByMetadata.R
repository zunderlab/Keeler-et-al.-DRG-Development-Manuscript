# Created by Austin Keeler & Ashley Hirt
# 15 July 2020
# Plot a single UMAP colored by a metadata variable

library(tidyverse)
library(scales)
library(data.table)
library(viridis)


# ONLY WORKS IF SAME # CELLS IN UMAP.LAYOUT AS CONCAT.TRANSFORMED.
# ALSO TO ENSURE YOU GET COLD TO HOT COLORING THE WAY YOU WANT, MAKE SURE 
# THE ROWS ARE ARRANGED SO COLDER LEVELS AT THE TOP (EG EARLY TIMEPOINTS) 
# AND HOTTER LEVELS ARE AT THE BOTTOM. THE ROW ORDERING WILL DETERMINE
# THE COLOR PROGRESSION.

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

# Input parameters ----
# make sure to set your working directory!!
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata.csv"
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
POINT.SIZE <- 0.0001
VAR <- "age" # make sure this matches your md column name exactly
OUTPUT.FILENAME <- paste0("Colored_by_", VAR, ".png") 
# ^^ this is just a default, you can change--supply a full string, eg "My_File_Name.png")


# Read files ----
# read metadata
md <- fread(paste0(INPUT.FOLDER, "/", METADATA.FILENAME), stringsAsFactors = F)
# assign file numbering
md <- cbind(md, "File" = 1:nrow(md)) 

# read concat.transformed (we're just using the File col)
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCATTRANSFORMED.FILENAME), 
                            stringsAsFactors = F)

# Read UMAP layout, has two columns, x and y, where each row is a cell
umap_layout <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
# this is terrible; combine concat File column w/ umap df
# this assumes same # of cells and same "ordering" of cells (first at the top)
plotting_df <- cbind(umap_layout, concat.transformed$File)
# Assign column names
colnames(plotting_df) <- c("umap_x", "umap_y", "File")


# Join and Plot ----

# join on File 
plotting.df <- merge(plotting_df, md)
set.seed(42)
random.plotting.df <- plotting.df[sample(nrow(plotting.df)),]

# save plot
ggsave(OUTPUT.FILENAME, 
       plot = {
         ggplot(data=random.plotting.df, aes_string(x="umap_x", y="umap_y", color=VAR)) + 
           # note color=VAR, that determines coloring
           geom_point(size = POINT.SIZE) + 
           scale_color_viridis_d(option = "plasma") + 
           theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")
           )
       },
       device = "png", height = 14,width = 21, units = "in", dpi = 600 
)