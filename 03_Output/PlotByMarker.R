#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ggfortify)
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

## Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
OUTPUT.FOLDER.NAME <- "_Markers"
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
PANEL.FILENAME <- "/panel.csv"
POINT.SIZE <- 0.1
PLOT.HEIGHT <- 7
PLOT.WIDTH <- 7
#vary values in line 44 scale_color_viridis_c(values=) to alter color expressions

## Read needed files
concat.transformed <- fread(paste0(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME), stringsAsFactors = F)
colnames(concat.transformed) <- gsub("-", "_", colnames(concat.transformed))
colnames(concat.transformed) <- gsub(".", "_", colnames(concat.transformed), fixed = TRUE)
concat.transformed <- as.data.frame(concat.transformed)
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
panel <- read.panel(INPUT.FOLDER,PANEL.FILENAME)

## Prep dataframe for plotting
plotting.vars <- get.plotting.annotate(panel)
plotting.df <- as.data.frame(cbind(concat.transformed[,plotting.vars],layout.in))

## Save plots colored by each marker
#make output folder
time.now <- Sys.time()
output.dir <- paste0(OUTPUT.FOLDER,"/",substr(time.now,start=1,stop=10),"_",
                     substr(time.now,start=12,stop=13),".",substr(time.now,start=15,stop=16),".",
                     substr(time.now,start=18,stop=19),OUTPUT.FOLDER.NAME)
dir.create(output.dir)
setwd(output.dir)
#loop through variables to plot
for (var.to.plot in plotting.vars){
  ggsave(paste0(var.to.plot,".png"),plot = ggplot(plotting.df,aes_string(x="umap_x",y="umap_y",
                                                                         color=var.to.plot)) + 
           geom_point(size = POINT.SIZE) + scale_color_viridis_c(values = c(0,0.1,0.3,1)) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")),
         height = PLOT.HEIGHT, width = PLOT.WIDTH)
}