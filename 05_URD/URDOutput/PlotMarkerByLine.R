#Corey Williams, University of Virginia
#20 Aug, 2020
#Make marker expression vs pseudotime plots
#INSTRUCTIONS: In panel file, add column, default name 'URD' and mark for inclusion in plot

rm(list = ls())
.libPaths( c( .libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.6") )

library(dplyr)
library(ggfortify)

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
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
OUTPUT.FOLDER.NAME <- "_pseudotime_MarkerByLine_"
URD_IN <- "/URD.RData"
PANEL_IN <- "panel_URD.csv"
URD_COLNAME <- "URD" #column name in panel for plotting markers
START_SEGMENT <-  #Root or parent segment
  END_SEGMENT <-  #tip
  PLOT_FILENAME <- "Markers_Pseudotime_.png"

## Load URD workspace =============================================================================
input.folder.temp <- INPUT.FOLDER
output.folder.temp <- OUTPUT.FOLDER
load(paste0(INPUT.FOLDER,URD_IN))
INPUT.FOLDER <- input.folder.temp
OUTPUT.FOLDER <- output.folder.temp
object <- URD_Object.tree

## Get cells on trajectory ========================================================================
#get order of segments from end to start
segments_trajectory <- END_SEGMENT
this_segment <- as.character(segments_trajectory)
while (this_segment %in% object@tree$segment.joins$child){
  this_segment <- dplyr::filter(object@tree$segment.joins,child == this_segment)$parent
  segments_trajectory <- c(segments_trajectory,as.integer(this_segment))
}
#Get cells in each segment
cells_trajectory <- unlist(lapply(segments_trajectory,function(this_segment){
  object@tree$cells.in.segment[[toString(this_segment)]]
}))

## Get markers for plotting =======================================================================
#read panel
panel <- read.panel(INPUT.FOLDER,PANEL.FILENAME)
#get markers labeled for plotting
urd_markers <- as.character(paste0(panel$Antigen[panel[,URD_COLNAME] == 1],"_",panel$Metal[panel[,URD_COLNAME] == 1]))

## Calculate splines ==============================================================================
#get pseudotime for each cell
pseudotime_trajectory <- object@pseudotime[cells_trajectory,"pseudotime"]
#get expression data for selected markers
exprs_trajectory <- t(as.matrix(object@logupx.data[urd_markers,cells_trajectory]))
#get splines for each marker
splines_trajectory <- lapply(urd_markers,function(this_marker){
  smooth.spline(pseudotime_trajectory,exprs_trajectory[,this_marker],df=5)
})

## Plot ===========================================================================================
#Make output directory
time.now <- Sys.time()
output.dir <- paste0(OUTPUT.FOLDER,"/",substr(time.now,start=1,stop=10),"_",
                     substr(time.now,start=12,stop=13),".",substr(time.now,start=15,stop=16),".",
                     substr(time.now,start=18,stop=19),OUTPUT.FOLDER.NAME)
dir.create(output.dir)
setwd(output.dir)

#Generate df for plotting splines
splines_markers <- lapply(seq_along(urd_markers),function(this_marker){data.frame(splines_trajectory[[this_marker]]$x,splines_trajectory[[this_marker]]$y,rep(urd_markers[this_marker],length(splines_trajectory[[this_marker]]$x)))})
plotting_df_1 <- do.call(rbind,splines_markers)
colnames(plotting_df_1) <- c("pseudotime","expression","marker")
#output all spline curves on one plot
ggsave("splines.png",ggplot(plotting_df_1) + geom_line(aes(x=pseudotime,y=expression,color=marker)) 
       + scale_color_brewer(type = "qual",palette = 2) + 
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"), 
               axis.text = element_blank(),axis.title.y = element_text(size = 8),
               axis.title.x = element_text(size = 8)))

#generate df for plotting points
plotting_df_2 <- data.frame(pseudotime_trajectory,as.matrix(exprs_trajectory))
colnames(plotting_df_2)[1] <- "pseudotime"
#output points and curve on plot
sapply(urd_markers,function(this_marker){
  ggsave(paste0(gsub("[.]","_",this_marker),".png"),ggplot() + geom_point(data = plotting_df_2,aes_string(x = "pseudotime",y = this_marker,alpha = 0.1)) + geom_line(data = dplyr::filter(plotting_df_1,marker == this_marker),aes(x = pseudotime,y=expression,size=2)) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                 legend.position = "none",axis.text = element_blank(),axis.title.y = element_text(size = 8),
                 axis.title.x = element_text(size = 8)))
})

#reset working directory
setwd("..")