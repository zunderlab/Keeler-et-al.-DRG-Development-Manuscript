#Corey Williams, University of Virginia
#20 Aug, 2020
#Make marker expression vs pseudotime plots
#INSTRUCTIONS: In panel file, add column, default name 'URD' and mark for inclusion in plot

rm(list = ls())
.libPaths( c( .libPaths(), "~/R/3.6.3_URD_Debug") )

library(dplyr)
library(reshape2)
library(ggfortify)
library(Matrix)

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
OUTPUT.FOLDER.NAME <- "_multiple_branches_"
URD_IN <- "/URD.RData"
PANEL_IN <- "panel_URD.csv"
URD_COLNAME <- "URD" #column name in panel for plotting markers
START_SEGMENT <-  #roots
  END_SEGMENT <- c() #multiple tips
PLOT_FILENAME <- "Markers_Pseudotime.png"
DF <- 5
LINE_WIDTH <- 2
DISCRETE_COLORS <- c()
PLOT_HEIGHT <- 7
PLOT_WIDTH <- 10

## Check for user inputs ==========================================================================
if (is.null(START_SEGMENT) | is.null(END_SEGMENT)){
  stop("Specify START_SEGMENT and END_SEGMENT")
}

## Load URD workspace =============================================================================
input.folder.temp <- INPUT.FOLDER
output.folder.temp <- OUTPUT.FOLDER
load(paste0(INPUT.FOLDER,URD_IN))
INPUT.FOLDER <- input.folder.temp
OUTPUT.FOLDER <- output.folder.temp
object <- URD_Object.tree

## Get cells on trajectories ======================================================================
segments_trajectory <- lapply(END_SEGMENT,function(this_tip){
  #get order of segments from end to start
  these_segments <- this_tip
  this_segment <- as.character(these_segments)
  do_segment <- TRUE
  while (do_segment == TRUE){
    this_segment <- dplyr::filter(object@tree$segment.joins,child == this_segment)$parent
    these_segments <- c(these_segments,as.integer(this_segment))
    if (this_segment == as.character(START_SEGMENT)){
      do_segment <- FALSE
    }
  }
  these_segments
}) %>%
  unlist() %>%
  unique() %>%
  sort()

#Get cells in each segment
cells_trajectory <- lapply(segments_trajectory,function(this_segment){
  object@tree$cells.in.segment[[toString(this_segment)]]
}) 

#attach segment IDs to cells
names(cells_trajectory) <- segments_trajectory

## Get markers for plotting =======================================================================
#read panel
panel <- read.panel(INPUT.FOLDER,PANEL.FILENAME)
#get markers labeled for plotting
urd_markers <- as.character(paste0(panel$Antigen[panel[,URD_COLNAME] == 1],"_",
                                   panel$Metal[panel[,URD_COLNAME] == 1]))

## Calculate splines ==============================================================================
#get pseudotime for each cell
cells_trajectory <- lapply(as.character(segments_trajectory),function(this_segment){
  df_out <- data.frame(t(as.matrix(object@logupx.data[urd_markers,
                                                      cells_trajectory[[this_segment]]])),
                       object@pseudotime[cells_trajectory[[this_segment]],"pseudotime"],
                       rep(this_segment,length(cells_trajectory[[this_segment]])))
  colnames(df_out)[(ncol(df_out)-1):ncol(df_out)] <- c("pseudotime","segment")
  df_out
})
#make df for plotting cells
names(cells_trajectory) <- as.character(segments_trajectory)
cells_plotting <- do.call(rbind,cells_trajectory)

#get splines for each marker
splines_trajectory <- lapply(urd_markers, function(this_marker){
  this_list <- lapply(as.character(segments_trajectory), function(this_segment){
    #check that there are cells in segment
    if(nrow(cells_trajectory[[this_segment]]) > 0){
      #Check if there is enough variation to calculate a spline
      if(IQR(cells_trajectory[[this_segment]]$pseudotime)>0){
        #calculate splines
        spline_out <- smooth.spline(x=cells_trajectory[[this_segment]]$pseudotime,
                                    y=cells_trajectory[[this_segment]][,this_marker],df=DF)
        spline_out <- data.frame(spline_out$x,spline_out$y,rep(this_segment,length(spline_out$x)))
        colnames(spline_out) <- c("pseudotime",this_marker,"segment")
        return(spline_out)
      } else {
        return(NULL)
      }
    }
  })
  #output list
  names(this_list) <- as.character(segments_trajectory)
  return(this_list)
})
#make spline list for plotting
splines_plotting <- lapply(splines_trajectory,function(this_marker){do.call(rbind,this_marker)})
names(splines_plotting) <- urd_markers

## Plot ===========================================================================================
#Make output directory
time.now <- Sys.time()
output.dir <- paste0(OUTPUT.FOLDER,"/",substr(time.now,start=1,stop=10),"_",
                     substr(time.now,start=12,stop=13),".",substr(time.now,start=15,stop=16),".",
                     substr(time.now,start=18,stop=19),OUTPUT.FOLDER.NAME)
dir.create(output.dir)
setwd(output.dir)

#plot
#plot just lines
lapply(urd_markers,function(this_marker){
  ggsave(paste0(this_marker,".png"),ggplot() + 
           geom_line(data = splines_plotting[[this_marker]],
                     aes_string(x="pseudotime",y=this_marker,color="segment"),
                     size = LINE_WIDTH) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                 axis.text = element_blank(),axis.title.y = element_text(size = 8),
                 axis.title.x = element_text(size = 8))) 
})
#plot lines and points
lapply(urd_markers,function(this_marker){
  ggsave(paste0(this_marker,"_with_lines.png"),ggplot() + 
           geom_point(data = cells_plotting,
                      aes_string(x = "pseudotime",y = this_marker,color = "segment"),
                      alpha = 0.1) + 
           geom_line(data = splines_plotting[[this_marker]],
                     aes_string(x="pseudotime",y=this_marker,color="segment"),
                     size = LINE_WIDTH) +
           scale_color_manual(values=DISCRETE_COLORS) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                 axis.text = element_blank(),axis.title.y = element_text(size = 8),
                 axis.title.x = element_text(size = 8)),
         height = PLOT_HEIGHT,width = PLOT_WIDTH)
})

#reset working directory
setwd("..")