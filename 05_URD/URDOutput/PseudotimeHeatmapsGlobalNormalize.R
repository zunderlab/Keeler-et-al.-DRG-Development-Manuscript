#Corey Williams, University of Virginia
#20 Aug, 2020
#Make marker expression vs pseudotime plots
#INSTRUCTIONS: In panel file, add column, default name 'URD' and mark for inclusion in plot
#Warning: Sometimes timepoint heatmap will be shifted 1-2 pixels out of alignment with expression

rm(list = ls())
#.libPaths( c( .libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.6") )

library(dplyr)
library(reshape2)
library(ggfortify)
library(egg)
library(ggpubr)

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
OUTPUT.FOLDER.NAME <- "_pseudotime_heatmap_all"
URD_IN <- "/URD.RData"
PANEL_IN <- "/panel_URD_Heatmaps_FULL.csv"
URD_COLNAME <- "URD" #column name in panel for plotting markers
START_SEGMENT <-  #root
  
  #----------------changed to run on all tips
  #END_SEGMENT <- 1
  END_SEGMENTS <- c() #all tips of interest
CELL_TYPES <- c() #label with names
#PLOT_FILENAME <- "Pseudotime_Heatmap_tip1_SmoothMuscle.png"
PLOT_FILE_BASENAME <- "Pseudotime_Heatmap_tip__X__Y.png"
#----------------changed to run on all tips

PSEUDOTIME_BINS <- 1000 #number of columns in heatmap
NORMALIZE_GLOBALLY <- TRUE
NORMALIZE_LOCALLY <- FALSE
#Plot parameters
PLOT_WIDTH <- 22.5
PLOT_HEIGHT <- 17.5
LEGEND_WIDTH_FACTOR <- 50
TEXT_SIZE <- 11
MARKER_ORDER <- NULL
#MARKER_ORDER <- c(32,38,41,10,21,7,2,22,6,3,17,18,19,25,8,24,5,15,16,33,34,36,31,28,14,39,1,4,29,30,40,11,13,20,35,23,26,9,12,27,37) #to reorder, run to line 67 then look at urd_markers then (cont. below)
#use a vector like c(1,3,2) to reorder. 
#Will hierarchically cluster markers when null

## Check for user inputs ==========================================================================
if (is.null(START_SEGMENT) | is.null(END_SEGMENTS)){
  stop("Specify START_SEGMENT and END_SEGMENT")
}

## Load URD workspace =============================================================================
input.folder.temp <- INPUT.FOLDER
output.folder.temp <- OUTPUT.FOLDER
load(paste0(INPUT.FOLDER,URD_IN))
INPUT.FOLDER <- input.folder.temp
OUTPUT.FOLDER <- output.folder.temp
object <- URD_Object.tree


## Get cells on trajectory ========================================================================
#get order of segments from end to start

#----------------changed to run on all tips
cells_trajectory_list <- list()
for (END_SEGMENT in END_SEGMENTS) {
  #----------------changed to run on all tips
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
  cells_trajectory_list[[END_SEGMENT]] <- cells_trajectory
}

## Get markers for plotting =======================================================================
#read panel
panel <- read.panel(INPUT.FOLDER,PANEL_IN)
#get markers labeled for plotting
urd_markers <- as.character(paste0(panel$Antigen[panel[,URD_COLNAME] == 1],"_",
                                   panel$Metal[panel[,URD_COLNAME] == 1]))

## Make plotting data frame =======================================================================

heatmap_list <- list()
heatmap_concat <- data.frame(row.names=sort(urd_markers))
timepoint_list <- list()
for (END_SEGMENT in END_SEGMENTS) {
  cells_traj_current <- cells_trajectory_list[[END_SEGMENT]]
  #get pseudotime for each cell
  pseudotime_trajectory <- object@pseudotime[cells_traj_current,"pseudotime"]
  #get pseudotime quantiles
  pseudotime_quantiles <- quantile(pseudotime_trajectory,seq(0,1,length.out=PSEUDOTIME_BINS))
  #get pseudotime bin info
  pseudotime_bins <- cut(pseudotime_trajectory,unique(pseudotime_quantiles[1:(PSEUDOTIME_BINS)]),
                         include.lowest=TRUE)
  #get cells for each bin
  cell_ids_in_bin <- lapply(sort(unique(pseudotime_bins)),
                            function(this_bin){which(pseudotime_bins==this_bin)})
  # #get expression info for each bin
  # exprs_bin <- lapply(urd_markers,function(this_marker){
  #   data.frame(sort(unique(pseudotime_bins)),
  #              sapply(cell_ids_in_bin,function(this_bin){
  #                mean(object@logupx.data[this_marker,cells_traj_current[this_bin]])}),
  #              rep(this_marker,times = length(unique(pseudotime_bins))))
  # })
  # heatmap_df <- do.call(rbind,exprs_bin)
  # colnames(heatmap_df) <- c("pseudotime","expression","marker")
  # #normalize expression of markers (need for hierarchical clustering or for plotting)
  # if (NORMALIZE_GLOBALLY == TRUE) {
  #   for(this_marker in urd_markers){
  #     these_bins <- which(heatmap_df$marker == this_marker)
  #     heatmap_df$expression[these_bins] <- heatmap_df$expression[these_bins]/
  #       max(object@logupx.data[this_marker,])
  #   }  
  # } else if (NORMALIZE_LOCALLY == TRUE){
  #   for(this_marker in urd_markers){
  #     these_bins <- which(heatmap_df$marker == this_marker)
  #     heatmap_df$expression[these_bins] <- heatmap_df$expression[these_bins]/
  #       max(heatmap_df$expression[these_bins])
  #   }
  # }
  # #normalize markers for plotting
  # if (NORMALIZE_LOCALLY == TRUE){
  #   heatmap_df <- heatmap_normalized
  # }
  # 
  # heatmap_ind <- dcast(heatmap_df,marker ~ pseudotime,value.var = "expression")
  # 
  # # Remove first column which is measurement parameter name for some reason!
  # # Add back as row names
  # row_names_params <- heatmap_ind[,1]
  # heatmap_ind <- heatmap_ind[,2:ncol(heatmap_ind)]
  # row.names(heatmap_ind) <- row_names_params
  # 
  # heatmap_list[[END_SEGMENT]] <- heatmap_df
  # heatmap_concat <- cbind(heatmap_concat, heatmap_ind)
  # cat("Concatenated heatmap for tip", END_SEGMENT, "trajectory.\n")
  
  #get timepoint info
  timepoints <- sapply(cell_ids_in_bin,function(this_bin){
    #get timepoint values from metadata
    these_timepoints <- object@meta[cells_trajectory[this_bin],"stage"]
    #get frequencies of each timepoint value
    timepoint_count <- sapply(unique(these_timepoints),function(this_timepoint){
      sum(these_timepoints==this_timepoint)
    })
    #get most common timepoint value
    unique(these_timepoints)[which.max(timepoint_count)]
  })
  #make timepoint data frame, adding column of ones called "stage" for ggplot y axis
  timepoint_df <- data.frame(sort(unique(pseudotime_bins)),timepoints,rep(1,length(timepoints)))
  colnames(timepoint_df) <- c("pseudotime","stage","timepoint")
  
  timepoint_list[[END_SEGMENT]] <- timepoint_df
  cat("Added timepoints to list for tip", END_SEGMENT, "trajectory.\n")
}

# Export because it takes so long to process
save(heatmap_list, file="heatmap_list.RData")
save(heatmap_concat, file="heatmap_concat.RData")
save(timepoint_list, file="timepoint_list.RData")

#load("heatmap_list.RData")
#load("heatmap_concat.RData")
#load("timepoint_list.RData")



#do hierarchical clustering whenever marker order is not manually specified
if (is.null(MARKER_ORDER) == TRUE){
  #do hierarchical clustering to get order of markers
  
  #hclust_df <- dcast(heatmap_df,marker ~ pseudotime,value.var = "expression")
  hclust_df <- heatmap_concat
  
  hclust_dist <- dist(hclust_df,method = "euclidean")
  hclust_out <- hclust(hclust_dist,method = "ward.D")
  
  marker_order <- sort(urd_markers)[hclust_out$order]
} else {
  marker_order <- urd_markers[rev(MARKER_ORDER)]
}

for (i in 1:length(END_SEGMENTS)) {
  
  end_segment <- END_SEGMENTS[i]
  
  heatmap_df <- heatmap_list[[end_segment]]
  timepoint_df <- timepoint_list[[end_segment]]
  
  #reorder levels of plotting df based on hclust
  heatmap_df$marker <- factor(heatmap_df$marker,levels=marker_order)
  
  ## Plot ===========================================================================================
  #Make output directory
  #time.now <- Sys.time()
  #output.dir <- paste0(OUTPUT.FOLDER,"/",substr(time.now,start=1,stop=10),"_",
  #                     substr(time.now,start=12,stop=13),".",substr(time.now,start=15,stop=16),".",
  #                     substr(time.now,start=18,stop=19),OUTPUT.FOLDER.NAME)
  #dir.create(output.dir)
  #setwd(output.dir)
  
  #make expression heatmap
  exprs_plot <- ggplot(heatmap_df) +
    geom_tile(aes(x = pseudotime,y = marker,fill = expression)) +
    scale_fill_viridis_c() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = TEXT_SIZE,color = "black",hjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),
          axis.ticks = element_blank(),legend.title = element_text(size = TEXT_SIZE),
          legend.key.width = unit(PLOT_WIDTH/LEGEND_WIDTH_FACTOR,"in"),plot.margin = margin())
  
  #make timepoint heatmap
  timepoint_plot <- ggplot(timepoint_df) +
    geom_tile(aes(x = pseudotime,y = timepoint,fill = stage)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text = element_blank(),
          axis.title.y = element_text(size=TEXT_SIZE,hjust=1,vjust=0.5,angle=0),
          axis.title.x = element_text(size = TEXT_SIZE,margin = margin()),
          axis.ticks = element_blank(),legend.position = "none",plot.margin = margin())
  
  egg::ggarrange(exprs_plot,timepoint_plot,ncol = 1,heights = c(length(urd_markers),1))
  
  outfile_name <- sub("__X", end_segment, PLOT_FILE_BASENAME)
  outfile_name <- sub("_Y", CELL_TYPES[i], outfile_name)
  
  #save plot
  ggsave(outfile_name,
         egg::ggarrange(exprs_plot,timepoint_plot,ncol = 1,heights = c(length(urd_markers)-1,1)),width = PLOT_WIDTH,height = PLOT_HEIGHT)
  cat("Output", outfile_name)
}

#reset working directory
#setwd("..")
