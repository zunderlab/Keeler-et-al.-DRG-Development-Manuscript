#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start PlotByCluster.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(ggfortify)
library(dplyr)

print("libraries loaded")

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
METADATA.FILENAME <- "/metadata.csv"
CONCAT.TRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
CLUSTER.PLOT.FILENAME <- "clusters_.png" #label for condition of interest
POINT.SIZE <- 0.1
GRAY.COLOR <- "#bdbdbd"
METADATA.TAG <- "genotype"#header of metadata column to be grayed out
METADATA.GRAY <- c("NA") #Which metadata conditions to gray. Use string or integer to index

print("input parameters loaded, reading needed files")

## Read needed files ==============================================================================
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)
metadata <- read.metadata(INPUT.FOLDER,METADATA.FILENAME)
#This is inefficient, but currently the best option I have to link file # to cells in other files
concat.transformed <- read.concat.transformed(INPUT.FOLDER,CONCAT.TRANSFORMED.FILENAME) 

print("needed files read, prepping data to plot")

## Prep dataframe for plotting ====================================================================
plotting.df <- as.data.frame(cbind(layout.in,clusters.in,concat.transformed$File))
colnames(plotting.df) <- c("umap_x","umap_y","cluster","File")
#Get file IDs for graying
if (is.character(METADATA.GRAY)==1){
  gray.files <- which(metadata[,METADATA.TAG] %in% METADATA.GRAY)
} else if (is.numeric(METADATA.GRAY) == 1){
  gray.files <- which(metadata[,METADATA.TAG] %in% unique(metadata[,METADATA.TAG])[METADATA.GRAY])
} else {
  print("METADATA.GRAY in unrecognized format")
  quit()
}

print("data ready to plot, plotting")

## Save plots colored by each marker ==============================================================
#set output folder
setwd(OUTPUT.FOLDER)
#loop through variables to plot
ggsave(CLUSTER.PLOT.FILENAME,plot = ggplot() + 
         geom_point(data = dplyr::filter(plotting.df,File %in% gray.files),aes(x=umap_x,y=umap_y),
                    color = GRAY.COLOR, size = POINT.SIZE) + 
         geom_point(data = dplyr::filter(plotting.df,!(File %in% gray.files)),
                    aes(x=umap_x,y=umap_y,color=factor(cluster)), size = POINT.SIZE) +
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black")) +
         guides(colour = guide_legend(override.aes = list(shape=15, size=8))),
       #^^should work for changing size/shape of legend elements... 
       #might have to tweak size per preference
       height = 7,width = 7)

print("data plotted and file outputted")
print("End PlotByClusterGrayMetadata.R")