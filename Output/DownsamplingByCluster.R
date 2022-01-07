#Corey Williams, University of Virginia
#14 May, 2020
#Downsample files to equally represent clusters

print("Start downsampleByCluster.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(data.table)
library(dplyr)

## Input parameters ===============================================================================
SUBSAMPLES <- 1000 #Number of cells per cluster. If a cluster contains less, all cells will be kept

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
DOWNSAMPLE.OUT.FILENAME <- "Concat_Transformed_.csv" #rename
SUBSAMPLE.ID.OUT.FILENAME <- "Subsample_ID_by_Cluster.csv" #Relative to concat_transformed file, rename
LAYOUT.OUT.FILENAME <- "UMAP_layout_Subsampled.csv" #rename
CLUSTERS.OUT.FILENAME <- "clusters_Subsampled.csv" #will be needed to cooperate with pipeline, rename

print("Finished reading input parameters, reading files")

## Read needed files ==============================================================================
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCATTRANSFORMED.FILENAME), stringsAsFactors = F)
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)

print("Finished reading needed files, downsampling")

## Perform cluster-based downsampling =============================================================
#Initialize lists for outputs
numcluster <- length(unique(clusters.in))
subsample.ids <- vector(mode = "list",length = numcluster)
expression.out <- subsample.ids
layout.out <- subsample.ids
clusters.out <- subsample.ids

#loop through each cluster
for (this_cluster in sort(unique(clusters.in))) {
  
  #Check if subsampling is needed or if all cells for cluster will be used
  if (sum(clusters.in == this_cluster) < SUBSAMPLES) {
    #get ids of cells in cluster file
    subsample.ids[[this_cluster]] <- which(clusters.in == this_cluster)
    #use subsample ids to subsample from expression info
    expression.out[[this_cluster]] <- concat.transformed[subsample.ids[[this_cluster]],]
    # store layout
    layout.out[[this_cluster]] <- layout.in[subsample.ids[[this_cluster]],]
    #store cluster id's
    clusters.out[[this_cluster]] <- clusters.in[subsample.ids[[this_cluster]]]
    
  } else {
    #Get vector of all ids for this cluster
    all.ids <- which(clusters.in == this_cluster)
    #randomly subsample from the selection of cells
    subsample.ids[[this_cluster]] <- sample(all.ids,size = SUBSAMPLES)
    #use subsample ids to subsample from expression info
    expression.out[[this_cluster]] <- concat.transformed[subsample.ids[[this_cluster]],]
    #store subsampled layout
    layout.out[[this_cluster]] <- layout.in[subsample.ids[[this_cluster]],]
    #store cluster id's
    clusters.out[[this_cluster]] <- clusters.in[subsample.ids[[this_cluster]]]
    
  }
}

#Move outputs from list format to vector or array to output
subsample.ids.out <- unlist(subsample.ids)
expression.out <- rbindlist(expression.out)
layout.out <- rbindlist(layout.out)
clusters.out <- unlist(clusters.out)

print("Finished downsampling, writing output file")

## Write output file ==============================================================================
fwrite(list(subsample.ids.out),file = SUBSAMPLE.ID.OUT.FILENAME,col.names = FALSE,sep = ",")
fwrite(expression.out,file = DOWNSAMPLE.OUT.FILENAME,row.names = FALSE)
fwrite(layout.out,file = LAYOUT.OUT.FILENAME)
fwrite(list(clusters.out),file = CLUSTERS.OUT.FILENAME,col.names = FALSE,sep = ",")

print(table(clusters.out))

print("Completed downsampleByCluster.R")
