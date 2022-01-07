#Austin Keeler, University of Virginia
#5 Dec, 2021
#subset cells based on segment

print("Start ExtractClusters.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(data.table)

## Input parameters ===============================================================================
SEGMENTs_KEEP <- c() #comma separated segments if more than one

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
SEGMENTS.FILENAME <- "/segments.csv"
PSEUDO.SEGMENTS.FILENAME <- "/pseudo_segments.csv"
PSEUDO.SEGMENTS.OUT.FILENAME <- "pseudo_segments_.csv" #rename to include segments kept

print("Finished reading input parameters, reading files")

## Read needed files ==============================================================================
segments.in <- fread(paste0(INPUT.FOLDER,SEGMENTS.FILENAME), stringsAsFactors = F)
segments.in <- as.data.frame(segments.in)
segments.in <- segments.in[which(rowSums(segments.in) > 0),]
pseudo.segments.in <- fread(paste0(INPUT.FOLDER,PSEUDO.SEGMENTS.FILENAME), stringsAsFactors = F)
pseudo.segments.in <- as.data.frame(pseudo.segments.in)
pseudo.segments.in <- pseudo.segments.in[which(rowSums(pseudo.segments.in) > 0),]

print("Finished reading needed files, extracting cells based on cluster id")

## Subset cells based on cluster id ===============================================================
segments_keep <- which(segments.in %in% SEGMENTs_KEEP)
pseudo.segments_out <- pseudo.segments.in[segments_keep,]
pseudo.segments_out <- as.data.frame(pseudo.segments_out)

print("Finished downsampling, writing output file")

## Write output file ==============================================================================
fwrite(pseudo.segments_out,file = PSEUDO.SEGMENTS.OUT.FILENAME,col.names = FALSE)

print("Completed ExtractClusters.R")
