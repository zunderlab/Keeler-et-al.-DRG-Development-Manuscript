#Austin Keeler, adapated from Corey Williams, University of Virginia
#01 August, 2020
#Plot cells over threshold MARKER EXPRESSION
#COMPLETED!!


print("Start PlotBySilhouetteScore_Threshold.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(ggfortify)
library(data.table)
library(dplyr)

print("libraries loaded")

## Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CLUSTERS.FILENAME <- "/clusters_FINAL.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CONCAT.TRANSFORMED <- "/Concat_Transformed_Extracted_Neurons_sans6_16.csv"
CONCAT.OUT.FILENAME <- "Concat_Transformed_TuJ1low.csv"
LAYOUT.OUT.FILENAME <- "UMAP_layout_TuJ1low.csv"
CLUSTERS.OUT.FILENAME <- "clusters_TuJ1low.csv" #will be needed to cooperate with pipeline
EXPRESSION.VALUE <- 2

print("input parameters loaded, reading needed files")

## Read needed files
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)
concat.transformed.in <- fread(paste0(INPUT.FOLDER,CONCAT.TRANSFORMED), stringsAsFactors = F)
colnames(concat.transformed.in) <- gsub("-", "_", colnames(concat.transformed.in))
colnames(concat.transformed.in) <- gsub(".", "_", colnames(concat.transformed.in), fixed = TRUE)
concat.transformed.in <- as.data.frame(concat.transformed.in)

#assign cell ID number to silhouettescore
cell_id <- 1:nrow(layout.in)
cell_id_df <- data.frame("cell_id" = cell_id)
layout.id <- cbind(layout.in,cell_id_df)
layout.clusters.id <- cbind(clusters.in,layout.id)

#assign cell ID number to concat_transform
cell_id <- 1:nrow(layout.in)
cell_id_df <- data.frame("cell_id" = cell_id)
layout.id <- cbind(layout.in,cell_id_df)
layout.concat.transform.id <- cbind(layout.id,concat.transformed.in)

print("needed files read, prepping data to threshold")

## Prep dataframe for thresholding
#thresholding.df <- merge(layout.clusters.id,layout.concat.transform.id)
#colnames(thresholding.df) <- c("cell_id","clusters","umap_x","umap_y","silhouettescore")
thresholding.ct.df <- merge(layout.clusters.id,layout.concat.transform.id)
thresholded.ct.df <- thresholding.ct.df[thresholding.ct.df$TuJ1_Y89Di <= EXPRESSION.VALUE, ]

print("data ready to write")

## Pull thresholded.ct.df apart into layout, cluster, and concat.transform

layout.out <- thresholded.ct.df[,c(1:2)]
clusters.out <- thresholded.ct.df[,c(4)]
concat.transformed.out <- thresholded.ct.df[,c(5:46)]

## Write output file ==============================================================================
fwrite(concat.transformed.out,file = CONCAT.OUT.FILENAME,row.names = FALSE)
fwrite(layout.out,file = LAYOUT.OUT.FILENAME)
fwrite(list(clusters.out),file = CLUSTERS.OUT.FILENAME,col.names = FALSE,sep = ",")
