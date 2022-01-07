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
#IMPORTANT - ADD Number_X rows to match the number of clusters in your parent data (before any was extracted).
#My data had 20 clusters, so there are 20 rows of NUMBER_1 through NUMBER_20.  Also I duplicated 20 downsampling
#sets below, one for each cluster. Modify how many there are to match the number of clusters.

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
DOWNSAMPLE.OUT.FILENAME <- "Concat_Transformed_.csv" #new file
SUBSAMPLE.ID.OUT.FILENAME <- "Subsample_ID_.csv" #new file
LAYOUT.OUT.FILENAME <- "UMAP_layout_.csv" #new file
CLUSTERS.OUT.FILENAME <- "clusters_.csv" #new file
NUMBER_1 <- 0
NUMBER_2 <- 0
NUMBER_3 <- 0
NUMBER_4 <- 0
NUMBER_5 <- 0
NUMBER_6 <- 0
NUMBER_7 <- 0
NUMBER_8 <- 0
NUMBER_9 <- 0
NUMBER_10 <- 0
NUMBER_11 <- 0
NUMBER_12 <- 0
NUMBER_13 <- 0
NUMBER_14 <- 0
NUMBER_15 <- 0
NUMBER_16 <- 0
NUMBER_17 <- 0
NUMBER_18 <- 0
NUMBER_19 <- 0
NUMBER_20 <- 0
NUMBER_21 <- 0
NUMBER_22 <- 0
NUMBER_23 <- 0
NUMBER_24 <- 0
NUMBER_25 <- 0
NUMBER_26 <- 0
NUMBER_27 <- 0
NUMBER_28 <- 0
NUMBER_29 <- 0
NUMBER_30 <- 0
NUMBER_31 <- 0
NUMBER_32 <- 0
NUMBER_33 <- 0
NUMBER_34 <- 0
NUMBER_35 <- 0
NUMBER_36 <- 0
NUMBER_37 <- 0
NUMBER_38 <- 0
NUMBER_39 <- 0
NUMBER_40 <- 0
NUMBER_41 <- 0
NUMBER_42 <- 0

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


#CLUSTER 1
    #Get vector of all ids for this cluster
all.ids <- which(clusters.in == 1)
    #randomly subsample from the selection of cells
subsample.ids[[1]] <- sample(all.ids,size = NUMBER_1)
    #use subsample ids to subsample from expression info
expression.out[[1]] <- concat.transformed[subsample.ids[[1]],]
    #store subsampled layout
layout.out[[1]] <- layout.in[subsample.ids[[1]],]
    #store cluster id's
clusters.out[[1]] <- clusters.in[subsample.ids[[1]]]

#CLUSTER 2
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 2)
#randomly subsample from the selection of cells
subsample.ids[[2]] <- sample(all.ids,size = NUMBER_2)
#use subsample ids to subsample from expression info
expression.out[[2]] <- concat.transformed[subsample.ids[[2]],]
#store subsampled layout
layout.out[[2]] <- layout.in[subsample.ids[[2]],]
#store cluster id's
clusters.out[[2]] <- clusters.in[subsample.ids[[2]]]

#CLUSTER 3
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 3)
#randomly subsample from the selection of cells
subsample.ids[[3]] <- sample(all.ids,size = NUMBER_3)
#use subsample ids to subsample from expression info
expression.out[[3]] <- concat.transformed[subsample.ids[[3]],]
#store subsampled layout
layout.out[[3]] <- layout.in[subsample.ids[[3]],]
#store cluster id's
clusters.out[[3]] <- clusters.in[subsample.ids[[3]]]

#CLUSTER 4
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 4)
#randomly subsample from the selection of cells
subsample.ids[[4]] <- sample(all.ids,size = NUMBER_4)
#use subsample ids to subsample from expression info
expression.out[[4]] <- concat.transformed[subsample.ids[[4]],]
#store subsampled layout
layout.out[[4]] <- layout.in[subsample.ids[[4]],]
#store cluster id's
clusters.out[[4]] <- clusters.in[subsample.ids[[4]]]

#CLUSTER 5
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 5)
#randomly subsample from the selection of cells
subsample.ids[[5]] <- sample(all.ids,size = NUMBER_5)
#use subsample ids to subsample from expression info
expression.out[[5]] <- concat.transformed[subsample.ids[[5]],]
#store subsampled layout
layout.out[[5]] <- layout.in[subsample.ids[[5]],]
#store cluster id's
clusters.out[[5]] <- clusters.in[subsample.ids[[5]]]

#CLUSTER 6
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 6)
#randomly subsample from the selection of cells
subsample.ids[[6]] <- sample(all.ids,size = NUMBER_6)
#use subsample ids to subsample from expression info
expression.out[[6]] <- concat.transformed[subsample.ids[[6]],]
#store subsampled layout
layout.out[[6]] <- layout.in[subsample.ids[[6]],]
#store cluster id's
clusters.out[[6]] <- clusters.in[subsample.ids[[6]]]

#CLUSTER 7
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 7)
#randomly subsample from the selection of cells
subsample.ids[[7]] <- sample(all.ids,size = NUMBER_7)
#use subsample ids to subsample from expression info
expression.out[[7]] <- concat.transformed[subsample.ids[[7]],]
#store subsampled layout
layout.out[[7]] <- layout.in[subsample.ids[[7]],]
#store cluster id's
clusters.out[[7]] <- clusters.in[subsample.ids[[7]]]

#CLUSTER 8
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 8)
#randomly subsample from the selection of cells
subsample.ids[[8]] <- sample(all.ids,size = NUMBER_8)
#use subsample ids to subsample from expression info
expression.out[[8]] <- concat.transformed[subsample.ids[[8]],]
#store subsampled layout
layout.out[[8]] <- layout.in[subsample.ids[[8]],]
#store cluster id's
clusters.out[[8]] <- clusters.in[subsample.ids[[8]]]

#CLUSTER 9
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 9)
#randomly subsample from the selection of cells
subsample.ids[[9]] <- sample(all.ids,size = NUMBER_9)
#use subsample ids to subsample from expression info
expression.out[[9]] <- concat.transformed[subsample.ids[[9]],]
#store subsampled layout
layout.out[[9]] <- layout.in[subsample.ids[[9]],]
#store cluster id's
clusters.out[[9]] <- clusters.in[subsample.ids[[9]]]

#CLUSTER 10
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 10)
#randomly subsample from the selection of cells
subsample.ids[[10]] <- sample(all.ids,size = NUMBER_10)
#use subsample ids to subsample from expression info
expression.out[[10]] <- concat.transformed[subsample.ids[[10]],]
#store subsampled layout
layout.out[[10]] <- layout.in[subsample.ids[[10]],]
#store cluster id's
clusters.out[[10]] <- clusters.in[subsample.ids[[10]]]

#CLUSTER 11
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 11)
#randomly subsample from the selection of cells
subsample.ids[[11]] <- sample(all.ids,size = NUMBER_11)
#use subsample ids to subsample from expression info
expression.out[[11]] <- concat.transformed[subsample.ids[[11]],]
#store subsampled layout
layout.out[[11]] <- layout.in[subsample.ids[[11]],]
#store cluster id's
clusters.out[[11]] <- clusters.in[subsample.ids[[11]]]

#CLUSTER 12
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 12)
#randomly subsample from the selection of cells
subsample.ids[[12]] <- sample(all.ids,size = NUMBER_12)
#use subsample ids to subsample from expression info
expression.out[[12]] <- concat.transformed[subsample.ids[[12]],]
#store subsampled layout
layout.out[[12]] <- layout.in[subsample.ids[[12]],]
#store cluster id's
clusters.out[[12]] <- clusters.in[subsample.ids[[12]]]

#CLUSTER 13
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 13)
#randomly subsample from the selection of cells
subsample.ids[[13]] <- sample(all.ids,size = NUMBER_13)
#use subsample ids to subsample from expression info
expression.out[[13]] <- concat.transformed[subsample.ids[[13]],]
#store subsampled layout
layout.out[[13]] <- layout.in[subsample.ids[[13]],]
#store cluster id's
clusters.out[[13]] <- clusters.in[subsample.ids[[13]]]

#CLUSTER 14
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 14)
#randomly subsample from the selection of cells
subsample.ids[[14]] <- sample(all.ids,size = NUMBER_14)
#use subsample ids to subsample from expression info
expression.out[[14]] <- concat.transformed[subsample.ids[[14]],]
#store subsampled layout
layout.out[[14]] <- layout.in[subsample.ids[[14]],]
#store cluster id's
clusters.out[[14]] <- clusters.in[subsample.ids[[14]]]
    
#CLUSTER 15
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 15)
#randomly subsample from the selection of cells
subsample.ids[[15]] <- sample(all.ids,size = NUMBER_15)
#use subsample ids to subsample from expression info
expression.out[[15]] <- concat.transformed[subsample.ids[[15]],]
#store subsampled layout
layout.out[[15]] <- layout.in[subsample.ids[[15]],]
#store cluster id's
clusters.out[[15]] <- clusters.in[subsample.ids[[15]]]

#CLUSTER 16
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 16)
#randomly subsample from the selection of cells
subsample.ids[[16]] <- sample(all.ids,size = NUMBER_16)
#use subsample ids to subsample from expression info
expression.out[[16]] <- concat.transformed[subsample.ids[[16]],]
#store subsampled layout
layout.out[[16]] <- layout.in[subsample.ids[[16]],]
#store cluster id's
clusters.out[[16]] <- clusters.in[subsample.ids[[16]]]

#CLUSTER 17
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 17)
#randomly subsample from the selection of cells
subsample.ids[[17]] <- sample(all.ids,size = NUMBER_17)
#use subsample ids to subsample from expression info
expression.out[[17]] <- concat.transformed[subsample.ids[[17]],]
#store subsampled layout
layout.out[[17]] <- layout.in[subsample.ids[[17]],]
#store cluster id's
clusters.out[[17]] <- clusters.in[subsample.ids[[17]]]

#CLUSTER 18
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 18)
#randomly subsample from the selection of cells
subsample.ids[[18]] <- sample(all.ids,size = NUMBER_18)
#use subsample ids to subsample from expression info
expression.out[[18]] <- concat.transformed[subsample.ids[[18]],]
#store subsampled layout
layout.out[[18]] <- layout.in[subsample.ids[[18]],]
#store cluster id's
clusters.out[[18]] <- clusters.in[subsample.ids[[18]]]

#CLUSTER 19
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 19)
#randomly subsample from the selection of cells
subsample.ids[[19]] <- sample(all.ids,size = NUMBER_19)
#use subsample ids to subsample from expression info
expression.out[[19]] <- concat.transformed[subsample.ids[[19]],]
#store subsampled layout
layout.out[[19]] <- layout.in[subsample.ids[[19]],]
#store cluster id's
clusters.out[[19]] <- clusters.in[subsample.ids[[19]]]

#CLUSTER 20
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 20)
#randomly subsample from the selection of cells
subsample.ids[[20]] <- sample(all.ids,size = NUMBER_20)
#use subsample ids to subsample from expression info
expression.out[[20]] <- concat.transformed[subsample.ids[[20]],]
#store subsampled layout
layout.out[[20]] <- layout.in[subsample.ids[[20]],]
#store cluster id's
clusters.out[[20]] <- clusters.in[subsample.ids[[20]]]

#CLUSTER 21
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 21)
#randomly subsample from the selection of cells
subsample.ids[[21]] <- sample(all.ids,size = NUMBER_21)
#use subsample ids to subsample from expression info
expression.out[[21]] <- concat.transformed[subsample.ids[[21]],]
#store subsampled layout
layout.out[[21]] <- layout.in[subsample.ids[[21]],]
#store cluster id's
clusters.out[[21]] <- clusters.in[subsample.ids[[21]]]

#CLUSTER 22
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 22)
#randomly subsample from the selection of cells
subsample.ids[[22]] <- sample(all.ids,size = NUMBER_22)
#use subsample ids to subsample from expression info
expression.out[[22]] <- concat.transformed[subsample.ids[[22]],]
#store subsampled layout
layout.out[[22]] <- layout.in[subsample.ids[[22]],]
#store cluster id's
clusters.out[[22]] <- clusters.in[subsample.ids[[22]]]

#CLUSTER 23
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 23)
#randomly subsample from the selection of cells
subsample.ids[[23]] <- sample(all.ids,size = NUMBER_23)
#use subsample ids to subsample from expression info
expression.out[[23]] <- concat.transformed[subsample.ids[[23]],]
#store subsampled layout
layout.out[[23]] <- layout.in[subsample.ids[[23]],]
#store cluster id's
clusters.out[[23]] <- clusters.in[subsample.ids[[23]]]

#CLUSTER 24
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 24)
#randomly subsample from the selection of cells
subsample.ids[[24]] <- sample(all.ids,size = NUMBER_24)
#use subsample ids to subsample from expression info
expression.out[[24]] <- concat.transformed[subsample.ids[[24]],]
#store subsampled layout
layout.out[[24]] <- layout.in[subsample.ids[[24]],]
#store cluster id's
clusters.out[[24]] <- clusters.in[subsample.ids[[24]]]

#CLUSTER 25
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 25)
#randomly subsample from the selection of cells
subsample.ids[[25]] <- sample(all.ids,size = NUMBER_25)
#use subsample ids to subsample from expression info
expression.out[[25]] <- concat.transformed[subsample.ids[[25]],]
#store subsampled layout
layout.out[[25]] <- layout.in[subsample.ids[[25]],]
#store cluster id's
clusters.out[[25]] <- clusters.in[subsample.ids[[25]]]

#CLUSTER 26
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 26)
#randomly subsample from the selection of cells
subsample.ids[[26]] <- sample(all.ids,size = NUMBER_26)
#use subsample ids to subsample from expression info
expression.out[[26]] <- concat.transformed[subsample.ids[[26]],]
#store subsampled layout
layout.out[[26]] <- layout.in[subsample.ids[[26]],]
#store cluster id's
clusters.out[[26]] <- clusters.in[subsample.ids[[26]]]

#CLUSTER 27
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 27)
#randomly subsample from the selection of cells
subsample.ids[[27]] <- sample(all.ids,size = NUMBER_27)
#use subsample ids to subsample from expression info
expression.out[[27]] <- concat.transformed[subsample.ids[[27]],]
#store subsampled layout
layout.out[[27]] <- layout.in[subsample.ids[[27]],]
#store cluster id's
clusters.out[[27]] <- clusters.in[subsample.ids[[27]]]

#CLUSTER 28
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 28)
#randomly subsample from the selection of cells
subsample.ids[[28]] <- sample(all.ids,size = NUMBER_28)
#use subsample ids to subsample from expression info
expression.out[[28]] <- concat.transformed[subsample.ids[[28]],]
#store subsampled layout
layout.out[[28]] <- layout.in[subsample.ids[[28]],]
#store cluster id's
clusters.out[[28]] <- clusters.in[subsample.ids[[28]]]

#CLUSTER 29
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 29)
#randomly subsample from the selection of cells
subsample.ids[[29]] <- sample(all.ids,size = NUMBER_29)
#use subsample ids to subsample from expression info
expression.out[[29]] <- concat.transformed[subsample.ids[[29]],]
#store subsampled layout
layout.out[[29]] <- layout.in[subsample.ids[[29]],]
#store cluster id's
clusters.out[[29]] <- clusters.in[subsample.ids[[29]]]

#CLUSTER 30
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 30)
#randomly subsample from the selection of cells
subsample.ids[[30]] <- sample(all.ids,size = NUMBER_30)
#use subsample ids to subsample from expression info
expression.out[[30]] <- concat.transformed[subsample.ids[[30]],]
#store subsampled layout
layout.out[[30]] <- layout.in[subsample.ids[[30]],]
#store cluster id's
clusters.out[[30]] <- clusters.in[subsample.ids[[30]]]

#CLUSTER 31
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 31)
#randomly subsample from the selection of cells
subsample.ids[[31]] <- sample(all.ids,size = NUMBER_31)
#use subsample ids to subsample from expression info
expression.out[[31]] <- concat.transformed[subsample.ids[[31]],]
#store subsampled layout
layout.out[[31]] <- layout.in[subsample.ids[[31]],]
#store cluster id's
clusters.out[[31]] <- clusters.in[subsample.ids[[31]]]

#CLUSTER 32
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 32)
#randomly subsample from the selection of cells
subsample.ids[[32]] <- sample(all.ids,size = NUMBER_32)
#use subsample ids to subsample from expression info
expression.out[[32]] <- concat.transformed[subsample.ids[[32]],]
#store subsampled layout
layout.out[[32]] <- layout.in[subsample.ids[[32]],]
#store cluster id's
clusters.out[[32]] <- clusters.in[subsample.ids[[32]]]

#CLUSTER 33
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 33)
#randomly subsample from the selection of cells
subsample.ids[[33]] <- sample(all.ids,size = NUMBER_33)
#use subsample ids to subsample from expression info
expression.out[[33]] <- concat.transformed[subsample.ids[[33]],]
#store subsampled layout
layout.out[[33]] <- layout.in[subsample.ids[[33]],]
#store cluster id's
clusters.out[[33]] <- clusters.in[subsample.ids[[33]]]

#CLUSTER 34
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 34)
#randomly subsample from the selection of cells
subsample.ids[[34]] <- sample(all.ids,size = NUMBER_34)
#use subsample ids to subsample from expression info
expression.out[[34]] <- concat.transformed[subsample.ids[[34]],]
#store subsampled layout
layout.out[[34]] <- layout.in[subsample.ids[[34]],]
#store cluster id's
clusters.out[[34]] <- clusters.in[subsample.ids[[34]]]

#CLUSTER 35
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 35)
#randomly subsample from the selection of cells
subsample.ids[[35]] <- sample(all.ids,size = NUMBER_35)
#use subsample ids to subsample from expression info
expression.out[[35]] <- concat.transformed[subsample.ids[[35]],]
#store subsampled layout
layout.out[[35]] <- layout.in[subsample.ids[[35]],]
#store cluster id's
clusters.out[[35]] <- clusters.in[subsample.ids[[35]]]

#CLUSTER 36
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 36)
#randomly subsample from the selection of cells
subsample.ids[[36]] <- sample(all.ids,size = NUMBER_36)
#use subsample ids to subsample from expression info
expression.out[[36]] <- concat.transformed[subsample.ids[[36]],]
#store subsampled layout
layout.out[[36]] <- layout.in[subsample.ids[[36]],]
#store cluster id's
clusters.out[[36]] <- clusters.in[subsample.ids[[36]]]

#CLUSTER 37
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 37)
#randomly subsample from the selection of cells
subsample.ids[[37]] <- sample(all.ids,size = NUMBER_37)
#use subsample ids to subsample from expression info
expression.out[[37]] <- concat.transformed[subsample.ids[[37]],]
#store subsampled layout
layout.out[[37]] <- layout.in[subsample.ids[[37]],]
#store cluster id's
clusters.out[[37]] <- clusters.in[subsample.ids[[37]]]

#CLUSTER 38
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 38)
#randomly subsample from the selection of cells
subsample.ids[[38]] <- sample(all.ids,size = NUMBER_38)
#use subsample ids to subsample from expression info
expression.out[[38]] <- concat.transformed[subsample.ids[[38]],]
#store subsampled layout
layout.out[[38]] <- layout.in[subsample.ids[[38]],]
#store cluster id's
clusters.out[[38]] <- clusters.in[subsample.ids[[38]]]

#CLUSTER 39
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 39)
#randomly subsample from the selection of cells
subsample.ids[[39]] <- sample(all.ids,size = NUMBER_39)
#use subsample ids to subsample from expression info
expression.out[[39]] <- concat.transformed[subsample.ids[[39]],]
#store subsampled layout
layout.out[[39]] <- layout.in[subsample.ids[[39]],]
#store cluster id's
clusters.out[[39]] <- clusters.in[subsample.ids[[39]]]

#CLUSTER 40
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 40)
#randomly subsample from the selection of cells
subsample.ids[[40]] <- sample(all.ids,size = NUMBER_40)
#use subsample ids to subsample from expression info
expression.out[[40]] <- concat.transformed[subsample.ids[[40]],]
#store subsampled layout
layout.out[[40]] <- layout.in[subsample.ids[[40]],]
#store cluster id's
clusters.out[[40]] <- clusters.in[subsample.ids[[40]]]

#CLUSTER 41
#Get vector of all ids for this cluster
all.ids <- which(clusters.in == 41)
#randomly subsample from the selection of cells
subsample.ids[[41]] <- sample(all.ids,size = NUMBER_41)
#use subsample ids to subsample from expression info
expression.out[[41]] <- concat.transformed[subsample.ids[[41]],]
#store subsampled layout
layout.out[[41]] <- layout.in[subsample.ids[[41]],]
#store cluster id's
clusters.out[[41]] <- clusters.in[subsample.ids[[41]]]


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
