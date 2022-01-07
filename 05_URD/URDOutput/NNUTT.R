#Originally created by Emily Puleo. Modified on 10/19/21 by Austin Keeler


rm(list = ls())

# Load packages
#library(tidyverse)
library(ZunderPipelineFunctions)
library(ggfortify)
library(data.table)
library(URD)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define constants
INPUT_FOLDER <- getwd()
WRITE_PATH<-getwd()
# these are the cells you want to map
RA_FILENAME1 <- "Concat_Transformed_1.csv" #ie, dataset1
#check segmentlist on line 88 and change as needed

#this is the URD you are using
IV_FILENAME <- "Concat_Transformed.csv" #ie, dataset2
URD_IN <- "/URD.RData"
#^open the URD workspace into the environment manually

# define distance function
euc_dist <- function(row1, row2) {
  row1 <- as.numeric(row1)
  row2 <- as.numeric(row2)
  # make sure lengths are equal
  if (length(row1) == length(row2)) { 
    num_params <- length(row1)
  } else {
    stop("Length of row1 does not equal length of row2.")
  }
  ind_sum <- vector(length = num_params)
  for (param_i in seq_len(num_params)) { 
    #iterate over parameters, store each component of dist in ind_sum
    ind_sum[param_i] <- (row1[param_i] - row2[param_i])^2
  }
  dist <- sqrt(sum(ind_sum))
  return(dist)
}
# Read in files =====================================================================================
load(paste0(INPUT_FOLDER,URD_IN))
dataset_RA <- fread(file.path(INPUT_FOLDER, RA_FILENAME1), 
                     stringsAsFactors = F,
                     data.table = F)

iv_dataset <- fread(file.path(INPUT_FOLDER, IV_FILENAME),
                    stringsAsFactors = F,
                    data.table = F)

##FILE 1 ======================================
#isolate file
expr_RA <- dataset_RA

# combine urd expr data w/ comp tSNE layout coords
iv_dataset <- cbind(iv_dataset, URD_Object@tsne.y$tSNE1,URD_Object@tsne.y$tSNE2)
iv_expr <- select(iv_dataset, -c("File","URD_Object@tsne.y$tSNE1","URD_Object@tsne.y$tSNE2"))

#make sure panel is the same
expr_RA<-expr_RA[,colnames(expr_RA)==colnames(iv_expr)]
iv_expr<-iv_expr[,colnames(expr_RA)==colnames(iv_expr)]
colnames(expr_RA)==colnames(iv_expr)
#^ check to make sure all TRUE

#make dists vector
dists <- vector(mode = "list", length = nrow(expr_RA))

#calc dist
for (i in 1:nrow(expr_RA)) {
  cell <- expr_RA[i,]
  dists[[i]] <- apply(iv_expr, 1, function(comp_cell) {
    euc_dist(comp_cell, cell)
  }
  )
}


tsne_ind<-rep(0,length(dists))
for (i in 1:length(dists)){
  tsne_ind[i]<-which(dists[[i]]==min(dists[[i]]))
}

segmentlist<-c() #list all segments in your data
new_index<-matrix(data=NA, nrow=length(tsne_ind),ncol=1)
for (k in 1:length(segmentlist)){
  for (i in 1:length(tsne_ind)){
    for (j in 1:length(URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]])){
      if (paste0("V",tsne_ind[i]) == URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]][j]){
        new_index[i]<-tsne_ind[i]
        #which(paste0("V",tsne_ind[i]) == URD_Object.tree@tree$cells.in.segment[[segmentlist[k]]][j])]
      }
    }}}

new_index<-data.frame(new_index)
new_index %>% drop_na()
new_index<-as.list(new_index)
#new_index<-as.matrix(new_index)
new_index<-as.numeric(unlist(new_index))

#plotURD
similar.cells1<-URD_Object@group.ids$init[new_index] #most similar points
#URD_Object.tree.new<-urdSubset(URD_Object, cells.keep=whichCells(URD_Object,"init",similar.cells))
png("NNUTT.png")
plotTree(URD_Object.tree,"Stage",cells.highlight=whichCells(URD_Object,"init",similar.cells1), cell.alpha=0, cells.highlight.alpha = 1, cells.highlight.size=2)
dev.off()
