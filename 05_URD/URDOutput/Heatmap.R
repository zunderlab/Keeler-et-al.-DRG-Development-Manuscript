#Austin Keeler, University of Virginia
#30 Nov, 2021
#Heatmap for URD marker expression

print("Start extracting marker expression by segment")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

#Try 2
library(ComplexHeatmap)
library(dendsort)
library(viridis)
library(circlize)
library(data.table)

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
HEATMAP.FILENAME <- "/heatmap.csv"
MAIN.COLORMAP <- "viridis"
RIGHT.DEND.WIDTH <- 50
TOP.DEND.HEIGHT <- 50
MAIN.FONTSIZE <- 20

OUTPUT.TYPE <- "pdf" # choose png or pdf
PNG.WIDTH <- 1000
PNG.HEIGHT <- 1000
PDF.WIDTH <- 20
PDF.HEIGHT <- 20

## Read needed files ==============================================================================
heatmap.in <- fread(paste0(INPUT.FOLDER,HEATMAP.FILENAME),header = TRUE, stringsAsFactors = F)
heatmap.in <- as.data.frame(heatmap.in) #maybe change the format
heatmap_in <- as.data.frame(heatmap.in[,2:36])

heatmap.in.seg <- as.data.frame(heatmap.in[,1])
segment.order <- as.character(heatmap.in[,1])

heatmap_df <- t(heatmap_in)
colnames(heatmap_df) <- c(segment.order)

#scaling
heatmap_scale <- as.data.frame(scale(heatmap.in[,2:36]))

#----------------Prepare Dendrograms
row_dend = dendsort(hclust(dist(t(heatmap_scale))))
row_dend = rev(row_dend) #flip to put high expressers on top
col_dend = dendsort(hclust(dist(heatmap_scale)))

#color
col_fun = colorRamp2(c(-3.4, 0, 6), c("blue", "white", "red"))

#----------------Output heatmap image
hm <- Heatmap(heatmap_df, col=col_fun,
              cluster_columns = col_dend,
              row_dend_width=unit(RIGHT.DEND.WIDTH, "mm"),
              column_dend_height=unit(TOP.DEND.HEIGHT, "mm"),
              row_names_gp=gpar(fontsize = MAIN.FONTSIZE),
              show_column_names = TRUE, column_names_side = "bottom",
              column_names_gp = gpar(fontsize = MAIN.FONTSIZE))
# cluster_rows = row_dend, #line removed to make function
# These steps make it take forever to run:
#clustering_distance_rows = robust_dist,
#clustering_distance_columns = robust_dist,

if (OUTPUT.TYPE == "png") {
  png("heatmap.png",width=PNG.WIDTH,height=PNG.HEIGHT)
}
if (OUTPUT.TYPE == "pdf") {
  pdf("heatmap.pdf",width=PDF.WIDTH,height=PDF.HEIGHT)
}

hm
dev.off()

