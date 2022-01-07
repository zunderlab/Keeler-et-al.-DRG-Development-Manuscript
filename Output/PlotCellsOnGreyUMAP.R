#Austin Keeler, University of Virginia
#2 Nov, 2021
#Highlight specified cells on grey UMAP

print("Start PlotCellsOnGreyUMAP.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ZunderPipelineFunctions)
library(ggfortify)

print("libraries loaded")

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
GREY.LAYOUT.FILENAME <- "/UMAP_layout.csv" #This should be the larger dataset, to provide a background umap
GREY.POINT.SIZE <- 0.1 #should be smaller point size, creating a backdrop for highlighted cells
GREY.COLOR <- "#bdbdbd"
HIGHLIGHT.LAYOUT.FILENAME <- "/UMAP_layout_highlight.csv" #smaller dataset of cells to be highlighted
HIGHLIGHT.POINT.SIZE <- 1 #should be larger point size so these cells pop, vary depending on density
HIGHLIGHT.COLOR <- "red" #some bright, obvious color
HIGHLIGHT.PLOT.FILENAME <- "Highlighted.png" #change name to reflect what is shown
HEIGHT <- 5
WIDTH <- 5

print("input parameters loaded, reading needed files")

## Read needed files ==============================================================================
grey.layout.in <- read.layout(INPUT.FOLDER,GREY.LAYOUT.FILENAME)
highlight.layout.in <- read.layout(INPUT.FOLDER,HIGHLIGHT.LAYOUT.FILENAME)

print("needed files read, prepping data to plot")

#plot
ggsave(HIGHLIGHT.PLOT.FILENAME,plot = ggplot() + 
         geom_point(data = grey.layout.in,
                    aes(x=umap_x,y=umap_y),color = GREY.COLOR, size = GREY.POINT.SIZE) + 
         geom_point(data = highlight.layout.in,
                    aes(x=umap_x,y=umap_y),color = HIGHLIGHT.COLOR, size = HIGHLIGHT.POINT.SIZE)+ 
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black")) +
         guides(colour = guide_legend(override.aes = list(shape=15, size=8))),
       height = HEIGHT, width = WIDTH, units = "in", dpi = 300)

print("data plotted and file outputted")
