# Created by Austin Keeler & Ashley Hirt
# 15 July 2020
# Plot a single UMAP colored by a metadata variable

library(tidyverse)
library(scales)
library(ZunderPipelineFunctions)
library(data.table)
library(viridis)


# ONLY WORKS IF SAME # CELLS IN UMAP.LAYOUT AS CONCAT.TRANSFORMED.
# ALSO TO ENSURE YOU GET COLD TO HOT COLORING THE WAY YOU WANT, MAKE SURE 
# THE ROWS ARE ARRANGED SO COLDER LEVELS AT THE TOP (EG EARLY TIMEPOINTS) 
# AND HOTTER LEVELS ARE AT THE BOTTOM. THE ROW ORDERING WILL DETERMINE
# THE COLOR PROGRESSION.



# Input parameters ----
# make sure to set your working directory!!
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata.csv"
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
POINT.SIZE <- 0.0001
VAR <- "age" # make sure this matches your md column name exactly
OUTPUT.FILENAME <- paste0("Colored_by_", VAR, ".png") 
# ^^ this is just a default, you can change--supply a full string, eg "My_File_Name.png")


# Read files ----
# read metadata
md <- fread(paste0(INPUT.FOLDER, "/", METADATA.FILENAME), stringsAsFactors = F)
# assign file numbering
md <- cbind(md, "File" = 1:nrow(md)) 

# read concat.transformed (we're just using the File col)
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCATTRANSFORMED.FILENAME), 
                            stringsAsFactors = F)

# Read UMAP layout, has two columns, x and y, where each row is a cell
umap_layout <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
# this is terrible; combine concat File column w/ umap df
# this assumes same # of cells and same "ordering" of cells (first at the top)
plotting_df <- cbind(umap_layout, concat.transformed$File)
# Assign column names
colnames(plotting_df) <- c("umap_x", "umap_y", "File")


# Join and Plot ----

# join on File 
plotting.df <- merge(plotting_df, md)
set.seed(42)
random.plotting.df <- plotting.df[sample(nrow(plotting.df)),]

# save plot
ggsave(OUTPUT.FILENAME, 
       plot = {
               ggplot(data=random.plotting.df, aes_string(x="umap_x", y="umap_y", color=VAR)) + 
                       # note color=VAR, that determines coloring
                       geom_point(size = POINT.SIZE) + 
                       scale_color_viridis_d(option = "plasma") + 
                       theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.background = element_blank(), 
                             axis.line = element_line(colour = "black")
                       )
       },
       device = "png", height = 14,width = 21, units = "in", dpi = 600 
)

