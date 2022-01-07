#Kristen Fread, University of Virginia
#20 Oct 2021
#code to save file numbers from concat_transformed separately so that we don't have to load full concat transformed everytime

rm(list = ls(all = TRUE))
#.libPaths("/project/zunderlab/R/3.6.3")

library(ZunderPipelineFunctions)
library(data.table)

## Input parameters
INPUT.FOLDER <- getwd()
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
Fil_name <- "filteredmatrix_file_numbers.csv"


## Read needed files
concat.transformed <- read.concat.transformed(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME)


#save file numbers
fil_nums_filtered <- as.data.frame(concat.transformed[,"File"])
colnames(fil_nums_filtered) <- "File" #renaming file column
fwrite(fil_nums_filtered,file = paste(INPUT.FOLDER,Fil_name,sep = '/'),row.names = FALSE)