rm(list = ls(all = TRUE))
.libPaths( c( .libPaths(), "/projects/zunderlab/R/3.6") )

INPUT.DIR <- "" #directory where FCS files are
print(INPUT.DIR)
source("BatchAdjust.R")

BatchAdjust(basedir=".", outdir=".", channelsFile = "ChannelsToAdjust.txt", batchKeyword="Set", anchorKeyword="UNIV",
            method="65p",transformation=FALSE,addExt="_65p",plotDiagnostics=TRUE)
