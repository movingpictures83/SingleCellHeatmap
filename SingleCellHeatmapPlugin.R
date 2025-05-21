library(dplyr)
library(Seurat)
library(patchwork)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
	pdf(outputfile)
	pbmc.markers <- readRDS(paste(pfix, parameters["rdsfile", 2], sep="/"))
	pbmc <- readRDS(paste(pfix, parameters["rdsfile2", 2], sep="/"))
topN <- as.integer(parameters["N", 2])

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_logFC > 1) %>%
    slice_head(n = topN) %>%
    ungroup() -> top10

#print(str(pbmc))
write.csv(pbmc@commands$RunPCA.RNA@params$features, paste(outputfile, "csv", sep="."))
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
}


