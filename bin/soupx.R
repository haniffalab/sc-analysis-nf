#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##Run SoupX to remove ambient RNA
#script 1 of the pipeline

##install packages 
#install.packages("SoupX")
#install.packages("matric")
#install.packages("ggplot2")
#install.packages("knitr")
#install.packages("Seurat")
#install.packages("remotes")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
########## weird install #########
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DropletUtils")

#list.of.packages <- c("Dropletutils")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) BiocManager::install("DropletUtils")
##################################
##load packages
library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)
library(sp)
library(SeuratObject)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input directory).n", call.=FALSE)
} 

data_dir <- file.path(args[1], 'output/Gene/filtered')
data_dir_raw <-file.path(args[1], 'output/Gene/raw')
filt.matrix <- Read10X(data_dir)
raw.matrix <- Read10X(data_dir_raw)
srat <- CreateSeuratObject(counts = filt.matrix)
srat <- SCTransform(srat, verbose = F)
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = T)
soup.channel <- SoupChannel(raw.matrix, filt.matrix)
meta <- srat@meta.data
umap <- srat@reductions$umap@cell.embeddings
soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel <- setDR(soup.channel, umap)
soup.channel <- autoEstCont(soup.channel)
########### save plot #############
png("soup.channel.png", width = 6.67, height = 6.67, units = "in", res = 300)
print(autoEstCont(soup.channel))
dev.off()
###################################
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
DropletUtils:::write10xCounts("adjusted_counts", adj.matrix)



