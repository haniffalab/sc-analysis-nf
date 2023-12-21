#!/usr/bin/env Rscript

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
##################################
##load packages
library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)

######## why use "=" here and then later change to "<-" ########
dir= '/Users/nlg143/projects/sc-analysis-nf/input/'
dir2= '/output/Gene/filtered/'
dir2_raw='/output/Gene/raw/'

sample = 'FCAImmP7241240'
###################################################
data_dir <- file.path(dir, sample, dir2)
data_dir_raw <-file.path(dir, sample, dir2_raw)
filt.matrix <- Read10X(data_dir )
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
out_dir= '/Users/nlg143/projects/sc-analysis-nf/output/'
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
DropletUtils:::write10xCounts('soup', adj.matrix)

