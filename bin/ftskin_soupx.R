##Run SoupX to remove ambient RNA

##install packages 
install.packages("SoupX")
install.packages("Matric")
install.packages("DoubletFinder")
##load packages
library(seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)

######## why use "=" here and then later change to "<-" ########
dir= '/lustre/scratch127/cellgen/cellgeni/tickets/tic-2250/data/fetal_skin'
dir2= '/output/Gene/filtered/'
dir2_raw='/output/Gene/raw/'
sample_id= read.csv("input/soup_samples.csv")

sample_list = list(sample_id)
###################################################

for (sample in sample_list){
  data_dir <- file.path(dir, sample, dir2)
  data_dir_raw <-file.path(dir, sample, dir2_raw)
  filt.matrix <- Read10X(data_dir )
  raw.matrix  <- Read10X(data_dir_raw)
  srat  <- CreateSeuratObject(counts = filt.matrix)
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  soup.channel  <- autoEstCont(soup.channel)
  out_dir= '/nfs/team298/ew17/fetal_skin/Reanalysis/data/SoupX/'
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  DropletUtils:::write10xCounts(paste0(out_dir, sample, 'soupX_filt'), adj.matrix)
}