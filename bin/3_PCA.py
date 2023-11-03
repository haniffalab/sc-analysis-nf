#!/Users/nlg143/anaconda3/envs/scanpy_38/bin/python

#This includes PCA, "Computing the neighborhood graph", "Embedding the neighborhood graph" and 
#"Clustering the neighborhood graph" sections of the Scanpy PBMC3 tutorial 

import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("adata_2pp_out.h5ad")

# Reduce the dimensionality of the data by running principal component analysis (PCA),
# which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata, svd_solver='arpack')

#We can make a scatter plot in the PCA coordinates, but we will not use that later on.
sc.pl.pca(adata, color='CST3', save= "_pcaplot.png", show=False)

#Let us inspect the contribution of single PCs to the total variance in the data. 
# This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, 
# e.g. used in the clustering function sc.tl.louvain() or tSNE sc.tl.tsne().
# In our experience, often a rough estimate of the number of PCs does fine.
sc.pl.pca_variance_ratio(adata, log=True, save= "_plot.png", show=False)

# Save the result.
adata.write(results_file)
###this will save into the directorty called write. 
adata 

##################### Computing the neighborhood graph #####################
# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix.
# You might simply use default values here. 
# For the sake of reproducing Seurat’s results, let’s take the following values.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

##################### Embedding the neighborhood graph #####################
#We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below.
#  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e.,
#  it better preserves trajectories. 
# In some ocassions, you might still observe disconnected clusters and similar connectivity violations. 
# They can usually be remedied by running:
#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#sc.tl.umap(adata, init_pos='paga')
sc.tl.umap(adata)

sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], save="_raw_embedding_neighbourhood.png", show=False)
#As we set the .raw attribute of adata, the previous plots showed the “raw” (normalized, logarithmized, but uncorrected) gene expression.
# You can also plot the scaled and corrected gene expression by explicitly stating that you don’t want to use .raw.

sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False, save="_embedding_neighbourhood.png", show=False)

##################### Clustering the neighborhood graph #####################
#As with Seurat & other frameworks, we recommend  Leiden graph-clustering method (community detection based on optimizing modularity) by Traag et al. (2018).
### Note!! that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
sc.tl.leiden(adata)
#Plot the clusters, which agree quite well with the result of Seurat:
#### UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'], save="_leiden_graph_clustering.png", show=False)

#Save the result.
#adata.write(results_file)
## make this adata file into output:
# To allow for input of the processed adata object to be utilised into further scripts.
#adata.write_h5ad("adata_2pp_out.h5ad") 

###################################################################################################################
# come back to this!! I need to save it to the write folder really but I haven't actually made that yet..         #
###################################################################################################################

