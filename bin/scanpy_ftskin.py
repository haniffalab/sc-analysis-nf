#!/usr/bin/env python3

#0.4.1 QC for mitochondrial contamination, gene counts, cell number and Co
#Merged object fetal skin used
#kernel scVelo

import scanpy as sc
import os
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import math
from scipy.stats import median_abs_deviation
##UserWarning: h5py is running against HDF5 1.14.3 when it was built against 1.14.2, this may cause problems

adata= sc.read('/sc-analysis-nf/output/FCAImmP7241240/Fskin_obj_0_3_1.h5ad')
adata.raw = adata
sample_list = ['FCAImmP7241240']

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
p1.savefig("qcmetrics_totalcounts.png")
# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt", save="_qcmetrics_pctcountsmt.png", show=False )
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_qcmetrics.png", show=False)

# Automatic threshholding based on median absolute deviation (MAD)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

samples_QC = {}
QC = {}
for sample in sample_list:
    adata_byid = adata[adata.obs.sanger_id == sample]
    adata_byid.obs["outlier"] = (
        is_outlier(adata_byid, "log1p_total_counts", 8)
        | is_outlier(adata_byid, "log1p_n_genes_by_counts", 8)
        )
    adata_byid.obs.outlier.value_counts()
    samples_QC[sample] = adata_byid
    QC[sample] = adata_byid
###### RuntimeWarning: Mean of empty slice.
adata = ad.concat(samples_QC, index_unique= None)
adata.obs.outlier.value_counts()

samples_QC = {}
Q = {}
for sample in sample_list:
    adata_byid = adata[adata.obs.sanger_id == sample]
    adata_byid.obs["mt_outlier"] = is_outlier(adata_byid, "pct_counts_mt", 5) | (
        adata_byid.obs["pct_counts_mt"] > 30
        )
    adata_byid.obs.mt_outlier.value_counts()
    samples_QC[sample] = adata_byid
    QC[sample] = adata_byid
#ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
adata = ad.concat(samples_QC, index_unique= None)
adata.obs.mt_outlier.value_counts()

###### would be good to print this in the html page ######
print(f"Total number of cells: {adata.n_obs}")
adataQC = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adataQC.n_obs}")

#### plots4/5/6
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True,
        groupby="sanger_id", save="_all_counts.png", show=False)
sc.pl.violin(adataQC, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True,
        groupby="sanger_id", save="_qc_all_counts.png", show=False)
sc.pl.highest_expr_genes(adata, n_top=20, save="_FCAImmP7241240.png", show=False )

conditions = [
    (adata.obs['predicted_doublets'] == True),
    (adata.obs['outlier'] == True),
        (adata.obs['mt_outlier'] == True),
    (adata.obs['n_genes_by_counts'] < 500),
   (adata.obs['n_genes_by_counts'] >= 500) & (adata.obs['predicted_doublets'] != True) & (adata.obs['mt_outlier'] != True) & (adata.obs['outlier'] != True) 
]

values = ['Doublet', 'Count_outlier','MT_outlier','Low_nFeature', 'Pass']
adata.obs['QC'] = np.select(conditions, values)
adata.obs['QC'] = adata.obs['QC'].astype('category')
adata.obs[['predicted_doublets','outlier', 'mt_outlier','n_genes_by_counts','QC']].head(10)

#adata.obs['QC'].value_counts()
#adata
adata.layers["counts"] = adata.X.copy()
adata.layers

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='sanger_id')
## plot 7 - filter_genes_dispersion_scatter_FCAImmP7241240.png
sc.pl.highly_variable_genes(adata, save="_scatter_FCAImmP7241240.png", show=False)

adata.layers["raw"] = adata.X.copy()
#adata

#adata = adata[:, adata.var.highly_variable]
#adata
unique_donors = list(adata.obs['donor'].unique())
print(unique_donors)

for gender in ['male', 'female']:
    gender_data = adata[adata.obs['gender'] == gender]
    unique_donors = gender_data.obs['donor'].nunique()
    unique_pcw = gender_data.obs['pcw'].unique()

    print(f"Gender: {gender}")
    print(f"Number of Unique Donors: {unique_donors}\n")
    print(f"Unique pcw values: {unique_pcw}\n")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
##plot 8
sc.pl.pca_scatter(adata, color='KRT4', save="_FCAImmP7241240.png", show=False)
##plot 9 
sc.pl.pca_variance_ratio(adata, log=True, save="_FCAImmP7241240.png", show=False)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
##plot 10
#####################################
#colours need to be customisable???
#####################################
sc.pl.umap(adata, color=['LYZ','MS4A1','CD3D'], save="_neighbour_FCAImmP7241240.png", show=False)

sc.tl.leiden(adata)
##needed to conda install -c conda-forge leidenalg

##############################
# turn these into csv files?? 
##############################
s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1',
           'UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7',
           'POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
           'USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8']
g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B',
             'MKI67','TMPO','CENPF','TACC3','PIMREG','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1',
             'KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','JPT1','CDC20','TTK','CDC25C',
             'KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
             'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
cell_cycle_genes = s_genes + g2m_genes

adata= adata.raw.to_adata()
adata
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
len(cell_cycle_genes)

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
#plot11
sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_doubletscores.png", show=False)
adata
adata = adata[adata.obs['QC'] == 'Pass']
adata

adata.layers["counts"] = adata.X.copy()
adata

#plot12 
sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_qc_doubletscores.png", show=False)
#plot13
sc.pl.umap(adata, color=['leiden'], save="_leiden.png", show=False) #UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
#plot14
sc.pl.umap(adata, color=['pct_counts_ribo','pct_counts_mt'], save="_pct_counts.png", show=False)
#plot15
with plt.rc_context({'figure.figsize': (5.5, 4)}): 
    sc.pl.violin(adata, ['pct_counts_ribo', 'pct_counts_mt'], groupby='leiden', stripplot = False, save="_leiden_pct_counts.png", show=False) #userwarning
#plot16 
sc.pl.umap(adata, color=['S_score','G2M_score','leiden'], save="_cellcycle.png", show=False)
#plot17
with plt.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['S_score','G2M_score'], groupby = 'leiden', stripplot = False, save="_cellcycle.png", show=False)
#plot18
with plt.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['outlier','mt_outlier'], groupby = 'leiden', stripplot = False, save="_outlier.png", show=False)
adata
#plot19
sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_plot19.png", show=False)
#plot20
sc.pl.umap(adata, color=['pct_counts_ribo','pct_counts_mt','leiden'], save="_plot20_leiden.png", show=False)
#plot21
with plt.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['pct_counts_ribo', 'pct_counts_mt'], groupby='leiden', stripplot = False, save="_pctcounts_leiden.png", show=False )
#plot22
with plt.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['outlier','mt_outlier'], groupby = 'leiden', stripplot = False, save="_outlier_leiden.png", show=False)
#plot23
sc.pl.umap(adata, color=['S_score','G2M_score','leiden'], save="_cellcycle_leiden.png", show=False)
adata.raw = adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, batch_key='sanger_id')
adata.var.highly_variable.value_counts()
adata

##rather than saving and re-reeading in the same file! made a copy
adata_normalised = adata.copy()
adata_normalised

sc.pp.scale(adata_normalised, max_value=10)

#plot24
sc.tl.pca(adata_normalised, svd_solver='arpack')
sc.pl.pca(adata_normalised, color='CST3', save="_plot24_normalised.png", show=False)
#plot25
sc.pl.umap(adata_normalised, color=['pcw', 'sanger_id', 'gender', 'sorting', 'chemistry_sorting'], save="_metadata.png", show=False)

adata.layers

###Why are they doing this??
#adata.write('/nfs/team298//ew17/fetal_skin/Reanalysis/dataOutputs/Objects/Raw_Objects_v3/Fskin_obj_0_4_1.h5ad')
#adata= sc.read('/nfs/team298//ew17/fetal_skin/Reanalysis/dataOutputs/Objects/Raw_Objects_v3/Fskin_obj_0_4_1.h5ad')


# Remove sample X due to bad quality
#plot26
sc.pl.violin(adata_normalised, ['n_genes_by_counts'],rotation=90, groupby="sanger_id", save="_n_genes_by_counts_sangerid.png", show=False)
#plot27
sc.pl.violin(adata_normalised, ['total_counts'],rotation=90, groupby="sanger_id", save="_total_counts_sangerid.png", show=False)
#plot28
sc.pl.violin(adata_normalised, [ 'pct_counts_mt'],rotation=90,groupby="sanger_id", save="_pct_counts_mt_sangerid.png", show=False)


sc.tl.umap(adata_normalised, min_dist=0.3, spread=1.0,init_pos='spectral') # Using the default distance scaling from umap-learn
adata_normalised.obsm['X_umap_disp_03'] = adata_normalised.obsm['X_umap']

adata.write('Fskin_obj_0_4_1.h5ad')

