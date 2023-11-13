#!/usr/bin/env python3

import fire
import scanpy as sc
import pandas as pd
import numpy as np


def marker_genes(mg_adata_input:str, mg_outfile_mid:str, mg_outfile_end:str, mg_outfile_withoutx:str ):
    #read in output file from third script 
    adata = sc.read_h5ad(mg_adata_input)
    results_file = 'pbmc3k.h5ad'  # the file that will store the analysis result
    ###################################### Finding marker genes #############################################
    #Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute 
    # of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_plot1.png",show=False)
    ########################################################################################################
    #                         RuntimeWarning: invalid value encountered in log2                            #
    #           an extra graph in my plot?? also is all of the script inplace=true?                        # 
    ########################################################################################################

    sc.settings.verbosity = 2  # reduce the verbosity
    # The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar. 
    # We recommend using the latter in publications, see e.g., Sonison & Robinson (2018). 
    # You might also consider much more powerful differential testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_plot_low_verbo.png", show=False)


    ##save the result 
    mg_outfile_mid = "pbmc3k_4.h5ad"
    adata.write(mg_outfile_mid)


    ##############################################################################################
    #       NOT ENTIRELY SURE OF WHAT THE ADATA OBJECT IS DOING? WHY ARE WE SAVING IT?
    #       SURELY JUST ASSIGNING IT AS A DIFFERENT OBJECT WITH A LOGICAL NAME MAKES MORE SENSE?
    ##############################################################################################   

    #As an alternative, let us rank genes using logistic regression.For instance, this has been suggested by Natranos et al. (2018). 
    #The essential difference is that here, we use a multi-variate appraoch whereas conventional differential tests are uni-variate. Clark et al. (2014) has more details. 
    sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_log_reg.png", show=False)
    ######## not sure if this will work was getting an iteration error?? #####

    # With exceptions of IL7R, which is only found by the t-test and FCER1A, which is only found by the other two appraoches, 
    # all marker genes are recovered in all approaches.

    # Define a list of marker genes for later reference.
    marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14','LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1','FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

    #Reload the object that has been save with the Wilcoxon Rank-Sum test result.
    adata = sc.read(results_file_4)
    #Show the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

    # Get a table with scroes and groups 
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(5)

    #Compare to a single cluster:
    sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20, save="_single_cluster.png", show=False)

    ## If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.
    sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8, save="1.png", show=False)


    #Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):
    adata = sc.read(mg_outfile_mid)

    sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8, save="2.png", show=False)

    #If you want to compare a certain gene across groups, use the following.

    sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden', save="_compare.png", show=False)
    #########################################################################################################
    #             will need to make it so these genes are not hard coded...                                 #
    #########################################################################################################

    #Actually mark the cell types.
    new_cluster_names = [
        'CD4 T', 'CD14 Monocytes',
        'B', 'CD8 T',
        'NK', 'FCGR3A Monocytes',
        'Dendritic', 'Megakaryocytes']
    adata.rename_categories('leiden', new_cluster_names)

    ####################    REASON FOR PDF RATHER THAN PNG/JPEG???    ######################
    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='_leiden_cell_types.pdf',show=False)

    ##Now that we annotated the cell types, let us visualize the marker genes.
    sc.pl.dotplot(adata, marker_genes, groupby='leiden', save='.png', show=False);

    ##There is also a very compact violin plot.
    #########################################################################################################
    #           getting errors when making this plot!!! asked Jake already                                  #
    #########################################################################################################
    #sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', save='.png', show=False);

    #During the course of this analysis, the AnnData accumlated the following annotations.
    mg_outfile_end = "pbmc3k_4_end.h5ad"
    adata.write(mg_outfile_end, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

    #If you want to share this file with people who merely want to use it for visualization, 
    # a simple way to reduce the file size is by removing the dense scaled and corrected data matrix.
    # The file still contains the raw data used in the visualizations in adata.raw.
    adata.raw.to_adata().write(mg_outfile_withoutx)

if __name__ == "__main__":
    fire.Fire(marker_genes)