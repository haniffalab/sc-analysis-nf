#!/usr/bin/env python3

import fire
import scanpy as sc
import pandas as pd
import numpy as np


def marker_genes(mg_adata_input:str, text_in:str, mg_outfile_end:str, mg_outfile_withoutx:str ):
    #read in output file from third script 
    adata = sc.read_h5ad(mg_adata_input)
    results_file = 'pbmc3k.h5ad'  # the file that will store the analysis result
    ###################################### Finding marker genes #############################################
    #Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute 
    # of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_4.png",show=False)
    ########################################################################################################
    #                         RuntimeWarning: invalid value encountered in log2                            #
    #           an extra graph in my plot?? also is all of the script inplace=true?                        # 
    ########################################################################################################

    sc.settings.verbosity = 2  # reduce the verbosity
    # The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar. 
    # We recommend using the latter in publications, see e.g., Sonison & Robinson (2018). 
    # You might also consider much more powerful differential testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_wilcoxon_4.png", show=False)


    ##save the result 
    #mg_outfile_mid = "pbmc3k_4.h5ad"
    #adata.write(mg_outfile_mid)
    results_file = adata.copy()
    results_file2 = adata.copy()  

    ##############################################################################################
    #       NOT ENTIRELY SURE OF WHAT THE ADATA OBJECT IS DOING? WHY ARE WE SAVING IT?
    #       SURELY JUST ASSIGNING IT AS A DIFFERENT OBJECT WITH A LOGICAL NAME MAKES MORE SENSE?
    ##############################################################################################   

    #As an alternative, let us rank genes using logistic regression.For instance, this has been suggested by Natranos et al. (2018). 
    #The essential difference is that here, we use a multi-variate appraoch whereas conventional differential tests are uni-variate. Clark et al. (2014) has more details. 
    sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_logistic_regression_4.png", show=False)
    ######## not sure if this will work was getting an iteration error?? #####

    # With exceptions of IL7R, which is only found by the t-test and FCER1A, which is only found by the other two appraoches, 
    # all marker genes are recovered in all approaches.

    # Define a list of marker genes for later reference.
    mgenes_list = pd.read_csv(text_in)
    marker_genes = list(mgenes_list) 



    #Reload the object that has been save with the Wilcoxon Rank-Sum test result.
    #adata = sc.read(mg_outfile_mid)
    ### instead of read in a saved file just use the object with wilcoxon results in 
    results_file

    #Show the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.
    pd.DataFrame(results_file.uns['rank_genes_groups']['names']).head(5)

    # Get a table with scores and groups 
    result = results_file.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(5)

    #Compare to a single cluster:
    sc.tl.rank_genes_groups(results_file, 'leiden', groups=['0'], reference='1', method='wilcoxon')
    sc.pl.rank_genes_groups(results_file, groups=['0'], n_genes=20, save="_single_cluster_4.png", show=False)

    ## If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.
    sc.pl.rank_genes_groups_violin(results_file, groups='0', n_genes=8, save="_violin_4.png", show=False)


    #Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):
    #adata = sc.read(mg_outfile_mid)
    ### instead reload object adata_wilcoxon and reassign back to adata 
    #results_file2

    sc.pl.rank_genes_groups_violin(results_file2, groups='0', n_genes=8, save="_violin_differential_expression_4.png", show=False)

    #If you want to compare a certain gene across groups, use the following.

    sc.pl.violin(results_file2, ['CST3', 'NKG7', 'PPBP'], groupby='leiden', save="_compare_across_groups_4.png", show=False)
    #########################################################################################################
    #             will need to make it so these genes are not hard coded...                                 #
    #########################################################################################################

    #Actually mark the cell types.
    new_cluster_names = [
        'CD4 T', 'CD14 Monocytes',
        'B', 'CD8 T',
        'NK', 'FCGR3A Monocytes',
        'Dendritic', 'Megakaryocytes']
    results_file2.rename_categories('leiden', new_cluster_names)

    ####################    REASON FOR PDF RATHER THAN PNG/JPEG???    ######################
    sc.pl.umap(results_file2, color='leiden', legend_loc='on data', title='', frameon=False, save='_leiden_cell_types_4.pdf',show=False)

    ##Now that we annotated the cell types, let us visualize the marker genes.
    sc.pl.dotplot(results_file2, marker_genes, groupby='leiden', save='leiden_marker_genes_4.png', show=False);

    ##There is also a very compact violin plot.
    #########################################################################################################
    #           getting errors when making this plot!!! asked Jake already                                  #
    #########################################################################################################
    sc.pl.stacked_violin(results_file2, marker_genes, groupby='leiden', save='leiden_marker_genes_4.png', show=False);

    #During the course of this analysis, the AnnData accumlated the following annotations.
    mg_outfile_end = "pbmc3k_4_end.h5ad"
    results_file2.write(mg_outfile_end, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

    #If you want to share this file with people who merely want to use it for visualization, 
    # a simple way to reduce the file size is by removing the dense scaled and corrected data matrix.
    # The file still contains the raw data used in the visualizations in adata.raw.
    results_file2.raw.to_adata().write(mg_outfile_withoutx)

if __name__ == "__main__":
    fire.Fire(marker_genes)