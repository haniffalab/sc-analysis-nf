#!/usr/bin/env python3
## script 2 of the pipeline

# Read in SoupX output
# kernel clone_20230628_ cellrank
import fire
import sys, os
import anndata as ad
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import numpy as np3
import pandas as pd
import scanpy as sc
import numpy.random as random
import scipy.sparse
import anndata as ad
import numpy as np

################## unsure why these are commented out yet ##################
# from collections import Counte
# from fitter import Fitter, get_common_distributions, get_distributions
# import utils as eu
############################################################################
import datetime


def load_data(load_in: str, load_raw_out: str, metadata_in: str, load_meta_out: str):
    # to set unique date for outputs
    calc_date = datetime.datetime.now()
    date = calc_date.strftime("%Y-%m-%d")
    # date=calc_date.strftime('%d-%m-%Y')
    date = date.replace("-", "")

    sc.settings.verbosity = 2

    # Note: package versions can be located in dependencies at the bottom of the script

    root_directory = "/sc-analysis-nf/output/"
    out_directory = "/sc-analysis-nf/output/"

    sample = "FCAImmP7241240"

    adata_list = []

    n = str(sample)
    print(f"Data to load: ${load_in}")
    data = sc.read_10x_mtx(
        "./", var_names="gene_symbols", make_unique=True, cache=False
    )
    # data = sc.read_10x_mtx(root_directory + n + 'soupX_filt', var_names = 'gene_symbols', make_unique=True, cache=False)
    # (root_directory + n + '/output/Gene/filtered/', Method='mtx')
    data.obs["sanger_id"] = n
    adata_list.append(data)
    print(
        "\033[1m"
        + f"Data {n} loaded successfully. Data shape is {data.shape}"
        + "\033[0m",
        "\n",
    )
    print("=" * 100, "\n")

    adata = sc.AnnData.concatenate(
        *adata_list,
        join="outer",
        batch_key=None,
        batch_categories=None,
        index_unique=None,
    )

    adata.obs_names_make_unique

    adata.obs["sanger_id"].unique()

    adata.obs["index1"] = adata.obs.index + "-1-" + adata.obs["sanger_id"]

    adata.write(load_raw_out)
    # output/'+n+'soupX_filt/
    adata.obs["index1"]

    ##############################################################################################
    ## Read in metadata & set index
    old_raw_metadata = pd.read_csv(metadata_in, low_memory=False)
    old_raw_metadata.set_index("index", inplace=True)
    old_raw_metadata
    # remove? checking structure
    old_raw_metadata.sanger_id.unique()

    # making a copy of metadata
    old_raw_metadata["index1"] = old_raw_metadata.index.values.copy()
    # renaming the columns
    raw_metadata = old_raw_metadata.drop(
        labels={
            "n_genes",
            "scrublet_score",
            "bh_doublet_pval",
            "cell_caller",
            "cluster_scrublet_score",
            "doublet_pval",
            "mt_prop",
            "n_counts",
            "n_genes",
            "annot",
            "hierarchy1",
        },
        axis="columns",
    )
    # remove? checking structure
    raw_metadata["sanger_id"].unique
    # remove?
    adata.obs

    # Convert "adata.obs" to a DataFrame
    adata_obs_df = adata.obs.reset_index()
    # Merge "raw_metadata" and "adata_obs_df" on the "index1" column
    merged_df = raw_metadata.merge(adata_obs_df, on="index1")

    # Update the "adata" object with the merged metadata
    adata.obs = adata.obs.merge(
        merged_df, left_on="index1", right_on="index1", how="left"
    )

    # Now, "adata" contains the merged metadata in AnnData format with the "index1" column preserved

    adata.obs = adata.obs.set_index("index1")

    # remove?
    len(adata.obs.index.unique())

    del adata.obs["sanger_id_y"]
    del adata.obs["sanger_id_x"]

    adata.obs.sorting.unique()
    adata.obs["sample"].unique()

    #############################################################################################################
    ########### this section doesn't seem to be needed but in jupyter notebook isnt returning an error ##########
    #############################################################################################################
    # Samples missing from metadata
    # 4834STDY7002871 4834STDY7002872
    # new samples
    #'HCA_rFSKI14449635', 'HCA_rFSKI14449636', 'HCA_rFSKI14449637', 'HCA_rFSKI14449638', 'HCA_rFSKI14449639', 'HCA_rFSKI14449640', 'HCA_rFSKI14449641', 'HCA_rFSKI14449642' 'WS_wEMB13942138', 'WS_wEMB13942139', 'WS_wEMB13942140', 'WS_wEMB13942141', 'WS_wEMB13942142', 'WS_wEMB13942143', 'WS_wEMB13942144', 'WS_wEMB13942145', 'WS_wEMB13942146', 'WS_wEMB13942147', 'WS_wEMB13942148', 'WS_wEMB13942149', 'WS_wEMB13942150', 'WS_wEMB13942151', 'WS_wEMB13942152', 'WS_wEMB13942153'
    # Fill up filling data on new samples

    # Find the indices of samples where sanger_id == "12343"
    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(
            [
                "HCA_rFSKI14449639",
                "HCA_rFSKI14449640",
                "HCA_rFSKI14449641",
                "HCA_rFSKI14449642",
                "HCA_rFSKI14449635",
                "HCA_rFSKI14449636",
                "HCA_rFSKI14449637",
                "HCA_rFSKI14449638",
            ]
        )
    ]

    # Update the 'X' column with "YXYX" for the selected samples
    adata.obs.loc[samples_to_fill_indices, "donor"] = "F220"
    adata.obs.loc[samples_to_fill_indices, "gender"] = "female"
    adata.obs.loc[samples_to_fill_indices, "pcw"] = "17"
    adata.obs.loc[samples_to_fill_indices, "chemistry"] = "SC5P-R2"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449635"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45P"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45P"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449636"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34N"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34N"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449637"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34N"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34N"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449638"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34P"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34P"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449639"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45P"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45P"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449640"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34N"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34N"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449641"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34N"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34N"

    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(["HCA_rFSKI14449642"])
    ]
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "CD45NCD34P"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_CD45NCD34P"

    #### Next donor round
    samples_to_fill_indices = adata.obs.index[
        adata.obs["sanger_id"].isin(
            [
                "WS_wEMB13942138",
                "WS_wEMB13942139",
                "WS_wEMB13942140",
                "WS_wEMB13942141",
                "WS_wEMB13942142",
                "WS_wEMB13942143",
                "WS_wEMB13942144",
                "WS_wEMB13942145",
                "WS_wEMB13942146",
                "WS_wEMB13942147",
                "WS_wEMB13942148",
                "WS_wEMB13942149",
                "WS_wEMB13942150",
                "WS_wEMB13942151",
                "WS_wEMB13942152",
                "WS_wEMB13942153",
            ]
        )
    ]
    # Update the 'X' column with "YXYX" for the selected samples
    adata.obs.loc[samples_to_fill_indices, "donor"] = "F217"
    adata.obs.loc[samples_to_fill_indices, "gender"] = "female"
    adata.obs.loc[samples_to_fill_indices, "pcw"] = "17"
    adata.obs.loc[samples_to_fill_indices, "chemistry"] = "SC5P-R2"
    adata.obs.loc[samples_to_fill_indices, "sorting"] = "Total"
    adata.obs.loc[samples_to_fill_indices, "chemistry_sorting"] = "SC5P-R2_Total"

    # Some cells have missing information in metadata column. We fill up information based on the samples with information

    # Some NA values, make labelling consistent
    # Define a list of columns to fill missing values in
    columns_to_fill = [
        "pcw",
        "chemistry",
        "sorting",
        "chemistry_sorting",
        "donor",
        "gender",
    ]  # Add more columns if needed
    # 'gender','sample', add back later
    # Replace 'nan' values with NaN in all specified columns
    for column in columns_to_fill:
        adata.obs[column] = adata.obs[column].replace("nan", np.nan)

    # Group your data by 'sanger_id' and fill NaN values with the mode
    for sanger_id, group in adata.obs.groupby("sanger_id"):
        for column in columns_to_fill:
            mode_value = group[column].mode().values[0]  # Calculate the mode
            adata.obs.loc[group.index, column] = adata.obs.loc[
                group.index, column
            ].fillna(mode_value)

    # Define a list of columns to fill missing values in
    columns_to_fill = ["pcw", "chemistry", "sorting", "chemistry_sorting", "donor"]
    # Add more columns if needed (e.g., 'gender', 'sample')

    # Replace 'nan' values with NaN in all specified columns
    for column in columns_to_fill:
        adata.obs[column] = adata.obs[column].replace("nan", np.nan)

    # Group your data by 'sanger_id' and fill NaN values with the mode
    for sanger_id, group in adata.obs.groupby("sanger_id"):
        for column in columns_to_fill:
            # Calculate the mode and check if it's empty
            mode_values = group[column].mode()
            if not mode_values.empty:
                mode_value = mode_values.values[0]
                adata.obs.loc[group.index, column] = adata.obs.loc[
                    group.index, column
                ].fillna(mode_value)
            else:
                print(
                    f"Mode is empty for '{column}' in group '{sanger_id}'. Handle this case as needed."
                )

    adata.obs["pcw"] = adata.obs["pcw"].astype(str)

    adata.write(load_meta_out)


if __name__ == "__main__":
    fire.Fire(load_data)
