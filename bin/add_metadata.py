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


def add_metadata(
    add_meta_in: str, metadata_in: str, samples: str, add_metadata_out: str
):
    ## Read in raw h5ad object
    adata = sc.read(add_meta_in)
    ## had to add these in from the first part of the load_data.py script ##
    sample = samples
    n = str(sample)
    adata.obs["sanger_id"] = n

    adata.obs_names_make_unique
    adata.obs["sanger_id"].unique()
    adata.obs["index1"] = adata.obs.index + "-1-" + adata.obs["sanger_id"]
    ########################################################################
    ## Read in metadata & set index
    old_raw_metadata = pd.read_csv(metadata_in, low_memory=False)
    old_raw_metadata.set_index("index", inplace=True)

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
    del adata.obs["sanger_id_y"]
    del adata.obs["sanger_id_x"]

    # remove?
    adata.obs.sorting.unique()
    adata.obs["sample"].unique()

    adata.obs["pcw"] = adata.obs["pcw"].astype(str)

    adata.write(add_metadata_out)


if __name__ == "__main__":
    fire.Fire(add_metadata)
