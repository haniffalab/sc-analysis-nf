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



def load_mtx():

    adata = sc.read_10x_mtx(
        "./", var_names="gene_symbols", make_unique=True, cache=False
    )

    adata.write("adata.h5ad")
    
if __name__ == "__main__":
    fire.Fire(load_mtx)
