#!/usr/bin/env python3

import fire
import scanpy as sc
import pandas as pd
import numpy as np

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")


def read_data(path: str, outfile: str):
    # results_file = 'pbmc3k.h5ad'  # the file that will store the analysis result

    # Read in the count matrix into an AnnData object, which holds many slots for annotations and different representations of the data.
    # It also comes with its own HDF5-based file format: .h5ad.
    adata = sc.read_10x_mtx(
        path,  # the directory with the `.mtx` file
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=True,
    )  # write a cache file for faster subsequent reading

    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    adata

    adata.write_h5ad(outfile)


if __name__ == "__main__":
    fire.Fire(read_data)
