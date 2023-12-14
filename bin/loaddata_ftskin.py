#!/usr/bin/env python3
## script 2 of the pipeline

#Read in SoupX output 
#kernel clone_20230628_ cellrank 

import sys,os
import scvi
import anndata as ad
import matplotlib
import seaborn as sns
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
#from collections import Counte
#from fitter import Fitter, get_common_distributions, get_distributions
#import utils as eu
############################################################################
import datetime

#to set unique date for outputs
calc_date = datetime.datetime.now()
date=calc_date.strftime('%Y-%m-%d')
#date=calc_date.strftime('%d-%m-%Y')
date = date.replace('-', '')

sc.settings.verbosity = 2

# Note: package versions can be located in dependencies at the bottom of the script

root_directory = '/Users/nlg143/projects/sc-analysis-nf/output/'
out_directory = '/Users/nlg143/projects/sc-analysis-nf/output/'

sample = 'FCAImmP7241240'

adata_list = []

n = str(sample)
print(f'Data to load: {n}')
data = sc.read_10x_mtx(root_directory + n + 'soupX_filt', var_names = 'gene_symbols', make_unique=True, cache=False)
#(root_directory + n + '/output/Gene/filtered/', Method='mtx')
data.obs['sanger_id'] = n
adata_list.append(data)
print('\033[1m' + f'Data {n} loaded successfully. Data shape is {data.shape}' + '\033[0m', "\n")
print("="*100,"\n")

adata = sc.AnnData.concatenate(*adata_list, join = 'outer', batch_key = None, batch_categories = None, index_unique = None)

adata.obs_names_make_unique

adata.obs['sanger_id'].unique()

adata.obs['index1'] = adata.obs.index +'-1-'+ adata.obs['sanger_id']

adata.write('Fskin_0.1.1_raw')
#output/'+n+'soupX_filt/



