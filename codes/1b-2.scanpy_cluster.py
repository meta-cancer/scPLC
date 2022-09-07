import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
# import anndata2ri
from rpy2.robjects import r
import re

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_header()
# sc.settings.set_figure_params(dpi=80, facecolor='white')


## input and output file 
# input_data = 'Seurat_tmp.counts.h5ad'
input_data = sys.argv[1]
adata = sc.read_h5ad(input_data)

coi = re.sub('.h5ad','',input_data)

n_hvg = int(sys.argv[2])
n_PC = int(sys.argv[3])
n_res = float(sys.argv[4])


# n_hvg = 2000
# n_PC = 10
# n_res = 1


results_file = coi+'.scanpy.hvg'+str(n_hvg)+'_PC'+str(n_PC)+'_res'+str(n_res)+'.h5ad' # the file that will store the analysis results
results_fig = coi+'.scanpy.plot_hvg'+str(n_hvg)+'_PC'+str(n_PC)+'_res'+str(n_res)+'.png'

## Step 0: preprocessing 

# sc.pl.highest_expr_genes(adata, n_top=20, )


# # base filtering 
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# 
# # compute some metrics 
# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
# 
# # violin plot 
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)
#             
# # scatter plot 
# sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
#
#
# # split adata 
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

## Step 1: scale data 
# Total-count normalize to 10000
sc.pp.normalize_total(adata, target_sum=1e4)


# Logarithmize the data.
sc.pp.log1p(adata)

# identify highly-viriable genes 
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

# sc.pl.highly_variable_genes(adata)

# get back an AnnData of the object in .raw by calling .raw.to_adata().
adata.raw = adata


# split hvg genes 
adata = adata[:, adata.var.highly_variable]


# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. 
# Scale the data to unit variance.
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# sc.pp.regress_out(adata, [ 'pct_counts_mt'])
sc.pp.regress_out(adata, [ 'percent.mt'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)


### Step 2: Principal component analysis
# PCA
sc.tl.pca(adata)

# # pca variance 
# sc.pl.pca_variance_ratio(adata, log=True)
# sc.pl.pca(adata, color='CST3')

### Step 3: compute the neighborhood graph 
# neighbors 
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_PC)

# # bbknn
# sc.external.pp.bbknn(adata, batch_key='batch')  # running bbknn 1.3.6


# adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat


### Step 4: Embedding the neighborhood graph
# # paga
# sc.tl.paga(adata)
# sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
# sc.tl.umap(adata, init_pos='paga')

# umap 
sc.tl.umap(adata)
# sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])


### Step 5: Clustering the neighborhood graph
# # leiden 
# sc.tl.leiden(adata)


# louvain
sc.tl.louvain(adata, resolution=n_res)

# AnnData object with n_obs × n_vars = 963276 × 2000
#     obs: 'batch', 'louvain'
#     uns: 'pca', 'neighbors', 'umap', 'louvain'
#     obsm: 'X_pca', 'X_umap'
#     varm: 'PCs'


# plot 
mpl.rcParams['figure.figsize'] = (10,10)

sc.pl.umap(adata, color=['louvain','Sample'], legend_loc='on data')

mpl.pyplot.savefig(results_fig,dpi=300)


### final : save file 
# write data in h5ad format, and read data by ReadH5AD function
adata.write_h5ad(results_file)
