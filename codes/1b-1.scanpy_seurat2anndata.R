library(Seurat)
library(sceasy)
library(reticulate)

### input and output file ### 
# input: seurat.rds 
args = commandArgs(T)
input_file <- args[1]

# output: anndata.h5ad 
res_file <- gsub('.rds', '.rds.h5ad', input_file)


# read file
combined_seu <- readRDS(input_file)

## Step 0: configure condaenv ------

use_python(python = '/zhanglab/home/caoqi/anaconda3/bin/python', required = T)
use_condaenv('EnvironmentName', conda = '/zhanglab/home/caoqi/anaconda3/bin/anaconda')

# Sys.getenv()
# Sys.setenv()
# py_config()
# py_available()
# py_module_available('loompy')
loompy <- reticulate::import('loompy')
anndata <- reticulate::import('anndata')
# loompy <- reticulate::import_from_path(module = 'loompy', path = '/zhanglab/home/caoqi/anaconda3/lib/python3.7/site-packages/')
# anndata <- reticulate::import_from_path(module = 'anndata', path = '/zhanglab/home/caoqi/anaconda3/lib/python3.7/site-packages/anndata/')

print(paste(Sys.time(),'condaenv is built'))


## Step 1: format transform -----

# Seurat to AnnData
sceasy::convertFormat(combined_seu, from="seurat", to="anndata", outFile=res_file)

print(paste(Sys.time(),'Seurat to AnnData completed'))


# # AnnData to Seurat
# sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
#                       outFile='filename.cds')
