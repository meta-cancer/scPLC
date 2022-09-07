library(reticulate)

# create a new environment 
conda_create("r-reticulate")


# install SciPy
conda_install("r-reticulate", "scipy")
# conda_install("r-reticulate", "anndata==0.6.19") # fail to find anndata==0.6.19 in bioconda
conda_install("r-reticulate", "loompy==2.0.17")

# import SciPy (it will be automatically discovered in "r-reticulate")
scipy <- import("scipy")

# indicate that we want to use a specific condaenv
# Sys.setenv(RETICULATE_PYTHON = "/zhanglab/home/caoqi/anaconda3/bin/python")
# use_python(python = '/zhanglab/home/caoqi/anaconda3/bin/python3', required = T)
use_condaenv("r-reticulate")

# import SciPy (will use "r-reticulate" as per call to use_condaenv)
scipy <- import("scipy")
anndata <- import("anndata")



library(Seurat)
library(sceasy)

args = commandArgs(T)

### input and output file ### 
# input: anndata.h5ad 
# input_file <- args[1]
# input_file <- 'Seurat_tmp.counts.rds'
input_file <- 'Seurat_tmp.counts.h5ad'
# input_file <- 'Seurat_tmp.counts.scanpy.hvg2000_PC10_res1.h5ad'


coi <- gsub('.h5ad', '', input_file)

setwd('~/Proj_01_scRNA/2.sce/')
# output: seurat.cds
res_file <- paste0(coi,'.cds')


# read file
# combined_seu <- readRDS(input_file)
# combined_seu <- ReadH5AD(input_file)

## Step 0: configure condaenv 

# use_python(python = '/zhanglab/home/caoqi/anaconda3/bin/python', required = T)
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

scanpy <- reticulate::import('scanpy')

print(paste(Sys.time(),'condaenv is built'))

## Step 1: read input 
combined_seu = anndata$read_h5ad(input_file)


## Step 2: format transform and output 


# Seurat to AnnData
# sceasy::convertFormat(combined_seu, from="seurat", to="anndata", outFile=res_file)
# 
# print(paste(Sys.time(),'Seurat to AnnData completed'))


# AnnData to Seurat
sceasy::convertFormat(combined_seu, from="anndata", to="seurat", outFile=res_file)


