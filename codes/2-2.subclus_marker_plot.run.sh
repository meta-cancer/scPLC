#!/usr/bin/env bash
# usage: bash script.sh [seurat|sce] file.rds
input_format=$1
input_file=$2
output=`basename $input_file`
dir=`pwd`
R -e "rmarkdown::render('/zhanglab/home/caoqi/Proj_01_scRNA/2.sce/04.subclus_marker_plot.Rmd',output_file=paste0('$dir/$output.',Sys.Date(),'.html'))" --args $input_format $dir/$input_file
