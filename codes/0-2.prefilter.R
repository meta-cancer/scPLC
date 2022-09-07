library(dplyr)
library(DropletUtils)
# library(BiocParallel)
library(SingleCellExperiment)
library(scater)
library(scRNAseq)
library(igraph)
library(scran)

args=commandArgs(T)

## coi and data file 
# coi <- 'A093_ICC'
coi <- args[1]



## Step 0: read 10X file ============
# data.dir <- paste0("/zhanglab/home/caoqi/Proj_01_scRNA/1.rename_cellranger//",coi,"/outs/filtered_feature_bc_matrix/")
# sce <- read10xCounts(data.dir, col.names = T)
# # rownames and colnames
# rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
# colnames(sce) <- paste0(coi,"_",1:ncol(sce))
# colData(sce)[,'Sample'] <- coi
# # save file
# saveRDS(sce, paste0('./',coi,'_raw.rds'))
# 
data.dir <- paste0("/zhanglab/home/caoqi/Proj_01_scRNA/2.sce/sce_raw_A124/")
sce <- readRDS(paste0(data.dir,coi,'_raw.rds'))

## creat log file
log.file <- data.frame('process'='sce_prefilter')
log.file$name <- coi
# log: cellranger count
log.file$cellranger_count <- ncol(sce)




## Step 1: remove empty droplet ===============

set.seed(100)
a <- colSums(counts(sce)); 
b <- sort(a)[2] # find the second least UMI count as the lower bound 
# emptyDrops
e.out <- emptyDrops(counts(sce), lower = b, BPPARAM=MulticoreParam())
# filter 
sce_noempty <- sce[,which(e.out$FDR <= 0.01)] # cutoff FDR <= 0.01
# summary empty droplet
tmp <- summary(e.out$FDR <= 0.01) # 0.01 is cutoff from qiming's cell paper; 0.001 is default cutoff  
tmp
# log: empty cell count
log.file$empty_count <- ncol(sce) - ncol(sce_noempty) 




## Step 2: primary filter ============
# cutoff: sum  < 30000, 500 < gene < 6000, mt_pct < 50 (2020-10-5)
# cutoff2: sum  < 30000, 100 < gene < 6000, mt_pct < 50 (2020-10-30)


# mitochodrial
gene_MT <- read.table("/genome/human.GRCh38/gene_name_type.MT.txt")
is.mito <- rownames(sce_noempty) %in% gene_MT[,2]
sce_noempty <- addPerCellQC(sce_noempty, subsets=list(Mito=is.mito))

df <- perCellQCMetrics(sce_noempty, subsets=list(Mito=is.mito)) # test QCMetrics
# filter isOutlier per cell with specific cutoff
qc.lib <- isOutlier(df$sum, log=TRUE, type="lower") | df$sum > 30000 
qc.nexprs <- isOutlier(df$detected, log=TRUE, type="lower") | df$detected < 500 | df$detected > 6000 ## minexp =500 !!!
qc.mito <- isOutlier(df$subsets_Mito_percent, type="higher") | df$subsets_Mito_percent > 50
attr(qc.lib, "thresholds")
attr(qc.nexprs, "thresholds")
attr(qc.mito, "thresholds")

discard <- qc.lib | qc.nexprs |  qc.mito
log.file[1,c("libSize",'NExprs','MitoProp','discard')] <- c(sum(qc.lib), sum(qc.nexprs), sum(qc.mito), sum(discard))

# add discard info to sce 
colData(sce_noempty)[,'discard'] <- discard

# rm discard cells
sce2 <- sce_noempty[,! sce_noempty$discard]

# plot checking 
sce_data <- sce_noempty
gridExtra::grid.arrange(
  plotColData(sce_data, y = 'sum', colour_by = 'discard')+ 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce_data, y = 'detected', colour_by = 'discard')+ 
    scale_y_log10() + ggtitle("Detected genes"),
  plotColData(sce_data, y = 'subsets_Mito_percent', colour_by = 'discard')+ 
    ggtitle("Mito percent"),
  ncol=2
)


plot(sce2$sum,sce2$subsets_Mito_percent)
plot(sce2$sum,sce2$detected)
plot(sce2$subsets_Mito_percent,sce2$detected)

# saveRDS(sce2, paste0(plc70[i,1],".filter.rds"))




## Step 3: rm doublet cells ========

# doubletCells score 
scores <- doubletCells(sce2)

# doublet cutoff 
cutoff <-  median(log10(scores+1)) + 3* mad( log10(scores+1) ) ## !!!!!!!

# cluster 
sce2log <- logNormCounts(sce2)
sce2log <- scater::runPCA(sce2log, scale_features=NULL, ntop = 1000, ncomponents=50)

g <- buildKNNGraph(sce2log, use.dimred="PCA", BPPARAM = MulticoreParam(), k=5)
clusters <- membership( cluster_louvain(g) )
# clusters <- igraph::cluster_fast_greedy(g)$membership


# outlier cluster
table(clusters)
boxplot(split(log10(scores+1), clusters))
abline(h=cutoff)

x <- split(log10(scores+1), clusters)
y <- as.data.frame( lapply(x, median) )

tmp <- names(y)[y > cutoff]
tmp <- sub("X", "", tmp)

# rm doublets 
sce_nodoublet <- sce2[,! clusters %in% tmp ]

# log: doubltes cells count 
log.file$doublet <- ncol(sce2) - ncol(sce_nodoublet)

log.file$retained_cells <- ncol(sce_nodoublet)



# Step 4: output and log file -=========

# save count file 
saveRDS(sce_nodoublet, paste0(coi,"_count.rds"))
write.table(log.file, paste0("log_prefilter.",coi,'.txt'))
