library(Seurat)
library(DoubletFinder)
# library(tidyverse)
library(metacell)
scdb_init('scdb', force_reinit = T)
wd = '/home/feshap/raid/proj/mmcortex'
data_path = file.path(wd, 'scrna_data_gz')


# mat = scdb_mat('cort5')
# sc_cort = colnames(mat@mat)
# 
# write(x = sc_cort, file = './data/sc_cort.txt')

dirs = list.dirs(data_path)

# sample_names = tail(head(tail(lapply(stringr::str_split(dirs, '\\/'), function(x) dplyr::last(x)), -2), -1), -2)
sample_names = sapply(stringr::str_split(grep('rep', list.dirs(data_path), v=T), '\\/'), function(x) x[[length(x)]])
sample_names

# mats_all = lapply(head(tail(dirs, -2), -1), Read10X)
mats_all = lapply(grep('rep', list.dirs(data_path), v=T), Seurat::Read10X)

dim_all = lapply(mats_all, dim)
dim_all

## Multiplet data from 
##        https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
##        https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php

df_10x = vroom::vroom('./data/df_10x_doublets.tsv')

doub_lm <- lm(prc ~ recovered, data=df_10x)
doub_lm

samp_n_df = data.frame(cbind(as.character(sample_names), as.numeric(purrr::map(dim_all, 2))))
colnames(samp_n_df) = c('sample', 'n')
samp_n_df$inf_prc = 0.01*(doub_lm$coefficients[[1]] + doub_lm$coefficients[[2]] * as.numeric(samp_n_df[,'n']))
samp_n_df

readr::write_tsv(samp_n_df, './data/samp_n_df_pl.tsv')

# mats_all = lapply(head(tail(dirs, -2), -1), Read10X)

remove_mito_small = function(mats_merge) {
    nms = rownames(mats_merge)
    mito_genes = grep('^mt-', nms)
    mito_umis = Matrix::colSums(mats_merge[mito_genes,])
    tot_umis = Matrix::colSums(mats_merge)
    mito_frac = mito_umis/tot_umis
    mats_no_mito = mats_merge[,mito_frac < 0.15]
    cs = Matrix::colSums(mats_no_mito)
    mats_filt_size = mats_no_mito[,cs >= 1.5e+03 & cs <= 3e+04]
    return(mats_filt_size)
}

mats_all = lapply(mats_all, remove_mito_small)

# mats_all = lapply(mats_all, function(x) x = x[,dimnames(x)[[2]] %in% sc_cort])

# mats_all = tail(mats_all, -2)

saveRDS(object = mats_all, file = './data/mats_all_pl.rds')

mats_all = readRDS('./data/mats_all_pl.rds')

# seu_mmc = purrr::map2(tail(mats_all, -2), sample_names, function(x,y) CreateSeuratObject(counts = x, project = y))
seu_mmc = purrr::map2(mats_all, sample_names, function(x,y) CreateSeuratObject(counts = x, project = y))

purrr::map2(seu_mmc, sample_names, function(x,y) saveRDS(x, file.path(wd, 'pl_doublet', glue('{y}_seu.rds'))))

get_bcmvn = function(mmc){
    sweep.res.list_mmc <- paramSweep_v3(mmc, PCs = 1:30, sct = FALSE)
    sweep.stats_mmc <- summarizeSweep(sweep.res.list_mmc, GT = FALSE)
    bcmvn_mmc <- find.pK(sweep.stats_mmc)
    return(bcmvn_mmc)
}

for (i in 1:length(seu_mmc)) {
    seu_mmc[[i]] = NormalizeData(seu_mmc[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
    seu_mmc[[i]] = FindVariableFeatures(seu_mmc[[i]], selection.method = "vst", nfeatures = 2000)
    seu_mmc[[i]] = ScaleData(seu_mmc[[i]])
    seu_mmc[[i]] = RunPCA(seu_mmc[[i]], features = VariableFeatures(object = seu_mmc[[i]]))
    seu_mmc[[i]] = RunUMAP(seu_mmc[[i]], features = VariableFeatures(object = seu_mmc[[i]]))
    seu_mmc[[i]] <- FindNeighbors(seu_mmc[[i]], dims = 1:30)
    seu_mmc[[i]] <- FindClusters(seu_mmc[[i]], resolution = 0.6)
}


bcmvn_vec = lapply(seu_mmc, get_bcmvn)

calculate_nExp = function(pbmc, inf_prc, i) {
    print(inf_prc[[i]])
    homotypic.prop <- modelHomotypic(pbmc[[i]]@meta.data$seurat_clusters)
    nExp_poi <- round(0.75 * inf_prc[[i]] * nrow(pbmc[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
    return(nExp_poi.adj)
}

saveRDS(object = seu_mmc, file = './data/seu_mmc_pl.rds')

saveRDS(object = bcmvn_vec, file = './data/bcmvn_vec_pl.rds')

seu_mmc = readRDS(file = './data/seu_mmc_pl.rds')
bcmvn_vec = readRDS(file.path('./data/bcmvn_vec_pl.rds'))


samp_n_df = vroom::vroom('./data/samp_n_df_pl.tsv')

nepa_vec = lapply(seq_along(seu_mmc), function(x, p, i) {calculate_nExp(x, p, i)}, x = seu_mmc, p = samp_n_df$inf_prc)

pk_max = purrr::map(seq_along(bcmvn_vec), 
                    function(i) bcmvn_vec[[i]][which.max(bcmvn_vec[[i]][2:(dim(bcmvn_vec[[i]])[[1]] - 1),'BCmetric'])+1,'pK'])
pk_max = purrr::map(pk_max, function(x) as.numeric(levels(x)[x]))
pk_max

for (i in 1:length(seu_mmc)) {
    print('sample i')
    seu_mmc[[i]] = doubletFinder_v3(seu_mmc[[i]], PCs = 1:30, pN = 0.25, pK = pk_max[[i]], 
                                    nExp = nepa_vec[[i]], reuse.pANN = FALSE, sct = FALSE)
}


            
# doublet_preds = lapply(seu_mmc, function(x) x@meta.data[, dim(x@meta.data)[[2]]]) %>% unlist
doublet_preds = lapply(seu_mmc, function(x) x@meta.data[, grep('DF', colnames(x@meta.data),v=T)]) %>% unlist
doublet_names = lapply(seu_mmc, function(x) rownames(x@meta.data)) %>% as.matrix %>% unlist

df_pred = data.frame(cbind(doublet_names, doublet_preds))
head(df_pred)

saveRDS(df_pred, file.path(wd, 'data', 'df_pred_pl.rds'))



