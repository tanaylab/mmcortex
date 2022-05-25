library(Seurat)
library(Matrix)
library(glue)
library(DoubletFinder)
library(tidyverse)

wd = '/home/feshap/raid/proj/mmcortex'
data_path = file.path(wd, 'scrna_data_gz')


dirs = list.dirs(data_path)

mats_all = lapply(head(tail(dirs, -2), -1), Read10X)

remove_mito_small = function(mats_merge) {
    nms = rownames(mats_merge)
    mito_genes = grep('^mt-', nms)
    mito_umis = Matrix::colSums(mats_merge[mito_genes,])
    tot_umis = Matrix::colSums(mats_merge)
    mito_frac = mito_umis/tot_umis
    mats_no_mito = mats_merge[,mito_frac < 0.15]
    mats_filt_small = mats_no_mito[,Matrix::colSums(mats_no_mito) >= 4000]
    return(mats_filt_small)
}

mats_all = lapply(mats_all, remove_mito_small)

sample_names = head(tail(lapply(stringr::str_split(dirs, '\\/'), function(x) dplyr::last(x)), -2), -1)

seu_mmc = purrr::map2(mats_all, sample_names, function(x,y) CreateSeuratObject(counts = x, project = y))

purrr::walk2(seu_mmc, sample_names, function(x,y) saveRDS(x, file.path(wd, 'doublet_seu', glue('{y}_seu.rds'))))

for (i in 1:length(seu_mmc)) {
    seu_mmc[[i]] = NormalizeData(seu_mmc[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
    seu_mmc[[i]] = FindVariableFeatures(seu_mmc[[i]], selection.method = "vst", nfeatures = 2000)
    seu_mmc[[i]] = ScaleData(seu_mmc[[i]])
    seu_mmc[[i]] = RunPCA(seu_mmc[[i]], features = VariableFeatures(object = seu_mmc[[i]]))
    seu_mmc[[i]] = RunUMAP(seu_mmc[[i]], features = VariableFeatures(object = seu_mmc[[i]]))
}

get_bcmvn = function(mmc){
    sweep.res.list_mmc <- paramSweep_v3(mmc, PCs = 1:10, sct = FALSE)
    sweep.stats_mmc <- summarizeSweep(sweep.res.list_mmc, GT = FALSE)
    bcmvn_mmc <- find.pK(sweep.stats_mmc)
    return(bcmvn_mmc)
}

bcmvn_vec = lapply(seu_mmc, get_bcmvn)

purrr::walk2(seu_mmc, sample_names, function(x,y) saveRDS(x, file.path(wd, 'doublet_seu', glue('{y}_seu.rds'))))

saveRDS(bcmvn_vec, file.path(wd, 'doublet_seu', 'bcmvn_vec.rds'))

bcmvn_vec = readRDS(file.path(wd, 'doublet_seu', 'bcmvn_vec.rds'))

seu_mmc = lapply(grep('seu', file.path(wd, list.files(file.path(wd, 'doublet_seu'))), v=T), readRDS)

purrr::map(seq_along(bcmvn_vec), function(i) plot(bcmvn_vec[[i]][,'pK'], bcmvn_vec[[i]][,'BCmetric'], main = sample_names[[i]]))

calculate_nExp = function(pbmc) {
    pbmc = pbmc
    pbmc <- FindNeighbors(pbmc, dims = 1:15)
    pbmc <- FindClusters(pbmc, resolution = 0.6)
    homotypic.prop <- modelHomotypic(pbmc@meta.data$seurat_clusters)
    nExp_poi <- round(0.075*nrow(pbmc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    return(nExp_poi.adj)
}

nepa_vec = lapply(seu_mmc, calculate_nExp)

pk_max = purrr::map(seq_along(bcmvn_vec), 
                    function(i) bcmvn_vec[[i]][which.max(bcmvn_vec[[i]][2:(dim(bcmvn_vec[[i]])[[1]] - 1),'BCmetric'])+1,'pK'])
pk_max = map(pk_max, function(x) as.numeric(levels(x)[x]))
# pk_max

for (i in 1:length(seu_mmc)) {
    print('sample i')
    seu_mmc[[i]] = doubletFinder_v3(seu_mmc[[i]], PCs = 1:15, pN = 0.25, pK = pk_max[[i]], 
                                    nExp = nepa_vec[[i]], reuse.pANN = FALSE, sct = FALSE)
}

purrr::walk2(seu_mmc, sample_names, function(x,y) saveRDS(x, file.path(wd, 'doublet_seu', glue('{y}_seu.rds'))))