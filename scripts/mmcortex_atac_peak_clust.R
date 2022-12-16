library(metacell)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
devtools::load_all("~/src/mcATAC/")
gset_genome("mm10")


atac_mc <- readRDS('./data/mmcortex_atac_mc.rds')

pk_cl <- gen_atac_peak_clust(atac_mc = atac_mc, k = 64)
saveRDS(pk_cl, './data/mmcortex_atac_filt_peak_clustering_k=64.rds')
atac_mc@peaks$pk_cl_k=64 <- pk_cl