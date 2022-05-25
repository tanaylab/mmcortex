library(parallel)
SEED = 1337
set.seed(SEED)
wd = "/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex"
setwd(wd)
pmc_mot_filt = readRDS('./data/peak_mc_motifs_200bp_filt.rds')
cols_motifs = grep('chrom|start|end$', colnames(pmc_mot_filt), inv=T)
enh_subsamp = readRDS('./data/enh_subsamp.rds')
pmc_mot_subsamp = pmc_mot_filt[enh_subsamp,cols_motifs]
kis = c(10, 25, seq(100, 500, 50))
km_vec = mclapply(kis, function(ki) {print(ki); return(tglkmeans::TGL_kmeans(t(pmc_mot_subsamp), k = ki))}, mc.cores = length(kis))
saveRDS(km_vec, './data/test_km_on_pmc_mot_subsamp.rds')