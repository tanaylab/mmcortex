library(prego)
library(misha.ext)
gset_genome('mm10')
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
# scdb_init(base_dir = './scdb', force_reinit = T)
# nm <- "pl_cort_feat_peaks"

# mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
# cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
#                       'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
# cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
#                                                                   rep(s, length(which(mcmd$cell_type == s))))))                         

####################
### Don't remove
## How to create feat_peak mat object
load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/multi_mmcortex.Rda",v=T)
# load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_mmcortex.Rda",v=T)
# load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_annot_mmcortex.Rda",v=T)

# #all intervals
# head(multi_model$atac_intervs)
# #those in clusters annotated as variable
# feat_peak = multi_model$atac_intervs[acn$f_var_peak,]
feat_peak = multi_model$atac_intervs
# scc <- scc_read('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads/')

# mcc <- scc_project_on_mc(sc_counts=scc, cell_to_metacell=c2mc)
# feat_peak$peak_name <- peak_names(feat_peak)
# mcc_write(mcc, './output/mcatac/mmcortex_mcc_feat_peak')
### Don't remove
####################

# mca <- mcc_to_mcatac(mc_counts=mcc, peaks=feat_peak,metadata=mcmd)
# saveRDS(mca, './output/mcatac/mmcortex_mcatac_feat_peaks.rds')
mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks.rds')
a_legc <- log2(1e-5 + mca@egc)

# km_a_legc <- tglkmeans::TGL_kmeans(as.matrix(a_legc), k = 80, seed = 1337)
# saveRDS(km_a_legc, './output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
km_a_legc <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')

new_clust_vec <- km_a_legc$cluster
seq_coords <- mca@peaks[,c('chrom', 'start', 'end')]

coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)

seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

prego::set_parallel(10)
aa <- getOption('prego.parallel.nc')
print(aa)
options(error = recover)

rm(list = ls()[!(ls() %in% c('seqs_all', 'new_clust_vec', "ALLGENOME", "GINTERVS", "GROOT", "GTRACKS", "GWD"))])
print(ls())

res <- prego::regress_pwm.clusters(sequences=seqs_all, clusters = new_clust_vec, 
                                        use_sge = T, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
print(warnings())
saveRDS(res, './output/sequence_modeling/mmcortex_mcatac_rna_match_feat_peak_res.rds')