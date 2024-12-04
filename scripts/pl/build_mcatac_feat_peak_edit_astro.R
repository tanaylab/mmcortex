library(metacell)
devtools::load_all("~/src/mcATAC/")
gset_genome('mm10')
SEED = 1337
set.seed(SEED)
wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
nm <- "pl_cort_feat_peaks"

mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))                         

####################
### Don't remove
## How to create feat_peak mat object
load("/net/mraid20/export/tgdata/users/atanay/proj/mmcortex/work0922/data/multi_mmcortex.Rda",v=T)
# load("/net/mraid20/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_mmcortex.Rda",v=T)
# load("/net/mraid20/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_annot_mmcortex.Rda",v=T)

# #all intervals
# head(multi_model$atac_intervs)
# #those in clusters annotated as variable
# feat_peak = multi_model$atac_intervs[acn$f_var_peak,]
feat_peak = multi_model$atac_intervs
scc <- scc_read('/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads_28122022/')
c2mc <- readr::read_tsv("./output/mcatac/c2mc.tsv")

## Move scATACs from fAMCs 602,603 (glia) to other fAMCs that AMCs 16_16, 16_20 flow to
flow_res_path <- file.path(wd, "output/mcatac/pl_cort_flow_mat.tsv")
flow_res_results <- readr::read_tsv(flow_res_path)
frr_mat <- as.matrix(dplyr::select(flow_res_results, -rowname))
rownames(frr_mat) <- flow_res_results$rowname
glia_ct <- c('Astrocytes', 'OPCs')
mc_glia <- which(mcmd$cell_type %in% glia_ct)
nnz_frcs <- apply(frr_mat[c('16_16', '16_20'),], 1, function(x) {sum_glia <- sum(x[mc_glia]); x[mc_glia] <- 0; nnz_mc <- which(x > 0); x_nnz_frc <- x[nnz_mc]/sum(x[nnz_mc]); return(setNames(x_nnz_frc, nnz_mc))})
nnz_cumsum <- lapply(nnz_frcs, function(x) cumsum(sort(x)))
flow_by_ct <- t(tgs_matrix_tapply(frr_mat, mcmd$cell_type, sum))
flow_by_ct_norm <- flow_by_ct/rowSums(flow_by_ct)

famc_new <- do.call('rbind', lapply(c(602, 603), function(x) {
    scah <- c2mc$cell[c2mc$metacell == x]; 
    amch <- rownames(frr_mat)[which(frr_mat[,x] > 0)];
    sca_rn <- runif(n = length(scah))
    famc_new <- as.numeric(sapply(seq_along(scah), function(i) mc_new <- names(nnz_cumsum[[amch]][nnz_cumsum[[amch]] >= sca_rn[[i]]])[[1]]))
    # print(cbind(sca_rn, famc_new))
    return(tibble::enframe(setNames(famc_new, scah), name = 'cell', value = 'metacell'))
}))

c2mc_new <- c2mc
c2mc_new$metacell[match(famc_new$cell, c2mc_new$cell)] <- famc_new$metacell


mcc <- scc_project_on_mc(sc_counts=scc, cell_to_metacell = dplyr::rename(c2mc_new, cell_id = cell))
feat_peak$peak_name <- peak_names(feat_peak, tad_based = F)
mcc_write(mcc, './output/mcatac/mmcortex_mcc_feat_peak_no_602_603', overwrite = T)

### Don't remove
####################


mca <- mcc_to_mcatac(mc_counts=mcc, peaks=feat_peak, metadata=mcmd)
saveRDS(mca, './output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds')
# mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks.rds')
# a_legc <- log2(1e-5 + mca@egc)

# km_a_legc <- tglkmeans::TGL_kmeans(as.matrix(a_legc), k = 80, seed = 1337)
# saveRDS(km_a_legc, './output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
# # km_a_legc <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')

# new_clust_vec <- km_a_legc$cluster
# seq_coords <- mca@peaks[,c('chrom', 'start', 'end')]

# coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)

# seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# prego::set_parallel(10)
# aa <- getOption('prego.parallel.nc')
# print(aa)
# options(error = recover)

# rm(list = ls()[!(ls() %in% c('seqs_all', 'new_clust_vec', "ALLGENOME", "GINTERVS", "GROOT", "GTRACKS", "GWD"))])
# print(ls())

# res <- prego::regress_pwm.clusters(sequences=seqs_all, clusters = new_clust_vec, 
#                                         use_sge = T, 
#                                         match_with_db=T, use_sample = T, 
#                                         parallel=T)
# print(warnings())
# saveRDS(res, './output/sequence_modeling/mmcortex_mcatac_rna_match_feat_peak_res.rds')