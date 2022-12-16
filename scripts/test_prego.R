library(metacell)
library(pheatmap)
devtools::load_all("~/src/prego/")
devtools::load_all("~/src/mcATAC/")
gset_genome('mm10')
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
nm <- "pl_cort_feat_peaks"

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))                         



####################
### Don't remove
## How to create feat_peak mat object
load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/multi_mmcortex.Rda",v=T)
load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_mmcortex.Rda",v=T)
load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_annot_mmcortex.Rda",v=T)

#all intervals
head(multi_model$atac_intervs)
#those in clusters annotated as variable
# feat_peak = multi_model$atac_intervs[acn$f_var_peak,]
feat_peak = multi_model$atac_intervs
scc <- scc_read('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads/')
mat_prom <- scdb_mat('pl_prom_cort')
scmat <- scc_extract(scc = scc, intervals = feat_peak)
scmd <- mat_prom@cell_metadata[colnames(scmat),]
mat_new <- scm_new_matrix(mat = scmat, cell_metadata=scmd, stat_type='umi')
scdb_add_mat(nm, mat_new)
### Don't remove
####################


mat_feat = scdb_mat(nm)


day_mat <- mcmd[,grep('^E', colnames(mcmd))]
colnames(day_mat) <- gsub('E', '', colnames(day_mat))

col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))

cs <- Matrix::colSums(mat_feat@mat)
mat_ds <- scm_downsamp(mat_feat@mat, quantile(cs, 0.1))

### 13:18 are days of samples
sc_res <- lapply(13:18, function(d) {
    cells_day <- colnames(mat_ds)[colnames(mat_ds) %in% rownames(mat_feat@cell_metadata)[mat_feat@cell_metadata$day == d]]
    
    mc_ord_day <- cust_mc_ord_st[cust_mc_ord_st %in% which(day_mat[,as.character(d)] >= 5)]
    
    sc_cor <- tgs_cor(log2(1+as.matrix(mat_ds[,cells_day])))
    return(list(sc_cor = sc_cor
               ))
})

sc_cor_kms <- parallel::mclapply(sc_res, function(x) tglkmeans::TGL_kmeans(x$sc_cor, k = 30, seed = 1337), mc.cores = 6)

mcl_days <- mapply(FUN = function(d, km) {
    cells_day <- colnames(mat_ds)[colnames(mat_ds) %in% rownames(mat_feat@cell_metadata)[mat_feat@cell_metadata$day == d]]
    mcls_day <- t(tgs_matrix_tapply(as.matrix(mat_ds[,cells_day]), km$cluster, sum))
    mcls_day <- as(mcls_day, 'dgCMatrix')
    return(mcls_day)
},  d= 13:18, km=  sc_cor_kms)

saveRDS(sc_cor_kms, '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/microcluster_assignment.RDS')

mcl_all <- do.call('cbind', mcl_days)
colnames(mcl_all) <- as.character(sapply(13:18, function(x) stringr::str_c(x, 1:30, sep = '_')))
old_clust_vec <- atac_clsts$km$cluster[atac_clsts$km$cluster %in% atac_clsts$vclst_nms]

a_legc <- log2(1e-5 + t(t(mcl_all/colSums(mcl_all))))
colnames(a_legc) <- colnames(mcl_all) 

aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
mcl_all <- aaa$mcl_all
old_clust_vec <- aaa$old_clust_vec
km_a_legc <- tglkmeans::TGL_kmeans(as.matrix(a_legc), k = 80, seed = 1337)
saveRDS(list(mcl_all = mcl_all, sc_res = sc_res, mat_ds = mat_ds, old_clust_vec = old_clust_vec, km_a_legc = km_a_legc),
             "/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
new_clust_vec <- km_a_legc$cluster
# clust_ind_mat <- sapply(sort(unique(old_clust_vec)), function(cl) ifelse(old_clust_vec == cl, 1, 0))
clust_ind_mat <- sapply(sort(unique(new_clust_vec)), function(cl) ifelse(new_clust_vec == cl, 1, 0))
rownames(clust_ind_mat) <- rownames(a_legc)

seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(a_legc))
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)
# seqs_all <- unlist(sapply(1:ceiling(nrow(seq_coords)/1000), function(i) gseq.extract(seq_coords[(1+(i-1)*1000):(min(i*1000, nrow(seq_coords))),])))
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# clust_ind_ls <- lapply(sort(unique(km_a_legc$cluster)), function(cl) ifelse(km_a_legc$cluster == cl, 1, 0))

# res <- regress_pwm.two_phase(seqs_all, clust_ind_mat[, 1], two_phase_sample_frac = c(0.1, 1), first_phase_metric = "ks")
resa <- sapply(1:ncol(clust_ind_mat), function(i) regress_pwm.sample(seqs_all, clust_ind_mat[, i], improve_epsilon = 0.01, match_with_db = TRUE))

res <- regress_pwm.clusters(sequences = seqs_all, clusters = new_clust_vec, use_sample = TRUE, match_with_db = TRUE, sample_frac = NULL, sample_ratio = 1, final_metric = "ks", parallel = TRUE)

# dir.create('./output/mcatac/figs')
dir.create('./output/mcatac/figs/pwm_regress_122k_seqs')
for (i in 1:length(resa)) {
    # png(glue::glue('./output/mcatac/figs/pwm_regress_test/cluster_{i}.png'), w = 800, h = 800)
    plot_regression_qc(resa[[i]])
    ggsave(glue::glue('./output/mcatac/figs/pwm_regress_122k_seqs/cluster_{i}.png'))
    # dev.off()
}

saveRDS(resa, "/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_pwm_regress_122k_seqs.rds")
# keys_to_take <- unique(unlist(sapply(c('POU3', 'NEURO', 'SOX2'), function(g) grep(g, pssms[[1]][[1]]$track))))
# names(keys_to_take) <- unique(unlist(sapply(c('POU3', 'NEURO', 'SOX2'), function(g) grep(g, pssms[[1]][[1]]$track,v=T))))
# pssms_to_take <- lapply(keys_to_take, function(k) pssms[[2]][[1]][pssms[[2]][[1]]$key == k,])

# samp_inds <- sample(1:length(clust_ind_ls), 10)
# pwms_neuro <- parallel::mclapply(samp_inds, function(i) prego::regress_pwm(sequences = seqs_all, 
#                                                                             response = clust_ind_ls[[i]], 
#                                                                             motif = pssms_to_take$NEUROG2), 
#                                                                             mc.cores = 10)
# pwms_null <- parallel::mclapply(clust_ind_ls, function(rs) prego::regress_pwm(sequences = seqs_all, response = rs), mc.cores = 16)