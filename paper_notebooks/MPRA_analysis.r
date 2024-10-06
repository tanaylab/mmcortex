library(tgstat)
library(dplyr)
library(matrixStats)
library(pheatmap)
library(MPRAnalyze)
library(data.table)
library(cowplot)
library(patchwork)
library(misha)
library(misha.ext)
library(xgboost)
library(metacell)

wd <- "/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex"
setwd(wd)


gsetroot('/home/aviezerl/mm10')

mca <- readRDS(file.path(wd, 'output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds'))
a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))
mcp <- dplyr::select(mca@peaks, chrom, start, end, peak_name)
scdb_init(file.path(wd, 'scdb'), f=T)
mc <- scdb_mc('pl_cort')
mat <- scdb_mat('pl_cort')

print('dim of mat@mat')
print(dim(mat@mat))
tss <- gintervals.load('intervs.global.tss')
tads <- gintervals.load('intervs.global.tad_names')

orig_lib <- readr::read_tsv(file.path(wd, 'data/st_and_temporal_and_e14_shadow_enh_11-1-22.tsv'))

ol_crs_df <- as.data.frame(cbind(do.call('rbind', stringr::str_split(orig_lib$rowname, '\\.|-')), orig_lib$rowname))
colnames(ol_crs_df) <- c('cell_type', 'pattern', 'number', 'crs_name')
head(ol_crs_df)

load(file.path(wd, 'output/mcatac/ct_peaks.rda'))

nei_orig_ct <- gintervals.neighbors(dplyr::relocate(orig_lib, rowname, .after = end), ct_peaks, mindist = 0, maxdist = 0, maxneighbors = 1e+3)

mcmd <- readr::read_tsv(file.path(wd, 'BonevCollab/mcmd_pl_cort.tsv'))
mcmd <- mcmd[!(as.numeric(mcmd$metacell) %in% 602:603),]

col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
col_key

# Import data
mpra_annot <- read.table(file.path(wd, 'MPRA_data/st_and_temporal_and_e14_shadow_enh_11-1-22.tsv'),header=T)
res <- read.table(file.path(wd, 'MPRA_data/MPRA_dev_090223.tsv'),header=T)

# Correct p-value (BH)
res <- res[,c(grep('names',colnames(res)),grep('statistic',colnames(res)),grep('mad.score',colnames(res)),grep('pval.mad',colnames(res)),grep('enh_type',colnames(res)))]
res[,grep('pval.mad',colnames(res))] <- apply(res[,grep('pval.mad',colnames(res))],2,function(x){
  return(p.adjust(x,method='BH'))
})

# df <- res[res$enh_type=='NSC',]
df <- res
df$temp <- gsub('seqs_','',as.data.frame(stringr::str_split(paste0(df$names), pattern ='\\.' , n = 3, simplify = TRUE))[,2])
df <- df[rowMins(as.matrix(df[,grep('NSC|IPC', grep('pval.mad',colnames(df), v=T),v=T)]), na.rm=T)<=0.1,]           #Select enhancers which are significant in at least one timepoint 

df_nsc <- df[,c('names', grep('NSC', grep('mad.score', colnames(df), v=T), v=T))]
rownames(df_nsc) <- df_nsc$names
df_ipc <- df[,c('names', grep('IPC', grep('mad.score', colnames(df), v=T), v=T))]
rownames(df_ipc) <- df_ipc$names
pltmt_nsc <- df_nsc[apply(df_nsc, 1, function(x) all(!is.na(x))),grep('_E1\\d', colnames(df_nsc))]
rn <- rownames(pltmt_nsc)
pltmt_nsc <- apply(pltmt_nsc, 2, function(x) log2(x - min(x) + 5e-1))              
rownames(pltmt_nsc) <- rn 
                   
                   
pltmt_ipc <- df_ipc[apply(df_ipc, 1, function(x) all(!is.na(x))),grep('_E1\\d', colnames(df_ipc))]
rn <- rownames(pltmt_ipc)
pltmt_ipc <- apply(pltmt_ipc, 2, function(x) log2(x - min(x) + 5e-1))              
rownames(pltmt_ipc) <- rn

orig_lib_f <- orig_lib[grep('chr', orig_lib$rowname, inv = T),]
e14_lib <- as.data.frame(stringr::str_split(gsub('E14_', '', grep('chr', orig_lib$rowname, v=T)), '-', simplify = T))
colnames(e14_lib) <- c('chrom', 'start', 'end')
e14_lib[,2:3] <- apply(e14_lib[,2:3], 2, as.numeric)


rownames(e14_lib) <- paste0('E14_', grep('chr', orig_lib$rowname, v=T))
shadow_lib <- as.data.frame(orig_lib[grep('chr', orig_lib$rowname),c('chrom', 'start', 'end')])
rownames(shadow_lib) <- paste0('shadow_', grep('chr', orig_lib$rowname, v=T))

orig_lib_a <- dplyr::bind_rows(orig_lib_f, tibble::rownames_to_column(e14_lib), tibble::rownames_to_column(shadow_lib))
# save(orig_lib_a, file = file.path(wd, 'output/MPRA/orig_lib_w_e14_and_shadow_coords.rda'))
seqs_tdf <- as.data.frame(stringr::str_split(orig_lib_a$rowname, '_|\\.', simplify = T))
seqs_tdf[seqs_tdf[,3] == '',3] <- NA
rownames(seqs_tdf) <- orig_lib_a$rowname
colnames(seqs_tdf) <- c('ct', 'seqs', 'temp', 'dummy1', 'dummy2')


monuc_track <- gtrack.ls('seq.G_or_C')

monuc_seq_coords <- gextract(monuc_track, intervals = as.data.frame(orig_lib_a[,c('chrom', 'start', 'end')]), colnames = gsub('seq\\.', '', monuc_track), iterator = 1)
monuc_mat <- subset(monuc_seq_coords, select = -c(chrom, start, end, intervalID))
monuc_per_peak <- as.data.frame(tgs_matrix_tapply(t(monuc_mat), monuc_seq_coords$intervalID, function(x) length(which(x == 1)))/266)
rownames(monuc_per_peak) <- orig_lib_a$rowname
monuc_per_peak$A_or_T <- 1 - monuc_per_peak$G_or_C
# save(monuc_per_peak, file = file.path(wd, 'output/MPRA/monuc_per_peak_orig_lib_w_e14_and_shadow.rda'))
# load(file = file.path(wd, 'output/MPRA/monuc_per_peak_orig_lib_w_e14_and_shadow.rda'))

dinuc_tracks <- gtrack.ls('seq\\.[ACGT]{2}$')

dinucs_seq_coords <- gextract(dinuc_tracks, intervals = as.data.frame(orig_lib_a[,c('chrom', 'start', 'end')]), colnames = gsub('seq\\.', '', dinuc_tracks), iterator = 1)

# saveRDS(object = dinucs_seq_coords, file = file.path(wd, 'output/MPRA/dinucs_orig_lib.rds'))
dinucs_mat <- subset(dinucs_seq_coords, select = -c(chrom, start, end, intervalID))

dinucs_per_peak <- tgs_matrix_tapply(t(dinucs_mat), dinucs_seq_coords$intervalID, function(x) length(which(x == 1)))/266
rownames(dinucs_per_peak) <- orig_lib_a$rowname
# save(dinucs_per_peak, file = file.path(wd, 'output/MPRA/dinucs_per_peak_orig_lib_w_e14_and_shadow.rda'))
# load(file.path(wd, 'output/MPRA/dinucs_per_peak_orig_lib_w_e14_and_shadow.rda'))

source(file.path(wd, 'scripts/util.r'))

load(file.path(wd, 'output/sequence_modeling/xgb_cv_res.rda'))
preds_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['contrib_meth']]))
colnames(preds_all)

agg_id <- readr::read_csv(file.path(wd, 'scatac_data//aggregation_id.csv'))
mca@cell_to_metacell$agg_id <- as.numeric(unlist(purrr::map(stringr::str_split(mca@cell_to_metacell$cell_id, '-'), 2)))
agg_id_day <- data.frame(cbind(1:14, rep(12:18, each = 2)))
colnames(agg_id_day) <- c('agg_id', 'day')
mca@cell_to_metacell$day <- agg_id_day$day[mca@cell_to_metacell$agg_id]
atac_mc_day <- table(mca@cell_to_metacell$metacell, mca@cell_to_metacell$day)
atac_mc_day_norm <- atac_mc_day/rowSums(atac_mc_day)
nsc_mcs <- which(mcmd$cell_type == 'NSC')
egc_by_day <- mca@egc[,nsc_mcs] %*% atac_mc_day_norm[nsc_mcs,]
egc_by_day_n <- t(t(egc_by_day)/colSums(egc_by_day))
colnames(egc_by_day_n) <- paste0('E', colnames(egc_by_day_n))
a_legc_by_day_n <- log2(1e-5 + egc_by_day_n)

load(file.path(wd, 'output/sequence_modeling/enhflow_model_n_motif=16_w_clustering.rda'))

rename_prego_motifs <- function(model, tfs_in) {
    amd <- prego::all_motif_datasets()
    amou <- unique(amd[,c('motif', 'motif_orig')])
    motif_clustering <- model@params$distilled_features
    prego_motifs_renamed <- sapply(names(model@motif_models), function(cni) {
        candidates <- motif_clustering$feat[motif_clustering$clust == cni]
        candidates <- setNames(candidates, unlist(purrr::map(stringr::str_split(amou$motif_orig[match(candidates, amou$motif)], '_'), 1)))
        candidates_all <- candidates[tolower(names(candidates)) %in% tolower(tfs_in)]
        # candidates_feats <- candidates[tolower(names(candidates)) %in% tolower(intersect(tfs_in, feats))]
        # if (length(candidates_feats) > 0) {
        #     return(candidates_feats[[1]])
        # } else
        if (length(candidates_all) > 0) {
            return(candidates_all[[1]])
        } else {
            return(candidates[[1]])
        }
    })
    return(prego_motifs_renamed) 
}
tfs_in <- readLines(file.path(wd, 'output/metacell_model/tfs_in.txt'))
prego_motifs <- do.call('rbind', lapply(seq_along(tm_w_add_feat@motif_models), function(i) dplyr::mutate(tm_w_add_feat@motif_models[[i]]$pssm, motif = names(tm_w_add_feat@motif_models)[[i]])))
orig_lib_prego_ie <- prego::gextract_pwm(intervals = orig_lib_a[,c('chrom', 'start', 'end', 'rowname')], motifs = unique(prego_motifs$motif), dataset = prego_motifs)
orig_lib_prego_ie_mat <- subset(orig_lib_prego_ie, select = -c(chrom, start, end, rowname))
rownames(orig_lib_prego_ie_mat) <- orig_lib_prego_ie$rowname
prego_motifs_renamed <- rename_prego_motifs(tm_w_add_feat, tfs_in = tfs_in)
colnames(orig_lib_prego_ie_mat) <- prego_motifs_renamed[colnames(orig_lib_prego_ie_mat)]

dtest_ole_mtt_meth_0_atac_1 <- cbind(orig_lib_prego_ie_mat, 'NSC_ATAC' = mean(quantile(rowMeans(a_legc_by_day_n[rownames(preds_all),]), c(0,1))),'methylation' = 0, 
                                        dinucs_per_peak[rownames(orig_lib_prego_ie_mat),intersect(colnames(dinucs_per_peak), colnames(preds_all))], 
                                        'GC content' = 0.05, 'prox_ATAC' = 10, 'prox_RNA' = 10)


pred01 <- setNames(predict(xgb_cv_res[[1]]$bstDMatrix_meth, xgb.DMatrix(as.matrix(dtest_ole_mtt_meth_0_atac_1)), predcontrib=FALSE, predinteraction = FALSE), rownames(dtest_ole_mtt_meth_0_atac_1))

inds_active_mpra <- multintersect(union(rownames(pltmt_nsc[rowMeans(pltmt_nsc) > quantile(rowMeans(pltmt_nsc), 0.75),]),
                             rownames(pltmt_ipc[rowMeans(pltmt_ipc) > quantile(rowMeans(pltmt_ipc), 0.75),])), rownames(pltmt_nsc), rownames(pltmt_ipc))

png(file.path(wd, 'output/MPRA/figs/test/MPRA_vs_xgb_pred_inds_active_by_day.png'), h = 400, w = 2000, res = 50)
par(mfrow = c(1,10), cex.axis = 3, cex.lab = 4, cex.main = 4, las = 2)
bins_en <- seq(-2, 2, l = 17)
par(col.main = col_key[['NSC']])
names_nc <- inds_active_mpra
vvv <- sapply(colnames(pltmt_nsc), function(cni) {
    if (cni == colnames(pltmt_nsc)[[1]]) {par(mar = c(12,8,3,3))}
    else {par(mar = c(12,2,3,3))}
    boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_nsc[names_nc,cni], nm = gsub('_mad.score', '', cni), bins = bins_en, 
                    ylab = '', ylim = c(-0.2,4.95), xlim = c(3.5,13), text_y_factor = 1., text_cex = 2, xaxt = 's', show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    # bin_and_ks(xvec = pred01[names_nc], yvec = rowMeans(pltmt_nsc[names_nc,]), bins_ks = bins_ks, ks.alternative = 'less')
    title(ylab = 'NSC MPRA', line = 4)
})
# par(mar = c(12,8,1,3))
par(col.main = col_key[['IPC']])
vvv <- sapply(colnames(pltmt_ipc), function(cni) {
    if (cni == colnames(pltmt_ipc)[[1]]) {par(mar = c(12,8,3,3))}
    else {par(mar = c(12,2,3,3))}
    boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_ipc[names_nc,cni], nm = gsub('_mad.score', '', cni), bins = bins_en, 
        ylab = '', ylim = c(0.04,5.95), xlim = c(3.5,13), text_y_factor = 1., text_cex = 2, show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    # bin_and_ks(xvec = pred01[names_nc], yvec = rowMeans(pltmt_ipc[names_nc,]), bins_ks = bins_ks, ks.alternative = 'greater')
    title(ylab = 'IPC MPRA', line = 4)
    title(xlab = '<- NSC  xgboost prediction  IPC ->', line = 17)
})
dev.off()

png(file.path(wd, 'output/MPRA/figs/test/astrocyte_specific_crs_vs_all_active_mpra_inds_ecdf.png'), h = 400, w = 400)

par(cex.lab = 2, cex.main = 1.5, cex.axis = 1.5, mar = c(5,5,3,1))
ksx <- apply(pltmt_nsc[intersect(rownames(pltmt_nsc), unique(nei_orig_ct$rowname[nei_orig_ct$type == 'astro_peak'])),], 1, function(x) sum((x - min(x))*(1:ncol(pltmt_nsc)))/sum((x - min(x))))
ksy <- apply(pltmt_nsc[setdiff(inds_active_mpra,unique(nei_orig_ct$rowname[nei_orig_ct$type == 'astro_peak'])),], 1, function(x) sum((x - min(x))*(1:ncol(pltmt_nsc)))/sum((x - min(x))))
plot(ecdf(ksx), do.points = F, col = col_key[['Astrocytes']], add = F, 
                xlab = 'Center of mass', ylab = 'ECDF', main = 'Astrocyte-specific CRSs are maximally\nactivated in later NSC MPRA time points', lwd = 2)
plot(ecdf(ksy), do.points =F, add = T, lwd = 2) 
legend('topleft', legend = c('Astrocyte-specific CRSs',paste0('n = ', length(ksx)), 'All active CRSs',paste0('n = ', length(ksy))), col = c(col_key[['Astrocytes']], NA, 'black', NA), lwd = 2)
dev.off()

ct_tracks <- c('mmcortex_nsc_marginal', 'mmcortex_ipc_marginal', 'mmcortex_cpn_l23_marginal','mmcortex_cpn_l56_marginal','mmcortex_cthpn_marginal','mmcortex_scpn_marginal')
ol_marg <- gextract(ct_tracks,intervals =  orig_lib_a[,c('chrom', 'start', 'end')], iterator = orig_lib_a[,c('chrom', 'start', 'end')])
mcp_trim <- dplyr::mutate(mcp[,c('chrom', 'start', 'end')],start = start + 17, end = end - 17)
ct_tracks_on_manifold <- gextract(ct_tracks,intervals =  mcp_trim, iterator = mcp_trim)
ol_marg$rowname <- orig_lib_a$rowname[ol_marg$intervalID]
ol_marg[,ct_tracks] <- apply(ol_marg[,ct_tracks], 2, function(x) {y <- x; y[is.na(x)] <- 0; return(y)})
ol_marg$log_nsc <- log2(1e-5 + ol_marg$mmcortex_nsc_marginal/sum(ol_marg$mmcortex_nsc_marginal))
ol_marg$log_ipc <- log2(1e-5 + ol_marg$mmcortex_ipc_marginal/sum(ol_marg$mmcortex_ipc_marginal))
ol_marg$log_scpn <- log2(1e-5 + ol_marg$mmcortex_scpn_marginal/sum(ol_marg$mmcortex_scpn_marginal))
ol_marg$log_cthpn <- log2(1e-5 + ol_marg$mmcortex_cthpn_marginal/sum(ol_marg$mmcortex_cthpn_marginal))
ol_marg$log_cpn_l23 <- log2(1e-5 + ol_marg$mmcortex_cpn_l23_marginal/sum(ol_marg$mmcortex_cpn_l23_marginal))
ol_marg$log_cpn_l56 <- log2(1e-5 + ol_marg$mmcortex_cpn_l56_marginal/sum(ol_marg$mmcortex_cpn_l56_marginal))

THRSH_A <- -12
THRSH_B <- -14.5
new_class <- setNames(rep(NA, nrow(ol_marg)), ol_marg$rowname)

new_class[abs(ol_marg$log_nsc - ol_marg$log_ipc) < 1 & (ol_marg$log_ipc >= (-ol_marg$log_nsc + 2*THRSH_A))] <- 'A'
new_class[abs(ol_marg$log_nsc - ol_marg$log_ipc) < 1 & (ol_marg$log_ipc < (-ol_marg$log_nsc + 2*THRSH_A))] <- 'B'
new_class[abs(ol_marg$log_nsc - ol_marg$log_ipc) < 1 & (ol_marg$log_ipc < (-ol_marg$log_nsc + 2*THRSH_B))] <- 'C'
new_class[abs(ol_marg$log_nsc - ol_marg$log_ipc) >= 1 & ol_marg$log_nsc > ol_marg$log_ipc] <- 'NSC'
new_class[abs(ol_marg$log_nsc - ol_marg$log_ipc) >= 1 & ol_marg$log_nsc < ol_marg$log_ipc] <- 'IPC'

col_key[c('A', 'B', 'C')] <- c('orange', 'darkgreen', 'magenta')
ol_marg[match(names(new_class), ol_marg$rowname),'new_class'] <- new_class

ol_marg$pred01 <- pred01[ol_marg$rowname]

ol_marg$pred01_cut <- cut(ol_marg$pred01, breaks = quantile(ol_marg$pred01, c(0,0.2,0.8,1), na.rm = T), labels = c('pred NSC', 'pred und.', 'pred IPC'))

rm_nsc <- rowMeans(pltmt_nsc[,])
rm_ipc <- rowMeans(pltmt_ipc[,])

ol_marg$rm_nsc <- rm_nsc[ol_marg$rowname]

ol_marg$rm_ipc <- rm_ipc[ol_marg$rowname]

ol_marg$cut_rm_nsc <- cut(ol_marg$rm_nsc, breaks = quantile(ol_marg$rm_nsc, na.rm = T, seq(0,1,l=4)), labels = c('NSC MPRA lo', 'NSC MPRA mid', 'NSC MPRA hi'))
ol_marg$cut_rm_ipc <- cut(ol_marg$rm_ipc, breaks = quantile(ol_marg$rm_ipc, na.rm = T, seq(0,1,l=4)), labels = c('IPC MPRA lo', 'IPC MPRA mid', 'IPC MPRA hi'))

neu_rn <- ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc <= -14 & ol_marg$pred01 >= quantile(ol_marg$pred01, 0.7, na.rm = T) & ol_marg$cut_rm_ipc %in% paste0('IPC MPRA ', c('hi')))]
ipc_pc_rn <- ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc >= -14 & ol_marg$pred01 >= quantile(ol_marg$pred01, 0.7, na.rm = T) & ol_marg$cut_rm_ipc %in% paste0('IPC MPRA ', c('hi')))]
anti_neu_rn <- ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc <= -14 & ol_marg$pred01 <= quantile(ol_marg$pred01, 0.3, na.rm = T) & ol_marg$cut_rm_nsc %in% paste0('NSC MPRA ', c('hi')))]
nsc_pc_rn <- ol_marg$rowname[which(ol_marg$log_nsc >= -14 & ol_marg$log_ipc <= -14 & ol_marg$pred01 <= quantile(ol_marg$pred01, 0.3, na.rm = T) & ol_marg$cut_rm_nsc %in% paste0('NSC MPRA ', c('hi')))]

hi_nsc_lo_ipc_rn <- ol_marg$rowname[which(ol_marg$log_ipc <= -14 & ol_marg$log_nsc >= -14 & ol_marg$pred01 <= quantile(ol_marg$pred01, 0.7, na.rm = T) & ol_marg$pred01 >= quantile(ol_marg$pred01, 0.3, na.rm = T))]
hi_ipc_lo_nsc_rn <- ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc >= -14 & ol_marg$pred01 >= quantile(ol_marg$pred01, 0.3, na.rm = T) & ol_marg$pred01 <= quantile(ol_marg$pred01, 0.7, na.rm = T))]
nc_neu_rn <- setdiff(ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc <= -14)], union(neu_rn, anti_neu_rn))
nc_hi_ipc_lo_nsc_rn <- setdiff(ol_marg$rowname[which(ol_marg$log_nsc <= -14 & ol_marg$log_ipc >= -14)], union(hi_ipc_lo_nsc_rn, hi_nsc_lo_ipc_rn))
nc_hi_nsc_lo_ipc_rn <- setdiff(ol_marg$rowname[which(ol_marg$log_ipc <= -14 & ol_marg$log_nsc >= -14)], union(hi_ipc_lo_nsc_rn, hi_nsc_lo_ipc_rn))

ctxt_rn_intervs <- ol_marg[,c('chrom','start', 'end', 'rowname')]
ctxt_rn_intervs$peak_name <- mcATAC::peak_names(ctxt_rn_intervs, tad_based = F)

# nei_orig_peaks <- gintervals.neighbors(orig_lib_a[,c('chrom', 'start', 'end', 'rowname')], 
#                                        as.data.frame(mcp),  
#                                        mindist = -30e+5, maxdist = 30e+5, maxneighbors = 1e+6)
# nei_orig_peaks <- nei_orig_peaks[abs(nei_orig_peaks$dist) >= 1e+3,]

# ipc_marginal_peaks_umis <- Matrix::rowSums(mca@mat[nei_orig_peaks$peak_name,which(mcmd$cell_type == 'IPC')])
# nsc_marginal_peaks_umis <- Matrix::rowSums(mca@mat[nei_orig_peaks$peak_name,which(mcmd$cell_type == 'NSC')])
# nei_orig_peaks$is_anti_neu <- ifelse(nei_orig_peaks$rowname %in% anti_neu_rn, TRUE, FALSE)
# nei_orig_peaks$is_ctxt_ipc <- ifelse(nei_orig_peaks$rowname %in% hi_ipc_lo_nsc_rn, TRUE, FALSE)
# nei_orig_peaks$is_ctxt_nsc <- ifelse(nei_orig_peaks$rowname %in% hi_nsc_lo_ipc_rn, TRUE, FALSE)
# nei_orig_peaks$is_neu <- ifelse(nei_orig_peaks$rowname %in% neu_rn, TRUE, FALSE)
# nei_orig_peaks$is_ipc_pc <- ifelse(nei_orig_peaks$rowname %in% ipc_pc_rn, TRUE, FALSE)
# nei_orig_peaks$is_nsc_pc <- ifelse(nei_orig_peaks$rowname %in% nsc_pc_rn, TRUE, FALSE)

# nei_orig_peaks$ipc_atac_umis <- ipc_marginal_peaks_umis[nei_orig_peaks$peak_name]
# nei_orig_peaks$nsc_atac_umis <- nsc_marginal_peaks_umis[nei_orig_peaks$peak_name]

# get_neighbor_mat <- function(order_metric_name, stat_type, nei_peak_df, crs_rownames, atac_umi_colname = 'ipc_atac_umis', bin_vec = seq(0, 30e+5, 1e+4)) {
#     mat_inds <- which(nei_orig_peaks$rowname %in% crs_rownames)
#     sks <- unique(nei_orig_peaks$rowname[mat_inds])
#     metric <- sapply(seq_along(sks), function(i) {
#         sk_inds <- which(nei_peak_df$rowname == sks[[i]])
#         if (order_metric_name == 'com_abs') {
#             y <- sum(abs(nei_peak_df$dist[sk_inds]) * nei_peak_df[sk_inds, atac_umi_colname])/sum(nei_peak_df[sk_inds, atac_umi_colname])
#         } else if (order_metric_name == 'com') {
#             y <- sum(nei_orig_peaks$dist[sk_inds] * nei_peak_df[sk_inds, atac_umi_colname])/sum(nei_peak_df[sk_inds, atac_umi_colname])
#         } else if (order_metric_name == 'entropy') {
#             dist_bins <- cut(nei_orig_peaks$dist[sk_inds], breaks = bin_vec)
#             pv <- tapply(nei_peak_df[sk_inds, atac_umi_colname], dist_bins, sum, na.rm = T)
#             pv[is.na(pv)] <- 0
#             y <- entropy::entropy(pv)
#         } else if (order_metric_name == 'sum') {
#             y <- sum(nei_peak_df[sk_inds, atac_umi_colname])
#         }
#         else {stop('order_metric_name must be either: "com", "com_abs", "entropy" or "sum"')}
#         return(y)
#     })
#     metric_ord <- order(metric)
#     neighbor_mat <- t(sapply(seq_along(metric_ord), function(i) {
#         sk_inds <- which(nei_peak_df$rowname == sks[[metric_ord[[i]]]])
#         dist_bins <- cut(nei_orig_peaks$dist[sk_inds], breaks = bin_vec)
#         if (stat_type == 'sum_umis') {
#             pv <- tapply(nei_peak_df[sk_inds, atac_umi_colname], dist_bins, sum, na.rm = T)   
#         } else if (stat_type == 'num_neighbors') {
#             pv <- table(dist_bins)
#         }
#         pv[is.na(pv)] <- 0
#         return(pv)
#     }))
#     rownames(neighbor_mat) <- sks[metric_ord]
#     colnames(neighbor_mat) <- head(bin_vec, -1)
#     return(neighbor_mat)
# }
# neu_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = neu_rn)
# print('Done mat 1')
# ipc_pc_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = ipc_pc_rn)
# print('Done mat 2')
# ipc_ctxt_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = hi_ipc_lo_nsc_rn)
# print('Done mat 3')
# neu_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = neu_rn)
# print('Done mat 4')
# ipc_pc_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = ipc_pc_rn)
# print('Done mat 5')
# ipc_ctxt_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = hi_ipc_lo_nsc_rn)
# print('Done mat 6')
# anti_neu_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 7')
# nsc_pc_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 8')
# nsc_ctxt_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 1')
# anti_neu_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 10')
# nsc_pc_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 11')
# nsc_ctxt_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')
# print('Done mat 12')


# png(file.path(wd, 'output/MPRA/figs/test/sum_nei_and_umis_over_50kbp_per_ipc_crs_categories.png'), h = 400, w = 800)
# par(mfrow = c(1,2))
# par(cex.lab = 2, mar = c(5,7,2,1))
# K <- 25
# LWD <- 3
# bg_col <- 'whitesmoke'
# x <- as.numeric(colnames(neu_num_neighbor_mat))
# plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(0.8,2e+0), ylab = glue::glue('Mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T), col = 'red', lwd = LWD)
# lines(x, zoo::rollmean(colMeans(ipc_pc_num_neighbor_mat), k = K, na.pad = T), col = col_key[['IPC']], lwd = LWD)
# lines(x, zoo::rollmean(colMeans(ipc_ctxt_num_neighbor_mat), k = K, na.pad = T), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# par(cex.lab = 2, mar = c(5,7,2,1))
# plot(x, zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(8e+2,2.3e+3), ylab = glue::glue('Mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T), col = 'red', lwd = LWD)
# lines(x, zoo::rollmean(colMeans(ipc_pc_umi_neighbor_mat), k = K, na.pad = T), col = col_key[['IPC']], lwd = LWD)
# lines(x, zoo::rollmean(colMeans(ipc_ctxt_umi_neighbor_mat), k = K, na.pad = T), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# dev.off()

# png(file.path(wd, 'output/MPRA/figs/test/sum_nei_and_umis_over_50kbp_per_nsc_crs_categories.png'), h = 400, w = 800)
# par(mfrow = c(1,2))
# par(cex.lab = 2, mar = c(5,7,2,1))
# K <- 25
# LWD <- 3
# bg_col <- 'whitesmoke'
# x <- as.numeric(colnames(neu_num_neighbor_mat))
# plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(0.8,1.62e+0), ylab = glue::glue('Mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, zoo::rollmean(colMeans(anti_neu_num_neighbor_mat), k = K, na.pad = T), col = 'magenta', lwd = LWD)
# lines(x, zoo::rollmean(colMeans(nsc_pc_num_neighbor_mat), k = K, na.pad = T), col = col_key[['NSC']], lwd = LWD)
# lines(x, zoo::rollmean(colMeans(nsc_ctxt_num_neighbor_mat), k = K, na.pad = T), col = 'green', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# par(cex.lab = 2, mar = c(5,7,2,1))
# plot(x, zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(8e+2,1.82e+3), ylab = glue::glue('Mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, zoo::rollmean(colMeans(anti_neu_umi_neighbor_mat), k = K, na.pad = T), col = 'magenta', lwd = LWD)
# lines(x, zoo::rollmean(colMeans(nsc_pc_umi_neighbor_mat), k = K, na.pad = T), col = col_key[['NSC']], lwd = LWD)
# lines(x, zoo::rollmean(colMeans(nsc_ctxt_umi_neighbor_mat), k = K, na.pad = T), col = 'green', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# dev.off()


source(file.path(wd, 'scripts/util.r'))
nsc_pca <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]), restrict_to_tads = F,
                                       mc_sel = which(mcmd$cell_type == 'NSC'), mca = mca , mat_rna = mat , mc_rna = mc, tss = tss, tads = tads
                           )

ipc_pca <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]),restrict_to_tads = F,
                                       mc_sel = which(mcmd$cell_type == 'IPC'), mca = mca , mat_rna = mat , mc_rna = mc, tss = tss, tads = tads
                           )
nsc_pca$rowname <- ctxt_rn_intervs$rowname[match(nsc_pca$peak_name, ctxt_rn_intervs$peak_name)]
nsc_pca$new_class <- ol_marg$new_class[match(nsc_pca$rowname, ol_marg$rowname)]
ipc_pca$rowname <- ctxt_rn_intervs$rowname[match(ipc_pca$peak_name, ctxt_rn_intervs$peak_name)]
ipc_pca$new_class <- ol_marg$new_class[match(ipc_pca$rowname, ol_marg$rowname)]
nsc_pca$is_ctxt <- ifelse(nsc_pca$rowname %in% hi_ipc_lo_nsc_rn, TRUE, FALSE)
ipc_pca$is_ctxt <- ifelse(ipc_pca$rowname %in% hi_nsc_lo_ipc_rn, TRUE, FALSE)
nsc_pca$is_neu <- ifelse(nsc_pca$rowname %in% neu_rn, TRUE, FALSE)
ipc_pca$is_neu <- ifelse(ipc_pca$rowname %in% neu_rn, TRUE, FALSE)
nsc_pca$is_anti_neu <- ifelse(nsc_pca$rowname %in% anti_neu_rn, TRUE, FALSE)
ipc_pca$is_anti_neu <- ifelse(ipc_pca$rowname %in% anti_neu_rn, TRUE, FALSE)


png(file.path(wd, 'output/MPRA/figs/test/atac_mpra_prox_atac_umi_of_potentially_silenced_crs_groups.png'), h = 700, w = 1550)

####### IPC plot ######
######################


nri <- ol_marg$rowname %in% neu_rn
ord <- order(as.numeric(nri))
layout(mat = rbind(c(1,1,1,2,2,3,3,4,4,5,5),c(6,6,6,7,7,8,8,9,9,10,10)))
par(cex.lab = 2.5, cex.axis = 1.5, mar = c(5,5,1,1))

CEX <- 0.75
plot(ol_marg$log_nsc[ord], ol_marg$log_ipc[ord], pch = 16, col = ifelse(ol_marg$rowname[ord] %in% neu_rn, 'red', 'black'), xlab = 'NSC ATAC', ylab = 'IPC ATAC', cex = ifelse(ol_marg$rowname[ord] %in% neu_rn, 2*CEX, CEX))
points(ol_marg$log_nsc[match(hi_ipc_lo_nsc_rn, ol_marg$rowname)], ol_marg$log_ipc[match(hi_ipc_lo_nsc_rn, ol_marg$rowname)], pch = 16, col = 'orange3', cex = 2*CEX)
points(ol_marg$log_nsc[match(ipc_pc_rn, ol_marg$rowname)], ol_marg$log_ipc[match(ipc_pc_rn, ol_marg$rowname)], pch = 16, col = col_key[['IPC']], cex = 2*CEX)
legend('bottomright', col = c('red', col_key[['IPC']], 'orange3'), legend = paste0('n = ', c(length(neu_rn), length(ipc_pc_rn), length(hi_ipc_lo_nsc_rn))), pch = 15, cex = 2)


par(mar = c(4,5,1,1))
locs <- seq(1,3.5,0.5)
boxplot(ol_marg$pred01[ol_marg$rowname %in% neu_rn], at = locs[[1]], col = 'red', add = F, xlim = c(0.75,2.25), ylim = c(-1.2,.51), ylab = 'NSC <-    xgb pred    -> IPC')
boxplot(ol_marg$pred01[ol_marg$rowname %in% ipc_pc_rn], at = locs[[2]], col = col_key[['IPC']], add = T)
boxplot(ol_marg$pred01[ol_marg$rowname %in% hi_ipc_lo_nsc_rn], at = locs[[3]], col = 'orange3', add = T)

par(mar = c(5,5,1,1))
locs <- seq(1,2.5,0.5)
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% neu_rn], at = locs[[1]], col = 'red', add = F, xlim = c(0.75,2.25), ylim = c(0.9,2.7), ylab = 'MPRA score')
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% ipc_pc_rn], at = locs[[2]], col = col_key[['IPC']], add = T)
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% hi_ipc_lo_nsc_rn], at = locs[[3]], col = 'orange3', add = T)

locs <- seq(1,3,l=5)
boxplot(ipc_pca$atac_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.62e+4), ylab = '', col = 'darkgray')
boxplot(ipc_pca$atac_umi[ipc_pca$is_neu == T], at = locs[[2]], add = T, col = 'red')
boxplot(ipc_pca$atac_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], at = locs[[3]], add = T, col = 'orange3')
boxplot(ipc_pca$atac_umi[match(ipc_pc_rn, ipc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['IPC']])

grid(col = 'lightblue', lwd = 2, lty = 2)
title(ylab = 'Proximal ATAC UMI (50kbp)', line = 3)

locs <- seq(1,5,l=9)

boxplot(ipc_pca$rna_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.262e+5), ylab = '', col = 'darkgray', yaxt = 'n')
boxplot(ipc_pca$rna_umi[ipc_pca$is_neu == T], at = locs[[2]], add = T, col = 'red', yaxt = 'n')
boxplot(ipc_pca$rna_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], at = locs[[3]], add = T, col = 'orange3', yaxt = 'n')
# boxplot(ipc_pca$atac_umi[ipc_pca$is_ctxt == T], at = locs[[4]], add = T, col = 'green')
boxplot(ipc_pca$rna_umi[match(ipc_pc_rn, ipc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['IPC']], yaxt = 'n')
axis(2, at = seq(0,1.5e+5, 5e+4), labels = signif(seq(0,1.5e+5, 5e+4), 2))
title(ylab = 'Proximal RNA UMI (500kbp)', line = 3)


####### NSC plot ######
######################
CEX <- 0.75
nri <- ol_marg$rowname %in% anti_neu_rn
ord <- order(as.numeric(nri))
plot(ol_marg$log_nsc[ord], ol_marg$log_ipc[ord], pch = 16, col = ifelse(ol_marg$rowname[ord] %in% anti_neu_rn, 'magenta', 'black'), xlab = 'NSC ATAC', ylab = 'IPC ATAC', cex = ifelse(ol_marg$rowname[ord] %in% anti_neu_rn, 2*CEX, CEX))
points(ol_marg$log_nsc[match(hi_nsc_lo_ipc_rn, ol_marg$rowname)], ol_marg$log_ipc[match(hi_nsc_lo_ipc_rn, ol_marg$rowname)], pch = 16, col = 'green', cex = 2*CEX)
points(ol_marg$log_nsc[match(nsc_pc_rn, ol_marg$rowname)], ol_marg$log_ipc[match(nsc_pc_rn, ol_marg$rowname)], pch = 16, col = col_key[['NSC']], cex = 2*CEX)
legend('bottomright', col = c('magenta', col_key[['NSC']], 'green'), legend = paste0('n = ', c(length(anti_neu_rn), length(nsc_pc_rn), length(hi_nsc_lo_ipc_rn))), pch = 15, cex = 2)


par(mar = c(4,5,1,1))
locs <- seq(1,3,0.5)
boxplot(ol_marg$pred01[ol_marg$rowname %in% anti_neu_rn], at = locs[[1]], col = 'magenta', add = F, xlim = c(0.75,2.25), ylim = c(-1.2,.51), ylab = 'NSC <-    xgb pred    -> IPC')
boxplot(ol_marg$pred01[ol_marg$rowname %in% nsc_pc_rn], at = locs[[2]], col = col_key[['NSC']], add = T)
boxplot(ol_marg$pred01[ol_marg$rowname %in% hi_nsc_lo_ipc_rn], at = locs[[3]], col = 'green', add = T)

par(mar = c(5,5,1,1))
locs <- seq(1,3,0.5)
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% anti_neu_rn], at = locs[[1]], col = 'magenta', add = F, xlim = c(0.75,2.35), ylim = c(0.9,2.7), ylab = 'MPRA score')
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% nsc_pc_rn], at = locs[[2]], col = col_key[['NSC']], add = T)
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% hi_nsc_lo_ipc_rn], at = locs[[3]], col = 'green', add = T)
text(c(1.5,2.75), rep(0.43,2), labels = paste0('in ', c('IPC', 'NSC')), xpd = T, col = 'darkgray', cex = 2)

locs <- seq(1,3,l=5)
boxplot(nsc_pca$atac_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.62e+4), ylab = '', col = 'darkgray')
boxplot(nsc_pca$atac_umi[nsc_pca$is_anti_neu == T], at = locs[[2]], add = T, col = 'magenta')
boxplot(nsc_pca$atac_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], at = locs[[3]], add = T, col = 'green')
boxplot(nsc_pca$atac_umi[match(nsc_pc_rn, nsc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['NSC']])
grid(col = 'lightblue', lwd = 2, lty = 2)
title(ylab = 'Proximal ATAC UMI (50kbp)', line = 3)

locs <- seq(1,3,l=5)

boxplot(nsc_pca$rna_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.262e+5), ylab = '', col = 'darkgray', yaxt = 'n')
boxplot(nsc_pca$rna_umi[nsc_pca$is_anti_neu == T], at = locs[[2]], add = T, col = 'magenta', yaxt = 'n')
boxplot(nsc_pca$rna_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], at = locs[[3]], add = T, col = 'green', yaxt = 'n')
boxplot(nsc_pca$rna_umi[match(nsc_pc_rn, nsc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['NSC']], yaxt = 'n')
axis(2, at = seq(0,1.5e+5, 5e+4), labels = signif(seq(0,1.5e+5, 5e+4), 2))
title(ylab = 'Proximal RNA UMI (500kbp)', line = 3)
dev.off()