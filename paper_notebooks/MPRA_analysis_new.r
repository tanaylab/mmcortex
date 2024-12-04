library(xgboost)
library(tgstat)
library(dplyr)
library(matrixStats)

wd <- "/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex"

setwd(wd)

library(metacell)
library(misha)
library(misha.ext)
gsetroot('/home/aviezerl/mm10')

scdb_init(file.path(wd, 'scdb'), f=T)
mat <- scdb_mat('pl_cort')

mc <- scdb_mc('pl_cort')

source('./scripts/util.r')

mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds')

a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))

mcp <- dplyr::select(mca@peaks, chrom, start, end, peak_name)
library(metacell)

scdb_init('scdb', f=T)

mc <- scdb_mc('pl_cort')
mat <- scdb_mat('pl_cort')

tss <- gintervals.load('intervs.global.tss')

tads <- gintervals.load('intervs.global.tad_names')

nei_mcp_tads <- gintervals.neighbors(as.data.frame(mcp), tads, mindist = 0, maxdist = 0, maxneighbors = 1e+6)

orig_lib <- readr::read_tsv('./data/st_and_temporal_and_e14_shadow_enh_11-1-22.tsv')

ol_crs_df <- as.data.frame(cbind(do.call('rbind', stringr::str_split(orig_lib$rowname, '\\.|-')), orig_lib$rowname))
colnames(ol_crs_df) <- c('cell_type', 'pattern', 'number', 'crs_name')
head(ol_crs_df)

load('./output/mcatac/ct_peaks.rda')

nei_orig_ct <- gintervals.neighbors(dplyr::relocate(orig_lib, rowname, .after = end), ct_peaks, mindist = 0, maxdist = 0, maxneighbors = 1e+3)

mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
mcmd <- mcmd[!(as.numeric(mcmd$metacell) %in% 602:603),]

col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
col_key

# Import data
mpra_annot <- read.table('./MPRA_data/st_and_temporal_and_e14_shadow_enh_11-1-22.tsv',header=T)
res <- read.table('./MPRA_data/MPRA_dev_090223.tsv',header=T)


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
save(orig_lib_a, file = './output/MPRA/orig_lib_w_e14_and_shadow_coords.rda')
seqs_tdf <- as.data.frame(stringr::str_split(orig_lib_a$rowname, '_|\\.', simplify = T))
seqs_tdf[seqs_tdf[,3] == '',3] <- NA
rownames(seqs_tdf) <- orig_lib_a$rowname
colnames(seqs_tdf) <- c('ct', 'seqs', 'temp', 'dummy1', 'dummy2')


# monuc_track <- gtrack.ls('seq.G_or_C')

# monuc_seq_coords <- gextract(monuc_track, intervals = as.data.frame(orig_lib_a[,c('chrom', 'start', 'end')]), colnames = gsub('seq\\.', '', monuc_track), iterator = 1)
# monuc_mat <- subset(monuc_seq_coords, select = -c(chrom, start, end, intervalID))
# monuc_per_peak <- as.data.frame(tgs_matrix_tapply(t(monuc_mat), monuc_seq_coords$intervalID, function(x) length(which(x == 1)))/266)
# rownames(monuc_per_peak) <- orig_lib_a$rowname
# monuc_per_peak$A_or_T <- 1 - monuc_per_peak$G_or_C
# save(monuc_per_peak, file = './output/MPRA/monuc_per_peak_orig_lib_w_e14_and_shadow.rda')

load(file = './output/MPRA/monuc_per_peak_orig_lib_w_e14_and_shadow.rda')

# dinuc_tracks <- gtrack.ls('seq\\.[ACGT]{2}$')

# dinucs_seq_coords <- gextract(dinuc_tracks, intervals = as.data.frame(orig_lib_a[,c('chrom', 'start', 'end')]), colnames = gsub('seq\\.', '', dinuc_tracks), iterator = 1)

# saveRDS(object = dinucs_seq_coords, file = './output/MPRA/dinucs_orig_lib.rds')
# dinucs_mat <- subset(dinucs_seq_coords, select = -c(chrom, start, end, intervalID))

# dinucs_per_peak <- tgs_matrix_tapply(t(dinucs_mat), dinucs_seq_coords$intervalID, function(x) length(which(x == 1)))/266
# rownames(dinucs_per_peak) <- orig_lib_a$rowname
# save(dinucs_per_peak, file = './output/MPRA/dinucs_per_peak_orig_lib_w_e14_and_shadow.rda')

load('./output/MPRA/dinucs_per_peak_orig_lib_w_e14_and_shadow.rda')

source('./scripts/util.r')

load('./output/sequence_modeling/xgb_cv_res.rda')
preds_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['contrib_meth']]))
colnames(preds_all)

agg_id <- readr::read_csv('./scatac_data//aggregation_id.csv')
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
                                        'GC content' = 0.47, 'prox_ATAC' = 10, 'prox_RNA' = 10)

pred01 <- setNames(predict(xgb_cv_res[[1]]$bstDMatrix_meth, xgb.DMatrix(as.matrix(dtest_ole_mtt_meth_0_atac_1)), predcontrib=FALSE, predinteraction = FALSE), rownames(dtest_ole_mtt_meth_0_atac_1))

inds_active_mpra <- multintersect(union(rownames(pltmt_nsc[rowMeans(pltmt_nsc) > quantile(rowMeans(pltmt_nsc), 0.75),]),
                             rownames(pltmt_ipc[rowMeans(pltmt_ipc) > quantile(rowMeans(pltmt_ipc), 0.75),])), rownames(pltmt_nsc), rownames(pltmt_ipc))


plot_asterisks <- function(x, y, npi) {
        xl <- sapply(seq_along(npi), function(i) {
            if (npi[[i]] > 1) {
                return(seq(x[[i]] - 0.5*npi[[i]], x[[i]] + 0.5*npi[[i]], l = npi[[i]]))
            } else if (npi[[i]] == 1) {return(x[[i]])}
        })
        points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = 3, col = 'red', lwd= 2)
    }

ks_on_boxplot <- function(xvec, yvec, bins_vec, alternative = 'two.sided') {
    xvcnf <- cut(xvec, bins_vec)
    xvc <- droplevels(xvcnf)
    ksh <- ks.test(yvec[xvc %in% head(levels(xvc), 3)], yvec[xvc %in% tail(levels(xvc), 3)], alternative = alternative)
    print(ksh$p.value)
    LBINX <- match(head(levels(xvc), 3), levels(xvcnf))
    RBINX <- match(tail(levels(xvc), 3), levels(xvcnf))
    Y_LOW_LINE <- 3.3
    Y_HIGH_LINE <- 3.5
    LWD <- 2
    lines(LBINX[c(1,length(LBINX))], rep(Y_LOW_LINE, 2), col = 'red', lwd = LWD)
    lines(RBINX[c(1,length(RBINX))], rep(Y_LOW_LINE, 2), col = 'red', lwd = LWD)
    lines(rep(mean(LBINX), 2), c(Y_LOW_LINE,Y_HIGH_LINE), col = 'red', lwd = LWD)
    lines(rep(mean(RBINX), 2), c(Y_LOW_LINE,Y_HIGH_LINE), col = 'red', lwd = LWD)
    lines(c(mean(LBINX), mean(RBINX)), rep(Y_HIGH_LINE, 2), col = 'red', lwd = LWD)
    powvec <- 5*10**seq(-4,0,1)
    npi <- 4 - which.max(powvec > ksh$p.value)
    npi <- ifelse(npi < 1, 0, npi)
    npi <- ifelse(npi > 3, 3, npi)
    # print('this has run')
    plot_asterisks(mean(c(LBINX, RBINX)), Y_HIGH_LINE + 0.1, npi)
}

ct_tracks <- c('mmcortex_nsc_marginal', 'mmcortex_ipc_marginal', 'mmcortex_cpn_l23_marginal','mmcortex_cpn_l56_marginal','mmcortex_cthpn_marginal','mmcortex_scpn_marginal')
ol_marg <- gextract(ct_tracks,intervals =  orig_lib_a[,c('chrom', 'start', 'end')], iterator = orig_lib_a[,c('chrom', 'start', 'end')])
len_ol_marg <- ol_marg$end - ol_marg$start
ol_marg$start[len_ol_marg < quantile(len_ol_marg, 0.5)] <- ol_marg$start[len_ol_marg < quantile(len_ol_marg, 0.5)] - 1
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

nei_ol_tads <- gintervals.neighbors(ol_marg, tads, mindist = 0, maxdist = 0, maxneighbors = 1e+6)

ol_marg$tad_name <- nei_ol_tads$tad_name[match(ol_marg$rowname, nei_ol_tads$rowname)]

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

nei_orig_peaks <- gintervals.neighbors(orig_lib_a[,c('chrom', 'start', 'end', 'rowname')], 
                                       as.data.frame(mcp),  
                                       mindist = -30e+5, maxdist = 30e+5, maxneighbors = 1e+6)
nei_orig_peaks <- nei_orig_peaks[abs(nei_orig_peaks$dist) >= 1e+3,]

nei_ol_mcp <- gintervals.neighbors(orig_lib_a[,c('chrom', 'start', 'end', 'rowname')], 
                                       as.data.frame(mcp),  
                                       mindist = -500, maxdist = 500, maxneighbors = 1e+6)

# ctxt_rn_intervs$peak_name <- nei_orig_peaks$peak_name[match(ctxt_rn_intervs$rowname, nei_ol_mcp$rowname)]

ipc_marginal_peaks_umis <- Matrix::rowSums(mca@mat[nei_orig_peaks$peak_name,which(mcmd$cell_type == 'IPC')])
nsc_marginal_peaks_umis <- Matrix::rowSums(mca@mat[nei_orig_peaks$peak_name,which(mcmd$cell_type == 'NSC')])
nei_orig_peaks$is_anti_neu <- ifelse(nei_orig_peaks$rowname %in% anti_neu_rn, TRUE, FALSE)
nei_orig_peaks$is_ctxt_ipc <- ifelse(nei_orig_peaks$rowname %in% hi_ipc_lo_nsc_rn, TRUE, FALSE)
nei_orig_peaks$is_ctxt_nsc <- ifelse(nei_orig_peaks$rowname %in% hi_nsc_lo_ipc_rn, TRUE, FALSE)
nei_orig_peaks$is_neu <- ifelse(nei_orig_peaks$rowname %in% neu_rn, TRUE, FALSE)
nei_orig_peaks$is_ipc_pc <- ifelse(nei_orig_peaks$rowname %in% ipc_pc_rn, TRUE, FALSE)
nei_orig_peaks$is_nsc_pc <- ifelse(nei_orig_peaks$rowname %in% nsc_pc_rn, TRUE, FALSE)

nei_orig_peaks$ipc_atac_umis <- ipc_marginal_peaks_umis[nei_orig_peaks$peak_name]
nei_orig_peaks$nsc_atac_umis <- nsc_marginal_peaks_umis[nei_orig_peaks$peak_name]

nei_orig_peaks$tad_name_ol <- nei_ol_tads$tad_name[match(nei_orig_peaks$rowname, nei_ol_tads$rowname)]

nei_orig_peaks$tad_name_mcp <- nei_mcp_tads$tad_name[match(nei_orig_peaks$peak_name, nei_mcp_tads$peak_name)]

get_neighbor_mat <- function(order_metric_name, stat_type, nei_peak_df, crs_rownames, atac_umi_colname = 'ipc_atac_umis', bin_vec = seq(0, 30e+5, 1e+4)) {
    mat_inds <- which(nei_peak_df$rowname %in% crs_rownames)
    sks <- unique(nei_peak_df$rowname[mat_inds])
    metric <- sapply(seq_along(sks), function(i) {
        sk_inds <- which(nei_peak_df$rowname == sks[[i]])
        if (order_metric_name == 'com_abs') {
            y <- sum(abs(nei_peak_df$dist[sk_inds]) * nei_peak_df[sk_inds, atac_umi_colname])/sum(nei_peak_df[sk_inds, atac_umi_colname])
        } else if (order_metric_name == 'com') {
            y <- sum(nei_peak_df$dist[sk_inds] * nei_peak_df[sk_inds, atac_umi_colname])/sum(nei_peak_df[sk_inds, atac_umi_colname])
        } else if (order_metric_name == 'entropy') {
            dist_bins <- cut(nei_peak_df$dist[sk_inds], breaks = bin_vec)
            pv <- tapply(nei_peak_df[sk_inds, atac_umi_colname], dist_bins, sum, na.rm = T)
            pv[is.na(pv)] <- 0
            y <- entropy::entropy(pv)
        } else if (order_metric_name == 'sum') {
            y <- sum(nei_peak_df[sk_inds, atac_umi_colname])
        }
        else {stop('order_metric_name must be either: "com", "com_abs", "entropy" or "sum"')}
        return(y)
    })
    metric_ord <- order(metric)
    # print(length(metric_ord))
    # print(head(sks))
    # print(head(metric_ord))
    # print(head(sks[metric_ord]))
    neighbor_mat <- t(sapply(seq_along(metric_ord), function(i) {
        sk_inds <- which(nei_peak_df$rowname == sks[[metric_ord[[i]]]])
        dist_bins <- cut(nei_peak_df$dist[sk_inds], breaks = bin_vec)
        if (stat_type == 'sum_umis') {
            pv <- tapply(nei_peak_df[sk_inds, atac_umi_colname], dist_bins, sum, na.rm = T)   
        } else if (stat_type == 'num_neighbors') {
            pv <- table(dist_bins)
        }
        pv[is.na(pv)] <- 0
        return(pv)
    }))
    rownames(neighbor_mat) <- sks[metric_ord]
    colnames(neighbor_mat) <- head(bin_vec, -1)
    return(neighbor_mat)
}

neu_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = neu_rn)
ipc_pc_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = ipc_pc_rn)
ipc_ctxt_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = hi_ipc_lo_nsc_rn)
neu_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = neu_rn)
ipc_pc_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = ipc_pc_rn)
ipc_ctxt_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = hi_ipc_lo_nsc_rn)
anti_neu_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
nsc_pc_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
nsc_ctxt_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')
anti_neu_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
nsc_pc_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
nsc_ctxt_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')

# neu_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = neu_rn)
# ipc_pc_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = ipc_pc_rn)
# ipc_ctxt_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = hi_ipc_lo_nsc_rn)
# neu_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = neu_rn)
# ipc_pc_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = ipc_pc_rn)
# ipc_ctxt_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = hi_ipc_lo_nsc_rn)
# anti_neu_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
# nsc_pc_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
# nsc_ctxt_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')
# anti_neu_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = anti_neu_rn, atac_umi_colname = 'nsc_atac_umis')
# nsc_pc_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = nsc_pc_rn, atac_umi_colname = 'nsc_atac_umis')
# nsc_ctxt_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = hi_nsc_lo_ipc_rn, atac_umi_colname = 'nsc_atac_umis')

rand_crs <- sample(ol_marg$rowname, round(mean(c(length(neu_rn),length(ipc_pc_rn),length(hi_ipc_lo_nsc_rn)))))

rand_umi_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = nei_orig_peaks, crs_rownames = rand_crs)
# rand_umi_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'sum_umis', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = rand_crs)

rand_num_neighbor_mat <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = nei_orig_peaks, crs_rownames = rand_crs)
# rand_num_neighbor_mat_tad <- get_neighbor_mat(order_metric_name = 'com_abs', stat_type = 'num_neighbors', nei_peak_df = dplyr::filter(nei_orig_peaks, tad_name_ol == tad_name_mcp), crs_rownames = rand_crs)


nsc_pca <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]), 
                                       restrict_to_tads = F,
                                       mc_sel = which(mcmd$cell_type == 'NSC'), mca = mca , mat_rna = mat , mc_rna = mc, tss = tss, tads = tads
                           )

ipc_pca <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]),
                                       restrict_to_tads = F,
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

# nsc_pca_tad <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]), 
#                                            restrict_to_tads = T,
#                                        mc_sel = which(mcmd$cell_type == 'NSC'), mca = mca , mat_rna = mat , mc_rna = mc, tss = tss, tads = tads
#                            )

# ipc_pca_tad <- proximal_chromatin_activity(peaks_of_interest = ctxt_rn_intervs[,c(1:3,5)],background_peaks = as.data.frame(mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]),
#                                            restrict_to_tads = T,
#                                        mc_sel = which(mcmd$cell_type == 'IPC'), mca = mca , mat_rna = mat , mc_rna = mc, tss = tss, tads = tads
#                            )
# nsc_pca_tad$rowname <- ctxt_rn_intervs$rowname[match(nsc_pca$peak_name, ctxt_rn_intervs$peak_name)]
# nsc_pca_tad$new_class <- ol_marg$new_class[match(nsc_pca$rowname, ol_marg$rowname)]
# ipc_pca_tad$rowname <- ctxt_rn_intervs$rowname[match(ipc_pca$peak_name, ctxt_rn_intervs$peak_name)]
# ipc_pca_tad$new_class <- ol_marg$new_class[match(ipc_pca$rowname, ol_marg$rowname)]
# nsc_pca_tad$is_ctxt <- ifelse(nsc_pca$rowname %in% hi_ipc_lo_nsc_rn, TRUE, FALSE)
# ipc_pca_tad$is_ctxt <- ifelse(ipc_pca$rowname %in% hi_nsc_lo_ipc_rn, TRUE, FALSE)
# nsc_pca_tad$is_neu <- ifelse(nsc_pca$rowname %in% neu_rn, TRUE, FALSE)
# ipc_pca_tad$is_neu <- ifelse(ipc_pca$rowname %in% neu_rn, TRUE, FALSE)
# nsc_pca_tad$is_anti_neu <- ifelse(nsc_pca$rowname %in% anti_neu_rn, TRUE, FALSE)
# ipc_pca_tad$is_anti_neu <- ifelse(ipc_pca$rowname %in% anti_neu_rn, TRUE, FALSE)

plot_asterisks <- function(x, y, npi) {
    xl <- sapply(seq_along(npi), function(i) {
        if (npi[[i]] > 1) {
            return(seq(x[[i]] - 0.0751*npi[[i]], x[[i]] + 0.0751*npi[[i]], l = npi[[i]]))
        } else if (npi[[i]] == 1) {return(x[[i]])}
    })
    points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = 3, col = 'red', lwd= 3)
}

add_ks_stars_to_boxplot <- function(vec1, x1, vec2, x2, y_line, delta_y_asterisk, alternative = 'two.sided') {
    ksr <- ks.test(vec1, vec2, alternative = alternative)
    ksp <- ksr$p.value
    print(ksr$statistic)
    powvec <- 5*10**seq(-4,0,1)
    npi <- 4 - which.max(powvec > ksp)
    npi <- ifelse(npi < 1, 0, npi)
    npi <- ifelse(npi > 3, 3, npi)
    lines(c(x1, x2), rep(y_line, 2), col = 'red', lwd = 2)
    plot_asterisks(mean(c(x1, x2)), y_line+delta_y_asterisk, npi)
}



save(pltmt_nsc, pltmt_ipc,
        pred01,
        inds_active_mpra,
        ks_on_boxplot,
        boxplot_vec,
        nei_orig_ct,
        neu_num_neighbor_mat,
        ipc_pc_num_neighbor_mat,
        ipc_ctxt_num_neighbor_mat,
        anti_neu_num_neighbor_mat,
        nsc_pc_num_neighbor_mat,
        nsc_ctxt_num_neighbor_mat,
        neu_umi_neighbor_mat,
        ipc_pc_umi_neighbor_mat,
        ipc_ctxt_umi_neighbor_mat,
        anti_neu_umi_neighbor_mat,
        nsc_pc_umi_neighbor_mat,
        nsc_ctxt_umi_neighbor_mat,
        rand_num_neighbor_mat,
        rand_umi_neighbor_mat,
        ol_marg,
        ipc_pca,
        nsc_pca,
        ipc_pc_rn,
        hi_ipc_lo_nsc_rn,
        neu_rn,
        nc_neu_rn, 
        nc_hi_ipc_lo_nsc_rn, 
        add_ks_stars_to_boxplot, 
        plot_asterisks,
        anti_neu_rn, 
        hi_nsc_lo_ipc_rn, 
        nc_hi_nsc_lo_ipc_rn,
        nsc_pc_rn,
        file = file.path(wd, 'output/MPRA/fig_7_data.rda')
        )








# nei_tads_peaks <- gintervals.neighbors(tads, as.data.frame(mcp), maxdist = 0, mindist = 0, maxneighbors = 1e+4)

# nei_genes_tads <- gintervals.neighbors(tss[!duplicated(tss$geneSymbol) & tss$geneSymbol %in% rownames(mc@e_gc),], tads,maxdist = 0, mindist = 0, maxneighbors = 1e+4)

# mc_sel <- mcmd$metacell[mcmd$cell_type == 'NSC']

# marginal_peaks_umis <- Matrix::rowSums(mca@mat[mcp$peak_name,mc_sel])

# marginal_genes_umis <- Matrix::rowSums(mat@mat[,names(mc@mc[mc@mc %in% mc_sel])])

# head(marginal_genes_umis)

# head(marginal_peaks_umis)

# proximal_chromatin_activity <- function(peaks_of_interest,
#                                         background_peaks,
#                                         mc_sel,
#                                         # day,
#                                         mca,
#                                         mat_rna,
#                                         mc_rna,
#                                         tss,
#                                         tads,
#                                         d_puncture = 1e+3,
#                                         peak_clusters = NULL,
#                                         pred = NULL,
#                                         d_proximity_atac = 5e+4, 
#                                         d_proximity_rna = 5e+5, 
#                                         eps = 1e-5,
#                                         restrict_to_tads = TRUE) {
#     if (!(tibble::has_name(peaks_of_interest, 'peak_name') & tibble::has_name(background_peaks, 'peak_name'))) {
#         stop('Peaks of interest or background peaks do not have peak name')
#     }
#     # if (length(which(!(peaks_of_interest$peak_name %in% background_peaks$peak_name))) > 0) {
#     #     stop('Not all peaks of interest are in background peak set')
#     # }
#     peaks_of_interest <- peaks_of_interest[!duplicated(peaks_of_interest$peak_name),]
#     background_peaks <- background_peaks[!duplicated(background_peaks$peak_name),]
#     nei_peaks_peaks <- gintervals.neighbors(peaks_of_interest, background_peaks, maxdist = d_proximity_atac, mindist = -d_proximity_atac, maxneighbors = 1e+6)
#     colnames(nei_peaks_peaks)[grep('peak_name', colnames(nei_peaks_peaks))] <- c('peak_name_1', 'peak_name_2')
#     nei_peaks_peaks <- nei_peaks_peaks[nei_peaks_peaks$peak_name_1 != nei_peaks_peaks$peak_name_2 & abs(nei_peaks_peaks$dist) >= d_puncture,]

#     nei_bg_peaks_tads <- gintervals.neighbors(tads, background_peaks, maxdist = 0, mindist = 0, maxneighbors = 1e+4)
#     nei_fg_peaks_tads <- gintervals.neighbors(tads, peaks_of_interest, maxdist = 0, mindist = 0, maxneighbors = 1e+4)
#     marginal_peaks_umis <- Matrix::rowSums(mca@mat[background_peaks$peak_name,mc_sel])
#     nei_peaks_peaks$tads1 <- nei_fg_peaks_tads$tad_name[match(nei_peaks_peaks[,4], nei_fg_peaks_tads$peak_name)]
#     nei_peaks_peaks$tads2 <- nei_bg_peaks_tads$tad_name[match(nei_peaks_peaks[,8], nei_bg_peaks_tads$peak_name)]
#     if (restrict_to_tads) {
#         nei_peaks_peaks <- nei_peaks_peaks[which(nei_peaks_peaks$tads1 == nei_peaks_peaks$tads2),]
#     }
#     nei_peaks_peaks$atac_umis_2 <- marginal_peaks_umis[match(nei_peaks_peaks$peak_name_2, names(marginal_peaks_umis))]
#     tad_atac_umis <- tapply(nei_peaks_peaks$atac_umis_2, nei_peaks_peaks$peak_name_1, sum)
    
#     npgf_sum_atac <- tad_atac_umis
    
#     npgf_sum_atac[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% names(npgf_sum_atac))]] <- 0
#     npgf_sum_atac[is.na(npgf_sum_atac)] <- 0
    
#     npgf_sum_atac <- npgf_sum_atac[peaks_of_interest$peak_name]
    

    
#     nei_peaks_genes <- gintervals.neighbors(peaks_of_interest, tss, maxdist = d_proximity_rna, mindist = -d_proximity_rna, maxneighbors = 1e+4)
    
#     # print(head(nei_peaks_peaks))
#     peaks_tads <- gintervals.neighbors(nei_peaks_genes[,1:4], tads, maxdist = 0, mindist = 0, maxneighbors = 1)
#     genes_tads <- gintervals.neighbors(dplyr::rename(nei_peaks_genes[,c(5:7, grep('geneSymbol', colnames(nei_peaks_genes)))], 
#                                             chrom = chrom1, start = start1, end = end1),
#                                        tads, maxdist = 0, mindist = 0, maxneighbors = 1)
#     nei_peaks_genes$peak_tad <- peaks_tads$tad_name
#     nei_peaks_genes$gene_tad <- genes_tads$tad_name
#     if (restrict_to_tads) {
#         npgf <- dplyr::select(nei_peaks_genes[which(nei_peaks_genes$peak_tad == nei_peaks_genes$gene_tad),], 
#                       c(1:4, grep('geneSymbol|peak_tad|gene_tad', colnames(nei_peaks_genes))))
#     } else {
#         npgf <- dplyr::select(nei_peaks_genes, 
#                       c(1:4, grep('geneSymbol|peak_tad|gene_tad', colnames(nei_peaks_genes))))
#     }
    
#     genes_nei_umis <- Matrix::rowSums(mat_rna@mat[unique(npgf$geneSymbol[npgf$geneSymbol %in% rownames(mat_rna@mat)]),
#                                                       names(mc_rna@mc[mc_rna@mc %in% mc_sel])])
#     npgf$rna_umis <- genes_nei_umis[npgf$geneSymbol]
#     npgf$rna_umis[is.na(npgf$rna_umis)] <- 0
#     npgf_sum_rna <- tapply(npgf$rna_umis, npgf$peak_name, function(x) sum(x))
#     npgf_sum_rna[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% npgf$peak_name)]] <- 0
#     npgf_sum_rna[is.na(npgf_sum_rna)] <- 0
#     npgf_sum_rna <- npgf_sum_rna[peaks_of_interest$peak_name]
    
    

# tads_atac_umis <- tapply(marginal_peaks_umis, nei_tads_peaks$tad_name[match(names(marginal_peaks_umis), nei_tads_peaks$peak_name)], sum)

# tads_rna_umis <- tapply(marginal_genes_umis, nei_genes_tads$tad_name[match(names(marginal_genes_umis), nei_genes_tads$geneSymbol)], sum)

# tb <- intersect(names(tads_atac_umis), names(tads_rna_umis))

# quantile(log2(1+tads_atac_umis[setdiff(names(tads_atac_umis), tb)]))

# setdiff(names(tads_rna_umis), names(tads_atac_umis))

# colnames(ol_marg)



# # tads_neu_rn <- intersect(ol_marg$tad_name[ol_marg$rowname %in% neu_rn], tb)
# # tads_ipc_pc_rn <- intersect(ol_marg$tad_name[ol_marg$rowname %in% ipc_pc_rn], tb)
# # tads_hi_ipc_lo_nsc <- intersect(ol_marg$tad_name[ol_marg$rowname %in% hi_ipc_lo_nsc_rn], tb)
# tads_neu_rn <- ol_marg$tad_name[ol_marg$rowname %in% neu_rn]
# tads_ipc_pc_rn <- ol_marg$tad_name[ol_marg$rowname %in% ipc_pc_rn]
# tads_hi_ipc_lo_nsc <- ol_marg$tad_name[ol_marg$rowname %in% hi_ipc_lo_nsc_rn]

# length(tads_neu_rn)
# length(intersect(tads_neu_rn, tb))

# length(tads_ipc_pc_rn)
# length(intersect(tads_ipc_pc_rn, tb))

# length(tads_hi_ipc_lo_nsc)
# length(intersect(tads_hi_ipc_lo_nsc, tb))

# quantile(log2(1+tads_atac_umis))

# # png('./output/MPRA/figs/umis_per_tads.png', h = 800, w = 1200, res = 100)
# pdf('./output/MPRA/figs/umis_per_tads.pdf', h = 8, w = 12)
# options(repr.plot.width = 18)
# options(repr.plot.height = 6)
# q_atac <- quantile(log2(1+tads_atac_umis[tb]), c(0.1,0.9))
# q_rna <- quantile(log2(1+tads_rna_umis[tb]), c(0.1,0.9))

# par(mfrow = c(2,3))
# brks <- seq(4,19,l=100)
# brks_rna <- seq(0,22,l=100)
# hist(log2(1+tads_atac_umis[tb]), breaks = brks, main = 'ATAC UMIs per TAD', xlab = 'Log2 ATAC UMIs')
# hist(log2(1+tads_atac_umis[setdiff(names(tads_atac_umis), tb)]), breaks = brks, col = 'lightblue', add = T)
# legend('topleft', legend = c('tads w RNA', 'tads w/o RNA'), col = c('gray', 'lightblue'), pch = 15)
# hist(log2(1+tads_rna_umis[tb]), , breaks = brks_rna, main = 'RNA UMIs per TAD', xlab = 'Log2 RNA UMIs')
# hist(log2(1+tads_rna_umis[setdiff(names(tads_rna_umis), tb)]), breaks = brks_rna, col = 'pink', add = T)
# legend('topleft', legend = c('tads w ATAC', 'tads w/o ATAC'), col = c('gray', 'pink'), pch = 15)

# plot(log2(1+tads_atac_umis[tb]), log2(1+tads_rna_umis[tb]), cex = 1, pch = 1, xlab = 'Log2 ATAC UMIs per TAD', ylab = 'Log2 RNA UMIs per TAD')
# points(log2(1+tads_atac_umis[tads_neu_rn]), log2(1+tads_rna_umis[tads_neu_rn]), cex = 1, pch = 16, col = 'red')
# points(log2(1+tads_atac_umis[tads_ipc_pc_rn]), log2(1+tads_rna_umis[tads_ipc_pc_rn]), cex = 1, pch = 16, col = 'blue')
# points(log2(1+tads_atac_umis[tads_hi_ipc_lo_nsc]), log2(1+tads_rna_umis[tads_hi_ipc_lo_nsc]), cex = 1, pch = 16, col = 'orange3')
# # dev.off()

# # png('./output/MPRA/figs/umis_per_tads_density.png', h = 400, w = 1200, res = 100)
# # options(repr.plot.width = 12)
# # options(repr.plot.height = 6)
# # par(mfrow = c(1,2), mar = c(5,4,1,1))
# par(mar = c(5,4,1,1))
# BW_ATAC <- .45
# den_neu_atac <- density(log2(1+tads_atac_umis[tads_neu_rn]), bw = BW_ATAC)
# den_ipc_pc_atac <- density(log2(1+tads_atac_umis[tads_ipc_pc_rn]), bw = BW_ATAC)
# den_hi_ipc_lo_nsc_atac <- density(log2(1+tads_atac_umis[tads_hi_ipc_lo_nsc]), bw = BW_ATAC)

# plot(den_neu_atac$x, den_neu_atac$y, col = 'red', type = 'l', lwd = 2, ylim = c(0,.41), xlab = 'log2 sum ATAC UMIs', ylab = 'Density')
# lines(den_ipc_pc_atac$x, den_ipc_pc_atac$y, lwd = 2, col = 'blue')
# lines(den_hi_ipc_lo_nsc_atac$x, den_hi_ipc_lo_nsc_atac$y, lwd = 2, col = 'orange3')

# filt_na <- function(x) {x[!is.na(x)]}

# BW_RNA <- .75
# den_neu_rna <- density(filt_na(log2(1+tads_rna_umis[tads_neu_rn])), bw = BW_RNA)
# den_ipc_pc_rna <- density(filt_na(log2(1+tads_rna_umis[tads_ipc_pc_rn])), bw = BW_RNA)
# den_hi_ipc_lo_nsc_rna <- density(filt_na(log2(1+tads_rna_umis[tads_hi_ipc_lo_nsc])), bw = BW_RNA)

# plot(den_neu_rna$x, den_neu_rna$y, col = 'red', type = 'l', lwd = 2, ylim = c(0,.22), xlab = 'log2 sum RNA UMIs', ylab = 'Density')
# lines(den_ipc_pc_rna$x, den_ipc_pc_rna$y, lwd = 2, col = 'blue')
# lines(den_hi_ipc_lo_nsc_rna$x, den_hi_ipc_lo_nsc_rna$y, lwd = 2, col = 'orange3')
# dev.off()



# library(iceqream)

# apply(dtest_ole_mtt_meth_0_atac_1[train_peaks,17:ncol(dtest_ole_mtt_meth_0_atac_1)], 2, function(x) class(x))
# apply(apply(as.matrix(ol_marg[match(train_peaks, ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x))), 2, function(x) class(x))
# apply(as.matrix(dtest_ole_mtt_meth_0_atac_1[train_peaks,1:16]), 2, function(x) class(x))

# save(all_peaks_here, ol_marg, dtest_ole_mtt_meth_0_atac_1, file = './output/MPRA/vars_for_blue_red_reg.rda')

# normalizePath('./output/MPRA/vars_for_blue_red_reg.rda')

# load('./output/MPRA/vars_for_blue_red_reg.rda')

# library(glmnet)

# devtools::load_all("~/src/iceqream")

# all_peaks_here <- multunion(neu_rn,ipc_pc_rn)

# all_peaks_here <- union(neu_rn,ipc_pc_rn, hi_ipc_lo_nsc_rn)
# set.seed(1337)
# test_peaks <- sample(all_peaks_here, size = round(0.1*length(all_peaks_here)))
# train_peaks <- setdiff(all_peaks_here, test_peaks)
# rvb_tm <- regress_trajectory_motifs(peak_intervals = dplyr::mutate(as.data.frame(ol_marg[ol_marg$rowname %in% train_peaks,]), const = F), seed = 1337,
#                                     additional_features = dtest_ole_mtt_meth_0_atac_1[train_peaks,17:ncol(dtest_ole_mtt_meth_0_atac_1)],
#                                     atac_scores = apply(as.matrix(ol_marg[match(train_peaks, ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x))), 
#                                     motif_energies = as.matrix(dtest_ole_mtt_meth_0_atac_1[train_peaks,1:16]), 
#                                     max_motif_num=16, n_prego_motifs = 4,spat_bin_size = 10)

# load('./output/MPRA/red_vs_blue_iq_model.rda')

# rvb_tm_test <- iceqream::infer_trajectory_motifs(traj_model = rvb_tm, 
#                                   peak_intervals = dplyr::mutate(as.data.frame(ol_marg[ol_marg$rowname %in% test_peaks,]), const = F), 
#                                     additional_features = dtest_ole_mtt_meth_0_atac_1[test_peaks,17:ncol(dtest_ole_mtt_meth_0_atac_1)],
#                                     atac_scores = apply(as.matrix(ol_marg[match(test_peaks, ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x))))
#                                     #                     , 
#                                     # motif_energies = as.matrix(dtest_ole_mtt_meth_0_atac_1[test_peaks,1:16]))

# cor(rvb_tm_test@predicted_diff_score,rvb_tm_test@diff_score, method = 'spearman')

# options(repr.plot.width = 12)
# options(repr.plot.height = 6)

# load(file = './output/MPRA/rvb_iq_model_cv_ls.rda')

# iq_ev_preds <- lapply(seq_along(iq_cv_ls), function(i) {
#     pred_i <- iceqream::infer_trajectory_motifs(traj_model = iq_cv_ls[[i]]$tm_i, 
#                                   peak_intervals = dplyr::mutate(as.data.frame(ol_marg[ol_marg$rowname %in% iq_cv_ls[[i]]$test_peaks,]), const = F), 
#                                     additional_features = dtest_ole_mtt_meth_0_atac_1[iq_cv_ls[[i]]$test_peaks,17:ncol(dtest_ole_mtt_meth_0_atac_1)],
#                                     atac_scores = apply(as.matrix(ol_marg[match(iq_cv_ls[[i]]$test_peaks, ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x))))
#     return(pred_i)
# })

# test_preds_all <- do.call('c', lapply(seq_along(iq_ev_preds), function(i) {
#     test_peaks <- iq_cv_ls[[i]]$test_peaks
#     test_preds <- iq_ev_preds[[i]]@predicted_diff_score[test_peaks]
#     return(test_preds)
# }))

# test_obs_all <- do.call('c', lapply(seq_along(iq_ev_preds), function(i) {
#     test_peaks <- iq_cv_ls[[i]]$test_peaks
#     test_preds <- setNames(iq_ev_preds[[i]]@diff_score[match(test_peaks, rownames(iq_ev_preds[[i]]@additional_features))], test_peaks)
#     return(test_preds)
# }))

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# y_test <- test_preds_all
# # x_test <- apply(apply(as.matrix(ol_marg[match(names(test_preds_all), ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x))), 1, diff)
# x_test <- test_obs_all

# plot(x_test, y_test, col = ifelse(names(test_preds_all) %in% ipc_pc_rn, 'blue', 'red'), pch = 16, main = paste0('Test data (5-fold CV), n = ', length(test_preds_all)), xlab = 'Observed', ylab = 'Predicted')
# text(-0.52,0.25, labels = paste0('r = ', signif(cor(x_test, y_test, method = 'spearman'), 2)))
# abline(0,1,col='red')

# mmm <- apply(as.matrix(ol_marg[match(iq_cv_ls[[i]]$test_peaks, ol_marg$rowname),c('log_nsc','log_ipc')]), 2, function(x) (x - min(x))/(max(x) - min(x)))

# rownames(mmm) <- test_peaks

# plot(mmm[,1], mmm[,2], col = ifelse(rownames(mmm) %in% ipc_pc_rn, 'blue', 'red'))

# # par(mfrow = c(1,2), mar = c(5,4,2,1))
# # plot(rvb_tm@predicted_diff_score,rvb_tm@diff_score, col = ifelse(rownames(rvb_tm@additional_features) %in% ipc_pc_rn, 'blue', 'red'), pch = 16, main = paste0('Training data, n = ', length(train_peaks)), xlab = 'Observed', ylab = 'Predicted')
# # text(-0.6,0.35, labels = paste0('r = ', signif(cor(rvb_tm@predicted_diff_score,rvb_tm@diff_score), 2)))
# # abline(0,1,col='red')
# # x_test <- rvb_tm_test@predicted_diff_score[test_peaks]
# # y_test <- rvb_tm_test@diff_score[match(test_peaks, rownames(rvb_tm_test@additional_features))]
# # plot(x_test, y_test, col = ifelse(test_peaks %in% ipc_pc_rn, 'blue', 'red'), pch = 16, main = paste0('Test data (no CV), n = ', length(test_peaks)), xlab = 'Observed', ylab = 'Predicted')
# # text(-0.52,0.35, labels = paste0('r = ', signif(cor(x_test, y_test), 2)))
# # abline(0,1,col='red')

# iceqream::plot_traj_model_report(rvb_tm, filename = './output/MPRA/red_vs_blue_iq_model_report.pdf')

# prego_pssms <- lapply()



# source('./scripts/util.r')

# ct_tracks



# cfupn_peaks_df <-amos_peak_calling_function(trk = "mmcortex_cthpn_marginal", canonical = F)
# cpnl23_peaks_df <-amos_peak_calling_function(trk = "mmcortex_cpn_l23_marginal", canonical = F)

# astro_peaks_df <-amos_peak_calling_function(trk = "mmcortex_astrocytes_marginal", canonical = F)
# oligo_peaks_df <-amos_peak_calling_function(trk = "mmcortex_OPCs_marginal", canonical = F)
# nsc_peaks_df <-amos_peak_calling_function(trk = "mmcortex_nsc_marginal", canonical = F)

# nei_cfupn_mcp <- gintervals.neighbors(cfupn_peaks_df, as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+6)
# nei_cpn_mcp <- gintervals.neighbors(cpnl23_peaks_df, as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+6)
# nei_astro_mcp <- gintervals.neighbors(astro_peaks_df, as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+6)
# nei_oligo_mcp <- gintervals.neighbors(oligo_peaks_df, as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+6)
# nei_nsc_mcp <- gintervals.neighbors(nsc_peaks_df, as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+6)

# cfupn_peaks_df$is_in <- cfupn_peaks_df$intervalID %in% nei_cfupn_mcp$intervalID
# cpnl23_peaks_df$is_in <- cpnl23_peaks_df$intervalID %in% nei_cpn_mcp$intervalID
# astro_peaks_df$is_in <- astro_peaks_df$intervalID %in% nei_astro_mcp$intervalID
# oligo_peaks_df$is_in <- oligo_peaks_df$intervalID %in% nei_oligo_mcp$intervalID

# nsc_peaks_df$is_in <- nsc_peaks_df$intervalID %in% nei_nsc_mcp$intervalID

# length(which(nsc_peaks_df$mmcortex_nsc_marginal[!nsc_peaks_df$is_in] > mean(nsc_peaks_df$mmcortex_nsc_marginal[nsc_peaks_df$is_in])))

# par(mfrow = c(1,5))
# boxplot(mmcortex_cthpn_marginal ~ is_in, data = cfupn_peaks_df, main = 'CthPN')
# boxplot(mmcortex_cpn_l23_marginal ~ is_in, data = cpnl23_peaks_df, main = 'CPN_L2-3')
# boxplot(mmcortex_astrocytes_marginal ~ is_in, data = astro_peaks_df, main = 'Astrocytes')
# boxplot(mmcortex_OPCs_marginal ~ is_in, data = oligo_peaks_df, main = 'OPCs')
# boxplot(mmcortex_nsc_marginal ~ is_in, data = nsc_peaks_df, main = 'OPCs')

# table(mcmd$cell_type)
# # /sum(table(mcmd$cell_type))

# nrow(cfupn_peaks_df)

# length(setdiff(nsc_peaks_df$intervalID, nei_nsc_mcp$intervalID))

# length(setdiff(cfupn_peaks_df$intervalID, nei_cfupn_mcp$intervalID))

# length(setdiff(cpnl23_peaks_df$intervalID, nei_cpn_mcp$intervalID))

# length(setdiff(astro_peaks_df$intervalID, nei_astro_mcp$intervalID))

# length(setdiff(oligo_peaks_df$intervalID, nei_oligo_mcp$intervalID))

# table(nsc_peaks_df$is_in)

# length(which(nsc_peaks_df$mmcortex_nsc_marginal[!nsc_peaks_df$is_in] > mean(nsc_peaks_df$mmcortex_nsc_marginal[nsc_peaks_df$is_in])))

# length(which(cfupn_peaks_df$mmcortex_cthpn_marginal[!cfupn_peaks_df$is_in] > mean(cfupn_peaks_df$mmcortex_cthpn_marginal[cfupn_peaks_df$is_in])))

# length(which(cpnl23_peaks_df$mmcortex_cpn_l23_marginal[!cpnl23_peaks_df$is_in] > mean(cpnl23_peaks_df$mmcortex_cpn_l23_marginal[cpnl23_peaks_df$is_in])))

# length(which(astro_peaks_df$mmcortex_astrocytes_marginal[!astro_peaks_df$is_in] > mean(astro_peaks_df$mmcortex_astrocytes_marginal[astro_peaks_df$is_in])))

# length(which(oligo_peaks_df$mmcortex_OPCs_marginal[!oligo_peaks_df$is_in] > mean(oligo_peaks_df$mmcortex_OPCs_marginal[oligo_peaks_df$is_in])))

# nei_astro_not_in_oligo_not_in <- gintervals.neighbors(astro_peaks_df[!astro_peaks_df$is_in,], oligo_peaks_df[!oligo_peaks_df$is_in,], mindist = 0, maxdist = 0, maxneighbors = 1e+6)

# nei_cfupn_not_in_cpn_not_in <- gintervals.neighbors(cfupn_peaks_df[!cfupn_peaks_df$is_in,], cpnl23_peaks_df[!cpnl23_peaks_df$is_in,], mindist = 0, maxdist = 0, maxneighbors = 1e+6)

# nrow(astro_peaks_df[!astro_peaks_df$is_in,])
# nrow(oligo_peaks_df[!oligo_peaks_df$is_in,])

# nrow(cfupn_peaks_df[!cfupn_peaks_df$is_in,])
# nrow(cpnl23_peaks_df[!cpnl23_peaks_df$is_in,])

# nrow(nei_astro_not_in_oligo_not_in)

# nrow(nei_cfupn_not_in_cpn_not_in)

# head(nei_astro_not_in_oligo_not_in)

# cor(nei_cfupn_not_in_cpn_not_in$mmcortex_cthpn_marginal, nei_cfupn_not_in_cpn_not_in$mmcortex_cpn_l23_marginal, method = 'spearman')
# plot(nei_cfupn_not_in_cpn_not_in$mmcortex_cthpn_marginal, nei_cfupn_not_in_cpn_not_in$mmcortex_cpn_l23_marginal)

# quantile(nei_cfupn_not_in_cpn_not_in$mmcortex_cthpn_marginal)
# quantile(nei_cfupn_not_in_cpn_not_in$mmcortex_cpn_l23_marginal)

# quantile(cfupn_peaks_df$mmcortex_cthpn_marginal)
# quantile(cpnl23_peaks_df$mmcortex_cpn_l23_marginal)

# cor(nei_astro_not_in_oligo_not_in$mmcortex_astrocytes_marginal, nei_astro_not_in_oligo_not_in$mmcortex_OPCs_marginal, method = 'spearman')
# plot(nei_astro_not_in_oligo_not_in$mmcortex_astrocytes_marginal, nei_astro_not_in_oligo_not_in$mmcortex_OPCs_marginal)

# gquantiles(expr = 'mmcortex_cthpn_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)

# gquantiles(expr = 'mmcortex_cpn_l23_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)

# gquantiles(expr = 'mmcortex_cpn_l56_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)

# gquantiles(expr = 'mmcortex_scpn_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)

# gquantiles(expr = 'mmcortex_ipc_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)

# gquantiles(expr = 'mmcortex_nsc_marginal', percentiles = c(seq(0,900,100), 950, 990, 999, 1000)/1000, iterator = 20)



# head(cfupn_peaks)

# cfupn_peaks_can <-amos_peak_calling_function(trk = "mmcortex_cthpn_marginal", q = .95, canonical = T)

# cpn_peaks <-amos_peak_calling_function(trk = "mmcortex_cthpn_marginal", q = .95, canonical = F)





# head(ct_tracks_on_manifold)

# ct_tracks <- c('mmcortex_nsc_marginal', 'mmcortex_ipc_marginal')
# ol_marg <- gextract(ct_tracks,intervals =  orig_lib_a[,c('chrom', 'start', 'end')], iterator = orig_lib_a[,c('chrom', 'start', 'end')])
# mcp_trim <- dplyr::mutate(mcp[,c('chrom', 'start', 'end')],start = start + 17, end = end - 17)

# neighboring_intervs <- dplyr::mutate(orig_lib_a[match(neu_rn, orig_lib_a$rowname),], start = start - 5e+4 + (end - start)/2, end = end + 5e+4 - (end - start)/2)

# head(neighboring_intervs)

# ct_tracks_on_neighboring_intervals <- gextract(ct_tracks,intervals =  neighboring_intervs[,2:4], iterator =neighboring_intervs[,2:4])

# head(ct_tracks_on_neighboring_intervals)

# plot(ct_tracks_on_neighboring_intervals$mmcortex_nsc_marginal, ct_tracks_on_neighboring_intervals$mmcortex_ipc_marginal)
# abline(0,1,col = 'red')

# plot(ct_tracks_on_manifold$mmcortex_nsc_marginal, ct_tracks_on_manifold$mmcortex_ipc_marginal, xlim = c(0,42), ylim = c(0,42))
# abline(0,1,col = 'red')

# head(nsc_pca)

# rnb <- intersect(nsc_pca$peak_name, ipc_pca$peak_name)
# plot(unlist(nsc_pca[match(rnb, nsc_pca$peak_name),'atac_umi']), unlist(ipc_pca[match(rnb, ipc_pca$peak_name),'atac_umi']))
#      # , xlim = c(0,42), ylim = c(0,42)
#     # )
# abline(0,1,col = 'red')

# filt_na <- function(x) {y <- x; return(y[!is.na(y)])}

# table(nsc_pca$peak_name %in% mcp$peak_name)

# rnb <- intersect(nsc_pca$peak_name, ipc_pca$peak_name)
# mcp_rnb <- filt_na(setNames(nei_orig_mcp$peak_name[match(nsc_pca$rowname[match(rnb, nsc_pca$peak_name)], nei_orig_mcp$rowname)], rnb))
# smoothScatter(log2(1+unlist(nsc_pca[match(rnb, nsc_pca$peak_name),'rna_umi'])), log2(1+unlist(ipc_pca[match(rnb, ipc_pca$peak_name),'rna_umi'])))
#      # , xlim = c(0,42), ylim = c(0,42)
#     # )
# abline(0,1,col = 'red')
# abline(1,1,col = 'green')
# abline(-1,1,col = 'green')

# rna_diff_vec <- setNames(log2(1+unlist(ipc_pca[match(rnb, ipc_pca$peak_name),'rna_umi'])) - log2(1+unlist(nsc_pca[match(rnb, nsc_pca$peak_name),'rna_umi'])), rnb)
# crs_diff_rna <- names(rna_diff_vec[abs(rna_diff_vec) >= 2])

# atac_diff_vec <- setNames(unlist(ipc_pca[match(rnb, ipc_pca$peak_name),'atac_umi']) - unlist(nsc_pca[match(rnb, nsc_pca$peak_name),'atac_umi']), rnb)
# crs_diff_atac <- names(atac_diff_vec[abs(atac_diff_vec) >= 5e+3])

# length(crs_diff_rna)
# length(crs_diff_atac)

# length(intersect(crs_diff_rna,crs_diff_atac))

# plot(atac_diff_vec[crs_diff_rna], rna_diff_vec[crs_diff_rna])

# a_legc_avg_ct <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, mean))

# delta_ipc_nsc <- setNames(a_legc_avg_ct[,'IPC'] - a_legc_avg_ct[,'NSC'], rownames(a_legc_avg_ct))

# # png(file.path(wd, 'output/sequence_modeling/figs/delta_atac_vs_delta_proximal_atac.png'), h = 500, w = 600)
# pdf(file.path(wd, 'output/sequence_modeling/figs/delta_atac_vs_delta_proximal_atac.pdf'), h = 500/71, w = 600/71)
# par(las = 2, mar = c(13, 5, 3,1), cex.lab = 1.5)
# boxplot_vec(xvec = atac_diff_vec[names(mcp_rnb)], yvec = delta_ipc_nsc[mcp_rnb], nm = 'Delta ATAC vs delta proximal ATAC', num_bins = 11, ylab = 'Delta ATAC')
# title(xlab = 'Delta proximal IPC-NSC ATAC UMIs (50kbp)', line = 11)
# text(1.5, 2, labels = paste0('cor = ', signif(cor(atac_diff_vec[names(mcp_rnb)], delta_ipc_nsc[mcp_rnb], method = 'spearman', use = 'pairwise.complete.obs'), 2)), cex = 1.5)
# dev.off()

# match_rn_nsc_pca <- match(c(neu_rn, nc_neu_rn, ipc_pc_rn, hi_ipc_lo_nsc_rn), nsc_pca$rowname)
# boxplot(atac_diff_vec[nsc_pca$peak_name[match_rn_nsc_pca]] ~ c(rep('neu', length(neu_rn)), rep('nc_neu', length(nc_neu_rn)), rep('ipc_pc', length(ipc_pc_rn)), rep('hi_ipc_lo_nsc', length(hi_ipc_lo_nsc_rn))), ylim = c(-2e+3, 4e+3), 
#         main = 'Delta proximal ATAC', col = rev(c('red', 'darkgray', col_key[['IPC']], 'orange3')))
# # boxplot(atac_diff_vec[nsc_pca$peak_name[match_rn_nsc_pca]] ~ c(rep('neu', length(neu_rn)), rep('nc_neu', length(nc_neu_rn))), ylim = c(-2e+3, 2e+3), xlim = c(0,6))

# match_rn_nsc_pca <- match(c(neu_rn, nc_neu_rn, ipc_pc_rn, hi_ipc_lo_nsc_rn), nsc_pca$rowname)
# boxplot(rna_diff_vec[nsc_pca$peak_name[match_rn_nsc_pca]] ~ c(rep('neu', length(neu_rn)), rep('nc_neu', length(nc_neu_rn)), rep('ipc_pc', length(ipc_pc_rn)), rep('hi_ipc_lo_nsc', length(hi_ipc_lo_nsc_rn))), ylim = c(-1.2,.4), 
#         main = 'Delta proximal log RNA', col = rev(c('red', 'darkgray', col_key[['IPC']], 'orange3')))
# # boxplot(atac_diff_vec[nsc_pca$peak_name[match_rn_nsc_pca]] ~ c(rep('neu', length(neu_rn)), rep('nc_neu', length(nc_neu_rn))), ylim = c(-2e+3, 2e+3), xlim = c(0,6))

# boxplot_vec(xvec = rna_diff_vec[names(mcp_rnb)], yvec = delta_ipc_nsc[mcp_rnb], nm = 'Delta ATAC vs delta proximal RNA', num_bins = 11)
# cor(rna_diff_vec[names(mcp_rnb)], delta_ipc_nsc[mcp_rnb], method = 'spearman')



# load('./output/sequence_modeling/enhflow_model_n_motif=16_w_clustering.rda')

# options(gmax.data.size = 1e+9)

# rename_prego_motifs <- function(model, tfs_in) {
#     amd <- prego::all_motif_datasets()
#     amou <- unique(amd[,c('motif', 'motif_orig')])
#     motif_clustering <- model@params$distilled_features
#     prego_motifs_renamed <- sapply(names(model@motif_models), function(cni) {
#         candidates <- motif_clustering$feat[motif_clustering$clust == cni]
#         candidates <- setNames(candidates, unlist(purrr::map(stringr::str_split(amou$motif_orig[match(candidates, amou$motif)], '_'), 1)))
#         candidates_all <- candidates[tolower(names(candidates)) %in% tolower(tfs_in)]
#         # candidates_feats <- candidates[tolower(names(candidates)) %in% tolower(intersect(tfs_in, feats))]
#         # if (length(candidates_feats) > 0) {
#         #     return(candidates_feats[[1]])
#         # } else
#         if (length(candidates_all) > 0) {
#             return(candidates_all[[1]])
#         } else {
#             return(candidates[[1]])
#         }
#     })
#     return(prego_motifs_renamed) 
# }

# prego_motifs <- do.call('rbind', lapply(seq_along(tm_w_add_feat@motif_models), function(i) dplyr::mutate(tm_w_add_feat@motif_models[[i]]$pssm, motif = names(tm_w_add_feat@motif_models)[[i]])))
# prego_ie <- prego::gextract_pwm(intervals = mcp, motifs = unique(prego_motifs$motif), dataset = prego_motifs)
# prego_ie_mat <- subset(prego_ie, select = -c(chrom, start, end, peak_name))
# rownames(prego_ie_mat) <- prego_ie$peak_name
# colnames(prego_ie_mat) <- prego_motifs_renamed[colnames(prego_ie_mat)]

# load('./output/sequence_modeling/dinucs_feat_peaks.rda')

# preds_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['contrib_meth']]))

# dtest_prego_ie_mat <- cbind(prego_ie_mat, 'NSC_ATAC' = mean(quantile(rowMeans(a_legc_by_day_n[rownames(preds_all),]), c(0,1))),'methylation' = 0, 
#                                         dinucs_per_peak[rownames(prego_ie_mat),intersect(colnames(dinucs_per_peak), colnames(preds_all))], 
#                                         'GC content' = 0.05, 'prox_ATAC' = 10, 'prox_RNA' = 10)

# pred_mcp <- setNames(predict(xgb_cv_res[[1]]$bstDMatrix_meth, xgb.DMatrix(as.matrix(dtest_prego_ie_mat)), predcontrib=FALSE, predinteraction = FALSE), rownames(dtest_prego_ie_mat))

# nei_orig_mcp <- gintervals.neighbors(dplyr::relocate(orig_lib_a, rowname, .after = end), as.data.frame(mcp), mindist = 0, maxdist = 0, maxneighbors = 1e+3)

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# nei_orig_peaks_f <- nei_orig_peaks[nei_orig_peaks$dist <= 5e+4,]

# match_rn_nei_mcp <- match(c(neu_rn, nc_neu_rn), nei_orig_mcp$rowname)
# boxplot(pred_mcp[nei_orig_mcp$peak_name[match_rn_nei_mcp]] ~ c(rep('neu', length(neu_rn)), rep('nc_neu', length(nc_neu_rn))), xlim = c(0,9), add = F, at = c(1,2))
# match_rn_nei_peaks_f <- match(c(neu_rn, nc_neu_rn), nei_orig_peaks_f$rowname)
# boxplot(pred_mcp[nei_orig_peaks_f$peak_name[match_rn_nei_peaks_f]] ~ c(rep('neu_nei', length(neu_rn)), rep('nc_neu_nei', length(nc_neu_rn))), add = T, at = c(3,4))

# match_rn_nei_peaks_f <- match(c(ipc_pc_rn, hi_ipc_lo_nsc_rn), nei_orig_peaks_f$rowname)
# boxplot(pred_mcp[nei_orig_peaks_f$peak_name[match_rn_nei_peaks_f]] ~ c(rep('IPC PC', length(ipc_pc_rn)), rep('hi_ipc_lo_nsc', length(hi_ipc_lo_nsc_rn))), add = T, at = c(5,6))
# match_rn_nei_peaks_f <- match(c(ipc_pc_rn, hi_ipc_lo_nsc_rn), nei_orig_mcp$rowname)
# boxplot(pred_mcp[nei_orig_mcp$peak_name[match_rn_nei_peaks_f]] ~ c(rep('IPC PC', length(ipc_pc_rn)), rep('hi_ipc_lo_nsc', length(hi_ipc_lo_nsc_rn))), add = T, at = c(7,8))

# colnames(prego_ie_mat)

# mtt <- c('HOMER.Olig2', 'HOCOMOCO.OLIG2_HUMAN.H11MO.0.B', 'JASPAR.NEUROG1')

# all_rn_list <- list(neu_rn, nc_neu_rn, ipc_pc_rn, hi_ipc_lo_nsc_rn)
# plot(0, col = 'white', xlim = c(0,14), ylim = c(-16,-13))
# ccc <- lapply(seq_along(all_rn_list), function(i) {
#     x <- all_rn_list[[i]]
#     pks_nei <- nei_orig_peaks_f$peak_name[match(x, nei_orig_peaks_f$rowname)]
#     print(i)
#     boxplot(prego_ie_mat[pks_nei,mtt], at = (i-1)+seq(1,length(all_rn_list)*length(mtt),length(all_rn_list)), add = T)
# })







# options(repr.plot.width = 15)
# options(repr.plot.height = 5)
