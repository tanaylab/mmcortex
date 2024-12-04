library(tgstat)
library(metacell)
library(misha)
library(misha.ext)
library(xgboost)
library(matrixStats)
# devtools::load_all("~/src/prego")
# devtools::load_all('~/src/iceqream/')

wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')
source('./scripts/util.r')
options(gmax.data.size = 1e+9)
scdb_init('scdb', f=T)

mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
mcmd <- mcmd[-c(602:603),]

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))

cust_st_ord2 <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN/CfuPN','iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st2 <- unlist(lapply(cust_st_ord2, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))

cts <- c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')

color_key <- unique(mcmd[,c('cell_type', 'color')])
color_key <- color_key[match(cust_st_ord, color_key$cell_type),]
col_key <- tibble::deframe(color_key)
col_annot <- as.data.frame(mcmd[,c('cell_type', 'mean_day')])
rownames(col_annot) <- mcmd$metacell

ann_colors <- list(cell_type = tibble::deframe(color_key),
                   mean_day = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green3', 'blue', 'purple'))(100), 
                                      seq(13,18,l=100)))

clrmp <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue1', 'blue4', 'purple3'))(1000)

clrmp_abs <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(1000)
brks_abs <- seq(-16.6,-10, l=length(clrmp_abs))

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=length(clrmp_rel))                             

nsc_mcs <- mcmd$metacell[mcmd$cell_type == 'NSC']
ipc_mcs <- mcmd$metacell[mcmd$cell_type %in%  c('IPC', 'IPC_cyc')]

mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds')

mca@peaks$peak_name_ntb <- mcATAC::peak_names(subset(mca@peaks, select = -peak_name), tad_based = F)

mc_rna <- scdb_mc('pl_cort')
mat_rna <- scdb_mat('pl_cort')

mca@rna_egc <- mc_rna@e_gc

tss <- gintervals.load('intervs.global.tss')

tss <- tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]

tss2 <- tss[!duplicated(tss$geneSymbol),]

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

tads <- gintervals.load('intervs.global.tad_names')

legc <- log2(1e-5 + mca@rna_egc)

a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(mca@rna_egc))
head(tfs)

tfs_in <- intersect(tfs, rownames(legc)[which(rowMaxs(legc) >= -13.5)])

# agg_id <- readr::read_csv('./scatac_data//aggregation_id.csv')
# mca@cell_to_metacell$agg_id <- as.numeric(unlist(purrr::map(stringr::str_split(mca@cell_to_metacell$cell_id, '-'), 2)))
# agg_id_day <- data.frame(cbind(1:14, rep(12:18, each = 2)))
# colnames(agg_id_day) <- c('agg_id', 'day')
# mca@cell_to_metacell$day <- agg_id_day$day[mca@cell_to_metacell$agg_id]

# atac_mc_day <- table(mca@cell_to_metacell$metacell, mca@cell_to_metacell$day)
# atac_mc_day_norm <- atac_mc_day/rowSums(atac_mc_day)
# rna_mc_day <- as.matrix(mcmd[,grep('E\\d\\d', colnames(mcmd))])
# rna_mc_day_norm <- rna_mc_day/rowSums(rna_mc_day)
# nsc_mcs <- which(mcmd$cell_type == 'NSC')
# ipc_mcs <- which(mcmd$cell_type %in% c('IPC', 'IPC_cyc'))
# egc_by_day <- mca@egc[,nsc_mcs] %*% atac_mc_day_norm[nsc_mcs,]
# egc_by_day_n <- t(t(egc_by_day)/colSums(egc_by_day))
# colnames(egc_by_day_n) <- paste0('E', colnames(egc_by_day_n))
# egc_by_day_ipc <- mca@egc[,ipc_mcs] %*% atac_mc_day_norm[ipc_mcs,]
# egc_by_day_n_ipc <- t(t(egc_by_day_ipc)/colSums(egc_by_day_ipc))
# colnames(egc_by_day_n_ipc) <- paste0('E', colnames(egc_by_day_n_ipc))

# a_legc_by_day_n <- log2(1e-5 + egc_by_day_n)

load('./output/metacell_model/cell_cycle_phase_data.rda')

load('./output/metacell_model/pcurve_nsc_vs_mat_neuro.rda')
ro <- rev(pcu$ord)

load('./output/mcatac/var_peaks_after_enh_prom_separation.rda')

# a_legc_avg_cl <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
a_legc_avg_ct <- t(tgs_matrix_tapply(a_legc[dist_peaks$peak_name,], mcmd$cell_type, mean))
# a_legc_avg_cl_ct <- t(tgs_matrix_tapply(tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean), mcmd$cell_type, mean))


load(file='./output/mcatac/new_nsc_ipc_transition_peak_clusters_to_plot.rda')



load('./output/mcatac/var_peaks_after_enh_prom_separation.rda')


load('./output/methylation//avg_meth_all.rda')

rm_ama <- rowMeans(avg_meth_all[,grep('E\\d\\d', colnames(avg_meth_all))], na.rm = T)


print('check 1')

pca_nsc_ipc <- proximal_chromatin_activity(peaks_of_interest = mcp, background_peaks = mcp, restrict_to_tads = F,
                                           d_puncture = 1e+3, d_proximity_atac =5e+4, d_proximity_rna = 5e+5, 
                                           mc_sel = nsc_mcs, mat_rna = mat_rna, mc_rna = mc_rna, tss = tss, tads = tads, mca = mca)

rna_umi_vec <- log2(1+tibble::deframe(pca_nsc_ipc[,c('peak_name', 'rna_umi')]))

atac_umi_vec <- log2(1+tibble::deframe(pca_nsc_ipc[,c('peak_name', 'atac_umi')]))

load('./output/sequence_modeling/nuc_feat_peaks.rda')

scp <- sort(cln_plot)

load('./output/sequence_modeling/dinucs_feat_peaks.rda')

# n_5mCpG <- 300*dinucs_per_peak[,'CG']*rm_ama[rownames(dinucs_per_peak)]
# n_nmCpG <- 300*dinucs_per_peak[,'CG']*(1-rm_ama[rownames(dinucs_per_peak)])

# load('./output/sequence_modeling/enhflow_output_delta_ipc_nsc_n_motifs=8.rda')

# load('./output/sequence_modeling/enhflow_model_w_add_feat_120923.rda')

# load('./output/sequence_modeling/enhflow_model_50_motifs_w_hc.rda')
load('./output/sequence_modeling/enhflow_model_n_motif=16_w_clustering.rda')

# load('./output/sequence_modeling/enhflow_model_astro_nsc_16_motifs.rda')

# pssms <- lapply(tm_w_nsc_e13_17@motif_models, function(x) x$pssm)

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

prego_motifs_renamed <- rename_prego_motifs(tm_w_add_feat, tfs_in = tfs_in)
# prego_motifs_renamed

prego_motifs_manual <- setNames(c('E_box_1', 'T_box_1', 'GATA_1', 'Fos_Jun_1', 'RFX_1', 'AAAG_1', 'AAAG_2', 'SOX', 'E_box_2', 'E_box_3', 'NR2', 'GACA', 'E_box_4',  'E_box_5', 'E_box_6', 'NR1'), names(prego_motifs_renamed))
# prego_motifs_manual

prego_motifs <- do.call('rbind', lapply(seq_along(tm_w_add_feat@motif_models), function(i) dplyr::mutate(tm_w_add_feat@motif_models[[i]]$pssm, motif = names(tm_w_add_feat@motif_models)[[i]])))
prego_ie <- prego::gextract_pwm(intervals = mcp, motifs = unique(prego_motifs$motif), dataset = prego_motifs)

prego_ie_mat <- subset(prego_ie, select = -c(chrom, start, end, peak_name))
rownames(prego_ie_mat) <- prego_ie$peak_name
colnames(prego_ie_mat) <- prego_motifs_renamed[colnames(prego_ie_mat)]
# motif_clustering <- tm_w_add_feat@params$distilled_features

colnames(prego_ie_mat) <- prego_motifs_manual[match(colnames(prego_ie_mat), prego_motifs_renamed)]

# amd <- prego::all_motif_datasets()

# amou <- unique(amd[,c('motif', 'motif_orig')])

# feats <- names(scdb_gset('pl_cort_dns')@gene_set)


# km_prom_a_legc <- tglkmeans::TGL_kmeans(a_legc[prom_peaks$peak_name,], k = 10, seed = 1337)

# km_enh_a_legc <- tglkmeans::TGL_kmeans(a_legc[dist_peaks$peak_name,], k = 60, seed = 1337)

# a_legc_avg_cl_enh <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
# a_legc_avg_cl_prom <- tgs_matrix_tapply(t(a_legc[prom_peaks$peak_name,]), km_prom_a_legc$cluster, mean)


xnsc <- a_legc_avg_ct[,'NSC']
xipc <- a_legc_avg_ct[,'IPC']
delta_ipc_nsc <- xipc - xnsc
y_all <- delta_ipc_nsc

pb <- intersect(setdiff(names(rm_ama), prom_peaks$peak_name), var_peaks)
length(pb)

## True variables
x_all <- cbind(prego_ie_mat[pb,], xnsc[pb], rm_ama[pb], 
                    dinucs_per_peak[pb,c('AA','AC','AG','AT','CA','CC','CG','CT')], 
                    nuc_per_peak[pb,'C'], 
                    atac_umi_vec[pb], 
                    rna_umi_vec[pb])
colnames(x_all) <- c(colnames(prego_ie_mat),'NSC_ATAC','methylation', 
                        c('AA','AC','AG','AT','CA','CC','CG','CT'), 
                    'GC content', 'prox_ATAC', 'prox_RNA')



x_all <- as.data.frame(x_all[apply(x_all, 1, function(x) all(!is.na(x))),])

test_peaks <- rownames(x_all)[sample(1:nrow(x_all), round(nrow(x_all)/20))]
train_peaks <- rownames(x_all)[!(rownames(x_all) %in% test_peaks)]
val_peaks <- sample(train_peaks, round(nrow(x_all)/10))
train_peaks <- train_peaks[!(train_peaks %in% val_peaks)]

x_train = x_all[train_peaks,]
x_test = x_all[val_peaks,]
y_train = delta_ipc_nsc[train_peaks]
y_test = delta_ipc_nsc[val_peaks]


eval_results <- function(true, predicted) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  SSE_nna <- sum((predicted - true)^2, na.rm = T)
  SST_nna <- sum((true - mean(true))^2, na.rm = T)
    
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/length(SSE))
  # Model performance metrics
  return(data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  ))
}

run_xgb_cv <- function(x_all, y_all, num_peak_groups = 3) {
    peak_groups <- setNames(sample(1:num_peak_groups, nrow(x_all), replace = T), rownames(x_all))
    pga <- sort(unique(peak_groups))
    test_grp <- tail(pga, 1)
    train_grps <- setdiff(pga, test_grp)
    xgb_res <- lapply(train_grps, function(ui) {
        tgh <- setdiff(train_grps, ui)
        train_peaks <- names(peak_groups[peak_groups %in% tgh])
        val_peaks <- names(peak_groups[peak_groups == ui])
        x_train = x_all[train_peaks,]
        x_test = x_all[val_peaks,]
        y_train = y_all[train_peaks]
        y_test = y_all[val_peaks]
        dtrain <- xgb.DMatrix(data = x_train, label = y_train)
        dtest <- xgb.DMatrix(data = x_test, label = y_test)
        xgb_params = list('md' = 3, 'eta' = 0.3, 'nr' = 250, obj = "reg:squarederror", em = 'rmse')
        bstDMatrix_meth <- xgboost(data = dtrain, max.depth = xgb_params$md, verbose = 0,
                          eta = xgb_params$eta, 
                          nthread = 50, 
                          nrounds = xgb_params$nr, 
                          objective = xgb_params$obj,
                          eval_metric = xgb_params$em,
                                  )
        pred_meth <- predict(bstDMatrix_meth, dtest, predcontrib=FALSE, predinteraction = FALSE)
        names(pred_meth) <- val_peaks
        
        contrib_meth <- predict(bstDMatrix_meth, dtest, predcontrib=TRUE, predinteraction = FALSE)
        rownames(contrib_meth) <- val_peaks
        ev_xgb_meth = eval_results(y_test, pred_meth)
        return(list(ev_xgb_meth = ev_xgb_meth, pred_meth = pred_meth, y_test = y_test, 
                    bstDMatrix_meth = bstDMatrix_meth, contrib_meth = contrib_meth))
    })
    return(xgb_res)
}

LL <- -16
UL <- -15
inc_peaks <- intersect(names(rm_ama)[which(rm_ama > 0.5)], rownames(x_all))
mid_peaks <- intersect(names(rm_ama)[which(rm_ama > 0.1 & rm_ama <= 0.5)], rownames(x_all))
dec_peaks <- intersect(names(rm_ama)[which(rm_ama <= 0.1)], rownames(x_all))

xgb_cv_res <- run_xgb_cv(as.matrix(x_all), y_all, num_peak_groups = 11)

ev_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['ev_xgb_meth']]))
ev_all

preds_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['contrib_meth']]))

pred_meth_all <- do.call('c', lapply(xgb_cv_res, function(x) x[['pred_meth']]))


# png('./output/sequence_modeling/figs/xgboost_predicted_vs_observed_delta_atac_peak_clust_var_new.png', h = 500, w = 600)



fs <- sort(apply(preds_all, 2, sd), decreasing = T)
fs <- fs[names(fs) != 'BIAS']
feats_sort <- factor(names(fs), levels = names(fs))


pvl <- tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(subset(preds_all,select=-BIAS))), 
                cols = colnames(subset(preds_all, select = -BIAS)))

pvl$lvl <- feats_sort[match(pvl$name, levels(feats_sort))]

x_all_long <- tidyr::pivot_longer(tibble::rownames_to_column(x_all), cols = colnames(x_all))

mtch_inds_x_all_preds <- match(apply(pvl[,1:2], 1, paste0, collapse = '-'), apply(x_all_long[,1:2], 1, paste0, collapse = '-'))

pvl$feature_value <- x_all_long$value[mtch_inds_x_all_preds] 

pvl <- pvl %>% dplyr::group_by(name) %>% 
            dplyr::mutate(val_lin = (feature_value  - min(feature_value))/(max(feature_value)-min(feature_value)),
                                                color = clrmp_rel[1+round(length(clrmp_rel)*val_lin)])

integral_over_feature_contrib <- colSums(abs(subset(preds_all, select = -BIAS)))
norm_integral_over_feature_contrib <- integral_over_feature_contrib/sum(integral_over_feature_contrib)

samp_df <- dplyr::sample_frac(pvl, 0.1)


meth_fct <- setNames(factor(c(rep('meth_low', length(dec_peaks)), rep('meth_mid',length(mid_peaks)),
                        rep('meth_high', length(inc_peaks))), 
                        levels = c('meth_low', 'meth_mid', 'meth_high')), 
                        c(dec_peaks, mid_peaks, inc_peaks))

NSC_0_THRESH <- -16.4
meth_fct_filt_nsc <- meth_fct[xnsc[names(meth_fct)] <= NSC_0_THRESH]

olig2_vec <- setNames(x_all[names(meth_fct_filt_nsc),'E_box_1'], names(meth_fct_filt_nsc))
olig2_cut <- setNames(cut(olig2_vec, quantile(olig2_vec, (0:3)/3), labels = c('E_box_1_low', 'E_box_1_mid','E_box_1_high')), names(meth_fct_filt_nsc))


feature_cut <- function(feature_vec, cells, name, cut_values = NULL, cut_quantiles = NULL) {
    if (is.null(names(feature_vec))) {
        stop('feature_vec should be a named vector')
    }
    if (!is.null(cut_values)) {
        feature_vec_cut <- setNames(cut(feature_vec, breaks = cut_values, labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
                               names(feature_vec))
    } else if (!is.null(cut_quantiles)) {
        feature_vec_cut <- setNames(cut(feature_vec, breaks = quantile(feature_vec, cut_quantiles), labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
                               names(feature_vec))
    }
    else {
        feature_vec_cut <- setNames(cut(feature_vec, breaks = c(quantile(feature_vec, 0)-1e-3, quantile(feature_vec, (1:3)/3)), labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
                               names(feature_vec))
    }
    # print(length(which(is.na(feature_vec_cut))))
    return(feature_vec_cut)
}

meth_cells <- pb
meth_cut <- feature_cut(feature_vec = rm_ama[meth_cells], cells = meth_cells, cut_values = c(-1e-2,0.25,0.75,1), name = 'meth')
prox_atac_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'prox_ATAC'], names(meth_cut)), cells = names(meth_cut), name = 'prox_ATAC')
prox_rna_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'prox_RNA'], names(meth_cut)), cells = names(meth_cut), name = 'prox_RNA')
olig2_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'E_box_1'], names(meth_cut)), cells = names(meth_cut), name = 'E_box_1')
sox_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'SOX'], names(meth_cut)), cells = names(meth_cut), name = 'SOX')
eomes_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'T_box_1'], names(meth_cut)), cells = names(meth_cut), name = 'T_box_1')
nsc_atac_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'NSC_ATAC'], names(meth_cut)), cells = names(meth_cut), name = 'NSC_ATAC')

feat_cut_df <- as.data.frame(meth_cut)
feat_cut_df$prox_atac_cut <- prox_atac_cut
feat_cut_df$prox_rna_cut <- prox_rna_cut
feat_cut_df$olig2_cut <- olig2_cut
feat_cut_df$sox_cut <- sox_cut
feat_cut_df$eomes_cut <- eomes_cut
feat_cut_df$nsc_atac_cut <- nsc_atac_cut
rownames(feat_cut_df) <- names(meth_cut)




save(delta_ipc_nsc, 
            rm_ama,
            pred_meth_all, 
            ev_all,
            samp_df,
            norm_integral_over_feature_contrib,
            olig2_vec,
            feat_cut_df,
            file = './output/sequence_modeling/fig_6_data.rda')
            
pca_nsc_ipc <- proximal_chromatin_activity(peaks_of_interest = mcp, background_peaks = mcp, restrict_to_tads = F,
                                           d_puncture = 1e+3, d_proximity_atac =5e+4, d_proximity_rna = 5e+5, 
                                           mc_sel = as.character(nsc_mcs), mat_rna = mat_rna, mc_rna = mc_rna, tss = tss, 
                                           tads = tads, mca = mca)

nsc_pca <- pca_nsc_ipc
ipc_pca <- proximal_chromatin_activity(peaks_of_interest = mcp, background_peaks = mcp, restrict_to_tads = F,
                                           d_puncture = 1e+3, d_proximity_atac =5e+4, d_proximity_rna = 5e+5, 
                                           mc_sel = as.character(ipc_mcs), mat_rna = mat_rna, mc_rna = mc_rna, tss = tss, 
                                           tads = tads, mca = mca)

rnb <- intersect(nsc_pca$peak_name, ipc_pca$peak_name)
atac_diff_vec <- setNames(unlist(ipc_pca[match(rnb, ipc_pca$peak_name),'atac_umi']) - unlist(nsc_pca[match(rnb, nsc_pca$peak_name),'atac_umi']), rnb)

y_all_01 <- (y_all - min(y_all))/(max(y_all) - min(y_all))
y_all_clvls <- setNames(clrmp_rel[1+round((length(clrmp_rel)-1)*y_all_01)], names(y_all))

save(tm_w_add_feat,
        prego_motifs_manual,
        x_all,
        preds_all,
        delta_ipc_nsc, 
        y_all_clvls,
        rnb, 
        atac_diff_vec,
        file = './output/sequence_modeling/fig_s6_data.rda')






# intervs_energy_new <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_prego_motif_energy.rds')
# ie_mat <- as.matrix(subset(intervs_energy_new, select = -c(chrom, start, end, peak_name, mmcortex.marginal,	intervalID,	peak_name, peak_name_ntb)))
# rownames(ie_mat) <- intervs_energy_new$peak_name

# raq98 <- unlist(plyr::llply(colnames(ie_mat), function(x) quantile(ie_mat[,x], probs = 0.98), .parallel = T))

# names(raq98) <- colnames(ie_mat)

# ra_98_bin_int <- t(plyr::laply(1:length(raq98), function(i) as.numeric(ie_mat[,i] >= raq98[[i]]), .parallel = T))
# save(ra_98_bin_int, file='./output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda')
# colnames(ra_98_bin_int) <- colnames(ie_mat)
# rownames(ra_98_bin_int) <- rownames(ie_mat)

# load('./output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda')
# ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int[dist_peaks$peak_name,])), km_enh_a_legc$cluster, sum)
# ra_98_lfc_int <- log2(1e-2+ra_98_sum_clust_int/(0.02*as.numeric(km_enh_a_legc$size)))
# colnames(ra_98_lfc_int) <- colnames(ra_98_bin_int)

# mtt <- c("JASPAR.EOMES", 
#          "JOLMA.MEIS2_mono_DBD_2",
#          'JASPAR.NEUROG1',
#          "JASPAR.NEUROD1",
#          'HOMER.Sox2',
#          "JASPAR.EMX1",
#          'JASPAR.POU3F2',
#          'JASPAR.NFIA',
#          "JASPAR.MEF2C",
#          "JASPAR.FOXP1",
#          "HOCOMOCO.MECP2_MOUSE.H11MO.0.C",
#          'HOMER.CTCF',
#          'HOMER.NRF1',
#          'HOCOMOCO.KLF3_MOUSE.H11MO.0.A',
#          'JOLMA.ETV1_mono_DBD',
#          'JASPAR.NFIB'
#         )

# options(repr.plot.width =10)
# options(repr.plot.height =16)

# p_motifs_lfc_mtt <- pheatmap::pheatmap(round(ra_98_lfc_int[enh_cl_ord,mtt], 3), breaks = seq(-3,3, l = length(clrmp_rel)), col = clrmp_rel, cluster_rows = F, fontsize = 12, treeheight_col = 0)

# mcATAC::save_pheatmap(p_motifs_lfc_mtt, './output/mcatac/figs/motifs_annot_mat_mtt.png', h = 1800, w = 850)

# # New version - 25/04/24
# inds_glia <- which(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes'))
# inds_no_glia <- which(!(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes')))


# pltmt <- a_legc_avg_cl_enh[enh_cl_ord, cust_mc_ord_st2]

# brks <- seq(-16.7,-14.75,l=100)

# col_ha1 <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_no_glia],'cell_type'], 
#                                                      col = ann_colors[['cell_type']],
#                                                      height =unit(2, 'cm')), 
#                              mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_no_glia],'mean_day'], axis_param = list(gp = gpar(fontsize = 0)),
#                                                     # col = circlize::colorRamp2(seq(13,18,1), c('red', 'orange', 'yellow', 'green2', 'blue', 'purple')), 
#                                                     height =unit(2, 'cm')), 
#                              annotation_name_gp = gpar(fontsize = 18),
#                              show_legend = F)
# # row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_a_legc$size[as.numeric(hcct2)]), ylim = c(800,2500), gp = gpar(fill = 'black',fontsize = 18), 
# row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_enh_a_legc$size[as.numeric(enh_cl_ord)]), ylim = c(800,8000), gp = gpar(fill = 'black',fontsize = 18), 
#                                                         axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90)), 
#                             # frac_prom = anno_barplot(as.numeric(t(perc_prom[peak_clust_var[hcct2]])), ylim = c(0,1), 
#                             # gp = gpar(fill = 'black',fontsize = 18), 
#                             # axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90),
#                             annotation_name_gp = gpar(fontsize = 18),
#                             annotation_name_rot = 0,
#                             # cl2 = anno_simple(x = ct_hc_a_legc[ord1], col = setNames(peak_cl2_color_key$color[match(ct_hc_a_legc[ord1], peak_cl2_color_key$cl2)], ct_hc_a_legc[ord1])),
#                             annotation_name_offset = unit(3, 'cm'),        
#                             which = 'row',
#                             width = unit(5, 'cm')
#                            )
# clrmp_rel2 <- circlize::colorRamp2(breaks = c(-3,0,3), colors = c('blue3','white','red3'))
# motifs_anno_mat <- round(ra_98_lfc_int[enh_cl_ord,mtt], 3)
# # motifs_anno_mat <- round(ra_98_lfc_int[hcct2,mtt], 3)
# # motifs_anno_mat <- round(prego_98_lfc[hcct2,], 3)
# hc_mtt <- hclust(dist(t(motifs_anno_mat)), method = 'ward.D2')
# colnames(motifs_anno_mat) <- unlist(purrr::map(stringr::str_split(colnames(motifs_anno_mat), '\\.'), 2))
# colnames(motifs_anno_mat) <- unlist(purrr::map(stringr::str_split(colnames(motifs_anno_mat), '_'), 1))
# mam_lin <- apply(motifs_anno_mat, 2, function(x) {y <- x; y[x < -3] <- -3; y[x > 3] <- 3; return((y + 3)/6)})

# # motif_ha <- rowAnnotation(SOX = anno_numeric(motifs_anno_mat[,grep('Sox2', colnames(motifs_anno_mat), ign = T)], 
# #                                                  labels_gp = gpar('fontsize' = 0),
# #                                                     # bg_gp = gpar(fill = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('sox2', colnames(mam_lin), ign = T)])], each = 2), alpha.f = 0.75), ncol = 2))),
# #                                              # bg_gp = gpar('fill' = head(matrix(bg_gp = gpar('fill' = head(matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('sox2', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2), 1))),, ncol = 2), 1))),
# #                                              bg_gp = gpar(fill = head(clrmp_rel[1+round(999*mam_lin[,grep('sox2', colnames(mam_lin), ign = T)])], -1))),
# #                               NEUROD1 = anno_numeric(motifs_anno_mat[,grep('NEUROD1', colnames(motifs_anno_mat), ign = T)], 
# #                                                      labels_gp = gpar('fontsize' = 0),
# #                                                     bg_gp = gpar('fill' = head(matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('neurod1', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2), 1))),
# #                               EOMES = anno_numeric(motifs_anno_mat[,grep('eomes', colnames(motifs_anno_mat), ign = T)], 
# #                                                     labels_gp = gpar('fontsize' = 0),
# #                                                     bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('eomes', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
# #                               POU3F2 = anno_numeric(x = motifs_anno_mat[,grep('pou3f2', colnames(motifs_anno_mat), ign = T)], 
# #                                                     labels_gp = gpar('fontsize' = 0),
# #                                                     bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('pou3f2', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
# #                               # which = 'row',
# #                             width = unit(9, 'cm')
# #                              )

# motif_ha <- rowAnnotation(`motif\nenrichment` = motifs_anno_mat[,hc_mtt$order], col = list(`motif\nenrichment` = clrmp_rel2),
#                             width = unit(9, 'cm')
#                              )



# ch <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_no_glia], 
#                               name = 'log2\nfraction\nATAC',
#                               # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
#                               col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), 
#                                 colors = c('white', 'orange', 'red', 'purple', 'black')),
#                               column_split = factor(names(cust_mc_ord_st2[inds_no_glia]), levels = cust_st_ord2), 
#                               column_gap = unit(2, 'mm'),
#                               column_title_gp = gpar(fontsize = 20),
#                               column_title_rot = 90,
#                               row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), 
#                               row_gap = unit(2, 'mm'),
#                                 row_title_gp = gpar(fontsize = 0),
#                               top_annotation = col_ha1, 
#                               # bottom_annotation = col_ha2,
#                               show_heatmap_legend = T,
#                               show_column_names = F,
#                               right_annotation = motif_ha,
#                               heatmap_legend_param = list(legend_height = unit(5, 'in'), legend_width = unit(5, 'in'), labels_gp = gpar(fontsize = 16)),
#                                 row_names_gp = gpar(fontsize = 14),
#                               # heatmap_width = unit(92, 'cm'), heatmap_height = unit(50, 'cm'),
#                               heatmap_width = unit(45, 'cm'), 
#                               heatmap_height = unit(32, 'cm'),
#                         # left_annotation = row_ha,
#                             cluster_columns = F, cluster_rows = F)

# options(repr.plot.width =20)
# options(repr.plot.height =16)

# draw(ch)

# col_ha1_glia <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_glia],'cell_type'], 
#                                                      col = ann_colors[['cell_type']],
#                                                      height =unit(2, 'cm')), 
#                              mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_glia],'mean_day'], axis_param = list(gp = gpar(fontsize = 16)), ylim = c(13,18),
#                                                     height =unit(2, 'cm')), 
#                              annotation_name_gp = gpar(fontsize = 0),
#                              show_legend = F)

# # ch_glia <- ComplexHeatmap::Heatmap(matrix = pltmt[hcct2,inds_glia], 
# ch_glia <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_glia], 
#                                 # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
#                               col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
#                                    column_split = factor(names(cust_mc_ord_st2[inds_glia]), levels = c('OPCs', 'Astrocytes')), column_gap = unit(3, 'mm'),
#                               column_title_gp = gpar(fontsize = 20),
#                               column_title_rot = 90,
#                               row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), 
#                                    row_gap = unit(2, 'mm'),
#                               row_title_gp = gpar(fontsize = 0),
#                               top_annotation = col_ha1_glia, 
#                                    left_annotation = row_ha,
#                               show_heatmap_legend = F,
#                               show_column_names = F,
#                                    show_row_names = F,
#                               heatmap_width = unit(10, 'cm'), 
#                                    heatmap_height = unit(32, 'cm'),
#                               row_names_gp = gpar(fontsize = 14),
#                         cluster_columns = F, cluster_rows = F)

# options(repr.plot.width = 14)
# options(repr.plot.height = 14)

# # draw(ch)

# # png('./output/mcatac/figs/mmcortex_famc_legc_cluster_size_frac_prom_annot_arb_met_0.45.png', w = 3200, h =1600)
# # png('./output/mcatac/figs/mmcortex_famc_legc_cluster_size_frac_prom_annot_left_arb_met_0.45.png', w = 2000, h =1600)
# # png('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc.png', w = 3000, h =1700)

# png('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc_no_glia.png', w = 1380, h =1170)
# draw(ch)
# dev.off()

# png('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc_glia.png', w = 380, h =1270)
# draw(ch_glia)
# dev.off()

# pdf('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc_no_glia.pdf', w = 1480/.75e+2, h =1170/.71e+2)
# draw(ch)
# dev.off()

# pdf('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc_glia.pdf', w = 320/.71e+2, h =1270/.71e+2)
# draw(ch_glia)
# dev.off()

# # cfupn clusters - 11, 19, 21, 22
# # cpn clusters - 23,25

# options(repr.plot.width =8)
# options(repr.plot.height =6)

# dir.create('./output/mcatac/figs/cpn_cfupn_clust_phm/')

# mcATAC::save_pheatmap

# plot_clust_pheatmap <- function(clj) {
#     pclj <- pheatmap::pheatmap(a_legc[dist_peaks$peak_name[km_enh_a_legc$cluster == clj],cust_mc_ord_st2], 
#                                fontsize = 16, treeheight_row = 0,
#                                main = paste0('cluster ', clj),
#                                show_rownames = F, show_colnames = F,
#                    annotation_col = col_annot, silent = T, color = clrmp_abs, annotation_colors = ann_colors, 
#                                cluster_cols =F, cluster_rows = T, clustering_method = 'ward.D2')
#     mcATAC::save_pheatmap(x = pclj, filename = paste0('./output/mcatac/figs/cpn_cfupn_clust_phm/', clj, '.png'),
#                          height = 1200, width = 1600)
# }
# ttt <- sapply(c(11, 12, 21, 22, 23, 24, 25), function(clj) plot_clust_pheatmap(clj))





# load(file = './output/methylation/avg_meth_all.rda')

# nsc_inc_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] >= 1)], rownames(avg_meth_all))
# nsc_inc_atac <- a_legc_by_day_n[nsc_inc_peaks,]
# nsc_dec_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] <= -1)], rownames(avg_meth_all))
# nsc_dec_atac <- a_legc_by_day_n[nsc_dec_peaks,]

# par(mfrow = c(1,2))
# # ord <- order(avg_meth_all[nsc_inc_peaks,'E17'] - avg_meth_all[nsc_inc_peaks,'E13'])
# ord <- order(rowMeans(a_legc_by_day_n[nsc_inc_peaks,]))
# p_inc_atac <- pheatmap::pheatmap(a_legc_by_day_n[nsc_inc_peaks[ord],], fontsize_col = 16,
#                                       # main = 'NSC increasing peaks',
#                                       fontsize = 16,
#                                       cluster_rows = F, 
#                                       clustering_method = 'ward.D2',
#                                       cluster_cols = F,
#                                       show_rownames = F,
#                                       treeheight_row = 0,
#                                       color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
#                                      breaks = seq(-16.6,-14,l = 1000))

# mcATAC::save_pheatmap(p_inc_atac, './output/mcatac//figs/phm_nsc_inc_atac.png', h = 600, w = 500)
# # ord <- order(avg_meth_all[nsc_dec_peaks,'E17'] - avg_meth_all[nsc_dec_peaks,'E13'], decreasing = T)

# ord <- order(rowMeans(a_legc_by_day_n[nsc_dec_peaks,]))
# p_dec_atac <- pheatmap::pheatmap(a_legc_by_day_n[nsc_dec_peaks[ord],], fontsize_col = 16,
#                                       # main = 'NSC decreasing peaks',
#                                       fontsize = 16,
#                                       cluster_rows = F, 
#                                       clustering_method = 'ward.D2',
#                                       cluster_cols = F,
#                                       show_rownames = F,
#                                       treeheight_row = 0,
#                                       color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
#                                      breaks = seq(-16.6,-14,l = 1000))

# mcATAC::save_pheatmap(p_dec_atac, './output/mcatac//figs/phm_nsc_dec_atac.png', h = 600, w = 500)

# source('./scripts/util.r')

# save_pheatmap_pdf(p_dec_atac, './output/mcatac//figs/phm_nsc_dec_atac.pdf', h = 600/91, w = 500/91)
# save_pheatmap_pdf(p_inc_atac, './output/mcatac//figs/phm_nsc_inc_atac.pdf', h = 600/91, w = 500/91)

# mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

# mpra_lib <- readr::read_tsv('./cpg_methylation//st_and_temporal_and_e14_shadow_enh_9-1-22.tsv')
# mpra_lib <- as.data.frame(dplyr::relocate(mpra_lib, rowname, .after = end))

# nei_mpra_mcp <- tidyr::drop_na(gintervals.neighbors1(mpra_lib, mcp, maxdist = 0, mindist = 0, maxneighbors = 3))

# nsc_asc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_asc')
# nsc_asc_atac <- a_legc_by_day_n[nsc_asc_nei$peak_name,]

# nsc_desc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_desc')
# nsc_desc_atac <- a_legc_by_day_n[nsc_desc_nei$peak_name,]


# vasc <- colMeans(a_legc[nsc_asc_nei$peak_name,])
# vdesc <- colMeans(a_legc[nsc_desc_nei$peak_name,])

# legc <- log2(1e-5 + mc_rna@e_gc)

# pltmt <- rbind(vasc,
#                vdesc,
#                legc[sort(grep('os|Sox10|Sox13|Sox7|Sox17|Sox15|Sox18|Sox30',
#                          # c(
#                              grep('Tet|Dnmt',rownames(mc_rna@e_gc), v=T) 
#                            # sort(grep('Eomes|Pou3f2|Neurog2|Sox|Hox', rownames(mc_rna@e_gc), v=T)))
#                         , inv=T, v=T)),mcmd$metacell]
#                )

# rownames(pltmt) <- c('inc_peaks', 'dec_peaks', tail(rownames(pltmt), -2))

# par(mfrow = c(2,3))
# cells_h <- as.character(mcmd$metacell[which(mcmd$cell_type == 'NSC')])
# lm_and_plot <- function(x, y, xlab = NULL, ylab = NULL, main = NULL) {
#     lm1 <- lm(y ~ x)
#     plot(x,y, xlab = xlab, ylab = ylab, main = glue::glue("cor = {round(cor(x,y,method = 'pearson'), 3)}, pval = {round(cor.test(x,y,method ='pearson')$p.value, 3)}"))
#     abline(lm1$coefficients[[1]], lm1$coefficients[[2]], col= 'red', lty =2)
# }

# png('./output/methylation/figs/scatter_cor_nsc_inc_dec_peaks_vs_meth_components.png', h = 600, w = 1000, res = 100)
# par(mfrow = c(2,3), cex.lab = 1.5, mar = c(5,5,2,1), cex.axis = 1.5, cex.main = 1.5)
# lm_and_plot(pltmt['Dnmt1',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt1', ylab = 'NSC dec_peaks ATAC')
# lm_and_plot(pltmt['Dnmt3a',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3a', ylab = 'NSC dec_peaks ATAC')
# lm_and_plot(pltmt['Dnmt3b',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3b', ylab = 'NSC dec_peaks ATAC')
# lm_and_plot(pltmt['Tet1',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet1', ylab = 'NSC inc_peaks ATAC')
# lm_and_plot(pltmt['Tet2',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet2', ylab = 'NSC inc_peaks ATAC')
# lm_and_plot(pltmt['Tet3',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet3', ylab = 'NSC inc_peaks ATAC')
# dev.off()



# options(repr.plot.width = 6)
# options(repr.plot.height = 36)

# boxplot_vec <- function(xvec, yvec, nm, num_bins = 7, ylab = '', xlab = '') {
#     xvec_rng <- c(min(xvec, na.rm = T),max(xvec, na.rm = T))
#     cvc <- seq(xvec_rng[[1]], xvec_rng[[2]], l = num_bins)
#     mtfc <- cut(xvec, breaks = cvc)
#     boxplot(yvec ~ mtfc, main = nm, ylab = ylab, xlab  = xlab)
#     text(1:length(levels(mtfc)), rep(quantile(yvec, 0.99), length(levels(mtfc))), labels = table(mtfc), col = 'red', cex = 1.2)
# }

# library(vioplot)

# vioplot_vec <- function(xvec, yvec, nm, num_bins = 7, ylab = '', xlab = '') {
#     xvec_rng <- c(min(xvec, na.rm = T),max(xvec, na.rm = T))
#     cvc <- seq(xvec_rng[[1]], xvec_rng[[2]], l = num_bins)
#     mtfc <- cut(xvec, breaks = cvc)
#     vioplot(yvec ~ mtfc, main = nm, ylab = ylab, xlab  = xlab)
#     text(1:length(levels(mtfc)), rep(quantile(yvec, 0.99), length(levels(mtfc))), labels = table(mtfc), col = 'red', cex = 1.2)
# }

# # png('./output/mcatac/figs/nsc_ipc_all_peaks_delta_atac_features_boxplots.png', h = 3000, w = 3000/8, res = 150)
# # par(mfrow = c(8,1), mar = c(5,5.5,2,1.5), cex.lab = 2, cex.main = 1.5, las = 2, cex.axis = 0.8)
# pks_here <- names(nsc_ipc_pks_ord)
# yall <- a_legc_avg_ct[pks_here,'IPC'] - a_legc_avg_ct[pks_here,'NSC']
# # boxplot_vec(rm_ama[pks_here], yall, 'Methylation', ylab = 'Delta ATAC')
# # boxplot_vec(scales[pks_here], yall, 'Scales', ylab = 'Delta ATAC')
# # boxplot_vec(nuc_per_peak[pks_here,'C'], yall, 'GC content', ylab = 'Delta ATAC')
# # uuu <- sapply(motifs_to_take[c(3,4,1,2,5)], function(mtfi) {
# #     boxplot_vec(nsc_ipc_peaks_motifs[pks_here,mtfi], yall, toupper(purrr::map(stringr::str_split(purrr::map(stringr::str_split(mtfi, '\\.'), 2), '\\.|_'), 1)), ylab = 'Delta ATAC')
# # })
# # dev.off()
# # ng2_rng <- c(min(nsc_ipc_peaks_motifs[pks_here,mtfi], na.rm = T),
# #             max(nsc_ipc_peaks_motifs[pks_here,mtfi], na.rm = T))
# #     cvc <- seq(ng2_rng[[1]], ng2_rng[[2]], l = 9)
# #     mtfc <- cut(nsc_ipc_peaks_motifs[pks_here,mtfi], breaks = cvc)
# #     boxplot(yall ~ mtfc, main = toupper(purrr::map(stringr::str_split(purrr::map(stringr::str_split(mtfi, '\\.'), 2), '\\.|_'), 1)), ylab = 'Delta ATAC', xlab  = '')

# col_function <- function(xvec, groupvec, clrmp) {
#     meds <- tapply(xvec, groupvec, mean, na.rm = T)
#     colvals <- clrmp[1+round((length(clrmp)-2)*(meds - min(xvec, na.rm = T))/(max(xvec, na.rm = T) - min(xvec, na.rm = T)))]
#     return(colvals)
# }

# clrmp_rel2 <- colorRampPalette(c('blue4','blue1', 'white', 'red1','red4'))(1000)

# nsc_ipc_pks_ord2 <- do.call('c', lapply(as.numeric(names(sort(a_legc_avg_cln_nsc, decreasing = T))), function(clj) nsc_ipc_pks_ord[nsc_ipc_pks_ord ==clj]))

# nsc_ipc_pks_ord2 <- factor(nsc_ipc_pks_ord, levels = rev(as.numeric(names(sort(a_legc_avg_cln_nsc, decreasing = T)))))
# # nsc_ipc_pks_ord2 <- factor(nsc_ipc_pks_ord, levels = as.numeric(names(sort(a_legc_avg_cln_nsc, decreasing = T))))

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# library(vioplot)

# png('./output/mcatac/figs/nsc_ipc_all_peaks_all_clusters_delta_atac_features_boxplots.png', h = 1200, w = 900, res = 100)
# par(mfrow = c(1,6), mar = c(8,3,2,1), cex.lab = 2, cex.axis = 2, cex.main = 2)
# ylimi <- c(0.85,length(cln_plot)+0.15)
# pks_here <- names(nsc_ipc_pks_ord)
# horiz <- TRUE
# vioplot(delta_ipc_nsc[pks_here] ~ nsc_ipc_pks_ord2, xlab = '', ylab = '', main = '',col = col_function(delta_ipc_nsc[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# title(xlab = 'Delta\nATAC', line = 5)
# # vioplot(a_legc_avg_ct[pks_here,'NSC'] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('NSC\nATAC'))), xlab = '', main = '',col = col_function(a_legc_avg_ct[pks_here,'NSC'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# # vioplot(a_legc_avg_ct[pks_here,'IPC'] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('IPC\nATAC'))), xlab = '', main = '',col = col_function(a_legc_avg_ct[pks_here,'IPC'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# # vioplot(manual_scale[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('Logistic\nscale'))), xlab = '', main = '',col = col_function(manual_scale[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# vioplot(rm_ama[pks_here] ~ nsc_ipc_pks_ord2, xlab = '', ylab = '', main = '',col = col_function(rm_ama[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# title(xlab = 'NSC\nMethylation', line = 5)
# # vioplot(scales[pks_here] ~ nsc_ipc_pks_ord2, ylab = 'Logistic\nscale', col = col_function(scales[pks_here], nsc_ipc_pks_ord2, clrmp_rel2))
# # vioplot(nuc_per_peak[pks_here,'C'] ~ nsc_ipc_pks_ord2, ylab = 'GC\ncontent',xlab = '', main = '',col = col_function(nuc_per_peak[pks_here,'C'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # vioplot(dinucs_per_peak[pks_here,'CG'] ~ nsc_ipc_pks_ord2, xlab = substitute(paste(bold('CpG \ncontent'))),ylab = '', main = '',ylim = c(0,0.05), col = col_function(nuc_per_peak[pks_here,'C'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# cns <- colnames(prego_ie_mat)
# cnlst <- unlist(purrr::map(stringr::str_split(string = cns, pattern = '\\.'), 2))
# yyy <- sapply(grep('sox|eomes|olig2$|neuro', cnlst, ign = T), function(i) {vioplot(prego_ie_mat[pks_here,cns[[i]]] ~ nsc_ipc_pks_ord2, ylab = '', main = '', xlab = '', 
#                                                                    col = col_function(prego_ie_mat[pks_here,cns[[i]]], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi); title(xlab = cnlst[[i]], line = 3)})
# # vioplot(atac_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, xlab = substitute(paste(bold('ATAC UMI\n(20kbp radius)'))),ylab = '',main = '', col = col_function(atac_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# # vioplot(rna_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, xlab = substitute(paste(bold('RNA UMI\n(500kbp radius)'))), ylab = '',main = '',col = col_function(rna_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = horiz, xlim = ylimi)
# dev.off()

# cnlst <- unlist(purrr::map(stringr::str_split(string = cns, pattern = '\\.'), 2))

# pks_here <- names(nsc_ipc_pks_ord2)

# options(repr.plot.width = 12)
# options(repr.plot.height = 12)

# # png('./output/mcatac/figs/nsc_ipc_all_peaks_all_clusters_delta_atac_features_boxplots_vert.png', w = 1600, h = 1600, res = 100)
# cnis <- grep('methyl|JASPAR|JOLMA|HOCOMOCO|HOMER|CG|AG|prox', colnames(x_all), v=T)
# par(mfrow = c(7+1,1), mar = c(3,8,2,1), cex.lab = 2, cex.axis = 2, cex.main = 2)
# ylimi <- c(0.85,length(cln_plot)+0.15)
# # pks_here <- intersect(pb, names(nsc_ipc_pks_ord))

# pks_here <- nsc_ipc_pks_ord2[names(nsc_ipc_pks_ord2) %in% rownames(x_all)]
# boxplot(delta_ipc_nsc[names(pks_here)] ~ pks_here, ylab = 'Delta\nATAC', xlab = '', main = '',col = col_function(delta_ipc_nsc[names(pks_here)], pks_here, clrmp_rel2), horizontal = F, xlim = ylimi)
# ttt <- sapply(cnis, function(cni) {
#     boxplot(x_all[names(pks_here),cni] ~ pks_here, xlab = '', ylab = cni, main = '',col = col_function(x_all[names(pks_here),cni], pks_here, clrmp_rel2), horizontal = F, xlim = ylimi, ylim = quantile(x_all[names(pks_here),cni], c(0.05,0.95)))
# })

# # vioplot(delta_ipc_nsc[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('Delta\nATAC'))), xlab = '', main = '',col = col_function(delta_ipc_nsc[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # # vioplot(manual_scale[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('Logistic\nscale'))), xlab = '', main = '',col = col_function(manual_scale[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # vioplot(rm_ama[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('Methylation'))),xlab = '', main = '',col = col_function(rm_ama[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # # vioplot(scales[pks_here] ~ nsc_ipc_pks_ord2, ylab = 'Logistic\nscale', col = col_function(scales[pks_here], nsc_ipc_pks_ord2, clrmp_rel2))
# # # vioplot(nuc_per_peak[pks_here,'C'] ~ nsc_ipc_pks_ord2, ylab = 'GC\ncontent',xlab = '', main = '',col = col_function(nuc_per_peak[pks_here,'C'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # vioplot(dinucs_per_peak[pks_here,'CG'] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('CpG \ncontent'))),xlab = '', main = '',ylim = c(0,0.05), col = col_function(nuc_per_peak[pks_here,'C'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # cns <- colnames(prego_ie_mat)
# # cnlst <- unlist(purrr::map(stringr::str_split(string = cns, pattern = '\\.'), 2))
# # yyy <- sapply(grep('ndf1|tcf3|tbx5', cnlst, ign = T), function(i) vioplot(prego_ie_mat[pks_here,cns[[i]]] ~ nsc_ipc_pks_ord2, xlab = '', main = '', ylab = cnlst[[i]], 
# #                                                                    col = col_function(prego_ie_mat[pks_here,cns[[i]]], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi))
# # vioplot(atac_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('ATAC UMI\n(50kbp radius)'))),xlab = '',main = '', col = col_function(atac_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # vioplot(rna_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, ylab = substitute(paste(bold('RNA UMI\n(50kbp radius)'))), xlab = '',main = '',col = col_function(rna_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = F, xlim = ylimi)
# # dev.off()



# ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')
# ipc_module
# astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')
# astro_module
# stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')
# stem_module

# tss <- gintervals.load('intervs.global.tss')

# tss <- tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]

# tss2 <- tss[!duplicated(tss$geneSymbol),]

# mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

# prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

# dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

# ipc_module_tss <- dplyr::mutate(tss2[tss2$geneSymbol %in% ipc_module,c('chrom', 'start', 'end', 'geneSymbol')], type = 'ipc_tss') %>% filter(!duplicated(.$geneSymbol))
# astro_module_tss <- dplyr::mutate(tss2[tss2$geneSymbol %in% astro_module,c('chrom', 'start', 'end', 'geneSymbol')], type = 'astro_tss') %>% filter(!duplicated(.$geneSymbol))
# stem_module_tss <- dplyr::mutate(tss2[tss2$geneSymbol %in% stem_module,c('chrom', 'start', 'end', 'geneSymbol')], type = 'stem_tss') %>% filter(!duplicated(.$geneSymbol))
# neuro_tss <- dplyr::mutate(tss2[tss2$geneSymbol %in% neuro_genes,c('chrom', 'start', 'end', 'geneSymbol')], type = 'neuro_tss') %>% filter(!duplicated(.$geneSymbol))

# ipc_module_peaks <- gintervals.neighbors(ipc_module_tss, mcp, mindist = 0, maxdist = 1e+6, maxneighbors = 1e+6) 
# ipc_module_peaks <- ipc_module_peaks[ipc_module_peaks$peak_name %in% intersect(var_peaks, dist_peaks$peak_name),]

# astro_module_peaks <- gintervals.neighbors(astro_module_tss, mcp, mindist = 0, maxdist = 1e+6, maxneighbors = 1e+6) 
# astro_module_peaks <- astro_module_peaks[astro_module_peaks$peak_name %in% intersect(var_peaks, dist_peaks$peak_name),]

# stem_module_peaks <- gintervals.neighbors(stem_module_tss, mcp, mindist = 0, maxdist = 1e+6, maxneighbors = 1e+6) 
# stem_module_peaks <- stem_module_peaks[stem_module_peaks$peak_name %in% intersect(var_peaks, dist_peaks$peak_name),]

# ipc_module_a_legc <- a_legc[ipc_module_peaks$peak_name,]
# astro_module_a_legc <- a_legc[astro_module_peaks$peak_name,]
# stem_module_a_legc <- a_legc[stem_module_peaks$peak_name,]

# km_ipc_module_a_legc <- tglkmeans::TGL_kmeans(ipc_module_a_legc[,nsc_mcs], seed = 1337, k = round(nrow(ipc_module_peaks)/1000))

# km_astro_module_a_legc <- tglkmeans::TGL_kmeans(astro_module_a_legc[,nsc_mcs], seed = 1337, k = round(nrow(astro_module_peaks)/500))
# km_stem_module_a_legc <- tglkmeans::TGL_kmeans(stem_module_a_legc[,nsc_mcs], seed = 1337, k = round(nrow(stem_module_peaks)/500))

# ipc_module_a_legc_avg_km <- tgs_matrix_tapply(t(ipc_module_a_legc), km_ipc_module_a_legc$cluster, mean)

# astro_module_a_legc_avg_km <- tgs_matrix_tapply(t(astro_module_a_legc), km_astro_module_a_legc$cluster, mean)
# stem_module_a_legc_avg_km <- tgs_matrix_tapply(t(stem_module_a_legc), km_stem_module_a_legc$cluster, mean)

# nsc_mcs_ord <- cust_mc_ord_st[names(cust_mc_ord_st) == 'NSC']

# pheatmap::pheatmap(ipc_module_a_legc_avg_km[,cust_mc_ord_st] - rowMeans(ipc_module_a_legc_avg_km), col = clrmp_rel, breaks = 0.3*brks_rel, cluster_cols = F, treeheight_row = 0, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F)

# pheatmap::pheatmap(astro_module_a_legc_avg_km[,cust_mc_ord_st] - rowMeans(astro_module_a_legc_avg_km), col = clrmp_rel, breaks = 0.3*brks_rel, cluster_cols = F, treeheight_row = 0, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F)

# pheatmap::pheatmap(stem_module_a_legc_avg_km[,cust_mc_ord_st] - rowMeans(stem_module_a_legc_avg_km), col = clrmp_rel, breaks = 0.3*brks_rel, cluster_cols = F, treeheight_row = 0, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F)

# sort(table(km_a_legc$cluster[match(stem_module_peaks[km_stem_module_a_legc$cluster == 1,'peak_name'], rownames(a_legc))]))

# uni <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('Astrocytes', 'NSC', 'IPC_cyc', 'IPC')]
# gci <- match(unique(names(uni)), names(uni))-1

# pheatmap::pheatmap(a_legc[stem_module_peaks[km_stem_module_a_legc$cluster == 1,'peak_name'],uni],gaps_col = gci, show_rownames = F, clustering_method = 'ward.D2',
#                    annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, color = clrmp_abs, breaks = seq(-16.6,-14.5,l=1000))

# v1 <- sort(a_legc_avg_cl_ct[,'iCfuPN'] - a_legc_avg_cl_ct[,'iCPN_late'])
# v2 <- sort(a_legc_avg_cl_ct[,'CthPN'] - a_legc_avg_cl_ct[,'CPN_L2-3'])
# v3 <- sort(a_legc_avg_cl_ct[,'NSC'] - a_legc_avg_cl_ct[,'IPC'])
# mod_in <- as.numeric(union(union(names(c(head(v1, 3), tail(v1, 3))), names(c(head(v2, 3), tail(v2, 3)))),
#                           names(c(head(v3, 3), tail(v3, 3)))))
# mod_in

# nsc_mcs <- which(mcmd$cell_type == 'NSC')

# ipc_cyc_mcs <- which(mcmd$cell_type == 'IPC_cyc')

# nsc_cl <- which(a_legc_avg_cl_ct[,'NSC'] - rowMaxs(subset(a_legc_avg_cl_ct, select = -NSC)) >= 0.25)

# options(repr.plot.height = 10)
# options(repr.plot.width = 13)

# col_annot$phase <- mc_median_phase

# pltmt <- a_legc_avg_cl[peak_clust_var,ro]

# ordn <- order(apply(pltmt, 1, which.max))

# peak_clust_var[ordn]

# pheatmap::pheatmap(pltmt[ordn,]  - rowMeans(pltmt[ordn,]), cluster_rows = F, cluster_cols =F, annotation_col = col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D',
#                    col = clrmp_rel, breaks = 0.5*brks_rel, fontsize = 18, annotation_legend = F, treeheight_row = 0, show_colnames = F)

# pheatmap::pheatmap(pltmt[ordn,]  - rowMeans(pltmt[ordn,]), cluster_rows = F, cluster_cols =F, annotation_col = col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D',
#                    col = clrmp_rel, breaks = 0.5*brks_rel, fontsize = 18, annotation_legend = F, treeheight_row = 0, show_colnames = F)

# options(repr.plot.height = 18)
# options(repr.plot.width = 13)

# pks_ordn <- do.call('c', lapply(ordn, function(x) setNames(rep(x, km_a_legc$size[[x]]), rownames(a_legc)[km_a_legc$cluster == x])))

# pltmt <- a_legc[names(pks_ordn),ro]

# p_single_peak_a_legc <- pheatmap::pheatmap(pltmt-rowMeans(pltmt), silent = T, cluster_rows = F, cluster_cols =F, annotation_col = col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D', show_rownames = F,
#                    col = clrmp_rel, breaks = 0.5*brks_rel, fontsize = 18, annotation_legend = F, treeheight_row = 0, show_colnames = F)

# a_legc_mrm <- a_legc_avg_cl - rowMeans(a_legc_avg_cl)

# mod_in <- mod_in[order(apply(a_legc_mrm[mod_in,cust_mc_ord_st], 1, function(x) sum(ifelse(x < 0, 0, x)*(1:length(x)))/sum(ifelse(x < 0, 0, x))))]
# mod_in

# options(repr.plot.height = 18)
# options(repr.plot.width = 13)

# length(mod_in)

# png('./output/mcatac/figs/differential_atac_clusters_nsc_ipc_scatter.png', h = 2000, w = 1000)
# mari <- rep(0.2, 4)
# par(mar = mari, mfrow = c(length(mod_in), 4), cex.lab = 6, cex.axis = 2)
# inds_nsc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'Astrocytes', 'OPCs')]
# inds_ipc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('IPC', 'IPC_cyc')]
# inds_imm <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('iCPN/CfuPN', 'iCfuPN', 'iCPN_early', 'iCPN_late')]
# inds_neu <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('CthPN', 'SCPN','CPN_L2-3', 'CPN_L5_6')]
# # inds_gli <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('Astrocytes', 'OPCs')]
# tvt <- sapply(head(mod_in, 14), function(md) {
#     mal <- max(a_legc_avg_cl[md,])
#     if (mal > -15) {
#         ylimi <- c(-16.7, mal)
#     } else {ylimi <- c(-16.7,-15)}
#     if (md == tail(mod_in, 1)) {
#         mari[[1]] <- 4
#         # mari <- c(4,0.2,0.2,0.2)
#         xaxti = 's'
#     } else {
#         mari[[1]] <- 0.2
#         # mari <- c(0.2,0.2,0.2,0.2)
#         xaxti = 'n'
#     }
#     # plot(mcmd$mean_day[inds_gli], legc_mrm[md,inds_gli], pch = 16, cex = 2, col = mcmd$color[inds_gli], 
#     #      xlim = c(15.5,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     # grid(col = 'pink')
#     tmpm <- mari[[2]]
#     mari[[2]] <- 8
#     par(mar = mari)
#     plot(mcmd$mean_day[inds_nsc], a_legc_avg_cl[md,inds_nsc], pch = 16, cex = 2, col = mcmd$color[inds_nsc], 
#          xlim = c(13,18), ylim = ylimi, ylab = md, xaxt =xaxti, xlab = '')
#     grid(col = 'pink')
#     mari[[2]] <- tmpm
#     par(mar = mari)
#     plot(mcmd$mean_day[inds_ipc], a_legc_avg_cl[md,inds_ipc], pch = 16, cex = 2, col = mcmd$color[inds_ipc], 
#          xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
#     grid(col = 'pink')
#     plot(mcmd$mean_day[inds_imm], a_legc_avg_cl[md,inds_imm], pch = 16, cex = 2, col = mcmd$color[inds_imm], 
#          xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
#     grid(col = 'pink')
#     plot(mcmd$mean_day[inds_neu], a_legc_avg_cl[md,inds_neu], pch = 16, cex = 2, col = mcmd$color[inds_neu], 
#          xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
#     grid(col = 'pink')
# })
# dev.off()

# genes_mods_in_sort <- setNames(unlist(sapply(mod_in, function(md) rep(md, length(which(km_a_legc$cluster == md))))), 
#                                              unlist(sapply(mod_in, function(md) rownames(a_legc)[km_a_legc$cluster == md])))
# head(genes_mods_in_sort)

# library(ComplexHeatmap)

# top_ha <- columnAnnotation(cell_type = anno_simple(x = mcmd$cell_type[cust_mc_ord_st2], col = tibble::deframe(color_key)),
#                           mean_day = anno_simple(x = mcmd$mean_day[cust_mc_ord_st2], 
#                                                  col = circlize::colorRamp2(breaks = seq(13,18,l=6), 
#                                                             colors = c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))))

# # p_mod_in_genes <- ComplexHeatmap::Heatmap(a_legc[names(genes_mods_in_sort),cust_mc_ord_st2] - rowMeans(a_legc[names(genes_mods_in_sort),]), 
# p_mod_in_genes <- ComplexHeatmap::Heatmap(a_legc[names(genes_mods_in_sort),cust_mc_ord_st2], 
#                                           row_title_gp = gpar(fontsize = 44),
#                                           show_heatmap_legend = F,show_row_names = F,
# #                                           row_names_gp = gpar(fontsize = 10),
                                          
#                                           split = factor(genes_mods_in_sort, levels = mod_in),
#                                           # fontsize_row = 5,
#                            # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
#                            cluster_rows = F,
#                                           show_column_names = F, 
#                                           # legend = T,
#          # annotation_row = rah,
#          cluster_columns = F, top_annotation = top_ha,
# #          col = circlize::colorRamp2(breaks = seq(-1,1,l=3), colors = c('blue3', 'white', 'red3')),
#         col = circlize::colorRamp2(breaks = seq(-16.6,-14.5,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
#         # breaks = seq(-3,3,l=100)
#                                          )

# png( './output/mcatac/figs/mod_in_by_peak.png', h = 2500, w = 1100)
# draw(p_mod_in_genes)
# dev.off()

# # md_cut <- cut(mcmd$mean_day, 13:18)
# # nsc_legc_t <- t(tgs_matrix_tapply(a_legc[,nsc_mcs], md_cut[nsc_mcs], mean))
# # pks_nsc_trans <- which(apply(nsc_legc_t, 1 , which.max) %in% 2:4 & rowMaxs(nsc_legc_t) - rowMins(nsc_legc_t) >= 0.75)

# # pltmt  <- t(apply(a_legc[pks_nsc_trans,nsc_mcs[order(mcmd$mean_day[nsc_mcs])]], 1, zoo::rollmean, k = 20, na.pad = F))
# # # hc <- hclust(dist(pltmt), method = 'ward.D2')

# # p1 <- pheatmap::pheatmap(pltmt, cluster_rows = T, cluster_cols = F, treeheight_row = 0, fontsize_row = 3, annotation_col = col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D2')

# # save_pheatmap(p1, './output/mcatac/figs/trans_peaks_in_IPC.png', h = 2600, w = 1400)

# # amd <- prego::all_motif_datasets()

# # tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
# # tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(legc))
# # # tfs

# # ## TF set which is more like "markers"

# # tf_motif_names <- unlist(sapply(tfs_hi, function(tf) unique(grep(tf, amd$motif, v=T))))
# # tf_motif_names

# # amd_in <- dplyr::filter(amd, motif %in% tf_motif_names)


# # readr::write_tsv(x = amd_in, file = './BonevCollab/pwms_of_tfs_in_data.tsv')

# rna_freqs <- table(mc_rna@mc)/length(mc_rna@mc)
# rna_ct_freqs <- tapply(as.numeric(rna_freqs), mcmd$cell_type, sum)
# atac_freqs <- table(as.numeric(mca@cell_to_metacell$metacell))/length(mca@cell_to_metacell$metacell)
# atac_ct_freqs <- tapply(as.numeric(atac_freqs), mcmd$cell_type, sum)
# atac_ct_freqs
# plot(as.numeric(rna_ct_freqs), as.numeric(atac_ct_freqs), col = color_key$color[match(names(atac_ct_freqs), color_key$cell_type)], pch = 16, cex = 1.5)
# abline(0,1,col='red')
# abline(-0.01,1,lty = 2)
# abline(0.01,1,lty = 2)
# par(cex.main = 2, cex.lab = 2, mar = c(5,5,4,2), cex.axis = 2)
# plot(as.numeric(rna_freqs), (as.numeric(atac_freqs) - as.numeric(rna_freqs))/pmax(as.numeric(rna_freqs)), col = mcmd$color, pch = 16, cex = 1.5, xlab = 'RMC frequency', ylab = '(f_ATAC - f_RNA)/f_RNA', main = 'Normalized deviation of fAMC\nfrequency from RMC frequency')
# abline(0,1,col='red')
# abline(-0.001,1,lty = 2)
# abline(0.001,1,lty = 2)

# par(las = 2, mar = c(12, 6, 4, 2), cex.main = 2, cex.axis = 1.5, cex.lab = 2,mgp =c(3,0,0))
# xbp <- (rna_ct_freqs - atac_ct_freqs)/rna_ct_freqs
# barplot(xbp[cust_st_ord], col = color_key$color[match(cust_st_ord, color_key$cell_type)],ylab = 'Relative deviation', main = 'Deviation of cell type\nproportions in ATAC vs. RNA')
# # grid(lwd = 2)

# # day_mcl_path <- file.path(wd, 'output/mcatac/pl_cort_day_mcls.rds')

# # day_mcls = readRDS(day_mcl_path)

# # dmcl_mat_all <- do.call('cbind', plyr::llply(day_mcls, function(x) amc_mat <- t(tgs_matrix_tapply(x$prom_mat, x$cor_km$cluster, sum)), .parallel = T))

# # dmcl_mat_norm <- log2(1e-5 + t(t(dmcl_mat_all)/colSums(dmcl_mat_all)))

# # markers <- scdb_gset('pl_cort_marks')
# # markers <- names(markers@gene_set)

# # sapply(seq(0.5,1,0.1), function(trs) {
# #     var_prom <- rownames(dmcl_mat_norm)[which(matrixStats::rowSds(dmcl_mat_norm) >= trs)]
# #     table(markers %in% var_prom)
# # })

# get_genes_specific_to_mcs <- function(legc, mc_pos = NULL, mc_neg = NULL, cl_vec = NULL) {
#     if (!is.null(mc_pos) && is.null(mc_neg)) {
#         cl_vec <- ifelse(1:ncol(legc) %in% mc_pos, 1, 0)
#     } else if (!is.null(mc_pos) && !is.null(mc_neg)) {
#         if (!(length(intersect(mc_pos, mc_neg)) == 0)) {
#             stop('mc_pos and mc_neg intersect')
#         }
#         legc <- legc[,c(mc_pos, mc_neg)]
#         cl_vec <- c(rep(1, length(mc_pos)), rep(0, length(mc_neg)))
#     }
#     legc_avg <- t(tgs_matrix_tapply(legc, cl_vec, mean))
#     if (ncol(legc_avg) == 2) {
#         diffs <- matrixStats::rowDiffs(legc_avg)
#         rownames(diffs) <- rownames(legc_avg)
#         return(diffs[order( diffs[,1], decreasing = T),])
#     } else {
#         diffs <- t(plyr::laply(1:ncol(legc_avg), function(i) legc_avg[,i] - matrixStats::rowMaxs(legc_avg[,-i]), .parallel = T))
#         colnames(diffs) <- colnames(legc_avg)
#         res <- lapply(1:ncol(diffs), function(i) {
#             df <- diffs[which(diffs[,i] > 0.1),]
# #             print(df)
#             if (!is.null(dim(df))) {
#                 return(df[order(df[,i], decreasing = T),])
#             } else {
#                 return(sort(df, decreasing = T))
#             }
            
#         })
#         names(res) <- colnames(legc_avg)
#         return(res)
#     }
# }

# # var_genes <- lapply(get_genes_specific_to_mcs(legc, cl_vec = mcmd$cell_type), function(x) head(rownames(x), 100))

# # names(var_genes) <- sort(unique(mcmd$cell_type))

# # load(file= './output/mcatac/feats_filt_from_cortex_matching.rda')

# # p_feats_filt_amc <- pheatmap(dmcl_mat_norm[feats_filt,], fontsize_row = 24, show_colnames = F, legend = T, labels_row = ifelse(feats_filt %in% union(tfs_hi, markers),feats_filt,''))

# # save_pheatmap(p_prom_var_gene_amc, './output/mcatac/figs/ct_specific_promoters_in_amcs.png', h = 2000, w = 1400, res = 200)

# # save_pheatmap(p_feats_filt_amc, './output/mcatac/figs/feats_filt_promoters_in_amcs.png', h = 4400, w = 3600, res = 200)

# # sc_cor_kms <- readRDS('./output/mcatac/microcluster_assignment.RDS')
# # lapply(seq_along(sc_cor_kms), function(i) length(sc_cor_kms[[i]]$cluster))

# seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(a_legc))

# tss <- dplyr::filter(gintervals.load('intervs.global.tss'), geneSymbol %in% rownames(mca@rna_egc))

# prom_peaks <- gintervals.neighbors(seq_coords, tss, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

# dist_peaks <- seq_coords[!(seq_coords$peak_name %in% prom_peaks$peak_name),]

# a_legc_prom <- a_legc[prom_peaks$peak_name,]

# a_legc_prom <- tgs_matrix_tapply(t(a_legc_prom), prom_peaks$geneSymbol, function(x) log2(mean(2**x)))

# # markers <- names(scdb_gset('pl_cort_marks_f')@gene_set)
# # markers <- markers[matrixStats::rowMaxs(legc[markers,]) < -10]

# goi = c('Pou3f1', 'Pou3f2', 'Cux1', 'Cux2', 'Neurod1', 'Neurog2', 'Id4',
#          'Eomes', 'Hes1', 'Apoe', 'Sox5', 'Tbr1', 'Foxp2', 'Foxp1', 'Nfia', 'Islr2', 
#          'Zbtb20', 'Bcl11b', 'Fezf2', 'Satb2', 'Mef2c', 'Nhlh1', 'Tle4',
#         'Rnd2',  'Runx1t1', 'Mapt', 'Mki67', 'Pcna',
#         'Fabp7', 'Olig1', 'Ldb2', 'Gadd45g', 'Syt4')
# markers <- goi

# pltmt <- legc[markers,cust_mc_ord_st]
# hcp <- hclust(dist(pltmt), method = 'ward.D2')
# genes_ord <- rownames(pltmt)[hcp$order]

# # par(las = 2)
# # prom_cor <- diag(tgs_cor(t(pltmt[genes_ord,]), t(a_legc_prom[genes_ord[genes_ord %in% rownames(a_legc_prom)],cust_mc_ord_st]), spearman = T))
# prom_cor <- diag(tgs_cor(t(pltmt[genes_ord,]), t(a_legc_prom[genes_ord,cust_mc_ord_st]), spearman = T))
# # barplot(prom_cor)

# clrmp <- colorRampPalette(c('white', 'pink', 'red3', 'purple3', 'blue3', 'black'))(100)
# brks <- seq(-16.7, -10, l=100)

# p_goi_rna_legc <- pheatmap::pheatmap(pltmt[genes_ord,], annotation_col = col_annot, annotation_colors = ann_colors, 
#                                      cluster_cols = F, show_colnames = F, cluster_rows = F, 
#                                      color = clrmp, breaks = brks, annotation_legend = F, fontsize_row = 18)
# # save_pheatmap(p_goi_rna_legc, './output/mcatac/figs/phm_goi_rna_legc.png', h = 1000, w = 2800, res = 150)

# # load('./output/mcatac/km_a_legc_k=80.rda')

# a_legc_prom <- a_legc[prom_peaks$peak_name,]

# dim(a_legc_prom)

# quantile(rowMeans(a_legc_prom))

# km_prom_a_legc <- tglkmeans::TGL_kmeans(a_legc[prom_peaks$peak_name,], k = 10, seed = 1337)

# km_enh_a_legc <- tglkmeans::TGL_kmeans(a_legc[dist_peaks$peak_name,], k = 60, seed = 1337)

# a_legc_avg_cl_prom <- tgs_matrix_tapply(t(a_legc[prom_peaks$peak_name,]), km_prom_a_legc$cluster, mean)

# a_legc_avg_cl_enh <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)

# a_legc_avg_cl_ct_enh <- t(tgs_matrix_tapply(a_legc_avg_cl_enh, mcmd$cell_type, mean))

# head(a_legc_avg_cl_ct_enh)

#     options(repr.plot.width=  12)
#     options(repr.plot.height=  6)
# par(mfrow = c(1,2))
# pv1 <- cor.test(a_legc_avg_cl_ct_enh[rownames(m2),'CthPN'], motifs_anno_mat[rownames(m2),'POU3F2'])$p.value
# plot(a_legc_avg_cl_ct_enh[rownames(m2),'CthPN'], motifs_anno_mat[rownames(m2),'POU3F2'], main = glue::glue("cor = {round(cor(a_legc_avg_cl_ct_enh[rownames(m2),'CthPN'], motifs_anno_mat[rownames(m2),'POU3F2']), 3)}, pv = {round(pv1, 3)}"))
# pv1 <- cor.test(a_legc_avg_cl_ct_enh[rownames(m2),'CPN_L2-3'], motifs_anno_mat[rownames(m2),'POU3F2'])$p.value
# plot(a_legc_avg_cl_ct_enh[rownames(m2),'CPN_L2-3'], motifs_anno_mat[rownames(m2),'POU3F2'], main =glue::glue("cor = {round(cor(a_legc_avg_cl_ct_enh[rownames(m2),'CPN_L2-3'], motifs_anno_mat[rownames(m2),'POU3F2']), 3)}, pv = {round(pv1, 3)}"))

# colnames(motifs_anno_mat)

#     options(repr.plot.width=  8)
#     options(repr.plot.height=  6)
# pheatmap::pheatmap(a_legc_avg_cl_prom[,cust_mc_ord_st2], annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, col = clrmp_abs, clustering_method = 'ward.D2')

# a_legc_avg_cl_enh <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
# cl_cvs <- rowSds(a_legc_avg_cl_enh)**2/(rowMeans(a_legc_avg_cl_enh) + 17.6)**2

# m1 <- a_legc_avg_cl_enh[head(names(sort(cl_cvs)),13),cust_mc_ord_st2]

# m1_hc <- hclust(dist(m1), method = 'ward.D2')
# m2 <- a_legc_avg_cl_enh[tail(names(sort(cl_cvs)),-13),cust_mc_ord_st2]
# m2_hc <- hclust(dist(m2), method = 'ward.D2')
# m2_ord <- order(apply(m2, 1, function(x) sum(x*1:length(x))/sum(x)))

# m2_hcct <- cutree(m2_hc, k = 16)

# m2_hcct_ord <- order(apply(tgs_matrix_tapply(t(m2), m2_hcct, mean), 1, function(x) sum(x*1:length(x))/sum(x)))
# m2_hcct_ord
# # m2_hcct[m2_hc$order]

# m2_cl_ord <- do.call('c', sapply(m2_hcct_ord, function(u) names(m2_hcct[m2_hcct == u])))

# enh_cl_ord <- c(rownames(m1)[m1_hc$order], rev(m2_cl_ord))

# var_peaks <- dist_peaks$peak_name[km_enh_a_legc$cluster %in% as.numeric(rownames(m2))]

# save(prom_peaks, dist_peaks, enh_cl_ord, km_enh_a_legc, km_prom_a_legc, m1, m2, var_peaks, file = './output/mcatac/var_peaks_after_enh_prom_separation.rda')

# length(var_peaks)

# nrow(dist_peaks)

# nrow(dist_peaks) - length(var_peaks)

#     options(repr.plot.width=  8)
#     options(repr.plot.height=  10)
#     # pheatmap::pheatmap(rbind(m1[rownames(m1)[m1_hc$order],], m2[rownames(m2)[m2_hc$order],]), cluster_rows = F, annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, col = clrmp_abs, clustering_method = 'ward.D')
#     pheatmap::pheatmap(rbind(m1[rownames(m1)[m1_hc$order],], m2[rev(m2_cl_ord),]), cluster_rows = F, annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, col = clrmp_abs, clustering_method = 'ward.D')

# prom_mxmns <- rowMaxmins(a_legc[prom_peaks$peak_name,])

# prom_sds <- setNames(rowSds(a_legc[prom_peaks$peak_name,]), prom_peaks$peak_name)

# Q_THRESH <- 0

# marks <- names(scdb_gset('pl_cort_marks_f')@gene_set)
# marks

# astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')

# ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')

# stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')

# proms_hi_var_astro <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% astro_module]
# proms_hi_var_ipc <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% ipc_module]
# proms_hi_var_stem <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% stem_module]
# proms_hi_var_marks <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% marks]

# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH)], astro_module))/length(astro_module)
# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH)], ipc_module))/length(ipc_module)
# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds,Q_THRESH)], stem_module))/length(stem_module)

# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH)], astro_module))/length(astro_module) + 
# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH)], ipc_module))/length(ipc_module) + 
# length(intersect(prom_peaks$geneSymbol[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds,Q_THRESH)], stem_module))/length(stem_module)

# pvals <- lapply(list(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem), function(x) {
#     p1 <- ks.test(prom_sds[x], prom_sds[setdiff(proms_hi_var_marks, unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem))))], alternative = 'less')$p.value
#     p2 <- ks.test(prom_sds[x], prom_sds[setdiff(names(prom_sds), unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem))))], alternative = 'less')$p.value
#     return(c(p1, p2))
#     })
    

# p.adjust(unlist(pvals))

# png('./output/mcatac/figs/tss-proximal_peak_variance_densities.png')
# par(cex.lab = 2, mar = c(5,5,2,1), cex.main = 2, cex.axis = 2)
# BW <- 0.1
# d1 <- density(prom_sds[setdiff(names(prom_sds), unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem, proms_hi_var_marks))))]**2, bw = BW)
# d2 <- density(prom_sds[proms_hi_var_astro]**2, bw = BW)
# d3 <- density(prom_sds[proms_hi_var_ipc]**2, bw = BW)
# d4 <- density(prom_sds[proms_hi_var_stem]**2, bw = BW)
# d5 <- density(prom_sds[setdiff(proms_hi_var_marks, unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem))))]**2, bw = BW)
# plot(d1$x, d1$y, col = 'black', type = 'l', lwd = 2, lty = 2, xlab = 'ATAC variance across metacells', ylab = glue::glue('Density, BW = {BW}'), 
#      main = 'TSS-proximal peak ATAC variance'
#     )
# lines(d2$x, d2$y, col = col_key[['Astrocytes']], lwd = 3, lty = 2)
# lines(d3$x, d3$y, col = col_key[['IPC']], lwd = 3, lty = 2)
# lines(d4$x, d4$y, col = col_key[['NSC']], lwd = 3, lty = 2)
# lines(d5$x, d5$y, col = 'darkorange', lwd = 3, lty = 2)
# legend('topright', legend = c('All TSS-proximal peaks', 'Near astro module TSSs', 'Near IPC module TSSs', 'Near stem module TSSs', 'Near marker gene TSSs'),
#        col = c('black', col_key[c('Astrocytes', 'IPC', 'NSC')], 'darkorange'), lty = rep(2,5), lwd = rep(2,5), cex = 1)
# dev.off()





# # km_a_legc <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
# load('./output/mcatac/var_peaks_after_enh_prom_separation.rda')

# a_legc_avg_cl <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)

# load('./output/metacell_model/pcurve_nsc_vs_mat_neuro.rda')

# pks_diff_cfupn_clpl23 <- rownames(a_legc_avg_ct)[a_legc_avg_ct[,'CthPN'] - a_legc_avg_ct[,'CPN_L2-3'] > 1]
# pks_diff_clpl23_cfupn <- rownames(a_legc_avg_ct)[a_legc_avg_ct[,'CthPN'] - a_legc_avg_ct[,'CPN_L2-3'] < -1]

# length(pks_diff_cfupn_clpl23)
# length(pks_diff_clpl23_cfupn)

# pks_cfupn_density <- 1/(sum(.misha$ALLGENOME[[1]]$end)/length(pks_diff_cfupn_clpl23))
# pks_cpn_density <- 1/(sum(.misha$ALLGENOME[[1]]$end)/length(pks_diff_clpl23_cfupn))

# pks_cfupn_density
# pks_cpn_density

# head(pks_diff_cfupn_clpl23)
# head(pks_diff_clpl23_cfupn)

# dim(legc)

# legc_avg_ct <- t(tgs_matrix_tapply(legc[,mcmd$metacell], mcmd$cell_type, mean))

# bias_cpnl23_vs_cthpn <- legc_avg_ct[,'CPN_L2-3'] - legc_avg_ct[,'CthPN']

# DDD <- 1e+5

# cor_mat_neuro_genes <- tgs_cor(t(legc), t(legc[c('Mef2c', 'Mapt', 'Runx1t1'),]))
# cor_mat_neuro_genes <- cor_mat_neuro_genes[order(rowMeans(cor_mat_neuro_genes)),]

# mat_neuro_genes <- sort(rownames(cor_mat_neuro_genes)[apply(cor_mat_neuro_genes, 1, max) >= 0.9])

# cpn_genes <- get_genes_specific_to_mcs(legc = legc, mc_pos = which(mcmd$cell_type == 'CPN_L2-3'), mc_neg = which(mcmd$cell_type == 'CthPN'))
# cfupn_genes <- sort(names(cpn_genes[cpn_genes <= -2]))
# cpn_genes <- sort(names(cpn_genes[cpn_genes >= 2]))

# cor_mat_neuro_genes <- tgs_cor(t(legc), t(legc[c('Mef2c', 'Mapt', 'Runx1t1'),]))
# cor_mat_neuro_genes <- cor_mat_neuro_genes[order(rowMeans(cor_mat_neuro_genes)),]

# mat_neuro_genes <- sort(rownames(cor_mat_neuro_genes)[apply(cor_mat_neuro_genes, 1, max) >= 0.85])

# length(mat_neuro_genes)

# mat_neuro_genes_nei_cfupn_peaks <- gintervals.neighbors(mcp[mcp$peak_name %in% pks_diff_cfupn_clpl23,], tss[tss$geneSymbol %in% mat_neuro_genes,], maxneighbors = 1e+6, maxdist = 1e+6, mindist = 0)

# mat_neuro_genes_nei_cpn_peaks <- gintervals.neighbors(mcp[mcp$peak_name %in% pks_diff_clpl23_cfupn,], tss[tss$geneSymbol %in% mat_neuro_genes,], maxneighbors = 1e+6, maxdist = 1e+6, mindist = 0)

# nei_cfupn_cpn_peaks <- gintervals.neighbors(mcp[mcp$peak_name %in% pks_diff_clpl23_cfupn,], mcp[mcp$peak_name %in% pks_diff_cfupn_clpl23,], maxneighbors = 1e+6, maxdist = 5e+5, mindist = 0)

# length(unique(nei_cfupn_cpn_peaks[,4]))
# length(unique(nei_cfupn_cpn_peaks[,8]))
# length(unique(nei_cfupn_cpn_peaks[,4]))/length(pks_diff_cfupn_clpl23)
# length(unique(nei_cfupn_cpn_peaks[,8]))/length(pks_diff_clpl23_cfupn)

# nei_cfupn_pks_tads <- gintervals.neighbors(mcp[mcp$peak_name %in% pks_diff_cfupn_clpl23,], tads, maxneighbors = 1e+6, maxdist = 1, mindist = 0)
# nei_cpn_pks_tads <- gintervals.neighbors(mcp[mcp$peak_name %in% pks_diff_clpl23_cfupn,], tads, maxneighbors = 1e+6, maxdist = 1, mindist = 0)

# tads_count_spec_peaks <- dplyr::mutate(tads, cfupn_peaks =  table(nei_cfupn_pks_tads$tad_name)[tads$tad_name], 
#                                        cpn_peaks =  table(nei_cpn_pks_tads$tad_name)[tads$tad_name],
#                                        tad_length = end - start,
#                                        bias_cpnl23_vs_cthpn = round(bias_cpnl23_vs_cthpn[geneSymbol], 2)
#                                       )
# tads_count_spec_peaks <- dplyr::mutate(tads_count_spec_peaks, spec_peak_total = rowSums(tads_count_spec_peaks[,c('cfupn_peaks', 'cpn_peaks')], na.rm = T),
#                                       spec_peak_density = spec_peak_total/tads_count_spec_peaks$tad_length,
#                                      cfupn_enr = cfupn_peaks/(pks_cfupn_density*tad_length),
#                                      cpn_enr = cpn_peaks/(pks_cfupn_density*tad_length))

# nei_peaks_tads <- gintervals.neighbors(mcp, tads, maxdist = 1, mindist = 0, maxneighbors = 1e+6)

# exons <- gintervals.load('intervs.global.exon')

# nei_peaks_exons <- gintervals.neighbors(mcp, exons, maxdist = 1, mindist = 0, maxneighbors = 1e+6)

# nei_peaks_tads$peak_type <- rep('general', nrow(nei_peaks_tads))
# nei_peaks_tads$peak_type[nei_peaks_tads$peak_name %in% pks_diff_cfupn_clpl23] <- 'CfuPN'
# nei_peaks_tads$peak_type[nei_peaks_tads$peak_name %in% pks_diff_clpl23_cfupn] <- 'CPN'

# nei_peaks_tads$is_exon <- ifelse(nei_peaks_tads$peak_name %in% nei_peaks_exons$peak_name, 'TRUE', 'FALSE')

# plot_peaks_tad <- function(tad, fdir) {
#     options(repr.plot.width=  12)
#     options(repr.plot.height=  8)
#     dff <- nei_peaks_tads[nei_peaks_tads$tad_name == tad,]
#     annh <- tibble::column_to_rownames(tibble(dff[,c('peak_name', 'peak_type', 'is_exon')]), 'peak_name')
#     # print(head(annh))
#     colnames(annh) <- c('peak_type', 'is_exon')
#     annc <- ann_colors
#     annc[['peak_type']] <- setNames(gplots::col2hex(c('blue', 'red', 'green')), c('general', 'CfuPN', 'CPN'))
#     annc[['is_exon']] <- setNames(gplots::col2hex(c('black', 'red')), c('FALSE', 'TRUE'))
#     # print(annc)
#     ppi <- pheatmap::pheatmap(a_legc[unique(dff$peak_name),cust_mc_ord_st], silent = T, show_colnames= F, annotation_legend = T, col = clrmp_abs,annotation_row = annh, annotation_col = col_annot, annotation_colors = annc, cluster_cols = F, cluster_rows = F)
#     save_pheatmap(ppi, paste0(fdir, '/', tad, '.png'), h = 1200, w = 1600)
# }

# # dir.create('./output/mcatac/figs/tads_with_ct_specific_peaks')

# plot_peaks_tad("3814_Cux1", './output/mcatac/figs/tads_with_ct_specific_peaks')

# plot_peaks_tad("1445_Fezf2", './output/mcatac/figs/tads_with_ct_specific_peaks')
# plot_peaks_tad("1733_Sybu", './output/mcatac/figs/tads_with_ct_specific_peaks')
# plot_peaks_tad("995_Nrcam", './output/mcatac/figs/tads_with_ct_specific_peaks')
# plot_peaks_tad("4301_Slc17a6", './output/mcatac/figs/tads_with_ct_specific_peaks')
# plot_peaks_tad("3031_Ptx3", './output/mcatac/figs/tads_with_ct_specific_peaks')
# plot_peaks_tad("3780_Cux2", './output/mcatac/figs/tads_with_ct_specific_peaks')

# gintervals.neighbors(tss[tss$geneSymbol == 'Sybu',c('chrom','start', 'end', 'geneSymbol')], tads)
# gintervals.neighbors(tads[tads$tad_name == '1733_Sybu',], tss[,c('chrom','start', 'end', 'geneSymbol')], maxdist = 0, maxneighbors = 1e+6)

# length(unique(mat_neuro_genes_nei_cfupn_peaks$geneSymbol))

# length(unique(mat_neuro_genes_nei_cpn_peaks$geneSymbol))

# length(intersect(mat_neuro_genes_nei_cfupn_peaks$geneSymbol, mat_neuro_genes_nei_cpn_peaks$geneSymbol))

# ct_genes_peaks_dist <- lapply(cust_st_ord, function(cti) {
#     ct_mcs <- which(mcmd$cell_type == cti)
#     ct_peaks_ord <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = ct_mcs)
#     ct_peaks <- names(ct_peaks_ord)[ct_peaks_ord >= 0.25]
#     ct_genes_ord <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = ct_mcs)
# })

# ct_peaks_all <- get_genes_specific_to_mcs(a_legc, cl_vec = mcmd$cell_type)

# npeaks <- get_genes_specific_to_mcs(a_legc, mc_pos = which(mcmd$cell_type == 'CPN_L2-3'))

# cti <- 'NSC'
# quantile(subset(a_legc_avg_ct, select = cti) - rowMaxs(subset(a_legc_avg_ct, select = -NSC)))

# png('./output/mcatac/figs/cpn_l23_vs_cthpn_specific_peaks.png')
# par(cex.lab = 2, cex.axis = 2, mar = c(5,5,1,1))
# colvec <- rep('black', nrow(a_legc_avg_ct))
# colvec[rownames(a_legc) %in% pks_diff_cfupn_clpl23] <- 'red'
# colvec[rownames(a_legc) %in% pks_diff_clpl23_cfupn] <- 'blue'
# plot(a_legc_avg_ct[,'CthPN'], a_legc_avg_ct[,'CPN_L2-3'], col = colvec, pch= 16, cex = 0.25, ylab = '', xlab = '')
# title(ylab = 'CPN_L2-3 ATAC', xlab = 'CthPN ATAC')
# abline(0,1,col='green')
# legend('topleft', legend = c('CthPN-specific', 'CPN_L2-3-specific'), col = c('red', 'blue'), pch = rep(16,2))
# grid()
# dev.off()

# options(repr.plot.width = 12
#        )
# options(repr.plot.height = 6)

# # png('./output/mcatac/figs/cpn_vs_cfupn_specific_atac_vs_diff_ord.png', h = 600, w = 1200)
# par(mfrow = c(1,2), cex.main = 3, mar = c(5,5,3,2), cex.lab = 2, cex.axis = 2)
# plot(1:ncol(a_legc), colMeans(a_legc[pks_diff_cfupn_clpl23,rev(pcu$ord)]), 
#      xlab = 'Differentiation order', 
#      pch =16, cex = 1.5,
#      col = mcmd$color[rev(pcu$ord)], xaxt = 'n',
#      main = paste0('CfuPN-specific peaks, n = ', length(pks_diff_cfupn_clpl23)), ylab = '', ylim = c(-16.7,-14.4))
# title(ylab = 'Mean expression', line = 3)
# grid(col = 'darkgray', lty = 2, lwd = 2)
# par(mar = c(5,1,3,1))
# plot(1:ncol(a_legc), colMeans(a_legc[pks_diff_clpl23_cfupn,rev(pcu$ord)]), 
#      xlab = 'Differentiation order', pch =16,  cex = 1.5,
#      col = mcmd$color[rev(pcu$ord)], 
#      main = paste0('CPN-specific peaks, n = ',  length(pks_diff_clpl23_cfupn)), 
#      ylab= '',xaxt = 'n', ylim = c(-16.7,-14.4))
# grid(col = 'darkgray', lty = 2, lwd = 2)
# # dev.off()

# options(repr.plot.width = 8)
# options(repr.plot.height = 8)

# options(repr.plot.width = 12
#        )
# options(repr.plot.height = 6)

# clrmp_atac <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
# # plot(1:100, rep(1,100), col = clrmp_atac)

# mat_cpn_vs_cth <- a_legc[pks_diff_clpl23_cfupn,]
# hc_cpn_vs_cth <- hclust(dist(mat_cpn_vs_cth), method = 'ward.D')
# mat_cth_vs_cpn <- a_legc[pks_diff_cfupn_clpl23,]
# hc_cth_vs_cpn <- hclust(dist(mat_cth_vs_cpn), method = 'ward.D')

# all(pcu$ord %in% as.numeric(mcmd$metacell))

# p_cpn_peaks <- pheatmap::pheatmap(mat_cpn_vs_cth[hc_cpn_vs_cth$order,cust_mc_ord_st2], cluster_rows = F, annotation_legend = F,
#                    col = clrmp_atac, breaks = seq(-16.6, -13, l = 100), cluster_cols = F, 
#                    annotation_col = subset(col_annot, select = -mean_day), annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
# mcATAC::save_pheatmap(p_cpn_peaks, './output/mcatac/figs/cpn_peaks_phm.png', h = 900, w = 800)
# save_pheatmap_pdf(p_cpn_peaks, './output/mcatac/figs/cpn_peaks_phm.pdf', h = 900/71, w = 800/71)

# p_cthpn_peaks <- pheatmap::pheatmap(mat_cth_vs_cpn[hc_cth_vs_cpn$order,cust_mc_ord_st2], cluster_rows = F, annotation_legend = F,
#                    col = clrmp_atac, breaks = seq(-16.6, -13, l = 100), 
#                    cluster_cols = F, annotation_col = subset(col_annot, select = -mean_day), 
#                    annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
# mcATAC::save_pheatmap(p_cthpn_peaks, './output/mcatac/figs/cthpn_peaks_phm.png', h = 900, w = 800)
# save_pheatmap_pdf(p_cthpn_peaks, './output/mcatac/figs/cthpn_peaks_phm.pdf', h = 900/71, w = 800/71)

# # p_cpn_peaks <- pheatmap::pheatmap(mat_cpn_vs_cth_sm[hc_cpn_vs_cth$order,], cluster_rows = F, annotation_legend = F,
# #                    col = clrmp_atac, breaks = seq(-16.6, -13, l = 100), cluster_cols = F, 
# #                    annotation_col = col_annot, annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
# # # save_pheatmap(p_cpn_peaks, './output/mcatac/figs/cpn_peaks_phm.png', h = 800, w = 800)

# # p_cthpn_peaks <- pheatmap::pheatmap(mat_cth_vs_cpn_sm[hc_cth_vs_cpn$order,], cluster_rows = F, annotation_legend = F,
# #                    col = clrmp_atac, breaks = seq(-16.6, -13, l = 100), 
# #                    cluster_cols = F, annotation_col = col_annot, 
# #                    annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
# # # save_pheatmap(p_cthpn_peaks, './output/mcatac/figs/cthpn_peaks_phm.png', h = 800, w = 800)

# mat_cpn_vs_cth_sm <- t(apply(mat_cpn_vs_cth[,rev(pcu$ord)], 1, zoo::rollmean, k = 10))

# png('./output/mcatac/figs/diff_order_index_of_max_acc_cthpn_vs_cpnl23.png')
# par(mar = c(8,4,5,1), cex.main = 2.5, cex.axis = 2)

# boxplot(c(apply(mat_cpn_vs_cth_sm,1,which.max),
#             apply(mat_cth_vs_cpn_sm,1,which.max)) ~ 
#                   c(rep('cpn', nrow(mat_cpn_vs_cth_sm)), rep('cfupn', nrow(mat_cth_vs_cpn_sm))), xaxt = 'n', xlab = '',ylab = '', main = 'Differentiation order index\nof maximal accessibility')
# axis(1,at =c(1,2), labels = c('CthPN\npeaks', 'CPN\npeaks'), col = NA,  line= 3)
# dev.off()

# # png('./output/mcatac/figs/diff_order_index_of_max_acc_cthpn_vs_cpnl23.png')
# par(mar = c(8,4,5,1), cex.main = 2.5, cex.axis = 2)

# vioplot(c(apply(mat_cpn_vs_cth_sm,1,which.max),
#             apply(mat_cth_vs_cpn_sm,1,which.max)) ~ 
#                   c(rep('cpn', nrow(mat_cpn_vs_cth_sm)), rep('cfupn', nrow(mat_cth_vs_cpn_sm))), xaxt = 'n', xlab = '',ylab = '', main = 'Differentiation order index\nof maximal accessibility')
# axis(1,at =c(1,2), labels = c('CthPN\npeaks', 'CPN\npeaks'), col = NA,  line= 3)
# # dev.off()

# # length(pks_diff_cfupn_clpl23)

# cor_tf_cl <- tgs_cor(t(a_legc_avg_cl[peak_clust_var,]), t(legc[tfs_hi,]), spearman = T)

# hcp_cor_mc <- hclust(dist(t(pltmt_cor)), method = 'ward.D2')

# ipc_legc_avg_atac_cl <- t(tgs_matrix_tapply(legc[,cust_mc_ord_st2[names(cust_mc_ord_st2) %in% ctis]], cutree(hcp_cor_mc, k = 3), mean))

# tail(ipc_legc_avg_atac_cl[order(ipc_legc_avg_atac_cl[,1] - rowMaxs(ipc_legc_avg_atac_cl[,2:3]), decreasing = T),], 20)
# tail(ipc_legc_avg_atac_cl[order(ipc_legc_avg_atac_cl[,2] - rowMaxs(ipc_legc_avg_atac_cl[,c(1,3)]), decreasing = T),], 20)
# tail(ipc_legc_avg_atac_cl[order(ipc_legc_avg_atac_cl[,3] - rowMaxs(ipc_legc_avg_atac_cl[,c(1,2)]), decreasing = T),], 20)

# a_legc_ct_max <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, max))


# a_legc_ct_min <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, min))

# a_legc_ct_sd <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, sd))

# options(repr.plot.width = 14)
# options(repr.plot.height = 14)

# ctis <- c('NSC','IPC', 'IPC_cyc', 'iCPN/CfuPN', 'iCPN_early', 'CthPN', 'CPN_L2-3')
# # ctis <- c('IPC', 'IPC_cyc')

# ipc_diff_peaks <- rownames(a_legc)[which(a_legc_ct_sd[,'IPC_cyc'] >= 0.65 & a_legc_ct_sd[,'IPC'] >= 0.65)]

# # ctis <- c('NSC','IPC', 'IPC_cyc', 'iCPN/CfuPN', 'iCPN_early')
# # ctis <- c('IPC', 'IPC_cyc')
# # pltmt_ipc_diff_peaks <- a_legc[rownames(a_legc)[which(a_legc_ct_sd[,'IPC_cyc'] >= 0.65 & a_legc_ct_sd[,'IPC'] >= 0.65)],
# #                     cust_mc_ord_st2[names(cust_mc_ord_st2) %in% ctis]]
# # hcp_ipc_dp <- hclust(dist(pltmt_ipc_diff_peaks), method = 'ward.D2')
# # hcp_ipc_dp_mc <- hclust(dist(t(pltmt_ipc_diff_peaks)), method = 'ward.D2')

# perc_prom <- tapply(rownames(a_legc), km_a_legc$cluster, function(x) length(which(x %in% prom_peaks$peak_name))/length(x))

# num_prom <- tapply(rownames(a_legc), km_a_legc$cluster, function(x) length(which(x %in% prom_peaks$peak_name)))

# # dir.create('./output/mcatac/figs/peak_clusters_across_metacells')

# options(repr.plot.width = 8)
# options(repr.plot.height = 8)

# # png('./output/mcatac/figs/peak_clust_var_filtering2.png')
# par(cex.lab = 2, mar = c(5,5,4,1), cex.axis = 1.5)
# plot(rowMeans(a_legc_avg_cl), (matrixStats::rowMaxs(a_legc_avg_cl) - matrixStats::rowMins(a_legc_avg_cl))/(2**rowMeans(a_legc_avg_cl)), 
#      pch = 16,
#      ylab = 'metric', xlab = 'log2 ATAC fraction', 
#      main = 'X = tgs_matrix_tapply(a_legc[var_peaks,],\nvar_peak_clusters, mean)\nmetric = (rowMaxs(X)-rowMins(X))/(2**rowMeans(X))')
# arb_met <- (matrixStats::rowMaxs(a_legc_avg_cl) - matrixStats::rowMins(a_legc_avg_cl))/(2**rowMeans(a_legc_avg_cl))
# abline(quantile(arb_met, 0.3), 0, col = 'red', lty = 2, lwd = 3)
# # dev.off()

# q_thresh <- quantile(pk_cvs[names(pk_cvs) %in% var_peaks], 0.9)
# q_thresh

# # png('./output/mcatac/figs/cumsum_num_prom_old_vs_new_peak_clust_var.png', h = 400, w = 800)
# par(mfrow = c(1,2), cex.main = 2, cex.lab = 2, mar = c(5,5,3,1))
# plot(1:length(num_prom), cumsum(sort(num_prom)), main = 'Old arbitrary metric', ylab = 'Cumulative num promoters')
# points(1:length(num_prom), cumsum(sort(num_prom)), pch = 16, col = ifelse(as.numeric(names(sort(num_prom))) %in% peak_clust_var, 'red', 'white'))
# legend('topleft', legend = c('variable peak cluster', 'non-variable peak cluster'), col = c('red', 'black'), pch = c(16,1), cex = 1.5)

# plot(1:length(num_prom), cumsum(sort(num_prom)), main = 'New manual selection', ylab = 'Cumulative num promoters')
# points(1:length(num_prom), cumsum(sort(num_prom)), pch = 16, col = ifelse(as.numeric(names(sort(num_prom))) %in% peak_clust_var_new, 'red', 'white'))
# legend('topleft', legend = c('variable peak cluster', 'non-variable peak cluster'), col = c('red', 'black'), pch = c(16,1), cex = 1.5)
# # dev.off()

# par(mfrow = c(1,2))
# plot(sort(arb_met), cumsum(num_prom[order(arb_met)]))
# points(sort(arb_met), cumsum(num_prom[order(arb_met)]), pch = 16, col = ifelse(as.numeric(names(sort(arb_met))) %in% peak_clust_var, 'red', 'white'))

# plot(sort(arb_met), cumsum(num_prom[order(arb_met)]))
# points(sort(arb_met), cumsum(num_prom[order(arb_met)]), pch = 16, col = ifelse(as.numeric(names(sort(arb_met))) %in% peak_clust_var_new, 'red', 'white'))

# par(mfrow = c(1,2))
# hist(rowSds(a_legc[prom_peaks$peak_name,]), 100)
# hist(rowMaxmins(a_legc[prom_peaks$peak_name,]), 100)

# v1 <- table(km_a_legc$cluster[match(prom_peaks$peak_name[rowSds(a_legc[prom_peaks$peak_name,]) >= 0.4], rownames(a_legc))])

# v1[as.character(1:80)[!(1:80 %in% as.numeric(names(v1)))]] <- 0

# v1 <- v1[order(as.numeric(names(v1)))]

# v2 <- table(km_a_legc$cluster[match(prom_peaks$peak_name, rownames(a_legc))])

# # quantile(arb_met)
# # hist(arb_met, 30)
# # lines(rep(quantile(arb_met, 1/3), 2), c(0, 10), col = 'red', lty = 2, lwd = 5)

# peak_clust_var <- which(arb_met >= quantile(arb_met, 0.35))

# pk_cvs <- rowSds(a_legc)**2/rowMeans(a_legc)**2

# peak_clust_var_new <- c(which(perc_prom < 0.1 | tapply(pk_cvs, km_a_legc$cluster, function(x) length(which(x >= q_thresh))/length(x)) > 0.1), 76)

# # peak_clust_var_new <- setdiff(as.numeric(rownames(a_legc_avg_cl)), c('45', '68', '72', '79'))
# peak_clust_var_new

# peak_clust_non_var <- setdiff(rownames(a_legc_avg_cl), peak_clust_var_new)

# save(km_a_legc, peak_clust_var, file = './output/mcatac/mmcortex_feat_peak_variable_peak_clusters.rda')

# load(file = './output/mcatac/mmcortex_feat_peak_variable_peak_clusters.rda')

# length(peak_clust_var)

# inds_glia <- which(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes'))
# inds_no_glia <- which(!(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes')))


# pltmt <- a_legc_avg_cl[peak_clust_var_new,cust_mc_ord_st2]


# hc_a_legc_avg_cl <- hclust(dist(pltmt - rowMeans(pltmt)), method = 'ward.D2')

# hcct <- cutree(hc_a_legc_avg_cl, k = 16)

# hcct_clr <- chameleon::distinct_colors(16)$name
# hcct_clr

# hcct_hc <- lapply(sort(unique(hcct)), function(cti) hclust(dist(pltmt[hcct == cti,]), method = 'ward.D2'))
# # ord_hcct <- c(7,8,9,5,6,3,2,1,4)
# # ord_hcct <- 1:16
# ord_hcct <- c(14,12,13,11,10,15,16,7,8,1,2,3,4,5,6,9)
# hcct2 <- unlist(sapply(ord_hcct, function(ctj) {if (length(hcct_hc[[ctj]]$order) < length(peak_clust_var_new)) {return(names(hcct)[which(hcct == ctj)][hcct_hc[[ctj]]$order])} else {return(which(hcct == ctj))}}))
# # ra <- as.data.frame(cbind(ecdf(arb_met)(arb_met[peak_clust_var_new]), as.character(hcct)))
# ra <- as.data.frame(as.character(hcct))
# # rownames(ra) <- 1:nrow(ra)
# rownames(ra) <- rownames(pltmt)
# # colnames(ra) <- c('arb_met', 'hcct_clr')
# colnames(ra) <- c('hcct_clr')
# # ra

# ann_colors[['arb_met']] <- setNames(colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100), seq(0,1, l=100))

# ann_colors[['hcct_clr']] <- setNames(c(hcct_clr, 'gray'), 1:16)

# pheatmap::pheatmap(pltmt[hcct2,], 
#                    # breaks= seq(-16.6,-14.75,l=100),
#                    annotation_row = ra,
# # pheatmap::pheatmap(pltmt[hc_a_legc_avg_cl$order,], breaks= seq(-16.6,-14.75,l=100), 
#                    # annotation_row = ra,
#                    show_colnames = F,
#                    color = colorRampPalette(c('white','orange', 'red','purple', 'black'))(100),
#                    annotation_col = col_annot, annotation_colors = ann_colors, 
#                    annotation_legend = F, cluster_cols = F, cluster_rows = F)

# pltmt_cor_tf_cl <- cor_tf_cl[hcct2,]
# hc_ctc <- hclust(dist(t(pltmt_cor_tf_cl)), method = 'ward.D')$order
# # p_cor_tf_cl <- pheatmap::pheatmap(pltmt_cor_tf_cl[,hc_ctc], fontsize = 14,
# #                    cluster_cols = F, cluster_rows = F, 
# #                    col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-1,1,l=100), 
# #                    treeheight_col = 0, treeheight_row = 0)
# # save_pheatmap(p_cor_tf_cl, './output/mcatac/figs/cor_tf_var_peak_cl_phm.png', h = 1600, w = 1600)

# pltmt_tf <- legc[tfs_hi[hc_ctc],cust_mc_ord_st2]
# # hc_ctc <- hclust(dist(t(pltmt_cor_tf_cl)), method = 'ward.D')$order
# p_tfs_hi_mc_rna <- pheatmap::pheatmap(pltmt_tf - rowMeans(pltmt_tf), gaps_col = tail(match(cust_st_ord2, names(cust_mc_ord_st2)) - 1, -1),
#                                       fontsize = 14, show_colnames = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F,
#                    cluster_cols = F, cluster_rows = F, 
#                    col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-4,4,l=100), 
#                    treeheight_col = 0, treeheight_row = 0)
# save_pheatmap(p_tfs_hi_mc_rna, './output/metacell_model//figs/tfs_hi_phm.png', h = 1600, w = 3000)

# # tfs_hi_motifs <- sapply(tfs_hi, function(tfi) grep(tfi, colnames(ie_mat), ign = T, v=T))
# # head(unlist(tfs_hi_motifs))

# # tfs_hi_motifs

# # length(unlist(tfs_hi_motifs))

# # library(matrixStats)

# # quantile(colMins(ie_mat))

# # ie_mat_avg_cl <- tgs_matrix_tapply(t(ie_mat), km_a_legc$cluster, mean)

# # ie_mat_avg_cl_tfs_hi <- ie_mat_avg_cl[,unlist(tfs_hi_motifs)]

# # dim(ie_mat_avg_cl_tfs_hi)

# # pltmt <- t(t(ie_mat_avg_cl_tfs_hi) - colMeans(ie_mat_avg_cl_tfs_hi))

# # pltmt_f <- pltmt[,colnames(pltmt)[colSds(pltmt) >= 0.6]]
# # hc_pf <- hclust(dist(t(pltmt_f)), method = 'ward.D')
# # ct_hc_pf <- cutree(hc_pf, 20)
# # pheatmap::pheatmap(pltmt_f[,hc_pf$order], cluster_cols = F, clustering_method = 'ward.D',col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100))
# # ct_hc_pf

# peak_clust_non_var <- setdiff(rownames(a_legc_avg_cl), peak_clust_var)

# options(repr.plot.width = 18)
# options(repr.plot.height = 8)

# png('./output/mcatac/figs/peak_cluster_sizes_and_frac_prom_barplots.png', w = 2200, h = 1200, res = 150)
# par(mfrow = c(2,2), mar = c(0,8,4,3), cex.axis = 2, cex.lab = 2, cex.main = 3, las = 2)
# mxi <- max(km_a_legc$size)
# barplot(sort(km_a_legc$size[peak_clust_var]), xaxt = 'n' ,main = 'Variable peak clusters', ylab = '', ylim = c(0,mxi))
# title(ylab = 'Count', line = 6)
# par(mar = c(0,1,4,0))
# barplot(sort(km_a_legc$size[peak_clust_non_var]), xaxt = 'n',main = 'Constitutive peak clusters', ylim = c(0,mxi))
# par(mar = c(4,8,2,3))
# barplot(perc_prom[peak_clust_var[order(km_a_legc$size[peak_clust_var])]], ylim = c(0,1), ylab = '', cex.names = 0.5)
# # axis(1, at = 1:length(peak_clust_var), labels = peak_clust_var[order(km_a_legc$size[peak_clust_var])], cex.axis = 1)
# title(ylab = 'Fraction TSS-proximal', line = 6)
# par(mar = c(4,1,2,0))
# barplot(perc_prom[peak_clust_non_var[order(km_a_legc$size[peak_clust_non_var])]], ylim = c(0,1), cex.names = 1.5)
# dev.off()

# png('./output/mcatac/figs/peak_cluster_sizes_and_frac_prom_barplots.png', w = 1800, h = 1200, res = 100)
# par(mfrow = c(2,1), mar = c(0,8,4,3), cex.axis = 2, cex.lab = 2, cex.main = 3, las = 2)
# mxi <- max(km_a_legc$size)
# x <- barplot(c(sort(km_a_legc$size[peak_clust_var]), 0, sort(km_a_legc$size[peak_clust_non_var])), cex.names = 1,
#              ylab = '', ylim = c(0,mxi))
# grid(lwd = 2, lty = 2)
# lines(rep(x[[length(peak_clust_var)]]+ 1.25, 2), c(0,1e+4), lwd = 2, lty = 2)
# text(x[[floor(length(peak_clust_var)/2)]], y = 6e+3, labels = c('Variable peak clusters'), cex = 2)
# text(x[[floor(length(peak_clust_var))*1.15]], y = 6e+3, labels = c('Constitutive\npeak clusters'), cex = 2)
# title(ylab = 'Cluster size', line = 6)

# par(mar = c(0,1,4,0))
# # barplot(, xaxt = 'n',main = 'Constitutive peak clusters', ylim = c(0,mxi))
# par(mar = c(4,8,2,3))
# barplot(c(sort(perc_prom[peak_clust_var]) ,0, sort(perc_prom[peak_clust_non_var])),
#         ylim = c(0,1), ylab = '', cex.names = 1)
# grid(lwd = 2, lty = 2)
# lines(rep(x[[length(peak_clust_var)]]+ 1.25, 2), c(0,1e+4), lwd = 2, lty = 2)
# text(x[[floor(length(peak_clust_var)/2)]], y = 0.8, labels = c('Variable peak clusters'), cex = 2)
# text(x[[floor(length(peak_clust_var))*1.15]], y = 0.8, labels = c('Constitutive\npeak clusters'), cex = 2)

# # axis(1, at = 1:length(peak_clust_var), labels = peak_clust_var[order(km_a_legc$size[peak_clust_var])], cex.axis = 1)
# title(ylab = 'Fraction TSS-proximal', line = 6)
# # par(mar = c(4,1,2,0))
# # barplot(, ylim = c(0,1), cex.names = 1.5)
# dev.off()

# png('./output/mcatac/figs/mmcortex_famc_legc_w_motif_q98_lfc_glia.png', w = 380, h =1270)
# draw(ch_glia)
# dev.off()

# plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6))
# {
#   if (!is.null(fig_fn)) {
# #     .plot_start(fig_fn, 400, 400)
#       png(fig_fn, 400, 400)
#   }
#   plot.new()
#   plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
#   rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
#   rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

# #   if (is.null(show_vals_ind)) {
# #     show_vals_ind = rep(T, length(cols))
# #   }
#   text(19, show_vals_ind,cex = 2, labels=round(vals[show_vals_ind], 3), pos=4)
# #   text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

#   if (!is.null(fig_fn)) {
#     dev.off()
#   }
# }
# shades = adjustcolor(colorRampPalette(c('darkblue','blue3', 'white', 'red','red4'))(100), alpha.f = 0.75)
# min_val = -3
# max_val = 3
# plot_color_bar(seq(min_val, max_val,l=101), shades, fig_fn = './output/mcatac/figs/atac_heatmap_motif_rgq_color_bar.png')

# pltmt2 <- a_legc_avg_cl_prom[,cust_mc_ord_st2]

# col_ha1 <- HeatmapAnnotation(cell_type = anno_simple(col_annot$cell_type[match(cust_mc_ord_st2, mcmd$metacell)], 
#                                                      col = ann_colors[['cell_type']],
#                                                      height =unit(1, 'cm')), 
#                              mean_day = anno_lines(x = col_annot$mean_day[match(cust_mc_ord_st2, mcmd$metacell)], axis_param = list(gp = gpar(fontsize = 16)),
#                                                     # col = circlize::colorRamp2(seq(13,18,1), c('red', 'orange', 'yellow', 'green2', 'blue', 'purple')), 
#                                                     height =unit(1, 'cm')), 
#                              annotation_name_gp = gpar(fontsize = 18),
#                              show_legend = F)

# row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_prom_a_legc$size), ylim = c(800,3000), gp = gpar(fill = 'black',fontsize = 18), 
#                                                         axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90)), 
#                             # frac_prom = anno_barplot(as.numeric(t(perc_prom[peak_clust_var[hcct2]])), ylim = c(0,1), 
#                             # gp = gpar(fill = 'black',fontsize = 18), 
#                             # axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90),
#                             annotation_name_gp = gpar(fontsize = 18),
#                             annotation_name_rot = 0,
#                             # cl2 = anno_simple(x = ct_hc_a_legc[ord1], col = setNames(peak_cl2_color_key$color[match(ct_hc_a_legc[ord1], peak_cl2_color_key$cl2)], ct_hc_a_legc[ord1])),
#                             annotation_name_offset = unit(3, 'cm'),        
#                             which = 'row',
#                             width = unit(5, 'cm')
#                            )


# ch2 <- ComplexHeatmap::Heatmap(matrix = pltmt2[,], name = 'log2\nfraction\nATAC', 
#                               col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
#                               column_split = factor(names(cust_mc_ord_st2), levels = cust_st_ord2), column_gap = unit(2, 'mm'),
#                               column_title_gp = gpar(fontsize = 0),
#                               top_annotation = col_ha1, 
#                               # bottom_annotation = col_ha2,
#                               show_heatmap_legend = T,
#                               show_column_names = F,
#                                show_row_dend = FALSE,
#                               column_title_rot = 90,
#                               # right_annotation = motif_ha,
#                               heatmap_legend_param = list(legend_height = unit(5, 'in'), legend_width = unit(5, 'in'), labels_gp = gpar(fontsize = 16)),
#                               # heatmap_width = unit(92, 'cm'), heatmap_height = unit(25, 'cm'),
#                                heatmap_width = unit(45, 'cm'), heatmap_height = unit(8, 'cm'),
#                         left_annotation = row_ha,
#                               row_names_gp = gpar(fontsize = 16),
#                         cluster_columns = F, cluster_rows = T,
#                               clustering_method_rows = 'ward.D2')

# # png('./output/mcatac/figs/mmcortex_famc_legc_non_var_peaks_cluster_size_frac_prom_annot.png', w = 2400, h = 1000)
# png('./output/mcatac/figs/mmcortex_famc_legc_non_var_peaks_cluster_size_frac_prom_annot_test.png', w = 1450, h = 350)
# draw(ch2)
# dev.off()

# # png('./output/mcatac/figs/mmcortex_famc_legc_non_var_peaks_cluster_size_frac_prom_annot.png', w = 2400, h = 1000)
# pdf('./output/mcatac/figs/mmcortex_famc_legc_non_var_peaks_cluster_size_frac_prom_annot_test.pdf', w = 1450/.71e+2, h = 350/.71e+2)
# draw(ch2)
# dev.off()

# set.seed(1337)
# samp_mc_ct <- sapply(tail(cust_st_ord, -2), function(cti) sample(which(names(cust_mc_ord_st) == cti), 1))

# pheatmap::pheatmap(pltmt[hcct2,c(inds_glia, samp_mc_ct)], cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, fontsize = 16)

# cor_pltmt <- tgs_cor(pltmt)

# # pheatmap::pheatmap(cor_pltmt[inds_glia,order(match(as.numeric(colnames(pltmt)), cust_mc_ord_st2))], annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F)
# pheatmap::pheatmap(cor_pltmt[inds_glia,which(as.numeric(colnames(pltmt)) %in% which(mcmd$cell_type %in% c('Astrocytes', 'OPCs', 'CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')))], annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, fontsize = 16)

# options(repr.plot.width = 18)
# options(repr.plot.height = 9)

# min_a_legc <- min(a_legc)
# var_peaks_n <- var_peaks[which(apply(a_legc[var_peaks,602:603], 1, function(x) all(x > min_a_legc)))]

# cor_mat <- tgs_cor(cbind(a_legc[var_peaks_n,which(mcmd$cell_type %in% c('Astrocytes', 'OPCs'))], a_legc_avg_ct[var_peaks_n,]), spearman = T)

# pheatmap::pheatmap(cor_mat[1:9,10:ncol(cor_mat)])

# par(mfrow = c(2,4))
# ttt <- sapply(c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3'), function(cti) {smoothScatter(a_legc[var_peaks,c(602)], a_legc[var_peaks,sample(x = which(mcmd$cell_type == cti), size = 1)], cex = .26, pch = 16); abline(0,1,col = 'red')})
# ttt <- sapply(c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3'), function(cti) {smoothScatter(a_legc[var_peaks,c(603)], a_legc[var_peaks,sample(x = which(mcmd$cell_type == cti), size = 1)], cex = .26, pch = 16); abline(0,1,col = 'red')})
# # ttt <- sapply(plot(as.data.frame(cbind(a_legc[var_peaks,c(602, 603)], a_legc_avg_ct[var_peaks,c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')])), cex = .16, pch = 16)

# # pltmt <- a_legc_avg_cl[peak_clust_var,cust_mc_ord_st]
# # ord <- hclust(tgs_dist(pltmt))$order
# var_peaks_ordered <- do.call('c', lapply(peak_clust_var[hcct2], function(x) rownames(a_legc)[which(km_a_legc$cluster == x)]))
# names(var_peaks_ordered) <- do.call('c', lapply(peak_clust_var[hcct2], function(x) rep(x, length(which(km_a_legc$cluster == x)))))

# a_legc_avg_cl_ct <- t(tgs_matrix_tapply(a_legc_avg_cl, mcmd$cell_type, mean))

# clusters_peaking_in_neurons <- rownames(a_legc_avg_cl_ct[peak_clust_var,])[apply(a_legc_avg_cl_ct[peak_clust_var,], 1, function(x) which.max(x) %in% match(cts, colnames(a_legc_avg_cl_ct)))]
# clusters_peaking_in_neurons                                                                                 
# #                                                                                  lapply(cts, function(ct) {
# #     rownames(a_legc_avg_cl_ct[peak_clust_var,])[apply(a_legc_avg_cl_ct[peak_clust_var,], 1, function(x) which.max(x) %in% cts)]
# # })
# # names(clusters_peaking_in_neurons) <- cts
# # clusters_peaking_in_neurons <- unlist(clusters_peaking_in_neurons)



# clusters_to_plot <- c(29,62,22,1)

# # clp <- rownames(a_legc_var_avg_pk_cl)[as.numeric(rownames(a_legc_var_avg_pk_cl)) >= 25]
# # clp <- rownames(a_legc)[as.numeric(rownames(a_legc)) >= 25]
# df_clp <- apply(expand.grid(clusters_to_plot, clusters_to_plot), 2, as.character)
# head(df_clp)

# options(repr.plot.width = 25)
# options(repr.plot.height = 25)

# png('./output/mcatac/figs/famc_averaged_peak_clusters_pairwise.png', h = 3000, w = 4500, res = 200)
# par(mfrow = c(4,4), cex.lab = 3, cex.main = 3, mar = c(6,6,1,1), cex.axis = 3)
# # par(mfrow = c(7,7), cex.lab = 3, cex.main = 3, mar = c(6,6,1,1))
# # tbnd <- sapply(1:16, function(n) {
# tbnd <- sapply(1:nrow(df_clp), function(n) {
#     plot(a_legc_avg_cl[df_clp[n,1],], 
#          a_legc_avg_cl[df_clp[n,2],], 
#          col = mcmd$color, 
#          pch = 16, cex = 1.5, 
#          xlab = paste0('cluster ', df_clp[n,1]),
#         ylab = paste0('cluster ', df_clp[n,2]),
# #         main = glue::glue('Clusters {df_clp[n,1]} vs {df_clp[n,2]}')
#         )
# })
# dev.off()

# pltmt <- t(a_legc[var_peaks_ordered, rev(cust_mc_ord_st)])

# p_var_peaks_all <- pheatmap::pheatmap(pltmt, silent = T,
#                                       annotation_colors = ann_colors, annotation_row = col_annot, annotation_legend = F,
#                                       show_colnames = F,  color = clrmp, breaks = brks,
#                         cluster_cols = F, cluster_rows = F, show_rownames = F)

# save_pheatmap(p_var_peaks_all, './output/mcatac/figs/famc_a_legc_var_peaks_all.png', h = 2000, w = 2000, res = 200)

# l23_peaks <- get_genes_specific_to_mcs(a_legc, mc_pos = which(mcmd$cell_type %in% c('CPN_L2-3')), mc_neg = which(mcmd$cell_type %in% c('CPN_L5_6')))
# l56_peaks <- get_genes_specific_to_mcs(a_legc, mc_pos = which(mcmd$cell_type %in% c('CPN_L5_6')), mc_neg = which(mcmd$cell_type %in% c('CPN_L2-3')))
# scpn_peaks <- get_genes_specific_to_mcs(a_legc, mc_pos = which(mcmd$cell_type == 'SCPN'), mc_neg = which(mcmd$cell_type %in% c('CthPN')))
# cthpn_peaks <- get_genes_specific_to_mcs(a_legc, mc_pos = which(mcmd$cell_type == 'CthPN'), mc_neg = which(mcmd$cell_type %in% c('SCPN')))
# mat_neuro_peaks <- unique(unlist(lapply(list(l23_peaks, l56_peaks, scpn_peaks, cthpn_peaks), function(x) head(names(x), 50))))

# tads <- gintervals.load('intervs.global.tad_names')

# # mat_neuro_peaks_all <- unique(unlist(lapply(list(l23_peaks, l56_peaks, scpn_peaks, cthpn_peaks), function(x) head(names(x), 1000))))

# options(repr.plot.height = 15)

# # pltmt <- a_legc_avg_cl[peak_clust_var,cust_mc_ord_st]
# # ord <- hclust(tgs_dist(pltmt))$order
# non_var_peaks_ordered <- do.call('c', lapply(peak_clust_non_var[ord2], function(x) rownames(a_legc)[which(km_a_legc$cluster == x)]))
# names(non_var_peaks_ordered) <- do.call('c', lapply(peak_clust_non_var[ord2], function(x) rep(x, length(which(km_a_legc$cluster == x)))))

# pltmt <- t(a_legc[non_var_peaks_ordered, rev(cust_mc_ord_st)])

# p_non_var_peaks_all <- pheatmap::pheatmap(pltmt, silent = T,
#                                       annotation_colors = ann_colors, annotation_row = col_annot, annotation_legend = F,
#                                       show_colnames = F,  color = clrmp, breaks = brks,
#                         cluster_cols = F, cluster_rows = F, show_rownames = F)

# save_pheatmap(p_non_var_peaks_all, './output/mcatac/figs/famc_a_legc_non_var_peaks_all.png', h = 2000, w = 2000, res = 200)

# rmal <- setNames(rowMeans(a_legc), rownames(a_legc))
# rsal <- setNames(rowSds(a_legc), rownames(a_legc))
# const_peaks <- rownames(a_legc)[rmal >= quantile(rmal, 0.95) & rsal < 0.6]

# png('./output/mcatac/figs/constitutive_peak_selection.png', h = 600, w  = 600)
# par(cex.axis = 1.5, cex.lab = 2, cex.main = 2, mar = c(4,6,2,1))
# plot(rmal, rsal, xlab = 'rowMeans(a_legc)', ylab = 'rowSds(a_legc)', main = 'Constitutive peak selection')
# points(rmal[const_peaks], rsal[const_peaks], col = 'red', pch = 16)
# dev.off()

# png('./output/mcatac/figs/constitutive_peaks_in_cell_types_vs_NSC.png', h = 1200, w  = 1600)
# par(mfrow = c(3,4), cex.axis = 2, mar = c(3,7,1,1), cex.lab = 3)
# vvv <- sapply(cust_st_ord[cust_st_ord != 'NSC'], function(cti) {
#     plot(a_legc_avg_ct[const_peaks,'NSC'], a_legc_avg_ct[const_peaks,cti], cex = 0.5, pch = 16, col = color_key$color[match(cti, color_key$cell_type)], ylab = cti, xlab = '')
#     abline(0,1,col='red')
# })
# dev.off()

# # motifs_to_take <- c('HOMER.Sox10', 'JOLMA.HOXD3_mono_DBD','JASPAR.NEUROG2.MA1642.1', 
# #                     'HOMER.Eomes','JOLMA.POU3F2_mono_DBD_1')
# motifs_to_take <- c('JASPAR.SOX2', 'JOLMA.EMX1_mono_DBD','JASPAR.NEUROG2.MA1642.1', 
#                     'HOMER.Eomes','JOLMA.POU3F2_mono_DBD_1')
#                     # names(head(sort(diff_energies, decreasing = T), 15)))

# intervs_energy <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')

# intervs_energy_new <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_prego_motif_energy.rds')

# # ie_mat <- as.matrix(subset(intervs_energy, select = -c(chrom, start, end, peak_name)))
# # rownames(ie_mat) <- peak_names(intervs_energy[,c('chrom', 'start', 'end')], tad_based =F)
# ie_mat <- as.matrix(subset(intervs_energy_new, select = -c(chrom, start, end, peak_name, mmcortex.marginal,	intervalID,	peak_name, peak_name_ntb)))
# # ie_mat <- as.matrix(intervs_energy_new[,colnames(prego_ie_mat)])
# rownames(ie_mat) <- intervs_energy_new$peak_name

# ie_mat_avg_cl_nsc_ipc <- tgs_matrix_tapply(t(ie_mat), km_a_legc_ipc_sm_pcu$cluster, mean, na.rm = T)

# ie_mat_avg_cl <- tgs_matrix_tapply(t(ie_mat[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean, na.rm = T)

# pheatmap::pheatmap(apply(ie_mat_avg_cl[,zbtb_sall_id4_motif_names], 2, function(x) {y <- x - max(x); y[y < -10] <- -10; return(y)}), col = clrmp_abs)

# ns_rg <- colSds(ie_mat_avg_cl)/abs(colMeans(ie_mat_avg_cl))

# ns_nsc_ipc <- colSds(ie_mat_avg_cl_nsc_ipc)/abs(colMeans(ie_mat_avg_cl_nsc_ipc))

# library(prego)

# amd <- prego::all_motif_datasets()

# mtfs <- lapply(colnames(ie_mat_avg_cl)[ns_rg >= quantile(ns_rg, 0.99)], function(x) dplyr::filter(amd, motif == x))

# library(matrixStats)

# tfs_in <- intersect(tfs, rownames(legc)[rowMaxs(legc) >= -13.5])

# write(tfs_in, file = './output/metacell_model/tfs_in.txt')

# length(tfs_in)

# head(tfs_in)

# amou <- unique(amd[,c('motif', 'motif_orig')])

# motifs_tfs_lst <- sapply(tfs_in, function(x) amou$motif[grep(x, amou$motif_orig, ign = T, v=F)])
# motifs_tfs_in <- intersect(unique(unlist(motifs_tfs_lst)), colnames(ie_mat))

# clrmp <- colorRampPalette(c('blue3', 'white', 'red3'))(100)
# brks <- seq(-16.6,-15.5, l = 100)
# nsc_clvls <- clrmp[1+round(99*(a_legc_avg_cl_ct[,'NSC'] + 16.6)/(-15+16.6))]
# ipc_clvls <- clrmp[1+round(99*(rowMeans(a_legc_avg_cl_ct[,c('IPC', 'IPC_cyc')]) + 16.6)/(-15+16.6))]
# y <- rowMeans(a_legc_avg_cl_ct[,c('IPC', 'IPC_cyc')]) - a_legc_avg_cl_ct[,'NSC']
# delta_ipc_nsc_clvls <- clrmp[1+round(99*(y - min(y))/(max(y) - min(y)))]
# cpn_clvls <- clrmp[1+round(99*(a_legc_avg_cl_ct[,'CPN_L2-3']) + 16.6)/(-15+16.6)]
# cfupn_clvls <- clrmp[1+round(99*(a_legc_avg_cl_ct[,'CthPN']) + 16.6)/(-15+16.6)]

# quantile(rowMeans(a_legc_avg_cl_ct[,c('IPC', 'IPC_cyc')]) - a_legc_avg_cl_ct[,'NSC'])

# par(mfrow = c(2,3))
# vvv <- sapply(motifs_to_take2, function(mtfi) {plot(cpg_cl[peak_clust_var], ie_mat_rgq_avg_cl[peak_clust_var,mtfi], col = delta_ipc_nsc_clvls[peak_clust_var], pch = 16, cex = 3, main = mtfi, ylab = 'motif energy', xlab = 'CpG quantile'); abline(0,-1)})

# par(mfrow = c(2,3))
# vvv <- sapply(motifs_to_take2, function(mtfi) {plot(cpg_cl[peak_clust_var], ie_mat_rgq_avg_cl[peak_clust_var,mtfi], col = cpn_clvls[peak_clust_var], pch = 16, cex = 3, main = mtfi, ylab = 'motif energy', xlab = 'CpG quantile'); abline(0,-1)})

# par(mfrow = c(2,3))
# vvv <- sapply(motifs_to_take2, function(mtfi) {plot(cpg_cl[peak_clust_var], ie_mat_rgq_avg_cl[peak_clust_var,mtfi], col = cfupn_clvls[peak_clust_var], pch = 16, cex = 3, main = mtfi, ylab = 'motif energy', xlab = 'CpG quantile'); abline(0,-1)})

# par(mfrow = c(2,3))
# vvv <- sapply(motifs_to_take2, function(mtfi) {plot(cpg_cl[peak_clust_var], ie_mat_rgq_avg_cl[peak_clust_var,mtfi], col = nsc_clvls[peak_clust_var], pch = 16, cex = 3, main = mtfi, ylab = 'motif energy', xlab = 'CpG quantile'); abline(0,-1)})

# par(mfrow = c(2,3))
# vvv <- sapply(motifs_to_take2, function(mtfi) {plot(cpg_cl[peak_clust_var], ie_mat_rgq_avg_cl[peak_clust_var,mtfi], col = ipc_clvls[peak_clust_var], pch = 16, cex = 3, main = mtfi, ylab = 'motif energy', xlab = 'CpG quantile'); abline(0,-1)})

# mod_in_lst <- list('NSC_' = c(29,36,31), 'IPC_' = c(64,19,16), 'CPN_' = c(4,2,1), 'Neuron_' = c(14,10,11), 'CfuPN_' = c(22,23))

# mod_in_lst_u <- unlist(mod_in_lst)
# mod_in_lst_u

# # bxp_df$cluster <- names(mod_in_lst_u)[match(bxp_df$cluster, mod_in_lst_u)]

# png('./output/sequence_modeling/figs/atac_and_motif_energies_select_peak_clusters_boxplots_take2.png', h = 1800, w = 730)
# par(mfrow = c(length(mod_in), 2), cex.lab = 3, cex.axis = 2)
# mari <- c(0.5,0.5,0.5,0.5)
# par(mar = mari, las = 2)
# vvv <- lapply(mod_in, function(clj) {
#     inds_cl <- which(bxp_df$cluster == clj)
#     if (clj == tail(mod_in, 1)) {
#         mari[[1]] <- 12
#         xaxti <- 's'
#     } else {
#         xaxti <- 'n'
#     }
#     mari[[2]] <- 8
#     mari[[4]] <- 0.5
#     par(mar = mari)
#     yaxti <- 's'
#     boxplot(bxp_df[inds_cl,grep('ATAC', colnames(bxp_df))], bty = 'n', xaxt = xaxti, yaxt = yaxti)
#     title(ylab = clj, line=  6)
#     mari[[2]] <- 0.5
#     mari[[4]] <- 5
#     par(mar = mari)
#     yaxti <- 'n'
#     boxplot(bxp_df[inds_cl,grep('ATAC|clu', colnames(bxp_df), inv=T)], bty = 'n', xaxt = xaxti, yaxt = yaxti)
#     axis(4, at = seq(0,1,l=5))
# })
# dev.off()

# options(repr.plot.height = 25)

# png('./output/sequence_modeling/figs/atac_and_motif_energies_select_peak_clusters_boxplots.png', h = 1250, w = 750)
# par(mfrow = c(9,1), cex.lab = 2, mar = c(4,8,.5,1), cex.axis = 2)
#     ttt <- sapply(colnames(bxp_df)[-5], function(cni) boxplot(bxp_df[,cni] ~ bxp_df$cluster, ylab = cni, col = clrsn, names = mod_in, xlab = ''))
# dev.off()

# library('vioplot')

# options(repr.plot.height = 20)
# options(repr.plot.width = 8)

# line_df <- cbind(c(3.5, 6.5, 9.5, 12.5), c(3.5, 6.5, 9.5, 12.5), rep(-100, 4), rep(100, 4))
# line_df

# clrmp <- colorRampPalette(c('blue3', 'white', 'red3'))(100)

# mari <- c(0.5,15,.25,.25)
# png('./output/sequence_modeling/figs/atac_and_motif_energies_select_peak_clusters_vioplots_mtt2.png', h = 2000, w = 900)
# par(mfrow = c(12,1), cex.lab = 4, cex.axis = 3)
#     ttt <- sapply(colnames(bxp_df)[-5], function(cni) {
#         if (cni == tail(colnames(bxp_df), 1)) {mari[[1]] <- 12; xaxti = 'n'} else {mari[[1]] <- 0.5; xaxti = 'n'}
#         par(mar = mari, las = 2)
#         medians <- tapply(bxp_df[,cni], bxp_df$cluster, median)
#         med_lin <- (medians - min(medians))/(max(medians) - min(medians))
#         clvls <- clrmp[1+round(99*med_lin)]
#         vioplot::vioplot(bxp_df[,cni] ~ bxp_df$cluster, col = clvls, names = mod_in, xlab = '', xaxt = xaxti, ylab = '');
#         uuu <- sapply(1:nrow(line_df), function(n) lines(line_df[n,1:2], line_df[n,3:4], lty = 1, lwd = 6))
#         title(ylab = gsub(' ', '\n', cni), line = 7)
#         if (cni == tail(colnames(bxp_df), 1)) {axis(1,at = 1:length(mod_in),labels = names(mod_in_lst_u)[match(mod_in, mod_in_lst_u)])}
        
#     })
# dev.off()

# png('./output/sequence_modeling/figs/atac_and_motif_energies_select_peak_clusters_boxplots.png', h = 1250, w = 750)
# par(mfrow = c(9,1), cex.lab = 2, mar = c(4,8,.5,1), cex.axis = 2)
#     ttt <- sapply(colnames(bxp_df)[-5], function(cni) boxplot(bxp_df[,cni] ~ bxp_df$cluster, ylab = cni, names = mod_in, xlab = ''))
# dev.off()

# top_1k_peaks_per_motif <- lapply(motifs_to_take, function(mi) head(rownames(ie_mat)[order(ie_mat[,mi], decreasing = T)], 500))

# names(top_1k_peaks_per_motif) <- motifs_to_take

# top_1k_peaks_per_motif_2 <- lapply(motifs_to_take2, function(mi) head(rownames(ie_mat)[order(ie_mat[,mi], decreasing = T)], 500))

# names(top_1k_peaks_per_motif_2) <- motifs_to_take2

# options(repr.plot.height = 20)
# options(repr.plot.width = 8)

# png('./output/sequence_modeling/figs/top_1k_energies_peaks_2_ATAC_in_NSC_by_t.png', h = 1200, w = 500)
# par(mfrow = c(6,1), cex.lab = 2, mar = c(3,6,1,1), cex.axis = 2)
# fff <- lapply(names(top_1k_peaks_per_motif_2), function(x) boxplot(log2(1e-5 + egc_by_day_n[top_1k_peaks_per_motif_2[[x]],]), ylim = c(-16.6, -14.5), ylab = x))
# dev.off()

# png('./output/sequence_modeling/figs/top_1k_energies_peaks_ATAC_in_NSC_by_t.png', h = 1200, w = 500)
# par(mfrow = c(5,1), cex.lab = 2, mar = c(3,6,1,1), cex.axis = 2)
# fff <- lapply(names(top_1k_peaks_per_motif), function(x) boxplot(log2(1e-5 + egc_by_day_n[top_1k_peaks_per_motif[[x]],]), ylim = c(-16.6, -14.5), ylab = x))
# dev.off()



# intervs_energy_new <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_prego_motif_energy.rds')

# # ie_mat <- as.matrix(subset(intervs_energy, select = -c(chrom, start, end, peak_name)))
# # rownames(ie_mat) <- peak_names(intervs_energy[,c('chrom', 'start', 'end')], tad_based =F)
# ie_mat <- as.matrix(subset(intervs_energy_new, select = -c(chrom, start, end, peak_name, mmcortex.marginal,	intervalID,	peak_name, peak_name_ntb)))
# # ie_mat <- as.matrix(intervs_energy_new[,colnames(prego_ie_mat)])
# rownames(ie_mat) <- intervs_energy_new$peak_name

# raq98 <- unlist(plyr::llply(colnames(ie_mat), function(x) quantile(ie_mat[,x], probs = 0.98), .parallel = T))

# names(raq98) <- colnames(ie_mat)

# ra_98_bin_int <- t(plyr::laply(1:length(raq98), function(i) as.numeric(ie_mat[,i] >= raq98[[i]]), .parallel = T))
# save(ra_98_bin_int, file='./output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda')
# colnames(ra_98_bin_int) <- colnames(ie_mat)
# rownames(ra_98_bin_int) <- rownames(ie_mat)

# load('./output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda')
# # ra_98_bin_int <- ra98bin_sp

# # ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int)), km_a_legc$cluster, sum)

# ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int[dist_peaks$peak_name,])), km_enh_a_legc$cluster, sum)

# # ra_98_lfc_int <- log2(1e-2+ra_98_sum_clust_int/(0.02*as.numeric(km_a_legc$size)))
# ra_98_lfc_int <- log2(1e-2+ra_98_sum_clust_int/(0.02*as.numeric(km_enh_a_legc$size)))
# colnames(ra_98_lfc_int) <- colnames(ra_98_bin_int)









# amd <- prego::all_motif_datasets()



# # save(ra98bin_sp, file='./output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda')

# load('./output/sequence_modeling/feat_peak_all_motif_q_98_lfc.rda')

# table(apply(ra_98_bin_int[,grep('CTCF', colnames(ra_98_bin_int))], 1, function(x) length(which(x == 1))))

# ctcf_sites <- rownames(ra_98_bin_int)[apply(ra_98_bin_int[,grep('CTCF', colnames(ra_98_bin_int))], 1, function(x) length(which(x == 1)) >= 2)]

# tbl_ctcf <- sort(table(km_a_legc$cluster[which(rownames(a_legc) %in% ctcf_sites)]))
# tbl_ctcf_norm <- tbl_ctcf/as.numeric(km_a_legc$size[as.numeric(names(tbl_ctcf))])

# round(sort(tbl_ctcf_norm), 3)

# length(round(sort(tbl_ctcf_norm), 3))

# length(intersect(as.numeric(names(tbl_ctcf_norm)), peak_clust_non_var))

# length(peak_clust_non_var)

# pheatmap::pheatmap(cbind(tbl_ctcf, ra_98_lfc_int[names(tbl_ctcf),grep('ctcf', colnames(ra_98_lfc_int), ign = T)]))

# ctcf_site_coords <- dplyr::filter(mcp, peak_name %in% ctcf_sites)

# head(tads)

# nrow(ctcf_site_coords)

# tad_borders <- do.call('rbind', lapply(unique(as.character(.misha$ALLGENOME[[1]]$chrom)), function(chri) {
#     starti <-  .misha$ALLGENOME[[1]]$end[.misha$ALLGENOME[[1]]$chrom == chri]
#     dplyr::filter(tads, chrom == chri) %>% select(chrom, start, end) %>% mutate(end = start + 1) %>% bind_rows(., as.data.frame(list(chrom = chri, start = starti- 1, end = starti)))
# }))

# nei_ctcf_tad <- gintervals.neighbors(ctcf_site_coords, tad_borders, maxdist = 1.5e+6, mindist = 0, maxneighbors = 1)

# peak_clust_non_var <- setdiff(1:80, peak_clust_var)

# nei_rand_peaks_tad <- gintervals.neighbors(dplyr::sample_n(mcp[mcp$peak_name %in% setdiff(rownames(a_legc)[km_a_legc$cluster %in% peak_clust_non_var], ctcf_sites),], size = nrow(ctcf_site_coords)), tad_borders, maxdist = 1.5e+6, mindist = 0, maxneighbors = 1)

# nrow(nei_ctcf_tad)

# nrow(nei_rand_peaks_tad)

# length(setdiff(rownames(a_legc)[km_a_legc$cluster %in% peak_clust_non_var], ctcf_sites))

# quantile(nei_ctcf_tad$dist, (0:20)/20)

# quantile(nei_rand_peaks_tad$dist, (0:20)/20)

# plot(ecdf(nei_ctcf_tad$dist), do.points = F)
# plot(ecdf(nei_rand_peaks_tad$dist), do.points = F, col = 'blue', add = T)

# vioplot::vioplot(c(nei_ctcf_tad$dist, nei_rand_peaks_tad$dist) ~ c(rep('ctcf', nrow(nei_ctcf_tad)), rep('rand', nrow(nei_rand_peaks_tad))), ylim = c(0,1e+6))

# ks.test(nei_ctcf_tad$dist, nei_rand_peaks_tad$dist)

# # tfs_hi_in <- unique(unlist(sapply(tfs_in, function(tfi) grep(paste0('\\.', tfi, '[_\\.$]'),colnames(ra_98_lfc_int), v=T, ign = T))))
# # tfs_hi_in <- unique(unlist(sapply(tfs_in, function(tfi) grep(paste0('\\.', tfi, '[_\\.$]?[^0-9]'),colnames(ra_98_lfc_int), v=T, ign = T))))

# grep('mga', tfs_hi_in, ign = T, v=T)

# grep('mga', rownames(legc), ign = T, v=T)

# tfs_orig <- unlist(purrr::map(stringr::str_split(tolower(amd$motif_orig), '_'), 1))

# # tfs_hi_in <- unique(unlist(sapply(tfs_in, function(tfi) grep(paste0('\\.', tfi, '[_\\.$]'),colnames(ra_98_lfc_int), v=T, ign = T))))
# tfs_hi_in <- unique(amd$motif[tfs_orig %in% tolower(tfs_in)])

# # tfs_hi_in <- amd$motif[amd$motif_orig]

# # tfs_hi_in <- unique(unlist(sapply(tfs_in, function(tfi) grep(paste0('\\.', tfi, '[_\\.$]'),colnames(ra_98_lfc_int), v=T, ign = T))))
# tfs_hi_out <- setdiff(unique(unlist(sapply(tfs_in, function(tfi) grep(tfi, colnames(ra_98_lfc_int), v=T, ign = T)))), tfs_hi_in)

# length(tfs_hi_in)

# length(tfs_hi_out)

# tfs_hi_in_var <- tfs_hi_in[apply(ra_98_lfc_int[,tfs_hi_in], 2, function(x) any(x >= log2(3)))]

# tfs_hi_in_var_f <- unlist(purrr::map(stringr::str_split(tfs_hi_in_var, '\\.'), 2))
# sort(tfs_hi_in_var_f)

# length(tfs_hi_in_var)

# ra_98_lfc_int[1:10,1:10]

# log2(3)



# apply(ra_98_lfc_int[,tfs_hi_in_var], 2, quantile)

# tfs_up_per_clust <- apply(ra_98_lfc_int[,tfs_hi_in_var], 1, function(x) tfs_hi_in_var[x > log2(3)])
# tfs_down_per_clust <- apply(ra_98_lfc_int[,tfs_hi_in_var], 1, function(x) tfs_hi_in_var[x < -log2(3)])
# # tfs_up_per_clust
# # tfs_down_per_clust

# options(repr.plot.width = 25)
# options(repr.plot.height = 6)

# f_tfs_hi_in_var <- factor(tfs_hi_in_var, levels = tfs_hi_in_var[order(apply(ra_98_lfc_int[,tfs_hi_in_var], 2, median))])

# length(unique(f_tfs_hi_in_var))

# ?interaction

# sort(f_tfs_hi_in_var)

# sort(f_tfs_hi_in_var)[1:14]

# par(las = 2, mar = c(10,5,2,1))
# boxplot(ra_98_lfc_int[,sort(f_tfs_hi_in_var)[1:14]], lex.order = F)

# par(las = 2, mar = c(10,5,2,1))
# boxplot(ra_98_lfc_int[,f_tfs_hi_in_var[1:14]], lex.order = F)





# tfs_hi_in_var

# length(tfs_hi_in_var_f)
# length(unique(tolower(tfs_hi_in_var_f)))

# x_up <- sapply(tfs_up_per_clust, length)
# x_down <- sapply(tfs_down_per_clust, length)
# boxplot(c(x_up, x_down) ~ c(rep('up', length(x_up)), rep('down', length(x_down))))

# quantile(x_up)
# quantile(x_down)

# length(tfs_down_per_clust)

# # tfs_up_per_clust

# grep('ctcf', tfs_hi_in_var_f, ign = T, v=T)

# hc_var <- hclust(dist(ra_98_lfc_int[peak_clust_var,tfs_hi_in_var]), method = 'ward.D')
# hc_non_var <- hclust(dist(ra_98_lfc_int[peak_clust_non_var,tfs_hi_in_var]), method = 'ward.D')

# hc_mtfs <- hclust(dist(t(ra_98_lfc_int[peak_clust_var,tfs_hi_in_var])), method = 'ward.D2')
# ct_ppp4 <- cutree(hc_mtfs, k = 30)

# unique(ct_ppp4[hc_mtfs$order])

# unique(ct_ppp4[hc_mtfs$order])[1:15]

# # km_mtfs <- tglkmeans::TGL_kmeans(t(ra_98_lfc_int[peak_clust_var,tfs_hi_in_var]), k = 16)
# # ct_ppp4 <- setNames(km_mtfs$cluster, tfs_hi_in_var)

# accl <- as.data.frame(tibble::column_to_rownames(tibble::enframe(ct_ppp4, name = 'motif', value = 'cluster'), 'motif'))
# clclrs <- list(cluster = setNames(chameleon::distinct_colors(n = length(unique(ct_ppp4)))$name, 1:length(unique(ct_ppp4))))
# # accl
# # clclrs

# options(repr.plot.width = 25)
# options(repr.plot.height = 12)

# ppp4 <- pheatmap::pheatmap(ra_98_lfc_int[c(peak_clust_var[hc_var$order], peak_clust_non_var[hc_non_var$order]),tfs_hi_in_var[hc_mtfs$order[ct_ppp4[hc_mtfs$order] %in% unique(ct_ppp4[hc_mtfs$order])]]], cluster_rows = F, cluster_cols = F,
# # ppp4 <- pheatmap::pheatmap(ra_98_lfc_int[c(peak_clust_var[hc_var$order], peak_clust_non_var[hc_non_var$order]),tfs_hi_in_var[order(ct_ppp4)]], cluster_rows = F, cluster_cols = F,
#                            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), fontsize_col = 6,
#                           annotation_col = accl, annotation_color = clclrs)

# a_legc_avg_cl_ct <- t(tgs_matrix_tapply(a_legc_avg_cl, mcmd$cell_type, mean))

# hc_a_legc_avg_cl_ct <- hclust(dist(a_legc_avg_cl_ct - rowMeans(a_legc_avg_cl_ct)), method = 'ward.D2')

# unique(ct_ppp4[hc_mtfs$order])

# # ct_ppp4_samp <- names(sapply(unique(ct_ppp4[hc_mtfs$order])[1:16], function(u) {
# ct_ppp4_samp <- names(sapply(unique(ct_ppp4[hc_mtfs$order]), function(u) {
#     if (length(which(ct_ppp4 == u)) >1) {
#         sample(ct_ppp4[ct_ppp4 == u], 1)
#     } else {
#         ct_ppp4[ct_ppp4 == u]
#     }
# }))
# ct_ppp4_samp

# lapply(unique(ct_ppp4[hc_mtfs$order]), function(x) tfs_hi_in_var_f[ct_ppp4 == x])

# mtt <- c("JASPAR.EOMES", 
#          "JOLMA.MEIS2_mono_DBD_2",
#          'JASPAR.NEUROG1',
#          "JASPAR.NEUROD1",
#          'HOMER.Sox2',
#          "JASPAR.EMX1",
#          'JASPAR.POU3F2',
#          'JASPAR.NFIA',
#          "JASPAR.MEF2C",
#          "JASPAR.FOXP1",
#          "HOCOMOCO.MECP2_MOUSE.H11MO.0.C",
#          # 'HOCOMOCO.MBD2_MOUSE.H11MO.0.B',
#          'HOMER.CTCF',
#          # 'JASPAR.ZBTB14',
#          'HOMER.NRF1',
#          # 'JASPAR.KLF7',
#          # 'JASPAR.SP1',
#          'HOCOMOCO.KLF3_MOUSE.H11MO.0.A',
#          'JOLMA.ETV1_mono_DBD',
#          # 'JASPAR.ZBTB33',
#          # 'HOCOMOCO.SP3_MOUSE.H11MO.0.B',
#          # 'JASPAR.Zfx',
#          # 'HOMER.E2F1',
#          # 'HOCOMOCO.EGR1_MOUSE.H11MO.0.A',
#          # 'HOMER.E2F7',
#          'JASPAR.NFIB'
#          # 'JASPAR.RFX3',
#          # 'JASPAR.Atf1',
#          # 'HOCOMOCO.MYC_MOUSE.H11MO.0.A'
#         )

# lapply(unique(ct_ppp4[hc_mtfs$order]), function(x) tfs_hi_in_var[ct_ppp4 == x])

# m_var <- cbind(
#     ra_98_lfc_int[peak_clust_var,ct_ppp4_samp],
#                                         a_legc_avg_cl_ct[peak_clust_var,] - rowMeans(a_legc_avg_cl_ct[peak_clust_var,]))
# m_non_var <- cbind(
#     ra_98_lfc_int[peak_clust_non_var,ct_ppp4_samp],
#                                    a_legc_avg_cl_ct[peak_clust_non_var,] - rowMeans(a_legc_avg_cl_ct[peak_clust_non_var,]))
# hc_together_var <- hclust(dist(m_var),
#                         method = 'ward.D2')
# hc_together_non_var <- hclust(dist(m_non_var),
#                                    method = 'average')

# mtt <- c("JASPAR.EOMES", 
#          "JOLMA.MEIS2_mono_DBD_2",
#          'JASPAR.NEUROG1',
#          "JASPAR.NEUROD1",
#          'HOMER.Sox2',
#          "JASPAR.EMX1",
#          'JASPAR.POU3F2',
#          'JASPAR.NFIA',
#          "JASPAR.MEF2C",
#          "JASPAR.FOXP1",
#          "HOCOMOCO.MECP2_MOUSE.H11MO.0.C",
#          'HOMER.CTCF',
#          'HOMER.NRF1',
#          'HOCOMOCO.KLF3_MOUSE.H11MO.0.A',
#          'JOLMA.ETV1_mono_DBD',
#          'JASPAR.NFIB'
#         )
# ra_98_lfc_int[,mtt]

# # peak_clust_ord <- c(peak_clust_var[hc_together_var$order], peak_clust_non_var[hc_together_non_var$order])

# # peak_clust_ord <- c(peak_clust_var[order(apply(a_legc_avg_cl_ct[peak_clust_var,cust_st_ord2], 1, function(x) sum(x*1:length(x))/sum(x)))],
# peak_clust_ord <- c(peak_clust_var[order(apply(a_legc_avg_cl_ct[peak_clust_var,cust_st_ord2], 1, function(x) which.max(x)))],
#                                                peak_clust_non_var[hc_together_non_var$order])

# options(repr.plot.width = 8)
# options(repr.plot.height = 16)



# # ppp5 <- pheatmap::pheatmap(ra_98_lfc_int[c(peak_clust_var[hc_var$order], peak_clust_non_var[hc_non_var$order]),ct_ppp4_samp],
# # ppp5 <- pheatmap::pheatmap(ra_98_lfc_int[hc_a_legc_avg_cl_ct$order,ct_ppp4_samp],  
# # ppp5 <- pheatmap::pheatmap(ra_98_lfc_int[peak_clust_ord,ct_ppp4_samp[c(1:13, grep('CTCF|SP', tfs_hi_in_var_f[match(ct_ppp4_samp, tfs_hi_in_var)]))]], treeheight_col = 0,
# ppp5 <- pheatmap::pheatmap(ra_98_lfc_int[peak_clust_ord,mtt], treeheight_col = 0,
#                            cluster_rows = F, cluster_cols = T, 
#                            labels_col = tfs_hi_in_var_f[match(mtt, tfs_hi_in_var)],
#                            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), clustering_method = 'ward.D2', fontsize_col = 16,
#                           # annotation_col = accl, annotation_color = clclrs
#                           )

# # pltmt <- a_legc_avg_cl_ct[hc_a_legc_avg_cl_ct$order,cust_st_ord]
# # pltmt <- a_legc_avg_cl_ct[c(peak_clust_var[hc_var$order], peak_clust_non_var[hc_non_var$order]),cust_st_ord]
# pltmt <- a_legc_avg_cl_ct[peak_clust_ord,cust_st_ord2]
# ppp6 <- pheatmap::pheatmap(pltmt - rowMeans(pltmt), cluster_rows = F, cluster_cols = F,
#                            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-1,1,l=100), clustering_method = 'ward.D', fontsize_col = 16,
#                           # annotation_col = accl, annotation_color = clclrs
#                           )

# save_pheatmap(ppp5, './output/sequence_modeling/figs/motif_lfc_over_q_98_per_cluster_k_motif=30_1-16.png', h = 1800, w = 1000)

# save_pheatmap(ppp6, './output/sequence_modeling/figs/a_legc_avg_cl_ct_minus_rowMeans.png', h = 1800, w = 1000)

# amd <- prego::all_motif_datasets()

# options(repr.plot.width = 4)
# options(repr.plot.height = 4)

# u <- 2
# xu <- names(ct_ppp4[ct_ppp4 == u])
# lni <- length(xu)
# par(mfrow = c(sqrt(lni), ceiling(sqrt(lni))))
# lapply(xu, function(mtfi) prego::plot_pssm_logo(dplyr::filter(amd, motif == mtfi), title = mtfi))

# lapply(unique(ct_ppp4[hc_mtfs$order]), function(u) names(ct_ppp4[ct_ppp4 == u]))

# umi <- unique(amd$motif[match(grep('^SOX11', unique(amd$motif_orig), v=T, ign = T), amd$motif_orig)])
# umi <- umi[umi %in% colnames(ra_98_lfc_int)]

# umi

# plot(ra_98_lfc_int[peak_clust_var[hcct2],'HOMER.Sox2'],
#     ra_98_lfc_int[peak_clust_var[hcct2],'JASPAR.Sox11'])
# abline(0,1,col='red')

# p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hcct2],umi]),
# # p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together$order],motifs_to_take]),                            
#                             # clustering_method = 'ward.D2', 
#                             cluster_rows = F,cluster_cols = F, silent = F,fontsize_row = 14,
#                 color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
#                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51))                           )
# # save_pheatmap(p_plt, './output/sequence_modeling/figs/var_peak_internal_motif_enrichment.png', h = 400, w = 2400, res = 150)

# p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hcct2],umi]),
# # p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together$order],motifs_to_take]),                            
#                             # clustering_method = 'ward.D2', 
#                             cluster_rows = F,cluster_cols = F, silent = F,fontsize_row = 14,
#                 color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
#                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51))                           )
# # save_pheatmap(p_plt, './output/sequence_modeling/figs/var_peak_internal_motif_enrichment.png', h = 400, w = 2400, res = 150)

# p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hcct2],umi]),
# # p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together$order],motifs_to_take]),                            
#                             # clustering_method = 'ward.D2', 
#                             cluster_rows = F,cluster_cols = F, silent = F,fontsize_row = 14,
#                 color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
#                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51))                           )
# # save_pheatmap(p_plt, './output/sequence_modeling/figs/var_peak_internal_motif_enrichment.png', h = 400, w = 2400, res = 150)





# pltmt1 <- t(ra_98_lfc_int[peak_clust_var[ord1],motifs_to_take])
# tp1 <- dist(t(pltmt1))

# hc_col_pltmt <- hclust(tp1, method = 'complete')

# colMeans(ra_98_lfc_int[,motifs_to_take])

# # mtg <- dplyr::bind_cols(ra_98_lfc_int[peak_clust_var,motifs_to_take], a_legc_avg_cl[peak_clust_var,])
# # head(mtg)
# # hc_together <- hclust(dist(mtg), method = 'ward.D2')

# p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together_rg$order],motifs_to_take]),
# # p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together$order],motifs_to_take]),                            
#                             # clustering_method = 'ward.D2', 
#                             cluster_rows = F,cluster_cols = F, silent = F,fontsize_row = 14,
#                 color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
#                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51))                           )
# save_pheatmap(p_plt, './output/sequence_modeling/figs/var_peak_internal_motif_enrichment.png', h = 400, w = 2400, res = 150)

# p_plt <- pheatmap::pheatmap(rbind(t(ra_98_lfc_int[,motifs_to_take]), ifelse(1:length(km_a_legc$size) %in% peak_clust_var, 1, -1)),
# # p_plt <- pheatmap::pheatmap(t(ra_98_lfc_int[peak_clust_var[hc_together$order],motifs_to_take]),                            
#                             # clustering_method = 'ward.D2', 
#                             cluster_rows = F,cluster_cols = T, silent = F,fontsize_row = 14,
#                 color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
#                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51))                           )
# # save_pheatmap(p_plt, './output/sequence_modeling/figs/var_peak_internal_motif_enrichment.png', h = 400, w = 2400, res = 150)

# pltmt2 <- a_legc_avg_cl[peak_clust_var,rev(cust_mc_ord_st)]
# brks <- seq(-16.7,-13.5,l=100)
# clrmp <- colorRampPalette(c('white', 'lightpink', 'red', 'black'))(100)
# col_ha <- HeatmapAnnotation(cell_type = anno_simple(col_annot$cell_type[rev(cust_mc_ord_st)], 
#                                             col = ann_colors[['cell_type']], height =unit(1.5, 'cm')), 
#                             show_legend = F, which = 'row')
# row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_a_legc$size[peak_clust_var[hc_together_rg$order]]), 
#                                                         gp = gpar(fill = 'black',fontsize = 24), 
#                                                         axis_param = list(gp = gpar(fontsize = 10), labels_rot = 90)), 
#                             frac_prom = anno_barplot(as.numeric(t(perc_prom[peak_clust_var[hc_together_rg$order]])), ylim = c(0,1), 
#                                 gp = gpar(fill = 'black',fontsize = 24), 
#                                     axis_param = list(gp = gpar(fontsize = 40), labels_rot = 90)),
#                                     which = 'column' 
#                             # width = unit(12, 'cm')
#                            )
# ch <- ComplexHeatmap::Heatmap(matrix = t(pltmt2[hc_together_rg$order,]), col = clrmp, 
#                               show_row_names = F,
#                               left_annotation = col_ha, 
#                               top_annotation = row_ha, 
#                               show_heatmap_legend = T, show_column_names = T,
#                               heatmap_legend_param = list(legend_height = unit(2, 'in'), legend_width = unit(2, 'in'), labels_gp = gpar(fontsize = 12)),
#                               heatmap_width = unit(25, 'npc'), heatmap_height = unit(30, 'npc'),
#                         row_names_gp = gpar(fontsize = 24),
#                         cluster_columns = F, cluster_rows = F)

# draw(ch)

# png('./output/sequence_modeling/figs/a_legc_avg_cl_order_motif_enrichment_test2.png', h = 800, w = 1800, res= 150)
# draw(ch)
# dev.off()

# meth_tracks <- gtrack.ls('NSC_meth_CpG.avg')

# meth_cov_tracks <- gtrack.ls('NSC_meth_CpG.cov')

# intervs_all <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(a_legc))

# intervs_all$peak_name <- peak_names(intervs_all, tad_based = F)

# avg_meth_all <- gextract(meth_tracks, intervals = intervs_all, iterator = intervs_all)

# options(gmax.data.size = 1e+9)

# cov_meth_all <- gextract(meth_cov_tracks, intervals = intervs_all, iterator = 1)

# cov_meth_all_peaks <- tgs_matrix_tapply(t(cov_meth_all[,grep('meth', colnames(cov_meth_all))]), cov_meth_all$intervalID, sum, na.rm = T)

# avg_meth_all <- avg_meth_all[which(rowMins(cov_meth_all_peaks) >= 20),]

# avg_meth_all$peak_name <- intervs_all$peak_name[avg_meth_all$intervalID]
# rownames(avg_meth_all) <- avg_meth_all$peak_name
# # avg_meth_ct_peaks[,'ct'] <- unlist(purrr::map(stringr::str_split(ct_peaks_df$peak_name1, '--'), 1))
# y <- colnames(avg_meth_all)[grep('E\\d\\d', colnames(avg_meth_all))]
# colnames(avg_meth_all)[grep('E\\d\\d', colnames(avg_meth_all))] <- unlist(purrr::map(stringr::str_split(purrr::map(stringr::str_split(y, '\\.'), 2), '_'), 1))

# nrow(avg_meth_all)

# avg_meth_all$cluster <- km_a_legc$cluster[match(rownames(avg_meth_all), rownames(a_legc))]

# head(avg_meth)

# all(names(avg_meth) == rownames(avg_meth_all))

# cl_avg_meth <- tgs_matrix_tapply(t(avg_meth_all[,grep('E\\d\\d', colnames(avg_meth_all))]), avg_meth_all$cluster, mean, na.rm = T)

# head(cl_avg_meth)

# meth_bins <- apply(cl_avg_meth, 2, do.call('rbind', tapply(avg_meth_non_prom[pb], k_amn$cluster[match(pb, rownames(all_mat_norm_non_prom))], function(x) table(cut(x, breaks = c(-1e-5,seq(0.1,1,0.1))))))

# avg_meth_prom <- avg_meth[names(avg_meth) %in% prom_peaks]

# avg_meth_non_prom <- avg_meth[setdiff(names(avg_meth), prom_peaks)]

# options(repr.plot.height = 5, repr.plot.width = 20)

# meth_bins <- do.call('rbind', tapply(avg_meth, avg_meth_all$cluster, function(x) table(cut(x, breaks = c(-1e-5,0.15,0.5,0.75,1)))))
# meth_bins <- meth_bins/rowSums(meth_bins)

# pheatmap::pheatmap(t(meth_bins[peak_clust_var[hc_together_rg$order],ncol(meth_bins):1]), cluster_rows = F, cluster_cols = F, col = colorRampPalette(c('white', 'blue4'))(100))

# meth_temp_bins <- cl_avg_meth

# mba <- sapply(colnames(meth_bins), function(cni) {
#                setNames(colorRampPalette(c('white', 'purple'))(100), seq(0,1,l=100))
#            }, simplify = F)

# ca <- as.data.frame(meth_bins)
# # colnames(ca) <- paste0('fraction_NSC_meth=', colnames(ca))
# # ca$cl_size <- k_amn$size

# # ac <- list('mean_NSC_meth' = setNames(colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(100), seq(0,1,l=100)))
# ac <- list()

# ac[names(mba)] <- mba

# ## Randome genome intervals

# #' Generate random genome motif PSSM matrix
# #'

# gen_random_genome_peak_motif_matrix <- function(num_peaks = 1e+5,
#                                                 peak_width = 2e+2,
#                                                 bp_from_chrom_edge_to_avoid = 3e+6) {
#     .misha$ALLGENOME[[1]] <- .misha$ALLGENOME[[1]][!grepl("_", .misha$ALLGENOME[[1]]$chrom), ]
#     chrom_lens <- apply(.misha$ALLGENOME[[1]][, 2:3], 1, diff)
#     chrom_fracs <- setNames(chrom_lens / sum(chrom_lens), .misha$ALLGENOME[[1]][, 1])
#     sample_seqs <- mapply(chrom_fracs, names(chrom_fracs), chrom_lens, FUN = function(x, y, z) {
#         chrom <- rep(y, round(x * num_peaks))
#         start <- sample.int(n = z, size = length(chrom))
#         end <- start + peak_width
#         return(as.data.frame(rbind(chrom, start, end)))
#     })
#     sample_seqs <- as.data.frame(do.call("rbind", lapply(sample_seqs, t)))
#     sample_seqs[, 2:3] <- apply(sample_seqs[, 2:3], 2, as.numeric)
#     sample_seqs <- sample_seqs[with(sample_seqs, order(chrom, start, end)), ]
#     end_shift <- .misha$ALLGENOME[[1]][match(sample_seqs$chrom, .misha$ALLGENOME[[1]][, 1]), 3] - bp_from_chrom_edge_to_avoid
#     sample_seqs <- sample_seqs[sample_seqs$start >= bp_from_chrom_edge_to_avoid & sample_seqs$end <= end_shift, ]
#     sample_seqs <- PeakIntervals(sample_seqs)
#     return(sample_seqs)
# }
# rg_coords <- as.data.frame(gen_random_genome_peak_motif_matrix(num_peaks = nrow(seq_coords), peak_width = 3e+2))
# rg_coords$peak_name <- peak_names(rg_coords, tad_based = F)
# rg_energy <- readRDS('./output/sequence_modeling/random_genome_motif_energy.rds')
# rg_energy$peak_name <- peak_names(rg_energy, tad_based = F)
# rg_mat <- subset(rg_energy, select = -c(chrom, start, end, peak_name))
# rownames(rg_mat) <- rg_energy$peak_name
# rgvm <- apply(as.matrix(rg_mat), 2, function(x) {y <- x; y[x < -1e+6] <- min(y[y > -1e+6]); return(y)})
# rg_q98 <- plyr::aaply(rgvm, 2, quantile, .98, .parallel = T)
# ra_98_bin_rg <- t(plyr::laply(names(rg_q98), function(i) as.numeric(ie_mat[,i] >= rg_q98[[i]]), .parallel = T))
# colnames(ra_98_bin_rg) <- colnames(ie_mat)
# rownames(ra_98_bin_rg) <- rownames(ie_mat)
# save(ra_98_bin_rg, file='./output/sequence_modeling/feat_peak_all_motif_vs_random_genome_q_98_binary_matrix.rda')

# rg_energy <- readRDS('./output/sequence_modeling/random_genome_motif_energy.rds')
# rg_energy$peak_name <- peak_names(rg_energy, tad_based = F)
# rg_mat <- subset(rg_energy, select = -c(chrom, start, end, peak_name))
# rownames(rg_mat) <- rg_energy$peak_name
# rgvm <- apply(as.matrix(rg_mat), 2, function(x) {y <- x; y[x < -1e+6] <- min(y[y > -1e+6]); return(y)})

# dim(rg_mat)

# rg_ecdf <- apply(rg_mat[,unlist(tfs_hi_motifs)], 2, ecdf)
# # rg_ecdf <- apply(rg_mat[,unlist(tfs_hi_motifs)], 2, ecdf)

# ie_mat_rgq <- sapply(unlist(tfs_hi_motifs), function(i) rg_ecdf[[i]](ie_mat[,i]))
# rownames(ie_mat_rgq) <- rownames(ie_mat)
# colnames(ie_mat_rgq) <- unlist(tfs_hi_motifs)

# ie_mat_rgq_bin_98 <- tgs_matrix_tapply(t(ie_mat_rgq), km_a_legc$cluster, function(x) length(which(x >= 0.98)))

# ie_mat_rgq_lfc_98 <- apply(ie_mat_rgq_bin_98, 2, function(x) log2(1e-4 + x/(0.02*km_a_legc$size)))

# ie_mat_rgq_lfc_98_diff <- t(t(ie_mat_rgq_lfc_98) - colMeans(ie_mat_rgq_lfc_98))

# options(repr.plot.height = 12)
# options(repr.plot.width = 16)

# pheatmap::pheatmap(ie_mat_rgq_lfc_98[peak_clust_var[hcct2],colSds(ie_mat_rgq_lfc_98[peak_clust_var,]) >= 1], 
#                    col = colorRampPalette(c('blue3', 'white', 'red3'))(100),
#                    fontsize_col = 15, breaks = seq(-3,3,l=100), cluster_rows = F)

# pheatmap::pheatmap(ie_mat_rgq_lfc_98_diff[peak_clust_var[hcct2],colSds(ie_mat_rgq_lfc_98[peak_clust_var,]) >= 0.5], fontsize_col = 7,breaks = seq(-3,3,l=100), cluster_rows = F)





# ie_mat_rgq2 <- ie_mat_rgq

# ie_mat_rgq_avg_cl <- tgs_matrix_tapply(t(ie_mat_rgq), km_a_legc$cluster, mean)

# dim(ie_mat_rgq_avg_cl)

# options(repr.plot.height = 5, repr.plot.width = 18)

# hcct2

# pheatmap::pheatmap(t(ie_mat_rgq_avg_cl[peak_clust_var[hcct2],motifs_to_take]), cluster_cols = F, col = colorRampPalette(c('blue4', 'white', 'red4'))(100), breaks =seq(0.1,.9,l=100), cluster_rows = F)

# options(repr.plot.height = 10, repr.plot.width = 18)

# pheatmap::pheatmap(t(a_legc_avg_cl[peak_clust_var[hcct2],rev(cust_mc_ord_st)]), show_rownames = F, 
#                    annotation_legend = F, annotation_row = col_annot, annotation_colors = ann_colors, cluster_rows = F, cluster_cols = F, col = clrmp_abs, breaks =0.2*brks_abs-13)

# pheatmap::pheatmap(ie_mat_rgq_avg_cl)

# cmx <- colMaxs(ie_mat_rgq_avg_cl)
# cmn <- colMins(ie_mat_rgq_avg_cl)
# cm <- colMeans(ie_mat_rgq_avg_cl)
# mtr <- (cmx-cmn)
# p_imrac <- pheatmap::pheatmap(ie_mat_rgq_avg_cl[peak_clust_var,head(order(colSds(ie_mat_rgq_avg_cl)/colMeans(ie_mat_rgq_avg_cl), decreasing = T), 80)], fontsize_col = 14, col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(0,1,l=100), cluster_rows = T)

# matc <- ie_mat_rgq_avg_cl[peak_clust_var,head(order(colSds(ie_mat_rgq_avg_cl)/colMeans(ie_mat_rgq_avg_cl), decreasing = T), 80)]
# hc_matc <- hclust(dist(t(matc)), method = 'ward.D')
# hc_matc_ct <- cutree(hc_matc, 12)

# cmx <- colMaxs(ie_mat_rgq_avg_cl)
# cmn <- colMins(ie_mat_rgq_avg_cl)
# cm <- colMeans(ie_mat_rgq_avg_cl)
# mtr <- (cmx-cmn)
# p_imrac <- pheatmap::pheatmap(ie_mat_rgq_avg_cl[peak_clust_var,], fontsize_col = 14, col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(0,1,l=100), cluster_rows = T)

# colnames(ie_mat_rgq2) <- c('SOX', 'HOX', 'NEUROG2', 'EOMES', 'POU3F2')

# clr_lvls <- c('red3', 'green2', 'blue3')

# png('./output/sequence_modeling/figs/conditional_motif_energy_densities.png', h = 800, w = 1200)
# par(mfrow = c(5,5), cex.lab = 3, cex.main = 3)
# mari = rep(0.5,4)
# ttt <- sapply(colnames(ie_mat_rgq2), function(cni) {
#     x <- ie_mat_rgq2[,cni]
#     xc <- cut(x, breaks = c(0,0.5,0.9,1), labels = paste0(cni, '_', c('lo', 'mid', 'hi')))
#     cnjs <- colnames(ie_mat_rgq2)
#     tttt <- sapply(cnjs, function(cnj) {
#         den_motif <- tapply(ie_mat_rgq2[,cnj], xc, density)
#         if (cni == cnjs[[1]]) {
#             mari[[3]] <- 4
#             maini <- cnj
#         } else {
#             mari[[3]] <- 1
#             maini <- ''
#         }
#         if (cnj == cnjs[[1]]) {
#             mari[[2]] <- 6
#             yaxti <- 's'
#         } else {
#             mari[[2]] <- 0.5
#             yaxti <- 'n'
#         }
#         if (cni == cnjs[[length(cnjs)]]) {
#             xaxti <- 's'
#         } else {xaxti <- 'n'}
#         par(mar = mari)
#         plot(0,0,col = 'white', main = maini, ylab = cni, xlim = c(0,1), ylim = c(0, max(sapply(den_motif, function(x) max(x$y)))), xlab = '', xaxt = xaxti, yaxt = yaxti)
#         vvv <- sapply(seq_along(den_motif), function(j) {
#             points(den_motif[[j]]$x, den_motif[[j]]$y, main = paste0(c(cni, cnj), collapse = '--'), type = 'l', col = clr_lvls[[j]])
#         })
#     })
# })
# dev.off()

# # ie_mat_rgq <- sapply(1:ncol(ie_mat), function(i) rg_ecdf[[i]](ie_mat[,i]))
# # rownames(ie_mat_rgq) <- rownames(ie_mat)
# # colnames(ie_mat_rgq) <- colnames(ie_mat)
# # ie_mat_rgq_prom <- ie_mat_rgq[prom_peaks$peak_name,]
# # ie_mat_rgq_non_prom <- ie_mat_rgq[setdiff(rownames(ie_mat), prom_peaks$peak_name),]

# # dim(ie_mat_rgq)
# # dim(ie_mat_rgq_non_prom)
# # dim(ie_mat_rgq_prom)

# # options(repr.plot.width = 24)
# # options(repr.plot.height = 9)

# # par(mfcol = c(2,5))
# # vtvn <- sapply(motifs_to_take, function(x) {hist(ie_mat_rgq_non_prom[,x], 50, xlim = c(0,1), ylim = c(0,7e+3),main = paste('Dist. of random genome energy quantiles\n in distal elements\n', x));
# #                                            hist(ie_mat_rgq_prom[,x], 50, xlim = c(0,1), ylim = c(0,8.3e+3), main = paste('Dist. of random genome energy quantiles\n in promoters\n', x));})

# png('./output/sequence_modeling/figs/random_genome_quantiles_of_motifs_tss-dist+prox.png', h = 500, w = 1250)
# par(mfcol = c(2,5), cex.lab = 2, cex.axis = 1.5, cex.main = 2)
# vtvn <- sapply(motifs_to_take, function(x) {y1 <- density(ie_mat_rgq_non_prom[,x]); 
#                                             par(mar = c(2,3,4,1))
#                                             # if (x == motifs_to_take[[1]]) {par(mar = c(2,4,2,1))}
#                                             plot(y1$x, y1$y, xlab = '', xlim = c(-0.2,1.2), main = x, type = 'l', ylab = '', xaxt = 'n');
#                                             axis(1,at = seq(0,1,0.2))
#                                             # if (x == motifs_to_take[[1]]) {title(ylab = 'TSS-distal peaks', line = 3)}
#                                             y2 <- density(ie_mat_rgq_prom[,x]); 
#                                             par(mar = c(5,3,2,1))
#                                             # if (x == motifs_to_take[[1]]) {par(mar = c(2,4,2,1))}
#                                             plot(y2$x, y2$y, xlab = '', xlim = c(-0.2,1.2), type = 'l', ylab = '', xaxt = 'n')
#                                             axis(1,at = seq(0,1,0.2), labels = seq(0,1,0.2))
#                                             # if (x == motifs_to_take[[1]]) {title(ylab = 'TSS-proximal peaks', line = 3)}
#                                             # density(y1) hist(ie_mat_rgq_non_prom[,x], 50, xlim = c(0,1), ylim = c(0,7e+3),main = paste('Dist. of random genome energy quantiles\n in distal elements\n', x));
#                                            # hist(ie_mat_rgq_prom[,x], 50, xlim = c(0,1), ylim = c(0,8.3e+3), main = paste('Dist. of random genome energy quantiles\n in promoters\n', x));
#                                            })
# dev.off()

# legc_avg_cl <- t(tgs_matrix_tapply(legc, mcmd$cell_type, mean))

# cpn_genes <- head(rownames(legc)[order(legc_avg_cl[,'CPN_L2-3'] - rowMaxs(subset(legc_avg_cl, select = -c(`CPN_L2-3`))), decreasing = T)], 1000)

# cthpn_genes <- head(rownames(legc)[order(legc_avg_cl[,'CthPN'] - rowMaxs(subset(legc_avg_cl, select = -c(CthPN))), decreasing = T)], 500)

# prom_peaks$geneSymbol[match(prom_peaks$peak_name, rownames(ie_mat_rgq_prom)[order(ie_mat_rgq_prom[,'JOLMA.POU3F2_mono_DBD_1'])])]

# length(which(is.na(prom_peaks$geneSymbol[match(prom_peaks$peak_name, rownames(ie_mat_rgq_prom)[order(ie_mat_rgq_prom[,'JOLMA.POU3F2_mono_DBD_1'])])])))

# plot(ecdf(ie_mat_rgq_prom[prom_peaks$peak_name[c %in% cpn_genes],'JOLMA.POU3F2_mono_DBD_1']), do.points = F)
# plot(ecdf(ie_mat_rgq_prom[prom_peaks$peak_name[!(prom_peaks$geneSymbol %in% cpn_genes)],'JOLMA.POU3F2_mono_DBD_1']), add = T, col = 'red', do.points = F)

# hist(ie_mat_rgq_prom[prom_peaks$peak_name[prom_peaks$geneSymbol %in% cpn_genes],'JOLMA.POU3F2_mono_DBD_1'])

# hist(ie_mat_rgq_prom[prom_peaks$peak_name[!(prom_peaks$geneSymbol %in% cpn_genes)],'JOLMA.POU3F2_mono_DBD_1'])

# # par(mfrow = c(1,2))
# vtvn <- sapply(motifs_to_take, function(x) {plot(ecdf(ie_mat_rgq_non_prom[,x]), add = F, col  = 'black', do.points = F, main = x);
#                                            plot(ecdf(ie_mat_rgq_prom[,x]), add = T, col  = 'red', do.points = F);})

# load(file='./output/sequence_modeling/feat_peak_all_motif_vs_random_genome_q_98_binary_matrix.rda')

# ra_98_sum_clust_rg <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_rg)), km_a_legc$cluster, sum)

# ra_98_lfc_rg <- log2(1e-1 + ra_98_sum_clust_rg/as.numeric(0.02*km_a_legc$size))

# dim(ra_98_lfc_rg)

# cl_lfc_rg <- tgs_matrix_tapply()

# inds <- which(matrixStats::colMaxs(ra_98_lfc_rg) - matrixStats::colMins(ra_98_lfc_rg) >= 4)

# length(inds)

# hc_cor_98 <- hclust(dist(t(ra_98_lfc_rg[peak_clust_var,inds])), method = 'ward.D2')

# # ct_hc_lfc <- cutree(hc_cor_98, k = 10)

# # motifs_new <- unlist(lapply(sort(unique(ct_hc_lfc)), function(x) {mi <- names(ct_hc_lfc[ct_hc_lfc == x]); return(mi[which.max(colSds(ra_98_lfc_rg[peak_clust_var,mi]))])}))

# ct_hc_lfc <- cutree(hc_cor_98, k = 16)

# table(ct_hc_lfc)

# motifs_new <- unlist(lapply(sort(unique(ct_hc_lfc)), function(x) {
#     mi <- names(ct_hc_lfc[ct_hc_lfc == x]); 
#     if(length(mi) > 1) {return(mi[which.max(colSds(ra_98_lfc_rg[peak_clust_var,mi]))])}
#     else {return(mi)}
#     }))

# motifs_new

# write(x = unlist(motifs_new), file = './output/sequence_modeling/motif_set_for_xgb.txt')



# ra_98_bin_int <- t(plyr::laply(1:length(raq98), function(i) as.numeric(ie_mat[,i] >= raq98[[i]]), .parallel = T))

# colnames(ra_98_bin_int) <- colnames(ie_mat)
# rownames(ra_98_bin_int) <- rownames(ie_mat)

# ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int)), km_a_legc$cluster, sum)

# # hc_pltmt_me_rg <- hclust(dist(t(pltmt3)), method = 'ward.D2')

# # # mtg_rg <- dplyr::bind_cols(cl_lfc_rg[peak_clust_var,motifs_to_take], a_legc_avg_cl_ct[peak_clust_var,cust_st_ord])
# # mtg_rg <- apply(dplyr::bind_cols(ra_98_lfc_rg[peak_clust_var,motifs_to_take], a_legc_avg_cl_ct[peak_clust_var,cust_st_ord]), 2, function(x) x - mean(x))
# # head(mtg_rg)

# # hc_together_rg <- hclust(dist(mtg_rg), method = 'ward.D')

# # pheatmap::pheatmap(t(mtg_rg[hc_together_rg$order,]), cluster_rows = F, cluster_cols = F)

# # ca2 <- cbind(ra_98_lfc_rg[,motifs_to_take], ca)

# # for (mot in motifs_to_take) {
# #     ac[[mot]] <- setNames(colorRampPalette(c('blue4', 'white', 'red4'))(100),
# #                 seq(-3,3,l=100))
# #     }

# # options(repr.plot.width = 18)
# # options(repr.plot.height = 5)

# # p_plt_rg <- pheatmap::pheatmap(t(ra_98_lfc_rg[peak_clust_var[hc_together_rg$order],motifs_to_take]), clustering_method = 'ward.D2', cluster_cols = F, cluster_rows = F, fontsize = 14, annotation_col = ca, annotation_colors = ac, annotation_legend = F,
# # # p_plt_rg <- pheatmap::pheatmap(t(cl_lfc_rg[peak_clust_var,]), clustering_method = 'ward.D2', cluster_cols = T, cluster_rows = T, fontsize = 14,
# # # p_plt_rg <- pheatmap::pheatmap(pltmt3[,hc_pltmt_me_rg$order], cluster_cols = F, cluster_rows = F, fontsize = 14,
# #                   color = colorRampPalette(c('blue4', 'white', 'red4'))(100),
# #                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)))

# # # save_pheatmap(p_plt_rg, './output/sequence_modeling/figs/motif_enrichment_var_peak_clusters_vs_random_genome.png',
# # #                  h = 400, w = 2400, res = 150)

# # dim(a_legc_avg_cl_ct)

# # mat_all <- cbind(a_legc_avg_cl_ct, ca2)

# # maz <- apply(mat_all, 2, function(x) (x - mean(x))/sd(x))

# # maz_hc <- hclust(dist(maz), method = 'ward.D2')$order

# # ac[['cell_type']] <- ann_colors[['cell_type']]

# # options(repr.plot.width = 18)
# # options(repr.plot.height = 12)

# # p_plt_rg_2 <- pheatmap::pheatmap(t(a_legc_avg_cl[peak_clust_var[hc_together_rg$order],rev(cust_mc_ord_st)]), cluster_cols = F, cluster_rows = F, fontsize = 14, annotation_col = ca2, annotation_colors = ac, annotation_row = col_annot, annotation_legend = T, show_rownames = F,
# # # p_plt_rg <- pheatmap::pheatmap(t(cl_lfc_rg[peak_clust_var,]), clustering_method = 'ward.D2', cluster_cols = T, cluster_rows = T, fontsize = 14,
# # # p_plt_rg <- pheatmap::pheatmap(pltmt3[,hc_pltmt_me_rg$order], cluster_cols = F, cluster_rows = F, fontsize = 14,
# #                   color = colorRampPalette(c('white', 'red4', 'black'))(100),
# #                 # breaks = seq(-16.6,-13.5)
# #                                 )

# # # save_pheatmap(p_plt_rg, './output/sequence_modeling/figs/motif_enrichment_var_peak_clusters_vs_random_genome.png',
# # #                  h = 400, w = 2400, res = 150)

# # p_plt_rg <- pheatmap::pheatmap(t(ra_98_lfc_rg[peak_clust_var,unlist(motifs_new)]), clustering_method = 'ward.D2', cluster_cols = T, cluster_rows = T, fontsize = 14,
# # # p_plt_rg <- pheatmap::pheatmap(t(cl_lfc_rg[peak_clust_var,]), clustering_method = 'ward.D2', cluster_cols = T, cluster_rows = T, fontsize = 14,
# # # p_plt_rg <- pheatmap::pheatmap(pltmt3[,hc_pltmt_me_rg$order], cluster_cols = F, cluster_rows = F, fontsize = 14,
# #                   color = colorRampPalette(c('blue4', 'white', 'red4'))(100),
# #                 breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)))

# # # save_pheatmap(p_plt_rg, './output/sequence_modeling/figs/motif_enrichment_var_peak_clusters_vs_random_genome.png',
# # #                  h = 400, w = 2400, res = 150)

# # # peak_clusters_per_motif <- apply(ra_98_lfc_rg[peak_clust_var,motifs_to_take], 2, function(x) peak_clust_var[which(x >= 2 | x >= quantile(x, 0.8))])
# # peak_clusters_per_motif <- apply(ra_98_lfc_rg[peak_clust_var,unlist(motifs_new)], 2, function(x) peak_clust_var[which(x >= 2 | x >= quantile(x, 0.8))])
# # peak_clusters_per_motif

# # cl_tbx5 <- as.numeric(tail(names(sort(ra_98_lfc_rg[,'JASPAR.TBX5']))))

# # apply(a_legc_avg_cl[cl_tbx5,], 1, which.max)

# # mcmd$cell_type[apply(a_legc_avg_cl[cl_tbx5,], 1, which.max)]

# # qg <- quantile(a_legc_avg_cl[unique(unlist(peak_clusters_per_motif)),], 0.98)

# # mcs_per_motif <- purrr::imap(peak_clusters_per_motif, function(.x,.y) {
# #     union(unique(as.numeric(unlist(
# #         apply(a_legc_avg_cl[.x,], 1, function(x) which(x >= qg))))),
# #             unlist(apply(a_legc_avg_cl[as.numeric(tail(names(sort(ra_98_lfc_rg[,'JASPAR.TBX5'])))),], 1, function(x) tail(order(x))))
# #                          )
# #               # )))
# #               # )
# # })

# # save(mcs_per_motif, file='./output/sequence_modeling/mcs_per_motif_a_legc_avg_cl_quantile_98.rda')

# # dim(pltmt2)



# # pltmt2 <- a_legc_avg_cl[peak_clust_var,rev(cust_mc_ord_st)]
# # brks <- seq(-16.7,-13.5,l=100)
# # clrmp <- colorRampPalette(c('white', 'lightpink', 'red', 'black'))(100)
# # col_ha <- HeatmapAnnotation(cell_type = anno_simple(col_annot$cell_type[rev(cust_mc_ord_st)], 
# #                                             col = ann_colors[['cell_type']], height =unit(1.5, 'cm')), 
# #                             show_legend = F, which = 'row')
# # col_ha2 <- HeatmapAnnotation(mean_day = anno_simple(col_annot$mean_day[rev(cust_mc_ord_st)], 
# #                                                     col = circlize::colorRamp2(seq(13,18,1), c('red', 'orange', 'yellow', 'green2', 'blue', 'purple')), 
# #                                                     height =unit(1.5, 'cm')), show_legend = F, which = 'row')
# # row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_a_legc$size[peak_clust_var[hc_together_rg$order]]), gp = gpar(fill = 'black',fontsize = 24), 
# #                                                         axis_param = list(gp = gpar(fontsize = 10), labels_rot = 90)), 
# #                             frac_prom = anno_barplot(as.numeric(t(perc_prom[peak_clust_var])), ylim = c(0,1), 
# #                                 gp = gpar(fill = 'black',fontsize = 6), 
# #                                     axis_param = list(gp = gpar(fontsize = 10), labels_rot = 90)),
                            
                            
# #                                     which = 'column'
# #                             # width = unit(12, 'cm')
# #                            )

# # ch <- ComplexHeatmap::Heatmap(matrix = t(pltmt2[hc_together_rg$order,]), col = clrmp, 
# #                               show_row_names = F,
# #                               left_annotation = col_ha, 
# #                               right_annotation = col_ha2,
# #                               top_annotation = row_ha, 
# #                               show_heatmap_legend = F, show_column_names = T,
# #                               # heatmap_legend_param = list(legend_height = unit(2, 'in'), legend_width = unit(2, 'in'), labels_gp = gpar(fontsize = 12)),
# #                               heatmap_width = unit(25, 'npc'), heatmap_height = unit(30, 'npc'),
# #                         row_names_gp = gpar(fontsize = 24),
# #                         cluster_columns = F, cluster_rows = F)

# # draw(ch)

# # png('./output/sequence_modeling/figs/a_legc_avg_cl_order_motif_enrichment_test.png', h = 800, w = 1800, res= 150)
# # draw(ch)
# # dev.off()

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# smoothScatter(log2(1e-5 + egc_by_day_n[,1]), log2(1e-5 + egc_by_day_n[,ncol(egc_by_day_n)]))

# plot(a_legc_avg_ct[,'Astrocytes'], log2(1e-5 + egc_by_day_n[rownames(a_legc_avg_ct),ncol(egc_by_day_n)]))

# smoothScatter(a_legc_avg_ct[,'Astrocytes'], log2(1e-5 + egc_by_day_n[rownames(a_legc_avg_ct),1]))

# a_legc_by_day_n_ipc <- log(5e-6 + egc_by_day_n_ipc)





# quantile(mcmd$mean_day[match(ipc_mcs, mcmd$metacell)])

# ipc_cl_avg_day <- t(tgs_matrix_tapply(a_legc_avg_cl[,ipc_mcs], cut(mcmd$mean_day[match(ipc_mcs, mcmd$metacell)], breaks = 13:18), mean))

# colnames(ipc_cl_avg_day) <- paste0('IPC_', colnames(ipc_cl_avg_day))

# ipc_cl_and_other_ct <- cbind(ipc_cl_avg_day[enh_cl_ord,], a_legc_avg_cl_ct[enh_cl_ord,tail(cust_st_ord2, -2)])

# # ord_here <- c(head(enh_cl_ord, 13), tail(enh_cl_ord, -13)[order(rowMeans(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),1:5]))])
# ord_here <- c(head(enh_cl_ord, 13), tail(enh_cl_ord, -13)[hclust(dist(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),]), method = 'ward.D2')$order])
# ord_here

# # ord_here <- c(head(enh_cl_ord, 13), tail(enh_cl_ord, -13)[order(rowMeans(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),1:5]))])
# # ord_here <- c(head(enh_cl_ord, 13), tail(enh_cl_ord, -13)[order(abs(apply(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),c(1,5)], 1, diff)))])
# ord_here <- c(head(enh_cl_ord, 13), tail(enh_cl_ord, -13)[order(rowMeans(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),11:13] - rowMeans(ipc_cl_and_other_ct[tail(enh_cl_ord, -13),14:16])))])
# ord_here

# cah <- tibble::column_to_rownames(as.data.frame(cbind(colnames(ipc_cl_and_other_ct), c(rep('IPC', 5), tail(cust_st_ord2, -2)))), 'V1') %>% dplyr::rename(cell_type = V2)
# cah

# options(repr.plot.width = 6)
# options(repr.plot.height = 9)

# p_avg_cl_ct_ipc <- pheatmap::pheatmap(ipc_cl_and_other_ct[ord_here,], col = clrmp_abs, gaps_row = 13, gaps_col = c(5,ncol(ipc_cl_and_other_ct)-7, ncol(ipc_cl_and_other_ct) - 3),
#                                       cluster_rows = F,clustering_method = 'complete', cluster_cols = F, treeheight_row = 0, annotation_col = cah, annotation_colors = ann_colors, annotation_legend = F)
# # save_pheatmap(p_avg_cl_ct_ipc, filename = './output/mcatac/figs/phm_avg_cl_ct_ipc_by_md.png', h = 1000, w = 700)

# p_avg_cl_ct_ipc$gtable$grobs[[5]]$gp

# ipc_peaks_all <- ct_peaks$peak_name[ct_peaks$type == 'ipc_peak']

# length(ipc_peaks_all)

# jacc2(ipc_peaks_all)

# jacc2(ipc_peaks_all, ipc_inc_peaks)
# jacc2(ipc_peaks_all, ipc_dec_peaks)



# jacc(ipc_inc_peaks, ipc_dec_peaks)

# length(ipc_inc_peaks)

# options(repr.plot.width = 14)
# options(repr.plot.height = 10)
# mi1 <- a_legc[ipc_inc_peaks,]
# mi1_hc <- hclust(tgs_dist(mi1), method = 'ward.D2')
# mi2 <- a_legc[setdiff(ipc_peaks_all, union(ipc_inc_peaks, ipc_dec_peaks)),]
# mi2_hc <- hclust(tgs_dist(mi2), method = 'ward.D2')
# mi3 <- a_legc[ipc_dec_peaks,]
# mi3_hc <- hclust(tgs_dist(mi3), method = 'ward.D2')


# mi1_avg_ct <- t(tgs_matrix_tapply(mi1, mcmd$cell_type, mean))
# mi2_avg_ct <- t(tgs_matrix_tapply(mi2, mcmd$cell_type, mean))
# mi3_avg_ct <- t(tgs_matrix_tapply(mi3, mcmd$cell_type, mean))

# col_key

# options(repr.plot.width = 7)
# options(repr.plot.height = 14)
# png('./output/mcatac/figs/temporal_ipc_peaks_by_ct_boxplots.png', h = 700, w = 400)
# par(las = 2, mar = c(1,5,3,1), mfrow = c(3,1), cex.axis = 2, cex.main = 2)
# boxplot(mi1_avg_ct[,cust_st_ord2], ylim = c(-16.6, -14.5), horizontal = F, col = col_key[cust_st_ord2], xaxt = 'n', main = 'Temporally activating peaks')
# par( mar = c(1,5,2,1))
# boxplot(mi2_avg_ct[,cust_st_ord2], ylim = c(-16.6, -14.5), horizontal = F, col = col_key[cust_st_ord2], xaxt = 'n', main = 'Generally accessible peaks')
# par( mar = c(10,5,2,1))
# boxplot(mi3_avg_ct[,cust_st_ord2], ylim = c(-16.6, -14.5), horizontal = F, col = col_key[cust_st_ord2], main = 'Temporally deactivating peaks')
# dev.off()



# day_mat <- as.matrix(mcmd[,grep('^E\\d\\d', colnames(mcmd))])
# rownames(day_mat) <- mcmd$metacell
# day_mat_norm <- day_mat/rowSums(day_mat)

# col_annot[,colnames(day_mat_norm)] <- day_mat_norm

# for (cni in colnames(day_mat_norm)) {
#     ann_colors[[cni]] <- setNames(colorRampPalette(c('white', 'black'))(2), c(0,1))
# }

# pltmt <- rbind(mi1[mi1_hc$order,cust_mc_ord_st2], 
#                          mi2[mi2_hc$order,cust_mc_ord_st2],
#                          mi3[mi3_hc$order,cust_mc_ord_st2])
# top_ha <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2,'cell_type'], 
#                                                      col = ann_colors[['cell_type']],
#                                                      # height =unit(2, 'cm')
#                                                     ), 
#                              mean_day = anno_lines(x = col_annot[cust_mc_ord_st2,'mean_day'], axis_param = list(gp = gpar(fontsize = 0)),
#                                                     # col = circlize::colorRamp2(seq(13,18,1), c('red', 'orange', 'yellow', 'green2', 'blue', 'purple')), 
#                                                     # height =unit(2, 'cm')
#                                                   ), 
#                              annotation_name_gp = gpar(fontsize = 18),
#                              show_legend = F)
# # top_ha <- HeatmapAnnotation(cell_type = anno_simple(x = col_annot[match(cust_mc_ord_st2, mcmd$metacell), 'cell_type'],
# #                             col = ann_colors[['cell_type']],
# #                                                      # height =unit(2, 'cm'),
# #                                                     )
# #                            mean_day = anno_lines(x = col_annot[match(cust_mc_ord_st2, mcmd$metacell), 'mean_day'] ),
# #                            which = 'column')

# ch_ipc_peaks <- ComplexHeatmap::Heatmap(pltmt, top_annotation = top_ha, show_column_names = F,show_row_names = F,
#                                         col = circlize::colorRamp2(breaks = seq(min(pltmt), max(pltmt), l = 5), colors = c('white', 'orange', 'red', 'purple', 'black')),
#                                 column_split = factor(names(cust_mc_ord_st2), levels = cust_st_ord2), 
#                                         # column_gap = unit(2, 'mm'),
#                               column_title_gp = gpar(fontsize = 20),
#                               column_title_rot = 90,
#                               row_split = factor(c(rep(1, nrow(mi1)), rep(2, nrow(mi2)), rep(3, nrow(mi3))), levels = c(1,2,3),
#                                 labels = c('Temporally activating', 'Generally accessible', 'Temporally deactivating')), 
#                                         # row_gap = unit(2, 'mm'),
#                                                             row_title_gp = gpar(fontsize = 20),
#                               show_heatmap_legend = T,
#                               # heatmap_legend_param = list(legend_height = unit(5, 'in'), legend_width = unit(5, 'in'), labels_gp = gpar(fontsize = 16)),
#                               heatmap_width = unit(35, 'cm'), heatmap_height = unit(35, 'cm'),
#                               # row_names_gp = gpar(fontsize = 14),                                                                                              
#                    # annotation_col = subset(col_annot, select = -mean_day), annotation_colors = ann_colors, 
#                                         cluster_columns = F,
#                                         # clustering_method = 'ward.D2', 
#                                         cluster_rows = F
#                                         # , col = clrmp_abs
#                                        )

# png('./output/mcatac/figs/ipc_atac_peaks_groups_legc_phm.png', h = 1600, w = 1400)
# draw(ch_ipc_peaks)
# dev.off()







# pltmt_avg_ct <- t(tgs_matrix_tapply(pltmt, names(cust_mc_ord_st2), mean))

# options(repr.plot.width = 9)
# options(repr.plot.height = 15)

# ppp <- ComplexHeatmap::Heatmap(pltmt_avg_ct[,cust_st_ord2], row_gap = unit(4, 'mm'),
#                                name = 'Cell type\nmean\nATAC', cluster_columns = F, cluster_rows = F, show_row_names = F, col = clrmp_abs, column_names_gp = gpar(fontsize = 20),
#                                                             row_split = factor(c(rep(1, nrow(mi1)), rep(2, nrow(mi2)), rep(3, nrow(mi3))), levels = c(1,2,3), 
#                                                                                labels = c('Temporally activating', 'Generally accessible', 'Temporally deactivating')),
#                               row_title_gp = gpar(fontsize = 20))

# draw(ppp)

# png('./output/mcatac/figs/ipc_atac_peaks_groups_legc_avg_ct_phm.png', h = 1000, w = 450)
# draw(ppp)
# dev.off()

# options(repr.plot.width = 15)
# options(repr.plot.height = 15)

# pheatmap::pheatmap(rbind(mi1[mi1_hc$order,as.character(ro[ro %in% as.numeric(mcmd$metacell)])], 
#                          mi2[mi2_hc$order,as.character(ro[ro %in% as.numeric(mcmd$metacell)])],
#                          mi3[mi3_hc$order,as.character(ro[ro %in% as.numeric(mcmd$metacell)])]), show_colnames = F, 
#                    annotation_col = col_annot, annotation_colors = ann_colors, show_rownames = F, cluster_cols = F,clustering_method = 'ward.D2', cluster_rows = F, col = clrmp_abs)

# quantile(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)]), (0:20)/20)

# length(which(abs(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)])) >= 1))

# ipc_inc_peaks <- names(which(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)]) >= 1))
# ipc_dec_peaks <- names(which(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)]) <= -1))

# options(repr.plot.width = 8)
# options(repr.plot.height = 8)
# quantile(rm_ama[intersect(names(rm_ama), ipc_inc_peaks)], (0:20)/20)
# quantile(rm_ama[intersect(names(rm_ama), ipc_dec_peaks)], (0:20)/20)
# hist(rm_ama[intersect(names(rm_ama), ipc_inc_peaks)], 20)
# hist(rm_ama[intersect(names(rm_ama), ipc_dec_peaks)], 20)

# length(ipc_inc_peaks)

# length(ipc_dec_peaks)

# length(intersect(ipc_inc_peaks, prom_peaks$peak_name))/length(ipc_inc_peaks)

# length(intersect(ipc_dec_peaks, prom_peaks$peak_name))/length(ipc_dec_peaks)

# length(intersect(ipc_inc_peaks, dist_peaks$peak_name))

# dim(a_legc[ipc_inc_peaks,])

# ipc_inc_hc <- hclust(tgs_dist(a_legc[intersect(ipc_inc_peaks, rownames(a_legc)),]), method = 'ward.D2')

# ipc_dec_hc <- hclust(tgs_dist(a_legc[intersect(ipc_dec_peaks, rownames(a_legc)),]), method = 'ward.D2')

# length(ipc_inc_peaks)
# length(ipc_dec_peaks)

# ipc_diff_peaks <- intersect(intersect(c(ipc_inc_peaks[ipc_inc_hc$order], ipc_dec_peaks[ipc_dec_hc$order]), which(dist_peaks$peak_name)

# dim(a_legc_by_day_n_ipc)

# # plot(rowMeans(a_legc_by_day_n_ipc[ipc_diff_peaks,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[ipc_diff_peaks,c(1:2)]), a_legc_avg_ct[ipc_diff_peaks,'CthPN'] - a_legc_avg_ct[ipc_diff_peaks,'CPN_L2-3'])
# plot(rowMeans(a_legc_by_day_n_ipc[intersect(ipc_inc_peaks, rownames(a_legc_avg_ct)),c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[intersect(ipc_inc_peaks, rownames(a_legc_avg_ct)),c(1:2)]),
#      a_legc_avg_ct[intersect(ipc_inc_peaks, rownames(a_legc_avg_ct)),'CthPN'] - a_legc_avg_ct[intersect(ipc_inc_peaks, rownames(a_legc_avg_ct)),'CPN_L2-3'])

# length(which(rowMeans(a_legc_by_day_n_ipc[peaks_here,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[peaks_here,c(1:2)]) >= 0))

# length(rownames(a_legc_by_day_n_ipc))

# peaks_here <- intersect(rownames(a_legc_by_day_n_ipc)[which(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)]) >= 0.5)], rownames(a_legc_avg_ct))

# xvec <- rowMeans(a_legc_by_day_n_ipc[peaks_here,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[peaks_here,c(1:2)])
# yvec <- a_legc_avg_ct[peaks_here,'CPN_L2-3'] - a_legc_avg_ct[peaks_here,'CthPN']
# boxplot_vec(xvec = xvec, yvec =yvec,
#             num_bins = 9,ylab = 'CPN_L2-3 - CthPN  ATAC', xlab = 'E18 - E13 IPC ATAC',
#             nm = 'CPN_L2-3 - CthPN delta ATAC vs ipc temporal trend')

# par(las = 2, mar = c(7, 5, 3,1))
# boxplot_vec(xvec = rowMeans(a_legc_by_day_n_ipc[intersect(rownames(a_legc_by_day_n_ipc),rownames(a_legc_avg_ct)),c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[intersect(rownames(a_legc_by_day_n_ipc),rownames(a_legc_avg_ct)),c(1:2)]), 
#             yvec = a_legc_avg_ct[,'CPN_L2-3'] - a_legc_avg_ct[,'CthPN'], 
#             num_bins = 13,ylab = 'CPN_L2-3 - CthPN  ATAC', xlab = 'E18 - E13 IPC ATAC',
#             nm = 'CPN_L2-3 - CthPN delta ATAC vs ipc temporal trend')

# head(which(abs(rowMeans(a_legc_by_day_n_ipc[,c(5:6)]) - rowMeans(a_legc_by_day_n_ipc[,c(1:2)])) >= 1))

# pheatmap::pheatmap(a_legc_avg_ct[intersect(c(ipc_inc_peaks[ipc_inc_hc$order], ipc_dec_peaks[ipc_dec_hc$order]), dist_peaks$peak_name),cust_st_ord2], show_rownames = F, cluster_cols = F, cluster_rows = F, col = clrmp_abs)

# ipc_diff_peaks <-intersect(c(ipc_inc_peaks, ipc_dec_peaks), dist_peaks$peak_name)
# hc_ipc_diff_peaks <- hclust(tgs_dist(a_legc[ipc_diff_peaks,]), method = 'ward.D2')

# options(repr.plot.width = 14)
# options(repr.plot.height = 10)
# pheatmap::pheatmap(a_legc[c(ipc_inc_peaks[ipc_inc_hc$order], ipc_dec_peaks[ipc_dec_hc$order]),cust_mc_ord_st2], show_colnames = F, 
#                    annotation_col = col_annot, annotation_colors = ann_colors, show_rownames = F, cluster_cols = F,clustering_method = 'ward.D2', cluster_rows = F, col = clrmp_abs)

# cor(colMeans(a_legc[intersect(pks_diff_cfupn_clpl23, ipc_inc_peaks),]),
#    colMeans(a_legc[intersect(pks_diff_clpl23_cfupn, ipc_inc_peaks),]), method = 'spearman')

# options(repr.plot.width = 12)
# options(repr.plot.height = 12)
# plot(colMeans(a_legc[intersect(pks_diff_cfupn_clpl23, ipc_inc_peaks),]),
#    colMeans(a_legc[intersect(pks_diff_clpl23_cfupn, ipc_inc_peaks),]), col = mcmd$color, pch = 16, cex = 2)

# options(repr.plot.width = 16)
# options(repr.plot.height = 8)
# par(mfrow = c(1,2))
# plot(colMeans(a_legc[pks_diff_cfupn_clpl23,]),
#    colMeans(a_legc[pks_diff_clpl23_cfupn,]), col = mcmd$color, pch = 16, cex = 2, xlab = 'Mean CfuPN peak accessibility', ylab = 'Mean CPN peak accessibility')
# md_v <- mcmd$mean_day
# md_clvls <- rainbow(100)[1+round(99*(md_v - min(md_v))/6)]
# plot(colMeans(a_legc[pks_diff_cfupn_clpl23,]),
#    colMeans(a_legc[pks_diff_clpl23_cfupn,]), col = md_clvls, pch = 16, cex = 2, xlab = 'Mean CfuPN peak accessibility', ylab = 'Mean CPN peak accessibility')

# options(repr.plot.width = 16)
# options(repr.plot.height = 8)
# par(mfrow = c(1,2))
# plot(md_v, colMeans(a_legc[pks_diff_cfupn_clpl23,]), col = mcmd$color, pch = 16, cex = 2, xlab = 'Mean day', ylab = 'Mean CfuPN peak accessibility')
# plot(md_v, colMeans(a_legc[pks_diff_clpl23_cfupn,]), col = mcmd$color, pch = 16, cex = 2, xlab = 'Mean day', ylab = 'Mean CPN peak accessibility')

# length(intersect(pks_diff_cfupn_clpl23, ipc_inc_peaks))
# length(intersect(pks_diff_cfupn_clpl23, ipc_dec_peaks))
# length(intersect(pks_diff_clpl23_cfupn, ipc_inc_peaks))
# length(intersect(pks_diff_clpl23_cfupn, ipc_dec_peaks))


# jacc2 <- function (x, y) 
# {
#     return(length(intersect(x, y))/min(c(length(x), length(y))))
# }

# jacc2(pks_diff_cfupn_clpl23, ipc_inc_peaks)
# jacc2(pks_diff_cfupn_clpl23, ipc_dec_peaks)
# jacc2(pks_diff_clpl23_cfupn, ipc_inc_peaks)
# jacc2(pks_diff_clpl23_cfupn, ipc_dec_peaks)


# load('./output/mcatac/ct_peaks.rda')

# intersect(ct_peaks$peak_name[ct_peaks$type == 'astro_peak'], prom_peaks$peak_name)

# sort(table(ct_peaks$type[ct_peaks$peak_name %in% ipc_inc_peaks]))

# sort(table(ct_peaks$type[ct_peaks$peak_name %in% ipc_dec_peaks]))

# length(pks_diff_cfupn_clpl23)

# options(repr.plot.width = 14)
# options(repr.plot.height = 10)
# pheatmap::pheatmap(a_legc[pks_diff_cfupn_clpl23,cust_mc_ord_st2], show_colnames = F, 
#                    annotation_col = col_annot, annotation_colors = ann_colors, show_rownames = F, cluster_cols = F, cluster_rows = T, clustering_method = 'ward.D2', col = clrmp_abs)

# quantile(rowMeans(egc_by_day_n_ipc[,c(5:6)]) - rowMeans(egc_by_day_n_ipc[,c(1:2)]), (0:20)/20)

# quantile(egc_by_day_n_ipc[,c(6)], (0:20)/20)
# quantile(egc_by_day_n_ipc[,c(1)], (0:20)/20)

# legc_by_day_n <- log2(1e-5 + egc_by_day_n)

# peaks_ord_nsc_t_diff <- sort(rowMeans(legc_by_day_n[,c('E17','E18')]) - rowMeans(legc_by_day_n[,c('E14','E13')]))

# delta_astro_nsc <- a_legc_avg_ct[,'Astrocytes'] - a_legc_avg_ct[,'NSC']

# delta_ipc_nsc <- a_legc_avg_ct[,'IPC'] - a_legc_avg_ct[,'NSC']

# nsc_pks_inc <- names(which(peaks_ord_nsc_t_diff >= 0.5))

# nsc_pks_dec <- names(which(peaks_ord_nsc_t_diff <= -0.5))

# astro_peaks <- names(which(delta_astro_nsc >= 0.5))

# source('./scripts/util.r')

# length(which(delta_ipc_nsc >= 1))
# length(which(delta_ipc_nsc <= -1))
# length(which(delta_ipc_nsc >= log2(1.5)))
# length(which(delta_ipc_nsc <= -log2(1.5)))

# 2**(1/2)

# log2(1/2)

# jacc(nsc_pks_inc, nsc_pks_dec)
# jacc(nsc_pks_inc, astro_peaks)
# jacc(nsc_pks_inc, names(nsc_ipc_pks_ord))
# jacc(nsc_pks_dec, names(nsc_ipc_pks_ord))

# length(intersect(nsc_pks_inc, union(nsc_pks_dec, union(astro_peaks, names(nsc_ipc_pks_ord)))))

# length(intersect(nsc_pks_inc, astro_peaks))
# min(length(nsc_pks_inc), length(astro_peaks))

# length(nsc_pks_inc)
# length(nsc_pks_dec)
# length(astro_peaks)
# length(nsc_ipc_pks_ord)
# length(nsc_pks_inc) +
# length(nsc_pks_dec) +
# length(astro_peaks) +
# length(nsc_ipc_pks_ord)

# km_a_nsc_pks_inc <- tglkmeans::TGL_kmeans(a_legc[nsc_pks_inc,mcs_here], k = 5)

# km_a_nsc_pks_dec <- tglkmeans::TGL_kmeans(a_legc[nsc_pks_dec,mcs_here], k = 5)

# km_a_nsc_pks_astro <- tglkmeans::TGL_kmeans(a_legc[astro_peaks,mcs_here], k = 3)

# options(repr.plot.width = 14)
# options(repr.plot.height = 10)

# mcs_here <- ro[ro %in% which(mcmd$cell_type %in% c('NSC', 'IPC', 'IPC_cyc', 'Astrocytes', 'OPCs'))]
# pks_here <- c(nsc_pks_inc[order(km_a_nsc_pks_inc$cluster)],
#                   nsc_pks_dec[order(km_a_nsc_pks_dec$cluster)],
#              astro_peaks[order(km_a_nsc_pks_astro$cluster)],
#              names(nsc_ipc_pks_ord))
# pltmt <- a_legc[pks_here,mcs_here]
# pltmt <- pltmt - rowMeans(pltmt)

# ppp4<- pheatmap::pheatmap(pltmt, show_rownames = F,show_colnames = F, annotation_legend = F,
#                    # clustering_method = 'ward.D',
#                    cluster_cols = F, 
#                    cluster_rows = F, annotation_col = col_annot, annotation_colors = ann_colors, 
#                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),
#                   breaks = seq(-2,2,l=100))

# # mcs_here <- ro[ro %in% which(mcmd$cell_type %in% c('NSC', 'IPC', 'IPC_cyc', 'Astrocytes', 'OPCs'))]
# mcs_here <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'IPC', 'IPC_cyc', 'Astrocytes', 'OPCs')]
# pks_here <- c(nsc_pks_inc[order(km_a_nsc_pks_inc$cluster)],
#                   nsc_pks_dec[order(km_a_nsc_pks_dec$cluster)],
#              astro_peaks[order(km_a_nsc_pks_astro$cluster)],
#              names(nsc_ipc_pks_ord))
# pltmt <- a_legc[pks_here,mcs_here]
# pltmt <- pltmt - rowMeans(pltmt)

# ppp4<- pheatmap::pheatmap(pltmt, show_rownames = F,show_colnames = F, annotation_legend = F,
#                    # clustering_method = 'ward.D',
#                    cluster_cols = F, 
#                    cluster_rows = F, annotation_col = col_annot, annotation_colors = ann_colors, 
#                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),
#                   breaks = seq(-2,2,l=100))



# # diff_ama <- apply(avg_meth_all[,c('E13','E17')], 1, diff)
# diff_ama <- rowMeans(avg_meth_all[,c('E16','E17')]) - rowMeans(avg_meth_all[,c('E13','E14')])

# head(diff_ama)

# boxplot(diff_ama[pks_here] ~ c(rep('inc', length(nsc_pks_inc)),
#                             rep('dec', length(nsc_pks_dec)),
#                             rep('astro', length(astro_peaks)),
#                             rep('ipc', length(nsc_ipc_pks_ord))))







# delta_egc_ipc_nsc <- log2(1e-5 + egc_by_day_n_ipc[,tail(colnames(egc_by_day_n_ipc), -2)]) - log2(1e-5 + egc_by_day_n[,head(colnames(egc_by_day_n), -2)])
# cor_delta_atac_meth <- tgs_cor(delta_egc_ipc_nsc[rownames(avg_meth_all),], t(apply(as.matrix(avg_meth_all[,grep('E\\d\\d', colnames(avg_meth_all))]), 1, diff)), spearman = T)
# meth_lvls <- apply(avg_meth_all[,grep('E\\d\\d', colnames(avg_meth_all))], 2, cut, breaks = seq(-1e-5,1,0.2))
# rownames(meth_lvls) <- rownames(avg_meth_all)

# delta_astro_nsc <- a_legc_avg_ct[,'Astrocytes'] - a_legc_avg_ct[,'NSC']

# quantile(delta_astro_nsc)

# quantile(delta_ipc_nsc)

# length(which(abs(delta_astro_nsc) >= 0.5))

# length(which(abs(delta_ipc_nsc) >= 0.5))

# quantile(delta_ipc_nsc)

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# smoothScatter(a_legc_avg_ct[which(km_a_legc$cluster %in% peak_clust_var),'Astrocytes'], a_legc_avg_ct[which(km_a_legc$cluster %in% peak_clust_var),'IPC'])

# plot(a_legc_avg_ct[which(km_a_legc$cluster %in% peak_clust_var),'Astrocytes'], a_legc_avg_ct[which(km_a_legc$cluster %in% peak_clust_var),'IPC'], pch = 16, cex = 0.1)

# a_legc_avg_ct <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, mean, na.rm = T))
# head(a_legc_avg_ct)

# xnsc <- a_legc_avg_ct[,'NSC']

# nsc_lvls <- cut(xnsc, breaks = quantile(xnsc, seq(0,1,1/3)), labels = paste0('NSC_', c('lo', 'mid', 'hi')))

# xipc <- a_legc_avg_ct[,'IPC']
# delta_ipc_nsc <- xipc - xnsc

# ipc_vs_nsc <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = ipc_mcs, mc_neg = nsc_mcs)

# astrocytes_vs_nsc <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = nsc_mcs, mc_neg = which(mcmd$cell_type == 'Astrocytes'))

# cpn_vs_ipc <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = which(mcmd$cell_type == 'CPN_L2-3'), mc_neg = ipc_mcs)

# cthpn_vs_ipc <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = which(mcmd$cell_type == 'CthPN'), mc_neg = ipc_mcs)

# cthpn_vs_cpn <- get_genes_specific_to_mcs(legc = a_legc, mc_pos = which(mcmd$cell_type == 'CthPN'), mc_neg = which(mcmd$cell_type == 'CPN_L2-3'))

# motifs_diff_between_peak_sets <- function(pks1, pks2, ie_mat) {
#     v1 <- colMeans(ie_mat[pks1,])
#     v2 <- colMeans(ie_mat[pks2,])
#     mtfs1 <- head(sort((v1 - v2)/(v1 + v2), decreasing = T), 100)
#     mtfs2 <- tail(sort((v1 - v2)/(v1 + v2), decreasing = T), 100)
#     return(list(pos_motifs = mtfs1, neg_motifs = mtfs2))
# }

# ipc_vs_nsc_motifs <- motifs_diff_between_peak_sets(names(ipc_vs_nsc)[ipc_vs_nsc >= 1], names(ipc_vs_nsc)[ipc_vs_nsc <= -1], ie_mat)

# astrocytes_vs_nsc_motifs <- motifs_diff_between_peak_sets(names(astrocytes_vs_nsc)[astrocytes_vs_nsc >= 1], names(astrocytes_vs_nsc)[astrocytes_vs_nsc <= -1], ie_mat)

# cpn_vs_ipc_motifs <- motifs_diff_between_peak_sets(names(cpn_vs_ipc)[cpn_vs_ipc >= 1], names(cpn_vs_ipc)[cpn_vs_ipc <= -1], ie_mat)

# cthpn_vs_ipc_motifs <- motifs_diff_between_peak_sets(names(cthpn_vs_ipc)[cthpn_vs_ipc >= 1], names(cthpn_vs_ipc)[cthpn_vs_ipc <= -1], ie_mat)

# cthpn_vs_cpn_motifs <- motifs_diff_between_peak_sets(names(cthpn_vs_cpn)[cthpn_vs_cpn >= 1], names(cthpn_vs_cpn)[cthpn_vs_cpn <= -1], ie_mat)

# motifs_to_take <- unique(unlist(lapply(list(ipc_vs_nsc_motifs, cpn_vs_ipc_motifs, astrocytes_vs_nsc_motifs, cthpn_vs_ipc_motifs, cthpn_vs_cpn_motifs), function(x) c(head(names(x[['pos_motifs']]), 15), tail(names(x[['neg_motifs']]), 15)))))



# din_lvls <- cut(delta_ipc_nsc, quantile(delta_ipc_nsc, c(0,0.1,0.5,0.9,1)), labels = paste0('delta{', c('--','-','+', '++'),'}'))

# options(repr.plot.width = 14)
# options(repr.plot.height = 12)

# # ttt <- lapply(motifs_to_take, function(mti) {
# png('./output/sequence_modeling/figs/neurog2_rg_quantiles_delta_ipc_nsc.png', h = 600, w = 1200)
# mti <- motifs_to_take[[3]]
#     databp = as.data.frame(list(hox = ie_mat_rgq[pb,mti], nsc_lvl = nsc_lvls, delta_lvl = din_lvls))

# par(las = 2, mar = c(18,9,4,1), cex.axis = 2, cex.lab = 2.5, cex.main = 2)
# boxplot(hox ~ ., data = databp, xlab = '', ylab = '', main = 'Neurog2 random genome energy quantile stratified by NSC ATAC and delta IPC-NSC ATAC')
# title(ylab = 'Neurog2 random genome\nenergy quantiles', line = 4)
# dev.off()

# options(repr.plot.width = 6)
# options(repr.plot.height = 12)

# # ttt <- lapply(motifs_to_take, function(mti) {
# png('./output/sequence_modeling/figs/neurog2_rg_quantiles_delta_ipc_nsc_horiz.png', w = 900, h = 1800)
# mti <- motifs_to_take[[3]]
#     databp = as.data.frame(list(hox = ie_mat_rgq[pb,mti], nsc_lvl = nsc_lvls, delta_lvl = din_lvls))

# par(las = 2, mar = c(15,25,3,5), cex.axis = 3, cex.lab = 4, cex.main = 3)
# boxplot(hox ~ ., horizontal = T, data = databp, xlab = '', ylab = ''
# #         , main = 'Neurog2 random genome energy quantile\nstratified by NSC ATAC and delta IPC-NSC ATAC'
#        )
# title(xlab = 'Neurog2 random genome\nenergy quantiles', line = 12, cex = 2)
# dev.off()

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
#     if (!(has_name(peaks_of_interest, 'peak_name') & has_name(background_peaks, 'peak_name'))) {
#         stop('Peaks of interest or background peaks do not have peak name')
#     }
#     if (length(which(!(peaks_of_interest$peak_name %in% background_peaks$peak_name))) > 0) {
#         stop('Not all peaks of interest are in background peak set')
#     }
#     nei_peaks_peaks <- gintervals.neighbors(peaks_of_interest, background_peaks, maxdist = d_proximity_atac, mindist = -d_proximity_atac, maxneighbors = 1e+4)
#     colnames(nei_peaks_peaks)[grep('peak_name', colnames(nei_peaks_peaks))] <- c('peak_name_1', 'peak_name_2')
#     nei_peaks_peaks <- nei_peaks_peaks[nei_peaks_peaks$peak_name_1 != nei_peaks_peaks$peak_name_2 & abs(nei_peaks_peaks$dist) >= d_puncture,]
#     nei_peaks_tads <- gintervals.neighbors(tads, background_peaks, maxdist = 0, mindist = 0, maxneighbors = 1e+4)
#     nei_peaks_genes <- gintervals.neighbors(background_peaks, tss, maxdist = d_proximity_rna, mindist = -d_proximity_rna, maxneighbors = 1e+4)
#     peaks_tads <- gintervals.neighbors(nei_peaks_genes[,1:4], tads, maxdist = 0, mindist = 0, maxneighbors = 1)
#     genes_tads <- gintervals.neighbors(nei_peaks_genes[,c(5:7, grep('geneSymbol', colnames(nei_peaks_genes)))], 
#                                        tads, maxdist = 0, mindist = 0, maxneighbors = 1)
#     nei_peaks_peaks$tads1 <- nei_peaks_tads$tad_name[match(nei_peaks_peaks[,4], nei_peaks_tads$peak_name)]
#     nei_peaks_peaks$tads2 <- nei_peaks_tads$tad_name[match(nei_peaks_peaks[,8], nei_peaks_tads$peak_name)]
#     if (restrict_to_tads) {
#         nei_peaks_peaks <- nei_peaks_peaks[which(nei_peaks_peaks$tads1 == nei_peaks_peaks$tads2),]
#     }
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
#     # genes_nei_umis <- log2(1e-5 + genes_nei_umis/sum(genes_nei_umis))
#     # genes_nei_umis <- rna_egc_by_day_n[,day]
#     npgf$rna_umis <- genes_nei_umis[npgf$geneSymbol]
#     npgf$rna_umis[is.na(npgf$rna_umis)] <- 0
#     print(head(npgf))
#     print(quantile(npgf[,'rna_umis']))
#     npgf_sum_rna <- tapply(npgf$rna_umis, npgf$peak_name, function(x) sum(x))
#     npgf_sum_rna[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% npgf$peak_name)]] <- 0
#     npgf_sum_rna[is.na(npgf_sum_rna)] <- 0
#     npgf_sum_rna <- npgf_sum_rna[peaks_of_interest$peak_name]
#     # npgf_sum_rna <- log2(1e-5 + npgf_sum_rna/sum(npgf_sum_rna))
#     print(head(npgf_sum_rna))
#     marginal_peaks_umis <- Matrix::rowSums(mca@mat[background_peaks$peak_name,mc_sel])
#     # nei_peaks_peaks$atac_umis_2 <- egc_by_day_n[nei_peaks_peaks$peak_name_2,day]
#     nei_peaks_peaks$atac_umis_2 <- marginal_peaks_umis[match(nei_peaks_peaks$peak_name_2, names(marginal_peaks_umis))]
#     # print(quantile(nei_peaks_peaks$atac_umis_2))
#     print(length(nei_peaks_peaks$atac_umis_2))
#     print(length(nei_peaks_peaks$peak_name_1))
#     tad_atac_umis <- tapply(nei_peaks_peaks$atac_umis_2, nei_peaks_peaks$peak_name_1, sum)
#     # print(quantile(tad_atac_umis))
#     npgf_sum_atac <- tad_atac_umis
#     npgf_sum_atac[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% names(npgf_sum_atac))]] <- 0
#     npgf_sum_atac[is.na(npgf_sum_atac)] <- 0
#     npgf_sum_atac <- npgf_sum_atac[peaks_of_interest$peak_name]
#     # npgf_sum_atac <- log2(1e-5 + npgf_sum_atac/sum(npgf_sum_atac))
#     sum_umi_df <- dplyr::left_join(tibble::enframe(npgf_sum_atac, name = 'peak_name', value = 'atac_umi'), 
#                                    tibble::enframe(npgf_sum_rna, name = 'peak_name', value = 'rna_umi'), by = 'peak_name')
#     # r_vec <- egc_by_day_n[,day]
#     r_vec <- Matrix::rowSums(mca@mat[,mc_sel])
#     r_vec <- log2(1e-5 + r_vec/sum(r_vec))
#     sum_umi_df$r <- r_vec[sum_umi_df$peak_name]
#     if(!is.null(peak_clusters)) {
#         sum_umi_df$cluster <- as.numeric(names(peak_clusters)[match(sum_umi_df$peak_name, peak_clusters)])
#     }
#     if (!is.null(pred)) {sum_umi_df$pred <- pred[sum_umi_df$peak_name]}
#     return(sum_umi_df)
# }


# nm <- 'pl_cort'
# mc_rna <- scdb_mc(nm)
# mat <- scdb_mat(nm)

# tads <- gintervals.load('intervs.global.tad_names')

# # a_legc_max_ct <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, max))
# # a_legc_min_ct <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, min))

# # options(repr.plot.width = 9)
# # options(repr.plot.height = 6)

# # cust_mc_ord_st[names(cust_mc_ord_st) == 'NSC']

# # pheatmap(t(apply(a_legc_avg_cl[,cust_mc_ord_st[names(cust_mc_ord_st) == 'NSC']], 1, zoo::rollmean, k = 10)), cluster_cols = F, col = colorRampPalette(c('white', 'red', 'black'))(100))

# # plot(rowMeans(a_legc[,which(mcmd$cell_type == 'NSC')]), rowSds(a_legc[,which(mcmd$cell_type == 'NSC')]))

# # hist(a_legc_max_ct[,'NSC'] - a_legc_min_ct[,'NSC'])

# # hist(a_legc_max_ct[,'IPC'] - a_legc_min_ct[,'IPC'])



# ds <- c(2e+3,5e+3, 5e+4, 1e+5, 2.5e+5)
# pks <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))
# j <- 2
# pcai <- proximal_chromatin_activity(peaks_of_interest = pks, 
#                                                                     d_puncture = ds[[j-1]],
#                                                                     mca = mca,
#                                                                    background_peaks = pks,
#                                                                    d_proximity = ds[[j]], 
#                                                                    mc_sel = nsc_mcs,
#                                                                         mat_rna = mat_rna,
#                                                                     mc_rna = mc_rna,
#                                                                    restrict_to_tads = T, 
#                                                                     tss = tss,
#                                                                    tads = tads)

# ds <- c(2e+3,5e+3, 5e+4, 1e+5, 2.5e+5)
# pks <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))
# days <- paste0('E', 13:17)
# sud_lst_nsc <- plyr::llply(tail(seq_along(ds), -1), 
#                             function(j) proximal_chromatin_activity(peaks_of_interest = pks, 
#                                                                     d_puncture = ds[[j-1]],
#                                                                     mca = mca,
#                                                                    background_peaks = pks,
#                                                                    d_proximity = ds[[j]], 
#                                                                    mc_sel = nsc_mcs,
#                                                                         mat_rna = mat_rna,
#                                                                     mc_rna = mc_rna,
#                                                                    restrict_to_tads = T, 
#                                                                     tss = tss,
#                                                                    tads = tads), .parallel = T)

# names(sud_lst_nsc) <- tail(ds, -1)

# sud_lst_nsc_df <- tibble::column_to_rownames(as.data.frame(cbind(sud_lst_nsc[[1]][,1], do.call('cbind', lapply(sud_lst_nsc, function(x) x[,2:3])), sud_lst_nsc[[1]][,4])), 'peak_name')

# quantile(a_legc_avg_ct[,'NSC'], (0:20)/20)

# multintersect <- function(...) {
#     args <- list(...)
#     if (length(args) < 2) {stop('Must supply at least two arguments')}
#     temp <- args[[1]]
#     for (j in 2:length(args)) {
#         temp <- intersect(temp, args[[j]])
#     }
#     print(glue::glue('Length of arguments was {paste0(unlist(lapply(args, length)), collapse = ", ")}'))
#     print(glue::glue('Length of intersection is {length(temp)}'))
#     return(temp)
# }

# load('./output/methylation/avg_meth_all.rda')

# head(avg_meth_all)

# pa <- multintersect(rownames(ie_mat), names(delta_ipc_nsc), rownames(avg_meth_all), rownames(sud_lst_nsc_df)
#                     # , rownames(a_legc_avg_ct)[a_legc_avg_ct[,'NSC'] < -16]
#                    )

# dfa <- as.data.frame(cbind(ie_mat[pa,motifs_to_take], rm_ama[pa],  
#                            subset(sud_lst_nsc_df[pa,], select = -c(r)), 
#                            cut(delta_ipc_nsc[pa], 
#                                breaks = quantile(delta_ipc_nsc[pa], c(0,0.5,0.75,0.95,1)), 
#                                # labels = paste0('delta_', c('neg', 'weak', 'pos'))
#                               )))

# colnames(dfa)[[ncol(dfa)]] <- 'delta_atac_lvl'

# splits_motifs <- split(dfa[,1:6],dfa$delta_atac_lvl)
# names(splits_motifs) <- levels(dfa$delta_atac_lvl)

# splits_prox <- split(dfa[,7:(ncol(dfa)-1)],dfa$delta_atac_lvl)
# names(splits_prox) <- levels(dfa$delta_atac_lvl)

# options(repr.plot.width = 24)

# par(las = 2, mar = c(10, 4, 2, 1), mfrow = c(1,4))
# vvv <- lapply(seq_along(splits_motifs), function(i) {boxplot(splits_motifs[[i]], main = names(splits_motifs)[[i]]); grid(col = 'pink', lwd = 1, lty= 2)})
# par(las = 2, mar = c(10, 4, 2, 1), mfrow = c(1,4))
# vvv <- lapply(seq_along(splits_prox), function(i) {boxplot(splits_prox[[i]], main = names(splits_prox)[[i]], ylim = c(-16.8,-14)); grid(col = 'pink', lwd = 2, lty = 2)})

# ng2_mtf <- grep('neurog2', motifs_to_take, v=T, ign = T)

# # train_inds <- train_peaks
# # test_inds <- val_peaks

# options(repr.plot.width = 4)
# options(repr.plot.height = 4)

# # par(mfrow = c(4,4))
# lapply(seq_along(tm_w_add_feat@motif_models), function(i) suppressWarnings(prego::plot_pssm_logo(tm_w_add_feat@motif_models[[i]]$pssm, title = colnames(prego_ie_mat)[names(tm_w_add_feat@motif_models)[[i]]])))

# devtools::load_all('~/src/enhflow/')
# # load('./output/sequence_modeling/enhflow_output_delta_ipc_nsc_n_motifs=8.rda')

# # load('./output/sequence_modeling/enhflow_model_w_add_feat_120923.rda')

# # load('./output/sequence_modeling/enhflow_model_50_motifs_w_hc.rda')
# load('./output/sequence_modeling/enhflow_model_n_motif=16_w_clustering.rda')

# prego_motifs <- do.call('rbind', lapply(seq_along(tm_w_add_feat@motif_models), function(i) dplyr::mutate(tm_w_add_feat@motif_models[[i]]$pssm, motif = names(tm_w_add_feat@motif_models)[[i]])))

# prego_ie <- prego::gextract_pwm(intervals = mcp, motifs = unique(prego_motifs$motif), dataset = prego_motifs)

# prego_ie$cluster <- ifelse(rownames(a_legc) %in% dist_peaks$peak_name, paste0('enh_', km_enh_a_legc$cluster[match(prego_ie$peak_name, rownames(a_legc[dist_peaks$peak_name,]))]),
#                            paste0('prom_', km_prom_a_legc$cluster[match(prego_ie$peak_name, rownames(a_legc[prom_peaks$peak_name,]))]))


# # prego_ie$cluster <- km_a_legc$cluster[match(prego_ie$peak_name, rownames(a_legc))]

# colnames(prego_ie)

# prego_ie_mat <- subset(prego_ie, select = -c(chrom, start, end, peak_name, cluster))

# rownames(prego_ie_mat) <- prego_ie$peak_name

# motif_clustering <- tm_w_add_feat@params$distilled_features

# amd <- prego::all_motif_datasets()

# amou <- unique(amd[,c('motif', 'motif_orig')])

# feats <- names(scdb_gset('pl_cort_dns')@gene_set)

# library(gridExtra)

# options(repr.plot.width = 20)
# options(repr.plot.height = 20)
# par(mfrow = c(4,4))
# do.call('grid.arrange', c(lapply(seq_along(tm_w_add_feat@motif_models), function(i) prego::plot_pssm_logo(tm_w_add_feat@motif_models[[i]]$pssm, title = names(tm_w_add_feat@motif_models)[[i]])), ncol = 4))

# prego_98_bin <- apply(prego_ie_mat, 2, function(x) as.numeric(x >= quantile(x, .98)))

# rownames(prego_98_bin) <- rownames(prego_ie_mat)

# prego_98_num_bin <- tgs_matrix_tapply(t(prego_98_bin[dist_peaks$peak_name,]), km_enh_a_legc$cluster, function(x) length(which(x  == 1)))

# prego_98_lfc <- tgs_matrix_tapply(t(prego_98_bin[dist_peaks$peak_name,]), km_enh_a_legc$cluster, function(x) log2((1+length(which(x == 1)))/(1 + 0.02*length(x))))


# options(repr.plot.width = 12)
# options(repr.plot.height = 12)

# pheatmap::pheatmap(prego_98_lfc[enh_cl_ord,], cluster_rows = F, cluster_cols = T, clustering_method = 'ward.D2', col = clrmp_rel, breaks = seq(-3,3,l=1000))

# ### Fitting using GLMs
# log_norm_feats <- iceqream::create_logist_features(tm_w_add_feat@normalized_energies)

# pc <- intersect(rownames(x_all), rownames(tm_w_add_feat@normalized_energies))
# dp <- delta_ipc_nsc[pc]
# dpl <- (dp - min(dp))/(max(dp) - min(dp))
# # glm_model_lin <- glmnet::glmnet(tm_w_add_feat@normalized_energies[pc,], dpl, binomial(link = "logit"), alpha = 1, lambda = 1e-5, parallel = T, seed = 60427)
# # glm_model_lin <- glmnet::glmnet(log_norm_feats[pc,], dpl, binomial(link = "logit"), alpha = 1, lambda = 1e-5, parallel = T, seed = 60427)
# glm_model_lin <- glmnet::glmnet(prego_ie_mat[pb,], delta_ipc_nsc[pb], alpha = 1, lambda = 1e-5, parallel = T, seed = 60427)

# glm_model_lin$a0

# glm_model_preds <- predict(object = glm_model_lin, newx = as.matrix(x_all[pb,colnames(prego_ie_mat)]))

# smoothScatter(glm_model_preds, delta_ipc_nsc[pb])

# eval_results(glm_model_preds, delta_ipc_nsc[pb])

# cor(glm_model_preds, delta_ipc_nsc[pb])**2

# ttt <- iceqream::filter_traj_model(tm_w_add_feat, 1e-3)

# ttt <- iceqream::filter_traj_model(tm_w_add_feat, 1e-3, parallel = TRUE)

# tm_w_add_feat@params$distilled_features[grep('neuro', tm_w_add_feat@params$distilled_features$feat, ign =T),]

# names(ttt@features_r2[ttt@features_r2 >= 1e-3])

# names(ttt@features_r2[ttt@features_r2 >= 1e-3])

# names(ttt@features_r2[ttt@features_r2 < 1e-3])

# spatial_freqs_300 <- iceqream::compute_traj_model_spatial_freq(tm_w_add_feat, size = 300, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1, parallel = T, atac_track = "mmcortex.marginal")

# spatial_freqs_1000 <- iceqream::compute_traj_model_spatial_freq(tm_w_add_feat, size = 1000, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1, parallel = T, atac_track = "mmcortex.marginal")

# options(repr.plot.width = 16)
# options(repr.plot.height = 16)

# # par(mfrow = c(4,4))
# # uuu <- sapply(unique(spatial_freqs$type), function(typei) {
# #     sapply(unique(spatial_freqs$motif), function(mtfi) {
# #         inds <- which(spatial_freqs$type == typei & spatial_freqs$motif == mtfi)
# #         plot(spatial_freqs$pos[inds], zoo::rollmean(spatial_freqs$freq[inds], k = 20, na.pad = T), main = paste(typei, '\n', mtfi))
# #     })
# # })

# par(mfrow = c(4,4))
# uuu <- sapply(unique(spatial_freqs_1000$type), function(typei) {
#     sapply(unique(spatial_freqs_1000$motif), function(mtfi) {
#         inds <- which(spatial_freqs_1000$type == typei & spatial_freqs_1000$motif == mtfi)
#         plot(spatial_freqs_1000$pos[inds], zoo::rollmean(spatial_freqs_1000$freq[inds], k = 20, na.pad = T), main = paste(typei, '\n', mtfi))
#     })
# })

# options(repr.plot.width = 22)
# options(repr.plot.height = 22)
# par(mfrow = c(5,5))
# uuu <- sapply(unique(spatial_freqs_1000$type), function(typei) {
#     sapply(head(unique(spatial_freqs_1000$motif),25), function(mtfi) {
#         inds <- which(spatial_freqs_1000$type == typei & spatial_freqs_1000$motif == mtfi)
#         plot(spatial_freqs_1000$pos[inds], zoo::rollmean(spatial_freqs_1000$freq[inds], k = 20, na.pad = T), main = paste(typei, '\n', mtfi))
#     })
# })

# plot_traj_model_report(tm_w_add_feat, "./output/sequence_modeling/report-k=16-1000bp/", k = 3, spatial_freqs = spatial_freqs_1000)

# eval_results(predicted = psn, true = tm_w_add_feat@diff_score)

# lapply(seq_along(tm_w_add_feat@motif_models), function(i) prego::plot_pssm_logo(tm_w_add_feat@motif_models[[i]]$pssm, title = names(tm_w_add_feat@motif_models)[[i]]))

# # plot(psn, tm_w_add_feat@diff_score, cex = 0.15, pch = 16)
# # abline(0,1,col = 'red')

# # head(tm_w_add_feat@diff_score)

# # head(psn)

# # all(names(psn) == names(tm_w_add_feat@diff_score))

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)

# # psn <- tm_w_add_feat@predicted_diff_score

# # psn <- (psn - min(psn))/(max(psn) - min(psn))

# # psn <- psn*(max(tm_w_add_feat@diff_score) - min(tm_w_add_feat@diff_score)) + min(tm_w_add_feat@diff_score)

# # head(tm_w_add_feat@diff_score)

# # eval_results(tm_w_add_feat@predicted_diff_score, tm_w_add_feat@diff_score)

# # eval_results(tm@predicted_diff_score, tm@diff_score)

# # cor(tm_w_add_feat@diff_score, tm_w_add_feat@predicted_diff_score,  method = 'spearman')

# # sqrt(cor(tm_w_add_feat@diff_score, tm_w_add_feat@predicted_diff_score,  method = 'spearman'))

# calc_q_enrichment <- function(mat, cl = NULL, q = 0.98) {
#     # mat is a motif energy matrix (peaks x motifs)
#     # cl is the output of tglkmeans::TGL_kmeans
#     # q is a quantile
#     qs <- unlist(plyr::llply(colnames(mat), function(x) quantile(mat[,x], probs = q), .parallel = T))
#     q_bin_mat <- t(plyr::laply(1:length(qs), function(i) as.numeric(mat[,i] >= qs[[i]]), .parallel = T))
#     colnames(q_bin_mat) <- colnames(mat)
#     rownames(q_bin_mat) <- rownames(mat)
#     q_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(q_bin_mat)), cl$cluster, sum)
#     q_lfc <- log2(1e-2+q_sum_clust_int/(0.02*as.numeric(cl$size)))
#     colnames(q_lfc) <- colnames(q_bin_mat)
#     return(q_lfc)
# }

# prego_lfc_cl <- calc_q_enrichment(prego_ie_mat, km_a_legc)

# options(repr.plot.width = 9)
# options(repr.plot.height = 14)

# pheatmap::pheatmap(prego_lfc_cl[peak_clust_var[hcct2],],
#                   cluster_rows = F, cluster_cols = F,
#                            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), fontsize_col = 10)

# pim_q98 <- apply(prego_ie_mat, 2, quantile, 0.98, na.rm = T)
# pim_98_bin <- t(plyr::laply(names(pim_q98), function(i) as.numeric(prego_ie_mat[,i] >= pim_q98[[i]]), .parallel = F))
# pim_98_sum_clust <- tgs_matrix_tapply(t(as.matrix(pim_98_bin)), km_a_legc$cluster, sum)
# pim_98_lfc <- log2(1e-1+pim_98_sum_clust/(0.02*as.numeric(km_a_legc$size)))

# colnames(pim_98_lfc) <- names(pim_q98)

# pltmtd <- cbind(pim_98_lfc, tapply(delta_ipc_nsc, km_a_legc$cluster, mean))
# ppp <- pheatmap::pheatmap(pim_98_lfc[peak_clust_var,], 
#                           annotation_color = list(value = setNames(colorRampPalette(c('blue3', 'white', 'red3'))(100), seq(-1,1,l=100))),
#                                                                                                 annotation_row = tibble::column_to_rownames(tibble::enframe(tapply(delta_ipc_nsc, km_a_legc$cluster, mean)), 'name'), treeheight_col = 0, treeheight_row = 0, clustering_method = 'ward.D', col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), fontsize_col = 16)

# options(repr.plot.width = 8)
# options(repr.plot.height = 8)

# sqrt(0.49)

# ppp <- pheatmap::pheatmap(pim_98_lfc[c(11,22,23,1,2,4),], col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), fontsize_col = 16)

# feats_plot <- importance$Feature[importance$Gain >= 0.01]

# x_all_new <- cbind(x_all, rep(NA, nrow(x_all)))
# rnin <- intersect(names(nsc_ipc_pks_ord2), rownames(x_all_new))
# x_all_new[rnin,ncol(x_all_new)] <- nsc_ipc_pks_ord2[rnin]

# cn <- colnames(x_all_new)
# cn[[length(cn)]] <- 'pk_cl'
# colnames(x_all_new) <- cn

# library(vioplot)

# mfp <- amou$motif_orig[match(feats_plot, amou$motif)]
# fp_titles <- ifelse(!is.na(mfp), unlist(purrr::map(stringr::str_split(mfp, '_'), 1)), feats_plot)
# fp_titles

# col_function <- function(xvec, groupvec, clrmp) {
#     meds <- tapply(xvec, groupvec, mean, na.rm = T)
#     colvals <- clrmp[1+round((length(clrmp)-2)*(meds - min(xvec, na.rm = T))/(max(xvec, na.rm = T) - min(xvec, na.rm = T)))]
#     return(colvals)
# }

# options(repr.plot.width = 14)
# options(repr.plot.height = 14)

# slfp <- sqrt(length(feats_plot))
# par(mfrow = c(round(slfp), ceiling(slfp)), cex.lab = 2, mar = c(6,6,1,1))
# pks_here <- intersect(rownames(x_all), names(delta_ipc_nsc))
# vvv <- sapply(feats_plot, function(fti) {
#     smoothScatter(x_all_new[pks_here,fti], delta_ipc_nsc[pks_here], xlab = fti, ylab = 'Delta ATAC')
#     grid(lwd = 2, lty = 2)
#     text(quantile(x_all_new[pks_here,fti], 0.), 
#         quantile(delta_ipc_nsc[pks_here], .999),
#         labels = paste('cor =', round(cor(delta_ipc_nsc[pks_here], x_all_new[pks_here,fti]), 2)), cex = 2, adj = c(0, 1),
#         col = 'red')

# })

# slfp <- sqrt(length(feats_plot))
# par(mfrow = c(round(slfp), ceiling(slfp)), cex.lab = 2, mar = c(6,6,1,1))
# pks_here <- intersect(rownames(x_all), names(delta_ipc_nsc))
# vvv <- sapply(feats_plot, function(fti) {
#     plot(x_all_new[pks_here,fti], delta_ipc_nsc[pks_here], xlab = fti, ylab = 'Delta ATAC', pch = 16, cex = 0.1)
#     grid(lwd = 2, lty = 2)
#     text(quantile(x_all_new[pks_here,fti], 0), 
#         quantile(delta_ipc_nsc[pks_here], .999),
#         labels = paste('cor =', round(cor(delta_ipc_nsc[pks_here], x_all_new[pks_here,fti]), 2)), cex = 2, adj = c(0, 1),
#         col = 'red')

# })

# png('./output/mcatac/figs/nsc_ipc_all_peaks_all_clusters_delta_atac_features_vioplots_test.png', h = 2000, w = 2000, res = 100)
# par(mfrow = c(1,12), mar = c(10,3,4,1), cex.lab = 2, cex.axis = 2, cex.main = 2)
# ylimi <- c(0.85,length(cln_plot)+0.15)
# pks_here <- intersect(names(nsc_ipc_pks_ord), rownames(x_all_new))

# vioplot(delta_ipc_nsc[pks_here] ~ nsc_ipc_pks_ord2[pks_here], xlab = '', ylab = '', main = '',col = col_function(delta_ipc_nsc[pks_here], nsc_ipc_pks_ord2[pks_here], clrmp_rel2), horizontal = T, xlim = ylimi)

# title(xlab = substitute(paste(bold('Delta\nATAC'))), line = 5)
# ttt <- sapply(seq_along(feats_plot), function(i) {
#     fti <- feats_plot[[i]]
#     vioplot(x_all_new[pks_here,fti] ~ x_all_new[pks_here,'pk_cl'], xlab = '',ylab = '', main = '',col = col_function(x_all_new[pks_here,fti], nsc_ipc_pks_ord2[pks_here], clrmp_rel2), horizontal = T, xlim = ylimi)
#     title(xlab = fp_titles[[i]], line = 5)
# })

# # vioplot(rm_ama[pks_here] ~ nsc_ipc_pks_ord2, xlab = '',ylab = '', main = '',col = col_function(rm_ama[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = T, xlim = ylimi)
# # title(xlab = substitute(paste(bold('Methylation'))), line = 5)
# # vioplot(dinucs_per_peak[pks_here,'CG'] ~ nsc_ipc_pks_ord2, xlab = '',ylab = '', main = '', col = col_function(nuc_per_peak[pks_here,'C'], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = T, xlim = ylimi)
# # title(xlab = substitute(paste(bold('CpG \ncontent'))), line = 5)
# # cns <- colnames(prego_ie_mat)
# # cnlst <- unlist(purrr::map(stringr::str_split(string = cns, pattern = '\\.'), 2))
# # yyy <- sapply(, function(i) {vioplot(prego_ie_mat[pks_here,cns[[i]]] ~ nsc_ipc_pks_ord2, ylab = '', main = '', xlab = '', 
# #                                                                    col = col_function(prego_ie_mat[pks_here,cns[[i]]], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = T, xlim = ylimi);
# #                                                                    title(xlab = cnlst[[i]], line = 5)})
# # vioplot(atac_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, xlab = '', ylab = '',main = '', col = col_function(atac_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = T, xlim = ylimi)
# # title(xlab = substitute(paste(bold('ATAC UMI\n(50kbp radius)'))), line = 5)
# # vioplot(rna_umi_vec[pks_here] ~ nsc_ipc_pks_ord2, xlab = '', ylab = '',main = '',col = col_function(rna_umi_vec[pks_here], nsc_ipc_pks_ord2, clrmp_rel2), horizontal = T, xlim = ylimi)
# # title(xlab =substitute(paste(bold('RNA UMI\n(50kbp radius)')))tm_w_, line = 5)
# dev.off()

# # plot(as.data.frame(x_all), pch = 16, cex = 0.25)

# options(repr.plot.width = 6)
# options(repr.plot.height = 6)
# scs_here <- intersect(names(rm_ama), rownames(preds_all))

# # par(mfrow = c(1,2), cex.lab = 2, mar = c(5,5,2,1))
# png('./output/sequence_modeling/figs/ipc_vs_nsc_atac_color_by_meth.png', h = 1000, w = 1000, res = 100)
# par(cex.lab = 4, cex.axis = 2, mar = c(7,7,2,1))
# # plot(xnsc[scs_here], xipc[scs_here], pch = 16, cex = 0.32, xlab = 'NSC ATAC', ylab = 'IPC ATAC')
# vh <- rm_ama[scs_here]
# plot(xnsc[scs_here], xipc[scs_here],  cex = 0.52, lwd = 0.5, xlab = '', ylab = '')
# points(xnsc[scs_here], xipc[scs_here],cex = 0.52,pch = 16,
#          col = clrmp_abs[1+round(999*(vh - min(vh))/(max(vh) - min(vh)))])
# title(xlab = 'NSC ATAC', line = 4)
# title(ylab = 'IPC ATAC', line = 4)
# dev.off()

# plot_color_bar(vals = seq(0,1, l = 101), cols = clrmp_abs[seq(1,length(clrmp_abs),l=100)], fig_fn = './output/sequence_modeling/figs/meth_color_bar.png', title = 'methylation')

# # options(repr.plot.width = 12)
# # options(repr.plot.height = 6)
# # par(mfrow = c(1,2), cex.lab = 2, mar = c(5,5,2,1))
# # plot(xnsc[scs_here], xipc[scs_here], pch = 16, cex = 0.32, xlab = 'NSC ATAC', ylab = 'IPC ATAC')
# # vh <- delta_ipc_nsc[scs_here]
# # plot(xnsc[scs_here], xipc[scs_here],  cex = 0.52, lwd = 0.5, xlab = 'NSC ATAC', ylab = 'IPC ATAC')
# # points(xnsc[scs_here], xipc[scs_here],cex = 0.52,pch = 16,
# #          col = clrmp_rel[1+round(999*(vh - min(vh))/(max(vh) - min(vh)))])

# # par(mfrow = c(1,2))
# # plot(xnsc[scs_here], rm_ama[scs_here], pch = 16, cex = 0.2)
# # smoothScatter(xnsc[scs_here], rm_ama[scs_here])
# # plot(xipc[scs_here], rm_ama[scs_here], pch = 16, cex = 0.2)
# # vh <- rm_ama[scs_here]
# # plot(xnsc[scs_here], xipc[scs_here], pch = 16, cex = 0.72,
# #          col = clrmp_abs[1+round(999*(vh - min(vh))/(max(vh) - min(vh)))])

# # options(repr.plot.width = 12)
# # options(repr.plot.height = 12)


# # mat_for_cor <- cbind(as.matrix(preds_all[scs_here,]), delta_ipc_nsc[scs_here], xnsc[scs_here], xipc[scs_here])
# # cn <- colnames(mat_for_cor)
# # cn[(length(cn)-2):length(cn)] <- c('delta_ipc_nsc', 'nsc_atac', 'ipc_atac')
# # colnames(mat_for_cor) <- cn
# # pheatmap::pheatmap(tgs_cor(mat_for_cor, spearman = T), col = clrmp_rel, breaks = seq(-1,1,l = length(clrmp_rel)),clustering_method = 'ward.D2')

# options(repr.plot.width = 12)
# options(repr.plot.height = 6)
# png('./output/sequence_modeling/figs/delta_atac_vs_meth_rollmean_shap_and_delta.png', h = 400, w = 900)
# par(mfrow = c(1,2), cex.lab = 2, cex.axis = 1.5, mar = c(5,5,1,1))
# plot(rm_ama[scs_here], delta_ipc_nsc[scs_here], pch = 16, cex = 0.2, ylim = 0.8*c(-1,1), xlab = 'methylation', ylab = 'Delta IPC-NSC ATAC')
# lines(sort(rm_ama[scs_here]), zoo::rollmean(delta_ipc_nsc[names(sort(rm_ama[scs_here]))], k = 100, na.pad = T), col = 'red')
# lines(sort(rm_ama[scs_here]), preds_all[names(sort(rm_ama[scs_here])),'methylation'], col = 'green')
# legend('bottomright', legend = c('methylation SHAP value', 'rollmean delta ATAC'), col = c('green', 'red'), lwd = rep(1,2))
# smoothScatter(rm_ama[scs_here], delta_ipc_nsc[scs_here], ylim = 0.8*c(-1,1), xlab = 'methylation', ylab = 'Delta IPC-NSC ATAC')
# lines(sort(rm_ama[scs_here]), zoo::rollmean(delta_ipc_nsc[names(sort(rm_ama[scs_here]))], k = 100, na.pad = T), col = 'red')
# lines(sort(rm_ama[scs_here]), preds_all[names(sort(rm_ama[scs_here])),'methylation'], col = 'green')
# legend('bottomright', legend = c('methylation SHAP value', 'rollmean delta ATAC'), col = c('green', 'red'), lwd = rep(1,2))
# dev.off()

# options(repr.plot.width = 12)
# options(repr.plot.height = 6)
# # png('./output/sequence_modeling/figs/delta_atac_vs_olig2_rollmean_shap_and_delta.png', h = 400, w = 900)
# par(mfrow = c(1,2), cex.lab = 2, cex.axis = 1.5, mar = c(5,5,1,1))
# feati <- 'HOMER.Olig2'

# plot(x_all[scs_here,feati], delta_ipc_nsc[scs_here], pch = 16, cex = 0.2, ylim = 0.8*c(-1,1), xlab = 'methylation', ylab = 'Delta IPC-NSC ATAC')
# lines(sort(x_all[scs_here,feati]), zoo::rollmean(delta_ipc_nsc[scs_here[order(x_all[scs_here,feati])]], k = 100, na.pad = T), col = 'red')
# lines(sort(x_all[scs_here,feati]), preds_all[scs_here[order(x_all[scs_here,feati])],feati], col = 'green')
# legend('bottomright', legend = c('HOMER.Olig2 SHAP value', 'rollmean delta ATAC'), col = c('green', 'red'), lwd = rep(1,2))
# smoothScatter(x_all[scs_here,feati], delta_ipc_nsc[scs_here], ylim = 0.8*c(-1,1), xlab = 'methylation', ylab = 'Delta IPC-NSC ATAC')
# lines(sort(x_all[scs_here,feati]), zoo::rollmean(delta_ipc_nsc[scs_here[order(x_all[scs_here,feati])]], k = 100, na.pad = T), col = 'red')
# lines(sort(x_all[scs_here,feati]), preds_all[scs_here[order(x_all[scs_here,feati])],feati], col = 'green')
# legend('bottomright', legend = c('HOMER.Olig2 SHAP value', 'rollmean delta ATAC'), col = c('green', 'red'), lwd = rep(1,2))
# # dev.off()

# # plot(rm_ama[scs_here], delta_ipc_nsc[scs_here], pch = 16, cex = 0.2, ylim = 0.5*c(-1,1))
# # points(sort(rm_ama[scs_here]), delta_ipc_nsc[names(sort(rm_ama[scs_here]))], col = 'red')
# # points(sort(rm_ama[scs_here]), preds_all[names(sort(rm_ama[scs_here])),'methylation'], col = 'green')
# # smoothScatter(rm_ama[scs_here], delta_ipc_nsc[scs_here], ylim = 0.5*c(-1,1))
# # points(sort(rm_ama[scs_here]), delta_ipc_nsc[names(sort(rm_ama[scs_here]))], col = 'red')
# # points(sort(rm_ama[scs_here]), preds_all[names(sort(rm_ama[scs_here])),'methylation'], col = 'green')

# # options(repr.plot.width = 18)
# # options(repr.plot.height = 6)

# # colnames(x_all)

# # head(apply(abs(preds_all), 1, function(x) which(x == max(x))))

# # tbl1 <- sort(table(colnames(preds_all)[apply(abs(preds_all), 1, function(x) which(x == max(x)))]))/nrow(preds_all)

# # tbl1[colnames(preds_all)[!(colnames(preds_all) %in% names(tbl1))]] <- 0

# # cumsum(rev(tbl1))

# # meth_max_inds <- which(colnames(preds_all)[apply(abs(preds_all), 1, function(x) which(x == max(x)))] == 'methylation')

# # preds_all_punctured_meth <- subset(preds_all[meth_max_inds,], select = -methylation)
# # meth_max_next_max <- apply(abs(preds_all_punctured_meth), 1, function(x) which(x == max(x)))

# # tbl2 <- sort(table(colnames(preds_all_punctured_meth)[meth_max_next_max]))/nrow(preds_all_punctured_meth)
# # # tbl2[colnames(preds_all)[!(colnames(preds_all) %in% names(tbl2))]] <- 0
# # tbl2

# # meth_max_next_val <- sapply(seq_along(meth_max_inds), function(i) preds_all_punctured_meth[i, meth_max_next_max[[i]]])

# # options(repr.plot.width = 14)
# # options(repr.plot.height = 7)


# # par(mfrow = c(1,2))
# # plot(preds_all[meth_max_inds, 'methylation'], meth_max_next_val)
# # plot(preds_all[, 'methylation'], preds_all[, 'HOMER.Olig2'])

# # par(mfrow = c(1,2))
# # smoothScatter(preds_all[meth_max_inds, 'methylation'], meth_max_next_val)
# # smoothScatter(preds_all[, 'methylation'], preds_all[, 'HOMER.Olig2'])

# # preds_all_meth_tbl_mat <- as.matrix(table(preds_all[meth_max_inds, 'methylation'] <= 0, meth_max_next_val <= 0))
# # preds_all_meth_tbl_mat

# # fisher.test(preds_all_meth_tbl_mat)

# # options(repr.plot.width = 12)
# # options(repr.plot.height = 12)


# # meth_pred_fct <- cut(x_all[rownames(preds_all),'methylation'],breaks= c(-1e-2, 0.1, 0.5, 1))

# # cor_mats_by_meth_stratum <- lapply(levels(meth_pred_fct), function(lvli) {tgs_cor(preds_all[which(meth_pred_fct == lvli),strong_feats], spearman = T)})

# # table(meth_pred_fct)

# # par(mfrow = c(1,3))
# # uuu <- sapply(levels(meth_pred_fct), function(lvli) {x <- preds_all[which(meth_pred_fct == lvli),'methylation'];
# #                                                      y <- preds_all[which(meth_pred_fct == lvli),'HOMER.Olig2'];
# #                                                      smoothScatter(x, y, xlim = c(-0.7,0.7), ylim = c(-0.7,0.7));
# #                                                      abline(0,0,col='red',lty=2);
# #                                                      lines(c(0,0),c(-10,10),col='red',lty=2);
# #                                                     text(-0.5, 0.5, labels = round(cor(x,y,method ='spearman'), 3), cex = 2)})

# options(repr.plot.width = 8)
# options(repr.plot.height = 8)


# p1 <- pheatmap::pheatmap(cor_mats_by_meth_stratum[[1]], fontsize = 16, col = clrmp_rel, breaks = seq(-1,1, l = length(clrmp_rel)), clustering_method = 'ward.D2')
# p2 <- pheatmap::pheatmap(cor_mats_by_meth_stratum[[2]][p1$tree_row$order, p1$tree_row$order], fontsize = 16, cluster_rows = F, cluster_cols = F, col = clrmp_rel, breaks = seq(-1,1, l = length(clrmp_rel)), clustering_method = 'ward.D2')
# p3 <- pheatmap::pheatmap(cor_mats_by_meth_stratum[[3]][p1$tree_row$order, p1$tree_row$order], fontsize = 16, cluster_rows = F, cluster_cols = F, col = clrmp_rel, breaks = seq(-1,1, l = length(clrmp_rel)), clustering_method = 'ward.D2')

# meth_pred_fct <- cut(preds_all[rownames(preds_all),'methylation'],breaks= c(-2, 0,2))

# table(meth_pred_fct)

# boxplot_df <- tibble::rownames_to_column(as.data.frame(subset(preds_all,select=-BIAS)))
# boxplot_df$meth_level <- meth_pred_fct
# boxplot_df <- tidyr::pivot_longer(boxplot_df, cols = colnames(subset(preds_all, select = -BIAS)))
# # boxplot_df <- tidyr::pivot_longer(data = boxplot_df, )

# boxplot_df <- tibble::rownames_to_column(as.data.frame(x_all[rownames(preds_all),]))
# boxplot_df$meth_level <- meth_pred_fct
# boxplot_df <- tidyr::pivot_longer(boxplot_df, cols = colnames(x_all))
# # boxplot_df <- tidyr::pivot_longer(data = boxplot_df, )


# # length(intersect(inc_peaks, var_peaks))
# # length(intersect(inc_peaks, non_var_peaks))
# # length(intersect(dec_peaks, var_peaks))
# # length(intersect(dec_peaks, non_var_peaks))
# # dec_peaks <- intersect(names(which(delta_ipc_nsc < 0)), rownames(x_all))

# xgb_cv_res <- run_xgb_cv(as.matrix(x_all), y_all, num_peak_groups = 11)

# save(xgb_cv_res, file = './output/sequence_modeling/xgb_cv_res.rda')

# xgb_cv_res_nsc <- run_xgb_cv(as.matrix(x_all_nsc_e13), y_all_nsc_e13, num_peak_groups = 5)

# preds_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['contrib_meth']]))

# preds_all_nsc_e13 <- do.call('rbind', lapply(xgb_cv_res_nsc, function(x) x[['contrib_meth']]))

# ev_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['ev_xgb_meth']]))
# ev_all

# ev_all <- do.call('rbind', lapply(xgb_cv_res, function(x) x[['ev_xgb_meth']]))
# ev_all

# xgb_cv_res <- run_xgb_cv(as.matrix(x_all), y_all, num_peak_groups = 11)

# pred_meth_all <- do.call('c', lapply(xgb_cv_res, function(x) x[['pred_meth']]))
# # png('./output/sequence_modeling/figs/xgboost_predicted_vs_observed_delta_atac_peak_clust_var_new.png', h = 500, w = 600)
# pdf('./output/sequence_modeling/figs/xgboost_predicted_vs_observed_delta_atac_peak_clust_var_new.pdf',
#          h = 500/71, w = 600/71)
# par(cex.lab = 2, cex.main = 2, cex.axis = 1.5, mar = c(5,5,2,1))
# plot(delta_ipc_nsc[names(pred_meth_all)], pred_meth_all, pch = 16, cex = 0.25,
#              ylab = 'Predicted delta ATAC', xlab = 'Observed delta ATAC', 
#                 main = 'xgboost prediction of IPC-NSC ATAC')
# abline(0,1,col = 'red', lwd = 2)
# text(-1.2, 1.2, labels = paste0('R^2 = ', round(mean(ev_all$Rsquare), 2)), cex = 2)
# dev.off()

# fs <- sort(apply(preds_all, 2, sd), decreasing = T)
# fs <- fs[names(fs) != 'BIAS']
# feats_sort <- factor(names(fs), levels = names(fs))
#                      # , levels = 1:length(fs))

# fs/sum(fs)

# pvl <- tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(subset(preds_all,select=-BIAS))), cols = colnames(subset(preds_all, select = -BIAS)))

# feats_sort

# pvl$lvl <- feats_sort[match(pvl$name, levels(feats_sort))]

# x_all_long <- tidyr::pivot_longer(tibble::rownames_to_column(x_all), cols = colnames(x_all))

# mtch_inds_x_all_preds <- match(apply(pvl[,1:2], 1, paste0, collapse = '-'), apply(x_all_long[,1:2], 1, paste0, collapse = '-'))

# pvl$feature_value <- x_all_long$value[mtch_inds_x_all_preds] 

# pvl <- pvl %>% dplyr::group_by(name) %>% 
#             mutate(val_lin = (feature_value  - min(feature_value))/(max(feature_value)-min(feature_value)),
#                                                 color = clrmp_rel[1+round(length(clrmp_rel)*val_lin)])
#             # mutate(val_lin = ifelse(feature_value >= 0, 0.5+0.5*feature_value/max(feature_value), 0.5*(feature_value - min(feature_value))/(-min(feature_value))),
#             #                                     color = clrmp_rel[1+round(length(clrmp_rel)*val_lin)])

# library(beeswarm)

# integral_over_feature_contrib <- colSums(abs(subset(preds_all, select = -BIAS)))
# norm_integral_over_feature_contrib <- integral_over_feature_contrib/sum(integral_over_feature_contrib)

# samp_df <- dplyr::sample_frac(pvl, 0.1)

# # png('./output/sequence_modeling/figs/xgb_cv_res_shap_per_feat_beeswarm.png', h = 1050, w = 2500)
# pdf('./output/sequence_modeling/figs/xgb_cv_res_shap_per_feat_beeswarm.pdf', h = 1050/71, w = 2500/71)
# par(las = 2, mar = c(14,12,2,1), cex.axis = 3, cex.lab = 6)
# # inds <- which(pvl$lvl %in% levels(pvl$lvl)[1:2])
# beeswarm(value ~ lvl, 
#          method = 'compactswarm',
#          spacing = 0.2,
#          corral = 'wrap',
#          data = samp_df,
#          pwcol = samp_df$color, 
#          ylim = c(-0.9,0.55), ylab = '', xlab = '')
# title( ylab = 'SHAP value', line = 7)
# grid(lwd = 5)
# v <- norm_integral_over_feature_contrib
# text(match(names(v), levels(pvl$lvl))-0.45, -0.25, col = 'black', labels = round(v,2), cex = 3, srt = 90)

# dev.off()

# meth_fct <- setNames(factor(c(rep('meth_low', length(dec_peaks)), rep('meth_mid',length(mid_peaks)),
#                         rep('meth_high', length(inc_peaks))), 
#                         levels = c('meth_low', 'meth_mid', 'meth_high')), 
#                         c(dec_peaks, mid_peaks, inc_peaks))

# NSC_0_THRESH <- -16.4
# meth_fct_filt_nsc <- meth_fct[xnsc[names(meth_fct)] <= NSC_0_THRESH]

# olig2_vec <- setNames(x_all[names(meth_fct_filt_nsc),'E_box_1'], names(meth_fct_filt_nsc))
# olig2_cut <- setNames(cut(olig2_vec, quantile(olig2_vec, (0:3)/3), labels = c('E_box_1_low', 'E_box_1_mid','E_box_1_high')), names(meth_fct_filt_nsc))

# options(repr.plot.width = 11)
# options(repr.plot.height = 6)
# # png('./output/sequence_modeling/figs/meth_vs_E_box_1_bpxlot_vec.png', h = 400, w = 650)
# pdf('./output/sequence_modeling/figs/meth_vs_E_box_1_bpxlot_vec.pdf', h = 400/71, w = 650/71)
# par(las = 2, mar = c(8, 5, 3, 1))
# boxplot_vec(xvec = olig2_vec, yvec = rm_ama[names(olig2_vec)], 
#             nm = 'methylation vs E_box_1 energy', num_bins = 10, 
#             ylab = 'NSC methylation', xlab = '')
# title(xlab = 'E_box_1 energy', line = 6)
# dev.off()

# feature_cut <- function(feature_vec, cells, name, cut_values = NULL, cut_quantiles = NULL) {
#     if (is.null(names(feature_vec))) {
#         stop('feature_vec should be a named vector')
#     }
#     if (!is.null(cut_values)) {
#         feature_vec_cut <- setNames(cut(feature_vec, breaks = cut_values, labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
#                                names(feature_vec))
#     } else if (!is.null(cut_quantiles)) {
#         feature_vec_cut <- setNames(cut(feature_vec, breaks = quantile(feature_vec, cut_quantiles), labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
#                                names(feature_vec))
#     }
#     else {
#         feature_vec_cut <- setNames(cut(feature_vec, breaks = c(quantile(feature_vec, 0)-1e-3, quantile(feature_vec, (1:3)/3)), labels = paste0(glue::glue('{name}_'), c('low', 'mid', 'high'))),
#                                names(feature_vec))
#     }
#     # print(length(which(is.na(feature_vec_cut))))
#     return(feature_vec_cut)
# }

# length(pb)

# meth_cells <- pb
# meth_cut <- feature_cut(feature_vec = rm_ama[meth_cells], cells = meth_cells, cut_values = c(-1e-2,0.25,0.75,1), name = 'meth')
# prox_atac_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'prox_ATAC'], names(meth_cut)), cells = names(meth_cut), name = 'prox_ATAC')
# prox_rna_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'prox_RNA'], names(meth_cut)), cells = names(meth_cut), name = 'prox_RNA')
# olig2_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'E_box_1'], names(meth_cut)), cells = names(meth_cut), name = 'E_box_1')
# sox_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'SOX'], names(meth_cut)), cells = names(meth_cut), name = 'SOX')
# eomes_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'T_box_1'], names(meth_cut)), cells = names(meth_cut), name = 'T_box_1')
# nsc_atac_cut <- feature_cut(feature_vec = setNames(x_all[names(meth_cut),'NSC_ATAC'], names(meth_cut)), cells = names(meth_cut), name = 'NSC_ATAC')

# # feat_cut_df <- as.data.frame(cbind(meth_cut, prox_atac_cut, prox_rna_cut, olig2_cut, sox_cut))
# feat_cut_df <- as.data.frame(meth_cut)
# feat_cut_df$prox_atac_cut <- prox_atac_cut
# feat_cut_df$prox_rna_cut <- prox_rna_cut
# feat_cut_df$olig2_cut <- olig2_cut
# feat_cut_df$sox_cut <- sox_cut
# feat_cut_df$eomes_cut <- eomes_cut
# feat_cut_df$nsc_atac_cut <- nsc_atac_cut
# rownames(feat_cut_df) <- names(meth_cut)

# boxplot(delta_ipc_nsc[rownames(feat_cut_df)] ~ feat_cut_df[,'prox_atac_cut']*feat_cut_df[,'olig2_cut'])

# boxplot_and_ks <- function(factor_df, factor1, factor2, y_vec, ks.alternative = 'greater', 
#                            main = NULL, legend_labels = NULL, col_axis = NULL, col_boxplot = NULL, legend_x = NULL, legend_y = NULL, ylab = NULL, ylim = NULL, xlim = NULL) {
#     if (!is.null(col_axis)) {xaxti = 'n'} else {xaxti = 's'}
#     if (!is.null(col_boxplot)) {colbx = col_boxplot} else {colbx = 'lightgray'}
#     filt_inds <- factor_df[,names(factor1)] %in% levels(factor_df[,names(factor1)])[factor1[[1]]] & factor_df[,names(factor2)] %in% levels(factor_df[,names(factor2)])[factor2[[1]]]
#     factor_df <- factor_df[filt_inds,]
#     y_vec <- y_vec[rownames(factor_df)]
#     boxplot(y_vec ~ factor_df[,names(factor1)]*factor_df[,names(factor2)], main = main, xlab = '', xaxt = xaxti, col = colbx, ylab = ylab, ylim = ylim, xlim = xlim, outcex = 0.5)
#     fct_eg <- expand.grid(levels(factor_df[,names(factor1)]),
#                       levels(feat_cut_df[,names(factor2)]))
#     print(fct_eg)
#     eg <- expand.grid(head(levels(factor_df[,names(factor1)])[factor1[[1]]], -1), tail(levels(factor_df[,names(factor1)])[factor1[[1]]], -1),
#                       levels(feat_cut_df[,names(factor2)])[factor2[[1]]])
#     eg <- eg[as.character(eg[,1]) != as.character(eg[,2]),]
#     # for (i in seq_along(col_axis)) {axis(1, at = i, labels = paste0(fct_eg[i,], collapse = '\\.'),col = col_axis[[i]])}
#     Map(axis, side=1, at=1:nrow(fct_eg), col.axis=col_axis, labels=fct_eg[,2], lwd=0, las=2)
#     axis(1,at=1:nrow(fct_eg),labels=FALSE)
#     eg_inds <- cbind(match(eg[,1], levels(factor_df[,names(factor1)])),
#                     match(eg[,2], levels(factor_df[,names(factor1)]))) + 
#                     (match(eg[,3], levels(factor_df[,names(factor2)])) - 1)*length(factor1[[1]])
#     eg$ks_p_res <- p.adjust(apply(eg, 1, function(x) {
#         fct1 <- factor_df[,names(factor1)]
#         fct2 <- factor_df[,names(factor2)]
#         pks1 <- rownames(factor_df)[fct1 == x[[1]][[1]] & fct2 == x[[3]][[1]]]
#         pks2 <- rownames(factor_df)[fct1 == x[[2]][[1]] & fct2 == x[[3]][[1]]]
#         ks_p <- round(ks.test(y_vec[pks1],y_vec[pks2], alternative = ks.alternative)$p.value, 3)
#         return(ks_p)
#     }))
#     if (is.null(legend_y)) {legend_y <- 2*quantile(max(y_vec, na.rm = T))}
#     if (is.null(legend_x)) {legend_x <- 1}
#     legend(x = legend_x, y = legend_y, xpd = T, legend = levels(factor_df[,names(factor1)]), cex = 1.5, col = unique(col_boxplot), pch = rep(15,3))
    
    
    
#     get_num_asterisks <- function(pv_vec) {
#         powvec <- 5*10**seq(-4,0,1)
#         print(pv_vec)
#         npi <- 4 - sapply(pv_vec, function(x) which.max(powvec > x))
#         print(npi)
#         npi <- ifelse(npi < 1, 0, npi)
#         npi <- ifelse(npi > 3, 3, npi) 
#     return(npi)
#                                       }

#     plot_asterisks <- function(x, y, npj) {
#         if (npj > 1) {
#             xl <- seq(from = x - 0.2*npj, to = x + 0.2*npj, length.out = npj)
#         } else if (npj == 1) {xl <- x}
#         else {xl <- NULL}
#         points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = 2.5, col = 'red', lwd= 1)
#     }

#     npi <- get_num_asterisks(eg$ks_p_res)
#     print(npi)
#     for (i in 1:nrow(eg_inds)) {
#         lines(eg_inds[i,1:2], rep(1.175+0.125*(i%%3),2), lwd = 3.5, col = 'red')
#         plot_asterisks(x = mean(eg_inds[i,1:2]), y = 1.25 +0.125*(i%%3), npj = npi[[i]])
#     }

#     print(eg)
#     # print(apply(eg, 2, function(x)  as.character))
# }

# options(repr.plot.width = 13)
# options(repr.plot.height = 13)
# # png('./output/sequence_modeling/figs/delta_atac_vs_homer.olig2_and_meth_and_nsc_atac_levels.png', h = 600, w = 850, res = 100)
# pdf('./output/sequence_modeling/figs/delta_atac_vs_homer.olig2_and_meth_and_nsc_atac_levels.pdf', h = 500/71, w = 850/71)
# par(mfrow = c(1,3),las = 1, mar = c(12,5,8,2), cex.lab = 1.2, cex.main = 1.32, cex.axis = 1.5)
# boxplot_and_ks(feat_cut_df, factor1 = list('meth_cut' = c(1:3)), factor2 = list('olig2_cut' = 2:3),ks.alternative = 'less', 
#                y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '', col_boxplot = rep(adjustcolor(c('lightblue', 'blue2', 'blue4'), alpha.f = 0.7), 3),
#                col_axis = rep(c('orange1', 'orange2', 'orange4'), each = 3), legend_x = 4, legend_y = 2.4, ylim = c(-1.5,1.5), ylab = 'IPC-NSC ATAC', xlim = c(3.5,9.5))
# # dev.off()


# # png('./output/sequence_modeling/figs/delta_atac_vs_homer.olig2_and_prox_atac_levels.png', h = 600, w = 300)
# par(las = 1, mar = c(12,2,8,2), cex.main = 1.72)
# boxplot_and_ks(feat_cut_df, factor1 = list('prox_atac_cut' = c(1:3)), factor2 = list('olig2_cut' = 2:3),ks.alternative = 'greater', y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '',
#                col_boxplot = rep(adjustcolor(c('lightgreen', 'green2', 'green4'), alpha.f = 0.7), 3),
#               col_axis = rep(c('orange1', 'orange2', 'orange4'), each = 3), legend_x = 4, legend_y = 2.4, ylim = c(-1.5,1.5), xlim = c(3.5,9.5))

# boxplot_and_ks(feat_cut_df, factor1 = list('prox_atac_cut' = c(1:3)), factor2 = list('nsc_atac_cut' = 1:3),ks.alternative = 'greater', y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '',
#               col_boxplot = rep(adjustcolor(c('lightgreen', 'green2', 'green4'), alpha.f = 0.7), 3),
#               col_axis = rep(c('red1', 'red3', 'red4'), each = 3), legend_y = 2.4, ylim = c(-1.5,1.5), xlim = c(0.5,9.65))


# dev.off()


