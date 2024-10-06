library(metacell)
library(misha)
library(misha.ext)
library(matrixStats)
library(ComplexHeatmap)
library(mcATAC)
library(tgstat)
library(plyr)
# devtools::load_all("/home/feshap/src/mcATAC")
devtools::load_all('~/src/iceqream/')

wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')
scdb_init(file.path(wd, 'scdb'), f=T)

mcmd <- readr::read_tsv(file.path(wd, 'BonevCollab/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]


options(gmax.data.size = 1e+9)
options(future.globals.maxSize = 8000 * 1024^2)

cust_st_ord <- c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))

cust_st_ord2 <- c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN/CfuPN','iCPN_early','iCPN_late',
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
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)
                             

c2mc_path <- file.path(wd, "output/mcatac/c2mc.tsv")
c2mc <- readr::read_tsv(c2mc_path)


mc_rna <- scdb_mc('pl_cort')

tss <- gintervals.load('intervs.global.tss')
tss <- tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]
tss2 <- tss[!duplicated(tss$geneSymbol),]
tads <- gintervals.load('intervs.global.tad_names')
nsc_mcs <- mcmd$metacell[mcmd$cell_type == 'NSC']
ipc_mcs <- mcmd$metacell[mcmd$cell_type %in%  c('IPC', 'IPC_cyc')]


# scc <- scc_read(path = file.path(wd, 'data/frag_reads_28122022/'))

# devtools::load_all("~/src/mcATAC/")

# scpeaks <- scc_to_peaks(sc_counts = scc, peaks = dplyr::select(mca@peaks, chrom, start, end, peak_name))

# scpeaks@peaks <- dplyr::mutate(scpeaks@peaks, peak_name_tb = peak_name, peak_name = peak_names(scpeaks@peaks[,1:3], tad_based = F))

# scpm <- scpeaks@mat

# rownames(scpm) <- scpeaks@peaks$peak_name

# astro_sc <- c2mc$cell[c2mc$metacell != 602 & c2mc$metacell %in% which(mcmd$cell_type == 'Astrocytes')]
# sc_602 <- c2mc$cell[c2mc$metacell == 602]

# astro_sc_mat <- scpm[,c(sc_602, astro_sc)]

# load(file.path(wd, 'output/mcatac/ct_peaks.rda'))

# ra <- ct_peaks[ct_peaks$type %in% c('astro_peak', 'neuro_peak'),] %>% select(peak_name, type)
# ra <- ra[!duplicated(ra$peak_name),]
# rownames(ra) <- ra$peak_name
# ra <- dplyr::select(ra, type)

# ac <- list(type = setNames(c('red', 'green'), c('neuro_peak', 'astro_peak')))

# table(ra$type)

# p_astro_peak_in_mc_602 <- pheatmap::pheatmap(as.matrix(astro_sc_mat[ct_peaks$peak_name[ct_peaks$type %in% c('astro_peak', 'neuro_peak')],]), cluster_cols = F, clustering_method = 'ward.D2', show_colnames = F, col = clrmp_abs,
#                    show_rownames = F,annotation_row = ra, annotation_colors = ac)
# save_pheatmap(p_astro_peak_in_mc_602, file.path(wd, 'output/mcatac/figs/test/astro_peaks_in_scatacs_from_mc_602.png'), h = 800, w = 600)


# ####################
# ### Don't remove
# ## How to create feat_peak mat object
# load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/multi_mmcortex.Rda",v=T)

# feat_peak = multi_model$atac_intervs
# scc <- scc_read('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads_28122022/')
# c2mc <- readr::read_tsv("./output/mcatac/c2mc.tsv")

# ## Move scATACs from fAMCs 602,603 (glia) to other fAMCs that AMCs 16_16, 16_20 flow to
# flow_res_path <- file.path(wd, "output/mcatac/pl_cort_flow_mat.tsv")
# flow_res_results <- readr::read_tsv(flow_res_path)
# frr_mat <- as.matrix(dplyr::select(flow_res_results, -rowname))
# rownames(frr_mat) <- flow_res_results$rowname
# glia_ct <- c('Astrocytes', 'Oligodendrocytes')
# mc_glia <- which(mcmd$cell_type %in% glia_ct)
# nnz_frcs <- apply(frr_mat[c('16_16', '16_20'),], 1, function(x) {sum_glia <- sum(x[mc_glia]); x[mc_glia] <- 0; nnz_mc <- which(x > 0); x_nnz_frc <- x[nnz_mc]/sum(x[nnz_mc]); return(setNames(x_nnz_frc, nnz_mc))})
# nnz_cumsum <- lapply(nnz_frcs, function(x) cumsum(sort(x)))
# flow_by_ct <- t(tgs_matrix_tapply(frr_mat, mcmd$cell_type, sum))
# flow_by_ct_norm <- flow_by_ct/rowSums(flow_by_ct)

# famc_new <- do.call('rbind', lapply(c(602, 603), function(x) {
#     scah <- c2mc$cell[c2mc$metacell == x]; 
#     amch <- rownames(frr_mat)[which(frr_mat[,x] > 0)];
#     sca_rn <- runif(n = length(scah))
#     famc_new <- as.numeric(sapply(seq_along(scah), function(i) mc_new <- names(nnz_cumsum[[amch]][nnz_cumsum[[amch]] >= sca_rn[[i]]])[[1]]))
#     # print(cbind(sca_rn, famc_new))
#     return(tibble::enframe(setNames(famc_new, scah), name = 'cell', value = 'metacell'))
# }))

# c2mc_new <- c2mc
# c2mc_new$metacell[match(famc_new$cell, c2mc_new$cell)] <- famc_new$metacell


# mcc <- scc_project_on_mc(sc_counts=scc, cell_to_metacell = dplyr::rename(c2mc_new, cell_id = cell))
# feat_peak$peak_name <- peak_names(feat_peak, tad_based = F)
# mcc_write(mcc, './output/mcatac/mmcortex_mcc_feat_peak_no_602_603', overwrite = T)



# mca <- mcc_to_mcatac(mc_counts=mcc, peaks=feat_peak, metadata=mcmd)
# saveRDS(mca, './output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds')

# ### Don't remove
# ####################


mca <- readRDS(file.path(wd, 'output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds'))
mca@rna_egc <- mc_rna@e_gc
mca@peaks$peak_name_ntb <- mcATAC::peak_names(subset(mca@peaks, select = -peak_name), tad_based = F)

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

# saveRDS(mca, file.path(wd, 'output/mcatac/mmcortex_mcatac_feat_peaks.rds'))

legc <- log2(1e-5 + mca@rna_egc)

a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))


# legc_tf <- legc[tfs_hi,]

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(mca@rna_egc))
head(tfs)

tfs_in <- intersect(tfs, rownames(legc)[which(rowMaxs(legc) >= -13.5)])

# tfs_hi <- tfs[(matrixStats::rowMaxs(legc[tfs,]) >= -13 & 
#               matrixStats::rowMaxs(legc[tfs,]) - matrixStats::rowMins(legc[tfs,]) >= 3) &
             # matrixStats::rowSds(legc[tfs,]) >= 0.6]
# tfs_hi <- tfs[which(log2(matrixStats::rowMaxs(mc_rna@mc_fp[tfs,])) - log2(matrixStats::rowMins(mc_rna@mc_fp[tfs,])) >= 3)]
# length(tfs_hi)


agg_id <- readr::read_csv(file.path(wd, 'scatac_data//aggregation_id.csv'))
mca@cell_to_metacell$agg_id <- as.numeric(unlist(purrr::map(stringr::str_split(mca@cell_to_metacell$cell_id, '-'), 2)))
agg_id_day <- data.frame(cbind(1:14, rep(12:18, each = 2)))
colnames(agg_id_day) <- c('agg_id', 'day')
mca@cell_to_metacell$day <- agg_id_day$day[mca@cell_to_metacell$agg_id]

atac_mc_day <- table(mca@cell_to_metacell$metacell, mca@cell_to_metacell$day)
atac_mc_day_norm <- atac_mc_day/rowSums(atac_mc_day)
rna_mc_day <- as.matrix(mcmd[,grep('E\\d\\d', colnames(mcmd))])
rna_mc_day_norm <- rna_mc_day/rowSums(rna_mc_day)
nsc_mcs <- which(mcmd$cell_type == 'NSC')
ipc_mcs <- which(mcmd$cell_type %in% c('IPC', 'IPC_cyc'))
egc_by_day <- mca@egc[,nsc_mcs] %*% atac_mc_day_norm[nsc_mcs,]
egc_by_day_n <- t(t(egc_by_day)/colSums(egc_by_day))
colnames(egc_by_day_n) <- paste0('E', colnames(egc_by_day_n))
egc_by_day_ipc <- mca@egc[,ipc_mcs] %*% atac_mc_day_norm[ipc_mcs,]
egc_by_day_n_ipc <- t(t(egc_by_day_ipc)/colSums(egc_by_day_ipc))
colnames(egc_by_day_n_ipc) <- paste0('E', colnames(egc_by_day_n_ipc))

a_legc_by_day_n <- log2(1e-5 + egc_by_day_n)

load(file.path(wd, 'output/metacell_model/cell_cycle_phase_data.rda'))

# load(file.path(wd, 'output/metacell_model/pcurve_nsc_vs_mat_neuro.rda'))
# ro <- rev(pcu$ord)
# # mcmd$cell_type[ro]

# load(file.path(wd, 'output/mcatac/mmcortex_feat_peak_variable_peak_clusters.rda'))
load(file.path(wd, 'output/mcatac/var_peaks_after_enh_prom_separation.rda'))

a_legc_avg_cl <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
a_legc_avg_ct <- t(tgs_matrix_tapply(a_legc[dist_peaks$peak_name,], mcmd$cell_type, mean))
a_legc_avg_cl_ct <- t(tgs_matrix_tapply(tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean), mcmd$cell_type, mean))

load(file=file.path(wd, 'output/mcatac/new_nsc_ipc_transition_peak_clusters_to_plot.rda'))

# options(repr.plot.width = 16)
# options(repr.plot.height = 16)

# scales <- setNames(unlist(sapply(lfits, function(x) x$estimate[['scale']])), names(lfits))

load(file.path(wd, 'output/methylation//avg_meth_all.rda'))

rm_ama <- rowMeans(avg_meth_all[,grep('E\\d\\d', colnames(avg_meth_all))], na.rm = T)

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

mat_rna <- scdb_mat('pl_cort')

source(file.path(wd, 'scripts/util.r'))

pca_nsc_ipc <- proximal_chromatin_activity(peaks_of_interest = mcp, background_peaks = mcp, 
                                           d_puncture = 1e+3, d_proximity_atac = sum(.misha$ALLGENOME[[1]]$end)/nrow(a_legc), d_proximity_rna = 5e+5, 
                                           mc_sel = nsc_mcs, mat_rna = mat_rna, mc_rna = mc_rna, tss = tss, tads = tads, mca = mca)

rna_umi_vec <- log2(1+tibble::deframe(pca_nsc_ipc[,c('peak_name', 'rna_umi')]))

atac_umi_vec <- log2(1+tibble::deframe(pca_nsc_ipc[,c('peak_name', 'atac_umi')]))

load(file.path(wd, 'output/sequence_modeling/nuc_feat_peaks.rda'))

scp <- sort(cln_plot)

load(file.path(wd, 'output/sequence_modeling/dinucs_feat_peaks.rda'))

# n_5mCpG <- 300*dinucs_per_peak[,'CG']*rm_ama[rownames(dinucs_per_peak)]
# n_nmCpG <- 300*dinucs_per_peak[,'CG']*(1-rm_ama[rownames(dinucs_per_peak)])

# load(file.path(wd, 'output/sequence_modeling/enhflow_output_delta_ipc_nsc_n_motifs=8.rda'))

# load(file.path(wd, 'output/sequence_modeling/enhflow_model_w_add_feat_120923.rda'))

# load(file.path(wd, 'output/sequence_modeling/enhflow_model_50_motifs_w_hc.rda'))
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

prego_motifs <- do.call('rbind', lapply(seq_along(tm_w_add_feat@motif_models), function(i) dplyr::mutate(tm_w_add_feat@motif_models[[i]]$pssm, motif = names(tm_w_add_feat@motif_models)[[i]])))
prego_ie <- prego::gextract_pwm(intervals = mcp, motifs = unique(prego_motifs$motif), dataset = prego_motifs)
prego_ie_mat <- subset(prego_ie, select = -c(chrom, start, end, peak_name))
rownames(prego_ie_mat) <- prego_ie$peak_name
prego_motifs_renamed <- rename_prego_motifs(tm_w_add_feat, tfs_in)
colnames(prego_ie_mat) <- prego_motifs_renamed[colnames(prego_ie_mat)]
motif_clustering <- tm_w_add_feat@params$distilled_features

amd <- prego::all_motif_datasets()

amou <- unique(amd[,c('motif', 'motif_orig')])

feats <- names(scdb_gset('pl_cort_dns')@gene_set)

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

km_prom_a_legc <- tglkmeans::TGL_kmeans(a_legc[prom_peaks$peak_name,], k = 10, seed = 1337)

km_enh_a_legc <- tglkmeans::TGL_kmeans(a_legc[dist_peaks$peak_name,], k = 60, seed = 1337)

a_legc_avg_cl_enh <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
cl_cvs <- rowSds(a_legc_avg_cl_enh)**2/(rowMeans(a_legc_avg_cl_enh) + 17.6)**2

m1 <- a_legc_avg_cl_enh[head(names(sort(cl_cvs)),13),cust_mc_ord_st2]

m1_hc <- hclust(dist(m1), method = 'ward.D2')
m2 <- a_legc_avg_cl_enh[tail(names(sort(cl_cvs)),-13),cust_mc_ord_st2]
m2_hc <- hclust(dist(m2), method = 'ward.D2')
m2_ord <- order(apply(m2, 1, function(x) sum(x*1:length(x))/sum(x)))

m2_hcct <- cutree(m2_hc, k = 16)

m2_hcct_ord <- order(apply(tgs_matrix_tapply(t(m2), m2_hcct, mean), 1, function(x) sum(x*1:length(x))/sum(x)))
m2_hcct_ord
# m2_hcct[m2_hc$order]

m2_cl_ord <- do.call('c', sapply(m2_hcct_ord, function(u) names(m2_hcct[m2_hcct == u])))

enh_cl_ord <- c(rownames(m1)[m1_hc$order], rev(m2_cl_ord))

var_peaks <- dist_peaks$peak_name[km_enh_a_legc$cluster %in% as.numeric(rownames(m2))]

# save(prom_peaks, dist_peaks, enh_cl_ord, km_enh_a_legc, km_prom_a_legc, m1, m2, var_peaks, file = file.path(wd, 'output/mcatac/var_peaks_after_enh_prom_separation.rda'))

intervs_energy_new <- readRDS(file.path(wd, 'output/sequence_modeling/mmcortex_feat_peak_prego_motif_energy.rds'))
ie_mat <- as.matrix(subset(intervs_energy_new, select = -c(chrom, start, end, peak_name, mmcortex.marginal,	intervalID,	peak_name, peak_name_ntb)))
rownames(ie_mat) <- intervs_energy_new$peak_name

raq98 <- unlist(plyr::llply(colnames(ie_mat), function(x) quantile(ie_mat[,x], probs = 0.98), .parallel = T))

names(raq98) <- colnames(ie_mat)

ra_98_bin_int <- t(plyr::laply(1:length(raq98), function(i) as.numeric(ie_mat[,i] >= raq98[[i]]), .parallel = T))
# save(ra_98_bin_int, file=file.path(wd, 'output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda'))
colnames(ra_98_bin_int) <- colnames(ie_mat)
rownames(ra_98_bin_int) <- rownames(ie_mat)

load(file.path(wd, 'output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda'))
ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int[dist_peaks$peak_name,])), km_enh_a_legc$cluster, sum)
ra_98_lfc_int <- log2(1e-2+ra_98_sum_clust_int/(0.02*as.numeric(km_enh_a_legc$size)))
colnames(ra_98_lfc_int) <- colnames(ra_98_bin_int)

mtt <- c("JASPAR.EOMES", 
         "JOLMA.MEIS2_mono_DBD_2",
         'JASPAR.NEUROG1',
         "JASPAR.NEUROD1",
         'HOMER.Sox2',
         "JASPAR.EMX1",
         'JASPAR.POU3F2',
         'JASPAR.NFIA',
         "JASPAR.MEF2C",
         "JASPAR.FOXP1",
         "HOCOMOCO.MECP2_MOUSE.H11MO.0.C",
         'HOMER.CTCF',
         'HOMER.NRF1',
         'HOCOMOCO.KLF3_MOUSE.H11MO.0.A',
         'JOLMA.ETV1_mono_DBD',
         'JASPAR.NFIB'
        )

options(repr.plot.width =10)
options(repr.plot.height =16)

p_motifs_lfc_mtt <- pheatmap::pheatmap(round(ra_98_lfc_int[enh_cl_ord,mtt], 3), breaks = seq(-3,3, l = length(clrmp_rel)), col = clrmp_rel, cluster_rows = F, fontsize = 12, treeheight_col = 0)
save_pheatmap(p_motifs_lfc_mtt, file.path(wd, 'output/mcatac/figs/test/motifs_annot_mat_mtt.png'), h = 1800, w = 850)

# New version - 25/04/24
inds_glia <- which(names(cust_mc_ord_st2) %in% c('Oligodendrocytes', 'Astrocytes'))
inds_no_glia <- which(!(names(cust_mc_ord_st2) %in% c('Oligodendrocytes', 'Astrocytes')))


pltmt <- a_legc_avg_cl_enh[enh_cl_ord, cust_mc_ord_st2]

brks <- seq(-16.7,-14.75,l=100)

col_ha1 <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_no_glia],'cell_type'], 
                                                     col = ann_colors[['cell_type']],
                                                     height =unit(2, 'cm')), 
                             mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_no_glia],'mean_day'], axis_param = list(gp = gpar(fontsize = 0)),
                                                    # col = circlize::colorRamp2(seq(13,18,1), c('red', 'orange', 'yellow', 'green2', 'blue', 'purple')), 
                                                    height =unit(2, 'cm')), 
                             annotation_name_gp = gpar(fontsize = 18),
                             show_legend = F)
# row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_a_legc$size[as.numeric(hcct2)]), ylim = c(800,2500), gp = gpar(fill = 'black',fontsize = 18), 
row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_enh_a_legc$size[as.numeric(enh_cl_ord)]), ylim = c(800,8000), gp = gpar(fill = 'black',fontsize = 18), 
                                                        axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90)), 
                            # frac_prom = anno_barplot(as.numeric(t(perc_prom[peak_clust_var[hcct2]])), ylim = c(0,1), 
                            # gp = gpar(fill = 'black',fontsize = 18), 
                            # axis_param = list(facing = 'inside', gp = gpar(fontsize = 20), labels_rot = -90),
                            annotation_name_gp = gpar(fontsize = 18),
                            annotation_name_rot = 0,
                            # cl2 = anno_simple(x = ct_hc_a_legc[ord1], col = setNames(peak_cl2_color_key$color[match(ct_hc_a_legc[ord1], peak_cl2_color_key$cl2)], ct_hc_a_legc[ord1])),
                            annotation_name_offset = unit(3, 'cm'),        
                            which = 'row',
                            width = unit(5, 'cm')
                           )
clrmp_rel2 <- circlize::colorRamp2(breaks = c(-3,0,3), colors = c('blue3','white','red3'))
motifs_anno_mat <- round(ra_98_lfc_int[enh_cl_ord,mtt], 3)
# motifs_anno_mat <- round(ra_98_lfc_int[hcct2,mtt], 3)
# motifs_anno_mat <- round(prego_98_lfc[hcct2,], 3)
hc_mtt <- hclust(dist(t(motifs_anno_mat)), method = 'ward.D2')
colnames(motifs_anno_mat) <- unlist(purrr::map(stringr::str_split(colnames(motifs_anno_mat), '\\.'), 2))
mtt_ord <- colnames(motifs_anno_mat)[hc_mtt$order]
mam_lin <- apply(motifs_anno_mat, 2, function(x) {y <- x; y[x < -3] <- -3; y[x > 3] <- 3; return((y + 3)/6)})


motif_ha <- HeatmapAnnotation(SOX = anno_numeric(motifs_anno_mat[,grep('Sox2', colnames(motifs_anno_mat), ign = T)], 
                                                 labels_gp = gpar('fontsize' = 0),
                                                    bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('sox2', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
                              NEUROD1 = anno_numeric(motifs_anno_mat[,grep('NEUROD1', colnames(motifs_anno_mat), ign = T)], 
                                                     labels_gp = gpar('fontsize' = 0),
                                                    bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('neurod1', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
                              EOMES = anno_numeric(motifs_anno_mat[,grep('eomes', colnames(motifs_anno_mat), ign = T)], 
                                                    labels_gp = gpar('fontsize' = 0),
                                                    bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('eomes', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
                              POU3F2 = anno_numeric(motifs_anno_mat[,grep('pou3f2', colnames(motifs_anno_mat), ign = T)], 
                                                    labels_gp = gpar('fontsize' = 0),
                                                    bg_gp = gpar('fill' = matrix(adjustcolor(rep(clrmp_rel[1+round(999*mam_lin[,grep('pou3f2', colnames(mam_lin), ign = T)])], 2), alpha.f = 0.75), ncol = 2))),
                              which = 'row',
                            width = unit(9, 'cm')
                             )
ch <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_no_glia], name = 'log2\nfraction\nATAC',
                              # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
                              col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
                              column_split = factor(names(cust_mc_ord_st2[inds_no_glia]), levels = cust_st_ord2), column_gap = unit(2, 'mm'),
                              column_title_gp = gpar(fontsize = 20),
                              column_title_rot = 90,
                              row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), row_gap = unit(2, 'mm'),
                                                            row_title_gp = gpar(fontsize = 0),
                              top_annotation = col_ha1, 
                              # bottom_annotation = col_ha2,
                              show_heatmap_legend = T,
                              show_column_names = F,
                              right_annotation = motif_ha,
                              heatmap_legend_param = list(legend_height = unit(5, 'in'), legend_width = unit(5, 'in'), labels_gp = gpar(fontsize = 16)),
                                row_names_gp = gpar(fontsize = 14),
                              # heatmap_width = unit(92, 'cm'), heatmap_height = unit(50, 'cm'),
                              heatmap_width = unit(45, 'cm'), heatmap_height = unit(32, 'cm'),
                        # left_annotation = row_ha,
                        cluster_columns = F, cluster_rows = F)

col_ha1_glia <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_glia],'cell_type'], 
                                                     col = ann_colors[['cell_type']],
                                                     height =unit(2.2, 'cm')), 
                             mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_glia],'mean_day'], axis_param = list(gp = gpar(fontsize = 16)), ylim = c(13,18),
                                                    height =unit(2.2, 'cm')), 
                             annotation_name_gp = gpar(fontsize = 0),
                             show_legend = F)

# ch_glia <- ComplexHeatmap::Heatmap(matrix = pltmt[hcct2,inds_glia], 
ch_glia <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_glia], 
                                # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
                              col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
                                   column_split = factor(names(cust_mc_ord_st2[inds_glia]), levels = c('Oligodendrocytes', 'Astrocytes')), column_gap = unit(3, 'mm'),
                              column_title_gp = gpar(fontsize = 20),
                              column_title_rot = 90,
                              row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), row_gap = unit(2, 'mm'),
                              row_title_gp = gpar(fontsize = 0),
                              top_annotation = col_ha1_glia, 
                                   left_annotation = row_ha,
                              show_heatmap_legend = F,
                              show_column_names = F,
                                   show_row_names = F,
                              heatmap_width = unit(10, 'cm'), heatmap_height = unit(36, 'cm'),
                              row_names_gp = gpar(fontsize = 14),
                        cluster_columns = F, cluster_rows = F)

options(repr.plot.width = 14)
options(repr.plot.height = 14)

# draw(ch)

# png(file.path(wd, 'output/mcatac/figs/test/mmcortex_famc_legc_cluster_size_frac_prom_annot_arb_met_0.45.png'), w = 3200, h =1600)
# png(file.path(wd, 'output/mcatac/figs/test/mmcortex_famc_legc_cluster_size_frac_prom_annot_left_arb_met_0.45.png'), w = 2000, h =1600)
# png(file.path(wd, 'output/mcatac/figs/test/mmcortex_famc_legc_w_motif_q98_lfc.png'), w = 3000, h =1700)
png(file.path(wd, 'output/mcatac/figs/test/mmcortex_famc_legc_w_motif_q98_lfc_no_glia.png'), w = 1380, h =1170)
draw(ch)
dev.off()

load(file = file.path(wd, 'output/methylation/avg_meth_all.rda'))

nsc_inc_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] >= 1)], rownames(avg_meth_all))
nsc_inc_atac <- a_legc_by_day_n[nsc_inc_peaks,]
nsc_dec_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] <= -1)], rownames(avg_meth_all))
nsc_dec_atac <- a_legc_by_day_n[nsc_dec_peaks,]

par(mfrow = c(1,2))
# ord <- order(avg_meth_all[nsc_inc_peaks,'E17'] - avg_meth_all[nsc_inc_peaks,'E13'])
ord <- order(rowMeans(a_legc_by_day_n[nsc_inc_peaks,]))
p_inc_atac <- pheatmap::pheatmap(a_legc_by_day_n[nsc_inc_peaks[ord],], fontsize_col = 16,
                                      main = 'NSC increasing peaks',
                                      fontsize = 10,
                                      cluster_rows = F, 
                                      clustering_method = 'ward.D2',
                                      cluster_cols = F,
                                      show_rownames = F,
                                      treeheight_row = 0,
                                      color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
                                     breaks = seq(-16.6,-14,l = 1000))
save_pheatmap(p_inc_atac, file.path(wd, 'output/mcatac/figs/test/phm_nsc_inc_atac.png'), h = 600, w = 500)
# ord <- order(avg_meth_all[nsc_dec_peaks,'E17'] - avg_meth_all[nsc_dec_peaks,'E13'], decreasing = T)
ord <- order(rowMeans(a_legc_by_day_n[nsc_dec_peaks,]))
p_dec_atac <- pheatmap::pheatmap(a_legc_by_day_n[nsc_dec_peaks[ord],], fontsize_col = 16,
                                      main = 'NSC decreasing peaks',
                                      fontsize = 10,
                                      cluster_rows = F, 
                                      clustering_method = 'ward.D2',
                                      cluster_cols = F,
                                      show_rownames = F,
                                      treeheight_row = 0,
                                      color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
                                     breaks = seq(-16.6,-14,l = 1000))
save_pheatmap(p_dec_atac, file.path(wd, 'output/mcatac/figs/test/phm_nsc_dec_atac.png'), h = 600, w = 500)

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

mpra_lib <- readr::read_tsv(file.path(wd, 'cpg_methylation//st_and_temporal_and_e14_shadow_enh_9-1-22.tsv'))
mpra_lib <- as.data.frame(dplyr::relocate(mpra_lib, rowname, .after = end))
ct_seq <- unlist(purrr::map(stringr::str_split(mpra_lib$rowname, '\\.'), 1))
ct_seq <- ifelse(ct_seq  %in%  c('IPC', 'NSC', 'Mature'), ct_seq, 'shadow')

dir_seq <- as.character(purrr::map(stringr::str_split(mpra_lib$rowname, '\\.'), 2))

mpra_lib$ct <- ct_seq
mpra_lib$dir <- dir_seq

nei_mpra_mcp <- tidyr::drop_na(gintervals.neighbors1(mpra_lib, mcp, maxdist = 0, mindist = 0, maxneighbors = 3))

nsc_asc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_asc')
nsc_asc_atac <- a_legc_by_day_n[nsc_asc_nei$peak_name,]

nsc_desc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_desc')
nsc_desc_atac <- a_legc_by_day_n[nsc_desc_nei$peak_name,]


vasc <- colMeans(a_legc[nsc_asc_nei$peak_name,])
vdesc <- colMeans(a_legc[nsc_desc_nei$peak_name,])

legc <- log2(1e-5 + mc_rna@e_gc)

pltmt <- rbind(vasc,
               vdesc,
               legc[sort(grep('os|Sox10|Sox13|Sox7|Sox17|Sox15|Sox18|Sox30',
                         # c(
                             grep('Tet|Dnmt',rownames(mc_rna@e_gc), v=T) 
                           # sort(grep('Eomes|Pou3f2|Neurog2|Sox|Hox', rownames(mc_rna@e_gc), v=T)))
                        , inv=T, v=T)),mcmd$metacell]
               )

rownames(pltmt) <- c('inc_peaks', 'dec_peaks', tail(rownames(pltmt), -2))

par(mfrow = c(2,3))
cells_h <- as.character(mcmd$metacell[which(mcmd$cell_type == 'NSC')])
lm_and_plot <- function(x, y, xlab = NULL, ylab = NULL, main = NULL) {
    lm1 <- lm(y ~ x)
    plot(x,y, xlab = xlab, ylab = ylab, main = glue::glue("cor = {round(cor(x,y,method = 'pearson'), 3)}, pval = {round(cor.test(x,y,method ='pearson')$p.value, 3)}"))
    abline(lm1$coefficients[[1]], lm1$coefficients[[2]], col= 'red', lty =2)
}

png(file.path(wd, 'output/methylation/figs/test/scatter_cor_nsc_inc_dec_peaks_vs_meth_components.png'), h = 600, w = 1000, res = 100)
par(mfrow = c(2,3), cex.lab = 1.5, mar = c(5,5,2,1), cex.axis = 1.5, cex.main = 1.5)
lm_and_plot(pltmt['Dnmt1',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt1', ylab = 'NSC dec_peaks ATAC')
lm_and_plot(pltmt['Dnmt3a',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3a', ylab = 'NSC dec_peaks ATAC')
lm_and_plot(pltmt['Dnmt3b',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3b', ylab = 'NSC dec_peaks ATAC')
lm_and_plot(pltmt['Tet1',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet1', ylab = 'NSC inc_peaks ATAC')
lm_and_plot(pltmt['Tet2',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet2', ylab = 'NSC inc_peaks ATAC')
lm_and_plot(pltmt['Tet3',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet3', ylab = 'NSC inc_peaks ATAC')
dev.off()