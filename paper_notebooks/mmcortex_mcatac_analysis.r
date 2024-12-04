library(metacell)
library(misha)
library(misha.ext)
library(matrixStats)
library(ComplexHeatmap)
# library(mcATAC)
library(tgstat)
library(plyr)
devtools::load_all("/home/feshap/src/mcATAC")
devtools::load_all('~/src/iceqream/')

wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')
scdb_init(file.path(wd, 'scdb'), f=T)

source(file.path(wd, 'scripts/util.r'))


mcmd <- readr::read_tsv(file.path(wd, 'BonevCollab/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]

options(gmax.data.size = 1e+9)
options(future.globals.maxSize = 8000 * 1024^2)

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
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)

mc_rna <- scdb_mc('pl_cort')
mat_rna <- scdb_mat('pl_cort')
marks <- names(scdb_gset('pl_cort_marks_f')@gene_set)

tss <- gintervals.load('intervs.global.tss')
tss <- tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]
tss2 <- tss[!duplicated(tss$geneSymbol),]
tads <- gintervals.load('intervs.global.tad_names')
nsc_mcs <- mcmd$metacell[mcmd$cell_type == 'NSC']
ipc_mcs <- mcmd$metacell[mcmd$cell_type %in%  c('IPC', 'IPC_cyc')]

mca <- readRDS(file.path(wd, 'output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds'))
mca@rna_egc <- mc_rna@e_gc
mca@peaks$peak_name_ntb <- mcATAC::peak_names(subset(mca@peaks, select = -peak_name), tad_based = F)

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

legc <- log2(1e-5 + mca@rna_egc)

a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(mca@rna_egc))
head(tfs)

tfs_in <- intersect(tfs, rownames(legc)[which(rowMaxs(legc) >= -13.5)])

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

a_legc_by_day_n <- log2(1e-5 + egc_by_day_n)

load(file.path(wd, 'output/metacell_model/cell_cycle_phase_data.rda'))

load(file.path(wd, 'output/metacell_model/pcurve_nsc_vs_mat_neuro.rda'))
ro <- rev(pcu$ord)

load(file.path(wd, 'output/methylation//avg_meth_all.rda'))

load(file.path(wd, 'output/sequence_modeling/enhflow_model_n_motif=16_w_clustering.rda'))

load('./output/metacell_model/diff_order_data.rda')

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

km_prom_a_legc <- tglkmeans::TGL_kmeans(a_legc[prom_peaks$peak_name,], k = 10, seed = 1337)

km_enh_a_legc <- tglkmeans::TGL_kmeans(a_legc[dist_peaks$peak_name,], k = 60, seed = 1337)

a_legc_avg_cl_enh <- tgs_matrix_tapply(t(a_legc[dist_peaks$peak_name,]), km_enh_a_legc$cluster, mean)
a_legc_avg_cl_prom <- tgs_matrix_tapply(t(a_legc[prom_peaks$peak_name,]), km_prom_a_legc$cluster, mean)
a_legc_avg_cl_ct <- t(tgs_matrix_tapply(a_legc_avg_cl_enh, mcmd$cell_type, mean))

cl_cvs <- rowSds(a_legc_avg_cl_enh)**2/(rowMeans(a_legc_avg_cl_enh) + 17.6)**2

m1 <- a_legc_avg_cl_enh[head(names(sort(cl_cvs)),13),cust_mc_ord_st2]

m1_hc <- hclust(dist(m1), method = 'ward.D2')
m2 <- a_legc_avg_cl_enh[tail(names(sort(cl_cvs)),-13),cust_mc_ord_st2]
m2_hc <- hclust(dist(m2), method = 'ward.D2')
m2_ord <- order(apply(m2, 1, function(x) sum(x*1:length(x))/sum(x)))

m2_hcct <- cutree(m2_hc, k = 16)

m2_hcct_ord <- order(apply(tgs_matrix_tapply(t(m2), m2_hcct, mean), 1, function(x) sum(x*1:length(x))/sum(x)))

m2_cl_ord <- do.call('c', sapply(m2_hcct_ord, function(u) names(m2_hcct[m2_hcct == u])))

enh_cl_ord <- c(rownames(m1)[m1_hc$order], rev(m2_cl_ord))

var_peaks <- dist_peaks$peak_name[km_enh_a_legc$cluster %in% as.numeric(rownames(m2))]

# save(prom_peaks, 
#     dist_peaks, 
#     enh_cl_ord, 
#     km_enh_a_legc, 
#     km_prom_a_legc, 
#     m1, 
#     m2, 
#     var_peaks, 
#             file = file.path(wd, 'output/mcatac/var_peaks_after_enh_prom_separation.rda'))

intervs_energy_new <- readRDS(file.path(wd, 'output/sequence_modeling/mmcortex_feat_peak_prego_motif_energy.rds'))
ie_mat <- as.matrix(subset(intervs_energy_new, 
            select = -c(chrom, start, end, peak_name, mmcortex.marginal, intervalID, peak_name, peak_name_ntb)))
rownames(ie_mat) <- intervs_energy_new$peak_name

raq98 <- unlist(plyr::llply(colnames(ie_mat), function(x) quantile(ie_mat[,x], probs = 0.98), .parallel = T))
names(raq98) <- colnames(ie_mat)

ra_98_bin_int <- t(plyr::laply(1:length(raq98), function(i) as.numeric(ie_mat[,i] >= raq98[[i]]), .parallel = T))
colnames(ra_98_bin_int) <- colnames(ie_mat)
rownames(ra_98_bin_int) <- rownames(ie_mat)

load(file.path(wd, 'output/sequence_modeling/feat_peak_all_motif_q_98_binary_matrix.rda'))
ra_98_sum_clust_int <- tgs_matrix_tapply(t(as.matrix(ra_98_bin_int[dist_peaks$peak_name,])), km_enh_a_legc$cluster, sum)
ra_98_lfc_int <- log2(1e-2+ra_98_sum_clust_int/(0.02*as.numeric(km_enh_a_legc$size)))
colnames(ra_98_lfc_int) <- colnames(ra_98_bin_int)

a_legc_avg_ct <- t(tgs_matrix_tapply(a_legc, mcmd$cell_type, mean, na.rm = T))

cond001 <- rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) > 1 & rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - a_legc_avg_ct[,'IPC'] > 1
cond010 <- rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) < -1 & rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) - a_legc_avg_ct[,'IPC'] > 1
cond110 <- rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) < -1 & abs(rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) - a_legc_avg_ct[,'IPC']) < 1
cond101 <- rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) > 1 & abs(rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - a_legc_avg_ct[,'IPC']) < 1
cond011 <- abs(rowMeans(a_legc_avg_ct[,c('CthPN', 'SCPN')]) - rowMeans(a_legc_avg_ct[,c('CPN_L2-3', 'CPN_L5_6')]) < 1)

pks011 <- rownames(a_legc_avg_ct[cond011,])
pks001 <- setdiff(rownames(a_legc_avg_ct[cond001,]), pks011)
pks010 <- setdiff(rownames(a_legc_avg_ct[cond010,]), pks011)
pks110 <- rownames(a_legc_avg_ct[cond110,])
pks101 <- rownames(a_legc_avg_ct[cond101,])

library(mcATAC)

mcacp <- mca
mcacp <- normalize_egc(mcacp, 'mmcortex.marginal')

mcacp@metadata <- subset(mcacp@metadata, subset = !(1:nrow(mcacp@metadata) %in% 602:603))

mcacp <- add_const_peaks(mcacp, -16)

mcacp <- normalize_to_prob(mcacp)



bins_nsc_neu_f <- bins_nsc_neu[as.numeric(mcmd$metacell)]

save(ra_98_lfc_int, 
    enh_cl_ord,
    m1,
    km_enh_a_legc,
    km_prom_a_legc,
    a_legc_avg_ct,
    a_legc_avg_cl_prom,
    a_legc_avg_cl_enh,
    pks001,
    pks010,
    pks101,
    pks110,
    mcacp,
    file = file.path(wd, 'output/mcatac/fig3_data.rda'))

nsc_inc_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] >= 1)], rownames(avg_meth_all))
nsc_inc_atac <- a_legc_by_day_n[nsc_inc_peaks,]
nsc_dec_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] <= -1)], rownames(avg_meth_all))
nsc_dec_atac <- a_legc_by_day_n[nsc_dec_peaks,]


astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')
ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')
stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')
Q_THRESH <- 0
prom_sds <- setNames(rowSds(a_legc[prom_peaks$peak_name,]), prom_peaks$peak_name)
prom_mxmns <- rowMaxmins(a_legc[prom_peaks$peak_name,])
proms_hi_var_astro <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% astro_module]
proms_hi_var_ipc <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% ipc_module]
proms_hi_var_stem <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% stem_module]
proms_hi_var_marks <- prom_peaks$peak_name[prom_mxmns >= quantile(prom_mxmns, 0) & prom_sds >= quantile(prom_sds, Q_THRESH) & prom_peaks$geneSymbol %in% marks]

save(prom_sds, 
        proms_hi_var_astro, 
        proms_hi_var_ipc, 
        proms_hi_var_stem, 
        proms_hi_var_marks,
        file = file.path(wd, 'output/mcatac/fig_s3_data.rda'))


save(nsc_inc_peaks, 
        nsc_dec_peaks, 
        a_legc_by_day_n,
        file = file.path(wd, 'output/mcatac/fig4_atac_data.rda'))