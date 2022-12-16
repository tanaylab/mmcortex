### Script to compare clusters obtained from scATAC-derived microclusters 
## and microclusters derived from scATAC-RNA metacell correlations
wd <- "/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/"
setwd(wd)
library(matrixStats)
library(pheatmap)
library(metacell)
scdb_init('scdb')
mat_prom <- scdb_mat('pl_prom_cort')
options(gmax.data.size = 1e+9)
devtools::load_all("~/src/mcATAC/")
library(prego)
gset_genome('mm10')
doMC::registerDoMC(cores = 70)
mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
colnames(mcmd)[colnames(mcmd) == 'st'] <- 'cell_type'
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))                         
col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))


aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
cl_raw <- aaa$km_a_legc
cl_rna <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
mcl_all <- aaa$mcl_all


mat_tbl <- as.matrix(table(cl_raw$cluster, cl_rna$cluster))
colnames(mat_tbl) <- paste0('raw_', colnames(mat_tbl))
rownames(mat_tbl) <- paste0('rna_', rownames(mat_tbl))
mat_tbl_norm_col <- t(t(mat_tbl)/colSums(mat_tbl))


ord <- order(apply(mat_tbl, 1, sum))
ord_col <- order(colSums(mat_tbl))
ord_row_norm <- order(rowSums(mat_tbl_norm_col))
ord_col_norm <- order(colSums(mat_tbl_norm_col))
ord_row_max <- order(apply(mat_tbl, 1, max))
ord_row_max_norm <- order(apply(mat_tbl_norm_col, 1, max))
ord_col_max <- order(apply(mat_tbl, 2, max))


# p_raw <- pheatmap(mat_tbl, color = colorRampPalette(c('white','red', 'black'))(100))
# p_norm <- pheatmap(mat_tbl_norm_col, color = colorRampPalette(c('white','red', 'black'))(100))

sc_cor_kms <- readRDS('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/microcluster_assignment.RDS')

day_mcl_prom <- t(do.call('rbind', 
                    lapply(1:length(sc_cor_kms), 
                            function(i) tgs_matrix_tapply(as.matrix(mat_prom@mat[,colnames(aaa$sc_res[[i]]$sc_cor)]), sc_cor_kms[[i]]$cluster, sum)
                            )))

nm <- "pl_prom_cort"

gstat = scdb_gstat(nm)

x = log(gstat$ds_mean)
init_filt = which(x >= -4)
y = gstat$ds_log_varmean[init_filt]
name_filt = gstat$name
x = x[init_filt]
xcut = cut(x, breaks = seq(min(x), max(x), l = 50))
top_q_inds = lapply(levels(xcut), function(l) {inds = which(xcut == l); 
                                  xfilt = y[inds]; 
                                  xtop = head(inds[order(xfilt, decreasing = T)], 30); 
                                  return(xtop)
                                 }
      )

names(top_q_inds) = levels(xcut)
mc_rna <- scdb_mc('pl_cort')
feats = intersect(rownames(gstat)[init_filt[unlist(top_q_inds)]], rownames(mc_rna@mc_fp))
cor_mcl_mc <- tgs_cor(day_mcl_prom[feats,], mc_rna@mc_fp[feats,], spearman=T)
rownames(cor_mcl_mc) <- sapply(13:18, function(x) paste0(x, '_', 1:30))

order_columns_hc_com <- function(mat, hc_cut_deg) {
    hcm <- hclust(dist(t(mat)), method = 'ward.D2')
    hc_ct <- cutree(hcm, hc_cut_deg)
    com_hc <- sapply(unique(hc_ct), function(ci) {
                        vec <- rowMeans(mat[,hc_ct == ci]); 
                        return(sum(vec*1:length(vec))/sum(vec))
                        })
    ord <- as.numeric(do.call('c', lapply(unique(hc_ct)[order(com_hc)], function(ci) which(hc_ct == ci))))
    return(ord)
}

amc_st_assn <- mcmd$cell_type[apply(cor_mcl_mc, 1, which.max)]
amc_cat_assn <- ifelse(amc_st_assn %in% c('NSC', 'IPC', 'IPC_cyc', 'Astrocytes', 'Oligodendrocytes'), amc_st_assn, 'Neuron')
amc_cat_assn <- ifelse(amc_st_assn %in% c('Astrocytes', 'Oligodendrocytes'), 'NSC', amc_cat_assn)
amc_ord <- do.call('c', lapply(c('NSC' ,'IPC', 'Neuron'), function(ct) {
    vec <- which(amc_cat_assn == ct);
    com_ord <- order_columns_hc_com(t(cor_mcl_mc[vec,]), 5);
    return(setNames(vec[com_ord], rep(ct, length(vec))))
    }))

row_annot <- as.data.frame(amc_cat_assn)
rownames(row_annot) <- colnames(mcl_all)
colnames(row_annot) <- 'cell_type'
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))
ann_colors[['cell_type']] <- c(ann_colors[['cell_type']], setNames('yellow', 'Neuron'))

save(row_annot, amc_ord, ann_colors, file = './output/mcatac/atac_metacell_order_annotation.rda')

pltmt <- cor_mcl_mc
rownames(pltmt) <- colnames(mcl_all)

# p_cor_mcl_mc_order <- pheatmap(pltmt[amc_ord,cust_mc_ord_st], 
#                                     color = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#                                     breaks = seq(-.3,.3, l = 101),
#                                     cluster_cols = F, cluster_rows = F,
#                                      annotation_col = col_annot, annotation_row = row_annot,
#                                      annotation_colors = ann_colors, 
#                                      fontsize_row = 4, fontsize_col = 4)
# save_pheatmap(p_cor_mcl_mc_order, './output/mcatac/figs/cor_amc_rmc.png', height=2500, width = 3500, res = 300)


mcl_all_norm <- t(t(mcl_all)/Matrix::colSums(mcl_all))
mcl_all_sum_atac <- t(log2(1e-5 + tgs_matrix_tapply(t(mcl_all_norm), cl_raw$cluster, mean)))
mcl_all_sum_rna <- t(log2(1e-5 + tgs_matrix_tapply(t(mcl_all_norm), cl_rna$cluster, mean)))
var_thresh_up <- -15
var_thresh_dn <- -16
varp_atac <- which(colMaxs(mcl_all_sum_atac) > var_thresh_up & 
                    colMins(mcl_all_sum_atac) < var_thresh_dn)
varp_rna <- which(colMaxs(mcl_all_sum_rna) > var_thresh_up & 
                            colMins(mcl_all_sum_rna) < var_thresh_dn)

pltmt_atac <- mcl_all_sum_atac[amc_ord,varp_atac]
pltmt_rna <- mcl_all_sum_rna[amc_ord,varp_rna]
ord_col_atac <- order_columns_hc_com(2**pltmt_atac, 4)
ord_col_rna <- order_columns_hc_com(2**pltmt_rna, 4)

# p_sum_atac <- pheatmap(pltmt_atac[,ord_col_atac], cluster_cols =F, cluster_rows = F, fontsize_row = 5, fontsize_col = 8, main = 'log2(UMI) - ATAC clusters',
#                         annotation_row = row_annot,annotation_colors = ann_colors)
# p_sum_rna <- pheatmap(pltmt_rna[,ord_col_rna], cluster_cols = F, cluster_rows = F, fontsize_row = 8, fontsize_col = 5, main = 'log2(UMI) - RNA clusters',
#                         annotation_row = row_annot,annotation_colors = ann_colors)
# save_pheatmap(p_sum_atac, './output/mcatac/figs/peak_clustering_sum_atac.png', height=2600, width = 2600, res = 300)
# save_pheatmap(p_sum_rna, './output/mcatac/figs/peak_clustering_sum_rna.png', height=2600, width = 2600, res = 300)

sp_varp_atac <- do.call('c', lapply(varp_atac[ord_col_atac], function(pcl) {vec <- rownames(mcl_all)[which(cl_raw$cluster == pcl)]; 
                                    setNames(vec, rep(pcl, length(vec)))}))
sp_varp_rna <- do.call('c', lapply(varp_rna[ord_col_rna], function(pcl) {vec <- rownames(mcl_all)[which(cl_rna$cluster == pcl)]; 
                                    setNames(vec, rep(pcl, length(vec)))}))
save(sp_varp_atac, file = './output/mcatac/peak_indices_var_atac_clust.rda')
save(sp_varp_rna, file = './output/mcatac/peak_indices_var_rna_clust.rda')
load('./output/mcatac/peak_indices_var_atac_clust.rda')
load('./output/mcatac/peak_indices_var_rna_clust.rda')


### prego regression 
# get sequences
seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(mcl_all))
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# save(seqs_all, file = './output/mcatac/feat_peak_seqs_500bp.rda')
load('./output/mcatac/feat_peak_seqs_500bp.rda')

# Identify promoter-proximal peaks from each group
tss <- gintervals.load('intervs.global.tss')
nei_sp_atac_prom <- gintervals.neighbors(seq_coords[match(sp_varp_atac, seq_coords$peak_name),], tss, maxneighbors = 1, mindist = -1e3, maxdist = 1e3)
nei_sp_rna_prom <- gintervals.neighbors(seq_coords[match(sp_varp_rna, seq_coords$peak_name),], tss, maxneighbors = 1, mindist = -1e3, maxdist = 1e3)
sp_varp_atac_prom <- sp_varp_atac[sp_varp_atac %in% nei_sp_atac_prom$peak_name]
sp_varp_rna_prom <- sp_varp_rna[sp_varp_rna %in% nei_sp_rna_prom$peak_name]
sp_varp_atac_dist <- sp_varp_atac[!(sp_varp_atac %in% nei_sp_atac_prom$peak_name)]
sp_varp_rna_dist <- sp_varp_rna[!(sp_varp_rna %in% nei_sp_rna_prom$peak_name)]

pk_col_annot_atac <- as.data.frame(tibble::enframe(sp_varp_atac))
rownames(pk_col_annot_atac) <- pk_col_annot_atac$value
pk_col_annot_atac <- dplyr::select(pk_col_annot_atac, name)
colnames(pk_col_annot_atac) <- 'cluster_atac'

pk_col_annot_rna <- as.data.frame(tibble::enframe(sp_varp_rna))
rownames(pk_col_annot_rna) <- pk_col_annot_rna$value
pk_col_annot_rna <- dplyr::select(pk_col_annot_rna, name)
colnames(pk_col_annot_rna) <- 'cluster_rna'
ann_colors[['cluster_atac']] <- setNames(gplots::col2hex(sample(grep('white|gray|grey', colors(), inv=T, v=T), 
                                    length(unique(pk_col_annot_atac$cluster_atac)))), 
                                    unique(pk_col_annot_atac$cluster_atac))
ann_colors[['cluster_rna']] <- setNames(gplots::col2hex(sample(grep('white|gray|grey', colors(), inv=T, v=T), 
                                    length(unique(pk_col_annot_rna$cluster_rna)))), 
                                    unique(pk_col_annot_rna$cluster_rna))

p_sp_atac <- pheatmap(log2(1e-5+t(mcl_all_norm[setdiff(sp_varp_atac, sp_varp_atac_prom),amc_ord])), 
                cluster_rows = F, 
                cluster_cols = F, 
                annotation_row = row_annot, 
                annotation_col = pk_col_annot_atac, 
                annotation_colors = ann_colors,
                show_colnames = F, 
                show_rownames = F,silent = F)
p_sp_rna <- pheatmap(log2(1e-5+t(mcl_all_norm[setdiff(sp_varp_rna, sp_varp_rna_prom),amc_ord])), 
                cluster_rows = F, cluster_cols = F, 
                annotation_row = row_annot, 
                annotation_col = pk_col_annot_rna, 
                annotation_colors = ann_colors,
                show_colnames = F, 
                show_rownames = F,silent = T)
p_sp_atac_prom <- pheatmap(log2(1e-5+t(mcl_all_norm[sp_varp_atac_prom,amc_ord])), 
                cluster_rows = F, 
                cluster_cols = F, 
                annotation_row = row_annot, 
                annotation_col = pk_col_annot_atac, 
                annotation_colors = ann_colors,
                show_colnames = F, 
                show_rownames = F,silent = F)
p_sp_rna_prom <- pheatmap(log2(1e-5+t(mcl_all_norm[sp_varp_rna_prom,amc_ord])), 
                cluster_rows = F, cluster_cols = F, 
                annotation_row = row_annot, 
                annotation_col = pk_col_annot_rna, 
                annotation_colors = ann_colors,
                show_colnames = F, 
                show_rownames = F,silent = F)
save_pheatmap(p_sp_atac, './output/mcatac/figs/amc_vs_single_peaks_variable_atac_clusters.png', height=3000, width = 4000, res = 300)
save_pheatmap(p_sp_rna, './output/mcatac/figs/amc_vs_single_peaks_variable_rna_clusters.png',  height=3000, width = 4000, res = 300)
save_pheatmap(p_sp_atac_prom, './output/mcatac/figs/amc_vs_single_promoter_peaks_variable_atac_clusters.png', height=3000, width = 2500, res = 300)
save_pheatmap(p_sp_rna_prom, './output/mcatac/figs/amc_vs_single_promoter_peaks_variable_rna_clusters.png',  height=3000, width =2500, res = 300)


## Integrate reik data
gset_genome('mm10')
load( "/net/mraid14/export/tgdata/users/atanay/proj/enhflow/scdb/rna_md.Rmd")
ct_neuro <- c('Epiblast', 'Rostral neural plate',
            'Neural crest', 'Definitive ectoderm', 'Forebrain/Midbrain/Hindbrain')
ct_neuro_mc <- lapply(ct_neuro, function(ct) {y <- md$metacell[md$cell_type == ct]; as.numeric(levels(y)[y])})
names(ct_neuro_mc) <- ct_neuro
rk_trk <- grep('marginal|DELETE|smoothed|R1|R2|grouped', gtrack.ls('wt_reik'), inv=T, v=T)
rk_trk <- rk_trk[order(as.numeric(gsub('wt_reik\\.mc', '', rk_trk)))]
# aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
peaks <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(aaa$mcl_all))
reik_all_peaks <- gextract(rk_trk[unlist(ct_neuro_mc)], peaks, iterator = peaks)
reik_neuro_peaks <- subset(reik_all_peaks, select = -c(chrom, start, end, intervalID))
rownames(reik_neuro_peaks) <- misha.ext::convert_misha_intervals_to_10x_peak_names(reik_all_peaks[,1:3])
reik_neuro_peaks[is.na(reik_neuro_peaks)] <- 0
# load('./output/mcatac/peak_indices_var_atac_clust.rda')
rnn <- t(t(reik_neuro_peaks)/colSums(reik_neuro_peaks))
rnnl <- log2(1e-5 + rnn)
colnames(rnnl) <- as.numeric(gsub('wt_reik\\.mc', '', colnames(rnnl)))
color_key <- unique(md[,c('cell_type', 'color')])
row_annot_reik <- tibble::column_to_rownames(md[,c('metacell','cell_type')], 'metacell')
row_annot_reik <- as.data.frame(as.character(row_annot_reik[order(as.numeric(rownames(row_annot_reik))),]))
colnames(row_annot_reik)  <- 'cell_type'
rownames(row_annot_reik)  <- as.numeric(rownames(row_annot_reik))
ann_colors[['cell_type']] <- c(ann_colors[['cell_type']], tibble::deframe(color_key))
ann_colors[['cluster_atac']] <- setNames(gplots::col2hex(sample(grep('white|gray|grey', colors(), inv=T, v=T),
                                    length(unique(pk_col_annot_atac$cluster_atac)))),
                                    unique(pk_col_annot_atac$cluster_atac))
ac_filt <- ann_colors[c('cell_type', 'cluster_atac')]
ac_filt[['cell_type']] <- ac_filt[['cell_type']][names(ac_filt[['cell_type']]) %in% 
                c("NSC", "IPC", "Neuron", ct_neuro)]

reik_mmc_merge <- as.matrix(bind_cols(rnnl, as.matrix(log2(1e-5+as.matrix(mcl_all_norm[,amc_ord])))))
rownames(reik_mmc_merge) <- rownames(mcl_all_norm)
ra_all <- bind_rows(row_annot_reik, row_annot)
prmm <- pheatmap(t(reik_mmc_merge[sp_varp_atac_dist,]), annotation_row = ra_all, annotation_colors = ac_filt, 
                cluster_cols = F, cluster_rows = F, 
                show_rownames=F, show_colnames=F, 
                annotation_col = pk_col_annot_atac)
save_pheatmap(prmm, './output/mcatac/figs/var_atac_peaks_mmcortex_and_reik.png', 
            height = 2000, width = 3500, res = 150)

# remove unnecessary variables
rm(list = ls()[!(ls() %in% c(grep('G', ls(), v=T), 'wd', 'sp_varp_atac', 'sp_varp_rna', 'seqs_all'))])

# regress on ATAC peak clusters
res <- prego::regress_pwm.clusters(sequences=seqs_all[sp_varp_atac], clusters = as.numeric(names(sp_varp_atac)), 
                                        min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        motif_num = 3,
                                        score_metric = "ks",
                                        final_metric = "ks",
                                        multi_kmers = T,
                                        use_sge = T, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
save(res, file = './output/sequence_modeling/prego_res_variable_atac_peak_clusters.rda')
# regress on RNA peak clusters
res <- prego::regress_pwm.clusters(sequences=seqs_all[sp_varp_rna], clusters = as.numeric(names(sp_varp_rna)), 
                                        min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        motif_num = 3,
                                        score_metric = "ks",
                                        final_metric = "ks",
                                        multi_kmers = T,
                                        use_sge = T, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
save(res, file = './output/sequence_modeling/prego_res_variable_rna_peak_clusters.rda')

## Regress only NSC clusters vs random genome
gen_random_genome_intervals <- function(num_peaks = 1e+5,
                                                peak_width = 2e+2,
                                                bp_from_chrom_edge_to_avoid = 3e+6,
                                                datasets_of_interest = NULL,
                                                motif_regex = NULL,
                                                motif_tracks = NULL,
                                                parallel = TRUE) {
    ALLGENOME[[1]] <- ALLGENOME[[1]][!grepl("_", ALLGENOME[[1]]$chrom), ]
    chrom_lens <- apply(ALLGENOME[[1]][, 2:3], 1, diff)
    chrom_fracs <- setNames(chrom_lens / sum(chrom_lens), ALLGENOME[[1]][, 1])
    sample_seqs <- mapply(chrom_fracs, names(chrom_fracs), chrom_lens, FUN = function(x, y, z) {
        chrom <- rep(y, round(x * num_peaks))
        start <- sample.int(n = z, size = length(chrom))
        end <- start + peak_width
        return(as.data.frame(rbind(chrom, start, end)))
    })
    sample_seqs <- as.data.frame(do.call("rbind", lapply(sample_seqs, t)))
    sample_seqs[, 2:3] <- apply(sample_seqs[, 2:3], 2, as.numeric)
    sample_seqs <- sample_seqs[with(sample_seqs, order(chrom, start, end)), ]
    end_shift <- ALLGENOME[[1]][match(sample_seqs$chrom, ALLGENOME[[1]][, 1]), 3] - bp_from_chrom_edge_to_avoid
    sample_seqs <- sample_seqs[sample_seqs$start >= bp_from_chrom_edge_to_avoid & sample_seqs$end <= end_shift, ]
    sample_seqs <- PeakIntervals(sample_seqs)
    return(sample_seqs)
}
rg_ints_tmp <- gen_random_genome_intervals(num_peaks = 5e+3, peak_width = 500)
seqs_rg <- gseq.extract(rg_ints_tmp)
load('./output/mcatac/peak_indices_var_atac_clust.rda')
seqs_cl_nsc <- seqs_all[sp_varp_atac[as.numeric(names(sp_varp_atac)) %in% 1:5]]
seqs_both <- c(seqs_cl_nsc, seqs_rg)
cl_both <- c(names(sp_varp_atac)[as.numeric(names(sp_varp_atac)) %in% 1:5], rep(6, length(seqs_rg)))
rm(list = ls()[!(ls() %in% c(grep('G', ls(), v=T), 'wd', 'seqs_both', 'cl_both'))])
res_nsc <- prego::regress_pwm.clusters(sequences = seqs_both, 
                                clusters = cl_both, 
                                min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        # motif_num = 3,
                                        score_metric = "ks",
                                        final_metric = "r2",
                                        multi_kmers = T,
                                        use_sge = F, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
save(res_nsc, file = './output/sequence_modeling/prego_res_nsc_atac_clusters_f_met_r2.rda')

# Load predictions

load('./output/sequence_modeling/prego_res_variable_atac_peak_clusters.rda')
resa <- res
load('./output/sequence_modeling/prego_res_variable_rna_peak_clusters.rda')
resr <- res

rownames(resa$pred_mat) <- rownames(mcl_all)[sp_varp_atac]
rownames(resr$pred_mat) <- rownames(mcl_all)[sp_varp_rna]


## Plot predictions of consensus motifs
resapm <- t(resa$pred_mat)
resapm <- resapm - rowMaxs(resapm)
resapm[resapm < -15] <- -15

rownames(resapm) <- apply(resa$stats[,c('cluster', 'consensus')], 1, paste, collapse = ': ')
p_resapm <- pheatmap(resapm, annotation_col = pk_col_annot_atac, fontsize_row = 24,
                    annotation_colors = ann_colors, 
                    cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
p_resapm_prom <- pheatmap(resapm[,sp_varp_atac_prom], annotation_col = pk_col_annot_atac, fontsize_row = 24,
                    annotation_colors = ann_colors, 
                    cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
save_pheatmap(p_resapm, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_atac_clusters_1motifs_per.png',  height=1600, width = 2800, res = 150)
save_pheatmap(p_resapm_prom, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_promoters_from_atac_clusters.png',  height=1600, width = 2800, res = 150)

resrpm <- t(resr$pred_mat)
resrpm <- resrpm - rowMaxs(resrpm)
resrpm[resrpm < -15] <- -15

rownames(resrpm) <- apply(resr$stats[,c('cluster', 'consensus')], 1, paste, collapse = ': ')
p_resrpm <- pheatmap(resrpm, annotation_col = pk_col_annot_rna, fontsize_row = 24,
                    annotation_colors = ann_colors, 
                    cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
p_resrpm_prom <- pheatmap(resrpm[,sp_varp_rna_prom], annotation_col = pk_col_annot_rna, fontsize_row = 24,
                    annotation_colors = ann_colors, 
                    cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
save_pheatmap(p_resrpm, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_rna_clusters_1motifs_per.png',  height=1600, width = 2800, res = 150)
save_pheatmap(p_resrpm_prom, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_promoters_from_RNA_clusters.png',  height=1600, width = 2800, res = 150)

## Count #q>0.98 hits
count_quantile_hits_in_clusters <- function(mat, cl, q_thresh = 0.98) {
    mat_q <- plyr::aaply(mat, 1, function(x) ecdf(x)(x), .parallel = T)
    # print(dim(mat_q))
    mat_count_q_thresh <- tgs_matrix_tapply(mat_q, cl, function(x) sum(x >= q_thresh)/(length(x)*(1-q_thresh)))    
    # print(dim(mat_count_q_thresh))
    return(mat_count_q_thresh)
}
# resapm_q <- plyr::aaply(resapm, 1, function(x) ecdf(x)(x), .parallel = T)
# resapm_q_count_98 <- tgs_matrix_tapply(resapm_q, names(sp_varp_atac), function(x) sum(x >= 0.98))
resapm_q_count_98 <- count_quantile_hits_in_clusters(t(resa$pred_mat), names(sp_varp_atac))
pheatmap(resapm_q_count_98, main = "Fold enrichment over random of #q>0.98 hits in clusters")
## Compute predictions - all models
preds_all_atac <- plyr::llply(resa$models, function(x) plyr::llply(x$models, function(y) prego::compute_pwm(seqs_all[sp_varp_atac], pssm = y$pssm, spat = y$spat, prior= 0.1), .parallel = T), .parallel = T)
preds_all_rna <- plyr::llply(resr$models, function(x) plyr::llply(x$models, function(y) prego::compute_pwm(seqs_all[sp_varp_rna], pssm = y$pssm, spat = y$spat, prior= 0.1), .parallel = T), .parallel = T)
# preds_all_atac <- plyr::llply(resa$models, function(x) plyr::llply(x$models, function(y) y$pred))
# preds_all_rna <- plyr::llply(resr$models, function(x) plyr::llply(x$models, function(y) y$pred))
paa_mat <- t(do.call('cbind', lapply(preds_all_atac, function(x) do.call('cbind', x))))
par_mat <- t(do.call('cbind', lapply(preds_all_rna, function(x) do.call('cbind', x))))
paa_mat <- paa_mat - rowMaxs(paa_mat)
par_mat <- par_mat - rowMaxs(par_mat)
paa_cons_list <- unlist(plyr::llply(resa$models, function(x) plyr::llply(x$models, function(y) y$consensus)))
par_cons_list <- unlist(plyr::llply(resr$models, function(x) plyr::llply(x$models, function(y) y$consensus)))
paa_cons_list <- ifelse(is.na(paa_cons_list),'None',paa_cons_list)
par_cons_list <- ifelse(is.na(par_cons_list),'None',par_cons_list)

colnames(paa_mat) <- rownames(mcl_all)[sp_varp_atac]
rownames(paa_mat) <- apply(cbind(rep(varp_atac[ord_col_atac], each=3),
                                rep(1:3, length(varp_atac)), 
                                paa_cons_list),1,
                            stringr::str_c, collapse = '_')
colnames(par_mat) <- rownames(mcl_all)[sp_varp_rna]
rownames(par_mat) <- apply(cbind(rep(varp_rna[ord_col_rna], each=3),
                            rep(1:3, length(varp_rna)),
                            par_cons_list),1,
                            stringr::str_c, collapse = '_')
paa_row_annot <- as.data.frame(rep(varp_atac[ord_col_atac], each=3))
rownames(paa_row_annot) <- rownames(paa_mat)
par_row_annot <- as.data.frame(rep(varp_rna[ord_col_rna], each=3))
rownames(par_row_annot) <- rownames(par_mat)
colnames(paa_row_annot) <- 'cluster_atac'
colnames(par_row_annot) <- 'cluster_rna'
p_paa <- pheatmap(paa_mat[,], annotation_col = pk_col_annot_atac,
                    annotation_row = paa_row_annot, annotation_colors = ann_colors, 
                    cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
p_par <- pheatmap(par_mat[,], annotation_col = pk_col_annot_rna, 
                        annotation_row = par_row_annot, annotation_colors = ann_colors, 
                        cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
save_pheatmap(p_paa, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_atac_clusters_3motifs_per.png',  height=1600, width = 2800, res = 150)
save_pheatmap(p_par, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_rna_clusters_3motifs_per.png',  height=1600, width = 2800, res = 150)

p_pred_atac <- pheatmap(2**(t(resa$pred_mat) - colMaxs(resa$pred_mat)), annotation_col = pk_col_annot_atac, 
                        # annotation_row = paa_row_annot,
                         annotation_colors = ann_colors, 
                        cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
p_pred_rna <- pheatmap(2**(t(resr$pred_mat) - colMaxs(resr$pred_mat)), annotation_col = pk_col_annot_rna, 
                        # annotation_row = par_row_annot,
                         annotation_colors = ann_colors, 
                        cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F)
save_pheatmap(p_pred_atac, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_atac_clusters_1motifs_per.png',  height=1600, width = 2800, res = 150)
save_pheatmap(p_pred_rna, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_rna_clusters_1motifs_per.png',  height=1600, width = 2800, res = 150)


## Get pwm energies of kmers
consensus_to_pssm <- function(cons) {
        cons_s <- unlist(stringr::str_split(cons, ''))
        pssm <- matrix(0, nrow = length(cons_s), ncol = 4)
        colnames(pssm) <- c("A", "C", "G", "T")
        for (i in 1:length(cons_s)) {
            # print(i)
            # print(cons_s[[i]])
            if (cons_s[[i]] %in% colnames(pssm)) {
                pssm[i, cons_s[[i]]] <- 1
            } else if (cons_s[[i]] == "*" | cons_s[[i]] == "N") {
                pssm[i, ] <- 0.25
            } else if (cons_s[[i]] %in% c('M', 'R', 'W', 'S', 'Y', 'K')) {
                switch(cons_s[[i]],
                       M = pssm[i, c("A", "C")] <- 0.5,
                       R = pssm[i, c("A", "G")] <- 0.5,
                       W = pssm[i, c("A", "T")] <- 0.5,
                       S = pssm[i, c("C", "G")] <- 0.5,
                       Y = pssm[i, c("C", "T")] <- 0.5,
                       K = pssm[i, c("G", "T")] <- 0.5
                )
            }
            else {
                cli_abort("Unknown character number {.val {i}} in consensus: {.val {cons_s[[i]]}}")
            }
        }
        return(pssm)
}

df_from_consensus_seqs <- function(cons_list) {
    cons_df <- lapply(1:length(cons_list), function(ci, i) {
        cmat <- as.data.frame(consensus_to_pssm(ci[[i]])) %>%
                dplyr::mutate(key = i, pos = 1:nrow(.)) %>%
                dplyr::relocate(key, pos, everything())
    }, ci = cons_list)
    return(do.call('rbind',cons_df))
}


c_df <- df_from_consensus_seqs(unique(unlist(lapply(resa$models, function(x) x$kmers))))
test_pwm_kmer <- plyr::laply(unique(c_df$key), function(i) prego::compute_pwm(sequences = seqs_all[sp_varp_atac], 
                                                            pssm = dplyr::filter(c_df, key == i), 
                                                            prior = 0.01), .parallel=T)

kmp_lq <- plyr::aaply(test_pwm_kmer, 1, function(x) {
    vec <- -log2(1 - ecdf(x)(x)); 
    vec[is.infinite(vec)] <- max(vec[!is.infinite(vec)]); 
    return(vec)
    }, .parallel=T)

p_kmp_lq <- pheatmap(kmp_lq, cluster_cols = F, show_colnames = F, annotation_col = pk_col_annot_atac, 
                     annotation_colors = ann_colors, cluster_rows = T, show_rownames = T)


### Plot motifs of expressed TFs in NSCs

### Get peaks
# mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks.rds')

# ### Get cluster assignments
mca_rna_km <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
scatac_km <- aaa$km_a_legc

intervs_energy <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')
rg_energy <- readRDS('./output/sequence_modeling/random_genome_25k_motif_energy.rds')

#' modified from mcATAC::calculate_d_stats
calculate_d_stats <- function(pssm_fg, pssm_bg, fg_clustering = NULL, parallel = getOption("mcatac.parallel"), alternative = "greater", nc = getOption("mcatac.parallel.nc")) {
    cols_fg <- grep("chrom|start|end$|interval|peak_name", colnames(pssm_fg), ignore.case = T, invert = T, value = T)
    cols_bg <- grep("chrom|start|end$|interval|peak_name", colnames(pssm_bg), ignore.case = T, invert = T, value = T)
    cols_both <- intersect(cols_fg, cols_bg)
    if (!parallel) {
        nc <- pmax(2, round(0.1 * nc))
    }
    if (!is.null(fg_clustering)) {
        ks_test_results <- parallel::mclapply(cols_both, FUN = function(x, i) {
            return(tapply(x[, i], fg_clustering, function(y) suppressWarnings(ks.test(y, pssm_bg[, i], alternative = alternative))))
        }, mc.cores = nc, x = pssm_fg)
        ks_d <- sapply(ks_test_results, function(x) sapply(x, function(y) y$statistic))
        colnames(ks_d) <- cols_both
    } else {
        ks_test_results <- parallel::mclapply(cols_both, function(x) {
            ks.test(x = pssm_fg[, x], y = pssm_bg[, x], alternative = alternative)
        },
        mc.cores = nc
        )
        ks_d <- setNames(sapply(ks_test_results, function(x) x$statistic), cols_both)
    }
    return(ks_d)
}


d_atac <- calculate_d_stats(pssm_fg = intervs_energy, pssm_bg = rg_energy, fg_clustering = scatac_km$cluster, alternative = 'less')
d_rna <- calculate_d_stats(pssm_fg = intervs_energy, pssm_bg = rg_energy, fg_clustering = mca_rna_km$cluster, alternative = 'less')

nsc_mc <- which(mcmd$cell_type == 'NSC')
inds_d <- which(colMaxs(d_atac) >= 0.4)

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- c('Bcl11b', tfs$Symbol)
tfs <- tfs[tfs %in% rownames(mc_rna@e_gc)]
nsc_tfs <- tfs[tfs %in% rownames(mc_rna@mc_fp)[rowMaxs(mc_rna@mc_fp[,nsc_mc]) > 3]]
fp_tfs <- tfs[tfs %in% rownames(mc_rna@mc_fp)[rowMaxs(mc_rna@mc_fp) > 3]]
# nsc_tf_motif <- unlist(plyr::llply(nsc_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(d_atac)[inds_d], v=T, ign=T), .parallel = T))
nsc_tf_motif <- unlist(plyr::llply(fp_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(d_atac)[inds_d], v=T, ign=T), .parallel = T))
nsc_motif_pwm <- t(as.matrix(intervs_energy[match(sp_varp_atac_dist, seq_coords$peak_name),nsc_tf_motif]))
colnames(nsc_motif_pwm) <- misha.ext::convert_misha_intervals_to_10x_peak_names(intervs_energy[,1:3])[sp_varp_atac_dist]
nsc_motif_pwm_prom <- t(as.matrix(intervs_energy[match(sp_varp_atac_prom, seq_coords$peak_name),nsc_tf_motif]))
colnames(nsc_motif_pwm_prom) <- misha.ext::convert_misha_intervals_to_10x_peak_names(intervs_energy[,1:3])[sp_varp_atac_prom]

count_quantile_hits_in_clusters <- function(mat, cl, q_thresh = 0.98) {
    mat_q <- plyr::aaply(mat, 1, function(x) ecdf(x)(x), .parallel = T)
    mat_count_q_thresh <- tgs_matrix_tapply(mat_q, cl, function(x) sum(x >= q_thresh)/(length(x)*(1-q_thresh)))    
    return(mat_count_q_thresh)
}
nmc_lfc <- log2(1e-2 + count_quantile_hits_in_clusters(nsc_motif_pwm, names(sp_varp_atac_dist)))
nmc_lfc_prom <- log2(1e-2 + count_quantile_hits_in_clusters(nsc_motif_pwm, names(sp_varp_atac_prom)))
brks <- c(seq(min(nmc_lfc), 0, l=50), seq(0.01,max(nmc_lfc),l=51))
clrmp <- colorRampPalette(c('blue4', 'white', 'red4'))(100)

p_fp_tf_enrich <- pheatmap(t(nmc_lfc[as.character(varp_atac[ord_col_atac]),colMaxs(abs(nmc_lfc)) >= 1.5]),
            cluster_cols=F, color = clrmp, breaks = brks, fontsize_row = 6)
p_fp_tf_enrich <- pheatmap(t(nmc_lfc[,colMaxs(abs(nmc_lfc)) >= 1]),
            cluster_cols=F, color = clrmp, breaks = brks, fontsize_row = 6)
save_pheatmap(p_fp_tf_enrich, './output/sequence_modeling/figs/motif_enrichment_of_high_fp_tfs.png',  height=2000, width = 4500, res = 300)


### Regress on enriched motifs
motifs_regress <- setNames(colnames(nmc_lfc)[apply(nmc_lfc, 1, which.max)], rownames(nmc_lfc))
amd <- all_motif_datasets()
res_reg <- plyr::llply(names(motifs_regress), function(i) {
    prego::regress_pwm(sequences=seqs_all[match(sp_varp_atac_dist, rownames(mcl_all))], 
                            response = ifelse(names(sp_varp_atac_dist) == i, 1, 0), 
                            motif = dplyr::filter(amd, motif == motifs_regress[[i]]),
                                        min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        score_metric = "r2",
                                        final_metric = "ks",
                                        use_sge = F, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
}, .parallel = T)
res_ks_df <- as.data.frame(list(input_motif = motifs_regress, 
                            ks_reg = sapply(res_reg, function(x) x$ks$statistic), 
                            ks_db_match = sapply(res_reg, function(x) x$db_match_ks$statistic),
                            db_match_motif = sapply(res_reg, function(x) x$db_match)))
# res_ks_df <- tibble::rownames_to_column(res_ks_df[order(match(rownames(res_ks_df), varp_atac[ord_col_atac])),], var='cluster')
# plot(res_ks_df$ks_db_match, res_ks_df$ks_reg, pch = 16, xlab = 'KS D of DB-matched motif', ylab = 'KS D of Regressed motif')
# abline(0,1,col='red',lty=2)

res_dn <- prego::regress_pwm.clusters(sequences=seqs_all[match(sp_varp_atac_dist, rownames(mcl_all))], 
                            clusters = names(sp_varp_atac_dist), 
                                        min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        score_metric = "r2",
                                        final_metric = "ks",
                                        use_sge = F, 
                                        match_with_db=T, 
                                        use_sample = T, 
                                        parallel=T)
res_dn_ks_df <- as.data.frame(list(input_motif = motifs_regress, 
                            ks_reg = sapply(res_dn, function(x) x$ks$statistic), 
                            ks_db_match = sapply(res_dn, function(x) x$db_match_ks$statistic),
                            db_match_motif = sapply(res_dn, function(x) x$db_match)))
# res_ks_df <- tibble::rownames_to_column(res_ks_df[order(match(rownames(res_ks_df), varp_atac[ord_col_atac])),], var='cluster')
# plot(res_dn_ks_df$ks_db_match, res_ks_df$ks_reg, pch = 16, xlab = 'KS D of DB-matched motif', ylab = 'KS D of Regressed motif')
# abline(0,1,col='red',lty=2)
library(gridExtra)
motif_plots_seed <- lapply(c(res_ks_df$input_motif, res_ks_df$db_match_motif), function(nm) prego::plot_pssm_logo(filter(amd, motif == nm), title = nm))
names(motif_plots_seed) <- c(res_ks_df$input_motif, res_ks_df$db_match_motif)
gall_seed <- marrangeGrob(grobs = motif_plots_seed, ncol =2, nrow = length(motif_plots_seed)/2)

motif_plots_dn <- lapply(res_dn$stats$cluster, function(nm) prego::plot_pssm_logo(filter(res_dn$motif_dataset, motif == nm), title = paste0('cluster_', nm)))
names(motif_plots_dn) <- paste0('cluster_', res_dn$stats$cluster)
gall_dn <- marrangeGrob(grobs = motif_plots_dn, ncol =1, nrow = length(motif_plots_dn))

ggsave('./output/sequence_modeling/figs/motifs_inferred_with_seed_per_cluster.png', gall_seed, height = 40, width = 16, units = 'in')
ggsave('./output/sequence_modeling/figs/motifs_inferred_de_novo_per_cluster.png', gall_dn, height = 40, width = 8, units = 'in')




nsc_motif_pwm_norm <- nsc_motif_pwm
nsc_motif_pwm_norm <- nsc_motif_pwm_norm - rowMaxs(nsc_motif_pwm_norm)
nsc_motif_pwm_norm[nsc_motif_pwm_norm < -15] <- -15

nsc_motif_pwm_rna <- t(as.matrix(intervs_energy[sp_varp_rna,nsc_tf_motif]))
colnames(nsc_motif_pwm_rna) <- misha.ext::convert_misha_intervals_to_10x_peak_names(intervs_energy[,1:3])[sp_varp_rna]
nsc_motif_pwm_rna_norm <- nsc_motif_pwm_rna
nsc_motif_pwm_rna_norm <- nsc_motif_pwm_rna_norm - rowMaxs(nsc_motif_pwm_rna_norm)
nsc_motif_pwm_rna_norm[nsc_motif_pwm_rna_norm < -15] <- -15
# nsc_motif_lq <- plyr::aaply(nsc_motif_pwm, 1, function(x) {
#     vec <- -log2(1 - ecdf(x)(x)); 
#     vec[is.infinite(vec)] <- max(vec[!is.infinite(vec)]); 
#     return(vec)}, 
#     .parallel=T)
# colnames(nsc_motif_lq) <- colnames(nsc_motif_pwm)

pk_col_annot_atac <- as.data.frame(tibble::enframe(sp_varp_atac))
rownames(pk_col_annot_atac) <- colnames(nsc_motif_pwm)
pk_col_annot_atac <- dplyr::select(pk_col_annot_atac, name)
colnames(pk_col_annot_atac) <- 'cluster_atac'
ann_colors[['cluster_atac']] <- setNames(gplots::col2hex(sample(grep('white|gray|grey', colors(), inv=T, v=T), 
                                    length(unique(pk_col_annot_atac$cluster_atac)))), 
                                    unique(pk_col_annot_atac$cluster_atac))

p_pwm_norm <- pheatmap::pheatmap(nsc_motif_pwm_norm, cluster_cols = F, show_colnames=F,fontsize_row = 6,
                annotation_col = pk_col_annot_atac, annotation_colors=ann_colors)
save_pheatmap(p_pwm_norm, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_motifs_atac_clusters.png',  height=1600, width = 2800, res = 150)
p_pwm_norm_prom <- pheatmap::pheatmap(nsc_motif_pwm_norm[,sp_varp_atac_prom], cluster_cols = F, show_colnames=F,fontsize_row = 6,
                annotation_col = pk_col_annot_atac, annotation_colors=ann_colors)
save_pheatmap(p_pwm_norm_prom, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_motifs_promoters_of_atac_clusters.png',  height=1600, width = 2800, res = 150)

p_pwm_rna_norm <- pheatmap::pheatmap(nsc_motif_pwm_rna_norm, cluster_cols = F, show_colnames=F,fontsize_row = 6,
                annotation_col = pk_col_annot_rna, annotation_colors=ann_colors)
save_pheatmap(p_pwm_rna_norm, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_motifs_rna_clusters.png',  height=1600, width = 2800, res = 150)
p_pwm_rna_norm_prom <- pheatmap::pheatmap(nsc_motif_pwm_rna_norm[,sp_varp_rna_prom], cluster_cols = F, show_colnames=F,fontsize_row = 6,
                annotation_col = pk_col_annot_rna, annotation_colors=ann_colors)
save_pheatmap(p_pwm_rna_norm_prom, './output/sequence_modeling/figs/cluster_comparison/predicted_energies_motifs_promoters_of_rna_clusters.png',  height=1600, width = 2800, res = 150)


###########



dir.create('./output/sequence_modeling/figs/cluster_comparison/var_peaks_atac')
dir.create('./output/sequence_modeling/figs/cluster_comparison/var_peaks_rna')
tt <- sapply(1:length(res$models), function(i) {
    mdl <- res$models[[i]]
    pathi <- glue::glue(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/var_peaks_atac/cl_{i}.png'))
    suppressWarnings(png(pathi, width = 3000, height = 3000, res = 150))
    print(plot_regression_qc_multi(mdl))
    dev.off()
})

ttt <- sapply(1:length(res$models), function(i) {
    mdl <- res$models[[i]]
    pathi <- glue::glue(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/var_peaks_rna/cl_{i}.png'))
    png(pathi, width = 3000, height = 3000, res = 150)
    print(suppressWarnings(plot_regression_qc_multi(mdl)))
    dev.off()
})

library(officer)

print_slide = function(i, folder_name, ppt) {
    path_hm = paste0(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/', folder_name), glue::glue('/cl_{i}.png'))
    print(path_hm)
    img_hm = external_img(src = path_hm, height = 7, width = 5)
#     ppt = add_slide(ppt)
    ppt = add_slide(ppt, layout = "Blank", master = "Office Theme")
    ppt = ph_with(x = ppt, value = img_hm,
                    location = ph_location(left = 0, 
                    top = 0, width = 9, height = 7
                    ))
    # ppt = ph_with(x = ppt, value = img_2d, location = ph_location(left = 0.25, 
                            # top = 2, width = 4.2, height = 4.2))
}

raw_cluster_files <- list.files('./output/sequence_modeling/figs/cluster_comparison/var_peaks_atac/', pattern = '.png')
rna_cluster_files <- list.files('./output/sequence_modeling/figs/cluster_comparison/var_peaks_rna/', pattern = '.png')
ppt = read_pptx()
bb <- lapply(seq_along(raw_cluster_files), function(i) print_slide(i, 'var_peaks_atac', ppt))
dir.create('output/sequence_modeling/ppt/')
print(ppt, target = file.path(wd, 'output/sequence_modeling/ppt/var_atac_peak_clusters.pptx'))
ppt = read_pptx()
cc <- lapply(seq_along(rna_cluster_files), function(i) print_slide(i, 'var_peaks_rna', ppt))
print(ppt, target = file.path(wd, 'output/sequence_modeling/ppt/var_rna_peak_clusters.pptx'))

mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks.rds')
c2mc <- mca@cell_to_metacell
amc_prom <- t(tgs_matrix_tapply(mat_prom@mat[,c2mc$cell_id], c2mc$metacell, sum))
a_legc_prom <- log2(1e-5 + t(t(as.matrix(amc_prom))/Matrix::colSums(amc_prom)))
cor_famc_rmc <- tgs_cor(mc_rna@mc_fp[feats,], a_legc_prom[feats,], spearman=T)
p_cor_famc_rmc <- pheatmap(t(cor_famc_rmc[cust_mc_ord_st, cust_mc_ord_st]), 
                    color = colorRampPalette(c('blue4', 'white', 'red4'))(100),
                                    breaks = seq(-.3,.3, l = 101),
                    cluster_cols = F, cluster_rows = F, 
                    show_rownames=F, show_colnames=F, 
                    annotation_row = col_annot, annotation_col = col_annot, 
                    annotation_colors=ann_colors)
save_pheatmap(p_cor_famc_rmc, './output/mcatac/figs/cor_famc_rmc.png', height=3000, width = 3500, res = 300)


egc_norm <- t(t(mca@egc)/colSums(mca@egc))
legc_cl_atac <- log2(1e-5 + tgs_matrix_tapply(t(egc_norm), cl_raw$cluster, mean))
legc_cl_rna <- log2(1e-5 + tgs_matrix_tapply(t(egc_norm), cl_rna$cluster, mean))
peak_rna_ord <- order(match(mcmd$cell_type[apply(legc_cl_rna, 1, which.max)], cust_st_ord))
peak_atac_ord <- order(match(mcmd$cell_type[apply(legc_cl_atac, 1, which.max)], cust_st_ord))
p_famc_atac_pk <- pheatmap(t(legc_cl_atac[peak_atac_ord,cust_mc_ord_st]), main = 'log2(UMI) - ATAC clusters', cluster_cols = T, cluster_rows = F, annotation_row = col_annot, annotation_colors = ann_colors)
p_famc_rna_pk <- pheatmap(t(legc_cl_rna[peak_rna_ord,cust_mc_ord_st]), main = 'log2(UMI) - RNA clusters', cluster_cols = T, cluster_rows = F, annotation_row = col_annot, annotation_colors = ann_colors)
save_pheatmap(p_famc_atac_pk, './output/mcatac/figs/famc_vs_atac_cl_peaks.png', height=3000, width = 3000, res = 300)
save_pheatmap(p_famc_rna_pk, './output/mcatac/figs/famc_vs_rna_cl_peaks.png', height=3000, width = 3000, res = 300)


mat_varp_atac <- as.matrix(mca@mat[sp_varp_atac,])
a_legc_varp_atac <- log2(1e-5 + t(t(mat_varp_atac)/colSums(mat_varp_atac)))
p_a_legc_varp_atac <- pheatmap(t(a_legc_varp_atac[,cust_mc_ord_st]), 
                            cluster_cols = F, cluster_rows=F, 
                            show_rownames=F, show_colnames=F, 
                            annotation_row = col_annot, annotation_colors = ann_colors)
save_pheatmap(p_sp_rna, './output/mcatac/figs/amc_vs_single_peaks_variable_rna_clusters.png',  height=3000, width = 3000, res = 300)





p_max <- pheatmap(mat_tbl[ord_row_max,ord_col_max], color = colorRampPalette(c('white','red', 'black'))(100), cluster_rows = F, cluster_cols = F)
p_max_norm <- pheatmap(mat_tbl_norm_col[ord_row_max_norm,ord_col_max], color = colorRampPalette(c('white','red', 'black'))(100), cluster_rows = F, cluster_cols = F)

p_raw_ord_rowsum <- pheatmap(mat_tbl[ord,ord_col], 
                    cluster_cols = F, 
                    cluster_rows=F, 
                    color = colorRampPalette(c('white','red', 'black'))(100))
p_norm_ord_rowsum <- pheatmap(mat_tbl_norm_col[ord_row_norm,ord_col], 
                    cluster_cols = F, 
                    cluster_rows=F, 
                    color = colorRampPalette(c('white','red', 'black'))(100))

save_pheatmap(p_max, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table_ord_rowmax.png', height=1600, width = 1600)
save_pheatmap(p_max_norm, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table_norm_cols_ord_rowmax.png', height=1600, width = 1600)

save_pheatmap(p_raw_ord_rowsum, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table_ord_rowsum.png', height=1600, width = 1600)
save_pheatmap(p_norm_ord_rowsum, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table_norm_cols_ord_rowsum.png', height=1600, width = 1600)
save_pheatmap(p_raw, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table.png', height=1600, width = 1600)
save_pheatmap(p_norm, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_joint_table_norm_cols.png', height=1600, width = 1600)

res_raw <- readRDS(file.path(wd, 'output/sequence_modeling/mmc_mcl_feat_peak_motif_reg.rds'))
res_rna <- readRDS(file.path(wd, 'output/sequence_modeling/mmcortex_mcatac_rna_match_feat_peak_res.rds'))

pssm_raw <- lapply(res_raw$models, function(x) x$pssm)
pssm_rna <- lapply(res_rna$models, function(x) x$pssm)
diff_mat <- sapply(pssm_raw, function(xi) sapply(pssm_rna, function(yi) prego::pssm_diff(xi, yi)))

sapply(seq(1,9,2), function(xi) sapply(seq(2,10,2), function(yi) xi*yi))

clrmp <- colorRampPalette(c('white', 'yellow', 'red', 'blue','black'))(100)

brks <- quantile(diff_mat, (0:100)/100)
p_kl <- pheatmap(diff_mat, main = 'KL divergence between inferred motifs', 
                color = clrmp, 
                breaks =brks)
save_pheatmap(p_kl, './output/sequence_modeling/figs/cluster_comparison/kl_diffs_between_inferred_motifs.png', 
                    height=1600, width = 1600)
ks_all <- ks.test(res_raw$stats$ks_D, res_rna$stats$ks_D)

ecdf_raw <- ecdf(res_raw$stats$ks_D)
ecdf_rna <- ecdf(res_rna$stats$ks_D)
png('./output/sequence_modeling/figs/cluster_comparison/ecdf_raw_vs_rna_clusters.png', h = 800, w = 800)
plot(ecdf_raw, col = 'red', verticals=T, do.points = F, xlab = 'KS_D')
plot(ecdf_rna, add=T, col=  'blue', verticals=T, do.points = F)
legend('left', legend = c('raw_clusters', 'rna_clusters'), col = c('red', 'blue'), lty = 1, lwd = 1)
dev.off()
# dir.create('./output/sequence_modeling/figs/')
# dir.create('./output/sequence_modeling/figs/cluster_comparison')
# dir.create('./output/sequence_modeling/figs/cluster_comparison/raw_atac')
# dir.create('./output/sequence_modeling/figs/cluster_comparison/rna_match')
colnames(res_raw$pred_mat) <- paste0('raw_', colnames(res_raw$pred_mat))
colnames(res_rna$pred_mat) <- paste0('rna_', colnames(res_rna$pred_mat))
cor_pred <- tgs_cor(res_raw$pred_mat, res_rna$pred_mat)
cor_pred_raw <- tgs_cor(res_raw$pred_mat)
cor_pred_rna <- tgs_cor(res_rna$pred_mat)

p_cor <- pheatmap(t(cor_pred), color = colorRampPalette(c('blue', 'green4', 'white','yellow4', 'red'))(100),
                    breaks = c(seq(-1,1,l=101)))
p_cor_raw <- pheatmap(t(cor_pred_raw), color = colorRampPalette(c('blue', 'green4', 'white','yellow4', 'red'))(100),
                    breaks = c(seq(-1,1,l=101)))
p_cor_rna <- pheatmap(t(cor_pred_rna), color = colorRampPalette(c('blue', 'green4', 'white','yellow4', 'red'))(100),
                    breaks = c(seq(-1,1,l=101)))

save_pheatmap(p_cor, './output/sequence_modeling/figs/cluster_comparison/correlation_between_inferred_motif_energies_both.png', 
                    height=1600, width = 1600)
save_pheatmap(p_cor_raw, './output/sequence_modeling/figs/cluster_comparison/correlation_between_inferred_motif_energies_raw.png', 
                    height=1600, width = 1600)
save_pheatmap(p_cor_rna, './output/sequence_modeling/figs/cluster_comparison/correlation_between_inferred_motif_energies_rna.png', 
                    height=1600, width = 1600)


tt <- sapply(1:length(res_raw$models), function(i) {
    mdl <- res_raw$models[[i]]
    pathi <- glue::glue(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/raw_atac/cl_{i}.png'))
    suppressWarnings(png(pathi, width = 1400, height = 1000))
    print(plot_regression_qc(mdl))
    dev.off()
})

ttt <- sapply(1:length(res_rna$models), function(i) {
    mdl <- res_rna$models[[i]]
    pathi <- glue::glue(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/rna_match/cl_{i}.png'))
    png(pathi, width = 1400, height = 1000)
    print(suppressWarnings(plot_regression_qc(mdl)))
    dev.off()
})


library(officer)

print_slide = function(i, folder_name, ppt) {
    path_hm = paste0(file.path(wd, 'output/sequence_modeling/figs/cluster_comparison/', folder_name), glue::glue('/cl_{i}.png'))
    print(path_hm)
    img_hm = external_img(src = path_hm, height = 7, width = 5)
#     ppt = add_slide(ppt)
    ppt = add_slide(ppt, layout = "Blank", master = "Office Theme")
    ppt = ph_with(x = ppt, value = img_hm,
                    location = ph_location(left = 1, 
                    top = 0.5, width = 9, height = 7
                    ))
    # ppt = ph_with(x = ppt, value = img_2d, location = ph_location(left = 0.25, 
                            # top = 2, width = 4.2, height = 4.2))
}

raw_cluster_files <- list.files('./output/sequence_modeling/figs/cluster_comparison/raw_atac/', pattern = '.png')
rna_cluster_files <- list.files('./output/sequence_modeling/figs/cluster_comparison/rna_match/', pattern = '.png')
ppt = read_pptx()
bb <- lapply(seq_along(raw_cluster_files), function(i) print_slide(i, 'raw_atac', ppt))
dir.create('output/sequence_modeling/ppt/')
print(ppt, target = file.path(wd, 'output/sequence_modeling/ppt/raw_clusters.pptx'))
ppt = read_pptx()
cc <- lapply(seq_along(rna_cluster_files), function(i) print_slide(i, 'rna_match', ppt))
print(ppt, target = file.path(wd, 'output/sequence_modeling/ppt/rna_clusters.pptx'))