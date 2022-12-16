library(prego)
library(pheatmap)
library(ComplexHeatmap)
library(misha.ext)
library(metacell)
library(matrixStats)
devtools::load_all("~/src/mcATAC/")
doMC::registerDoMC(cores = 70)
misha.ext::gset_genome('mm10')

SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
options(gmax.data.size = 1e+9)
mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]],
                                                                  rep(st, length(which(mcmd$cell_type == st))))))


col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))

mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks_test.rds')
mca@peaks$peak_name <- peak_names(subset(mca@peaks, select = -peak_name), tad_based = F)
rownames(mca@egc) <- mca@peaks$peak_name
a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))
mc_rna <- scdb_mc('pl_cort')
load('./output/mcatac/km_a_legc_k=80.rda')


seq_coords <- mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]
seqs_all <- unlist(sapply(1:ceiling(nrow(seq_coords)/1000), function(i) gseq.extract(seq_coords[(1+(i-1)*1000):(min(i*1000, nrow(seq_coords))),])))
names(seqs_all) <- mca@peaks$peak_name
tss <- dplyr::filter(gintervals.load('intervs.global.tss'), geneSymbol %in% rownames(mc_rna@e_gc))
nei_seqs_tss <- gintervals.neighbors(as.data.frame(seq_coords), tss, mindist = -2e+3, maxdist = 2e+3, maxneighbors = 1)



a_legc_avg_cl <- tgs_matrix_tapply(t(a_legc), km_a_legc$cluster, mean)
peak_clust_var <- which(matrixStats::rowMaxs(a_legc_avg_cl) - matrixStats::rowMins(a_legc_avg_cl) >= 1)
pltmt <- a_legc_avg_cl[peak_clust_var,cust_mc_ord_st]
ord <- hclust(tgs_dist(pltmt))$order
var_peaks_ordered <- do.call('c', lapply(peak_clust_var[ord], function(x) rownames(a_legc)[which(km_a_legc$cluster == x)]))
names(var_peaks_ordered) <- do.call('c', lapply(peak_clust_var[ord], function(x) rep(x, length(which(km_a_legc$cluster == x)))))

var_peaks_prom <- var_peaks_ordered[var_peaks_ordered %in% nei_seqs_tss$peak_name]
var_peaks_dist <- var_peaks_ordered[!(var_peaks_ordered %in% nei_seqs_tss$peak_name)]

# res_clusters <- prego::regress_pwm.clusters(sequences = seqs_all[var_peaks_ordered], clusters = names(var_peaks_ordered))
# save(res_clusters, file='./output/sequence_modeling/prego_on_mcatac_a_legc_var_peak_clusters.rda')

load('./output/sequence_modeling/prego_on_mcatac_a_legc_var_peak_clusters.rda')

inferred_motif_df <- dplyr::relocate(do.call('rbind', lapply(names(res_clusters$models), function(x) {
                            y <- res_clusters$models[[x]]$pssm; y$motif <- x; as.data.frame(y)}
                            )),motif, .before = 'pos')

vm <- gextract_pwm(intervals = misha.ext::convert_10x_peak_names_to_misha_intervals(var_peaks_ordered), dataset = inferred_motif_df)


intervs_energy <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')
ie_mat <- subset(intervs_energy, select = -c(chrom, start, end, peak_name))
rownames(ie_mat) <- peak_names(intervs_energy[,c('chrom', 'start', 'end')], tad_based =F)

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- c('Bcl11b', tfs$Symbol)
tfs <- tfs[tfs %in% rownames(mc_rna@e_gc)]
fp_tfs <- tfs[tfs %in% rownames(mc_rna@mc_fp)[rowMaxs(mc_rna@mc_fp) > 3]]
fp_tf_motif <- unlist(plyr::llply(fp_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(ie_mat), v=T, ign=T), .parallel = T))
calculate_d_in_cluster <- function(vec, cl_vec,
                #  parallel = getOption("mcatac.parallel"),
                parallel = TRUE,
                 alternative = "less"
                 ) {
    # print(head(vec))
    # print(head(cl_vec))
    return(plyr::laply(unique(cl_vec), function(ci) ks.test(vec[cl_vec == ci], vec[cl_vec != ci], alternative = alternative)$statistic, .parallel = parallel))
}
iemf <- ie_mat[var_peaks_ordered,]

res <- plyr::laply(colnames(iemf), function(x) calculate_d_in_cluster(iemf[,x], names(var_peaks_ordered)), .parallel = TRUE)
rownames(res) <- colnames(iemf)
save(res, file = './output/sequence_modeling/a_legc_var_peak_clusters_all_motif_ks.rda')
load('./output/sequence_modeling/a_legc_var_peak_clusters_all_motif_ks.rda')

colnames(res) <- unique(names(var_peaks_ordered))
matf <- res[unique(apply(res, 2, which.max)),]


raq98 <- apply(iemf, 2, quantile, probs = 0.98)
ra_98_bin <- sapply(1:length(raq98), function(i) as.numeric(iemf[,i] >= raq98[[i]]))
colnames(ra_98_bin) <- colnames(iemf)
rownames(ra_98_bin) <- rownames(iemf)
ra_98_sum_clust <- tgs_matrix_tapply(t(ra_98_bin), names(var_peaks_ordered), sum)
ra_98_lfc <- log2(0.1+ra_98_sum_clust/(0.02*as.numeric(table(names(var_peaks_ordered)))))
p_r9l <- pheatmap(ra_98_lfc, color = colorRampPalette(c('blue4', 'white', 'red4'))(100), show_colnames = F,
                breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)))
inds <- which(colMaxs(ra_98_lfc) - colMins(ra_98_lfc) >= 4)
cor_mot <- tgs_cor(ra_98_lfc[,inds], spearman = T)
hccm <- hclust(tgs_dist(cor_mot))
ctmot <- setNames(cutree(hccm, k = 20), colnames(ra_98_lfc)[inds])

mot_max_ks <- colMaxs(ra_98_lfc) - colMins(ra_98_lfc)
names(mot_max_ks) <- colnames(iemf)
motifs_to_take <- sapply(sort(unique(ctmot)), function(k) {
                k_inds <- which(ctmot == k); 
                ctmk <- ctmot[k_inds]; 
                return(ctmk[which.max(mot_max_ks[names(ctmk)])])
                })

pk_col_annot_atac <- tibble::column_to_rownames(as.data.frame(tibble::enframe(var_peaks_ordered, name = 'cluster_atac', value = 'peak_name')), 'peak_name')



ann_colors[['cluster_atac']] <- setNames(sample(grep('white|gray|grey', colors(), inv=T, v=T), length(unique(names(var_peaks_ordered)))), unique(names(var_peaks_ordered)))

# ra_98_sum_clust <- tgs_matrix_tapply(t(ra_98_bin[var_peaks_prom,]), names(var_peaks_prom), sum)
# ra_98_lfc <- log2(0.1+ra_98_sum_clust/(0.02*as.numeric(table(names(var_peaks_prom)))))
pltmt <- t(ra_98_lfc[order(match(rownames(ra_98_lfc), names(var_peaks_ordered))),names(motifs_to_take)])


p_plt <- pheatmap(pltmt, cluster_rows = T,cluster_cols = T, silent = F,
                color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)))
svat <- table(names(var_peaks_prom))
pltmt_rep <- do.call('cbind', lapply(colnames(pltmt), function(i) matrix(rep(pltmt[,i], svat[[i]]), nrow = nrow(pltmt))))
rownames(pltmt_rep) <- rownames(pltmt)
colnames(pltmt_rep) <- var_peaks_prom
pprep_prom <- pheatmap(pltmt_rep[rev(p_plt$tree_row$order),], 
                annotation_col = pk_col_annot_atac, 
                annotation_colors = ann_colors,
                color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                breaks = c(seq(min(pltmt_rep),0,l=50), seq(0.01,max(pltmt_rep),l=51)), 
                show_colnames = F,
                cluster_rows = F,
                cluster_cols = F,
                fontsize_row = 8)
save_pheatmap(pprep_prom, './output/mcatac/figs/mmcortex_fold_change_of_motifs_across_promoters_in_var_peak_clusters.png', 
                w = 3600, h = 2400, res = 200)

# ra_98_sum_clust <- tgs_matrix_tapply(t(ra_98_bin[var_peaks_dist,]), names(var_peaks_dist), sum)
# ra_98_lfc <- log2(0.1+ra_98_sum_clust/(0.02*as.numeric(table(names(var_peaks_dist)))))
pltmt <- t(ra_98_lfc[order(match(rownames(ra_98_lfc), names(var_peaks_dist))),names(motifs_to_take)])
p_plt <- pheatmap(pltmt, cluster_rows = T,cluster_cols = F, silent = F)
ch_plt <- Heatmap(pltmt, cluster_rows = F, cluster_columns = F)
dev.new(); draw(ch_plt)


svat <- table(names(var_peaks_dist))
pltmt_rep <- do.call('cbind', lapply(colnames(pltmt), function(i) matrix(rep(pltmt[,i], svat[[i]]), nrow = nrow(pltmt))))
rownames(pltmt_rep) <- rownames(pltmt)
colnames(pltmt_rep) <- var_peaks_dist
pprep_dist <- pheatmap(pltmt_rep[rev(p_plt$tree_row$order),], 
                annotation_col = pk_col_annot_atac, 
                annotation_colors = ann_colors,
                color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                breaks = c(seq(min(pltmt_rep),0,l=50), seq(0.01,max(pltmt_rep),l=51)), 
                show_colnames = F,
                cluster_rows = F,
                cluster_cols = F,
                fontsize_row = 8)

save_pheatmap(pprep_dist, './output/mcatac/figs/mmcortex_fold_change_of_motifs_across_distal_peaks_in_var_peak_clusters.png', 
                w = 3600, h = 2400, res = 200)

motif_plots <- lapply(names(motifs_to_take)[rev(p_plt$tree_row$order)], 
                        function(nm) plot_pssm_logo(pssm=dplyr::filter(amd, motif == nm), 
                            title = nm))
gall3 <- marrangeGrob(grobs = motif_plots, ncol =1, nrow = length(motif_plots))
ggsave('./output/sequence_modeling/figs/log_fold_q98_enriched_motifs.png', gall3, height = 40, width = 5, units = 'in')



vm <- gextract_pwm(intervals = misha.ext::convert_10x_peak_names_to_misha_intervals(var_peaks_ordered), 
                    motifs = names(motifs_to_take))
vmm <- as.matrix(subset(vm, select = -c(chrom, start, end, peak_name)))
vmm <- apply(vmm, 2, function(x) ifelse(x < -1e6, min(x[x > -1e6]), x))
vmmn <- t(vmm)
vmmn <- vmmn  - rowMaxs(vmmn)

p_acc_prom <- pheatmap(t(a_legc[var_peaks_prom,cust_mc_ord_st]), 
            annotation_row = col_annot, annotation_col = pk_col_annot_atac, 
            annotation_colors = ann_colors, annotation_legend = F,
            cluster_rows = F, cluster_cols = F, 
            show_colnames =F, show_rownames = F)
save_pheatmap(p_acc_prom, './output/mcatac/figs/mmcortex_a_legc_var_peak_promoters_annot_peaks.png', 
                w = 3600, h = 2400, res = 200)

p_acc_dist <- pheatmap(t(a_legc[var_peaks_dist,cust_mc_ord_st]), 
            annotation_row = col_annot, annotation_col = pk_col_annot_atac, 
            annotation_colors = ann_colors, annotation_legend = F,
            cluster_rows = F, cluster_cols = F, 
            show_colnames =F, show_rownames = F)
save_pheatmap(p_acc_dist, './output/mcatac/figs/mmcortex_a_legc_var_distal_peak_clusters_annot_peaks.png', 
                w = 3600, h = 2400, res = 200)

alac_z <- t(apply(a_legc_avg_cl, 1, function(x) (x - mean(x))/sd(x)))
peak_clust_ct <- peak_clust_var[apply(alac_z[as.character(peak_clust_var),], 2, which.max)]
peak_clust_gr_louv <- lapply(gr_louv, function(x) setNames(x, peak_clust_ct[as.numeric(x)]))
peak_clust_gr_louv_max <- lapply(peak_clust_gr_louv, function(x) head(names(sort(table(names(x)), decreasing=T)), 1))
peak_clust_gr_louv_motif <- lapply(peak_clust_gr_louv_max, function(x) names(which.max(matf[,x])))




cor_ra <- tgs_cor(t(res_all), spearman=T)
k_c_ra <- tglkmeans::TGL_kmeans(cor_ra, k = 60, seed = 1337)

best_motif_per_cluster <- unlist(lapply(sort(unique(k_c_ra$cluster)), 
                            function(cl) {y <- res_all[which(k_c_ra$cluster == cl),]; 
                                    return(rownames(y)[which.max(rowMaxs(y))])}))

vm <- gextract_pwm(intervals = misha.ext::convert_10x_peak_names_to_misha_intervals(var_peaks_ordered), 
                    motif = best_motif_per_cluster)
vmm <- as.matrix(subset(vm, select = -c(chrom, start, end, peak_name)))
vmm <- apply(vmm, 2, function(x) ifelse(x < -1e6, min(x[x > -1e6]), x))
vmmn <- t(vmm)
vmmn <- vmmn  - rowMaxs(vmmn)


png('./output/sequence_modeling/figs/mcatac_legc_var_peak_clusters_prego_boxplots.png', width = 3000, height = 2400, res = 150)
par(mfrow = c(5,6),mar = c(4,2,4,1), cex.axis = 2, cex.main = 3)
purrr::walk(as.character(sort(as.numeric(colnames(res_clusters$pred_mat)))), 
                    function(i) boxplot(acc ~ mot, data = as.data.frame(list(acc = res_clusters$pred_mat[,i], 
                                                                        mot = ifelse(names(var_peaks_ordered) == i, 1, 0))), 
                                                                        main = paste0('cl ', i, ' ks_d = ', 
                                                                            round(res_clusters$stats$ks_D[res_clusters$stats$cluster == i], 2))))
dev.off()

motif_plots <- lapply(names(res_clusters$models), 
                        function(nm) plot_pssm_logo(pssm=res_clusters$models[[nm]]$pssm, 
                        title = paste0('cluster ', nm)))
gall3 <- marrangeGrob(grobs = motif_plots, ncol =6, nrow = 5)
ggsave('./output/sequence_modeling/figs/mcatac_legc_var_peak_clusters_prego_motifs.png', gall3, height = 25, width = 30, units = 'in')

inf_pssm_df <- dplyr::relocate(do.call('rbind', lapply(names(res_clusters$models), function(x) {
                            y <- res_clusters$models[[x]]$pssm; y$motif <- x; as.data.frame(y)}
                            )), motif, .before = 'pos')

vm <- gextract_pwm(intervals = misha.ext::convert_10x_peak_names_to_misha_intervals(var_peaks_ordered), dataset = inf_pssm_df)
vmm <- subset(vm, select = -c(chrom, start, end, peak_name))

scdb_init(file.path(wd, 'scdb'))
mg <- scdb_mgraph('pl_cort')
mg <- mg@mgraph
library(igraph)
mg_ig <- igraph::graph_from_data_frame(mg, directed = F) 
mg_louv <- igraph::cluster_louvain(graph = mg_ig, weights = 1/mg$dist)
gr_louv <- igraph::communities(mg_louv)
ct_gr_louv <- lapply(gr_louv, function(x) setNames(x, mcmd$cell_type[as.numeric(x)]))
louv_cl_vec <- do.call('c', lapply(sort(as.numeric(names(ct_gr_louv))), 
                                    function(x) setNames(as.numeric(ct_gr_louv[[as.character(x)]]), 
                                        rep(x, length(ct_gr_louv[[as.character(x)]])))))
sample_gr <- ct_gr_louv[sample.int(length(gr_louv), 5)]

peak_clust_ct <- as.character(peak_clust_var)[apply(alac_z[as.character(peak_clust_var),], 2, which.max)]
peak_clust_gr_louv <- lapply(gr_louv, function(x) setNames(x, peak_clust_ct[as.numeric(x)]))
peak_clust_gr_louv_max <- lapply(peak_clust_gr_louv, function(x) head(names(sort(table(names(x)), decreasing=T)), 1))
peak_clust_gr_louv_motif <- lapply(peak_clust_gr_louv_max, function(x) names(which.max(matf[,x])))

amd <- all_motif_datasets()

prego_on_gr_louv <- function(gr, gr_louv, peak_names=NULL, mot = NULL) {
    if (!(gr %in% names(gr_louv))) {
        stop(glue::glue('Group {gr} not in names of gr_louv'))
    }
    mcs <- as.numeric(gr_louv[[gr]])
    r <- mca@egc[,mcs]
    if (length(mcs) > 1) {
        r <- rowSums(r)
    }
    r <- log2(1e-5 + r/sum(r))
    rn <- (r - mean(r))/sd(r)
    pssm <- dplyr::filter(amd, motif == unlist(mot))
    return(prego::regress_pwm(sequences=sa, response = rn[peak_names], 
                                        motif = pssm,
                                        match_with_db=T, 
                                        unif_prior = 0.01,
                                        # min_kmer_cor = 1e-5,
                                        # final_metric = 'ks',
                                        # use_sample = T, 
                                        parallel=T))
}
sa <- seqs_all[var_peaks_ordered]
res3 <-  plyr::llply(seq_along(peak_clust_gr_louv_motif), .f = function(gr, mot, i) {
    return(prego_on_gr_louv(gr[[i]], gr_louv=gr_louv, peak_names=var_peaks_ordered, mot = mot[[i]]))
        }, gr = names(peak_clust_gr_louv_motif), mot = peak_clust_gr_louv_motif, .parallel = T)
res_prego_on_seed_motifs <- res
save(res_prego_on_seed_motifs, file = './output/sequence_modeling/res_prego_on_seed_motifs.rda')
load('./output/sequence_modeling/res_prego_on_seed_motifs.rda')
res <- res_prego_on_seed_motifs
pred_mat <- apply(do.call('cbind',lapply(res, function(x) unlist(x['pred']))), 2, as.numeric)
pred_mat <- apply(pred_mat, 2, function(x) ifelse(x < -1e6, min(x[x>-1e6]),x))
pred_mat <- 2**pred_mat
r_mat <- apply(do.call('cbind',lapply(res, function(x) unlist(x['response']))), 2, as.numeric)
cor_pred_r <- diag(tgs_cor(as.matrix(pred_mat), as.matrix(r_mat), spearman=T))

m_pred <- colMeans(pred_mat)
sd_pred <- colSds(pred_mat)
m_r <- colMeans(r_mat)
sd_r <- colSds(r_mat)
pred_z <- t((t(pred_mat) - m_pred)/sd_pred)
r_z <- t((t(r_mat) - m_r)/sd_r)
res_z <- pred_z - r_z
res_in_pred_units <- t(t(res_z) * sd_pred + m_pred)

doMC::registerDoMC(cores = 25)
results_prego_residuals_min_kcor_001 <- plyr::llply(1:ncol(res_in_pred_units), function(i) prego::regress_pwm(sequences=sa, response = res_in_pred_units[,i], 
                                        match_with_db=T, 
                                        unif_prior = 0.01,
                                        min_kmer_cor = 1e-2,
                                        final_metric = 'r2',
                                        parallel=F), .parallel = T)
save(results_prego_residuals,file= './output/sequence_modeling/results_prego_residuals_on_louv_group_mcatac.rda')

load('./output/sequence_modeling/res_prego_on_seed_motifs.rda')
load( './output/sequence_modeling/results_prego_residuals_on_louv_group_mcatac.rda')
res1 <- res_prego_on_seed_motifs
res2 <- results_prego_residuals
r2_1_vec <- do.call('c', lapply(res1, function(x) x$r2))
r2_2_vec <- do.call('c', lapply(res2, function(x) x$r2))
sum_r2 <- r2_1_vec + r2_2_vec

make_regression_and_residal_plots <- function(inds, res1, res2) {
    motif_plots <- sapply(inds, function(i) {
                        sapply(1:4, function(j) {
                            if (j == 1) {
                                return(plot_pssm_logo(pssm=res1[[i]]$pssm, title = paste0('MC cl. ', i, ' - regressed from seed motif')))
                            } else if (j == 3) {
                                return(plot_pssm_logo(pssm=res2[[i]]$pssm, title = paste0('MC cl. ', i, ' - regressed from residuals')))
                            } else if (j == 2) {
                                return(plot_regression_prediction(pred = res1[[i]]$pred, 
                                    r = res1[[i]]$response))
                            } else if (j == 4) {
                                return(plot_regression_prediction(pred = res2[[i]]$pred, 
                                    r = res2[[i]]$response))
                            }})})
    ml <- matrix(1:length(motif_plots), ncol = 4, nrow = length(inds), byrow = T)
    grobs_arranged <- marrangeGrob(grobs = motif_plots, ncol =4, nrow = length(inds), layout_matrix = ml)
    return(grobs_arranged)
}

library(gridExtra)

ggsave('./output/sequence_modeling/figs/regression_on_residual_qc.png', gall4, height = 24, width = 16, units = 'in')



ggp <- prego::plot_regression_qc(reg = ttt[[2]])
ggsave(filename='./output/sequence_modeling/figs/prego_repressor_example.png',plot=ggp,  w = 10, h= 7)
# res <- prego::regress_pwm(sequences=seqs_all, clusters = new_clust_vec, 
#                                         use_sge = T, 
#                                         match_with_db=T, use_sample = T, 
#                                         parallel=T)
# print(warnings())
# saveRDS(res, './output/sequence_modeling/mmcortex_mcatac_rna_match_feat_peak_res.rds')