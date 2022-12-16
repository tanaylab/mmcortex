### Script to compare clusters obtained from scATAC-derived microclusters 
## and microclusters derived from scATAC-RNA metacell correlations
wd <- "/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/"
setwd(wd)

library(pheatmap)
library(metacell)
scdb_init('scdb')
mat_prom <- scdb_mat('pl_prom_cort')

devtools::load_all("~/src/mcATAC/")
library(prego)

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))
cust_mc_ord_st_day = unlist(lapply(cust_st_ord, function(s) {ord_md <- order(mcmd$mean_day[mcmd$cell_type == s]);
                                                            setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))[ord_md]}))
col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))


aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
cl_raw <- aaa$km_a_legc
cl_rna <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')

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


p_raw <- pheatmap(mat_tbl, color = colorRampPalette(c('white','red', 'black'))(100))
p_norm <- pheatmap(mat_tbl_norm_col, color = colorRampPalette(c('white','red', 'black'))(100))
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

mcl_all_sum_atac <- t(tgs_matrix_tapply(t(as.matrix(aaa$mcl_all)), cl_raw$cluster, sum))
mcl_all_sum_rna <- t(tgs_matrix_tapply(t(as.matrix(aaa$mcl_all)), cl_rna$cluster, sum))
p_sum_atac <- pheatmap(log2(mcl_all_sum_atac[p_norm$tree_col$order,]), cluster_cols = F, fontsize_row = 5, fontsize_col = 8, main = 'log2(UMI) - ATAC clusters')
p_sum_rna <- pheatmap(log2(t(mcl_all_sum_rna[p_norm$tree_row$order,])), cluster_rows = F, fontsize_row = 8, fontsize_col = 5, main = 'log2(UMI) - RNA clusters')
save_pheatmap(p_sum_atac, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_sum_atac.png', height=1600, width = 1600)
save_pheatmap(p_sum_rna, './output/sequence_modeling/figs/cluster_comparison/peak_clustering_sum_rna.png', height=1600, width = 1600)

# sc_cor_kms <- parallel::mclapply(aaa$sc_res, function(x) tglkmeans::TGL_kmeans(x$sc_cor, k = 30), mc.cores = 6)
sc_cor_kms <- readRDS('./output/mcatac/microcluster_assignment.RDS')

day_mcl_prom <- t(do.call('rbind', 
                    lapply(1:length(sc_cor_kms), 
                            function(i) {ret_mat <- tgs_matrix_tapply(as.matrix(mat_prom@mat[,colnames(aaa$sc_res[[i]]$sc_cor)]), sc_cor_kms[[i]]$cluster, sum);
                            rownames(ret_mat) <- paste0(12+i, '_', rownames(ret_mat));
                            return(ret_mat)})))

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
p_cor_mcl_mc_cl_row <- pheatmap(cor_mcl_mc[,cust_mc_ord_st], 
                            cluster_rows =T, 
                            cluster_cols = T, 
                            annotation_col = col_annot, 
                            annotation_colors = ann_colors, 
                            color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                            breaks = seq(-0.3,0.3,l=101), 
                            fontsize_row = 5, 
                            fontsize_col = 3)
p_cor_mcl_mc_cl_row_no_cl_col <- pheatmap(cor_mcl_mc[,cust_mc_ord_st_day], 
                            cluster_rows =T, 
                            cluster_cols = F, 
                            annotation_col = col_annot, 
                            annotation_colors = ann_colors, 
                            color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                            breaks = seq(-0.3,0.3,l=101), 
                            fontsize_row = 5, 
                            fontsize_col = 3)
p_cor_mcl_mc_no_cl_row <- pheatmap(cor_mcl_mc[,cust_mc_ord_st_day], 
                            cluster_rows =F, 
                            cluster_cols = F, 
                            annotation_col = col_annot, 
                            annotation_colors = ann_colors, 
                            color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                            breaks = seq(-0.3,0.3,l=101), 
                            fontsize_row = 5, 
                            fontsize_col = 3)
save_pheatmap(p_cor_mcl_mc_cl_row, './output/mcatac/figs/cor_prom_mcl_mc_sort_day.png', height=1600, width = 2600)
save_pheatmap(p_cor_mcl_mc_no_cl_row, './output/mcatac/figs/cor_prom_mcl_mc_no_sort_day.png', height=1600, width = 2600)
save_pheatmap(p_cor_mcl_mc_cl_row_no_cl_col, './output/mcatac/figs/cor_prom_mcl_mc_no_sort_col.png', height=1600, width = 2600)

cor_mcl_all <- tgs_cor(as.matrix(aaa$mcl_all), spearman = T)
p_cor_mcl_all <- pheatmap(cor_mcl_all[p_cor_mcl_mc_cl_row$tree_row$order,p_cor_mcl_mc_cl_row$tree_row$order], 
                    fontsize_row = 4, fontsize_col = 5,
                    cluster_rows = F,
                    cluster_cols = F,
                            color = colorRampPalette(c('white', 'yellow', 'red4', 'black'))(100), 
                            breaks = seq(0,1,l=101))
save_pheatmap(p_cor_mcl_all, './output/mcatac/figs/cor_mcl_all.png', height=1600, width = 1600)

p_sum_atac_sort <- pheatmap(log2(mcl_all_sum_atac[p_cor_mcl_mc_cl_row$tree_row$order,]), cluster_rows = F, fontsize_row = 5, fontsize_col = 8, main = 'log2(UMI) - ATAC clusters')
p_sum_rna_sort <- pheatmap(log2(mcl_all_sum_rna[p_cor_mcl_mc_cl_row$tree_row$order,]), cluster_rows = F, fontsize_row = 5, fontsize_col = 8, main = 'log2(UMI) - RNA clusters')
save_pheatmap(p_sum_atac_sort, './output/mcatac/figs/peak_clustering_sum_atac_sort.png', height=4600, width = 3000, res = 300)
save_pheatmap(p_sum_rna_sort, './output/mcatac/figs/peak_clustering_sum_rna_sort.png', height=4600, width = 3000, res = 300)

# p_max <- pheatmap(mat_tbl[ord_row_max,ord_col_max], color = colorRampPalette(c('white','red', 'black'))(100), cluster_rows = F, cluster_cols = F)
# p_max_norm <- pheatmap(mat_tbl_norm_col[ord_row_max_norm,ord_col_max], color = colorRampPalette(c('white','red', 'black'))(100), cluster_rows = F, cluster_cols = F)


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