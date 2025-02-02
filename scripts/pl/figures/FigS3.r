wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))
library(metacell)
scdb_init('scdb',f=T)

mc_rna <- scdb_mc('pl_cort')


mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]


cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))

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

load(file.path(wd, 'output/mcatac/fig3_data.rda'))
load(file.path(wd, 'output/metacell_model/diff_order_data.rda'))
load(file.path(wd, 'output/mcatac/fig_s3_data.rda'))
load(file = file.path(wd, 'output/mcatac/fig_s3_atac_qc_data.rda'))

dir.create(file.path(wd, 'output/paper_figs/FigS3/'))
device <- 'pdf'
fig_s3a_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3A.{device}'))
fig_s3b_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3B.{device}'))
fig_s3c_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3C.{device}'))
fig_s3d_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3D.{device}'))
fig_s3e_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3E.{device}'))
fig_s3f_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3F.{device}'))
fig_s3g_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3G.{device}'))
fig_s3h_path <- glue::glue(file.path(wd, 'output/paper_figs/FigS3/FigS3H.{device}'))

pdf(fig_s3a_path, h = 500/71, w = 500/71)
# png('./output/mcatac/figs/total_frags_per_cells_per_sample.png', h = 500, w = 500)
par(las = 2, mar = c(7, 7, 2,1))
boxplot(total_frags_per_cells ~ mca@cell_to_metacell$sample_name[match(names(total_frags_per_cells), mca@cell_to_metacell$cell_id)],
               ylab = '', xlab = '')
title(xlab = 'Sample', line = 5)
title(ylab = 'Fragments per cell', line = 5)
dev.off()

pdf(fig_s3b_path, h = 500/71, w = 500/71)
# png('./output/mcatac/figs/fragments_in_peaks_per_cell_per_sample.png', h = 500, w = 500)
par(las = 2, mar = c(7, 5, 2,1))
boxplot(frags_in_peaks_per_cells/total_frags_per_cells ~ mca@cell_to_metacell$sample_name[match(names(total_frags_per_cells), mca@cell_to_metacell$cell_id)],
               ylab = '', xlab = '')
title(xlab = 'Sample', line = 5)
title(ylab = 'Fragments in peaks per cell', line = 3)
dev.off()


pdf(fig_s3c_path, h = 400/71, w = 400/71)
# png('./output/mcatac/figs/tss_enrichment_score.png', h = 400, w = 400)
plot(as.numeric(names(fc_vec)), fc_vec, ylab = 'TSS enrichment score', xlab = 'Distance from TSS')
lines(c(-1000,1000), rep(1,2), lwd = 2, lty = 2, col = 'red')
dev.off()


## Fig S3D
# mat_feat_prom_amc <- do.call('cbind', lapply(day_mcls, function(x) x[[2]][feats,]))

load(file.path(wd, 'output/mcatac/mat_feat_prom_amc.rda'))
marks <- names(scdb_gset('pl_cort_marks')@gene_set)

pltmt <- mat_feat_prom_amc - rowMeans(mat_feat_prom_amc)
hc_pltmt <- hclust(dist(pltmt), method = 'ward.D')
pp <- pheatmap::pheatmap(pltmt[hc_pltmt$order,], cluster_rows = F, 
            labels_row = ifelse(rownames(pltmt[hc_pltmt$order,]) %in% marks, rownames(pltmt[hc_pltmt$order,]) ,''), 
            show_colnames=F, treeheight_col = 0, treeheight_row = 0,fontsize = 12, clustering_method = 'ward.D',
            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-0.5,0.5,l=101))
# save_pheatmap_png(pp, './output/mcatac/figs/atac_cortex_matching_promoter_activity.png', h = 3200, w = 1600)
save_pheatmap_pdf(pp, fig_s3d_path, h = 800/71, w = 400/71)

## Fig S3E
ord_pp <- rownames(pltmt[hc_pltmt$order,])
legc <- log2(1e-5 + mc_rna@e_gc)
pltmt2 <- legc[ord_pp,] - rowMeans(legc[ord_pp,])
pp2 <- pheatmap::pheatmap(pltmt2[,cust_mc_ord_st], cluster_rows = F, cluster_cols = F, 
    annotation_col = col_annot, annotation_colors = ann_colors, labels_row = ifelse(ord_pp %in% marks, ord_pp ,''),
    show_colnames=F, treeheight_col = 0, treeheight_row = 0,fontsize = 12, clustering_method = 'ward.D',
    col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=101))
# save_pheatmap_png(pp2, './output/mcatac/figs/atac_cortex_matching_rna_feature_activity.png', h = 3200, w = 2000)
save_pheatmap_pdf(pp2, fig_s3e_path, h = 800/71, w = 400/71)



## Fig S3F
# png('./output/mcatac/figs/tss-proximal_peak_variance_densities.png')
pdf(fig_s3f_path)
par(cex.lab = 2, mar = c(5,5,2,1), cex.main = 2, cex.axis = 2)
BW <- 0.1
d1 <- density(prom_sds[setdiff(names(prom_sds), unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem, proms_hi_var_marks))))]**2, bw = BW)
d2 <- density(prom_sds[proms_hi_var_astro]**2, bw = BW)
d3 <- density(prom_sds[proms_hi_var_ipc]**2, bw = BW)
d4 <- density(prom_sds[proms_hi_var_stem]**2, bw = BW)
d5 <- density(prom_sds[setdiff(proms_hi_var_marks, unique(unlist(c(proms_hi_var_astro, proms_hi_var_ipc, proms_hi_var_stem))))]**2, bw = BW)
plot(d1$x, d1$y, col = 'black', type = 'l', lwd = 2, lty = 2, xlab = 'ATAC variance across metacells', ylab = glue::glue('Density, BW = {BW}'), 
     main = 'TSS-proximal peak ATAC variance'
    )
lines(d2$x, d2$y, col = col_key[['Astrocytes']], lwd = 3, lty = 2)
lines(d3$x, d3$y, col = col_key[['IPC']], lwd = 3, lty = 2)
lines(d4$x, d4$y, col = col_key[['NSC']], lwd = 3, lty = 2)
lines(d5$x, d5$y, col = 'darkorange', lwd = 3, lty = 2)
legend('topright', legend = c('All TSS-proximal peaks', 'Near astro module TSSs', 'Near IPC module TSSs', 'Near stem module TSSs', 'Near marker gene TSSs'),
       col = c('black', col_key[c('Astrocytes', 'IPC', 'NSC')], 'darkorange'), lty = rep(2,5), lwd = rep(2,5), cex = 1)
dev.off()


## Fig S3G

cust_grid <- function(ylab = 'Normalized ATAC') {
    LWD <- 2
    RNG <- 1e+6
    Y_LINE <- 0.1
    Y_LINE2 <- 0.5
    X_LINE <- 10.5
    COL_Y <- 'lightblue'
    COL_X <- 'lightblue'
    lines(RNG*c(-1,1), rep(0.1,2), col = COL_Y, lty = 2, lwd = LWD)
    lines(RNG*c(-1,1), rep(0.5,2), col = COL_Y, lty = 2, lwd = LWD)
    lines(RNG*c(-1,1), rep(0.9,2), col = COL_Y, lty = 2, lwd = LWD)
    lines(rep(X_LINE, 2), RNG*c(-1,1), col = COL_X, lty = 2, lwd = LWD)
    axis(2, at = c(0.1,0.5,0.9), labels = c(0.1,0.5,0.9))
    title(ylab = ylab, line = 2)
}

diff_bin_clvls <- clrmp_abs[round(seq(1,length(clrmp_abs), l = 20))]

vec_layout <- c(rep(c(1,2), 2), rep(c(3,4), 1), rep(c(5,6), 2), rep(c(7,8), 2), rep(c(9,10), 1))
mat_layout <- matrix(vec_layout, nrow = 8, ncol = 2, byrow = T)


pdf(fig_s3g_path, h = 600/71, w = 600/71)

MAR_LL <- 4
MAR_LR <- 1

layout(mat = mat_layout)

par(mar = c(0.5,MAR_LL,3,1), cex.lab = 1.52, cex.main = 1.5)

ctsi <- c('CthPN', 'SCPN', 'CPN_L2-3', 'CPN_L5_6', 'IPC')
boxplot(a_legc_avg_ct[pks101,ctsi], col = col_key[ctsi], 
                main = glue::glue('CfuPN+, IPC+ peaks, n = {length(pks101)}'), ylab = '',
                 xaxt = 'n', yaxt = 's', ylim = c(-16.6, -14))
title(ylab = 'Mean accessibility', line = 2)

par(mar = c(0.5,MAR_LR,3,1))
boxplot(a_legc_avg_ct[pks110,ctsi], col = col_key[ctsi], 
        main = glue::glue('CPN+, IPC+ peaks, n = {length(pks110)}'), xaxt = 'n', yaxt = 'n', ylim = c(-16.6, -14))
par(mar = c(1.5,MAR_LL,1.5,1))
barplot(t(tbl_vnn_ct_cfupn_norm), col = col_key[colnames(tbl_vnn_ct_cfupn_norm)], 
        xaxt = 'n', yaxt = 'n', bty = 'n')
par(mar = c(1.5,MAR_LR,1.5,1))
barplot(t(tbl_vnn_ct_cfupn_norm), col = col_key[colnames(tbl_vnn_ct_cfupn_norm)], xaxt = 'n', yaxt = 'n', bty = 'n')

mcsh <- as.character(mcmd$metacell[mcmd$cell_type %in% colnames(tbl_vnn_ct_cfupn_norm)])
par(mar = c(0.5,MAR_LL,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks101,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls, xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid()
par(mar = c(0.5,MAR_LR,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks110,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid(ylab = '')

mcsh <- as.character(mcmd$metacell[mcmd$cell_type %in% colnames(tbl_vnn_ct_cpn_norm)])
par(mar = c(0.5,MAR_LL,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks101,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid()
par(mar = c(0.5,MAR_LR,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks110,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid(ylab = '')

par(mar = c(4,MAR_LL,0.5,1))
barplot(t(tbl_vnn_ct_cpn_norm), col = col_key[colnames(tbl_vnn_ct_cpn_norm)], xaxt = 's', yaxt = 'n', bty = 'n')
title(xlab = 'Differentiation bin', line = 3)
par(mar = c(4,MAR_LR,0.5,1))
barplot(t(tbl_vnn_ct_cpn_norm), col = col_key[colnames(tbl_vnn_ct_cpn_norm)], xaxt = 's', yaxt = 'n', bty = 'n')
title(xlab = 'Differentiation bin', line = 3)
dev.off()




## Fig S3H


nsc_late <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day > 16.5])
nsc_early <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day < 14.5])


pdf(fig_s3h_path, h = 500/71, w = 1000/71)
# png('./output/mcatac/figs/astro_and_early_nsc_vs_late_nsc_atac.png', h = 500, w = 1000)
par(mfrow = c(1,2), cex.lab = 1.52, cex.main = 1.5)
plot(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), a_legc_avg_ct[,'Astrocytes'], pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) ATAC', ylab = 'Astrocytes ATAC', main = 'Astro vs late NSC - ATAC')
abline(a =-1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '0.5 LFC','1 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), a_legc_avg_ct[,'Astrocytes'], method = 'pearson')
text(-15.5,-13, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
plot(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_early]), pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) ATAC', ylab = 'NSC (mean day < 14.5) ATAC', main = 'Early vs late NSC - ATAC')
abline(a =-1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '0.5 LFC','1 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_early]), method = 'pearson')
text(-15.5,-13, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
dev.off()