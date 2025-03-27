# wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
wd <- '.'
setwd(wd)
library(grid)
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
par(las = 2, mar = c(7, 7, 2,1))
boxplot(total_frags_per_cells ~ mca@cell_to_metacell$sample_name[match(names(total_frags_per_cells), mca@cell_to_metacell$cell_id)],
               ylab = '', xlab = '')
title(xlab = 'Sample', line = 5)
title(ylab = 'Fragments per cell', line = 5)
dev.off()

pdf(fig_s3b_path, h = 500/71, w = 500/71)
par(las = 2, mar = c(7, 5, 2,1))
boxplot(frags_in_peaks_per_cells/total_frags_per_cells ~ mca@cell_to_metacell$sample_name[match(names(total_frags_per_cells), mca@cell_to_metacell$cell_id)],
               ylab = '', xlab = '')
title(xlab = 'Sample', line = 5)
title(ylab = 'Fragments in CREs per cell', line = 3)
dev.off()


pdf(fig_s3c_path, h = 400/71, w = 400/71)
plot(as.numeric(names(fc_vec)), fc_vec, ylab = 'TSS enrichment score', xlab = 'Distance from TSS')
lines(c(-1000,1000), rep(1,2), lwd = 2, lty = 2, col = 'red')
dev.off()


## Fig S3D

load(file.path(wd, 'output/mcatac/mat_feat_prom_amc.rda'))
marks <- names(scdb_gset('pl_cort_marks')@gene_set)

pltmt <- mat_feat_prom_amc - rowMeans(mat_feat_prom_amc)
hc_pltmt <- hclust(dist(pltmt), method = 'ward.D')
pp <- pheatmap::pheatmap(pltmt[hc_pltmt$order,], cluster_rows = F, 
            labels_row = ifelse(rownames(pltmt[hc_pltmt$order,]) %in% marks, rownames(pltmt[hc_pltmt$order,]) ,''), 
            show_colnames=F, treeheight_col = 0, treeheight_row = 0,fontsize = 12, clustering_method = 'ward.D',
            col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-0.5,0.5,l=101))

lr <- ifelse(rownames(pltmt) %in% marks, rownames(pltmt), '')

pp$gtable <- add.flag(pheatmap = pp, kept.labels = lr, repel.degree = 0)

save_pheatmap_pdf(pp, fig_s3d_path, h = 800/71, w = 400/71)

## Fig S3E
ord_pp <- rownames(pltmt[hc_pltmt$order,])
legc <- log2(1e-5 + mc_rna@e_gc)
pltmt2 <- legc[ord_pp,] - rowMeans(legc[ord_pp,])
pp2 <- pheatmap::pheatmap(pltmt2[,cust_mc_ord_st], cluster_rows = F, cluster_cols = F, 
    annotation_col = col_annot, annotation_colors = ann_colors, labels_row = ifelse(ord_pp %in% marks, ord_pp ,''),
    show_colnames=F, treeheight_col = 0, treeheight_row = 0,fontsize = 12, clustering_method = 'ward.D',
    col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=101))

lr <- ifelse(rownames(pltmt2) %in% marks, rownames(pltmt2), '')

pp2$gtable <- add.flag(pheatmap = pp2, kept.labels = lr, repel.degree = 0)

save_pheatmap_pdf(pp2, fig_s3e_path, h = 800/71, w = 400/71)