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

tss <- gintervals.load('intervs.global.tss')
tss <- tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]
tss2 <- tss[!duplicated(tss$geneSymbol),]
tads <- gintervals.load('intervs.global.tad_names')
nsc_mcs <- mcmd$metacell[mcmd$cell_type == 'NSC']
ipc_mcs <- mcmd$metacell[mcmd$cell_type %in%  c('IPC', 'IPC_cyc')]



dir.create(file.path(wd, 'output/paper_figs/Fig3/'))
device <- 'pdf'
fig_3a1_path <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/Fig3A1.{device}'))
fig_3a2_path <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/Fig3A2.{device}'))
fig_3b_path <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/Fig3B.{device}'))
fig_3cdef_path <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/Fig3CDEF.{device}'))
fig_cfupn_traj <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/cfupn_traj.{device}'))
fig_cpn_traj <- glue::glue(file.path(wd, 'output/paper_figs/Fig3/cpn_traj.{device}'))

load(file.path(wd, 'output/mcatac/fig3_data.rda'))
load(file.path(wd, 'output/metacell_model/diff_order_data.rda'))

# New version - 25/04/24

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

inds_glia <- which(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes'))
inds_no_glia <- which(!(names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes')))


pltmt <- a_legc_avg_cl_enh[enh_cl_ord, cust_mc_ord_st2]

brks <- seq(-16.7,-14.75,l=100)

col_ha1 <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_no_glia],'cell_type'], 
                                                     col = ann_colors[['cell_type']],
                                                     height =unit(2, 'cm')), 
                                                     mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_no_glia],'mean_day'], 
                                                        axis_param = list(gp = gpar(fontsize = 0)),
                                                        height =unit(2, 'cm')), 
                                                     annotation_name_gp = gpar(fontsize = 18),
                                                     show_legend = F)

row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_enh_a_legc$size[as.numeric(enh_cl_ord)]), 
                                                ylim = c(800,8000), gp = gpar(fill = 'black',fontsize = 18), 
                                                axis_param = list(facing = 'inside', 
                                                gp = gpar(fontsize = 20), labels_rot = -90)),
                                            annotation_name_gp = gpar(fontsize = 18),
                                            annotation_name_rot = 0,
                                            annotation_name_offset = unit(3, 'cm'),        
                                            which = 'row',
                                            width = unit(5, 'cm')
                           )
clrmp_rel2 <- circlize::colorRamp2(breaks = c(-3,0,3), colors = c('blue3','white','red3'))
motifs_anno_mat <- round(ra_98_lfc_int[enh_cl_ord,mtt], 3)

hc_mtt <- hclust(dist(t(motifs_anno_mat)), method = 'ward.D2')
colnames(motifs_anno_mat) <- unlist(purrr::map(stringr::str_split(colnames(motifs_anno_mat), '\\.'), 2))
colnames(motifs_anno_mat) <- unlist(purrr::map(stringr::str_split(colnames(motifs_anno_mat), '_'), 1))
mam_lin <- apply(motifs_anno_mat, 2, function(x) {y <- x; y[x < -3] <- -3; y[x > 3] <- 3; return((y + 3)/6)})

motif_ha <- rowAnnotation(`motif\nenrichment` = motifs_anno_mat[,hc_mtt$order], col = list(`motif\nenrichment` = clrmp_rel2),
                            width = unit(9, 'cm')
                             )



ch <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_no_glia], 
                              name = 'log2\nfraction\nATAC',
                              col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), 
                                colors = c('white', 'orange', 'red', 'purple', 'black')),
                              column_split = factor(names(cust_mc_ord_st2[inds_no_glia]), levels = cust_st_ord2), 
                              column_gap = unit(2, 'mm'),
                              column_title_gp = gpar(fontsize = 20),
                              column_title_rot = 90,
                              row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), 
                              row_gap = unit(2, 'mm'),
                              row_title_gp = gpar(fontsize = 0),
                              top_annotation = col_ha1, 
                              show_heatmap_legend = T,
                              show_column_names = F,
                              right_annotation = motif_ha,
                              heatmap_legend_param = list(legend_height = unit(5, 'in'), 
                              legend_width = unit(5, 'in'), 
                              labels_gp = gpar(fontsize = 16)),
                              row_names_gp = gpar(fontsize = 14),
                              heatmap_width = unit(45, 'cm'), 
                              heatmap_height = unit(32, 'cm'),
                              cluster_columns = F, 
                              cluster_rows = F)


col_ha1_glia <- HeatmapAnnotation(cell_type = anno_simple(col_annot[cust_mc_ord_st2[inds_glia],'cell_type'], 
                                                     col = ann_colors[['cell_type']],
                                                     height =unit(2, 'cm')), 
                             mean_day = anno_lines(x = col_annot[cust_mc_ord_st2[inds_glia],'mean_day'], axis_param = list(gp = gpar(fontsize = 16)), ylim = c(13,18),
                                                    height =unit(2, 'cm')), 
                             annotation_name_gp = gpar(fontsize = 0),
                             show_legend = F)


ch_glia <- ComplexHeatmap::Heatmap(matrix = pltmt[,inds_glia], 
                              col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
                                   column_split = factor(names(cust_mc_ord_st2[inds_glia]), levels = c('OPCs', 'Astrocytes')), column_gap = unit(3, 'mm'),
                              column_title_gp = gpar(fontsize = 20),
                              column_title_rot = 90,
                              row_split = factor(ifelse(rownames(pltmt) %in% rownames(m1), 1, 2), levels = c(1,2)), 
                                   row_gap = unit(2, 'mm'),
                              row_title_gp = gpar(fontsize = 0),
                              top_annotation = col_ha1_glia, 
                                   left_annotation = row_ha,
                              show_heatmap_legend = F,
                              show_column_names = F,
                                   show_row_names = F,
                              heatmap_width = unit(10, 'cm'), 
                                   heatmap_height = unit(32, 'cm'),
                              row_names_gp = gpar(fontsize = 14),
                        cluster_columns = F, cluster_rows = F)


pdf(fig_3a1_path, w = 1480/.75e+2, h =1170/.71e+2)
draw(ch)
dev.off()

pdf(fig_3a2_path, w = 320/.71e+2, h =1270/.71e+2)
draw(ch_glia)
dev.off()


## Fig 3B

pltmt2 <- a_legc_avg_cl_prom[,cust_mc_ord_st2]

col_ha1 <- HeatmapAnnotation(cell_type = anno_simple(col_annot$cell_type[match(cust_mc_ord_st2, mcmd$metacell)], 
                                                     col = ann_colors[['cell_type']],
                                                     height =unit(1, 'cm')), 
                             mean_day = anno_lines(x = col_annot$mean_day[match(cust_mc_ord_st2, mcmd$metacell)], 
                                                    axis_param = list(gp = gpar(fontsize = 16)),
                                                    height =unit(1, 'cm')), 
                             annotation_name_gp = gpar(fontsize = 18),
                             show_legend = F)

row_ha <- HeatmapAnnotation(cluster_size = anno_barplot(as.numeric(km_prom_a_legc$size), ylim = c(800,3000), 
                                                        gp = gpar(fill = 'black',fontsize = 18), 
                                                        axis_param = list(facing = 'inside', 
                                                        gp = gpar(fontsize = 20), 
                                                        labels_rot = -90)), 
                            annotation_name_gp = gpar(fontsize = 18),
                            annotation_name_rot = 0,
                            annotation_name_offset = unit(3, 'cm'),        
                            which = 'row',
                            width = unit(5, 'cm')
                           )


ch2 <- ComplexHeatmap::Heatmap(matrix = pltmt2[,], name = 'log2\nfraction\nATAC', 
                              col = circlize::colorRamp2(breaks =  seq(-16.6,-14,l=5), colors = c('white', 'orange', 'red', 'purple', 'black')),
                              column_split = factor(names(cust_mc_ord_st2), levels = cust_st_ord2), column_gap = unit(2, 'mm'),
                              column_title_gp = gpar(fontsize = 0),
                              top_annotation = col_ha1, 
                              show_heatmap_legend = T,
                              show_column_names = F,
                              show_row_dend = FALSE,
                              column_title_rot = 90,
                              heatmap_legend_param = list(legend_height = unit(5, 'in'), 
                              legend_width = unit(5, 'in'), 
                              labels_gp = gpar(fontsize = 16)),
                              heatmap_width = unit(45, 'cm'), 
                              heatmap_height = unit(8, 'cm'),
                              left_annotation = row_ha,
                              row_names_gp = gpar(fontsize = 16),
                              cluster_columns = F, cluster_rows = T,
                              clustering_method_rows = 'ward.D2')


pdf(fig_3b_path, w = 1450/71, h = 350/71)
draw(ch2)
dev.off()


## Fig 3CDEF - CthPN-specific vs CPN-specific peak trajectories

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

pdf(fig_3cdef_path, h = 600/71, w = 600/71)

MAR_LL <- 4
MAR_LR <- 1

layout(mat = mat_layout)

par(mar = c(0.5,MAR_LL,3,1), cex.lab = 1.52, cex.main = 1.5)

ctsi <- c('CthPN', 'SCPN', 'CPN_L2-3', 'CPN_L5_6', 'IPC')
boxplot(a_legc_avg_ct[pks001,ctsi], col = col_key[ctsi], ylab = '',
        main = glue::glue('CfuPN-specific peaks, n = {length(pks001)}'), xaxt = 'n', yaxt = 's', ylim = c(-16.6, -14))
title(ylab = 'Mean accessibility', line = 2)
par(mar = c(0.5,MAR_LR,3,1))
boxplot(a_legc_avg_ct[pks010,ctsi], col = col_key[ctsi], 
                    main = glue::glue('CPN-specific peaks, n = {length(pks010)}'), 
                    xaxt = 'n', yaxt = 'n', ylim = c(-16.6, -14), ylab  = '')

par(mar = c(1.5,MAR_LL,1.5,1))
barplot(t(tbl_vnn_ct_cfupn_norm), col = col_key[colnames(tbl_vnn_ct_cfupn_norm)], xaxt = 'n', yaxt = 'n', bty = 'n')
par(mar = c(1.5,MAR_LR,1.5,1))
barplot(t(tbl_vnn_ct_cfupn_norm), col = col_key[colnames(tbl_vnn_ct_cfupn_norm)], 
               xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '')

mcsh <- as.character(mcmd$metacell[mcmd$cell_type %in% colnames(tbl_vnn_ct_cfupn_norm)])
par(mar = c(0.5,MAR_LL,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks001,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls, 
               xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid()
par(mar = c(0.5,MAR_LR,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks010,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,
               xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid(ylab = '')

mcsh <- as.character(mcmd$metacell[mcmd$cell_type %in% colnames(tbl_vnn_ct_cpn_norm)])
par(mar = c(0.5,MAR_LL,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks001,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,
               xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid()
par(mar = c(0.5,MAR_LR,0.5,1))
boxplot(Matrix::colMeans(mcacp@egc[pks010,mcsh]) ~ bins_nsc_neu[mcsh], col = diff_bin_clvls,
               xaxt = 'n', ylab = '', ylim = c(0,1), yaxt = 'n')
cust_grid(ylab = '')

par(mar = c(4,MAR_LL,0.5,1))

barplot(t(tbl_vnn_ct_cpn_norm), col = col_key[colnames(tbl_vnn_ct_cpn_norm)], xaxt = 's', yaxt = 'n', bty = 'n')
title(xlab = 'Differentiation bin', line = 3)
par(mar = c(4,MAR_LR,0.5,1))
barplot(t(tbl_vnn_ct_cpn_norm), col = col_key[colnames(tbl_vnn_ct_cpn_norm)], xaxt = 's', yaxt = 'n', bty = 'n')
title(xlab = 'Differentiation bin', line = 3)
dev.off()


mc2d <- scdb_mc2d('pl_cort_not_cor_cc')
CEX <- 2
LWD <- 4

pdf(fig_cfupn_traj, h = 500/71, w = 500/71)

plot(mc2d@mc_x, mc2d@mc_y, pch = 1, cex = CEX, col = 'black', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
points(mc2d@mc_x, mc2d@mc_y, pch = 16, cex = CEX, col = mcmd$color, lwd = LWD)
lines(pcu_cfupn$s[pcu_cfupn$ord,1], pcu_cfupn$s[pcu_cfupn$ord,2], lwd = 4, col = )
arrows(x0 = pcu_cfupn$s[pcu_cfupn$ord[[10]],1], y0 = pcu_cfupn$s[pcu_cfupn$ord[[10]],2], 
       x1 = pcu_cfupn$s[pcu_cfupn$ord[[1]],1], y1 = pcu_cfupn$s[pcu_cfupn$ord[[1]],2], lwd = 4)
dev.off()

pdf(fig_cpn_traj, h = 500/71, w = 500/71)
plot(mc2d@mc_x, mc2d@mc_y, pch = 1, cex = CEX, col = 'black', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
points(mc2d@mc_x, mc2d@mc_y, pch = 16, cex = CEX, col = mcmd$color, lwd = LWD)
lines(pcu_cpn$s[pcu_cpn$ord,1], pcu_cpn$s[pcu_cpn$ord,2], lwd = 4, col = 'black')
arrows(x0 = pcu_cpn$s[pcu_cpn$ord[[length(pcu_cpn$ord) - 50]],1], y0 = pcu_cpn$s[pcu_cpn$ord[[length(pcu_cpn$ord) - 50]],2], 
       x1 = pcu_cpn$s[pcu_cpn$ord[[length(pcu_cpn$ord)]],1], y1 = pcu_cpn$s[pcu_cpn$ord[[length(pcu_cpn$ord)]],2], lwd = 4)
dev.off()