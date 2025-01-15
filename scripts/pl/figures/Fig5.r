
wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
library(shaman)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))


load(file = file.path(wd, 'output/methylation/avg_meth_all.rda'))

load(file = file.path(wd, 'output/hic/fig5_data.rda'))

mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]

color_key <- unique(mcmd[,c('cell_type', 'color')])
col_key <- tibble::deframe(color_key)

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))


dir.create('./output/paper_figs/Fig5/')
device <- 'pdf'
fig_5a_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/Fig5A.{device}'))
fig_5b_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/Fig5B.{device}'))
fig_5c_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/'))
fig_5d_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/Fig5D.{device}'))
fig_5e_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/Fig5E.{device}'))
fig_5f_path <- file.path(wd, glue::glue('./output/paper_figs/Fig5/Fig5F.{device}'))


## Fig 5A

p_tad_borders_all <- pheatmap::pheatmap(pltmt - rowMeans(pltmt, na.rm = T), gaps_row = nrow(hvmm_inc),
                                        cluster_rows = F, cluster_cols = F, 
                                        col = colorRampPalette(c('blue3', 'white', 'red3'))(100), 
                                        show_rownames = F, fontsize_col =  20,
                   breaks = seq(-0.35,0.35,l=100))
save_pheatmap_pdf(p_tad_borders_all, fig_5a_path, h = 800/100, w = 450/100)

## Fig 5B

pdf(fig_5b_path, h = 600/71, w = 600/71)
par(mar = c(5,5,3,1), cex.lab = 2, cex.axis = 1.5, cex.main = 2)
LWD <- 3
plot(1:length(tbl1_rm), tbl1_rm, type = 'l', col = 'purple', ylim = c(0,6), lwd = 3, main = 'Number of neighbors per region', xlab = 'Distance [kbp]', ylab = 'Rolling mean (k = 25)')
lines(c(-100,600), rep(mean(tbl1_rm, na.rm = T), 2), col = 'purple', lwd = 3, lty = 2)
lines(1:length(tbl2_rm), tbl2_rm, col = 'orange', type= 'l', lwd = 3)
lines(c(-100,600), rep(mean(tbl2_rm, na.rm = T), 2), col = 'orange', lwd = 3, lty = 2)
lines(c(-100,600), rep(0, 2), col = 'black', lwd = 1, lty = 2)
lines(rep(0, 2), c(-10,10), col = 'black', lwd = 1, lty = 2)
legend('bottomright', legend = c('Shallowing', 'Deepening', 'mean S', 'mean D'), bg = 'white', col = rep(c('purple', 'orange'), 2), lty = c(1,1,2,2), lwd = rep(LWD,4), cex = 1.5)
dev.off()

## Fig 5C -- Insulation/SHAMAN examples

clrmp_trks <- colorRampPalette(c('yellow', 'blue'))(5)

plot_insulation_and_shaman_series <- function(i, hvmi) {
    inti_nm <- paste0(hvm[hvmi,1:3], collapse = '_')
    
    tadi_fld <- fig_5c_path
    if (!dir.exists(tadi_fld)) {dir.create(tadi_fld)}
    SHIFT <- 50e+4
    gints1d_p <- dplyr::mutate(hvm[hvmi,], start = start -SHIFT, end = end + SHIFT)
    gints2d_p <- dplyr::mutate(hvm_2d[hvmi,], start1 = start1 -SHIFT, end1 = end1 + SHIFT,
                              start2 = start2 -SHIFT, end2 = end2 + SHIFT)

    nei_gints1d_ins_prc_all <- gintervals.neighbors(gints1d_p, ins_prc_all, mindist = 0, maxdist = 0, maxneighbors = 1e+7)
    inds <- nei_gints1d_ins_prc_all$intervalID
    
    pdf(glue::glue('{tadi_fld}/Fig5C{i}a.pdf'), h = 750/71, w = 1250/71)
    par(mar = c(7,10,1,1), cex.lab = 5, cex.main = 5, cex.axis = 5)
    yall <- -1 * ins_prc_all[inds,grep('E\\d\\d', colnames(ins_prc_all))]
    plot(ins_prc_all[inds,2], -1 * ins_prc_all[inds,grep('E13', colnames(ins_prc_all))], ylim = c(min(yall, na.rm = T), max(yall, na.rm = T)), bg = 'white',
         type = 'l', 
         lwd = 6,
         xaxt = 'n',
         col = clrmp_trks[[1]], pch = 16, cex = 1.5, 
         xlab = '',
         main = '',
         ylab = ''
         )
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E14', colnames(ins_prc_all))], col = clrmp_trks[[2]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E15', colnames(ins_prc_all))], col = clrmp_trks[[3]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E16', colnames(ins_prc_all))], col = clrmp_trks[[4]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E17', colnames(ins_prc_all))], col = clrmp_trks[[5]], pch = 16, lwd = 6)
    title(xlab = unique(as.character(nei_gints1d_ins_prc_all$chrom)), line = 6)
    title(ylab = 'Insulation', line = 5)
    axis(1, padj = 0.5, hadj = -0.25)
    legend('topleft', cex = 4, col = clrmp_trks, legend = paste0('E', 13:17), lwd = 6)

    dev.off()

    ps13 = gextract(score_tracks[[1]], intervals = gints2d_p, colnames="score")
    ps14 = gextract(score_tracks[[2]], intervals = gints2d_p, colnames="score")
    ps15 = gextract(score_tracks[[3]], intervals = gints2d_p, colnames="score")
    ps16 = gextract(score_tracks[[4]], intervals = gints2d_p, colnames="score")
    ps17 = gextract(score_tracks[[5]], intervals = gints2d_p, colnames="score")

    pdf(glue::glue('{tadi_fld}/Fig5C{i}b.pdf'), h = 800/71, w = 1000/71)
    shaman_plot_map_score_with_annotations(genome = 'mm10',add_genes = F, points_score = ps13, point_size = 2, interval_range = gints1d_p, add_ideogram = F)
    dev.off()

    pdf(glue::glue('{tadi_fld}/Fig5C{i}c.pdf'), h = 800/71, w = 1000/71)
    shaman_plot_map_score_with_annotations(genome = 'mm10',add_genes = F, points_score = ps14, point_size = 2, interval_range = gints1d_p, add_ideogram = F)
    dev.off()

    pdf(glue::glue('{tadi_fld}/Fig5C{i}d.pdf'), h = 800/71, w = 1000/71)
    shaman_plot_map_score_with_annotations(genome = 'mm10',add_genes = F, points_score = ps15, point_size = 2, interval_range = gints1d_p, add_ideogram = F)
    dev.off()

    pdf(glue::glue('{tadi_fld}/Fig5C{i}e.pdf'), h = 800/71, w = 1000/71)
    shaman_plot_map_score_with_annotations(genome = 'mm10',add_genes = F, points_score = ps16, point_size = 2, interval_range = gints1d_p, add_ideogram = F)
    dev.off()

    pdf(glue::glue('{tadi_fld}/Fig5C{i}f.pdf'), h = 800/71, w = 1000/71)
    shaman_plot_map_score_with_annotations(genome = 'mm10',add_genes = F, points_score = ps17, point_size = 2, interval_range = gints1d_p, add_ideogram = F)
    dev.off()
}

hvm_ls <- c('86', '214')

## Fig 5C1 -- Hapln1
plot_insulation_and_shaman_series(1, hvm_ls[[1]])

## Fig 5C2 -- Tnc
plot_insulation_and_shaman_series(2, hvm_ls[[2]])


## Fig 5D


# png('./output/hic/figs/ct_peaks_near_astro_nsc_tsss_new.png', h = 450, w = 1350, res = 100)
pdf(fig_5d_path, h = 450/71, w = 1350/71)
par(mfrow = c(1,3), cex.axis = 2, cex.lab = 2, mar = c(5,5,4,1), cex.main = 2)
bins <- seq(-0.5e+6, 0.5e+6, l = 500)
plot(0, xlim = c(-4.5e+5, max(bins)), col = 'white', xaxt = 'n',
     ylim = c(-0.6,1), 
     main = 'type-specific peak density\nnear Astrocyte TSSs', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from TSS')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'astro_tss', 'astro_peaks', 'Astrocytes', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'astro_tss', 'nsc_peaks', 'NSC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'astro_tss', 'ipc_peaks','IPC', bins = bins)
legend('topleft', legend = c(paste0('IPC peaks\nn = ', nrow(ipc_peaks_n), '\n'), paste0('NSC peaks\nn = ', nrow(nsc_peaks_n), '\n'), 
                              paste0('Astro peaks\nn = ', nrow(astro_peaks_n), '\n')), 
                              col = col_key[c('IPC', 'NSC', 'Astrocytes')], lwd = 2, lty = 1, cex = 1.75)
text(-4e+5, 2.83, labels = paste0('n_TSSs = ', nrow(astro_tss)), adj = c(0,0), cex = 1.5)
axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))


plot(0, xlim = c(-4.5e+5, max(bins)), col = 'white', xaxt = 'n',
     ylim = c(-0.6,1), 
     main = 'type-specific peak density\nnear NSC TSSs', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from TSS')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'nsc_tss', 'astro_peaks', 'Astrocytes', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'nsc_tss', 'nsc_peaks', 'NSC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'nsc_tss', 'ipc_peaks','IPC', bins = bins)

axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))
text(-4e+5, 2.3, labels = paste0('n_TSSs = ', nrow(nsc_tss)), adj = c(0,0), cex = 1.5)

plot(0, xlim = c(-4.5e+5, max(bins)), col = 'white', xaxt = 'n',
     ylim = c(-0.6,1),
     main = 'type-specific peak density\nnear IPC TSSs', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from TSS')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'ipc_tss', 'astro_peaks', 'Astrocytes', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'ipc_tss', 'nsc_peaks', 'NSC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_all, 'ipc_tss', 'ipc_peaks', 'IPC', bins = bins)
axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))
text(-4e+5, 2.3, labels = paste0('n_TSSs = ', nrow(ipc_tss)), adj = c(0,0), cex = 1.5)
dev.off()


## Fig 5E

pdf(fig_5e_path, h = 450/71, w = 1350/71)
par(mfrow = c(1,3), cex.axis = 2, cex.lab = 2, mar = c(5,5,4,1), cex.main = 2)
bins <- seq(-0.5e+6, 0.5e+6, l = 500)
plot(0, xlim = c(-5.5e+5, max(bins)), col = 'white', , xaxt = 'n',
     ylim = c(-0.6,1), 
     main = 'type-specific TSS density\nnear Astrocyte peaks', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from peak')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_astro, 'ipc_tss', 'astro_peaks', 'IPC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_astro, 'nsc_tss', 'astro_peaks', 'NSC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_astro, 'astro_tss', 'astro_peaks', 'Astrocytes', bins = bins)
legend('topleft', legend = c(paste0('IPC TSSs\nn = ', nrow(ipc_tss), '\n'), paste0('NSC TSSs\nn = ', nrow(nsc_tss), '\n'), 
                             paste0('Astro TSSs\nn = ', nrow(astro_tss), '\n')),
       col = col_key[c('IPC', 'NSC', 'Astrocytes')], lwd = 2, lty = 1, cex = 1.75)
axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))
text(-4e+5, 2.83, labels = paste0('n_peaks = ', nrow(astro_peaks_n)), adj = c(0,0), cex = 1.5)


plot(0, xlim = c(min(bins), max(bins)), col = 'white', , xaxt = 'n',
     ylim = c(-0.6,1),
     main = 'type-specific TSS density\nnear NSC peaks', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from peak')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_nsc, 'ipc_tss', 'nsc_peaks', 'IPC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_nsc, 'astro_tss', 'nsc_peaks', 'Astrocytes', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_nsc, 'nsc_tss', 'nsc_peaks', 'NSC', bins = bins)

axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))
text(-4e+5, 2.3, labels = paste0('n_peaks = ', nrow(nsc_peaks_n)), adj = c(0,0), cex = 1.5)


plot(0, xlim = c(min(bins), max(bins)), col = 'white', , xaxt = 'n',
     ylim = c(-0.6,1), 
     main = 'type-specific TSS density\nnear IPC peaks', ylab = 'Smoothed log-fold enrichment', xlab = 'bp from peak')
grid()
plot_ct_nei_distance_distribution(nei_peaks_tss_ipc, 'ipc_tss', 'ipc_peaks', 'IPC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_ipc, 'nsc_tss', 'ipc_peaks','NSC', bins = bins)
plot_ct_nei_distance_distribution(nei_peaks_tss_ipc, 'astro_tss', 'ipc_peaks','Astrocytes', bins = bins)

axis(1, at = c(-4e+5, 0, 4e+5), labels = c('-400kbp', '0', '400kbp'))
text(-4e+5, 2.3, labels = paste0('n_peaks = ', nrow(ipc_peaks_n)), adj = c(0,0), cex = 1.5)
dev.off()


## Fig 5F

pdf(fig_5f_path, h =400/71, w = 500/71)

par(cex.axis = 1.5, cex.lab = 2, mar = c(12,7,2,1), las = 2)
boxplot(dist ~ interaction, data = nei_ct_peaks_ct_peaks_1M, 
       ylab = '', xlab = '')
title(ylab = 'Distance [bp]', line = 5)
title(xlab = 'Peak type pairs', line = 9)
ttt <- sapply(seq(4.5,4.5+2*4,4), function(x) lines(rep(x,2), c(0,1e+7), lwd = 3))
dev.off()

