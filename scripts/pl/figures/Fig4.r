
wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))


load(file = file.path(wd, 'output/mcatac/fig4_atac_data.rda'))
load(file = file.path(wd, 'output/methylation/fig4_meth_data.rda'))


mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]

color_key <- unique(mcmd[,c('cell_type', 'color')])
col_key <- tibble::deframe(color_key)

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))


dir.create('./output/paper_figs/Fig4/')
device <- 'pdf'
fig_4a_path <- glue::glue('./output/paper_figs/Fig4/Fig4A.{device}')
fig_4b_path <- glue::glue('./output/paper_figs/Fig4/Fig4B.{device}')
fig_4c_path <- glue::glue('./output/paper_figs/Fig4/Fig4C.{device}')
fig_4d_path <- glue::glue('./output/paper_figs/Fig4/Fig4D.{device}')
fig_4e_path <- glue::glue('./output/paper_figs/Fig4/Fig4E.{device}')
fig_4f_path <- glue::glue('./output/paper_figs/Fig4/Fig4F.{device}')
fig_4g_path <- glue::glue('./output/paper_figs/Fig4/Fig4G.{device}')
fig_4h_path <- glue::glue('./output/paper_figs/Fig4/Fig4H.{device}')


## Fig 4A

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
save_pheatmap_pdf(p_inc_atac, file.path(wd, fig_4a_path), h = 600/71, w = 500/71)

## Fig 4B
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
save_pheatmap_pdf(p_dec_atac, file.path(wd, fig_4b_path), h = 600/71, w = 500/71)


## Fig 4C

pdf(fig_4c_path, w = 600/111,h = 500/61)
par(cex.lab = 1.82, cex.axis = 2, cex.main = 2.5, las =2, cex.lab = 2)

layout(mat = matrix(c(1, 1, 2, 2, 2), nrow = 5, ncol = 1))
xaxti <- 'n'
par(mar = c(2,8,3,1))
i <- 'inc_peaks'
boxplot(pltmt[i,] ~ factor(mcmd$cell_type, levels = cust_st_ord), 
                ylim = quantile(pltmt[i,], c(0.02,1)), 
                col = col_key[cust_st_ord], yaxt = 's',xaxt = 'n', xlab  = '', ylab = '',main = 'NSC increasing peaks', horizontal = F)
title(ylab = 'log2 fraction', line = 6)
par(mar = c(16,8,2,1))
i <- 'dec_peaks'
boxplot(pltmt[i,] ~ factor(mcmd$cell_type, levels = cust_st_ord), 
                ylim = quantile(pltmt[i,], c(0.02,1)), 
                col = col_key[cust_st_ord], yaxt = 's',xaxt = 's', xlab  = '', ylab = '',main = 'NSC decreasing peaks', horizontal = F)
title(ylab = 'log2 fraction', line = 6)
dev.off()



## Fig 4D
# Schematic


## Fig 4E

ord <- order(rowMeans(avg_meth_all[nsc_inc_peaks,grep('E\\d\\d', colnames(avg_meth_all))]))
p_meth_inc_atac <- pheatmap::pheatmap(avg_meth_all[nsc_inc_peaks[ord],grep('E\\d\\d', colnames(avg_meth_all))], fontsize_col = 16,
                                      # main = 'NSC methylation of\nincreasing peaks',
                                      fontsize = 16,
                                      cluster_rows = F, 
                                      clustering_method = 'ward.D2',
                                      cluster_cols = F,
                                      show_rownames = F,
                                      treeheight_row = 0,
                                      color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
                                     breaks = seq(-1e-4,1,l = 1000))
save_pheatmap_pdf(p_meth_inc_atac, file.path(wd, fig_4e_path), h = 900/101, w = 450/71)

## Fig 4F
ord <- order(rowMeans(avg_meth_all[nsc_dec_peaks,grep('E\\d\\d', colnames(avg_meth_all))]))
p_meth_dec_atac <- pheatmap::pheatmap(avg_meth_all[nsc_dec_peaks[ord],grep('E\\d\\d', colnames(avg_meth_all))], fontsize_col = 16,
                                      # main = 'NSC methylation of\ndecreasing peaks',
                                      fontsize = 16,
                                      cluster_rows = F, 
                                      clustering_method = 'ward.D2',
                                      cluster_cols = F,
                                      show_rownames = F,
                                      treeheight_row = 0,
                                      color = colorRampPalette(c('white', 'orange', 'red3', 'purple', 'black'))(1000),
                                     breaks = seq(-1e-4,1,l = 1000))
save_pheatmap_pdf(p_meth_dec_atac, file.path(wd, fig_4f_path), h = 900/101, w = 450/71)

## Fig 4G

pdf(file.path(wd, fig_4g_path), h = 900/81, w = 350/81)
par(mfrow = c(4,1), las = 2, mar = c(7,4,2,1), cex.axis = 1.5, cex.main = 2)

t13_nsc <- table(cut(avg_meth_all[nsc_peaks,grep('E13', colnames(avg_meth_all))], breaks = meth_qs))
t17_nsc <- table(cut(avg_meth_all[nsc_peaks,grep('E17', colnames(avg_meth_all))], breaks = meth_qs))
axis_labels <- names(t13_nsc)
axis_labels[[1]] <- '[0,0.01]'

barplot(do.call('c', lapply(seq_along(t13_nsc), function(i) c(t13_nsc[[i]], t17_nsc[[i]]))), space = rep(c(2, 0.8), length(t13_nsc)), main = 'NSC peaks', 
                            col = rep(c(col_key[['NSC']], adjustcolor(col_key[['NSC']], alpha.f = 0.2)), length(t13_nsc)))
axis(1, at=seq(2,5*(length(t13_nsc)-1),4.8), labels = axis_labels)
legend('topright', col = c(col_key[['NSC']], adjustcolor(col_key[['NSC']], alpha.f = 0.2)), cex = 2, pch = 15, legend = c('E13','E17'))
                
t13_astro <- table(cut(avg_meth_all[astro_peaks,grep('E13', colnames(avg_meth_all))], breaks = meth_qs))
t17_astro <- table(cut(avg_meth_all[astro_peaks,grep('E17', colnames(avg_meth_all))], breaks = meth_qs))

barplot(do.call('c', lapply(seq_along(t13_astro), function(i) c(t13_astro[[i]], t17_astro[[i]]))), space = rep(c(2, 0.8), length(t13_astro)), main = 'Astro peaks', 
                            col = rep(c(col_key[['Astrocytes']], adjustcolor(col_key[['Astrocytes']], alpha.f = 0.2)), length(t13_nsc)))
axis(1, at=seq(2,5*(length(t13_nsc)-1),4.8), labels = axis_labels)
legend('topright', col = c(col_key[['Astrocytes']], adjustcolor(col_key[['Astrocytes']], alpha.f = 0.2)),cex = 2, pch = 15, legend = c('E13','E17'))
                
t13_ipc <- table(cut(avg_meth_all[ipc_peaks,grep('E13', colnames(avg_meth_all))], breaks = meth_qs))
t17_ipc <- table(cut(avg_meth_all[ipc_peaks,grep('E17', colnames(avg_meth_all))], breaks = meth_qs))
barplot(do.call('c', lapply(seq_along(t13_ipc), function(i) c(t13_ipc[[i]], t17_ipc[[i]]))), space = rep(c(2, 0.8), length(t13_ipc)), main = 'IPC peaks', 
                            col = rep(c(col_key[['IPC']], adjustcolor(col_key[['IPC']], alpha.f = 0.2)), length(t13_nsc)))
axis(1, at=seq(2,5*(length(t13_nsc)-1),4.8), labels = axis_labels)
legend('topright', col = c(col_key[['IPC']], adjustcolor(col_key[['IPC']], alpha.f = 0.2)),cex = 2, pch = 15, legend = c('E13','E17'))

t13_neuro <- table(cut(avg_meth_all[neuro_peaks,grep('E13', colnames(avg_meth_all))], breaks = meth_qs))
t17_neuro <- table(cut(avg_meth_all[neuro_peaks,grep('E17', colnames(avg_meth_all))], breaks = meth_qs))
barplot(do.call('c', lapply(seq_along(t13_neuro), function(i) c(t13_neuro[[i]], t17_neuro[[i]]))), space = rep(c(2, 0.8), length(t13_neuro)), main = 'Neuronal peaks', 
                            col = rep(c(col_key[['CthPN']], adjustcolor(col_key[['CthPN']], alpha.f = 0.2)), length(t13_nsc)))
axis(1, at=seq(2,5*(length(t13_nsc)-1),4.8), labels = axis_labels)
legend('topleft', col = c(col_key[['CthPN']], adjustcolor(col_key[['CthPN']], alpha.f = 0.2)), cex = 2,pch = 15, legend = c('E13','E17'))
dev.off()



## Fig 4H

pdf(file.path(wd, fig_4h_path), h = 1000/91, w = 750/91)
mari <- c(4.3,3,3,1)
par(mfcol = c(3,2), mar = mari, cex.lab = 1.82, cex.axis = 1.52, cex.main = 2, las = 2)
plt_lst <- c('Dnmt1','Dnmt3a','Dnmt3b','Tet1','Tet2','Tet3')
iii <- sapply(plt_lst, function(i) {
        if (i %in% c('Dnmt1','Dnmt3a','Dnmt3b')) {
            mari[[2]] <- 5
            ylabi <- 'log2 fraction'
        } else {
            mari[[2]] <- 3
            ylabi <- ''
        }
        if (i %in% c('Dnmt3b','Tet3')) {
            mari[[1]] <- 12
            xlabi <- cust_st_ord
            xaxti <- 's'
        } else {
            mari[[1]] <- 1
            xaxti <- 'n'
        }
        par(mar = mari)
        boxplot(pltmt[i,] ~ factor(mcmd$cell_type, levels = cust_st_ord), 
                ylim = c(-16.6, -12), 
                # ylim = quantile(pltmt[i,], c(0.02,0.98)), 
                # ylim = c(-16.6, -12), 
                col = col_key[cust_st_ord], yaxt = 's',xaxt = xaxti, xlab  = '', ylab = ylabi,main = i, horizontal = F)
    })
dev.off()

