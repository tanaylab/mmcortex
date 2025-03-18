wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
library(metacell)
db_path <- file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))

load(file = file.path(wd, 'output/methylation/avg_meth_all.rda'))
load(file = file.path(wd, 'output/methylation/fig4_meth_data.rda'))
load(file = file.path(wd, 'output/mcatac/fig4_atac_data.rda'))
mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]
color_key <- unique(mcmd[,c('cell_type', 'color')])
col_key <- tibble::deframe(color_key)

mc <- scdb_mc('pl_cort')
legc <- log2(1e-5 + mc@e_gc)

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))


dir.create('./output/paper_figs/FigS4/')
device <- 'pdf'
fig_s4a_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/FigS4A.{device}'))
fig_s4b_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/FigS4B.{device}'))
fig_s4c_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/FigS4C.{device}'))
fig_s4d_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/FigS4D.{device}'))
fig_s4e_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/FigS4E.{device}'))
fig_s4f_path <- file.path(wd, glue::glue('output/paper_figs/FigS4/Figs4F.{device}'))

## Fig S4A

BW <-  0.08
v1 <- density(rowMeans(avg_meth_all[intersect(prom_peaks,rownames(avg_meth_all)),grep('E\\d\\d', colnames(avg_meth_all))]), bw = BW)
v21 <- density(rowMeans(avg_meth_all[yv,grep('E\\d\\d', colnames(avg_meth_all))]), bw = BW)
v31 <- density(rowMeans(avg_meth_all[yv2,grep('E\\d\\d', colnames(avg_meth_all))]), bw = BW)

pdf(fig_s4a_path, h = 500/71, w = 500/71)
par(cex.axis = 1.5, mar = c(5,5,5,1), cex.lab = 2, cex.main = 1.62)
plot(v1$x, v1$y/sum(v1$y), type = 'l', lwd = 2, col = 'orange', xlab = 'Mean NSC methylation', ylab = glue::glue('Density (bw = {BW})'), main = 'NSC Methylation in CREs by TSS proximity')
grid(lwd = 2, lty = 2)
lines(v21$x, v21$y/sum(v21$y), type = 'l', lwd = 2, col = 'green')
lines(v31$x, v31$y/sum(v31$y), type = 'l', lwd = 2, col = 'blue')
legend('topright', legend = c('TSSs', 'Constitutive CREs', 'Variable CREs'), 
       col = c('orange', 'blue', 'green'),
       cex = 1.5, lwd = rep(2,3), lty = rep(1,3))
dev.off()

## Fig S4B

pdf(fig_s4b_path, h=  400/71, w = 900/71)
par(cex.main = 2, mfrow = c(1,2), cex.lab = 2, mar = c(5,5,4,1), cex.axis =1.5)
hist(avg_meth_all[dist_peaks_acc,'E17'] - avg_meth_all[dist_peaks_acc,'E13'], 50,
     xlim = c(-1,1),
     xlab = 'E17 - E13 methylation', ylab = 'Count', main = 'E17 - E13 methylation\nin TSSs')
grid(lwd = 2, lty = 2)
hist(avg_meth_all[prom_peaks_acc,'E17'] - avg_meth_all[prom_peaks_acc,'E13'], 50,
     xlim = c(-1,1),
     xlab = 'E17 - E13 methylation', ylab = 'Count', main = 'E17 - E13 methylation\nin CREs')
grid(lwd = 2, lty = 2)
dev.off()

## Fig S4C

x1 <- rowMeans(a_legc_by_day_n[prom_peaks,])
y1 <- rowMeans(avg_meth_all[prom_peaks,grep('E\\d\\d', colnames(avg_meth_all))])
x2 <- rowMeans(a_legc_by_day_n[dist_peaks,])
y2 <- rowMeans(avg_meth_all[dist_peaks,grep('E\\d\\d', colnames(avg_meth_all))])
cor1 <- cor(x1, y1, method = 'spearman', use = 'pairwise.complete.obs')
cor2 <- cor(x2, y2, method = 'spearman', use = 'pairwise.complete.obs')


pdf(fig_s4c_path, h = 400/71, w = 800/71)
par(mfrow = c(1,2), mar = c(5,5,3,1), cex.main = 2, cex.lab = 2, cex.axis = 1.5)
plot(x1, y1, pch =16, cex = 0.1, main = 'TSSs', ylab = 'Mean NSC methylation', xlab = 'Mean NSC accessibility', 
     ylim = c(0,1), xlim = c(-16.6, -13.5)
    )
text(-15, 0.8, labels = glue::glue('Spearman = {round(cor1, 2)}'), adj = c(0,1), cex = 1.5)
par(mar = c(5,2,3,1), cex.main = 2, cex.lab = 2, cex.axis = 1.5)
plot(x2, y2, pch =16, cex = 0.1, main = 'CREs', ylab = '', xlab = 'Mean NSC accessibility', 
     ylim = c(0,1), xlim = c(-16.6, -13.5)
    )
text(-15, 0.8, labels = glue::glue('Spearman = {round(cor2, 2)}'), adj = c(0,1),cex = 1.5)
dev.off()

## Fig S4D

cells_h <- as.character(mcmd$metacell[which(mcmd$cell_type == 'NSC')])
lm_and_plot <- function(x, y, xlab = NULL, ylab = NULL, main = NULL) {
    lm1 <- lm(y ~ x)
    plot(x,y, xlab = xlab, ylab = ylab, 
                main = glue::glue("cor = {round(cor(x,y,method = 'pearson'), 3)}, pval = {round(cor.test(x,y,method ='pearson')$p.value, 3)}"))
#     abline(lm1$coefficients[[1]], lm1$coefficients[[2]], col= 'red', lty =2)
}

pdf(fig_s4d_path, h = 600/71, w = 1000/71)
par(mfrow = c(2,3), cex.lab = 1.5, mar = c(5,5,2,1), cex.axis = 1.5, cex.main = 1.5)
lm_and_plot(pltmt['Dnmt1',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt1', ylab = 'NSC deactivating CREs ATAC')
lm_and_plot(pltmt['Dnmt3a',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3a', ylab = 'NSC deactivating CREs ATAC')
lm_and_plot(pltmt['Dnmt3b',cells_h], pltmt['dec_peaks',cells_h], xlab = 'Dnmt3b', ylab = 'NSC deactivating CREs ATAC')
lm_and_plot(pltmt['Tet1',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet1', ylab = 'NSC activating CREs ATAC')
lm_and_plot(pltmt['Tet2',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet2', ylab = 'NSC activating CREs ATAC')
lm_and_plot(pltmt['Tet3',cells_h], pltmt['inc_peaks',cells_h], xlab = 'Tet3', ylab = 'NSC activating CREs ATAC')
dev.off()

## Fig S4E



nsc_late <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day > 16.5])
nsc_early <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day < 14.5])


pdf(fig_s4e_path, h = 500/71, w = 1000/71)
# png('./output/mcatac/figs/astro_and_early_nsc_vs_late_nsc_atac.png', h = 500, w = 1000)
par(mfrow = c(1,2), cex.lab = 1.52, cex.main = 1.5)
plot(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), a_legc_avg_ct[,'Astrocytes'], 
          pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) ATAC', ylab = 'Astrocytes ATAC', 
          main = 'Astro vs late NSC - ATAC')
abline(a =-1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '0.5 LFC','1 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), a_legc_avg_ct[,'Astrocytes'], method = 'pearson')
text(-15.5,-13, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
plot(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]),
           rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_early]), 
           pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) ATAC',
            ylab = 'NSC (mean day < 14.5) ATAC', main = 'Early vs late NSC - ATAC')
abline(a =-1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0.5,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '0.5 LFC','1 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_late]), rowMeans(a_legc[rownames(a_legc_avg_ct),nsc_early]), method = 'pearson')
text(-15.5,-13, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
dev.off()


## Fig S4F

mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))

nsc_glia_mcs <- mcmd$metacell[mcmd$cell_type %in% c('NSC', 'Astrocytes', 'OPCs')]

g <- 'Sall3'
MAR <- c(5,5,1,1)
CEX <- 2
W <- 450
H <- 350
pdf(fig_s4f_path, h = H/71, w = W/71)
par(cex.lab = 2, cex.axis = 1.5, mar = MAR)
plot(mcmd$mean_day[nsc_glia_mcs], legc[g,nsc_glia_mcs], col = mcmd$color[nsc_glia_mcs], cex = CEX,
            pch = 16, xlab = 'Metacell mean day', ylab =paste0(g, ' expression'))
# title(ylab = paste0(g, ' expression'), line = 4)
dev.off()