# wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
wd <- '.'
# setwd(wd)
library(vioplot)
source(file.path(wd, 'scripts/util.r'))


load(file = file.path(wd, 'output/methylation/fig4_meth_data.rda'))

load(file = file.path(wd, 'output/hic/figS5_data.rda'))

mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
mcmd <- mcmd[-c(602:603),]

color_key <- unique(mcmd[,c('cell_type', 'color')])
col_key <- tibble::deframe(color_key)

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))


dir.create('./output/paper_figs/FigS5/')
device <- 'pdf'
fig_s5a_path <- glue::glue('./output/paper_figs/FigS5/FigS5A.{device}')
fig_s5b_path <- glue::glue('./output/paper_figs/FigS5/FigS5B.{device}')
fig_s5c_path <- glue::glue('./output/paper_figs/FigS5/FigS5C.{device}')



## Fig S5A

# png('./output/hic/figs/cor_mat_ipa_by_cni.png', h = 1000, w = 1000)
pdf(fig_s5a_path, h = 1000/71, w = 1000/71)
mari <- c(1,1,1,1)
par(mfrow = c(5,5), mar = mari, cex.lab = 2)
for (cni in colnames(mat_ipa)) {
    for (cnj in colnames(mat_ipa)) {
        if (cnj == head(colnames(mat_ipa), 1)) {ylabi <- stringr::str_extract(cni, 'E\\d\\d'); mari[[2]] <- 5} else {ylabi <- ''; mari[[2]] <- 1}
        if (cni == tail(colnames(mat_ipa), 1)) {xlabi <- stringr::str_extract(cnj, 'E\\d\\d'); mari[[1]] <- 5} else {xlabi <- ''; mari[[1]] <- 1}
        par(mar = mari)
        smoothScatter(mat_ipa[,cnj], mat_ipa[,cni],
                      xlab = xlabi, ylab = ylabi)
        text(quantile(mat_ipa[,cni], 0.0001, na.rm = T), quantile(mat_ipa[,cnj], 0.9995, na.rm = T), labels = paste0('R^2 = ', round(cor_mat_ipa[cni,cnj]**2, 2)), cex = 1.2, adj = c(0,0.5))
        
    }
}
dev.off()

## Fig S5B

# png('./output/hic/figs/interval_type_num_nei_100kbp_boxplot.png', h = 400, w = 500)
pdf(fig_s5b_path, h = 400/71, w = 500/71)
par(mfrow = c(1,1), cex.lab = 1.5, cex.axis = 1.5, mar = c(6,5,2,1))
vioplot(intervalID ~ type, data = dfv, xlab = '', names = rep('', 3), ylab = 'Number of neighbors within 100kbp', col = c('purple','orange', 'darkgray'))
axis(1, at = 1:3,padj = rep(0.5, 3), labels = c('Deinsulating', 'Insulating', 'Random'))
title(xlab = 'Interval type', line = 4)
dev.off()
