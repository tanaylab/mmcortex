
wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
library(vioplot)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))


load(file = file.path(wd, 'output/methylation/avg_meth_all.rda'))

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

# ## Fig S5C

# pdf(fig_s5c_path, h = 350/71, w = 650/71)

# mari <- c(7,7,1,1)
# par(mfcol = c(1,3), mar = mari, cex.lab = 2, cex.axis = 2)
# Y_DFL <- 7
# DELTA_Y <- 2
# ncl <- setNames(unique(pval_sig_df$cre_type), names(score_lst))
# gtf <- setNames(c("NSC TSS", 'Astro TSS', 'IPC TSS'), c('nsc_gene', 'astro_gene', 'ipc_gene'))

# sss <- lapply(seq_along(score_lst), function(i) {
#     y <- score_lst[[i]]
#     nm <- names(score_lst)[[i]]
#     lapply(c('ipc_gene'), function(gt) {
#     #  lapply(c('nsc_gene', 'astro_gene'), function(gt) {
#         f <- y$gene_type == gt
#         pltmt <- y[f,grep('score', colnames(y))] - rowMeans(as.matrix(y[f,grep('score', colnames(y))]))
#         pv_df_h <- pval_sig_df[pval_sig_df$cre_type == ncl[names(score_lst)[[i]]] & pval_sig_df$tss_type == gt,]
#         if (i == 1) {
#             mari[[2]] <- 7
#         } else {
#             mari[[2]] <- 2
#         }
#         if (gt == 'ipc_gene') {
#             mari[[1]] <- 7
#         } else {
#             mari[[1]] <- 2
#         }
#         par(mar = mari)
#         boxplot(pltmt, ylim = c(-10, 23), 
#                 xaxt = 'n',
#                 ylab = '',
#                 xlab = '', main = '')
#         axis(1, at = 1:5, labels = paste0('E', 13:17))
#         if (i == 1) {
#             title(ylab = paste0(c(gtf[[gt]], 'delta SHAMAN D from mean'), collapse = '\n'), line = 3)
#         }
#         if (gt == 'ipc_gene') {
#             title(xlab =  ncl[[i]], line = 5)
#         }
#         lines(c(-5,10), rep(0,2), col = 'red')
#         if (nrow(pv_df_h) >= 1) {
#             npi <- get_num_asterisks(as.numeric(pv_df_h[,'pval']))
            
#             for (i in 1:nrow(pv_df_h)) {
#                 xh <- as.numeric(pv_df_h[i,c('x1','x2')])

#                 lines(xh, rep(Y_DFL+DELTA_Y*i - 1,2), lwd = 3.5, col = 'darkgray')

#                 plot_asterisks(x = mean(xh), y = Y_DFL +DELTA_Y*i,  spacing_factor = 0.075, npj = npi[[i]])
                
#             }
#         }
#     })
# })


# dev.off()
