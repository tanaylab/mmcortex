
wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
gsetroot('/home/aviezerl/mm10')
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



load(file = './output/MPRA/fig_7_data.rda')



plot_asterisks <- function(x, y, npi, scale_factor = 1) {
        xl <- sapply(seq_along(npi), function(i) {
            if (npi[[i]] > 1) {
                return(seq(x[[i]] - scale_factor*npi[[i]], x[[i]] + scale_factor*npi[[i]], l = npi[[i]]))
            } else if (npi[[i]] == 1) {return(x[[i]])}
        })
        points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = 3, col = 'red', lwd= 2)
    }

ks_on_boxplot <- function(xvec, yvec, bins_vec, alternative = 'two.sided') {
    xvcnf <- cut(xvec, bins_vec)
    xvc <- droplevels(xvcnf)
    ksh <- ks.test(yvec[xvc %in% head(levels(xvc), 3)], 
                yvec[xvc %in% tail(levels(xvc), 3)], alternative = alternative)
    print(ksh$p.value)
    LBINX <- match(head(levels(xvc), 3), levels(xvcnf))
    RBINX <- match(tail(levels(xvc), 3), levels(xvcnf))
    Y_LOW_LINE <- 3.3
    Y_HIGH_LINE <- 3.5
    LWD <- 2
    lines(LBINX[c(1,length(LBINX))], rep(Y_LOW_LINE, 2), col = 'red', lwd = LWD)
    lines(RBINX[c(1,length(RBINX))], rep(Y_LOW_LINE, 2), col = 'red', lwd = LWD)
    lines(rep(mean(LBINX), 2), c(Y_LOW_LINE,Y_HIGH_LINE), col = 'red', lwd = LWD)
    lines(rep(mean(RBINX), 2), c(Y_LOW_LINE,Y_HIGH_LINE), col = 'red', lwd = LWD)
    lines(c(mean(LBINX), mean(RBINX)), rep(Y_HIGH_LINE, 2), col = 'red', lwd = LWD)
    powvec <- 5*10**seq(-4,0,1)
    npi <- 4 - which.max(powvec > ksh$p.value)
    npi <- ifelse(npi < 1, 0, npi)
    npi <- ifelse(npi > 3, 3, npi)
    # print('this has run')
    plot_asterisks(mean(c(LBINX, RBINX)), Y_HIGH_LINE + 0.1, npi, scale_factor = 0.6)
}

dir.create('./output/paper_figs/Fig7/')
device <- 'pdf'
fig_7ab_path <- file.path(wd, glue::glue('./output/paper_figs/Fig7/Fig7A.{device}'))
fig_7c_to_7p_path <- file.path(wd, glue::glue('./output/paper_figs/Fig7/Fig7CHKP.{device}'))
# fig_7k_to_7p_path <- file.path(wd, glue::glue('./output/paper_figs/Fig7/Fig7KP.{device}'))
fig_7ij_path <- file.path(wd, glue::glue('./output/paper_figs/Fig7/Fig7IJ.{device}'))
fig_7qr_path <- file.path(wd, glue::glue('./output/paper_figs/Fig7/Fig7QR.{device}'))




pdf(fig_7ab_path, h = 400/71, w = 2000/71)
par(mfrow = c(1,10), cex.axis = 2, cex.lab = 2, cex.main = 4, las = 2)
bins_en <- seq(-2, 2, l = 17)
par(col.main = col_key[['NSC']])
names_nc <- inds_active_mpra
marb <- 12
mart <- 3
marr <- 1
vvv <- sapply(colnames(pltmt_nsc), function(cni) {
    if (cni == colnames(pltmt_nsc)[[1]]) {par(mar = c(marb,5,mart,marr), cex.axis = 1.5)}
    else {par(mar = c(marb,2,mart,marr), cex.axis = 2)}
    boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_nsc[names_nc,cni], 
                    nm = gsub('_mad.score', '', cni), bins = bins_en, 
                    ylab = '', ylim = c(-0.2,4.95), xlim = c(3.5,13), 
                    text_y_factor = 1., text_cex = 2, xaxt = 's', show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    ks_on_boxplot(xvec = pred01[names_nc], yvec = pltmt_nsc[names_nc,cni], 
        bins_vec = bins_en, alternative = 'less')
    # bin_and_ks(xvec = pred01[names_nc], yvec = rowMeans(pltmt_nsc[names_nc,]), bins_ks = bins_ks, ks.alternative = 'less')
    title(ylab = 'NSC MPRA', line = 3)
})
# par(mar = c(12,8,1,3))
par(col.main = col_key[['IPC']])
vvv <- sapply(colnames(pltmt_ipc), function(cni) {
    if (cni == colnames(pltmt_ipc)[[1]]) {par(mar = c(marb,5,mart,marr), cex.axis = 1.5)}
    else {par(mar =  c(marb,2,mart,marr), cex.axis = 2)}
    boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_ipc[names_nc,cni], 
    nm = gsub('_mad.score', '', cni), bins = bins_en, 
        ylab = '', ylim = c(0.04,5.95), xlim = c(3.5,13), text_y_factor = 1., text_cex = 2, show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    ks_on_boxplot(xvec = pred01[names_nc], yvec = pltmt_ipc[names_nc,cni], bins_vec = bins_en, 
    alternative = 'greater')
    # bin_and_ks(xvec = pred01[names_nc], yvec = rowMeans(pltmt_ipc[names_nc,]), bins_ks = bins_ks, ks.alternative = 'greater')
    title(ylab = 'IPC MPRA', line = 3)
    title(xlab = '<- NSC  xgboost prediction  IPC ->', line = 17)
})
dev.off()

# png('./output/MPRA/figs/test/astrocyte_specific_crs_vs_all_active_mpra_inds_ecdf.png', h = 400, w = 400)
# pdf('./output/MPRA/figs/test/astrocyte_specific_crs_vs_all_active_mpra_inds_ecdf.pdf', h = 400/71, w = 400/71)


# par(cex.lab = 2, cex.main = 1.5, cex.axis = 1.5, mar = c(5,5,3,1))
# ksx <- apply(pltmt_nsc[intersect(rownames(pltmt_nsc), 
#                 unique(nei_orig_ct$rowname[nei_orig_ct$type == 'astro_peak'])),], 1, 
#                         function(x) sum((x - min(x))*(1:ncol(pltmt_nsc)))/sum((x - min(x))))
# ksy <- apply(pltmt_nsc[setdiff(inds_active_mpra,unique(nei_orig_ct$rowname[nei_orig_ct$type == 'astro_peak'])),],
#                  1, function(x) sum((x - min(x))*(1:ncol(pltmt_nsc)))/sum((x - min(x))))
# plot(ecdf(ksx), do.points = F, col = col_key[['Astrocytes']], add = F, 
#                 xlab = 'Center of mass', ylab = 'ECDF', 
#                 main = 'Astrocyte-specific CRSs are maximally\nactivated in later NSC MPRA time points', lwd = 2)
# plot(ecdf(ksy), do.points =F, add = T, lwd = 2) 
# legend('topleft', legend = c('Astrocyte-specific CRSs',
#                     paste0('n = ', length(ksx)), 'All active CRSs',paste0('n = ', length(ksy))), 
#                         col = c(col_key[['Astrocytes']], NA, 'black', NA), lwd = 2)
# dev.off()




# pdf('./output/MPRA/figs/sum_nei_and_umis_over_50kbp_per_ipc_crs_categories_tad.pdf', h = 400/71, w = 800/71)
# par(mfrow = c(1,2))
# par(cex.lab = 2, mar = c(5,7,2,1))
# K <- 25
# LWD <- 3
# bg_col <- 'whitesmoke'
# x <- as.numeric(colnames(neu_num_neighbor_mat))
# plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat_tad), k = K, na.pad = T), col = 'white', 
#             ylim = c(-10,.4), ylab = glue::glue('Log mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(neu_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'red', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_pc_num_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# par(cex.lab = 2, mar = c(5,7,2,1))
# plot(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'white', 
#         ylim = c(0,16.3e+0), ylab = glue::glue('Log mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'red', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_pc_umi_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# dev.off()



# pdf('./output/MPRA/figs/sum_nei_and_umis_over_50kbp_per_ipc_crs_categories_tad.pdf', h = 400/71, w = 800/71)
# par(mfrow = c(1,2))
# par(cex.lab = 2, mar = c(5,9,2,1))
# K <- 25
# LWD <- 3
# bg_col <- 'whitesmoke'
# x <- as.numeric(colnames(neu_num_neighbor_mat))
# plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat_tad), k = K, na.pad = T), xlim = c(0,3e+6), col = 'white', ylim = c(-2,2), ylab = glue::glue('Log enrichment over random CRSs\nmean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(neu_num_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'red', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_pc_num_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_num_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# par(cex.lab = 2, mar = c(5,9,2,1))
# plot(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'white', ylim = c(-5,5), ylab = glue::glue('Log enrichment over random CRSs\nmean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'red', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_pc_umi_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_umi_neighbor_mat_tad), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# dev.off()


####### IPC plot ######
######################

pdf(fig_7c_to_7p_path, h = 750/71, w = 1550/71)


nri <- ol_marg$rowname %in% neu_rn
ord <- order(as.numeric(nri))
layout(mat = rbind(c(1,1,1,2,2,3,3,4,4,5,5),c(6,6,6,7,7,8,8,9,9,10,10)))
par(cex.lab = 2.5, cex.axis = 1.5, mar = c(5,5,1,1))

CEX <- 0.75
plot(ol_marg$log_nsc[ord], ol_marg$log_ipc[ord], pch = 16, col = ifelse(ol_marg$rowname[ord] %in% neu_rn, 
                'red', 'black'), xlab = 'NSC ATAC', ylab = 'IPC ATAC', 
                cex = ifelse(ol_marg$rowname[ord] %in% neu_rn, 2*CEX, CEX))
points(ol_marg$log_nsc[match(hi_ipc_lo_nsc_rn, ol_marg$rowname)], 
            ol_marg$log_ipc[match(hi_ipc_lo_nsc_rn, ol_marg$rowname)], 
            pch = 16, col = 'orange3', cex = 2*CEX)
points(ol_marg$log_nsc[match(ipc_pc_rn, ol_marg$rowname)], 
        ol_marg$log_ipc[match(ipc_pc_rn, ol_marg$rowname)], 
        pch = 16, col = col_key[['IPC']], cex = 2*CEX)
legend('bottomright', col = c('red', col_key[['IPC']], 'orange3'), 
        legend = paste0('n = ', c(length(neu_rn), length(ipc_pc_rn), length(hi_ipc_lo_nsc_rn))), 
        pch = 15, cex = 2)


par(mar = c(4,5,1,1))
locs <- seq(1,3.5,0.5)
boxplot(ol_marg$pred01[ol_marg$rowname %in% neu_rn], at = locs[[1]], col = 'red', add = F, xlim = c(0.75,2.25), ylim = c(-1.2,.51), ylab = 'NSC <-    xgb pred    -> IPC')
boxplot(ol_marg$pred01[ol_marg$rowname %in% ipc_pc_rn], at = locs[[2]], col = col_key[['IPC']], add = T)
boxplot(ol_marg$pred01[ol_marg$rowname %in% hi_ipc_lo_nsc_rn], at = locs[[3]], col = 'orange3', add = T)

par(mar = c(5,5,1,1))
locs <- seq(1,2.5,0.5)
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% neu_rn], at = locs[[1]], col = 'red', add = F, xlim = c(0.75,2.25), ylim = c(0.9,2.7), ylab = 'MPRA score')
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% ipc_pc_rn], at = locs[[2]], col = col_key[['IPC']], add = T)
boxplot(ol_marg$rm_ipc[ol_marg$rowname %in% hi_ipc_lo_nsc_rn], at = locs[[3]], col = 'orange3', add = T)

locs <- seq(1,3,l=5)
boxplot(ipc_pca$atac_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.75e+4), ylab = '', col = 'darkgray')
boxplot(ipc_pca$atac_umi[ipc_pca$is_neu == T], at = locs[[2]], add = T, col = 'red')
boxplot(ipc_pca$atac_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], at = locs[[3]], add = T, col = 'orange3')
boxplot(ipc_pca$atac_umi[match(ipc_pc_rn, ipc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['IPC']])
add_ks_stars_to_boxplot(vec1 = ipc_pca$atac_umi[ipc_pca$is_neu == T], x1 = locs[[2]], vec2 = ipc_pca$atac_umi[match(ipc_pc_rn, ipc_pca$rowname)], x2 = locs[[4]],y_line = 16.5e+3, delta_y_asterisk = 8e+2)
add_ks_stars_to_boxplot(vec1 = ipc_pca$atac_umi[ipc_pca$is_neu == T], x1 = locs[[2]], vec2 = ipc_pca$atac_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], x2 = locs[[3]],y_line = 14.8e+3, delta_y_asterisk = 8e+2)


grid(col = 'lightblue', lwd = 2, lty = 2)
title(ylab = 'Proximal ATAC UMI (50kbp)', line = 3)

locs <- seq(1,5,l=9)

boxplot(ipc_pca$rna_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), 
        ylim = c(0,1.462e+5), ylab = '', col = 'darkgray', yaxt = 'n')
boxplot(ipc_pca$rna_umi[ipc_pca$is_neu == T], at = locs[[2]], add = T, col = 'red', yaxt = 'n')
boxplot(ipc_pca$rna_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], at = locs[[3]], 
        add = T, col = 'orange3', yaxt = 'n')
# boxplot(ipc_pca$atac_umi[ipc_pca$is_ctxt == T], at = locs[[4]], add = T, col = 'green')
boxplot(ipc_pca$rna_umi[match(ipc_pc_rn, ipc_pca$rowname)], at = locs[[4]], add = T, 
        col = col_key[['IPC']], yaxt = 'n')
add_ks_stars_to_boxplot(vec1 = ipc_pca$rna_umi[ipc_pca$is_neu == T], x1 = locs[[2]], 
        vec2 = ipc_pca$rna_umi[match(ipc_pc_rn, ipc_pca$rowname)], x2 = locs[[4]],
        y_line = 1.35e+5, delta_y_asterisk = 7e+3)
add_ks_stars_to_boxplot(vec1 = ipc_pca$rna_umi[ipc_pca$is_neu == T], x1 = locs[[2]], 
        vec2 = ipc_pca$rna_umi[ipc_pca$rowname %in% nc_hi_ipc_lo_nsc_rn], x2 = locs[[3]],
        y_line = 1.18e+5, delta_y_asterisk = 7e+3)
axis(2, at = seq(0,1.5e+5, 5e+4), labels = signif(seq(0,1.5e+5, 5e+4), 2))
title(ylab = 'Proximal RNA UMI (500kbp)', line = 3)


####### NSC plot ######
######################

# pdf(fig_7k_to_7p_path, h = 750/71, w = 1550/71)

CEX <- 0.75
nri <- ol_marg$rowname %in% anti_neu_rn
ord <- order(as.numeric(nri))
plot(ol_marg$log_nsc[ord], ol_marg$log_ipc[ord], pch = 16, 
        col = ifelse(ol_marg$rowname[ord] %in% anti_neu_rn, 'magenta', 'black'), 
        xlab = 'NSC ATAC', ylab = 'IPC ATAC', cex = ifelse(ol_marg$rowname[ord] %in% anti_neu_rn, 2*CEX, CEX))
points(ol_marg$log_nsc[match(hi_nsc_lo_ipc_rn, ol_marg$rowname)], 
        ol_marg$log_ipc[match(hi_nsc_lo_ipc_rn, ol_marg$rowname)], 
        pch = 16, col = 'green', cex = 2*CEX)
points(ol_marg$log_nsc[match(nsc_pc_rn, ol_marg$rowname)],
        ol_marg$log_ipc[match(nsc_pc_rn, ol_marg$rowname)], 
        pch = 16, col = col_key[['NSC']], cex = 2*CEX)
legend('bottomright', col = c('magenta', col_key[['NSC']], 'green'), 
            legend = paste0('n = ', c(length(anti_neu_rn), length(nsc_pc_rn), length(hi_nsc_lo_ipc_rn))),
             pch = 15, cex = 2)


par(mar = c(4,5,1,1))
locs <- seq(1,3,0.5)
boxplot(ol_marg$pred01[ol_marg$rowname %in% anti_neu_rn], at = locs[[1]], col = 'magenta', add = F, 
        xlim = c(0.75,2.25), ylim = c(-1.2,.51), ylab = 'NSC <-    xgb pred    -> IPC')
boxplot(ol_marg$pred01[ol_marg$rowname %in% nsc_pc_rn], at = locs[[2]], col = col_key[['NSC']], add = T)
boxplot(ol_marg$pred01[ol_marg$rowname %in% hi_nsc_lo_ipc_rn], at = locs[[3]], col = 'green', add = T)

par(mar = c(5,5,1,1))
locs <- seq(1,3,0.5)
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% anti_neu_rn], at = locs[[1]], col = 'magenta', add = F, xlim = c(0.75,2.35), ylim = c(0.9,2.7), ylab = 'MPRA score')
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% nsc_pc_rn], at = locs[[2]], col = col_key[['NSC']], add = T)
boxplot(ol_marg$rm_nsc[ol_marg$rowname %in% hi_nsc_lo_ipc_rn], at = locs[[3]], col = 'green', add = T)
text(c(1.5,2.75), rep(0.43,2), labels = paste0('in ', c('IPC', 'NSC')), xpd = T, col = 'darkgray', cex = 2)

locs <- seq(1,3,l=5)
boxplot(nsc_pca$atac_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, 
        xlim = c(0.5,2.75), ylim = c(0,1.62e+4), ylab = '', col = 'darkgray')
boxplot(nsc_pca$atac_umi[nsc_pca$is_anti_neu == T], at = locs[[2]], add = T, col = 'magenta')
boxplot(nsc_pca$atac_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], at = locs[[3]], add = T, col = 'green')
boxplot(nsc_pca$atac_umi[match(nsc_pc_rn, nsc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['NSC']])
grid(col = 'lightblue', lwd = 2, lty = 2)
add_ks_stars_to_boxplot(vec1 = nsc_pca$atac_umi[nsc_pca$is_anti_neu == T], x1 = locs[[2]], 
        vec2 = nsc_pca$atac_umi[match(nsc_pc_rn, ipc_pca$rowname)], 
        x2 = locs[[4]],y_line = 14.5e+3, delta_y_asterisk = 8e+2)
add_ks_stars_to_boxplot(vec1 = nsc_pca$atac_umi[nsc_pca$is_anti_neu == T], 
        x1 = locs[[2]], vec2 = nsc_pca$atac_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], 
        x2 = locs[[3]],y_line = 12.8e+3, delta_y_asterisk = 8e+2)


title(ylab = 'Proximal ATAC UMI (50kbp)', line = 3)

locs <- seq(1,3,l=5)

boxplot(nsc_pca$rna_umi[ol_marg$rowname %in% nc_neu_rn], at = locs[[1]], add = F, xlim = c(0.5,2.75), ylim = c(0,1.262e+5), ylab = '', col = 'darkgray', yaxt = 'n')
boxplot(nsc_pca$rna_umi[nsc_pca$is_anti_neu == T], at = locs[[2]], add = T, col = 'magenta', yaxt = 'n')
boxplot(nsc_pca$rna_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], at = locs[[3]], add = T, col = 'green', yaxt = 'n')
boxplot(nsc_pca$rna_umi[match(nsc_pc_rn, nsc_pca$rowname)], at = locs[[4]], add = T, col = col_key[['NSC']], yaxt = 'n')
add_ks_stars_to_boxplot(vec1 = nsc_pca$rna_umi[nsc_pca$is_anti_neu == T], x1 = locs[[2]], 
            vec2 = nsc_pca$rna_umi[match(nsc_pc_rn, ipc_pca$rowname)], x2 = locs[[4]],y_line = 1.15e+5, delta_y_asterisk = 2e+3)
add_ks_stars_to_boxplot(vec1 = nsc_pca$rna_umi[nsc_pca$is_anti_neu == T], x1 = locs[[2]], 
            vec2 = nsc_pca$rna_umi[nsc_pca$rowname %in% nc_hi_nsc_lo_ipc_rn], 
            x2 = locs[[3]],y_line = 1.018e+5, delta_y_asterisk = 2e+3)
axis(2, at = seq(0,1.5e+5, 5e+4), labels = signif(seq(0,1.5e+5, 5e+4), 2))
title(ylab = 'Proximal RNA UMI (500kbp)', line = 3)
dev.off()



pdf(fig_7ij_path, h = 400/71, w = 800/71)


par(mfrow = c(1,2))
par(cex.lab = 2, mar = c(5,7,2,1))
K <- 25
LWD <- 3
bg_col <- 'whitesmoke'
x <- as.numeric(colnames(neu_num_neighbor_mat))
plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(-1,1), ylab = glue::glue('Mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = bg_col)
lines(x, log2(zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat), k = K, na.pad = T)), col = 'red', lwd = LWD)
lines(x, log2(zoo::rollmean(colMeans(ipc_pc_num_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_num_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_num_neighbor_mat), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
par(cex.lab = 2, mar = c(5,7,2,1))
plot(x, zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T), col = 'white', ylim = c(-1,1), ylab = glue::glue('Mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = bg_col)
lines(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat), k = K, na.pad = T)), col = 'red', lwd = LWD)
lines(x, log2(zoo::rollmean(colMeans(ipc_pc_umi_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat), k = K, na.pad = T)), col = col_key[['IPC']], lwd = LWD)
lines(x, log2(zoo::rollmean(colMeans(ipc_ctxt_umi_neighbor_mat), k = K, na.pad = T)/zoo::rollmean(colMeans(rand_umi_neighbor_mat), k = K, na.pad = T)), col = 'orange3', lwd = LWD)
# axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
dev.off()



pdf(fig_7qr_path, h = 400/71, w = 800/71)
par(mfrow = c(1,2))
par(cex.lab = 2, mar = c(5,7,2,1))
K <- 25
LWD <- 3
bg_col <- 'whitesmoke'
x <- as.numeric(colnames(neu_num_neighbor_mat))
plot(x, zoo::rollmean(colMeans(neu_num_neighbor_mat), k = K, na.pad = T), col = 'white', 
            ylim = c(0.8,1.62e+0), ylab = glue::glue('Mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = bg_col)
lines(x, zoo::rollmean(colMeans(anti_neu_num_neighbor_mat), k = K, na.pad = T), col = 'magenta', lwd = LWD)
lines(x, zoo::rollmean(colMeans(nsc_pc_num_neighbor_mat), k = K, na.pad = T), col = col_key[['NSC']], lwd = LWD)
lines(x, zoo::rollmean(colMeans(nsc_ctxt_num_neighbor_mat), k = K, na.pad = T), col = 'green', lwd = LWD)
# axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
par(cex.lab = 2, mar = c(5,7,2,1))
plot(x, zoo::rollmean(colMeans(neu_umi_neighbor_mat), k = K, na.pad = T), col = 'white', 
            ylim = c(8e+2,1.82e+3), ylab = glue::glue('Mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = bg_col)
lines(x, zoo::rollmean(colMeans(anti_neu_umi_neighbor_mat), k = K, na.pad = T), col = 'magenta', lwd = LWD)
lines(x, zoo::rollmean(colMeans(nsc_pc_umi_neighbor_mat), k = K, na.pad = T), col = col_key[['NSC']], lwd = LWD)
lines(x, zoo::rollmean(colMeans(nsc_ctxt_umi_neighbor_mat), k = K, na.pad = T), col = 'green', lwd = LWD)
# axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
dev.off()



# pdf('./output/MPRA/figs/test/sum_nei_and_umis_over_50kbp_per_nsc_crs_categories_tad.pdf', h = 400/71, w = 800/71)
# par(mfrow = c(1,2))
# par(cex.lab = 2, mar = c(5,7,2,1))
# K <- 25
# LWD <- 3
# bg_col <- 'whitesmoke'
# x <- as.numeric(colnames(neu_num_neighbor_mat))
# plot(x, log2(zoo::rollmean(colMeans(neu_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'white', 
#         ylim = c(-6,.1), 
#      ylab = glue::glue('Log mean num neighbors\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(anti_neu_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'magenta', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(nsc_pc_num_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['NSC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(nsc_ctxt_num_neighbor_mat_tad), k = K, na.pad = T)), col = 'green', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# par(cex.lab = 2, mar = c(5,7,2,1))
# plot(x, log2(zoo::rollmean(colMeans(neu_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'white', 
#     ylim = c(3.5,10.5), 
#      ylab = glue::glue('Log mean UMIs\nrollmean over {10*K}kbp'), xlab = 'bp')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = bg_col)
# lines(x, log2(zoo::rollmean(colMeans(anti_neu_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'magenta', lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(nsc_pc_umi_neighbor_mat_tad), k = K, na.pad = T)), col = col_key[['NSC']], lwd = LWD)
# lines(x, log2(zoo::rollmean(colMeans(nsc_ctxt_umi_neighbor_mat_tad), k = K, na.pad = T)), col = 'green', lwd = LWD)
# # axis(1, at = 1:100, labels = seq(-5e+4-1, 5e+4, l=100))
# dev.off()








# png('./output/MPRA/figs/umis_per_tads.png', h = 800, w = 1200, res = 100)
# pdf('./output/MPRA/figs/umis_per_tads.pdf', h = 8, w = 12)
# options(repr.plot.width = 18)
# options(repr.plot.height = 6)
# q_atac <- quantile(log2(1+tads_atac_umis[tb]), c(0.1,0.9))
# q_rna <- quantile(log2(1+tads_rna_umis[tb]), c(0.1,0.9))

# par(mfrow = c(2,3))
# brks <- seq(4,19,l=100)
# brks_rna <- seq(0,22,l=100)
# hist(log2(1+tads_atac_umis[tb]), breaks = brks, main = 'ATAC UMIs per TAD', xlab = 'Log2 ATAC UMIs')
# hist(log2(1+tads_atac_umis[setdiff(names(tads_atac_umis), tb)]), breaks = brks, col = 'lightblue', add = T)
# legend('topleft', legend = c('tads w RNA', 'tads w/o RNA'), col = c('gray', 'lightblue'), pch = 15)
# hist(log2(1+tads_rna_umis[tb]), , breaks = brks_rna, main = 'RNA UMIs per TAD', xlab = 'Log2 RNA UMIs')
# hist(log2(1+tads_rna_umis[setdiff(names(tads_rna_umis), tb)]), breaks = brks_rna, col = 'pink', add = T)
# legend('topleft', legend = c('tads w ATAC', 'tads w/o ATAC'), col = c('gray', 'pink'), pch = 15)

# plot(log2(1+tads_atac_umis[tb]), log2(1+tads_rna_umis[tb]), cex = 1, pch = 1, xlab = 'Log2 ATAC UMIs per TAD', ylab = 'Log2 RNA UMIs per TAD')
# points(log2(1+tads_atac_umis[tads_neu_rn]), log2(1+tads_rna_umis[tads_neu_rn]), cex = 1, pch = 16, col = 'red')
# points(log2(1+tads_atac_umis[tads_ipc_pc_rn]), log2(1+tads_rna_umis[tads_ipc_pc_rn]), cex = 1, pch = 16, col = 'blue')
# points(log2(1+tads_atac_umis[tads_hi_ipc_lo_nsc]), log2(1+tads_rna_umis[tads_hi_ipc_lo_nsc]), cex = 1, pch = 16, col = 'orange3')
# # dev.off()

# # png('./output/MPRA/figs/umis_per_tads_density.png', h = 400, w = 1200, res = 100)
# # options(repr.plot.width = 12)
# # options(repr.plot.height = 6)
# # par(mfrow = c(1,2), mar = c(5,4,1,1))
# par(mar = c(5,4,1,1))
# BW_ATAC <- .45
# den_neu_atac <- density(log2(1+tads_atac_umis[tads_neu_rn]), bw = BW_ATAC)
# den_ipc_pc_atac <- density(log2(1+tads_atac_umis[tads_ipc_pc_rn]), bw = BW_ATAC)
# den_hi_ipc_lo_nsc_atac <- density(log2(1+tads_atac_umis[tads_hi_ipc_lo_nsc]), bw = BW_ATAC)

# plot(den_neu_atac$x, den_neu_atac$y, col = 'red', type = 'l', lwd = 2, ylim = c(0,.41), xlab = 'log2 sum ATAC UMIs', ylab = 'Density')
# lines(den_ipc_pc_atac$x, den_ipc_pc_atac$y, lwd = 2, col = 'blue')
# lines(den_hi_ipc_lo_nsc_atac$x, den_hi_ipc_lo_nsc_atac$y, lwd = 2, col = 'orange3')

# filt_na <- function(x) {x[!is.na(x)]}

# BW_RNA <- .75
# den_neu_rna <- density(filt_na(log2(1+tads_rna_umis[tads_neu_rn])), bw = BW_RNA)
# den_ipc_pc_rna <- density(filt_na(log2(1+tads_rna_umis[tads_ipc_pc_rn])), bw = BW_RNA)
# den_hi_ipc_lo_nsc_rna <- density(filt_na(log2(1+tads_rna_umis[tads_hi_ipc_lo_nsc])), bw = BW_RNA)

# plot(den_neu_rna$x, den_neu_rna$y, col = 'red', type = 'l', lwd = 2, ylim = c(0,.22), xlab = 'log2 sum RNA UMIs', ylab = 'Density')
# lines(den_ipc_pc_rna$x, den_ipc_pc_rna$y, lwd = 2, col = 'blue')
# lines(den_hi_ipc_lo_nsc_rna$x, den_hi_ipc_lo_nsc_rna$y, lwd = 2, col = 'orange3')
# dev.off()

