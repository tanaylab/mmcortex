library(gridExtra)

# wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
wd <- ''
setwd(wd)
source(file.path(wd, 'scripts/util.r'))

load(file = './output/sequence_modeling/fig_s6_data.rda')


dir.create('./output/paper_figs/FigS6/')
device <- 'pdf'
fig_s6a_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6A.{device}'))
fig_s6b_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6B.{device}'))
fig_s6c_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6C.{device}'))
fig_s6d_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6D.{device}'))
fig_s6e_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6E.{device}'))
# fig_6def_path <- file.path(wd, glue::glue('./output/paper_figs/FigS6/FigS6DEF.{device}'))


## Fig S6A
uuu <- lapply(seq_along(tm_w_add_feat@motif_models), function(i) {
    pssmi <- tm_w_add_feat@motif_models[[i]]$pssm
    if (nrow(pssmi) > 0) {
        
        p <- suppressWarnings(prego::plot_pssm_logo(pssmi, 
                title = prego_motifs_manual[names(tm_w_add_feat@motif_models)[[i]]]))
    } else {
        p <- NULL
    }
    return(p)
})

# par(mfrow = c(4,4))
# uuu <- do.call('grid.arrange', c(lapply(seq_along(tm_w_add_feat@motif_models), 
#                 function(i) prego::plot_pssm_logo(tm_w_add_feat@motif_models[[i]]$pssm, 
#                    title = prego_motifs_manual[names(tm_w_add_feat@motif_models)[[i]]])), ncol = 4))
ggplot2::ggsave(file = fig_s6a_path, 
            arrangeGrob(grobs = uuu, ncol = 4), width = 16, height = 16)  ## save plot



## Fig S6B
pdf(fig_s6b_path, 
        h = 500/71, w = 500/71)
plot(delta_els[setdiff(names(pred_score_new), rownames(x_all))], 
        pred_score_new[setdiff(names(pred_score_new), rownames(x_all))],
        pch = 1, cex = 0.42, col = adjustcolor('black', alpha.f = 0.65), 
        main = 'IPC - NSC accessibility', xlab = 'Observed', ylab = 'Predicted')
points(delta_els[intersect(names(pred_score_new), rownames(x_all))], 
        pred_score_new[intersect(names(pred_score_new), rownames(x_all))], 
        pch = 1, cex = 0.42, col = adjustcolor('red', alpha.f = 0.25))
abline(0,1,col = 'green', lty = 2, lwd = 3)
legend('bottomright', legend = c('ENCODE SCREEN pELS/dELS', 'Cortex CREs'), 
        pch = 1, col = c('black', 'red'), cex = 1, pt.cex = 0.42, bg = 'white')
dev.off()



## Fig S6C

delta_atac_inc_peaks <- names(delta_ipc_nsc[delta_ipc_nsc > 0])
delta_atac_dec_peaks <- names(delta_ipc_nsc[delta_ipc_nsc <= 0])
# feats_here <- c('NSC_ATAC','E_box_1','T_box_1','methylation', 'E_box_4', 'MYB', 'SOX')
feats_here <- c('NSC_ATAC','E_box_1','T_box_1','methylation', 'E_box_4', 'Fos_Jun_1', 'SOX')

pdf(fig_s6b_path, 
            h = 1000/71, w = 3500/71)
par(mfcol = c(2,7), cex.lab = 3, cex.axis = 2)
for (i in 1:7) {
    if (i == 1) {
        ylabi = 'SHAP value'
        ylimi <- c(-1.5, 0.65)
        mari = c(10,6,1,2)
    } else {
        ylabi = ''
        ylimi <- c(-0.5, 0.65)
        mari = c(10,3,1,2)
    }
    par(mar = mari)
    pks_here <- intersect(rownames(preds_all), delta_atac_inc_peaks)
    plot(x_all[pks_here,feats_here[[i]]], preds_all[pks_here,feats_here[[i]]], 
         col = y_all_clvls[pks_here],
         panel.first = grid(5,5), xlab = '', ylab = ylabi, cex = 1, pch= 16, 
         xlim = quantile(x_all[,feats_here[[i]]])[c(1,5)], ylim = ylimi)
    title(xlab = paste0(feats_here[[i]], '\n', 'feature value'), line = 8)
    
    lines(c(-100,100),c(0,0), lty = 2, col = 'red', lwd = 1)
    
    pks_here <- intersect(rownames(preds_all), delta_atac_dec_peaks)
    plot(x_all[pks_here,feats_here[[i]]], preds_all[pks_here,feats_here[[i]]], 
         col = y_all_clvls[pks_here],
         panel.first = grid(5,5), xlab = '', ylab = ylabi, cex = 1, pch= 16, 
         xlim = quantile(x_all[,feats_here[[i]]])[c(1,5)], ylim = ylimi)
    title(xlab = paste0(feats_here[[i]], '\n', 'feature value'), line = 8)
    
    lines(c(-100,100),c(0,0), lty = 2, col = 'red', lwd = 1)
}
dev.off()



# par(mfcol = c(4,4), cex.lab = 3, cex.axis = 2)
# for (i in 1:16) {
#     if (i%%4 == 1) {
#         ylabi = 'SHAP value'
#         ylimi <- c(-1.5, 0.65)
#         mari = c(10,6,1,2)
#     } else {
#         ylabi = ''
#         ylimi <- c(-0.5, 0.65)
#         mari = c(10,3,1,2)
#     }
#     par(mar = mari)
#     # pks_here <- intersect(rownames(preds_all), delta_atac_inc_peaks)
#     # plot(x_all[pks_here,feats_here[[i]]], preds_all[pks_here,feats_here[[i]]], 
#     #      col = y_all_clvls[pks_here],
#     #      panel.first = grid(5,5), xlab = '', ylab = ylabi, cex = 1, pch= 16, 
#     #      xlim = quantile(x_all[,feats_here[[i]]])[c(1,5)], ylim = ylimi)
#     # title(xlab = paste0(feats_here[[i]], '\n', 'feature value'), line = 8)
    
#     # lines(c(-100,100),c(0,0), lty = 2, col = 'red', lwd = 1)
    
#     pks_here <- intersect(rownames(preds_all), delta_atac_dec_peaks)
#     plot(x_all[pks_here,colnames(preds_all)[[i]]], preds_all[pks_here,colnames(preds_all)[[i]]], 
#          col = y_all_clvls[pks_here],
#          panel.first = grid(5,5), xlab = '', ylab = ylabi, cex = 1, pch= 16, 
#          xlim = quantile(x_all[,colnames(preds_all)[[i]]])[c(1,5)], ylim = ylimi)
#     title(xlab = paste0(colnames(preds_all)[[i]], '\n', 'feature value'), line = 8)
    
#     lines(c(-100,100),c(0,0), lty = 2, col = 'red', lwd = 1)
# }



# png(file.path(wd, 'output/sequence_modeling/figs/delta_atac_vs_delta_proximal_atac.png'), h = 500, w = 600)

## Fig S6C

# pb <- multintersect(names(atac_diff_vec), names(delta_ipc_nsc), rownames(x_all))
# pdf(fig_s6c_path, h = 500/71, w = 600/71)
# par(las = 2, mar = c(13, 5, 3,1), cex.lab = 1.5)
# boxplot_vec(xvec = atac_diff_vec[pb], yvec = delta_ipc_nsc[pb], 
#         nm = 'Delta ATAC vs delta proximal ATAC', num_bins = 11, ylab = 'Delta ATAC')
# title(xlab = 'Delta proximal IPC-NSC ATAC UMIs (50kbp)', line = 11)
# text(1.5, 2, labels = paste0('cor = ', 
#                             signif(cor(atac_diff_vec[pb], 
#                                     delta_ipc_nsc[pb], 
#                                     method = 'spearman', 
#                                     use = 'pairwise.complete.obs'), 2)), 
#                                     cex = 1.5)
# dev.off()


## Fig S6D

pdf(fig_s6d_path, h = 400/71, w = 650/71)
par(las = 2, mar = c(8, 5, 3, 1))
boxplot_vec(xvec = olig2_vec, yvec = rm_ama[names(olig2_vec)], 
            nm = 'methylation vs E_box_1 energy', num_bins = 10, 
            ylab = 'NSC methylation', xlab = '')
title(xlab = 'E_box_1 energy', line = 6)
dev.off()

## Fig S6E


# pdf(file.path(wd, 'output/sequence_modeling/figs/xgboost_predicted_vs_observed_ENCODE_dELS_pELS_mmcortex.pdf'), 

## Fig S6E

# pdf(file.path(wd, 'output/sequence_modeling/figs/model_r2_vs_prox_d.pdf'), h = 250/71, w = 250/71)
pdf(fig_s6e_path, h = 250/71, w = 250/71)
par(las = 2, mar = c(6,6,2,1))
barplot(unlist(r2_add), xlab = '', ylab = '')
title(xlab = 'Neighbor counting radius [bp]', line = 4)
title(ylab = 'Improvement on model R^2', line = 4)
axis(1, at = 1.2*(1:length(unlist(r2_add)))-0.5, labels = ds)
dev.off()
