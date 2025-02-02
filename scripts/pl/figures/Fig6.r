library(beeswarm)

wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(misha)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))

load(file = './output/sequence_modeling/fig_6_data.rda')

dir.create('./output/paper_figs/Fig6/')
device <- 'pdf'
fig_6a_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6A.{device}'))
fig_6b_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6B.{device}'))
fig_6c_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6C.{device}'))
fig_6def_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6DEF.{device}'))



pdf(fig_6a_path, h = 500/71, w = 600/71)
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5, mar = c(5,5,2,1))
plot(delta_ipc_nsc[names(pred_meth_all)], pred_meth_all, pch = 16, cex = 0.25,
             ylab = 'Predicted delta ATAC', xlab = 'Observed delta ATAC', 
                main = 'xgboost prediction of IPC-NSC ATAC')
abline(0,1,col = 'red', lwd = 2)
text(-1.2, 1.2, labels = paste0('R^2 = ', round(mean(ev_all$Rsquare), 2)), cex = 2)
dev.off()


pdf(fig_6b_path, h = 1050/71, w = 2500/71)
par(las = 2, mar = c(14,12,2,1), cex.axis = 3, cex.lab = 6)
# inds <- which(pvl$lvl %in% levels(pvl$lvl)[1:2])
beeswarm(value ~ lvl, 
         method = 'compactswarm',
         spacing = 0.2,
         corral = 'wrap',
         data = samp_df,
         pwcol = samp_df$color, 
         ylim = c(-0.9,0.55), ylab = '', xlab = '')
title( ylab = 'SHAP value', line = 7)
grid(lwd = 5)
v <- norm_integral_over_feature_contrib
text(match(names(v), levels(samp_df$lvl)), -0.45, -0.25, col = 'black', labels = round(v,2), cex = 3, srt = 90)

dev.off()


pdf(fig_6c_path, h = 400/71, w = 650/71)
par(las = 2, mar = c(8, 5, 3, 1))
boxplot_vec(xvec = olig2_vec, yvec = rm_ama[names(olig2_vec)], 
            nm = 'methylation vs E_box_1 energy', num_bins = 10, 
            ylab = 'NSC methylation', xlab = '')
title(xlab = 'E_box_1 energy', line = 6)
dev.off()

boxplot_and_ks <- function(factor_df, factor1, factor2, y_vec, ks.alternative = 'greater', 
                           main = NULL, legend_labels = NULL, col_axis = NULL, col_boxplot = NULL, legend_x = NULL, legend_y = NULL, ylab = NULL, ylim = NULL, xlim = NULL) {
    if (!is.null(col_axis)) {xaxti = 'n'} else {xaxti = 's'}
    if (!is.null(col_boxplot)) {colbx = col_boxplot} else {colbx = 'lightgray'}
    filt_inds <- factor_df[,names(factor1)] %in% levels(factor_df[,names(factor1)])[factor1[[1]]] & 
                    factor_df[,names(factor2)] %in% levels(factor_df[,names(factor2)])[factor2[[1]]]
    factor_df <- factor_df[filt_inds,]
    y_vec <- y_vec[rownames(factor_df)]
    boxplot(y_vec ~ factor_df[,names(factor1)]*factor_df[,names(factor2)], main = main, 
                xlab = '', xaxt = xaxti, col = colbx, ylab = ylab, ylim = ylim, xlim = xlim, outcex = 0.5)
    fct_eg <- expand.grid(levels(factor_df[,names(factor1)]),
                      levels(feat_cut_df[,names(factor2)]))
    print(fct_eg)
    eg <- expand.grid(head(levels(factor_df[,names(factor1)])[factor1[[1]]], -1), 
                    tail(levels(factor_df[,names(factor1)])[factor1[[1]]], -1),
                      levels(feat_cut_df[,names(factor2)])[factor2[[1]]])
    eg <- eg[as.character(eg[,1]) != as.character(eg[,2]),]
    # for (i in seq_along(col_axis)) {axis(1, at = i, labels = paste0(fct_eg[i,], collapse = '\\.'),col = col_axis[[i]])}
    Map(axis, side=1, at=1:nrow(fct_eg), col.axis=col_axis, labels=fct_eg[,2], lwd=0, las=2)
    axis(1,at=1:nrow(fct_eg),labels=FALSE)
    eg_inds <- cbind(match(eg[,1], levels(factor_df[,names(factor1)])),
                    match(eg[,2], levels(factor_df[,names(factor1)]))) + 
                    (match(eg[,3], levels(factor_df[,names(factor2)])) - 1)*length(factor1[[1]])
    eg$ks_p_res <- p.adjust(apply(eg, 1, function(x) {
        fct1 <- factor_df[,names(factor1)]
        fct2 <- factor_df[,names(factor2)]
        pks1 <- rownames(factor_df)[fct1 == x[[1]][[1]] & fct2 == x[[3]][[1]]]
        pks2 <- rownames(factor_df)[fct1 == x[[2]][[1]] & fct2 == x[[3]][[1]]]
        ks_p <- round(ks.test(y_vec[pks1],y_vec[pks2], alternative = ks.alternative)$p.value, 3)
        return(ks_p)
    }))
    if (is.null(legend_y)) {legend_y <- 2*quantile(max(y_vec, na.rm = T))}
    if (is.null(legend_x)) {legend_x <- 1}
    legend(x = legend_x, y = legend_y, xpd = T, legend = levels(factor_df[,names(factor1)]), 
                    cex = 1.5, col = unique(col_boxplot), pch = rep(15,3))
    
    get_num_asterisks <- function(pv_vec) {
        powvec <- 5*10**seq(-4,0,1)
        print(pv_vec)
        npi <- 4 - sapply(pv_vec, function(x) which.max(powvec > x))
        print(npi)
        npi <- ifelse(npi < 1, 0, npi)
        npi <- ifelse(npi > 3, 3, npi) 
    return(npi)}

    plot_asterisks <- function(x, y, npj, ast_cex = 2.5, spacing_factor = 0.2) {
        if (npj > 1) {
            xl <- seq(from = x - spacing_factor*npj, to = x + spacing_factor*npj, length.out = npj)
        } else if (npj == 1) {xl <- x}
        else {xl <- NULL}
        points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = ast_cex, col = 'red', lwd= 1)
    }

    npi <- get_num_asterisks(eg$ks_p_res)
    print(npi)
    for (i in 1:nrow(eg_inds)) {
        lines(eg_inds[i,1:2], rep(1.175+0.125*(i%%3),2), lwd = 3.5, col = 'red')
        plot_asterisks(x = mean(eg_inds[i,1:2]), y = 1.25 +0.125*(i%%3), npj = npi[[i]])
    }

    print(eg)
    # print(apply(eg, 2, function(x)  as.character))
}


pdf(fig_6def_path, h = 500/71, w = 850/71)
par(mfrow = c(1,3),las = 1, mar = c(12,5,8,2), cex.lab = 1.2, cex.main = 1.32, cex.axis = 1.5)
boxplot_and_ks(feat_cut_df, factor1 = list('meth_cut' = c(1:3)), factor2 = list('olig2_cut' = 2:3),
                ks.alternative = 'less',y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '', 
                col_boxplot = rep(adjustcolor(c('lightblue', 'blue2', 'blue4'), alpha.f = 0.7), 3),
               col_axis = rep(c('orange1', 'orange2', 'orange4'), each = 3), legend_x = 4, 
               legend_y = 2.4, ylim = c(-1.5,1.5), ylab = 'IPC-NSC ATAC', xlim = c(3.5,9.5))

par(las = 1, mar = c(12,2,8,2), cex.main = 1.72)
boxplot_and_ks(feat_cut_df, factor1 = list('prox_atac_cut' = c(1:3)), factor2 = list('olig2_cut' = 2:3),
                ks.alternative = 'greater', y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '',
               col_boxplot = rep(adjustcolor(c('lightgreen', 'green2', 'green4'), alpha.f = 0.7), 3),
              col_axis = rep(c('orange1', 'orange2', 'orange4'), each = 3), legend_x = 4, 
              legend_y = 2.4, ylim = c(-1.5,1.5), xlim = c(3.5,9.5))

boxplot_and_ks(feat_cut_df, factor1 = list('prox_atac_cut' = c(1:3)), factor2 = list('nsc_atac_cut' = 1:3),
                ks.alternative = 'greater', y_vec = delta_ipc_nsc[rownames(feat_cut_df)], main = '',
              col_boxplot = rep(adjustcolor(c('lightgreen', 'green2', 'green4'), alpha.f = 0.7), 3),
              col_axis = rep(c('red1', 'red3', 'red4'), each = 3), legend_y = 2.4, ylim = c(-1.5,1.5), xlim = c(0.5,9.65))

dev.off()

