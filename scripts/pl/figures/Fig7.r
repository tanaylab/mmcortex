
# wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
wd <- '.'
# setwd(wd)
library(tgstat)
source(file.path(wd, 'scripts/util.r'))


mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
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


clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-.3,.3, l=length(clrmp_rel))                             



load(file = './output/MPRA/fig_7_data.rda')

plot_asterisks <- function(x, y, npi, delta_points = 0.2, cex = 3) {
    xl <- sapply(seq_along(npi), function(i) {
        if (npi[[i]] > 1) {
            return(seq(x[[i]] - delta_points*npi[[i]], x[[i]] + delta_points*npi[[i]], l = npi[[i]]))
        } else if (npi[[i]] == 1) {return(x[[i]])}
    })
    points(unlist(xl), rep(y, length(unlist(xl))), pch = 8, cex = cex, col = 'red', lwd= 3)
}

ks_on_boxplot <- function(xvec, yvec, bins_vec, alternative = 'two.sided', y_low_line = 2.3, y_high_line = 2.5) {
    xvcnf <- cut(xvec, bins_vec)
    xvc <- droplevels(xvcnf)
    ksh <- ks.test(yvec[xvc %in% head(levels(xvc), 3)], 
                yvec[xvc %in% tail(levels(xvc), 3)], alternative = alternative)
    print(ksh$p.value)
    LBINX <- match(head(levels(xvc), 3), levels(xvcnf))
    RBINX <- match(tail(levels(xvc), 3), levels(xvcnf))
    LWD <- 2
    lines(LBINX[c(1,length(LBINX))], rep(y_low_line, 2), col = 'red', lwd = LWD)
    lines(RBINX[c(1,length(RBINX))], rep(y_low_line, 2), col = 'red', lwd = LWD)
    lines(rep(mean(LBINX), 2), c(y_low_line,y_high_line), col = 'red', lwd = LWD)
    lines(rep(mean(RBINX), 2), c(y_low_line,y_high_line), col = 'red', lwd = LWD)
    lines(c(mean(LBINX), mean(RBINX)), rep(y_high_line, 2), col = 'red', lwd = LWD)
    powvec <- 5*10**seq(-4,0,1)
    npi <- 4 - which.max(powvec > ksh$p.value)
    npi <- ifelse(npi < 1, 0, npi)
    npi <- ifelse(npi > 3, 3, npi)
    plot_asterisks(mean(c(LBINX, RBINX)), y_high_line + 0.1, npi, cex = 2, delta_points = 0.4)
}

dir.create('./output/paper_figs/Fig7/')
device <- 'pdf'
fig_7b_path <- file.path(wd, glue::glue('output/paper_figs/Fig7/Fig7B.{device}'))
fig_7c_path <- file.path(wd, glue::glue('output/paper_figs/Fig7/Fig7C.{device}'))
fig_7d_path <- file.path(wd, glue::glue('output/paper_figs/Fig7/Fig7D.{device}'))
fig_7e_path <- file.path(wd, glue::glue('output/paper_figs/Fig7/Fig7E.{device}'))

## Fig 7B
pheatmap::pheatmap(tgs_cor(tbl_mpra, s = T, p = T), cluster_rows = F, cluster_cols = F, fontsize = 20, 
    col = clrmp_rel, breaks = brks_rel, 
    filename = fig_7b_path)




## Fig 7C
pheatmap::pheatmap(cor_mpra_atac, cluster_rows = F, cluster_cols = F, fontsize = 20, 
    col = clrmp_rel, breaks = brks_rel, 
    filename = fig_7c_path)


## Fig 7D
barplot_vec <- lapply(colnames(tbl_mpra), function(cni) {
    ct <- stringr::str_extract(cni, 'NSC|IPC')
    day <- stringr::str_extract(cni, 'E\\d\\d')
    signh <- ifelse(ct == 'NSC', -1, 1)
    ph <- signh * pred01[rownames(tbl_mpra)]
    matx <- cbind(tbl_atac[,grep(ct, grep(day, colnames(tbl_atac), v = T, ign = T), v = T)], 
                a_legc_avg_ct[rownames(tbl_atac),ct], ph)
    tgs_cor(matx, as.matrix(tbl_mpra[,cni]), spearman = T, p = T)
})



## Fig 7D
pdf(fig_7d_path, width = 650/71, height = 450/71)
par(las = 2, mar = c(10,5,2,1), cex.axis = 0.75)
barplot(unlist(barplot_vec), 
                width = 0.8, space = rep(c(1,0,0), 
                ncol(tbl_mpra)), 
                col = c(rep(col_key[['NSC']], 3*length(grep('NSC', colnames(tbl_mpra)))),
                        rep(col_key[['IPC']], 3*length(grep('IPC', colnames(tbl_mpra))))),
                 ylab = 'Correlation (Spearman)')
ax_labels_df <- do.call('c', lapply(colnames(tbl_mpra), function(x) {
    xh <- stringr::str_extract(x, '(NSC|IPC)_E\\d\\d')
    return(c(paste0(xh, ' ATAC'), paste0(stringr::str_extract(xh, '(NSC|IPC)'), ' mean ATAC'), 
                    paste0(ifelse(grepl('NSC', x), 'negative ', ''),  'pred score')))
    }))
axis(1, at = rep((1:3)*0.748, ncol(tbl_mpra)) + 3.25*rep(0:7, each = 3) + 0.65, labels = ax_labels_df)
text(x = 3.25*(1:8)-1.25, y = -.123, labels = gsub('_', '\n', 
            stringr::str_extract(colnames(tbl_mpra), '(NSC|IPC)_E\\d\\d')), xpd = T, cex = 1)
dev.off()



## Fig 7E
pdf(fig_7e_path, h = 1000/71, w = 1000/71)
par(mfrow = c(2,4), cex.axis = 2, cex.lab = 2, cex.main = 4, las = 2)
bins_en <- seq(-2, 2, l = 17)
par(col.main = col_key[['NSC']])
names_nc <- inds_active_mpra
marb <- 12
mart <- 3
marr <- 2
vvv <- sapply(grep('e13', colnames(pltmt_nsc_norm), ign = T, inv = T, v=T), function(cni) {
    if (cni == colnames(pltmt_nsc_norm)[[2]]) {par(mar = c(marb,5,mart,marr), cex.axis = 1.5)}
    else {par(mar = c(marb,3,mart,marr), cex.axis = 2)}
        boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_nsc_norm[names_nc,cni], 
        col = col_key[['NSC']], nm = '',
                    bins = bins_en, 
                ylab = '', ylim = c(-1,2), xlim = c(3.5,13), 
                    text_y_factor = 1., text_cex = 2, xaxt = 's', show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    text(8.5,1.9, labels = gsub('NSC_', '', gsub('_mad.score', '', cni)), cex = 2)
        ks_on_boxplot(xvec = pred01[names_nc], yvec = pltmt_nsc_norm[names_nc,cni], y_low_line = 1.2, y_high_line = 1.4,
                bins_vec = bins_en, alternative = 'less')
    title(ylab = 'NSC MPRA', line = 3)
})

par(col.main = col_key[['IPC']])
vvv <- sapply(grep('e16', colnames(pltmt_ipc_norm), ign = T, v= T, inv = T), function(cni) {
    if (cni == colnames(pltmt_ipc_norm)[[1]]) {par(mar = c(marb,5,mart,marr), cex.axis = 1.5)}
    else {par(mar =  c(marb,3,mart,marr), cex.axis = 2)}
        boxplot_vec(xvec = pred01[names_nc], yvec = pltmt_ipc_norm[names_nc,cni], 
        col = col_key[['IPC']], nm = '',
    bins = bins_en, 
        ylab = '', ylim = c(-1,3.5), xlim = c(3.5,13), text_y_factor = 1., text_cex = 2, show_text = F)
    grid(lwd = 2, col = 'lightblue', lty = 2)
    text(8.5,3.2, labels = gsub('IPC_', '', gsub('_mad.score', '', cni)), cex = 2)
        ks_on_boxplot(xvec = pred01[names_nc], yvec = pltmt_ipc_norm[names_nc,cni], bins_vec = bins_en, 
    alternative = 'greater')
    title(ylab = 'IPC MPRA', line = 3)
    title(xlab = '<- NSC  xgboost prediction  IPC ->', line = 17)
})
dev.off()