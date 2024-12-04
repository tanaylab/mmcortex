# devtools::load_all('~/src/metacell')
library(metacell)
devtools::load_all('~/src/metacell.flow')
library(pheatmap)
wd = '/home/feshap/raid/proj/mmcortex'
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scdb_flow_init()
SEED = 1337
K = 16
set.seed(SEED)
scfigs_init("figs/")
doMC::registerDoMC(60)
nm = 'pl_cort'

mc = scdb_mc(nm)
mat = scdb_mat(nm)
mct = scdb_mctnetwork(nm)
mcf = scdb_mctnetflow(nm)
mgraph <- scdb_mgraph(nm)
mc2d <- scdb_mc2d(nm)


save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

library(matrixStats)

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         
cust_st_ord2 = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st2 = unlist(lapply(cust_st_ord2, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         


goi = c('Pou3f1', 'Pou3f2', 'Cux1', 'Cux2', 'Neurod1', 'Neurog2', 'Id4',
         'Eomes', 'Hes1', 'Apoe', 'Sox5', 'Tbr1', 'Foxp2', 'Foxp1', 'Nfia', 'Islr2', 
         'Zbtb20', 'Bcl11b', 'Fezf2', 'Satb2', 'Mef2c', 'Nhlh1', 'Tle4',
        'Rnd2',  'Runx1t1', 'Mapt', 'Mki67', 'Pcna',
        'Fabp7', 'Olig1', 'Ldb2', 'Gadd45g', 'Syt4')
marks_filt = goi

m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")
cc_genes <- union(m_genes, s_genes)

col_annot = mcmd[,c('metacell', 'cell_type', 'mean_day')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')

ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

legc <- log2(1e-05 + mc@e_gc)

mg_bon_marks <- as.data.frame(t(sapply(apply(readr::read_csv('./BonevCollab//marker_genes.tsv'), 1, 
                                             stringr::str_split, ' '), function(x) c(x[[1]][[1]], x[[1]][[length(x[[1]])]]))))

colnames(mg_bon_marks) <- c('cell_type', 'marks')

mbm_lst <- lapply(1:nrow(mg_bon_marks), function(n) stringr::str_split(mg_bon_marks$marks[[n]], ',')[[1]])
names(mbm_lst) <- mg_bon_marks$cell_type

nsc_score <- colSums(legc[mbm_lst[['NSC']],])
mat_neu_score <- colSums(legc[c('Mapt', 'Mef2c', 'Runx1t1'),])

library(princurve)

pcu <- princurve::principal_curve(cbind(nsc_score, mat_neu_score))

options(repr.plot.width = 6)
options(repr.plot.height = 6)

mbm_lst

mg_bon_marks

pltmt <- legc - rowMeans(legc)

cor_pltmt <- tgs_cor(pltmt)

options(repr.plot.width = 16)
options(repr.plot.height = 6)

pheatmap::pheatmap(cor_pltmt[which(mcmd$cell_type %in% c('Astrocytes', 'OPCs')),cust_mc_ord_st2], annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F)

pheatmap::pheatmap(pltmt[unique(unlist(mbm_lst)),which(mcmd$cell_type == 'Astrocytes')])

pheatmap::pheatmap(pltmt[unique(unlist(mbm_lst)),which(mcmd$cell_type == 'OPCs')])

legc_avg_ct <- t(tgs_matrix_tapply(legc, mcmd$cell_type, mean))



length(which(legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')] > sqrt(4)))
length(which(legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')] < - sqrt(4)))

sum(length(which(legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')] > sqrt(4))) + 
length(which(legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')] < - sqrt(4))))

diffg <- rownames(legc_avg_ct)[which((legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')]) > sqrt(4) | 
            (legc_avg_ct[,c('CthPN')] - 
             legc_avg_ct[,c('CPN_L2-3')]) < -sqrt(4))]

diffg <- rownames(legc_avg_ct)[which((legc_avg_ct[,c('SCPN')] - 
             legc_avg_ct[,c('CPN_L2-3')]) > sqrt(4) | 
            (legc_avg_ct[,c('SCPN')] - 
             legc_avg_ct[,c('CPN_L2-3')]) < -sqrt(4))]

length(diffg)

ct_neu <- c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')

sapply(ct_neu, function(x) sapply(ct_neu, function(y) c(length(which(legc_avg_ct[,x] - 
             legc_avg_ct[,y] > sqrt(4))),
            length(which(legc_avg_ct[,x] - 
             legc_avg_ct[,y] < - sqrt(4))))))

options(repr.plot.width = 7)
options(repr.plot.height = 7)

diffg_icpn_cfupn <- rownames(legc)[rowSds(legc[,which(mcmd$cell_type == 'iCPN/CfuPN')]) >= 0.58]

hc_diffg <- hclust(dist(legc[diffg_icpn_cfupn,which(mcmd$cell_type == 'iCPN/CfuPN')] - rowMeans(legc[diffg_icpn_cfupn,which(mcmd$cell_type == 'iCPN/CfuPN')])), method = 'ward.D2')

hc_diffg <- order(apply(legc[diffg_icpn_cfupn,cust_mc_ord_st[names(cust_mc_ord_st) == 'iCPN/CfuPN']], 1, function(x) sum(2**x*1:length(x))/sum(2**x)))

# p_legc_nsc_mat_neu_scores_avg_ct <- pheatmap(legc[diffg_icpn_cfupn[hc_diffg$order],cust_mc_ord_st] - rowMeans(legc[diffg_icpn_cfupn[hc_diffg$order],]), 
p_legc_nsc_mat_neu_scores_avg_ct <- pheatmap(legc[diffg_icpn_cfupn[hc_diffg],cust_mc_ord_st] - rowMeans(legc[diffg_icpn_cfupn[hc_diffg],]), 
         annotation_col = col_annot, 
                                             annotation_colors = ann_colors,
         cluster_cols = F,cluster_rows = F,
                           clustering_method = 'ward.D2',annotation_legend = F, 
                                             treeheight_row = 0, fontsize_col = 14, fontsize_row = 14, 
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-3,3,l=100))

p_legc_nsc_mat_neu_scores_avg_ct <- pheatmap(legc[diffg,cust_mc_ord_st] - rowMeans(legc[diffg,]), 
         annotation_col = col_annot, 
                                             annotation_colors = ann_colors,
         cluster_cols = F,cluster_rows = T,
                           clustering_method = 'ward.D2',annotation_legend = F, 
                                             treeheight_row = 0, fontsize_col = 14, fontsize_row = 14, 
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-3,3,l=100))

p_legc_nsc_mat_neu_scores_avg_ct <- pheatmap(legc[diffg,cust_mc_ord_st] - rowMeans(legc[diffg,]), 
         annotation_col = col_annot, 
                                             annotation_colors = ann_colors,
         cluster_cols = F,cluster_rows = T,
                           clustering_method = 'ward.D2',annotation_legend = F, 
                                             treeheight_row = 0, fontsize_col = 14, fontsize_row = 14, 
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-3,3,l=100))

options(repr.plot.width = 7)
options(repr.plot.height = 7)

sqi <- seq(-17,-5,l=100)
errb <- 0.5*sqrt(2**sqi)

plot(legc_avg_ct[,c('CPN_L2-3')],  legc_avg_ct[,c('CthPN')])
lines(sqi, sqi-2**errb, col = 'red')
lines(sqi, sqi+2**errb, col = 'red')
abline(-2,1)
abline(2,1)

options(repr.plot.width = 14)
options(repr.plot.height = 14)

plot(as.data.frame(legc_avg_ct[,c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')]))

pltmt <- legc_avg_ct[c(mbm_lst[['NSC']],'Mapt', 'Mef2c', 'Runx1t1'),cust_st_ord]

pltmt

ach <- as.data.frame(colnames(pltmt))
rownames(ach) <- ach[,1]
colnames(ach) <- 'cell_type'



p_legc_nsc_mat_neu_scores_avg_ct <- pheatmap(pltmt - rowMeans(pltmt), 
         annotation_col = ach, 
                                             annotation_colors = list(cell_type = col_key),
         cluster_cols = F,cluster_rows = F,
                           clustering_method = 'ward.D2',annotation_legend = F, treeheight_row = 0, fontsize_col = 14, fontsize_row = 14, 
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-3,3,l=100))

save_pheatmap_png(p_legc_nsc_mat_neu_scores_avg_ct, './output/metacell_model/figs/pcu_genes_avg_ct.png', h = 600, w = 600)

pltmt <- legc[c(mbm_lst[['NSC']],'Mapt', 'Mef2c', 'Runx1t1'),cust_mc_ord_st2]
pheatmap(pltmt - rowMeans(pltmt), annotation_col = col_annot, annotation_colors = ann_colors,cluster_cols = F,cluster_rows = F,
                           clustering_method = 'ward.D2',annotation_legend = F, treeheight_row = 0, fontsize_col = 12, fontsize_row = 12, 
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-3,3,l=100))

options(repr.plot.width = 6)
options(repr.plot.height = 8)

png('./output/metacell_model/figs/pcurve.png')
par(cex.main = 2, cex.lab = 2, mar = c(4,5,2,1))
plot(nsc_score, mat_neu_score, col = mcmd$color, cex = 2, pch = 16, main = 'Principal curve', ylab = 'Neuron score', xlab = 'NSC score')
points(pcu$s[,1], pcu$s[,2], pch = 1)
points(pcu$s[,1], pcu$s[,2], col = 'green', pch = 16)
dev.off()

mc_no_glia <- as.numeric(mcmd$metacell[!(mcmd$cell_type %in% c('Astrocytes', 'OPCs'))])

png('./output/metacell_model/figs/differentiation_order_by_ct_no_glia.png', h = 500, w = 350)
par(las = 2, mar = c(6,10,2,1), cex.lab = 2, cex.axis = 1.5)
# plot(1:length(pcu$ord[pcu$ord %in% mc_no_glia]), match(mcmd$cell_type[rev(pcu$ord[pcu$ord %in% mc_no_glia])], tail(cust_st_ord, -2)), 
#      col = mcmd$color[rev(pcu$ord[pcu$ord %in% mc_no_glia])], pch = 16, xlab = '', 
#      yaxt = 'n', ylab = '', cex = 2)
plot(1:length(pcu$ord), match(mcmd$cell_type[rev(pcu$ord)], cust_st_ord), 
     col = mcmd$color[rev(pcu$ord)], pch = 16, xlab = '', 
     yaxt = 'n', ylab = '', cex = 2)
# axis(2, at = 1:length(tail(cust_st_ord, -2)), labels = tail(cust_st_ord, -2))
axis(2, at = 1:length(cust_st_ord), labels = cust_st_ord)
# grid(nx = 6, ny = 0)
title(ylab = 'Cell type', line = 8)
title(xlab = 'Differentiation index', line = 4.5)
dev.off()

col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))

library(vioplot)

png('./output/metacell_model/figs/pcu_diff_order_vioplot.png', h = 650, w = 450)
par(las = 2, mar  = c(8,11,2,2), cex.lab = 2, cex.axis = 1.5)
vioplot(rev(1:length(pcu$ord)) ~ factor(mcmd$cell_type[pcu$ord], levels = cust_st_ord), col = col_key[cust_st_ord], horizontal = T, ylab = '',xlab = '')
title(xlab = 'Differentiation order', line =6)
dev.off()

# png('./output/metacell_model/figs/differentiation_order_by_ct_no_glia.png', h = 500, w = 350)
par(las = 2, mar = c(10,6,2,1), cex.lab = 2, cex.axis = 1.5)
# ind_rev <- length(pcu$ord) - pcu$ord

boxplot(pcu$ord[pcu$ord %in% mc_no_glia] ~ factor(mcmd$cell_type[pcu$ord[pcu$ord %in% mc_no_glia]], levels = tail(cust_st_ord, -2)), 
     col = col_key[tail(cust_st_ord, -2)], pch = 16, ylab = '', 
     xaxt = 'n', xlab = '', cex = 2)
axis(1, at = 1:length(tail(cust_st_ord, -2)), labels = tail(cust_st_ord, -2))
# grid(nx = 6, ny = 0)
# title(ylab = 'Cell type', line = 8)
# title(xlab = 'Differentiation index', line = 4.5)
# dev.off()

mcmd$pcu_ord <- pcu$ord

save(pcu, file = './output/metacell_model/pcurve_nsc_vs_mat_neuro.rda')

options(repr.plot.width = 20)
options(repr.plot.height = 20)

dir.create('./output/metacell_model/figs/metacell_explanation')

png('./output/metacell_model/figs/metacell_explanation/single_cells_only.png', h = 1600, w = 1600)
par(cex.lab = 3, mar = c(5,7,4,1), cex.main = 5)
plot(mc2d@sc_x, mc2d@sc_y, pch = 16, cex = 0.3, xlab = 'UMAP dim. 1', ylab = 'UMAP dim. 2', main = 'Single cells', xaxt = 'n', yaxt = 'n', bty = 'n')
dev.off()

png('./output/metacell_model/figs/metacell_explanation/metacell_assignment.png', h = 1600, w = 1600)
par(cex.lab = 3, mar = c(5,7,4,1), cex.main = 5)
plot(mc2d@sc_x, mc2d@sc_y, pch = 16, cex = 0.3, xlab = 'UMAP dim. 1', ylab = 'UMAP dim. 2', main = 'Metacell assignment', xaxt = 'n', yaxt = 'n', bty = 'n')
arrows(x0 = mc2d@sc_x, 
       y0 = mc2d@sc_y, 
       x1 = mc2d@mc_x[mc@mc[names(mc2d@sc_x)]], 
       y1 = mc2d@mc_y[mc@mc[names(mc2d@sc_y)]], lwd = 0.2, cex = 0.1)
points(mc2d@mc_x, mc2d@mc_y, pch = 1, cex = 3, col = 'pink', lwd = 2)
dev.off()


png('./output/metacell_model/figs/metacell_explanation/metacells.png', h = 1600, w = 1600)
par(cex.lab = 3, mar = c(5,7,4,1), cex.main = 5)
# plot(mc2d@sc_x, mc2d@sc_y, pch = 16, cex = 0.3, xlab = 'UMAP dim. 1', ylab = 'UMAP dim. 2', main = 'Metacell assignment', xaxt = 'n', yaxt = 'n', bty = 'n')
plot(mc2d@mc_x, mc2d@mc_y, pch = 1, cex = 3, col = 'black', lwd = 2, xlab = 'UMAP dim. 1', ylab = 'UMAP dim. 2', main = 'Metacells', xaxt = 'n', yaxt = 'n', bty = 'n')
points(mc2d@mc_x, mc2d@mc_y, pch = 16, cex = 3, col = 'pink', lwd = 2)
dev.off()


nsc_mcs <- which(mcmd$cell_type == 'NSC')

plot(mcmd$mean_day[nsc_mcs], rev(pcu$ord[pcu$ord %in% nsc_mcs]))

ct_day_mat <- tgs_matrix_tapply(t(mcmd[,grep('E\\d\\d', colnames(mcmd))]), mcmd$cell_type, sum)

ct_day_mat_norm <- t(t(ct_day_mat)/colSums(ct_day_mat))

round(ct_day_mat_norm, 3)

delta_ct_day_mat_norm <- t(apply(cbind(rep(0, nrow(ct_day_mat_norm)), ct_day_mat_norm), 1, diff))
colnames(delta_ct_day_mat_norm) <- apply(cbind(c('E13', head(colnames(ct_day_mat_norm), -1)), c('E13', tail(colnames(ct_day_mat_norm), -1))), 1, paste0, collapse  = '-->')

delta_ct_day_mat_norm <- delta_ct_day_mat_norm[cust_st_ord,]

par(mfcol = c(6,1), mar = c(1,4,3,1), las = 2, cex.main = 2)
tbtn <- sapply(1:ncol(delta_ct_day_mat_norm), function(i) barplot(delta_ct_day_mat_norm[,i], col = color_key$color[match(rownames(delta_ct_day_mat_norm), color_key$cell_type)], main = colnames(delta_ct_day_mat_norm)[[i]], ylim = c(-.1,.1), xaxt = 'n'))

options(repr.plot.width = 8)
options(repr.plot.height = 8)

png('./output/metacell_model/figs/cell_type_fraction_by_time_point.png', h = 1000, w = 600)
par(cex.axis = 2.5, cex.main = 3, mar = c(5,5,4,1))
barplot(ct_day_mat_norm[cust_st_ord,], col = color_key$color[match(cust_st_ord, color_key$cell_type)], yaxt = 'n', main = 'Cell type fraction by time point')
dev.off()

legc <- log2(1e-05 + mc@e_gc)

st_legc <- as.data.frame(t(tgs_matrix_tapply(legc, mcmd$cell_type, mean)))

feats = scdb_gset('pl_filt_lat')
feats = names(feats@gene_set)
feats = feats[feats %in% rownames(mc@e_gc)]

feats_all <- scdb_gset('pl')

nsc_mcs <- which(mcmd$cell_type == 'NSC')
fc_legc_nsc <- rowMaxs(legc[,nsc_mcs]) - rowMins(legc[,nsc_mcs])

g_cor_cc <- lapply(c('Pcna', 'Mki67', 'Top2a', 'Mcm5'), function(g) {cor_legc_g <- tgs_cor(t(legc[,nsc_mcs]), as.matrix(legc[g,nsc_mcs]), spearman = F);
                                                        return(rownames(cor_legc_g)[cor_legc_g >= 0.5 & fc_legc_nsc >= 3])})

length(unique(unlist(g_cor_cc)))

length(unlist(g_cor_cc))

genes_cor_cc <- c(names(feats_all@gene_set[feats_all@gene_set %in% as.numeric(names(table(feats_all@gene_set[cc_genes])))]),
                  unique(unlist(g_cor_cc)))

# all_cc_genes <- union(genes_cor_cc, unlist(g_cor_cc))

legc_st = t(tgs_matrix_tapply(legc[feats,], mcmd$cell_type, mean))

legc_feat_st = t(tgs_matrix_tapply(legc[feats,], mcmd$cell_type, mean))

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(legc))
# tfs

nsc_ipc_mcs <- which(mcmd$cell_type %in% c('NSC', 'IPC', 'IPC_cyc'))

ordh <- rev(pcu$ord[pcu$ord %in% nsc_ipc_mcs])
tfs_cor_pcu_nsc_ipc <- tgs_cor(t(legc[tfs,ordh]), as.matrix(1:length(nsc_ipc_mcs)), spearman = T)

tfs_cor_pcu_nsc_ipc

q_tfs_nsc_ipc <- apply(legc[tfs,nsc_ipc_mcs], 1, quantile, c(0.05,0.95))

tfs_nsc_ipc_diff <- colnames(q_tfs_nsc_ipc)[which(apply(q_tfs_nsc_ipc, 2, diff) >= 3)]

hist(tfs_cor_pcu_nsc_ipc)

tfs_cor_pcu_nsc_ipc <- rownames(tfs_cor_pcu_nsc_ipc)[abs(tfs_cor_pcu_nsc_ipc) >= 0.5 & rownames(tfs_cor_pcu_nsc_ipc) %in% tfs_nsc_ipc_diff]
# length(rownames(tfs_cor_pcu_nsc_ipc)[tfs_cor_pcu_nsc_ipc >= 0.5 & rownames(tfs_cor_pcu_nsc_ipc) %in% tfs_nsc_ipc_diff])



options(repr.plot.width = 10)
options(repr.plot.height = 16)

p_tfs_nsc_ipc_ord_pcu <- pheatmap::pheatmap(legc[tfs_cor_pcu_nsc_ipc,ordh], cluster_cols = F, annotation_col=  col_annot, annotation_colors = ann_colors, 
                   clustering_method = 'ward.D2',annotation_legend = F, treeheight_row = 0, fontsize_row = 8, show_colnames = F,
                  col = colorRampPalette(c('white', 'orange', 'red', 'purple3', 'black'))(100))

tgs_cor(t(legc[tfs_cor_pcu_nsc_ipc,ordh]), spearman = T)

options(repr.plot.width = 15)
options(repr.plot.height = 15)

# p_tfs_nsc_ipc_ord_pcu <- 
pheatmap::pheatmap(tgs_cor(t(legc[tfs_cor_pcu_nsc_ipc,which(mcmd$cell_type %in% c('IPC','IPC_cyc'))]), spearman = T), cluster_cols = T,
                   # annotation_col=  col_annot, annotation_colors = ann_colors, 
                   clustering_method = 'ward.D2',annotation_legend = F, treeheight_row = 0, fontsize_col = 12, fontsize_row = 12, show_colnames = T,
                  col = colorRampPalette(c('blue3','white','red3'))(100), breaks = seq(-1,1,l=100)
                                           )

save_pheatmap_png(p_tfs_nsc_ipc_ord_pcu, './output/metacell_model/figs/tfs_in_nsc_ipc_transition_ord_pcu.png', h = 1200, w = 950)

pheatmap::pheatmap(legc[tfs_cor_pcu_nsc_ipc,cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'IPC', 'IPC_cyca')]], cluster_cols = F, annotation_col=  col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D2',annotation_legend = F,
                  col = colorRampPalette(c('white', 'orange', 'red', 'purple3', 'black'))(100))

calc_cor_gene_legc_w_md_in_st <- function(legc, cell_type, mcmd) {
    mcs_st <- which(mcmd$cell_type == cell_type)
    cor_genes_md <- tgs_cor(t(legc[,mcs_st]), as.matrix(mcmd$mean_day[mcs_st]), spearman=T)
    names(cor_genes_md) <- rownames(legc)
    return(cor_genes_md)
}

calc_fc_gene_legc_in_st <- function(legc, cell_type, mcmd) {
    mcs_st <- which(mcmd$cell_type == cell_type)
    fc_genes <- rowMaxs(legc[,mcs_st]) - rowMins(legc[,mcs_st])
    names(fc_genes) <- rownames(legc)
    return(fc_genes)
}

cor_gene_legc_w_md_in_st <- do.call('cbind', lapply(unique(mcmd$cell_type), function(ct) calc_cor_gene_legc_w_md_in_st(legc, ct, mcmd)))

colnames(cor_gene_legc_w_md_in_st) <- unique(mcmd$cell_type)

fc_genes_legc_in_st <- do.call('cbind', lapply(unique(mcmd$cell_type), function(ct) calc_fc_gene_legc_in_st(legc, ct, mcmd)))

colnames(fc_genes_legc_in_st) <- unique(mcmd$cell_type)

nsc_dyn_tfs <- names(which(abs(cor_gene_legc_w_md_in_st[tfs,'NSC']) >= 0.5 & fc_genes_legc_in_st[tfs,'NSC'] >= log2(3)))

options(repr.plot.width = 10)
options(repr.plot.height = 8)

# cor_palmd <- tgs_cor(t(legc[,which(mcmd$cell_type == 'NSC')]), as.matrix(legc['Palmd',which(mcmd$cell_type == 'NSC')]))

# cor_palmd <- tgs_cor(t(legc[,nsc_mcs_ord]), as.matrix(zoo::rollmean(legc['Palmd',nsc_mcs_ord], k = 25, na.pad = T)), pairwise.complete.obs = T, spearman = T)

# ttt <- sapply(names(head(cor_palmd[order(cor_palmd, decreasing = T),], 20)), function(g) {
#     plot(mcmd$mean_day[nsc_mcs_ord], legc[g,nsc_mcs_ord], main = g)
#     # mcell_mc_plot_gg('pl_cort', 'Palmd', g, mc_filt = nsc_mcs_ord, use_egc = T)
# })

# ## TF set for understanding TF modules
# tfs_hi <- tfs[matrixStats::rowMaxs(legc[tfs,]) >= -13 & 
#               matrixStats::rowMaxs(legc[tfs,]) - matrixStats::rowMins(legc[tfs,]) >= 3 & 
#              matrixStats::rowSds(legc[tfs,]) >= 0.5]
# tfs_pos_cor <- lapply(cust_st_ord, function(ct) which(cor_gene_legc_w_md_in_st[tfs_hi,ct] >= 0.5))
# tfs_neg_cor <- lapply(cust_st_ord, function(ct) which(cor_gene_legc_w_md_in_st[tfs_hi,ct] >= 0.5))

# names(tfs_pos_cor) <- cust_st_ord
# names(tfs_neg_cor) <- cust_st_ord

# ## TF set which is more like "markers"
# tfs_hi <- tfs[which(log2(matrixStats::rowMaxs(mc@mc_fp[tfs,])) - log2(matrixStats::rowMins(mc@mc_fp[tfs,])) >= 3)]
# length(tfs_hi)

# legc_tf_st = t(tgs_matrix_tapply(legc[tfs,], mcmd$cell_type, mean))

get_genes_specific_to_mcs <- function(legc, mc_pos = NULL, mc_neg = NULL, cl_vec = NULL) {
    if (!is.null(mc_pos) && is.null(mc_neg)) {
        cl_vec <- ifelse(1:ncol(legc) %in% mc_pos, 1, 0)
    } else if (!is.null(mc_pos) && !is.null(mc_neg)) {
        if (!(length(intersect(mc_pos, mc_neg)) == 0)) {
            stop('mc_pos and mc_neg intersect')
        }
        legc <- legc[,c(mc_pos, mc_neg)]
        cl_vec <- c(rep(1, length(mc_pos)), rep(0, length(mc_neg)))
    }
    legc_avg <- t(tgs_matrix_tapply(legc, cl_vec, mean))
    if (ncol(legc_avg) == 2) {
        diffs <- matrixStats::rowDiffs(legc_avg)
        rownames(diffs) <- rownames(legc_avg)
        return(diffs[order( diffs[,1], decreasing = T),])
    } else {
        diffs <- t(plyr::laply(1:ncol(legc_avg), function(i) legc_avg[,i] - matrixStats::rowMaxs(legc_avg[,-i]), .parallel = T))
        colnames(diffs) <- colnames(legc_avg)
        return(lapply(1:ncol(diffs), function(i) {
            df <- diffs[which(diffs[,i] > 0.1),]
            return(df[order(df[,i], decreasing = T),])
        }))
    }
}

all_ct_tfs <- get_genes_specific_to_mcs(legc[tfs,], cl_vec = mcmd$cell_type)

cthpn_vs_cpnl23_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'CthPN'), mc_neg = which(mcmd$cell_type == 'CPN_L2-3'))

cthpn_vs_cpnl23_genes_filt <- cthpn_vs_cpnl23_genes[abs(cthpn_vs_cpnl23_genes) >= 2]

legcf <- legc[names(cthpn_vs_cpnl23_genes_filt), cust_mc_ord_st[names(cust_mc_ord_st) %in% c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')]]
legcf <-legcf - rowMeans(legcf)

library(princurve)

pcu_legcf <- principal_curve(x = t(legcf))

dim(legcf)

length(pcu_legcf$ord)

options(repr.plot.width = 8)
options(repr.plot.height = 14)

p_neuro_markers <- pheatmap::pheatmap(legcf[,pcu_legcf$ord], cluster_cols = F, col = colorRampPalette(c('blue3', 'white', 'red3'))(100), treeheight_row = 0, show_colnames = F, annotation_legend = F,
                   breaks = seq(-3,3,l=101), annotation_col = col_annot, annotation_colors = ann_colors, clustering_method = 'ward.D')

save_pheatmap_png(p_neuro_markers, './output/metacell_model/figs/neuro_markers_gradient_phm.png', h = 3000, w = 2000)

options(repr.plot.width = 7)
options(repr.plot.height = 7)

ct_neuro <- c('CPN_L2-3', 'CPN_L5_6', 'SCPN', 'CthPN')

png('./output/metacell_model/figs/cthpn_vs_cpnl23_score.png', h = 500, w = 350)
par(cex.lab = 2, mar = c(10,5,5,2), cex.axis = 1.5, cex.main = 2, las = 2)
boxplot(colSums(legcf[,pcu_legcf$ord]) ~ factor(mcmd$cell_type[as.numeric(colnames(legcf)[pcu_legcf$ord])], levels = ct_neuro), col = col_key[ct_neuro], 
        ylab = '', xlab = '')
title(ylab = 'Score', line = 3)
title(main = 'Corticofugal vs callosal\nbranch score')
title(xlab = 'Cell type', line = 8)
dev.off()

# png('./output/metacell_model/figs/cthpn_vs_cpnl23_score.png')
par(cex.lab = 2, mar = c(5,5,3,1))
plot(1:length(pcu_legcf$ord), colSums(legcf[,pcu_legcf$ord]), ylab = '', xlab = '',
     col = mcmd$color[as.numeric(colnames(legcf)[pcu_legcf$ord])], 
     pch = 16, cex = 2)
title(ylab = 'Neuron branch score', line = 3)
title(xlab = 'Principal curve order', line = 3)
# dev.off()

names(all_ct_tfs) <- sort(unique(mcmd$cell_type))

options(repr.plot.width = 8)
options(repr.plot.height = 12)

l23_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type %in% c('CPN_L2-3') & mcmd$mean_day >= 16.5), mc_neg = which(mcmd$cell_type %in% c('CPN_L5_6', 'SCPN', 'CthPN')))
l56_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type %in% c('CPN_L5_6') & mcmd$mean_day >= 16.5), mc_neg = which(mcmd$cell_type %in% c('CPN_L2-3', 'SCPN', 'CthPN')))
scpn_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type == 'SCPN' & mcmd$mean_day >= 16.5), mc_neg = which(mcmd$cell_type %in% c('CPN_L5_6', 'CPN_L2-3', 'CthPN')))
cthpn_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type == 'CthPN' & mcmd$mean_day >= 16.5), mc_neg = which(mcmd$cell_type %in% c('CPN_L5_6', 'SCPN', 'CPN_L2-3')))
mat_neuro_tfs <- unique(unlist(lapply(list(l23_tfs, l56_tfs, scpn_tfs, cthpn_tfs), function(x) head(names(x), 10))))
pheatmap(legc[mat_neuro_tfs,cust_mc_ord_st[cust_mc_ord_st %in% which(mcmd$cell_type %in% c('CPN_L5_6', 'CPN_L2-3', 'SCPN', 'CthPN'))]], cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F,
        col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100))

late_neuro_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type %in% c('CPN_L5_6', 'CPN_L2-3', 'SCPN', 'CthPN')))

ipc_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc')), mc_neg = which(mcmd$cell_type %in% c('NSC', 'CPN_L2-3', 'CthPN')))
ipc_tfs <- ipc_tfs[ipc_tfs > 1]

nsc_tfs <- get_genes_specific_to_mcs(legc[tfs,], mc_pos = which(mcmd$cell_type %in% c('NSC')), mc_neg = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc', 'OPCs', 'Astrocytes')))
nsc_tfs <- nsc_tfs[nsc_tfs > 1]

nsc_tfs <- tfs[st_legc[tfs,'NSC'] >= -14]

ipc_spec_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc')), mc_neg = which(mcmd$cell_type %in% c('NSC', 'CPN_L2-3', 'CthPN')))
ipc_spec_genes <- ipc_spec_genes[ipc_spec_genes > 2]

nsc_spec_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('NSC')), 
                                            mc_neg = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc', 'OPCs', 'Astrocytes'))
                                           )
nsc_spec_genes <- nsc_spec_genes[nsc_spec_genes > 2]

nsc_ipc_spec_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('NSC', 'IPC' ,'IPC_cyc')), 
                                            mc_neg = which(mcmd$cell_type %in% c('CthPN', 'CPN_L2-3', 'OPCs', 'Astrocytes'))
                                           )
nsc_ipc_spec_genes <- nsc_ipc_spec_genes[nsc_ipc_spec_genes > 2]

prc_genes_df <- as.data.frame(readr::read_csv('~/raid/data_other/PRC_genes.csv'))

prc_genes <- prc_genes_df$gene[1:(match('Auts2', prc_genes_df$gene)-1)]

prc_genes_in <- unlist(sapply(prc_genes, function(g) lapply(unlist(stringr::str_split(g, ';')), function(gi) grep(paste0('^', gi, '$'), rownames(mc@e_gc), ign = T, v=T))))

meth_genes <- grep('os', grep('dnmt|tet', rownames(mc@e_gc), v=T, ign = T), inv=T, v=T)
meth_genes

options(repr.plot.width = 8)
options(repr.plot.height = 8)

plot(rowMeans(legc[,nsc_mcs]), rowSds(legc[,nsc_mcs]))

options(repr.plot.width = 12)
options(repr.plot.height = 16)

nsc_legc_gq <- t(apply(legc[,nsc_mcs], 1, quantile, c(0,0.01,0.05,0.1,0.9, 0.95,0.99,1)))

length(which(apply(nsc_legc_gq[,c(3,6)], 1, diff)>= 2 & nsc_legc_gq[,6] >= -14))

nsc_dyn_genes <- rownames(nsc_legc_gq)[which(apply(nsc_legc_gq[,c(3,6)], 1, diff)>= 2 & nsc_legc_gq[,6] >= -14)]
plot(nsc_legc_gq[,3], nsc_legc_gq[,6])
points(nsc_legc_gq[nsc_dyn_genes,3], nsc_legc_gq[nsc_dyn_genes,6], col ='red', pch = 16)

nsc_asc_genes <- intersect(nsc_dyn_genes, rownames(fc_genes_legc_in_st)[cor_gene_legc_w_md_in_st[,'NSC'] >= 0.4])
nsc_desc_genes <- intersect(nsc_dyn_genes, rownames(fc_genes_legc_in_st)[cor_gene_legc_w_md_in_st[,'NSC'] <= -0.4])
# nsc_dyn_genes <- rownames(fc_genes_legc_in_st)[which(fc_genes_legc_in_st[,'NSC'] >= 3)]
nsc_dyn_genes_no_asc_desc <- setdiff(nsc_dyn_genes, union(nsc_asc_genes, nsc_desc_genes))
nsc_dyn_genes_no_asc_desc_no_cc <- setdiff(nsc_dyn_genes_no_asc_desc, unique(unlist(g_cor_cc)))

cor_nsc_legc_dyn_genes <- tgs_cor(t(legc[nsc_dyn_genes,nsc_mcs]), spearman = T)

hc_cor_nsc <- hclust(dist(cor_nsc_legc_dyn_genes), method = 'ward.D2')

ct_hc_cor_nsc <- cutree(hc_cor_nsc, k = 8)

ra <- as.data.frame(tibble::column_to_rownames(tibble::enframe(ct_hc_cor_nsc, name = 'gene', value = 'cluster'), 'gene'))

ac <- list(cluster = setNames(chameleon::distinct_colors(8)$name, 1:8))

p_cor_dyn_genes_nsc <- pheatmap::pheatmap(cor_nsc_legc_dyn_genes[hc_cor_nsc$order,hc_cor_nsc$order], show_rownames = F, show_colnames = F, 
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-1,1,l=100),
                   cluster_cols = F, cluster_rows = F, 
                   annotation_col = ra, annotation_row = ra, annotation_colors = ac)
save_pheatmap_png(p_cor_dyn_genes_nsc, './output/metacell_model/figs/nsc_gene_module_analysis/cor_dyn_genes_nsc.png', h = 1600, w = 1600)

hc_nsc_legc_dyn <- hclust(dist(legc[nsc_dyn_genes,nsc_mcs] - rowMeans(legc[nsc_dyn_genes,nsc_mcs])), method = 'ward.D')

ct_hc_nsc_legc_dyn <- cutree(hc_nsc_legc_dyn, k = 4)

p_dyn_genes_nsc <- pheatmap::pheatmap(legc[names(ct_hc_cor_nsc[hc_cor_nsc$order]),nsc_mcs] - rowMeans(legc[names(ct_hc_cor_nsc[hc_cor_nsc$order]),nsc_mcs]),
                                      fontsize_row = 6, show_rownames = T, show_colnames = F, treeheight_col = 0, annotation_legend = F,
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-3,3,l=100),
                   cluster_cols = T, cluster_rows = F, 
                   annotation_col = col_annot, annotation_row = ra, annotation_colors = c(ac, ann_colors))
save_pheatmap_png(p_dyn_genes_nsc, './output/metacell_model/figs/nsc_gene_module_analysis/all_nsc_dynamic_genes.png', h = 5200, w = 2600)

legc_nsc_dyn_genes_avg_cl <- tgs_matrix_tapply(t(legc[nsc_dyn_genes,nsc_mcs]), ct_hc_cor_nsc, mean)

options(repr.plot.height = 6)
options(repr.plot.width = 18)

png('./output/metacell_model/figs/nsc_gene_module_analysis/cc_module_mc_scatter.png', w = 1200, h = 400)
par(mfrow = c(1,3), cex.lab = 3, mar = c(5,5,2,1), cex.axis = 2)
plot(legc_nsc_dyn_genes_avg_cl[1,], legc_nsc_dyn_genes_avg_cl[6,], xlab = 'Cluster 1 MC', ylab = 'Cluster 6 MC', pch = 16, cex = 1.5)
plot(legc_nsc_dyn_genes_avg_cl[7,], legc_nsc_dyn_genes_avg_cl[6,], xlab = 'Cluster 7 MC', ylab = 'Cluster 6 MC', pch = 16, cex = 1.5)
plot(legc_nsc_dyn_genes_avg_cl[7,], legc_nsc_dyn_genes_avg_cl[1,], xlab = 'Cluster 7 MC', ylab = 'Cluster 1 MC', pch = 16, cex = 1.5)
dev.off()

quantile(colSums(mat@mat), (0:20)/20)

head(which(colSums(mat@mat) >= 2.5e+3))

mat_ds <- scm_downsamp(mat@mat, 3000)

head(colnames(mat_ds))

head(rownames(mat_ds))

nsc_sc <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'NSC')]), colnames(mat_ds))
cl1_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 1])
cl6_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 6])
cl7_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 7])
cl1_sc <- colMeans(mat_ds[cl1_genes,nsc_sc])
cl6_sc <- colMeans(mat_ds[cl6_genes,nsc_sc])
cl7_sc <- colMeans(mat_ds[cl7_genes,nsc_sc])

length(nsc_sc)

png('./output/metacell_model/figs/nsc_gene_module_analysis/cc_module_sc_scatter.png', w = 1200, h = 400)
par(mfrow = c(1,3), cex.lab = 3, mar = c(5,5,2,1), cex.axis = 2)
plot(cl1_sc, cl6_sc, pch = 16, cex = 0.5)
plot(cl1_sc, cl7_sc, pch = 16, cex = 0.5)
plot(cl7_sc, cl6_sc, pch = 16, cex = 0.5)
dev.off()

col_key <- tibble::deframe(color_key)

options(repr.plot.height = 6)
options(repr.plot.width = 12)

dir.create('./output/metacell_model/figs/nsc_gene_module_analysis')
dir.create('./output/metacell_model/figs/nsc_gene_module_analysis/pheatmaps')
dir.create('./output/metacell_model/figs/nsc_gene_module_analysis/boxplots')


vvv <- lapply(sort(unique(ct_hc_cor_nsc)), function(clj) {
    gnj <- names(ct_hc_cor_nsc)[ct_hc_cor_nsc == clj]
    png(glue::glue('./output/metacell_model/figs/nsc_gene_module_analysis/boxplots/cluster_{clj}.png'), h = 600, w = 800)
    par(las = 2, cex.main = 2, cex.lab = 3,cex.axis = 2,mar = c(14,4,4,0.5))
    boxplot(st_legc[gnj,cust_st_ord2], col = col_key[cust_st_ord2], main = paste('Cluster', clj))
    dev.off()
})

clrmp_rna <- colorRampPalette(c('white', 'orange', 'red','purple', 'black'))(100)
clrmp_rel <- colorRampPalette(c('blue3', 'white', 'red3'))(100)

options(repr.plot.height = 8)
options(repr.plot.width = 16)

vvv <- lapply(sort(unique(ct_hc_cor_nsc)), function(clj) {
    gnj <- names(ct_hc_cor_nsc)[ct_hc_cor_nsc == clj]
    pltmtj <- legc[gnj,cust_mc_ord_st2] - rowMeans(legc[gnj,nsc_mcs])
    p_clj <- pheatmap::pheatmap(pltmtj, color = clrmp_rel, treeheight_row = 0, fontsize_row = 8, show_colnames = F,
                       # breaks = seq(-16.6,-8,l=100), 
                       breaks = seq(-3,3,l=100), 
                       annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, silent = T)
    save_pheatmap_png(p_clj, glue::glue('./output/metacell_model/figs/nsc_gene_module_analysis/pheatmaps/cluster_{clj}.png'), h = 16*nrow(pltmtj), w = 2000)
})

hist(fc_genes_legc_in_st[,'NSC'], 100)



length(nsc_dyn_genes)

length(unlist(g_cor_cc))





ipc_asc_genes <- rownames(fc_genes_legc_in_st)[which(fc_genes_legc_in_st[,'IPC'] >= 1 & cor_gene_legc_w_md_in_st[,'IPC'] >= 0.4)]
ipc_desc_genes <- rownames(fc_genes_legc_in_st)[which(fc_genes_legc_in_st[,'IPC'] >= 1 & cor_gene_legc_w_md_in_st[,'IPC'] <= -0.4)]
ipc_dyn_genes <- rownames(fc_genes_legc_in_st)[which(fc_genes_legc_in_st[,'IPC'] >= 3)]
ipc_dyn_genes_no_asc_desc <- setdiff(ipc_dyn_genes, union(ipc_asc_genes, ipc_desc_genes))
# ipc_dyn_genes_no_asc_desc_no_cc <- setdiff(ipc_dyn_genes_no_asc_desc, unlist(g_cor_cc))`

length(ipc_asc_genes)
length(ipc_desc_genes)
length(ipc_dyn_genes)
length(ipc_dyn_genes_no_asc_desc)

col_key <- tibble::deframe(color_key)

mcs_in <-cust_mc_ord_st2[names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes', 'NSC', 'IPC', 'IPC_cyc')]

nsc_mcs_ord <- cust_mc_ord_st[names(cust_mc_ord_st) == 'NSC']

ipc_mcs_ord <- cust_mc_ord_st[names(cust_mc_ord_st) == 'IPC']

marks <- scdb_gset('pl_cort_marks_f')
marks <- names(marks@gene_set)

pp3 <- pheatmap::pheatmap(legc[nsc_asc_genes,nsc_mcs_ord] - rowMeans(legc[nsc_asc_genes,nsc_mcs_ord]), legend = F,
                          fontsize = 26, main = 'Time-upregulated genes', fontsize_row = 14,labels_row = ifelse(rownames(legc[nsc_asc_genes,]) %in% marks, rownames(legc[nsc_asc_genes,]), ''),
                          treeheight_row = F,show_colnames = F,cluster_rows = T, cluster_cols = F, 
                          annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
save_pheatmap_png(pp3, './output/metacell_model/figs/NSC_asc_genes_phm.png', h = 2800, w = 1800)

pp4 <- pheatmap::pheatmap(legc[nsc_desc_genes,nsc_mcs_ord] - rowMeans(legc[nsc_desc_genes,nsc_mcs_ord]), legend = T,
                          fontsize = 26, main = 'Time-downregulated genes', fontsize_row = 14,
                          labels_row = ifelse(rownames(legc[nsc_desc_genes,]) %in% marks, rownames(legc[nsc_desc_genes,]), ''),
                          treeheight_row = F,show_colnames = F,cluster_rows = T, cluster_cols = F, 
                          annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
save_pheatmap_png(pp4, './output/metacell_model/figs/NSC_desc_genes_phm.png', h = 2800, w = 1800)

nsc_mc_clust_cor_cc <- hclust(dist(t(legc[unique(unlist(g_cor_cc)),nsc_mcs_ord])), method = 'ward.D2')
cutree_nsc <- cutree(nsc_mc_clust_cor_cc, k = 3)
table(cutree_nsc)
legc_cor_cc <- legc[unique(unlist(g_cor_cc)),nsc_mcs_ord[nsc_mc_clust_cor_cc$order]] - rowMeans(legc[unique(unlist(g_cor_cc)),nsc_mcs_ord])
hc_legc_cor_cc <- hclust(dist(legc_cor_cc), method = 'ward.D2')
hc_legc_cor_cc_cols <- hclust(dist(t(legc_cor_cc)), method = 'ward.D2')
legc_cc <- legc[cc_genes_in,nsc_mcs_ord[nsc_mc_clust_cor_cc$order]] - rowMeans(legc[cc_genes_in,nsc_mcs_ord])
hc_legc_cc <- hclust(dist(legc_cc), method = 'ward.D2')

pp_cc4 <- pheatmap::pheatmap(rbind(legc_cor_cc[hc_legc_cor_cc$ord,hc_legc_cor_cc_cols$order],legc_cc[hc_legc_cc$ord,hc_legc_cor_cc_cols$order]), 
                             gaps_row = length(unique(unlist(g_cor_cc))),
                             fontsize = 26, main = 'C.C.-correlated genes', fontsize_row = 14,
                             labels_row = ifelse(rownames(legc[unique(unlist(g_cor_cc)),]) %in% marks, rownames(legc[unique(unlist(g_cor_cc)),]), ''),
                             treeheight_row = 0, treeheight_col = 0,show_colnames = F,cluster_rows = F, cluster_cols = F, clustering_method = 'ward.D2',
                             annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, 
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-3,3,l=100))
# save_pheatmap_png(pp_cc4, './output/metacell_model/figs/NSC_genes_cor_cc_phm.png', h = 3200, w = 1800)

mat1 <- legc[nsc_dyn_genes_no_asc_desc_no_cc,nsc_mcs_ord] - rowMeans(legc[nsc_dyn_genes_no_asc_desc_no_cc,nsc_mcs_ord])
mat2 <- legc[unlist(g_cor_cc),nsc_mcs_ord] - rowMeans(legc[unlist(g_cor_cc),nsc_mcs_ord])
mat3 <- legc[nsc_asc_genes,nsc_mcs_ord] - rowMeans(legc[nsc_asc_genes,nsc_mcs_ord])
mat4 <- legc[nsc_desc_genes,nsc_mcs_ord] - rowMeans(legc[nsc_desc_genes,nsc_mcs_ord])

options(repr.plot.width = 6)
options(repr.plot.height = 6)


png('./output/metacell_model/figs/std_vs_gene_groups_in_nsc.png')
par(cex.lab = 2, cex.axis = 2, mar = c(8,5,3,1), las = 2, cex.main = 2)
v1 <- c(rowSds(mat1),rowSds(mat2),rowSds(mat3),rowSds(mat4))
v2 <- c(rep('rest', nrow(mat1)),
       rep('cor_cc', nrow(mat2)),
       rep('asc', nrow(mat3)),
       rep('desc', nrow(mat4)))
boxplot(v1 ~ v2, main = 'Standard deviation per\ngene group in NSC', xlab = '', ylab = '')
dev.off()

ks.test(rowSds(mat1), rowSds(mat2), alternative = 'greater')
ks.test(rowSds(mat1), rowSds(mat3), alternative = 'greater')
ks.test(rowSds(mat1), rowSds(mat4), alternative = 'greater')

options(repr.plot.width = 10)
options(repr.plot.height = 16)

pltmt <- legc[nsc_dyn_genes_no_asc_desc_no_cc,nsc_mcs_ord] - rowMeans(legc[nsc_dyn_genes_no_asc_desc_no_cc,nsc_mcs_ord])
pp_cc5 <- pheatmap::pheatmap(pltmt[,hc_legc_cor_cc_cols$order],fontsize = 26, main = 'All other NSC LFC>=3 genes ', fontsize_row = 14,
                             legend = F,labels_row = ifelse(rownames(pltmt) %in% marks, rownames(pltmt), ''),
                             treeheight_row = 0, treeheight_col = 0,show_colnames = F,
                          cluster_rows = T, cluster_cols = F, clustering_method = 'ward.D2',
                             annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, 
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
save_pheatmap_png(pp_cc5, './output/metacell_model/figs/NSC_genes_all_other_hi_var_genes_phm.png', h = 2800, w = 1800)

pp3 <- pheatmap::pheatmap(legc[ipc_asc_genes,ipc_mcs_ord] - rowMeans(legc[ipc_asc_genes,ipc_mcs_ord]), legend = F,
                          fontsize = 26, main = 'Time-upregulated genes', fontsize_row = 14,labels_row = ifelse(rownames(legc[ipc_asc_genes,]) %in% marks, rownames(legc[ipc_asc_genes,]), ''),
                          treeheight_row = F,show_colnames = F,cluster_rows = T, cluster_cols = F, 
                          annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
# save_pheatmap_png(pp3, './output/metacell_model/figs/ipc_asc_genes_phm.png', h = 2800, w = 1800)

pp4 <- pheatmap::pheatmap(legc[ipc_desc_genes,ipc_mcs_ord] - rowMeans(legc[ipc_desc_genes,ipc_mcs_ord]), legend = T,
                          fontsize = 26, main = 'Time-downregulated genes', fontsize_row = 14,
                          labels_row = ifelse(rownames(legc[ipc_desc_genes,]) %in% marks, rownames(legc[ipc_desc_genes,]), ''),
                          treeheight_row = F,show_colnames = F,cluster_rows = T, cluster_cols = F, 
                          annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
# save_pheatmap_png(pp4, './output/metacell_model/figs/ipc_desc_genes_phm.png', h = 2800, w = 1800)

ipc_mc_clust_cor_cc <- hclust(dist(t(legc[unlist(g_cor_cc),ipc_mcs_ord])), method = 'ward.D2')
cutree_ipc <- cutree(ipc_mc_clust_cor_cc, k = 3)
table(cutree_ipc)
legc_cor_cc <- legc[unlist(g_cor_cc),ipc_mcs_ord[ipc_mc_clust_cor_cc$order]] - rowMeans(legc[unlist(g_cor_cc),ipc_mcs_ord])
hc_legc_cor_cc <- hclust(dist(legc_cor_cc), method = 'ward.D2')
hc_legc_cor_cc_cols <- hclust(dist(t(legc_cor_cc)), method = 'ward.D2')
legc_cc <- legc[cc_genes_in,ipc_mcs_ord[ipc_mc_clust_cor_cc$order]] - rowMeans(legc[cc_genes_in,ipc_mcs_ord])
hc_legc_cc <- hclust(dist(legc_cc), method = 'ward.D2')

pp_cc4 <- pheatmap::pheatmap(rbind(legc_cor_cc[hc_legc_cor_cc$ord,hc_legc_cor_cc_cols$order],legc_cc[hc_legc_cc$ord,hc_legc_cor_cc_cols$order]), gaps_row = length(unlist(g_cor_cc)),
                             fontsize = 26, main = 'C.C.-correlated genes', fontsize_row = 14,
                             labels_row = ifelse(rownames(legc[unlist(g_cor_cc),]) %in% marks, rownames(legc[unlist(g_cor_cc),]), ''),
                             treeheight_row = 0, treeheight_col = 0,show_colnames = F,cluster_rows = F, cluster_cols = F, clustering_method = 'ward.D2',
                             annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, 
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-3,3,l=100))
save_pheatmap_png(pp_cc4, './output/metacell_model/figs/ipc_genes_cor_cc_phm.png', h = 3200, w = 1800)

pltmt <- legc[ipc_dyn_genes_no_asc_desc,ipc_mcs_ord] - rowMeans(legc[ipc_dyn_genes_no_asc_desc,ipc_mcs_ord])
pp_cc5 <- pheatmap::pheatmap(pltmt,fontsize = 26, main = 'All other IPC LFC>=3 genes ', fontsize_row = 14,
                          labels_row = ifelse(rownames(pltmt) %in% marks, rownames(pltmt), ''),
                             treeheight_row = 0, treeheight_col = 0,show_colnames = F,
                          cluster_rows = T, cluster_cols = T, clustering_method = 'ward.D2',
                             annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, 
         col = colorRampPalette(c('blue3','white','red3'))(100),breaks = seq(-1.5,1.5,l=100))
# save_pheatmap_png(pp_cc5, './output/metacell_model/figs/ipc_genes_all_other_hi_var_genes_phm.png', h = 2800, w = 1800)

# pp_cc4 <- pheatmap::pheatmap(legc[union(union(nsc_asc_genes, nsc_desc_genes), unlist(g_cor_cc)),nsc_mcs_ord[nsc_mc_clust_cor_cc$order]] - 
#                              rowMeans(legc[union(union(nsc_asc_genes, nsc_desc_genes), unlist(g_cor_cc)),nsc_mcs_ord]), 
#                           treeheight_row = F,
#                           # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
#                           show_colnames = F,
#                           cluster_rows = T, cluster_cols = T, clustering_method = 'ward.D2',
#                              annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, 
#         # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
#          col = colorRampPalette(c('blue3','white','red3'))(100),
#          breaks = seq(-3,3,l=100)
#                             )

cc_genes_in <- cc_genes[cc_genes %in% rownames(legc)]

# pp_cc5 <- pheatmap::pheatmap(legc[unlist(g_cor_cc),mcs_in] - rowMeans(legc[unlist(g_cor_cc),mcs_in]), 
pp_cc5 <- pheatmap::pheatmap(legc[cc_genes_in,mcs_in] - rowMeans(legc[cc_genes_in,mcs_in]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         # breaks = seq(-2,2,l=100)
                            )

plot()

pp_cc4$tree_col$order

pp_cc1 <- pheatmap::pheatmap(legc[cc_genes[cc_genes %in% rownames(legc)],nsc_mcs_ord] - rowMeans(legc[all_cc_genes,nsc_mcs_ord]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-1.5,1.5,l=100))
save

hc_ipc_legc_mc <- hclust(dist(t(ipc_legc[ipc_dyn_genes,])), method = "ward.D2")
hc_ipc_legc_genes <- hclust(dist(ipc_legc[ipc_dyn_genes,]), method = "ward.D2")

options(repr.plot.width = 12)
options(repr.plot.height = 22)

table(cutree(hc_ipc_legc_mc, k = 3), cut(mcmd$mean_day[which(mcmd$cell_type == 'IPC')], breaks = c(12.9,14:18)))

ipc_legc_avg_cl <- t(tgs_matrix_tapply(ipc_legc, cutree(hc_ipc_legc_mc, k = 3), mean))

head(ipc_legc_avg_cl[order(ipc_legc_avg_cl[,1] - rowMaxs(ipc_legc_avg_cl[,2:3]), decreasing = T),], 20)
head(ipc_legc_avg_cl[order(ipc_legc_avg_cl[,2] - rowMaxs(ipc_legc_avg_cl[,c(1,3)]), decreasing = T),], 20)
head(ipc_legc_avg_cl[order(ipc_legc_avg_cl[,3] - rowMaxs(ipc_legc_avg_cl[,c(1,2)]), decreasing = T),], 20)

dim(ipc_legc_avg_cl)

pp3 <- pheatmap::pheatmap(ipc_legc[ipc_dyn_genes[hc_ipc_legc_genes$order],hc_ipc_legc_mc$order] - rowMeans(ipc_legc[ipc_dyn_genes[hc_ipc_legc_genes$order],hc_ipc_legc_mc$order]),
                             # legend = F,
                          # fontsize = 26, main = 'Time-upregulated genes', 
                          fontsize_row = 14,
                          labels_row = ifelse(rownames(legc[nsc_asc_genes,]) %in% marks, rownames(legc[nsc_asc_genes,]), ''),
                          # show_rownames = F,
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = F, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100)
                         )
# save_pheatmap_png(pp3, './output/metacell_model/figs/NSC_asc_genes_phm.png', h = 2800, w = 1800)

mci <- nsc_mcs[[1]]
# plot(legc[all_cc_genes,mci], rowMeans(legc[all_cc_genes,nsc_mcs])*cc_score[[1]])
plot(legc[unlist(g_cor_cc),mci], rowMeans(legc[unlist(g_cor_cc),nsc_mcs]))

options(repr.plot.width = 16)
options(repr.plot.height = 6)

# # png('./output/metacell_model/figs/nsc_modules_scatter.png', h = 400, w = 1200)

# par(mfrow = c(1,3), cex.main = 2, cex.lab = 3, mar = c(6,8,3,2.5), cex.axis = 2)
# plot(mcmd$mean_day[nsc_mcs], colMeans(legc[unlist(g_cor_cc), nsc_mcs]), cex = 1.5,pch =16, col = color_key$color[color_key$cell_type == 'NSC'],
#     main = 'Cell-cycle correlated genes', xlab = 'Mean day',
#          ylab = 'Mean expression in metacell')
# par(mar = c(6,0.5,3,2.5))
# plot(mcmd$mean_day[nsc_mcs], colMeans(legc[nsc_asc_genes, nsc_mcs]), pch =16,  cex = 1.5,col = color_key$color[color_key$cell_type == 'NSC'], ylab = '',
#         main = 'Time-upregulated genes', xlab = 'Mean day',)
# plot(mcmd$mean_day[nsc_mcs], colMeans(legc[nsc_desc_genes, nsc_mcs]), pch =16,  cex = 1.5,col = color_key$color[color_key$cell_type == 'NSC'], ylab = '',
#         main = 'Time-downregulated genes', xlab = 'Mean day',)
# # plot(mcmd$mean_day[nsc_mcs], colMeans(legc[names(ipc_genes), nsc_mcs]), pch =16, col = color_key$color[color_key$cell_type == 'NSC'])
# # plot(mcmd$mean_day[nsc_mcs], colMeans(legc[names(glia_genes), nsc_mcs]), pch =16, col = color_key$color[color_key$cell_type == 'NSC'])
# # dev.off()

nsc_km <- tglkmeans::TGL_kmeans(t(legc[nsc_dyn_genes_no_asc_desc_no_cc,which(mcmd$cell_type == 'NSC')]), k = 3, seed = 1337)

nsc_avg_cl <- t(tgs_matrix_tapply(legc[,which(mcmd$cell_type == 'NSC')], nsc_km$cluster, mean))

options(repr.plot.width = 16)
options(repr.plot.height = 16)

f <- rownames(nsc_avg_cl) %in% union(cc_genes, union(genes_cor_cc, g_cor_cc))
nsc_scatter_df <- as.data.frame(nsc_avg_cl)
colnames(nsc_scatter_df) <- paste0('NSC metacell\ncluster ', colnames(nsc_scatter_df))

col_vec <- setNames(rep('lightgray', nrow(nsc_scatter_df)), rownames(nsc_scatter_df))
col_vec[names(ipc_genes)[ipc_genes > 1]] <- color_key$color[match('IPC', color_key$cell_type)]
col_vec[names(glia_genes)[glia_genes > 1]] <- color_key$color[match('Astrocytes', color_key$cell_type)]
col_vec[f] <- 'red'
head(col_vec)

unique(col_vec)

my_trellis <- function(df, ...) {
    nc <- length(colnames(df))
    par(mfrow = c(nc,nc))
    for (cni in colnames(df)) {
        for (cnj in colnames(df)) {
            if (cni == cnj) {
                plot(0,0,col = 'white', yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty = 'n')
                text(0,0,labels = cni, cex = 4)
            } else {
                plot(df[,cni], df[,cnj], xlab = '', ylab = '', ...)
                legend('topleft', legend = c('uncategorized', 'Glia-specific', 'IPC-specific', 'CC-correlated'),
                      col = unique(col_vec), pch = rep(16, 4), cex = 1)
            }
        }
    }
}

png('./output/metacell_model/figs/nsc_clusters_scatter.png', h = 1000, w = 1000)
par(cex.axis = 1.5
    # , mar = c(5,5,4,1)
   )
my_trellis(nsc_scatter_df, pch = 16, cex = 0.75, col = col_vec)
dev.off()

plot(nsc_scatter_df[nsc_dyn_genes_no_asc_desc_no_cc,], pch =16,cex = 0.75
     # col = col_vec, pch = 16, cex = ifelse(col_vec == 'black', 0.5, 1.5)
    )



# plot(nsc_scatter_df, col = 'black', pch = 1, cex = 1
#      # ifelse(col_vec == 'black', 0.5, 1.5)
#     )
par(cex.lab = 2)
plot(nsc_scatter_df, col = col_vec, pch = 16, cex = 1
     # ifelse(col_vec == 'black', 0.5, 1.5)
    )

mcs_in <-cust_mc_ord_st2[names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes', 'NSC', 'IPC', 'IPC_cyc')]

options(repr.plot.width = 12)
options(repr.plot.height = 6)

pltmt <- 

pp3 <- pheatmap::pheatmap(legc_avg_gt[,cust_mc_ord_st2] - rowMeans(legc_avg_gt), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), 
                          fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-2,2,l=100))
save_pheatmap_png(pp3, './output/metacell_model/figs/NSC_gene_module_across_time.png', h = 500, w = 1200)

options(repr.plot.width = 10)
options(repr.plot.height = 6)

pp3 <- pheatmap::pheatmap(legc_avg_gt[,cust_mc_ord_st2] - rowMeans(legc_avg_gt), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-2,2,l=100))

pp3 <- pheatmap::pheatmap(legc_avg_gt[,mcs_in] - rowMeans(legc_avg_gt), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-2,2,l=100))

options(repr.plot.width = 12)
options(repr.plot.height = 12)

gf <- gene_type_vec == 'NSC_hi-var'
# pp3 <- pheatmap::pheatmap(legc[gf,mcs_in] - rowMeans(legc[gf,mcs_in]), 
pp3 <- pheatmap::pheatmap(legc[gf,mcs_in] - rowMeans(legc[gf,nsc_mcs]), show_rownames = F, main = 'High-variance genes, excluding\ntime- and cell cycle-correlated', fontsize = 20,
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))
save_pheatmap_png(pp3, './output/metacell_model/figs/NSC_hi_var_genes_no_asc_desc_cc.png', h = 1800, w = 1000)

pp3 <- pheatmap::pheatmap(legc[nsc_asc_genes,mcs_in] - rowMeans(legc[nsc_asc_genes,mcs_in]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

pp3 <- pheatmap::pheatmap(legc[nsc_desc_genes,mcs_in] - rowMeans(legc[nsc_desc_genes,mcs_in]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

options(repr.plot.width = 12)
options(repr.plot.height = 20)

nsc_mcs <- which(mcmd$cell_type == 'NSC')

ml <- legc[nsc_dyn_genes_no_asc_desc_no_cc,nsc_mcs]
gg <- head(rownames(ml)[order(rowSds(ml), decreasing = T)], 50)

pheatmap::pheatmap(legc[gg,cust_mc_ord_st2] - rowMeans(legc[gg,]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

pp3 <- pheatmap::pheatmap(legc[nsc_dyn_genes_no_asc_desc_no_cc,mcs_in] - rowMeans(legc[nsc_dyn_genes_no_asc_desc_no_cc,mcs_in]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = T, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

pp3 <- pheatmap::pheatmap(legc[nsc_dyn_genes_no_asc_desc_no_cc,mcs_in] - rowMeans(legc[nsc_dyn_genes_no_asc_desc_no_cc,mcs_in]), 
                          treeheight_row = F,
                          # gaps_row = c(4, 4+length(nsc_dyn_tfs), 4+length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, 
                          show_colnames = F,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

n_astro_g <- 50
n_ipc_g <- 50

g1 <- !(colnames(nsc_avg_cl) %in% union(head(names(ipc_genes), n_ipc_g), head(names(glia_genes), n_astro_g)))
plot(nsc_avg_cl[1,g1], nsc_avg_cl[2,g1], pch = 16)
points(nsc_avg_cl[1,head(names(ipc_genes), n_ipc_g)], 
       nsc_avg_cl[2,head(names(ipc_genes), n_ipc_g)], pch = 16, cex = 2,
       col = color_key$color[match('IPC', color_key$cell_type)])
points(nsc_avg_cl[1,head(names(glia_genes), n_astro_g)], nsc_avg_cl[2,head(names(glia_genes), n_astro_g)], 
       pch = 16,cex = 2,
       col = color_key$color[match('Astrocytes', color_key$cell_type)])

options(repr.plot.width = 12)
options(repr.plot.height = 12)

f <- colnames(nsc_avg_cl) %in% union(cc_genes, genes_cor_cc)
nsc_scatter_df <- as.data.frame(t(nsc_avg_cl))
colnames(nsc_scatter_df) <- paste0('NSC metacell cluster ', colnames(nsc_scatter_df))
plot(nsc_scatter_df, col = ifelse(f, 'red', 'black'), pch = 16, cex = ifelse(f, 1, 0.5))

mcs_in <- cust_mc_ord_st2[names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes', 'NSC', 'IPC', 'IPC_cyc')]

mcs_in <- cust_mc_ord_st2

m4 <- legc[c('Top2a', 'Mki67', 'Pcna', 'Mcm5'),]
m11 <- legc[intersect(tfs, nsc_asc_genes),]
m12 <- legc[intersect(tfs, nsc_desc_genes),]
m2 <- legc[meth_genes,]
m3 <- legc[prc_genes_in,]
hc11 <- hclust(dist(m11[,mcs_in] - rowMeans(m11[,mcs_in])), method = 'ward.D2')$order
hc12 <- hclust(dist(m12[,mcs_in] - rowMeans(m12[,mcs_in])), method = 'ward.D2')$order
hc2 <- hclust(dist(m2[,mcs_in] - rowMeans(m2[,mcs_in])), method = 'ward.D2')$order
hc3 <- hclust(dist(m3[,mcs_in] - rowMeans(m3[,mcs_in])), method = 'ward.D2')$order
hc4 <- hclust(dist(m4[,mcs_in] - rowMeans(m4[,mcs_in])), method = 'ward.D2')$order
ma <- rbind(m11[hc11,],m12[hc12,], m2[hc2,], m3[hc3,])

ls()

pp3 <- pheatmap::pheatmap(ma[,cust_mc_ord_st2] - rowMeans(ma), gaps_row = cumsum(c(nrow(m11), nrow(m12), nrow(m2), nrow(m3))), main = 'Regulatory genes differentially\nexpressed in NSCs', fontsize = 20,
                          fontsize_row = 8, show_colnames = F,
                          cluster_rows = F, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))
save_pheatmap_png(pp3, './output/metacell_model/figs/tfs_meth_prc.png', h = 2800, w = 1600)

hc_cols <- hclust(dist(t(ma[!(rownames(ma) %in% cc_genes),])), method = 'ward.D2')$order

pp4 <- pheatmap::pheatmap(ma[,hc_cols] - rowMeans(ma), gaps_row = c(length(nsc_dyn_tfs), length(nsc_dyn_tfs) + length(meth_genes)), fontsize_row = 16, show_colnames = F,
                          cluster_rows = F, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('blue3','white','red3'))(100),
         breaks = seq(-3,3,l=100))

save_pheatmap_png(pp3, './output/metacell_model/nsc_dyn_tf_meth_prc_genes.png', h = 2200, w = 1500)

pp3 <- pheatmap::pheatmap(ma[,cust_mc_ord_st2[names(cust_mc_ord_st2) %in% c('OPCs', 'Astrocytes', 'NSC', 'IPC', 'IPC_cyc')]], fontsize_row = 14,
                          cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
        # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
         col = colorRampPalette(c('white', 'orange','red3', 'purple', 'black'))(100),
         # breaks = seq(-3,3,l=100)
                         )

# putative_gliogenic <- names(which(cor_gene_legc_w_md_in_st[,'NSC'] >= 0.5 & 
#                                   (st_legc[,'Astrocytes'] >= -15 | st_legc[,'OPCs'] >= -15) & 
#                                    (st_legc[,'CthPN'] <= -15 & st_legc[,'CPN_L2-3'] <= -15)))
# length(putative_gliogenic)

# nsc_asc_no_glia <- names(which(cor_gene_legc_w_md_in_st[,'NSC'] >= 0.5 & 
#                                   (st_legc[,'Astrocytes'] <= -15 & st_legc[,'OPCs'] >= -15)
#                                    # (st_legc[,'CthPN'] <= -15 | st_legc[,'CPN_L2-3'] <= -15))
#                         ))
# length(nsc_asc_no_glia)

# putative_neurogenic <- names(which(cor_gene_legc_w_md_in_st[,'NSC'] <= -0.5 & 
#                                    (st_legc[,'IPC'] >= -15 | st_legc[,'IPC_cyc'] >= -15) & 
#                                    (st_legc[,'Astrocytes'] <= -15 & st_legc[,'OPCs'] <= -15)))
# length(putative_neurogenic)

# pltmt <- legc[c(putative_gliogenic, putative_neurogenic),]
# p_nsc_clusters <- pheatmap(pltmt[,mcs_ipc] - rowMeans(pltmt[,mcs_ipc]),cluster_rows = F, show_colnames = F, legend = F,treeheight_row = F,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
#          col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         breaks = seq(-3,3,l=100))

# mcs_nsc <- cust_mc_ord_st[cust_mc_ord_st %in% which(mcmd$cell_type %in% c('Astrocytes', 'NSC', 'IPC', 'IPC_cyc'))]
# x <- sort(cor_gene_legc_w_md_in_st[,'NSC'])
# y1 <- names(head(x, 150))
# y2 <- names(tail(x, 150))
# x1 <- setdiff(union(names(which(cor_gene_legc_w_md_in_st[,'NSC'] <= -0.7 & st_legc[,'IPC'] >= -15)), names(nsc_spec_genes)), union(c(y1, y2), cc_genes))
# m1 <- legc[x1,mcs_nsc]
# hc1 <- hclust(dist(m1 - rowMeans(m1)), method = 'ward.D2')
# m2 <- legc[y1,mcs_nsc]
# hc2 <- hclust(dist(m2 - rowMeans(m2)), method = 'ward.D2')
# m3 <- legc[y2,mcs_nsc]
# hc3 <- hclust(dist(m3 - rowMeans(m3)), method = 'ward.D2')
# mall <- rbind(m1[hc1$order,],m2[hc2$order,],m3[hc3$order,])

# nsc_mcs_ord <- cust_mc_ord_st[cust_mc_ord_st %in% which(mcmd$cell_type == 'NSC')]

# mat_nsc <- legc[,nsc_mcs_ord] - rowMeans(legc[,nsc_mcs_ord])

# km_nsc <- tglkmeans::TGL_kmeans(mat_nsc, k = 32, seed = 1337)

# mat_avg_km_n <- tgs_matrix_tapply(t(mat_nsc), km_nsc$cluster, mean)

# mat_all_avg_km_n <- tgs_matrix_tapply(t(legc), km_nsc$cluster, mean)

mmc_gm <- readr::read_csv('./data/mmcortex_gene_modules-2023-04-01.csv')
mmc_gm <- mmc_gm[mmc_gm$gene %in% rownames(legc),]

mmc_gm_mat <- legc[mmc_gm$gene,]
mmc_gm_avg_module <- tgs_matrix_tapply(t(legc[mmc_gm$gene,]), mmc_gm$module, mean)

options(repr.plot.width = 16)
options(repr.plot.height = 8)

p_mod_avg_hm <- pheatmap::pheatmap(mmc_gm_avg_module[,cust_mc_ord_st2] - rowMeans(mmc_gm_avg_module), cluster_rows = T, show_colnames = F, legend = T,treeheight_row = F, 
         # annotation_row = rah,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))
save_pheatmap_png(p_mod_avg_hm, './output/metacell_model/figs/avg_mod_heatmap.png', h = 1000, w = 1200)

# delta_legc_mrm <- setNames(rowMaxs(legc_mrm) - rowMins(legc_mrm), rownames(legc_mrm))
# delta_legc_mrm

# mod_in <- which(delta_legc_mrm >= 3.5)
# mod_out <- which(delta_legc_mrm < 3.5)

# length(mod_in)

# length(mod_out)

# options(repr.plot.width = 10)
# options(repr.plot.height = 35)

# png('./output/metacell_model/figs/differential_gene_modules_nsc_ipc_scatter.png', h = 3500, w = 750)
# par(mfrow = c(length(mod_in), 2), cex.lab = 2, cex.axis = 2)
# inds_nsc <- cust_mc_ord_st[names(cust_mc_ord_st) == 'NSC']
# inds_ipc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('IPC', 'IPC_cyc')]
# legc_mrm <- mmc_gm_avg_module - rowMeans(mmc_gm_avg_module)
# tvt <- sapply(names(mod_in), function(md) {
#     if (md == tail(names(mod_in), 1)) {
#         par(mar = c(2,6,0.5,0))
#     } else {
#         par(mar = c(.51,6,.51,0))
#     }
#     mal <- max(abs(legc_mrm[md,c(inds_ipc,inds_nsc)]))
#     if (mal > 3) {
#         ylimi <- c(-mal, mal)
#     } else {ylimi <- c(-3,3)}
#     plot(mcmd$mean_day[inds_nsc], legc_mrm[md,inds_nsc], pch = 16, cex = 2, col = mcmd$color[inds_nsc], 
#          xlim = c(13,18), ylim = ylimi, ylab = md, xaxt ='n', xlab = '')
#     grid(col = 'pink')
#     if (md == tail(names(mod_in), 1)) {
#         axis(1, at = 13:18, labels = 13:18)
#     }
#     if (md == tail(names(mod_in), 1)) {
#         par(mar = c(2,1,0.5,0))
#     } else {
#         par(mar = c(.51,1,.51,0))
#     }
#     plot(mcmd$mean_day[inds_ipc], legc_mrm[md,inds_ipc], pch = 16, cex = 2, col = mcmd$color[inds_ipc], 
#          xlim = c(13,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     grid(col = 'pink')
#     if (md == tail(names(mod_in), 1)) {
#         axis(1, at = 13:18, labels = 13:18)
#     }
# })
# dev.off()

# unique(mcmd$cell_type[inds_mc])

# options(repr.plot.width = 10)
# options(repr.plot.height = 35)

# png('./output/metacell_model/figs/differential_gene_modules_nsc_ipc_boxplot.png', h = 3500, w = 750)
# par(mfrow = c(length(mod_in), 1), mar = c(0.5,6,0.5,0.5), cex.lab = 2)
# # ctout <- c('NSC', 'IPC', 'IPC_cyc')
# ctout <- c()
# inds_mc <- cust_mc_ord_st[!(names(cust_mc_ord_st) %in% ctout)]
# ctv <- factor(mcmd$cell_type[inds_mc], levels = cust_st_ord[!(cust_st_ord %in% ctout)])
# tbt <- sapply(names(mod_in), function(md) {
#     if (md != tail(mod_in, 1)) {
#         par(mar = c(0.5,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], xaxt = 'n', ylab = md)
#     } else {
#         par(las = 2, mar = c(3,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], ylab = md)
#     }
# })
# dev.off()

# png('./output/metacell_model/figs/differential_gene_modules_nsc_ipc_scatter.png', h = 3500, w = 1750)
# par(mfrow = c(length(mod_in), 4), cex.lab = 2, cex.axis = 2)
# inds_nsc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'Astrocytes', 'OPCs')]
# inds_ipc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('IPC', 'IPC_cyc')]
# inds_imm <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('iCPN/CfuPN', 'iCfuPN', 'iCPN_early', 'iCPN_late')]
# inds_neu <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('CthPN', 'SCPN','CPN_L2-3', 'CPN_L5_6')]
# # inds_gli <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('Astrocytes', 'OPCs')]
# legc_mrm <- mmc_gm_avg_module - rowMeans(mmc_gm_avg_module)
# tvt <- sapply(names(mod_in), function(md) {
#     mal <- max(abs(legc_mrm[md,]))
#     if (mal > 3) {
#         ylimi <- c(-mal, mal)
#     } else {ylimi <- c(-3,3)}
#         if (md == tail(names(mod_in), 1)) {
#         par(mar = c(4,0.2,0.2,0.2))
#     } else {par(mar = c(0.2,0.2,0.2,0.2))}
#     # plot(mcmd$mean_day[inds_gli], legc_mrm[md,inds_gli], pch = 16, cex = 2, col = mcmd$color[inds_gli], 
#     #      xlim = c(15.5,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     # grid(col = 'pink')
#     plot(mcmd$mean_day[inds_nsc], legc_mrm[md,inds_nsc], pch = 16, cex = 2, col = mcmd$color[inds_nsc], 
#          xlim = c(13,18), ylim = ylimi, ylab = md, xaxt ='n', xlab = '')
#     grid(col = 'pink')
#     plot(mcmd$mean_day[inds_ipc], legc_mrm[md,inds_ipc], pch = 16, cex = 2, col = mcmd$color[inds_ipc], 
#          xlim = c(13,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     grid(col = 'pink')
#     plot(mcmd$mean_day[inds_imm], legc_mrm[md,inds_imm], pch = 16, cex = 2, col = mcmd$color[inds_imm], 
#          xlim = c(13,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     grid(col = 'pink')
#     plot(mcmd$mean_day[inds_neu], legc_mrm[md,inds_neu], pch = 16, cex = 2, col = mcmd$color[inds_neu], 
#          xlim = c(13,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
#     grid(col = 'pink')
# })
# dev.off()

mod_out <- grep("AY036|clspn|Olig1|tfap2d|xist|rik|apoe|cenpf|dbi|gm35|id3|eomes|gucy1a1|mafb|neurog2|satb2|smc2|vim", rownames(mmc_gm_avg_module), v=T, ign = T)
mod_out

mod_out <- mod_out[order(apply(mmc_gm_avg_module[mod_out,cust_mc_ord_st] - rowMeans(mmc_gm_avg_module[mod_out,]), 
                               1, function(x) sum(ifelse(x < 0, 0, x)*(1:length(x)))/sum(ifelse(x < 0, 0, x))))]
mod_out

mod_in <- rownames(mmc_gm_avg_module)[!(rownames(mmc_gm_avg_module) %in% mod_out)]
mod_in

mod_in <- mod_in[order(apply(mmc_gm_avg_module[mod_in,cust_mc_ord_st] - rowMeans(mmc_gm_avg_module[mod_in,]), 
                             1, function(x) sum(ifelse(x < 0, 0, x)*(1:length(x)))/sum(ifelse(x < 0, 0, x))))]
mod_in

# mod_in_ord1 <- sapply(c('lhx2', 'nrgn', 'ptpr', 'igfbp', 'lhx9', 'nhlh', 'sema3'), function(x) grep(x, mod_in, v=T, ign=T))
# mod_in <- c(mod_in_ord1, mod_in[!(mod_in %in% mod_in_ord1)])
# mod_in

png('./output/metacell_model/figs/non_differential_gene_modules_nsc_ipc_scatter.png', h = 4000*length(mod_in)/(length(mod_in) + length(mod_out)), w = 750)
par(mfrow = c(length(mod_out), 4), cex.lab = 2, cex.axis = 2)
inds_nsc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'Astrocytes', 'OPCs')]
inds_ipc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('IPC', 'IPC_cyc')]
inds_imm <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('iCPN/CfuPN', 'iCfuPN', 'iCPN_early', 'iCPN_late')]
inds_neu <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('CthPN', 'SCPN','CPN_L2-3', 'CPN_L5_6')]
# inds_gli <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('Astrocytes', 'OPCs')]
legc_mrm <- mmc_gm_avg_module - rowMeans(mmc_gm_avg_module)
tvt <- sapply(mod_out, function(md) {
    mal <- max(abs(legc_mrm[md,]))
    if (mal > 3) {
        ylimi <- c(-mal, mal)
    } else {ylimi <- c(-3,3)}
    if (md == tail(mod_in, 1)) {
        mari <- c(4,0.2,0.2,0.2)
        xaxti = 's'
    } else {
        mari <- c(0.2,0.2,0.2,0.2)
        xaxti = 'n'
    }
    # plot(mcmd$mean_day[inds_gli], legc_mrm[md,inds_gli], pch = 16, cex = 2, col = mcmd$color[inds_gli], 
    #      xlim = c(15.5,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
    # grid(col = 'pink')
    tmpm <- mari[[2]]
    mari[[2]] <- 6
    par(mar = mari)
    plot(mcmd$mean_day[inds_nsc], legc_mrm[md,inds_nsc], pch = 16, cex = 2, col = mcmd$color[inds_nsc], 
         xlim = c(13,18), ylim = ylimi, ylab = md, xaxt =xaxti, xlab = '')
    grid(col = 'pink')
    mari[[2]] <- tmpm
    par(mar = mari)
    plot(mcmd$mean_day[inds_ipc], legc_mrm[md,inds_ipc], pch = 16, cex = 2, col = mcmd$color[inds_ipc], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
    plot(mcmd$mean_day[inds_imm], legc_mrm[md,inds_imm], pch = 16, cex = 2, col = mcmd$color[inds_imm], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
    plot(mcmd$mean_day[inds_neu], legc_mrm[md,inds_neu], pch = 16, cex = 2, col = mcmd$color[inds_neu], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
})
dev.off()

png('./output/metacell_model/figs/differential_gene_modules_nsc_ipc_scatter.png', h = 4000*length(mod_in)/(length(mod_in) + length(mod_out)), w = 750)
par(mfrow = c(length(mod_in), 4), cex.lab = 2, cex.axis = 2)
inds_nsc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('NSC', 'Astrocytes', 'OPCs')]
inds_ipc <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('IPC', 'IPC_cyc')]
inds_imm <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('iCPN/CfuPN', 'iCfuPN', 'iCPN_early', 'iCPN_late')]
inds_neu <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('CthPN', 'SCPN','CPN_L2-3', 'CPN_L5_6')]
# inds_gli <- cust_mc_ord_st[names(cust_mc_ord_st) %in% c('Astrocytes', 'OPCs')]
legc_mrm <- mmc_gm_avg_module - rowMeans(mmc_gm_avg_module)
tvt <- sapply(mod_in, function(md) {
    mal <- max(abs(legc_mrm[md,]))
    if (mal > 3) {
        ylimi <- c(-mal, mal)
    } else {ylimi <- c(-3,3)}
    if (md == tail(mod_in, 1)) {
        mari <- c(4,0.2,0.2,0.2)
        xaxti = 's'
    } else {
        mari <- c(0.2,0.2,0.2,0.2)
        xaxti = 'n'
    }
    # plot(mcmd$mean_day[inds_gli], legc_mrm[md,inds_gli], pch = 16, cex = 2, col = mcmd$color[inds_gli], 
    #      xlim = c(15.5,18), ylim = ylimi, xaxt ='n', yaxt = 'n', xlab = '')
    # grid(col = 'pink')
    tmpm <- mari[[2]]
    mari[[2]] <- 6
    par(mar = mari)
    plot(mcmd$mean_day[inds_nsc], legc_mrm[md,inds_nsc], pch = 16, cex = 2, col = mcmd$color[inds_nsc], 
         xlim = c(13,18), ylim = ylimi, ylab = md, xaxt =xaxti, xlab = '')
    grid(col = 'pink')
    mari[[2]] <- tmpm
    par(mar = mari)
    plot(mcmd$mean_day[inds_ipc], legc_mrm[md,inds_ipc], pch = 16, cex = 2, col = mcmd$color[inds_ipc], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
    plot(mcmd$mean_day[inds_imm], legc_mrm[md,inds_imm], pch = 16, cex = 2, col = mcmd$color[inds_imm], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
    plot(mcmd$mean_day[inds_neu], legc_mrm[md,inds_neu], pch = 16, cex = 2, col = mcmd$color[inds_neu], 
         xlim = c(13,18), ylim = ylimi, xaxt =xaxti, yaxt = 'n', xlab = '')
    grid(col = 'pink')
})
dev.off()

genes_mods_in_sort <- setNames(unlist(sapply(mod_in, function(md) rep(md, length(mmc_gm$gene[mmc_gm$module == md])))), 
                                             unlist(sapply(mod_in, function(md) mmc_gm$gene[mmc_gm$module == md])))
# genes_mods_in_sort

genes_mods_out_sort <- setNames(unlist(sapply(mod_out, function(md) rep(md, length(mmc_gm$gene[mmc_gm$module == md])))), 
                                             unlist(sapply(mod_out, function(md) mmc_gm$gene[mmc_gm$module == md])))
# genes_mods_in_sort

tfs_mods_in_sort <- genes_mods_in_sort[names(genes_mods_in_sort) %in% tfs]
# genes_mods_in_sort

tfs_mods_in_sort

tfs_mods_out_sort <- genes_mods_out_sort[names(genes_mods_out_sort) %in% tfs]
# genes_mods_in_sort

tfs_mods_out_sort

# p_mod_in_genes <- pheatmap(legc[names(genes_mods_in_sort),cust_mc_ord_st] - rowMeans(legc[names(genes_mods_in_sort),]), silent = T, fontsize_row = 5,
#                            gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
#                            cluster_rows = F, show_colnames = F, legend = T,treeheight_row = F, 
#          # annotation_row = rah,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D',
#          col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         breaks = seq(-3,3,l=100))
# save_pheatmap_png(p_mod_in_genes, './output/metacell_model/figs/mod_in_by_gene.png', h = 6000*length(mod_in)/(length(mod_in) + length(mod_out)), w = 1100)

library(ComplexHeatmap)

top_ha <- columnAnnotation(cell_type = anno_simple(x = mcmd$cell_type[cust_mc_ord_st], col = tibble::deframe(color_key)),
                          mean_day = anno_simple(x = mcmd$mean_day[cust_mc_ord_st], 
                                                 col = circlize::colorRamp2(breaks = seq(13,18,l=6), 
                                                            colors = c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))))

p_mod_in_genes <- ComplexHeatmap::Heatmap(legc[names(genes_mods_in_sort),cust_mc_ord_st] - rowMeans(legc[names(genes_mods_in_sort),]), show_heatmap_legend = F,
                                          row_names_gp = gpar(fontsize = 10),
                                          # row_title = unique(genes_mods_in_sort),
                                          split = factor(genes_mods_in_sort, levels = mod_in),
                                          # fontsize_row = 5,
                           # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
                           cluster_rows = F,
                                          show_column_names = F, 
                                          # legend = T,
         # annotation_row = rah,
         cluster_columns = F, top_annotation = top_ha,
         col = circlize::colorRamp2(breaks = seq(-4,4,l=3), colors = c('blue3', 'white', 'red3')),
        # breaks = seq(-3,3,l=100)
                                         )

png( './output/metacell_model/figs/mod_in_by_gene.png', h = 6000*length(mod_in)/(length(mod_in) + length(mod_out)), w = 1100)
draw(p_mod_in_genes)
dev.off()

p_mod_in_tfs <- ComplexHeatmap::Heatmap(legc[names(tfs_mods_in_sort),cust_mc_ord_st] - rowMeans(legc[names(tfs_mods_in_sort),]), show_heatmap_legend = F, row_title_gp = gpar(fontsize = 25),
                                          row_names_gp = gpar(fontsize = 20),
                                          # row_title = unique(genes_mods_in_sort),
                                          split = factor(tfs_mods_in_sort, levels = mod_in),
                                          # fontsize_row = 5,
                           # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
                           cluster_rows = F,
                                          show_column_names = F, 
                                          # legend = T,
         # annotation_row = rah,
         cluster_columns = F, top_annotation = top_ha,
         col = circlize::colorRamp2(breaks = seq(-4,4,l=3), colors = c('blue3', 'white', 'red3')),
        # breaks = seq(-3,3,l=100)
                                         )

png( './output/metacell_model/figs/mod_in_by_tf.png', h = 1500, w = 850)
draw(p_mod_in_tfs)
dev.off()

p_mod_out_tfs <- ComplexHeatmap::Heatmap(legc[names(tfs_mods_out_sort),cust_mc_ord_st] - rowMeans(legc[names(tfs_mods_out_sort),]), show_heatmap_legend = F, row_title_gp = gpar(fontsize = 25),
                                          row_names_gp = gpar(fontsize = 20),
                                          # row_title = unique(genes_mods_in_sort),
                                          split = factor(tfs_mods_out_sort, levels = mod_out),
                                          # fontsize_row = 5,
                           # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
                           cluster_rows = F,
                                          show_column_names = F, 
                                          # legend = T,
         # annotation_row = rah,
         cluster_columns = F, top_annotation = top_ha,
         col = circlize::colorRamp2(breaks = seq(-4,4,l=3), colors = c('blue3', 'white', 'red3')),
        # breaks = seq(-3,3,l=100)
                                         )

png( './output/metacell_model/figs/mod_out_by_tf.png', h = 2750, w = 850)
draw(p_mod_out_tfs)
dev.off()

p_mod_out_genes <- ComplexHeatmap::Heatmap(legc[names(genes_mods_out_sort),cust_mc_ord_st] - rowMeans(legc[names(genes_mods_out_sort),]), show_heatmap_legend = F,
                                          row_names_gp = gpar(fontsize = 8),
                                          # row_title = unique(genes_mods_in_sort),
                                          split = factor(genes_mods_out_sort, levels = mod_out),
                                          # fontsize_row = 5,
                           # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
                           cluster_rows = F,
                                          show_column_names = F, 
                                          # legend = T,
         # annotation_row = rah,
         cluster_columns = F, top_annotation = top_ha,
         col = circlize::colorRamp2(breaks = seq(-4,4,l=3), colors = c('blue3', 'white', 'red3')),
        # breaks = seq(-3,3,l=100)
                                         )

png( './output/metacell_model/figs/mod_out_by_gene.png', h = 8000*length(mod_out)/(length(mod_in) + length(mod_out)), w = 1400)
draw(p_mod_out_genes)
dev.off()

mmc_gm$gene[grep('neurog2', mmc_gm$module, ign=T)]

gi <- mmc_gm$gene[grep('neurog2|eomes', mmc_gm$module, ign=T)]
pheatmap::pheatmap(legc[gi,cust_mc_ord_st] - rowMeans(legc[gi,]), annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F,
                                      # show_heatmap_legend = F,
                                          # row_names_gp = gpar(fontsize = 8),
                                          # row_title = unique(genes_mods_in_sort),
                                          # split = factor(genes_mods_out_sort, levels = mod_out),
                                          # fontsize_row = 5,
                           # gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
                           cluster_rows = T,
                                          # show_column_names = F, 
                                          # legend = T,
         # annotation_row = rah,
         cluster_columns = F, 
                                      # top_annotation = top_ha,
         col = colorRampPalette(colors = c('blue3', 'white', 'red3'))(100),
        # breaks = seq(-3,3,l=100)
                                         )

# p_mod_in_genes <- pheatmap(legc[names(genes_mods_in_sort),cust_mc_ord_st] - rowMeans(legc[names(genes_mods_in_sort),]), silent = T, fontsize_row = 5,
#                            gaps_row = which(diff(as.numeric(factor(genes_mods_in_sort))) == 1),
#                            cluster_rows = F, show_colnames = F, legend = T,treeheight_row = F, 
#          # annotation_row = rah,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D',
#          col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         breaks = seq(-3,3,l=100))
# save_pheatmap_png(p_mod_in_genes, './output/metacell_model/figs/mod_in_by_gene.png', h = 6000*length(mod_in)/(length(mod_in) + length(mod_out)), w = 1100)

# png('./output/metacell_model/figs/non_differential_gene_modules_nsc_ipc_boxplot.png', h = 3500, w = 750)
# par(mfrow = c(length(mod_out), 1), mar = c(0.5,6,0.5,0.5), cex.lab = 2)
# # ctout <- c('NSC', 'IPC', 'IPC_cyc')
# ctout <- c()
# inds_mc <- cust_mc_ord_st[!(names(cust_mc_ord_st) %in% ctout)]
# ctv <- factor(mcmd$cell_type[inds_mc], levels = cust_st_ord[!(cust_st_ord %in% ctout)])
# tbt <- sapply(names(mod_out), function(md) {
#     if (md != tail(mod_out, 1)) {
#         par(mar = c(0.5,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], xaxt = 'n', ylab = md)
#     } else {
#         par(las = 2, mar = c(3,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], ylab = md)
#     }
# })
# dev.off()

options(repr.plot.width = 8)
options(repr.plot.height = 35)

# par(mfrow = c(length(mod_out), 1), mar = c(0.5,6,0.5,0.5), cex.lab = 2)
# # ctout <- c('NSC', 'IPC', 'IPC_cyc')
# ctout <- c()
# inds_mc <- cust_mc_ord_st[!(names(cust_mc_ord_st) %in% ctout)]
# ctv <- factor(mcmd$cell_type[inds_mc], levels = cust_st_ord[!(cust_st_ord %in% ctout)])
# tbt <- sapply(names(mod_out), function(md) {
#     mal <- max(abs(legc_mrm[md,inds_mc]))
#     if (mal > 3) {
#         ylimi <- c(-mal, mal)
#     } else {ylimi <- c(-3,3)}
#     if (md != tail(mod_in, 1)) {
#         par(mar = c(0.5,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], xaxt = 'n', ylab = md, ylim = ylimi)
#     } else {
#         par(las = 2, mar = c(3,6,0.5,0.5))
#         boxplot(legc_mrm[md,inds_mc] ~ ctv, col = color_key$color[match(unique(mcmd$cell_type[inds_mc]), color_key$cell_type)], ylab = md, ylim = ylimi)
#     }
#     grid(col = 'pink', lwd = 1, lty = 2)
# })

pheatmap(mmc_gm_avg_module[,cust_mc_ord_st] - rowMeans(mmc_gm_avg_module),cluster_rows = T, show_colnames = F, legend = T,treeheight_row = F, 
         # annotation_row = rah,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = T, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))





nsc_pk_cl <- c(2:4,17, 21,28)

# inds <- which(km_nsc$cluster %in% which(rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n) >= 2))
inds <- which(km_nsc$cluster %in% nsc_pk_cl)
inds <- inds[order(km_nsc$cluster[inds])]
rah <- as.data.frame(list(cluster = as.character(km_nsc$cluster[inds])))
rownames(rah) <- rownames(mat_nsc)[inds]
ann_colors[['cluster']] <- setNames((chameleon::distinct_colors(n = length(unique(rah$cluster))))$name, unique(rah$cluster))

pheatmap(mat_nsc[inds,],cluster_rows = F, show_colnames = F, legend = T,treeheight_row = F, annotation_row = rah,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = T, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

options(repr.plot.width = 12)
options(repr.plot.height = 20)

mcs_ipc <- cust_mc_ord_st[cust_mc_ord_st %in% which(mcmd$cell_type %in% c('Astrocytes','NSC', 'IPC', 'IPC_cyc', 'iCPN_early', 'iCPN/CfuPN', 'CPN_L2-3', 'CthPN'))]

library(matrixStats)

# inds <- which(km_nsc$cluster %in% which(rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n) >= 2))
mata <- do.call('rbind', lapply(nsc_pk_cl, function(i) {
    inds <- which(km_nsc$cluster %in% i)
    # inds <- inds[order(km_nsc$cluster[inds])]
    inds <- inds[rowMaxs(legc[inds,mcs_ipc]) - rowMins(legc[inds,mcs_ipc]) >= 4]
    hc <- hclust(dist(legc[inds,mcs_ipc] - rowMeans(legc[inds,mcs_ipc])), method = 'ward.D2')
    mati <- legc[inds[hc$order],]
    return(mati)
}))

options(repr.plot.width = 8)
options(repr.plot.height = 13)

dim(mata)

dim(mat_nsc[inds,])

pheatmap(mata[,mcs_ipc] - rowMeans(mata),cluster_rows = F, show_colnames = F, legend = F,treeheight_row = F, annotation_row = rah,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = T, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

options(repr.plot.width = 18)
options(repr.plot.height = 18)

pheatmap(mata[,mcs_ipc] - rowMeans(mata),cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F, annotation_row = rah,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = T, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

# # inds <- which(km_nsc$cluster %in% which(rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n) >= 2))
# inds <- which(km_nsc$cluster %in% c(2:4,21))
# inds <- inds[order(km_nsc$cluster[inds])]
# pheatmap(legc[inds,mcs_ipc],cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F, annotation_row = rah,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = T, clustering_method = 'ward.D2',
#          col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100)
#         #  col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         # breaks = seq(-3,3,l=100)
#         )

km_nsc$size

pheatmap(mat_avg_km_n[order(km_nsc$size),],cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

pheatmap(mat_all_avg_km_n[order(km_nsc$size),mcs_ipc] - rowMeans(mat_all_avg_km_n[order(km_nsc$size),mcs_ipc]),cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

hist(rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n))

plot(rowSds(mat_avg_km_n), rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n))
abline(0,1,col='red')

which(rowMaxs(mat_avg_km_n) - rowMins(mat_avg_km_n) >= 2)

inds <- which(km_nsc$cluster %in% c(13:15,1,26,9:11))
inds <- inds[order(km_nsc$cluster[inds])]
pheatmap(mat_nsc[inds,],cluster_rows = F, show_colnames = F, legend = F,treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

pheatmap(t(apply(mat_avg_km_n, 1, zoo::rollmean, k = 20, na.pad = T)),cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

dim(mat_nsc)
km_m_n <- tglkmeans::TGL_kmeans(mat_nsc - rowMeans(mat_nsc), k = 8)

km_m_n <- tglkmeans::TGL_kmeans(mat_nsc - rowMeans(mat_nsc), k = 8)

mat_avg_km_n <- tgs_matrix_tapply(t(mat_nsc - rowMeans(mat_nsc)), km_m_n$cluster, mean)

options(repr.plot.width = 8)
options(repr.plot.height = 5)

p_nsc_clusters <- pheatmap(mat_avg_km_n,cluster_rows = T, show_colnames = F, legend = F,treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

save_pheatmap_png(p_nsc_clusters, './output/metacell_model/figs/nsc_clusters.png', h = 600, w = 1000)

options(repr.plot.width = 8)
options(repr.plot.height = 15)

# p_nsc_genes <- pheatmap(mat_nsc[order(km_m_n$cluster),] - rowMeans(mat_nsc[order(km_m_n$cluster),]), cluster_rows = F, show_colnames = F, legend = F,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
#          col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         breaks = seq(-3,3,l=100)
# )
# save_pheatmap_png(p_nsc_genes, './output/metacell_model/figs/nsc_genes.png', h =1800, w = 1000)

p_nsc_genes <- pheatmap(mat_nsc[order(match(km_m_n$cluster, p_nsc_clusters$tree_row$order)),] - 
                        rowMeans(mat_nsc[order(match(km_m_n$cluster, p_nsc_clusters$tree_row$order)),]), 
                        fontsize_row = 6, cluster_rows = F, show_colnames = F, legend = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100)
)
save_pheatmap_png(p_nsc_genes, './output/metacell_model/figs/nsc_genes.png', h =3000, w = 1200)



# p_nsc_genes <- pheatmap(mall - rowMeans(mall), cluster_rows = F, show_colnames = F, legend = F, fontsize_row = 8,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
#         # col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
#         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         breaks = seq(-3,3,l=100)
#                                   )
# save_pheatmap_png(p_nsc_genes, './output/metacell_model/nsc_genes.png', h =2400, w = 1100)

# mall <- legc[which(cor_gene_legc_w_md_in_st[,'NSC'] <= -0.7 & st_legc[,'IPC'] >= -15),mcs_nsc]
# p_nsc_genes <- pheatmap(mall, cluster_rows = T, show_colnames = F, legend = F,
#          cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
#         col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100),
#         # col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
#         # breaks = seq(-3,3,l=100)
#                                   )


mcs_ipc <- cust_mc_ord_st[cust_mc_ord_st %in% which(mcmd$cell_type %in% c('Astrocytes','NSC', 'IPC', 'IPC_cyc', 'iCPN_early', 'iCPN/CfuPN', 'CPN_L2-3', 'CthPN'))]
x <- sort(cor_gene_legc_w_md_in_st[,'IPC'])
y1 <- names(head(x, 150))
y2 <- names(tail(x, 150))
x1 <- setdiff(names(ipc_spec_genes), union(c(y1, y2), cc_genes))
m1 <- legc[x1,mcs_ipc]
hc1 <- hclust(dist(m1 - rowMeans(m1)), method = 'ward.D2')
m2 <- legc[y1,mcs_ipc]
hc2 <- hclust(dist(m2 - rowMeans(m2)), method = 'ward.D2')
m3 <- legc[y2,mcs_ipc]
hc3 <- hclust(dist(m3 - rowMeans(m3)), method = 'ward.D2')
mall <- rbind(m1[hc1$order,],m2[hc2$order,],m3[hc3$order,])

mat_ipc <- legc[unique(c(x1, y1, y2)),mcs_ipc]
dim(mat_nsc)

km_m_i <- tglkmeans::TGL_kmeans(mat_ipc - rowMeans(mat_ipc), k = 14)

mat_avg_km_i <- tgs_matrix_tapply(t(mat_ipc - rowMeans(mat_ipc)), km_m_i$cluster, mean)

options(repr.plot.width = 8)
options(repr.plot.height = 5)

p_ipc_clusters <- pheatmap(mat_avg_km_i,cluster_rows = T, show_colnames = F, legend = F, treeheight_row = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

save_pheatmap_png(p_ipc_clusters, './output/metacell_model/figs/ipc_clusters.png', h = 600, w = 1000)

options(repr.plot.width = 8)
options(repr.plot.height = 15)

p_ipc_genes <- pheatmap::pheatmap(mat_ipc[order(match(km_m_i$cluster, p_ipc_clusters$tree_row$order)),] - 
                                  rowMeans(mat_ipc[order(match(km_m_i$cluster, p_ipc_clusters$tree_row$order)),]), fontsize_row = 6,
                                  cluster_rows = F, show_colnames = F, legend = F,
         cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F, clustering_method = 'ward.D2',
         col = colorRampPalette(c('blue4', 'white', 'red4'))(100),
        breaks = seq(-3,3,l=100))

save_pheatmap_png(p_ipc_genes, './output/metacell_model/figs/ipc_genes.png', h =3000, w = 1200)

cor_tfs_scpn_cthpn <- tgs_cor(t(legc[union(names(scpn_tfs[scpn_tfs > 0.1]),names(cthpn_tfs[cthpn_tfs > 0.1])),which(mcmd$cell_type %in% c('SCPN', 'CthPN'))]), spearman = T)

diag(cor_tfs_scpn_cthpn) <- 0

options(repr.plot.width = 18)
options(repr.plot.height = 18)

sdrow_cor_tf_scpn_cthpn <- matrixStats::rowSds(cor_tfs_scpn_cthpn)
genes_cor <- rownames(cor_tfs_scpn_cthpn)[which(sdrow_cor_tf_scpn_cthpn >= quantile(sdrow_cor_tf_scpn_cthpn, 0.8))]

all_genes_to_plot <- lapply(list(nsc_tfs, ipc_tfs, late_neuro_tfs, scpn_tfs, cthpn_tfs, l56_tfs, l23_tfs), function(l) head(names(l), 40))
for (i in 2:length(all_genes_to_plot)) {
    all_genes_to_plot[[i]] <- all_genes_to_plot[[i]][!(all_genes_to_plot[[i]] %in% unlist(sapply(all_genes_to_plot[1:(i-1)], head)))]
}
all_genes_select <- lapply(all_genes_to_plot, function(x) head(x,2))

names(all_genes_to_plot) <- c('NSC', 'IPC', 'Mature neuron', 'SCPN', 'CthPN', 'CPN_L5_6', 'CPN_L2-3')

all_genes_to_plot

all_genes_vec <- unique(do.call('c', purrr::imap(all_genes_to_plot, .f =function(.x,.y) setNames(.x, rep(.y, length(.x))))))
names(all_genes_vec) <- gsub('\\..*', '', names(all_genes_vec))
all_genes_vec

ag_wm <- apply(legc[all_genes_vec,cust_mc_ord_st], 1, which.max)

ag_com <- apply(legc[all_genes_vec,cust_mc_ord_st], 1, function(x) sum(x*1:length(x))/sum(x))

library(pheatmap)

pheatmap(legc[all_genes_vec[order(ag_com, decreasing = T)],cust_mc_ord_st], cluster_rows = F, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors)

options(repr.plot.width = 18)
options(repr.plot.height = 18)

pheatmap(legc[all_genes_vec[order(ag_wm, decreasing = F)],cust_mc_ord_st], 
         cluster_rows = F, cluster_cols = F, 
         annotation_col = col_annot, annotation_colors = ann_colors)

pheatmap(legc[all_genes_vec[order(ag_wm)],cust_mc_ord_st], cluster_rows = T, cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors)

library(pheatmap)

options(repr.plot.width = 18)
options(repr.plot.height = 12)

gp <- unique(unlist(all_genes_to_plot[4:7]))

# gpo <- gp[order(apply(legc[gp,cust_mc_ord_st], 1, function(x) sum(x*(1:length(x)))/sum(x)))]
gpo <- gp[hclust(tgs_dist(legc[gp,cust_mc_ord_st]), method = 'ward.D2')$order]

eps <- 1
pltmt <- t(apply(legc[gpo,cust_mc_ord_st], 1, function(x) (eps + x - min(x))/median(eps + x - min(x))))

pltmt <- log2(pltmt)

ord <- order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x)))

clrs <- colorRampPalette(c('blue3', 'white', 'red3'))(100)
brks <- c(seq(min(pltmt), 0, l = 50), seq(1e-5, max(pltmt), l = 51))

p_tf_hm <- pheatmap(pltmt[ord,],annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, cluster_rows = T, color = clrs, breaks =brks, main = 'log2((1 + x - min(x))/median(1 + x - min(x)))', fontsize = 12)

col_annot$mean_day <- mcmd$mean_day
ann_colors[['mean_day']] <- setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green3', 'blue2', 'purple'))(100), seq(13,18,l=100))

p_tf_hm_neurons <- pheatmap(pltmt[ord,which(cust_mc_ord_st %in% which(mcmd$cell_type %in% cts))],annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, cluster_rows = T, color = clrs, breaks =brks, main = 'log2((1 + x - min(x))/median(1 + x - min(x)))', fontsize = 12)

save_pheatmap_png(p_tf_hm, './output/metacell_model/figs/tf_log_fc_heatmap.png', w = 2400, h = 2000, res = 150)

save_pheatmap_png(p_tf_hm_neurons, './output/metacell_model/figs/tf_log_fc_heatmap_neurons.png', w = 1800, h = 1800, res = 150)

options(repr.plot.width = 18)
options(repr.plot.height = 18)

dir.create('./output//metacell_model/figs/')


pltmt = mc@mc_fp[goi,cust_mc_ord_st]
# pltmt = pltmt[order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x))),]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp = pheatmap::pheatmap(pltmt, color = colorRampPalette(pltt)(l), breaks = brks, fontsize = 14, legend = T, annotation_legend = F,
                        cluster_cols = F, cluster_rows = F, 
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors, 
                        show_colnames = F, fontsize_row = 14)
save_pheatmap_png(pp, './output//metacell_model/figs/marker_heatmap.png', w = 2400, h = 1200, res = 200)

goi <- setdiff(unique(unlist(sapply(mg_bon_marks[mg_bon_marks$cell_type %in% unique(mcmd$cell_type),'marks'], function(x) unlist(stringr::str_split(x,','))))),
                                    c('Gfap'))


pltmt = mc@mc_fp[goi,cust_mc_ord_st2]
# pltmt = pltmt[order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x))),]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp = pheatmap::pheatmap(pltmt, color = colorRampPalette(pltt)(l), breaks = brks, fontsize = 14, legend = T, annotation_legend = F,
                        cluster_cols = F, cluster_rows = F, 
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors, 
                        show_colnames = F, fontsize_row = 14)
# save_pheatmap_png(pp, './output//metacell_model/figs/marker_heatmap.png', w = 2400, h = 1200, res = 200)

pltmt = mc@mc_fp[goi,cust_mc_ord_st2]
# pltmt = pltmt[order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x))),]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))

library(ComplexHeatmap)


pltmt = mc@mc_fp[goi,cust_mc_ord_st2]
# pltmt = pltmt[order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x))),]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp = pheatmap::pheatmap(pltmt, color = colorRampPalette(pltt)(l), breaks = brks, fontsize = 14, legend = T, annotation_legend = F,
                        cluster_cols = F, cluster_rows = F, 
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors, 
                        show_colnames = F, fontsize_row = 14)
# save_pheatmap_png(pp, './output//metacell_model/figs/marker_heatmap.png', w = 2400, h = 1200, res = 200)

top_ha <- anno_simple(x = mcmd$cell_type[cust_mc_ord_st2], col = ann_colors[['cell_type']])

md_col_fun <- circlize::colorRamp2(breaks = 13:18, 
                                    colors = rev(c('purple', 'blue', 'green', 'yellow', 'orange', 'red')))
bottom_ha <- anno_simple(x = mcmd$mean_day[cust_mc_ord_st], 
                         col = md_col_fun)
lgd_MD = Legend(title = "mean day", col_fun = md_col_fun, at = 13:18, 
    labels = 13:18)

ch <- Heatmap(pltmt, name = "marker\nheatmap",
              top_annotation = HeatmapAnnotation(`cell type` = top_ha), 
              bottom_annotation = HeatmapAnnotation(`mean day` = bottom_ha, show_legend = T), 
              col = circlize::colorRamp2(breaks = c(minv, 0, maxv/2, maxv), colors = pltt), show_column_names = F, row_names_gp = gpar(fontsize = 14),
                             # legend = T, 
                             # annotation_legend = F,
                        cluster_columns = F, cluster_rows = F
                        # annotation_col = col_annot, 
                        # annotation_colors = ann_colors, 
                        )
# save_pheatmap_png(pp, './output//metacell_model/figs/marker_heatmap.png', w = 2400, h = 1200, res = 200)

png('./output//metacell_model/figs/marker_ComplexHeatmap.png', w = 2400, h =1200, res = 200)
draw(ch, annotation_legend_list = list(lgd_MD))
dev.off()
# draw(ch, annotation_legend_list = list(lgd_MD))

?par

# color_key = color_key[c('Astrocytes','early_NSC?','OPC','NSC','nNSC?','early_nNSC?','IPC','iCfuPN','iCPN_L2-3',
#               'iCPN_L5-6','CPN_L2-3','CPN_L5-6','SCPN','CthPN','Stellate_L4'),]
df = data.frame(color_key[order(match(color_key$cell_type, cust_st_ord), decreasing = F),])

# st_for_plot = df$st
# st_for_plot['iCPN_L2-3'] = 'immature CPN L2-3'
# st_for_plot['iCfuPN'] = 'immature CPN L2-3'


l = nrow(df)
scale_y = 1
# svg('./figs/legend_pl_rev.svg', width = 1200, height = 900)
png('./figs/legend_pl.png', width = 750, height = 2000, res = 250)
par(mar = c(4,1,4,0), bty = 'n')
plot(rep(0.9,l), scale_y*seq(l,1,-1), pch = 16, cex = 5, col = df$color, ylim = c(0,scale_y*l+1),
     xlim = c(-1,60),
    xlab = '', 
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(rep(6,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 1.81, df$cell_type)
dev.off()

dev.off()

source('./scripts/paper_scripts/util.r')

options(repr.plot.width = 14)
options(repr.plot.height = 14)

mc2d_bord_x_r <- 0.75
mc2d_bord_y_d <- 0.4
legend_bord_x_l <- mc2d_bord_x_r*0.8
par(fig = c(mc2d_bord_y_d,mc2d_bord_x_r, 0.02,0.95), mar = rep(0.2,4))
bb <- my_mcell_mc2d_plot('pl_cort')
# svg('./figs/legend_pl.svg', width = 1200, height = 900)
par(fig = c(legend_bord_x,mc2d_bord_x, 0.62,0.95), new = T)
df = data.frame(color_key[order(match(color_key$cell_type, cust_st_ord)),])
l = nrow(df)
scale_y = 2
plot(rep(0.93,l), scale_y*seq(l,1,-1), pch = 16, cex = 2, col = df$color, xlim = c(0.92, 0.97), ylim = c(0.5,scale_y*l+1),
    xlab = '', 
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(rep(0.94,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 2, df$cell_type)



# dev.off()

ct_per_day <- tgs_matrix_tapply(t(matrix(mc_ag, nrow = nrow(mc_ag), ncol = ncol(mc_ag))), mcmd$cell_type, sum)
ct_per_day <- t(t(ct_per_day)/colSums(ct_per_day))

colnames(ct_per_day) <- paste0('E', 13:18)

ct_per_day_n <- ct_per_day/rowSums(ct_per_day)

ra <- as.data.frame(tibble::column_to_rownames(color_key, 'cell_type')) %>% rename(cell_type = color)

ra$cell_type <- rownames(ra)

p_ct_per_day <- pheatmap(ct_per_day[cust_st_ord,], cluster_cols = F, cluster_rows = F, fontsize = 26, 
                         annotation_row = ra, annotation_colors = ann_colors, annotation_legend = F
         # color = colorRampPalette(c('white', 'yellow', 'red3', 'black'))(100)
        )
save_pheatmap_png(p_ct_per_day, './figs/ct_per_day_unnorm.png', h = 1200, w = 900, res = 100)

p_ct_per_day <- pheatmap(ct_per_day_n[cust_st_ord,], cluster_cols = F, cluster_rows = F, fontsize = 26, 
                         annotation_row = ra, annotation_colors = ann_colors, annotation_legend = F
         # color = colorRampPalette(c('white', 'yellow', 'red3', 'black'))(100)
        )
save_pheatmap_png(p_ct_per_day, './figs/ct_per_day.png', h = 1200, w = 900, res = 100)

load('./data/pl_cort_cc_score.rda')

cti <- c('Astrocytes', 'NSC', 'OPCs', 'IPC', 'IPC_cyc')
ct_mc <- which(mcmd$cell_type %in%  cti)

mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"day"])
mc_ag_n = mc_ag/rowSums(mc_ag)
mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)

get_egc_trajectory <- function(mc_final, mcf, t_final = NULL, t_initial = NULL, mc_initial = NULL, ct_initial = NULL, ct_final = NULL) {
    p_mc = rep(0,ncol(mc@e_gc))
    p_mc[mc_final] = mc_ag_c[mc_final,t_final]
    mc_pf_f <- mctnetflow_propogate_from_t(mcf = mcf,t = t_final,mc_p = p_mc)
    mpfpn <- t(t(mc_pf_f$probs)/colSums(mc_pf_f$probs))
    egc_t = mc@e_gc %*% mpfpn
    colnames(egc_t) <- paste0('E', 13:18)
    return(list(egc_t = egc_t, mc_pf_f = mc_pf_f))
}

# plot trajectory (just visualization)
backprop_and_plot_traj = function(mc, mct, mcf, color, net_id, name = NULL, tt, mc_ord = NULL, color_ord = NULL, mc_pf=0.1, w=2000, h=1000, edge_w_scale=5e-4, max_lwd=20){
    mc_p = rep(0,ncol(mc@e_gc))
    mc_p[mc@colors == color] = mct@mc_t[mc@colors == color,tt]
    mc_p = mc_pf*mc_p/sum(mc_p)

    card_prop = mctnetflow_propogate_from_t(mcf = mcf, t = tt, mc_p = mc_p)
    mm_mctnetwork_plot_net(w = w,
        h = h,
        mct_id = net_id,
        flow_id = net_id,
        mc_ord = mc_ord,
        propogate = card_prop$step_m,
        fn = paste0(paste0('figs/flow_tracebacks/',name),'_flow.png'),
#         colors_ordered = color_ord,
        edge_w_scale = edge_w_scale,  # 
        max_lwd = max_lwd)
    return(card_prop)
}

ct_t_df <- expand.grid(c('SCPN', 'CthPN', 'CPN_L5_6', 'CPN_L2-3'), 4:6, stringsAsFactors = F)
colnames(ct_t_df) <- c('cti', 'tti')

dir.create('./figs/flow_tracebacks')

bp_traj <- lapply(1:nrow(ct_t_df), function(n) backprop_and_plot_traj(mc, mct, mcf, name = paste0(ct_t_df$cti[[n]], '_', ct_t_df$tti[[n]]), net_id = 'pl_cort', color = color_key$color[color_key$cell_type == ct_t_df$cti[[n]]], tt = ct_t_df$tti[[n]], mc_ord = cust_mc_ord_st))

names(bp_traj) <- apply(ct_t_df, 1, function(x) paste0(x, collapse = '_'))

apply(ct_t_df, 1, function(x) paste0(x, collapse = '_'))

egc_t_traj <- lapply(bp_traj, function(x) mc@e_gc %*% x$probs)

names(egc_t_traj)

cor_tf_md <- lapply(tail(egc_t_traj, 4), function(x) tgs_cor(t(x[tfs_hi,]), as.matrix(1:ncol(x)), spearman = T))

lapply(cor_tf_md, function(x) c(names(head(x[order(x),])), names(tail(x[order(x),])), 'Bcl11b', 'Satb2', 'Rnd2', 'Sstr2'))

gct_df <- as.data.frame(tidyr::pivot_longer(as.data.frame(lapply(cor_tf_md, function(x) c(names(head(x[order(x),])), names(tail(x[order(x),])), c('Bcl11b', 'Satb2', 'Rnd2', 'Sstr2', 'Mef2c')))), cols = everything(), names_to = 'cell_type', values_to = 'gene'))

gct_df[,1] <- gsub('_6$', '', gct_df[,1])
gct_df[,1] <- gsub('\\.', '-', gct_df[,1])

head(gct_df)

gs <- c('Pou3f2', 'Pou3f1', 'Pou3f3', 'Zeb2', 'Ybx1', 'Zbtb20', 'Nhlh1', 'Nfia', 'Nfib', "Nfix")

cts <- c('SCPN', 'CthPN', 'CPN_L5_6', 'CPN_L2-3')

ct_g_max <- lapply(cts, function(cti) {
    v <- st_legc[tfs_hi,cti] - apply(st_legc[tfs_hi,which(colnames(st_legc) != cti)], 1, max)
    ordi <- order(v, decreasing = T)
    c(head(rownames(st_legc[tfs_hi,])[ordi], 5), tail(rownames(st_legc[tfs_hi,])[ordi], 5))
})
names(ct_g_max) <- cts

all_tfs_to_plot <- unlist(all_genes_select)

all_tfs_to_plot

all_tfs_to_plot[all_tfs_to_plot == 'Bcl11b'] <- 'Fezf2'
all_tfs_to_plot[all_tfs_to_plot == 'Tshz3'] <- 'Nfib'
all_tfs_to_plot[all_tfs_to_plot == 'Rorb'] <- 'Nfia'
all_tfs_to_plot[all_tfs_to_plot == 'Satb2'] <- 'Foxp1'
all_tfs_to_plot[all_tfs_to_plot == 'Neurod2'] <- 'Pou3f2'

all_tfs_to_plot

gct_df <- apply(expand.grid(cts, all_tfs_to_plot), 2, as.character)

# clrs <- c('blue', 'green', 'red')
clrs <- c('red')

unlist(all_tfs_to_plot)

nsc_flow_out <- do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgs_matrix_tapply(x[which(mcmd$cell_type == 'NSC'),], mcmd$cell_type, sum))/sum(colSums(x))))

ipc_flow_out <- do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgs_matrix_tapply(x[which(mcmd$cell_type == 'IPC'),], mcmd$cell_type, sum))/sum(colSums(x))))

ipcc_flow_out <- do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgs_matrix_tapply(x[which(mcmd$cell_type == 'IPC_cyc'),], mcmd$cell_type, sum))/sum(colSums(x))))

icpn_cfupn_flow_out <- do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgs_matrix_tapply(x[which(mcmd$cell_type == 'iCPN/CfuPN'),], mcmd$cell_type, sum))/sum(colSums(x))))

ct_flow_out_ls <- lapply(cust_st_ord, function(cti) do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgs_matrix_tapply(x[which(mcmd$cell_type == cti),], mcmd$cell_type, sum, na.rm = T), na.rm = T)/sum(colSums(x, na.rm = T), na.rm = T))))

ct_flow_out_ls_n <- lapply(ct_flow_out_ls, function(x) {y <- x/rowSums(x, na.rm = T); y[is.na(y)] <- 0; return(y)})
# ipcc_flow_out

names(ct_flow_out_ls_n) <- cust_st_ord

ct_flow_in_ls <- lapply(cust_st_ord, function(cti) do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgs_matrix_tapply(t(x[,which(mcmd$cell_type == cti)]), mcmd$cell_type, sum, na.rm = T), na.rm = T)/sum(rowSums(x, na.rm = T), na.rm = T))))

ct_flow_in_ls_n <- lapply(ct_flow_in_ls, function(x) {y <- x/rowSums(x, na.rm = T); y[is.na(y)] <- 0; return(y)})
# ipcc_flow_out

names(ct_flow_in_ls_n) <- cust_st_ord

length(ct_flow_in_ls)

length(ct_flow_out_ls)

save(ct_flow_in_ls, file = './output/metacell_flow/flows_in_by_ct.rda')

save(ct_flow_out_ls, file = './output/metacell_flow/flows_out_by_ct.rda')

options(repr.plot.width = 10)
options(repr.plot.height = 26)

png('./output/metacell_flow/figs/total_flows_by_ct.png', h = 2600, w = 700)
par(mfrow = c(length(cust_st_ord),2), las = 2, mar = c(1,8,2,1), cex.main = 2, cex.axis = 2, cex.lab = 2)
for (ct in cust_st_ord) {
    i <- match(ct, names(ct_flow_out_ls_n))
# ct <- names(ct_flow_out_ls_n)[[i]]
    xi <- ct_flow_in_ls[[i]]
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,8,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, 
         main = glue::glue('{ct} incoming flows'), ylab = '', xaxt = 'n', xlab = '')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Total flow', line = 6, cex.lab = 2)
    sususu <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
    xi <- ct_flow_out_ls[[i]]
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,8,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, 
         main = glue::glue('{ct} outgoing flows'), ylab = '', xaxt = 'n', xlab = '')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Total flow', line = 6, cex.lab = 2)
    sususu <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
}
dev.off()

png('./output/metacell_flow/figs/relative_flows_by_ct.png', h = 2600, w = 700)
par(mfrow = c(length(cust_st_ord),2), las = 2, mar = c(1,6,2,1), cex.main = 2, cex.axis = 2, cex.lab = 2)
for (ct in cust_st_ord) {
    i <- match(ct, names(ct_flow_out_ls_n))
# ct <- names(ct_flow_out_ls_n)[[i]]
    xi <- ct_flow_in_ls_n[[i]]
    # print(ct)
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,6,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, 
         main = glue::glue('{ct} incoming flows'), ylab = '', xaxt = 'n', xlab = '')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Relative flow', line = 4, cex.lab = 2)
    sususu <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
    xi <- ct_flow_out_ls_n[[i]]
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,6,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, 
         main = glue::glue('{ct} outgoing flows'), ylab = '', xaxt = 'n', xlab = '')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Relative flow', line = 4, cex.lab = 2)
    sususu <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
}
dev.off()



nsc_flow_in <- do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgs_matrix_tapply(t(x[,which(mcmd$cell_type == 'NSC')]), mcmd$cell_type, sum))/sum(rowSums(x))))

ipc_flow_in <- do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgs_matrix_tapply(t(x[,which(mcmd$cell_type == 'IPC')]), mcmd$cell_type, sum))/sum(rowSums(x))))

ipcc_flow_in <- do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgs_matrix_tapply(t(x[,which(mcmd$cell_type == 'IPC_cyc')]), mcmd$cell_type, sum))/sum(rowSums(x))))

astro_flow_in <- do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgs_matrix_tapply(t(x[,which(mcmd$cell_type == 'Astrocytes')]), mcmd$cell_type, sum))/sum(rowSums(x))))

aa <- paste0('E', 13:17)
bb <- paste0('E', 14:18)
xlabs <- apply(cbind(aa, rep('->', length(aa)), bb), 1, paste, collapse = ' ')

nsc_flow_out_n <- nsc_flow_out/rowSums(nsc_flow_out)

ipc_flow_out_n  <- ipc_flow_out/rowSums(ipc_flow_out)
# ipc_flow_out

icpn_cfupn_flow_out_n <- icpn_cfupn_flow_out/rowSums(icpn_cfupn_flow_out)
# ipcc_flow_out

ipcc_flow_out_n <- ipcc_flow_out/rowSums(ipcc_flow_out)
# ipcc_flow_out

nsc_flow_in_n <- nsc_flow_in/rowSums(nsc_flow_in)
# nsc_flow_in

ipc_flow_in_n  <- ipc_flow_in/rowSums(ipc_flow_in)
# ipc_flow_in

ipcc_flow_in_n <- ipcc_flow_in/rowSums(ipcc_flow_in)
# ipcc_flow_in

astro_flow_in_n <- astro_flow_in/rowSums(astro_flow_in)
# ipcc_flow_in

# png('./output/metacell_flow/figs/flows_outgoing_from_nsc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(nsc_flow_out_n)), lwd = 3, 
     main = 'Target types of flows outgoing from NSCs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
sususu <- sapply(colnames(nsc_flow_out_n), function(ct) {
        lines(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, ipc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(ipc_flow_out_n)), lwd = 3, 
     main = 'Target types of flows\noutgoing from IPCs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, ipc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(ipc_flow_out_n), function(ct) {
        lines(1:5, ipc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, ipc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc_cyc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, ipcc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(ipcc_flow_out_n)), lwd = 3, 
     main = 'Target types of flows\noutgoing from IPC_cycs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, ipcc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(ipcc_flow_out_n), function(ct) {
        lines(1:5, ipcc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, ipcc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc_cyc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, icpn_cfupn_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(icpn_cfupn_flow_out_n)), lwd = 3, 
     main = 'Target types of flows\noutgoing from iCPN/CfuPN', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, icpn_cfupn_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(icpn_cfupn_flow_out_n), function(ct) {
        lines(1:5, icpn_cfupn_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, icpn_cfupn_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_nsc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, nsc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(nsc_flow_in_n)), lwd = 3, 
     main = 'Source types of flows\nincoming to NSCs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
sususu <- sapply(colnames(nsc_flow_in_n), function(ct) {
        lines(1:5, nsc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, nsc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, ipc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(ipc_flow_in_n)), lwd = 3, 
     main = 'Source types of flows\nincoming to IPCs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, ipc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(ipc_flow_in_n), function(ct) {
        lines(1:5, ipc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, ipc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc_cyc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, ipcc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(ipcc_flow_in_n)), lwd = 3, 
     main = 'Source types of flows\nincoming to IPC_cycs', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, ipcc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(ipcc_flow_out_n), function(ct) {
        lines(1:5, ipcc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, ipcc_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

astro_flow_in_n[is.na(astro_flow_in_n)] <- 0

length(intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'Astrocytes')]), rownames(mat@cell_metadata[mat@cell_metadata$day == 'E16',])))

table(mcmd$cell_type[mc@mc], mat@cell_metadata[names(mc@mc),'day'])

# png('./output/metacell_flow/figs/flows_outgoing_from_ipc_cyc.png', h = 800, w = 800, res = 100)
par(las = 2, mar = c(12,8,5,5), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, astro_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(astro_flow_in_n)), lwd = 3, 
     main = 'Source types of flows\nincoming to Astrocytes', ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
points(1:5, astro_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
sususu <- sapply(colnames(astro_flow_in_n), function(ct) {
        lines(1:5, astro_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, astro_flow_in_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
# dev.off()

png('./figs/tf_tracebacks/ct_tfs_in_terminal_cts_2.png', h = 4000, w = 3500, res = 150)
par(mfrow = c(14,4), cex.main = 1, mar = c(1,1,1,1), cex.axis = 3)
vvv <- sapply(1:nrow(gct_df), function(n) {
    g <- gct_df[n,2] ##'Gadd45g'
    cti <- gct_df[n,1] ##'CPN_L2-3'
#     print(g)
#     print(cti)
#     ctt_traj_ls <- lapply(4:6, function(tf) {
    ctt_traj_ls <- lapply(6, function(tf) {
        mcsi <- which(mc_ag_cn[,tf] >= 0.4 & mcmd$cell_type == cti)
#         print(mcmd$cell_type[mcsi])
#         print(mcmd$mean_day[mcsi])
        return(get_egc_trajectory(mc_final = mcsi, mcf = mcf, t_final = tf, t_initial = 1))
    })
    ctt_traj <- lapply(ctt_traj_ls, function(x) x$egc_t)
#     print(head(ctt_traj[[1]]))
#     print(cti)
#     print(g)
#     print(ctt_traj[[1]][g,])
#     print(ctt_traj_ls[[1]]$step_m[1:4,1:4])
    mcs_in <- lapply(ctt_traj_ls, function(x) lapply(x$mc_pf_f$step_m, function(y) which(rowSums(y) > quantile(rowSums(y)[rowSums(y) > 0], 0.2))))
    maxy <- max(sapply(seq_along(ctt_traj), function(i) max(log2(1e-5+ctt_traj[[i]][g,]))))
    print(maxy)
#     print(mcs_in)
#     png(glue::glue('./figs/tf_tracebacks/{cti}_{g}.png'))
#     ctt_all <- log2(1e-5+do.call('rbind', lapply(ctt_traj, function(x) x[g,])))
#     print(ctt_all)
    if (maxy > -10) {
        ylm <- c(-17, maxy)
    } else {
        ylm <- c(-17, -10)
    }
    if (cti == 'SCPN') {
        yaxti <- NULL
        mari = c(1,4,1,1)
    } else {
        yaxti <- 'n'
        mari = c(1,1,1,1)
    }
    if (!(n %in% (nrow(gct_df)-3):(nrow(gct_df)))) {
        par(mar = mari)
        plot(0, xlim = c(13,18), ylim = ylm, col = 'white', ylab = '', xlab = 'day', main = paste(cti, g), xaxt = 'n', yaxt = yaxti)
    } else {
        mari[[1]] <- 3
        par(mar = mari)
        plot(0, xlim = c(13,18), ylim = ylm, col = 'white', ylab = '', xlab = 'day', main = paste(cti, g), yaxt = yaxti)
    }
    # if (!(n %in% (nrow(gct_df)-3):(nrow(gct_df)))) {
    #     plot(0, xlim = c(13,18), ylim = c(-16.7,-7), col = 'white', ylab = 'log2 RNA', xlab = 'day', main = paste(cti, g), xaxt = 'n')
    # } else {
    #     plot(0, xlim = c(13,18), ylim = c(-16.7,-7), col = 'white', ylab = 'log2 RNA', xlab = 'day', main = paste(cti, g))
    # }
    sapply(seq_along(ctt_traj), function(i) {
#         print(head(ctt_traj[[i]]))
        
#         points(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
#         lines(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
        sapply(seq_along(mcs_in[[i]]), function(j) {
#             print(mcs_in[[i]][[j]])
#             points(rep(12+j, length(mcs_in[[i]][[j]]))+runif(n = length(mcs_in[[i]][[j]]), min = -0.1,max = 0.1), legc[g,mcs_in[[i]][[j]]], col = clrs[[i]])
#             points(mcmd$mean_day[mcs_in[[i]][[j]]], legc[g,mcs_in[[i]][[j]]], col = clrs[[i]], pch = 16, cex = 1)
            points(mcmd$mean_day[mcs_in[[i]][[j]]], legc[g,mcs_in[[i]][[j]]], col = mcmd$color[mcs_in[[i]][[j]]], pch = 16, cex = 1)
        })
        clrmp_points <- colorRampPalette(c('blue', 'red'))(100)
        col_brks <- round(seq(-17, -8, l=100), 2)
        points_y <- log2(1e-5+ctt_traj[[i]][g,])
        points_y_cut <- as.numeric(cut(points_y, col_brks))
        lines(13:18, points_y, col = 'black', lwd = 6, lty = 2)
        points(13:18, points_y, col = clrmp_points[points_y_cut],pch = 16, cex = 5)
        # points(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = 'red', cex = 4, pch = 16)
        # lines(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = 'black', lwd = 6, lty = 2)
    })
    grid(lwd = 2)
#     legend('topright', legend = paste0("t_final = E", 16:18), col = clrs, pch = rep(1,3), lwd = rep(1,3))
#     dev.off()
})
dev.off()

# top_mc_fp_gene = apply(mc@mc_fp, 1, which.max)
# top_fp_gene = apply(mc@mc_fp, 1, max)
# top_mc_fp_gene = apply(legc, 1, which.max)
# top_fp_gene = apply(legc, 1, max)
# max_genes_by_st = lapply(cust_st_ord, function(st) setNames(top_fp_gene[top_mc_fp_gene %in% mcmd$metacell[mcmd$cell_type == st]],
#                                                            names(top_fp_gene)[top_mc_fp_gene %in% mcmd$metacell[mcmd$cell_type == st]]))
fp_by_st = t(tgs_matrix_tapply(mc@mc_fp, mcmd$cell_type, mean))

top_fp_gene = apply(fp_by_st, 1, max)

top_st_fp_gene = colnames(fp_by_st)[apply(fp_by_st, 1, which.max)]

head(top_fp_gene)

head(top_st_fp_gene)

max_genes_by_st = lapply(cust_st_ord, function(st) setNames(top_fp_gene[top_st_fp_gene == st],
                                                           names(top_fp_gene)[top_st_fp_gene == st])
                         )
names(max_genes_by_st) = cust_st_ord

n = 4
top_genes_by_st = lapply(max_genes_by_st, function(x) head(names(x)[order(x, decreasing = T)], n))
top_genes_by_st
top_genes_by_st_vec = unlist(top_genes_by_st)
top_genes_by_st_vec = top_genes_by_st_vec[!(top_genes_by_st_vec %in% c('Nts', 'Gm26917'))]
# top_genes_by_st_vec = top_genes_by_st_vec[top_genes_by_st_vec %in% marks]

# pltmt = t(tgs_matrix_tapply(mc@mc_fp[apply(mc@mc_fp, 1, max) >= 4,], mcmd$cell_type, mean))
pltmt = t(tgs_matrix_tapply(mc@mc_fp[top_genes_by_st_vec,], mcmd$cell_type, mean))

com_pltmt = apply(pltmt[,cust_st_ord[cust_st_ord %in% colnames(pltmt)]], 1, function(x) sum(x*1:length(x))/sum(x))

pltmt = log2(pltmt)
                  
# rownames(pltmt)[order(com_pltmt)]

# options(repr.plot.height = 20, repr.plot.width = 10)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp3 = pheatmap::pheatmap(pltmt[order(com_pltmt),cust_st_ord[cust_st_ord %in% colnames(pltmt)]], breaks = brks,
                   cluster_rows=F,cluster_cols=F, color = colorRampPalette(pltt)(l), fontsize_col = 24,
#                    annotation_col = col_annot, 
#                    annotation_colors = ann_colors, 
                   fontsize_row = 18)

save_pheatmap_png(pp3, './figs/pl_cort_top_genes_by_st.png')

n = 4
top_genes_by_st = lapply(max_genes_by_st, function(x) head(names(x)[order(x, decreasing = T)], n))
# top_genes_by_st
top_genes_by_st_vec = unlist(top_genes_by_st)
top_genes_by_st_vec = top_genes_by_st_vec[!(top_genes_by_st_vec %in% c('Nts', 'Gm26917'))]
# top_genes_by_st_vec = top_genes_by_st_vec[top_genes_by_st_vec %in% marks]

# pltmt = t(tgs_matrix_tapply(mc@mc_fp[apply(mc@mc_fp, 1, max) >= 4,], mcmd$cell_type, mean))
pltmt = mc@mc_fp[top_genes_by_st_vec,]

com_pltmt = apply(pltmt[,cust_st_ord[cust_st_ord %in% colnames(pltmt)]], 1, function(x) sum(x*1:length(x))/sum(x))

pltmt = log2(pltmt)
                  
# rownames(pltmt)[order(com_pltmt)]

# options(repr.plot.height = 20, repr.plot.width = 10)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))

pp4 = pheatmap::pheatmap(pltmt[order(com_pltmt),cust_mc_ord_st], breaks = brks,show_colnames = F, annotation_legend = F,
                   cluster_rows=F,cluster_cols=F, color = colorRampPalette(pltt)(l), fontsize_col = 24,
                   annotation_col = col_annot, 
                   annotation_colors = ann_colors, 
                   fontsize_row = 10)

save_pheatmap_png(pp4, './output/metacell_model/figs/pl_cort_phm_top_genes_by_st.png', h = 1200, w = 2400, res = 200)

mat <- scdb_mat('pl_cort')
mc2d <- scdb_mc2d('pl_cort')

colnames(mc@e_gc) <- 1:ncol(mc@e_gc)

colnames(mc@mc_fp) <- 1:ncol(mc@mc_fp)

scdb_add_mc(id = 'pl_cort', mc)

# mcell_mc2d_plot_by_factor('pl_cort', 'pl_cort', meta_field = 'day', single_plot = T)

# mcell_mc2d_plot_by_factor('pl_cort', 'pl_cort', meta_field = 'day', single_plot = F)

save_pheatmap_png(pp, './figs/pl_cort_goi_hm_for_rp_new.png', h = 2200, w = 3800,r=300)

st_fp <- as.data.frame(t(tgs_matrix_tapply(mc@mc_fp, mcmd$cell_type, mean)))
st_legc <- as.data.frame(t(tgs_matrix_tapply(legc, mcmd$cell_type, mean)))
# cs <- Matrix::colSums(mat@mat)

# summary(mc_umis)
# summary(mcn_umis)

# gs <- scdb_gset('pl_filt_lat')

# c('Pcna', 'Mki67', 'Top2a') %in% names(gs@gene_set)
# gs@gene_set[c('Pcna', 'Mki67', 'Top2a')]

# dror_path <- "/net/mraid20/export/tgdata/users/drorba/ambient/out/mmcortex_new/"

# options(repr.plot.width = 20)
# options(repr.plot.height = 8)
# par(mfrow = c(1,2))
# hist(mc_umis, 50)
# hist(mcn_umis, 50)
# par(mfrow = c(1,2))
# hist(as.numeric(table(mc@mc)), 50)
# hist(as.numeric(table(mcn@mc)), 50)
# options(repr.plot.width = 8)

# mcn <- scdb_mc('pl_cort_test')

# mcn_umis <- tapply(cs[names(mcn@mc)], mcn@mc, sum)

# tgconfig::set_param('mcp_heatmap_height', 2500, 'metacell')
# tgconfig::set_param('mcp_heatmap_width', 2500, 'metacell')

# mcell_mc_plot_marks(mc_id = 'pl_cort_test', gset_id = 'pl_cort_test_marks_f', mat_id = 'pl_cort')

# mc_umis <- tapply(cs[names(mc@mc)], mc@mc, sum)

# boxplot(sz ~ ct, data = as.data.frame(list(sz = mc_umis, ct = mcmd$cell_type)), col = color_key$color[order(color_key$cell_type)])

# boxplot(sz ~ ct, data = as.data.frame(list(sz = as.numeric(table(mc@mc)), ct = mcmd$cell_type)), col = color_key$color[order(color_key$cell_type)])

# pairs <- as.data.frame(list(ct1 = c('SCPN', 'SCPN', 'iCfuPN', 'iCPN/CfuPN', 'iCPN/CfuPN', 'iCPN/CfuPN', 'iCPN/CfuPN', 'NSC'), 
#                             ct2 = c('CPN_L5_6', 'CthPN', 'iCPN/CfuPN', 'iCPN_early', 'iCPN_late', 'SCPN', 'CPN_L5_6', 'Astrocytes')))
# pairs

# for (i in 1:nrow(pairs)) {
#     g1m2 <- head(rownames(st_fp[order(st_fp[,pairs$ct1[[i]]] - st_fp[,pairs$ct2[[i]]], decreasing = T),]))
#     g2m1 <- head(rownames(st_fp[order(st_fp[,pairs$ct2[[i]]] - st_fp[,pairs$ct1[[i]]], decreasing = T),]))
#     plot(colSums(legc[g1m2,]), colSums(legc[g2m1,]), col = mcmd$color, pch = 16, cex = 3, xlab = pairs$ct1[[i]], ylab = pairs$ct2[[i]])
#     text(colSums(legc[g1m2,]), -.5+colSums(legc[g2m1,]), labels = mcmd$metacell)
# }

# # scpn_vs_cpnl5_6_genes <- head(rownames(st_fp[order(st_fp$SCPN - st_fp$CPN_L5_6, decreasing = T),]))
# cpnl5_6_vs_scpn_genes <- head(rownames(st_fp[order(st_fp$CPN_L5_6 - st_fp$SCPN, decreasing = T),]))

# scpn_vs_cthpn <- head(rownames(st_fp[order(st_fp$SCPN - st_fp$CthPN, decreasing = T),]))
# cthpn_vs_scpn_genes <- head(rownames(st_fp[order(st_fp$CthPN - st_fp$SCPN, decreasing = T),]))

# plot(colSums(legc[scpn_vs_cthpn,]), colSums(legc[cthpn_vs_scpn_genes,]), col = mcmd$color, pch = 16, cex = 3)
# text(colSums(legc[scpn_vs_cthpn,]), -1+colSums(legc[cthpn_vs_scpn_genes,]), labels = mcmd$metacell)

# plot(colSums(legc[scpn_vs_cpnl5_6_genes,]), colSums(legc[cpnl5_6_vs_scpn_genes,]), col = mcmd$color, pch = 16, cex = 3)
# text(colSums(legc[scpn_vs_cpnl5_6_genes,]), -1+colSums(legc[cpnl5_6_vs_scpn_genes,]), labels = mcmd$metacell)

scfigs_init('figs')

mcell_mc2d_plot_by_factor(mc2d_id = 'pl_cort', single_plot = T, meta_field = 'day', mat_id = 'pl_cort')

st_fp_ord <- st_fp[order(st_fp$CthPN - st_fp$`CPN_L2-3`),]

st_legc_ord <- st_legc[order(st_legc$CthPN - st_legc$`CPN_L2-3`),]

cust_st_ord

genes_bif <- union(head(rownames(st_legc_ord), 10), tail(rownames(st_legc_ord), 10))
genes_bif

mcs_p <- which(mcmd$cell_type %in% c(
    'IPC',
    'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN'))

options(repr.plot.width = 8)

my_mcell_mc_plot_gg <- function(mc_id, g1, g2) {
    plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
}

png('./figs/Syt4_vs_Sox5.png', h = 800, w = 800) 
par(cex.lab = 3, cex.axis = 1, mar = c(6,6,4,1))
mcell_mc_plot_gg('pl_cort', 'Sox5', 'Syt4', cex = 3, use_egc = T, )
dev.off()
png('./figs/Pou3f1_vs_Nfib.png', h = 800, w = 800) 
par(cex.lab = 3, cex.axis = 1, mar = c(6,6,4,1))
mcell_mc_plot_gg('pl_cort', 'Nfib', 'Pou3f1',  cex = 3, use_egc = T, )
dev.off()

mcell_mc_plot_gg('pl_cort', 'Sox5', 'Syt4', use_egc = T)

ind <- 7
# g <- genes_bif[[ind]]
par(mfrow = c(1,2))
g <- 'Bcl11b'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
g <- 'Satb2'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
par(mfrow = c(1,2))
g <- 'Fezf2'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
g <- 'Cux1'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
par(mfrow = c(1,2))
g <- 'Runx1t1'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)
g <- 'Rnd2'
plot(mcmd$mean_day[mcs_p], legc[g,mcs_p], col = mcmd$color[mcs_p], pch = 16, cex = 2, main = g)

cor_md_by_st = apply(head(legc), 1, function(x) purrr::map2(.x = x, .y = mcmd$mean_day, .f = tapply, ... = mcmd$cell_type))

calc_cor_gene_legc_w_md_in_st <- function(legc, cell_type, mcmd) {
    mcs_st <- which(mcmd$cell_type == cell_type)
    cor_genes_md <- tgs_cor(t(legc[,mcs_st]), as.matrix(mcmd$mean_day[mcs_st]), spearman=T)
    names(cor_genes_md) <- rownames(legc)
#     vec <- sort(cor_genes_md, decreasing = T)
#     return(c(head(vec, 10), tail(vec, 10)))
    return(cor_genes_md)
}

ct_mcs <- lapply(cust_st_ord, function(ct) which(mcmd$cell_type == ct))
names(ct_mcs) <- cust_st_ord

library(ppcor)

my_pcor <- function(x, y, z, spearman = T) {
#     print(dim(x))
#     print(dim(y))
#     print(length(z))
    if (ncol(x) != ncol(y) ||  ncol(x) != nrow(z)) {
        stop('Check that all input data have compatible dimensions (x and y should by N*M matrices and z a M*1 matrix)')
    }
    cor_xz <- tgs_cor(t(x), z, spearman = spearman)
    cor_yz <- tgs_cor(t(y), z, spearman = spearman)
    cor_xy <- tgs_cor(t(x), t(y), spearman = spearman)
    print(dim(cor_xz))
    print(dim(cor_yz))
    res <- (cor_xy - cor_xz%*%t(cor_yz))/(sqrt(1-cor_xz**2)%*%t(sqrt(1-cor_yz**2)))
    return(res)
}

N = 100
K = 1
Z <- rnorm(N)
X <- Z + matrix(rnorm(N*K, sd = 1), nrow = K, ncol = N)
Y <- Z + matrix(rnorm(N*K, sd = 1), nrow = K, ncol = N)
ptt <- my_pcor(X,Y,as.matrix(Z))
ptt
cor.test(X,Y)
pcor.test(X,Y,Z, method = 'spearman')

ctall <- cor.test(X,Y)

ctall$conf.int

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(legc))
# tfs

calc_cor_gene_legc_w_md_in_st <- function(legc, cell_type, mcmd) {
    mcs_st <- which(mcmd$cell_type == cell_type)
    cor_genes_md <- tgs_cor(t(legc[,mcs_st]), as.matrix(mcmd$mean_day[mcs_st]), spearman=T)
    names(cor_genes_md) <- rownames(legc)
#     vec <- sort(cor_genes_md, decreasing = T)
#     return(c(head(vec, 10), tail(vec, 10)))
    return(cor_genes_md)
}

cor_gene_legc_w_md_in_st <- do.call('cbind', lapply(unique(mcmd$cell_type), function(ct) calc_cor_gene_legc_w_md_in_st(legc, ct, mcmd)))

colnames(cor_gene_legc_w_md_in_st) <- unique(mcmd$cell_type)

## TF set for understanding TF modules
tfs_hi <- tfs[matrixStats::rowMaxs(legc[tfs,]) >= -13 & 
              matrixStats::rowMaxs(legc[tfs,]) - matrixStats::rowMins(legc[tfs,]) >= 3 & 
             matrixStats::rowSds(legc[tfs,]) >= 0.5]
tfs_pos_cor <- lapply(cust_st_ord, function(ct) which(cor_gene_legc_w_md_in_st[tfs_hi,ct] >= 0.5))
tfs_neg_cor <- lapply(cust_st_ord, function(ct) which(cor_gene_legc_w_md_in_st[tfs_hi,ct] >= 0.5))

names(tfs_pos_cor) <- cust_st_ord
names(tfs_neg_cor) <- cust_st_ord

## TF set which is more like "markers"
tfs_hi <- tfs[which(log2(matrixStats::rowMaxs(mc@mc_fp[tfs,])) - log2(matrixStats::rowMins(mc@mc_fp[tfs,])) >= 3)]
length(tfs_hi)

cti <- 'Astrocytes'

# pcor_legc_md <- my_pcor(legc[tfs_hi,ct_mcs[[cti]]], legc[tfs_hi,ct_mcs[[cti]]], as.matrix(mcmd$mean_day[ct_mcs[[cti]]]))
pcor_legc_md <- my_pcor(legc[tfs_hi,], legc[tfs_hi,], as.matrix(mcmd$mean_day))
diag(pcor_legc_md) <- 0

# cor_g_md <- tgs_cor(t(legc[tfs_hi,ct_mcs[[cti]]]), as.matrix(mcmd$mean_day[ct_mcs[[cti]]]), spearman = T)
cor_g_md <- tgs_cor(t(legc[tfs_hi,]), as.matrix(mcmd$mean_day), spearman = T)

length(cor_g_md)

library(pheatmap)

options(repr.plot.width = 22)
options(repr.plot.height = 22)

# options(repr.plot.width = 30)
# options(repr.plot.height = 30)

# craw <- tgs_cor(t(legc[tfs_hi,ct_mcs[[cti]]]))
craw <- tgs_cor(t(legc[tfs_hi,]))
diag(craw) <- 0
dim(craw)

hc_raw <- hclust(tgs_dist(craw))
hc_pcor <- hclust(tgs_dist(pcor_legc_md))

acp <- list(cluster = setNames(c('blue', 'yellow', 'red', 'green'), 1:3), 
            cor_g_md = setNames(c('blue2', 'white', 'red2'), c(-1,0,1)),
           mean_day = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100), seq(13,18,l=100)))

ctp <- cutree(hc_pcor, k = 4)

rap <- tibble::column_to_rownames(tibble::enframe(ctp, name = 'gene', value = 'cluster'), 'gene')
rap$cor_g_md <- cor_g_md

phm <- pheatmap(legc[tfs_hi,ct_mcs[['Astrocytes']][order(mcmd$mean_day[ct_mcs[['Astrocytes']]])]], 
                cluster_cols = T,
                annotation_col = tibble::column_to_rownames(mcmd[,c('metacell', 'mean_day')], 'metacell'),
                annotation_row = rap, annotation_colors = acp)
p_raw <- pheatmap(craw[hc_raw$order,hc_raw$order], 
                  cluster_cols = F, cluster_rows = F, 
                  annotation_col = rap, annotation_row = rap, annotation_colors = acp, 
                  fontsize_row = 8, color = colorRampPalette(c('blue4', 'lightblue', 'white', 'yellow3', 'red3'))(100), 
                  breaks = seq(-1,1,l=101))
p_pcor <- pheatmap(pcor_legc_md[hc_pcor$order,hc_pcor$order],
                   cluster_cols = F, cluster_rows = F, 
                   annotation_col = rap, annotation_row = rap, annotation_colors = acp, 
                   fontsize_row = 8, color = colorRampPalette(c('blue4', 'lightblue','white', 'yellow3','red3'))(100), 
                   breaks = seq(-1,1,l=101))

sort(names(ctp)[ctp == 1])

options(repr.plot.width = 24)
options(repr.plot.height = 24)

cl_i <- 1
pheatmap(pcor_legc_md[names(ctp)[ctp == cl_i],names(ctp)[ctp == cl_i]])

length(tfs_hi)

'Zfhx4' %in% tfs_hi

# tfs_pos_cor <- lapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = T), 20))
tfs_pos_cor <- lapply(cust_st_ord, function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = T), 20))
# names(tfs_pos_cor) <- c('NSC', 'IPC', 'IPC_cyc')
# tfs_neg_cor <- lapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = F), 20))
tfs_neg_cor <- lapply(cust_st_ord, function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = F), 20))
# names(tfs_neg_cor) <- c('NSC', 'IPC', 'IPC_cyc')

library(pheatmap)

# mcs_ord <- unlist(sapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) {
mcs_ord <- unlist(sapply(cust_st_ord, function(ct) {
    f <- mcmd$cell_type == ct; 
    return(which(f)[order(mcmd$mean_day[f])])}))

# pltmt <- legc[union(names(unlist(tfs_pos_cor)), names(unlist(tfs_neg_cor))),mcs_ord]
pltmt <- legc[tfs_hi,mcs_ord]
ord <- order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x)))

hcp <- hclust(tgs_dist(pltmt), method = 'ward.D2')

cthcp <- cutree(hcp, 10)

cthcp[hcp$order]

?base::subset

ppp1 <- pheatmap(subset(pltmt, subset = !(rownames(pltmt) %in% c('Sox4', 'Sox11'))), annotation_legend = F,
                color = colorRampPalette(c('white', 'gold', 'red2','black'))(100),
         cluster_cols = F, cluster_rows = T, show_colnames =F, fontsize_row =8,
         annotation_col = col_annot, 
         annotation_colors = ann_colors)

save_pheatmap_png(ppp1, './figs/tf_pheatmap.png', w = 2000, h = 3000, res = 300)

ppp2 <- pheatmap(legc[union(names(unlist(tfs_pos_cor)), names(unlist(tfs_neg_cor))),mcs_ord], 
  
                color = colorRampPalette(c('white', 'gold', 'red2','black'))(100),
         cluster_cols = T, cluster_rows = T,
         annotation_col = col_annot, 
         annotation_colors = ann_colors)

ppp <- pheatmap(legc[unique(c(names(unlist(tfs_pos_cor)), names(unlist(tfs_neg_cor)))),mcs_ord], 
                color = colorRampPalette(c('white', 'gold', 'red2','black'))(100),
         cluster_cols = F, cluster_rows = T,
         annotation_col = col_annot, 
         annotation_colors = ann_colors)

pltmt = log2(mc@mc_fp[tfs_hi,mcs_ord])
pltmt = subset(pltmt, subset = !(rownames(pltmt) %in% c('Sox4', 'Sox11')))
minv = min(pltmt)
maxv = max(pltmt)
rng = maxv-minv
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))+1),
         seq(0+rng/l,maxv/2,l=round(l/(length(pltt)-1))),
        seq(maxv/2+rng/l,maxv,l=round(l/(length(pltt)-1))))

dim(pltmt)

p3 = pheatmap::pheatmap(pltmt[,], show_colnames = F, fontsize_row = 12, treeheight_row = 20, breaks = brks, annotation_legend = F,
                               annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, cluster_rows = T,
#                                color = topo.colors(pal_len, rev = T)
#                                color = colorRampPalette(rev(c("white", "yellow", "red", "black")))(pal_len),
                               color = colorRampPalette(pltt)(l),
#                                breaks = c(seq(min(fp_mc_norm), quantile(fp_mc_norm, 0.75), l=75*pal_len/100),
#                                           seq(quantile(fp_mc_norm, 0.76), quantile(fp_mc_norm, 0.9),l=14*pal_len/100),
#                                           seq(quantile(fp_mc_norm, 0.91), quantile(fp_mc_norm, 0.99), l=8*pal_len/100),
#                                           seq(quantile(fp_mc_norm, 0.991), max(fp_mc_norm), l=pal_len/100))
                               )
save_pheatmap_png(p, './figs/pl_cort_tf_mc_hm.png', width = 2600, height = 2600, res = 300)

fp_mc_tf <- mc@mc_fp[tfs_hi,]

cor_mc_tf = tgs_cor(t(pltmt), spearman = T)

head(cor_mc_tf)

km_tf = tglkmeans::TGL_kmeans(cor_mc_tf, k=12, seed = SEED)

tf_cl = setNames(km_tf$cluster, rownames(cor_mc_tf))
tf_cl[order(tf_cl)]

sort(names(tf_cl))

gaps_vec = tail(unlist(sapply(
    tapply(setNames(1:length(tf_cl), tf_cl[order(tf_cl)]), tf_cl[order(tf_cl)], identity)
    , function(x) head(x, 1))), -1) - 1
gaps_vec

p3$tree_row$order

options(repr.plot.width = 15, repr.plot.height = 15)
pal_len = 50
# p = pheatmap::pheatmap(cor_mc_tf[order(km_tf$cluster, decreasing = T),order(km_tf$cluster, decreasing = T)], 
p = pheatmap::pheatmap(cor_mc_tf[p3$tree_row$order,p3$tree_row$order], 
#                        gaps_col = gaps_vec, gaps_row = gaps_vec,
                       fontsize_col = 14, fontsize_row = 14, cluster_rows = F, cluster_cols = F,
                       silent = F, treeheight_row = 20, treeheight_col = 20, 
                       color = colorRampPalette(c("blue","skyblue","white", "yellow", "red"))(pal_len),
                       breaks = c(seq(-1, 0, length.out = pal_len/2), seq(1/(pal_len/2), 1, l=pal_len/2))
                      )
save_pheatmap_png(p, './figs/pl_cort_tf_cor_clust.png', w = 4200, h = 4200, r = 300)

ls()

quantile(matrixStats::rowMaxs(legc[tfs,]))

save_pheatmap_png(ppp, './figs/hi_var_tfs_nsc_ipc.png')

cor_legc_md <- tgs_cor(t(legc), as.matrix(mcmd$mean_day), spearman = T)

abs_diff_cor_ct_cor_all <- abs(rowMeans(subset(cor_gene_legc_w_md_in_st, select = -c(Astrocytes, OPCs))) - cor_legc_md)
abs_diff_cor_ct_cor_all <- setNames(abs_diff_cor_ct_cor_all[,1], rownames(abs_diff_cor_ct_cor_all))

names(cor_gene_legc_w_md_in_st) <- unique(mcmd$cell_type)

lapply(max_genes_by_st, head)

int_st = c('IPC_cyc','IPC','IPC_late', 'iCPN/CfuPN','iCfuPN','iCPN_early', 'iCPN_late')

max_genes_by_st_int = max_genes_by_st[int_st]

lat = scdb_gset('pl_lateral')

n = 8
top_genes_by_st_int = lapply(max_genes_by_st_int, function(x) head(names(x)[order(x, decreasing = T)], n))
top_genes_by_st_int

top_genes_by_st_int_vec = unlist(top_genes_by_st_int)
top_genes_by_st_int_vec = top_genes_by_st_int_vec[!(top_genes_by_st_int_vec %in% 
                                                    union(names(lat@gene_set), c('Nts', 'Gm26917')))]
# top_genes_by_st_vec = top_genes_by_st_vec[top_genes_by_st_vec %in% marks]

# pltmt = t(tgs_matrix_tapply(mc@mc_fp[apply(mc@mc_fp, 1, max) >= 4,], mcmd$cell_type, mean))
pltmt = t(tgs_matrix_tapply(mc@mc_fp[top_genes_by_st_int_vec,mcmd$metacell[mcmd$cell_type %in% int_st]], 
                            mcmd$cell_type[mcmd$cell_type %in% int_st], mean))
pltmt = pltmt[,int_st]
com_pltmt = apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x))

pltmt = log2(pltmt)
                  
# rownames(pltmt)[order(com_pltmt)]

# options(repr.plot.height = 20, repr.plot.width = 10)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp3 = pheatmap::pheatmap(pltmt[order(com_pltmt),], breaks = brks,
                   cluster_rows=F,cluster_cols=F, color = colorRampPalette(pltt)(l), fontsize_col = 18,
#                    annotation_col = col_annot, 
#                    annotation_colors = ann_colors, 
                   fontsize_row = 10)

save_pheatmap_png(pp3, './figs/pl_cort_int_st_hm.png', h = 1800, w = 800)

library(Matrix)
scdb_init("scdb/", force_reinit=T)
nm = 'pl_cort'

get_mc_cc = function() {
  ## function from Markus
  mat_id = nm
  mc_id = nm

  tag = nm
  m_0 = 0.0025
  s_0 = 0.001
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))
  mc2d = scdb_mc2d('pl_cort')

  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  png("figs/pl_cort_mc2d_cc.png", w=800, h=800)
  plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.4, col=ifelse(f[names(mc2d@sc_x)], "lightgray","black"))
  points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=2.5, bg=shades[101 - mc_cc$cc_score])
  dev.off()
  return(mc_cc)
}

mc_cc = get_mc_cc()

png('./figs/Bcll11b_vs_Satb2.png', h = 800, w = 800) 
par(cex.lab = 3, cex.axis = 1, mar = c(6,6,4,1))
mcell_mc_plot_gg('pl_cort', 'Bcl11b', 'Satb2', cex = 3, use_egc = T, )
dev.off()

png('./output/metacell_model/figs/cell_cycle_score_per_ct_vioplot.png', h = 350, w = 700)
par(las = 2, mar = c(12,7,1,1), cex.axis= 1.5, cex.lab = 2)
vioplot(100 - mc_cc$cc_score ~ factor(mcmd$cell_type, levels = cust_st_ord), col = col_key[cust_st_ord], ylab = '', xlab = '')
title(xlab= 'Cell type', line = 9)
title(ylab= 'Cell cycle score', line = 4)
dev.off()

options(repr.plot.width = 10)
options(repr.plot.height = 5)

cust_st_ord

col_key

mc_cc$cc_score

png('./output/metacell_model/figs/cell_cycle_score_per_ct_barplot.png', h = 350, w = 700)
par(las = 2, mar = c(12,7,1,1), cex.axis= 1.5, cex.lab = 2)
boxplot(100 - mc_cc$cc_score ~ factor(mcmd$cell_type, levels = cust_st_ord), col = col_key[cust_st_ord], ylab = '', xlab = '')
title(xlab= 'Cell type', line = 9)
title(ylab= 'Cell cycle score', line = 4)
dev.off()

plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6))
{
  if (!is.null(fig_fn)) {
#     .plot_start(fig_fn, 400, 400)
      png(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

#   if (is.null(show_vals_ind)) {
#     show_vals_ind = rep(T, length(cols))
#   }
  text(19, show_vals_ind,cex = 2, labels=round(vals[show_vals_ind], 3), pos=4)
#   text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
}
shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
min_val = min(mc_cc$cc_score)
max_val = max(mc_cc$cc_score)
plot_color_bar(seq(min_val-1, max_val,l=101), shades, fig_fn = './figs/pl_cort_cc_colorbar.png')

mc2d <- scdb_mc2d('pl_cort')

mat <- scdb_mat('pl_cort')

clrs <- c('red3', 'orange3', 'yellow3', 'green4', 'blue3', 'purple2')

days <- tail(unique(mat@cell_metadata$day), -1)

days

mc2d <- scdb_mc2d('pl_cort_not_cor_cc')

scdb_ls('mc2d')

length(mc2d@)

options(repr.plot.width = 16)
options(repr.plot.height = 16)

png('./output/metacell_model/figs/single_cell_day_of_origin.png', h = 2500, w = 2500, res = 150)
par(cex.main = 6, cex.lab = 2, mar = c(5,5,4,3))
plot(mc2d@sc_x, mc2d@sc_y, col = 'white', main = 'Cells by day of origin', xaxt = 'n', yaxt  = 'n', bty = 'n', xlab = '', ylab = '')
sapply(seq_along(days), function(i) {
    cells_d <- rownames(mat@cell_metadata)[mat@cell_metadata$day == days[[i]] & rownames(mat@cell_metadata) %in% names(mc2d@sc_x)]
    print(head(cells_d))
    points(mc2d@sc_x[cells_d], mc2d@sc_y[cells_d], col = clrs[[i]], pch = 16, cex = 1)
    })
legend('topright', legend = days, col = clrs, pch = rep(16, length(clrs)), cex = 5, bty = 'n')
dev.off()

load('./data/gene_modules_mcmd_pl_cort.Rda')

gene_modules

# ar <- as.data.frame(tibble::enframe(do.call('c', lapply(cust_st_ord, function(x) setNames(gene_modules[[x]], rep(x, 10)))), name = 'cell_type', value = 'gene'))
# rownames(ar) <- ar$gene
# ar <- ar[,c('cell_type')]
# rownames(ar) <- 1:nrow(ar)
# colnames(ar) <- 'cell_type'

all(

cumsum(table(mcmd$cell_type)[cust_st_ord])

gene_modules

pltmt <- legc[unlist(gene_modules[cust_st_ord]),cust_mc_ord_st]
pltmt <- pltmt - rowMeans(pltmt)
p_gene_module_phm <- pheatmap::pheatmap(pltmt, cluster_rows = F, cluster_cols = F, gaps_row = seq(10,10*length(gene_modules), 10), gaps_col = cumsum(table(mcmd$cell_type)[cust_st_ord]),
                   # annotation_row = ar, 
                   show_colnames = F,
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=101), annotation_col = col_annot, annotation_colors = ann_colors, annotation_legend = F)

save_pheatmap_png(p_gene_module_phm, './output/metacell_model/figs/gene_module_heatmap.png', h = 2800, w = 4200, res = 200)

            
st_by_day_lst = sapply(unique(mcmd$cell_type), function(st) {table(mat@cell_metadata$t[
    rownames(mat@cell_metadata) %in% names(mc@mc)[mc@mc %in% mcmd$metacell[mcmd$cell_type == st]]])})
st_by_day_mat = matrix(0, ncol = length(13:18), nrow = length(unique(mcmd$cell_type)))
colnames(st_by_day_mat) = 13:18
rownames(st_by_day_mat) = names(st_by_day_lst)

# st_by_day_mat

for (i in 1:length(st_by_day_lst)) {
    st_by_day_mat[names(st_by_day_lst)[[i]],names(st_by_day_lst[[i]])] = st_by_day_lst[[i]]
}

st_by_day_mat

st_by_day_mat_norm = apply(st_by_day_mat, 2, function(x) round(x/sum(x), 3))
st_by_day_mat_norm

st_by_day_mat_norm_st = t(apply(st_by_day_mat, 1, function(x) round(x/sum(x), 3)))
st_by_day_mat_norm_st

colSums(st_by_day_mat)

# lmo = lm(log2(st_by_day_mat_norm['NSC',]) ~ as.numeric(colnames(st_by_day_mat_norm)))
# lmo

# lmo$coefficients[[2]]

# y = st_by_day_mat_norm['NSC',]
# x = as.numeric(colnames(st_by_day_mat_norm))
# plot(x,log2(y))
# lines(x,lmo$coefficients[[1]] + lmo$coefficients[[2]]*x)

pp1 = pheatmap::pheatmap(st_by_day_mat_norm[cust_st_ord[cust_st_ord %in% rownames(st_by_day_mat_norm)],], 
                        fontsize = 24,
                        color = colorRampPalette(c('white', 'yellow', 'red', 'black'))(100), cluster_cols = F, cluster_rows = F)
pp2 = pheatmap::pheatmap(st_by_day_mat_norm_st[cust_st_ord[cust_st_ord %in% rownames(st_by_day_mat_norm_st)],], 
                        fontsize = 24,
                        color = colorRampPalette(c('white', 'yellow', 'red', 'black'))(100), cluster_cols = F, cluster_rows = F)

save_pheatmap_png(pp1, './figs/pl_cort_st_by_day_norm_day.png', h=800,w=650,res=100)
save_pheatmap_png(pp2, './figs/pl_cort_st_by_day_norm_st.png', h=800,w=650,res=100)

plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6))
{
  if (!is.null(fig_fn)) {
#     .plot_start(fig_fn, 400, 400)
      png(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

#   if (is.null(show_vals_ind)) {
#     show_vals_ind = rep(T, length(cols))
#   }
  text(19, show_vals_ind,cex = 2, labels=round(vals[show_vals_ind], 3), pos=4)
#   text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
}
clrmp = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'violet'))(100)

plot_color_bar(seq(13,18,l=101), clrmp, fig_fn = './figs/pl_cort_mean_day_colorbar.png')

flow_res_path <- file.path(wd, "output/mcatac/pl_cort_flow_mat.tsv")
flow_mat <- as.matrix(tibble::column_to_rownames(readr::read_tsv(flow_res_path), 'rowname'))

mc_from_mcl_flow = flow_res$mc_from_mcl_flow

min_cov_g = rowSums(mc_from_mcl_flow) > 1
cs = colSums(mc_from_mcl_flow)
gene_folds = t(t(mc_from_mcl_flow[min_cov_g,])/cs)

fp_reg = 1e-05
flow_fp = (gene_folds + fp_reg)/(apply(gene_folds, 1, median) + fp_reg)



flow_com = apply(flow_mat[,cust_mc_ord_st], 1, function(x) sum(x*(1:length(x)))/sum(x))

# flow_com

p = pheatmap::pheatmap(flow_mat[order(flow_com),cust_mc_ord_st], 
#                        border_color = 'black',
                       show_colnames = F, 
                       show_rownames = T, 
                       cluster_cols = F, 
                       cluster_rows = F,
                       gaps_col = tail(sapply(cust_st_ord, function(u) min(which(names(cust_mc_ord_st) == u))), -1) - 1,
                    color = colorRampPalette(c('black', 'red', 'orange','yellow', 'white'))(100),
                   annotation_col = col_annot, 
                   annotation_colors = ann_colors)
save_pheatmap_png(p, './figs/pl_cort_flow_mat.png')



# options(repr.plot.height = 18, repr.plot.width = 14)
par(mfrow = c(1,2))
plot_mat = log2(mc@mc_fp[marks_filt,cust_mc_ord_st])
min_val = min(plot_mat)
max_val = max(plot_mat)
plot_com_ord = order(apply(plot_mat, 1, function(x) sum(x*1:length(x))/sum(x)))
p = pheatmap::pheatmap(plot_mat[plot_com_ord,], cluster_rows = F, cluster_cols = F, 
                   annotation_col = col_annot, show_rownames = T, show_colnames = F,
                  annotation_colors = ann_colors, fontsize_row = 16,
                       breaks = c(seq(min_val, 0,l=50), seq(0.01,max_val,l=51)),
                      color = colorRampPalette(c('blue','white', 'red4'))(100),
#                        color = colorRampPalette(c('white', 'yellow', 'blue', 'brown', 'black'))(1000) 
                      )
# options(repr.plot.height = 6, repr.plot.width = 6)
save_pheatmap_png(p, './figs/pl_cort_rna_mc_hm.png', h = 1600, w = 4000)

# genes = 
# genes_ord_by_rna = p$tree_row$labels[p$tree_row$order]

# options(repr.plot.height = 18, repr.plot.width = 14)
pltmt = flow_fp[marks_filt,cust_mc_ord_st]
# min_val = min(plot_mat)
# max_val = max(plot_mat)

# plot_com_ord = order(apply(pltmt, 1, function(x) sum(x*1:length(x))/sum(x)))
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
p2 = pheatmap::pheatmap(pltmt[plot_com_ord,], cluster_rows = F, cluster_cols = F, 
                    annotation_col = col_annot, show_rownames = T, show_colnames = F,
                    annotation_colors = ann_colors, fontsize_row = 16,
#                     breaks = c(seq(min_val, 0,l=50), seq(0.01,max_val,l=51)),
                       breaks = brks,
                       color = colorRampPalette(pltt)(l),
#                     color = colorRampPalette(c('blue','white', 'red'))(100),
#                       color = colorRampPalette(c('white', 'red', 'black'))(100),
                      )
save_pheatmap_png(p2, './figs/pl_cort_atac_mc_hm.png', h = 1600, w = 4000)
options(repr.plot.height = 6, repr.plot.width = 6)

gstat = scdb_gstat('pl_prom_cort')
x = log(gstat$ds_mean)
init_filt = which(x >= -5)
y = gstat$ds_log_varmean[init_filt]
name_filt = gstat$name
x = x[init_filt]

xcut = cut(x, breaks = seq(min(x), max(x), l = 40))
top_q_inds = lapply(levels(xcut), function(l) {inds = which(xcut == l); 
                                  xfilt = y[inds]; 
                                  xtop = head(inds[order(xfilt, decreasing = T)], 40); 
                                  return(xtop)
                                 }
      )

names(top_q_inds) = levels(xcut)
feats = intersect(rownames(gstat)[init_filt[unlist(top_q_inds)]], rownames(mc@mc_fp))
ifn1_genes = c('Isg15', 'Wars', 'Ifit1')
cell_cyc = c('Mki67', 'Pcna', 'Smc4', 'Mcm3', 'Top2a')
stress = c('Fos', 'Hsp90ab1', 'Hspa1a', 'Hif1a')
misc = c('Xist', 'Tsix')
star_genes = c(ifn1_genes, cell_cyc, stress, misc)
star_genes = c(star_genes, c(grep('Mmc', feats, v=T),
                            grep('^Smc\\d', feats, v=T),
                            grep('^Cdk', feats, v=T),
                            grep('^Ccn', feats, v=T),
                             grep('^Ube', feats, v=T),
                             grep('^Rpl', feats, v=T),
                              grep('^Rps', feats, v=T)
                            )
              )

feats_filt = unique(feats[!(feats %in% star_genes)])


# mc_ord = mcs[[i]][order(mcmd$cell_type[mcmd$metacell %in% mcs[[i]]])];
# mc_ord = cust_mc_ord_st[cust_mc_ord_st %in% mc_ord]
# new_df = data.frame(setNames(annotation_col[mc_ord,], mc_ord))
# colnames(new_df) = 'CellType'
cor_mat = tgs_cor(mc_from_mcl_flow[feats_filt,], log2(mc@e_gc[feats_filt,] + 1e-7), spearman = T)
# rownames(cor_mat) = mc_ord
# colnames(cor_mat) = mc_ord
p = pheatmap::pheatmap(cor_mat[cust_mc_ord_st,cust_mc_ord_st], cluster_cols=F, cluster_rows = F, 
                    show_rownames = F,
                    show_colnames = F,
                    annotation_row = col_annot,
                    annotation_col = col_annot, 
                    annotation_colors = ann_colors,
#                         silent = T,
                    color = colorRampPalette(c('blue','green','yellow','red'))(1000),
                  );
save_pheatmap_png(p, glue::glue('./figs/pl_cort_cor_mcatac_mcrna.png'))

# cor_mat_marks = tgs_cor(t(mc_from_mcl_flow[marks_filt,]), t(log2(mc@e_gc[marks_filt,] + 1e-7)), spearman = T)
# sort(diag(cor_mat_marks))
# pheatmap::pheatmap(cor_mat_marks)
# cor_mat_marks

peak_mc = readRDS('./data/pl_cort_peak_mc.rds')

peak_st = t(tgs_matrix_tapply(peak_mc, mcmd$cell_type, mean))

peak_st_norm = t(apply(peak_st, 1, function(x) x/sum(x)))

icpn_st = c('iCPN/CfuPN', 'iCfuPN','iCPN_early', 'iCPN_late')
icpn_peak3 = which(apply(peak_st_norm, 1, max) >= 0.12 & 
                   colnames(peak_st)[apply(peak_st_norm, 1, which.max)] %in% c('CPN_L2-3', 'CPN_L5_6', icpn_st))

st_annot = tibble::column_to_rownames(data.frame(cbind(colnames(peak_st_norm), colnames(peak_st_norm))), 'X1')
colnames(st_annot)=  'st'

options(repr.plot.width = 8, repr.plot.height = 8)

pp = pheatmap::pheatmap(peak_st_norm[icpn_peak3,cust_st_ord], annotation_col = st_annot, show_rownames = F, fontsize = 14,
                   cluster_cols = F, annotation_colors = ann_colors)

save_pheatmap_png(pp, './figs/pl_cort_icpn_enh_hm.png', h = 1600, w=1400)

tfoi = unlist(read.delim('./data/tfoi.txt'))

distt = as.matrix(dist(legc_st))

options(repr.plot.width = 8)
options(repr.plot.height = 8)

pheatmap::pheatmap(distt[cust_st_ord,cust_st_ord], cluster_cols = F, cluster_rows = F)

# heatmap_genes_along_trajectory = function(mc_prob,mc,fig_dir,tag,min_max_fold =0.5,min_thr = -15,main = "",show_rownames = T,
#                                           highlighted_genes = NULL) {
  
  
#   legc = log2(mc@e_gc + 1e-5)
  
#   mc_prob = t(t(mc_prob)/colSums(mc_prob))
#   legc_traj = legc %*% mc_prob
#   min_traj = apply(legc_traj,1,min)
#   max_traj = apply(legc_traj,1,max)
#   cond1 <- max_traj - min_traj > min_max_fold
#   cond2 <- max_traj > min_thr
#   cond3 <- rownames(legc_traj) %in% tfoi
#     cond <- cond1 & cond2 & cond3
#   wc <- which(cond)
#     gnms_f = head(rownames(legc_traj)[wc[order((max_traj - min_traj)[wc], decreasing = T)]], 20)
#   legc_traj_n = (legc_traj[gnms_f,] - min_traj[gnms_f])/(max_traj[gnms_f]-min_traj[gnms_f])
  
#   legc_traj_nn = legc_traj_n/rowSums(legc_traj_n)
  
#   gene_mean_time = (legc_traj_nn %*% c(1:6))[,1]
#   gene_mode = apply(legc_traj_nn,1,which.max)
#   gene_rank = 10*gene_mode + gene_mean_time
    
#   gene_ord = rev(order(gene_rank))
#   print(dim(legc_traj_n))
#   legc_traj_n = legc_traj_n[gene_ord,]
#   legc_traj_n = legc_traj_n[c(nrow(legc_traj_n):1),]
  
#   nm_to_pos = seq(0,1,length.out = nrow(legc_traj_n))
#   names(nm_to_pos) = rownames(legc_traj_n)
#   colnames(legc_traj_n) = c(1:ncol(legc_traj_n))+12
  
#   shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9,name = "PuBu"))(1000)
#     p = pheatmap::pheatmap(legc_traj_n, color = shades, cluster_rows = F, cluster_cols = F, legend = F,
#                            fontsize_row = 14,
# #                            fontsize_row = 18 - nrow(legc_traj_n)/3, 
#                            main = paste(tag, mcmd$cell_type[as.numeric(tag)], collapse = '--'),
# #                            nrow(legc_traj_n)/30,
#                            fontsize_col = 20)
#     st = mcmd$cell_type[as.numeric(tag)]
#     save_pheatmap_png(p, paste0('./figs/ct_traj/', st, '_', tag, '.png'), h=300+50*nrow(legc_traj_n), w = 600)

# }

# plot_heatmap_trajectory = function() {
  
#   mc = scdb_mc(nm)
#   mat = scdb_mat(nm)
#   mct_id = nm
#   mct = scdb_mctnetwork(mct_id)
#   mcf = scdb_mctnetflow(mct_id)
#   marks = scdb_gset('pl_cort_marks_f')
#   marks = names(marks@gene_set)
  
#   mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"t"])
#   mc_ag_n = mc_ag/rowSums(mc_ag)
#   mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
#   mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
#   late_st = c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')
#   st_mcs = lapply(late_st, function(x) mcmd$metacell[mcmd$cell_type %in% x])
#   late_mcs = lapply(st_mcs, function(x) x[mc_ag_cn[x,ncol(mc_ag_cn)] > 0.3])
#   mcs_sel <- as.numeric(sapply(late_mcs, function(m) {
#     p_mc = rep(0,ncol(mc@e_gc))
#     p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
#     probs = mctnetflow_propogate_from_t(mcf, ncol(mc_ag_c), p_mc)
#     mc_prob = t(t(probs$probs)/colSums(probs$probs))
#     csi <- Matrix::colSums(probs$step_m[[length(probs$step_m)]][,m])
#     mci <- names(csi)[which.max(csi)]
#     return(mci)
#   }))
# #   print(st_mcs)
# #     print(late_mcs)
#   df_param = data.frame(mc = mcs_sel,
#                           edge_w_scale = c(5e-5,5e-5,5e-5,5e-5),
#                           fr_scale = c(1,1,1,1),
#                           max_lwd = c(20,20,20,20),
#                           lfp_thr = c(1,1.5,1.5,1.5),
#                           min_thr = c(-15,-14,-14,-14))
#     print(df_param)
  
#   for (i in 1:nrow(df_param)) {
#     m = df_param$mc[i]
#     lfp_thr = df_param$lfp_thr[i]
#     min_thr = df_param$min_thr[i]
    
#     p_mc = rep(0,ncol(mc@e_gc))
#     p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
#     probs = mctnetflow_propogate_from_t(mcf, ncol(mc_ag_c), p_mc)
#     mc_prob = t(t(probs$probs)/colSums(probs$probs))
#     markers = marks[marks %in% rownames(mc@e_gc)]
#     heatmap_genes_along_trajectory(mc_prob = mc_prob,mc = mc,
#                                    fig_dir = fig_dir,
#                                    tag = as.character(m),
#                                    min_max_fold = lfp_thr,
#                                    main = sprintf("Metacell %d",m),
#                                    min_thr = min_thr,show_rownames = F,
#                                    highlighted_genes = markers)
    
    
#   }
  
# }

# if (!dir.exists('./figs/ct_traj')) {dir.create('./figs/ct_traj')}

# plot_heatmap_trajectory()

# # g <- 'Fezf2'
# # cti <- 'CPN_L2-3'
# # cti <- 'SCPN'
# for (i in 1:nrow(gct_df)) {    
#     g <- as.character(gct_df[i,2])
#     cti <- gct_df[i,1]
#     ctt_traj_ls <- lapply(4:6, function(tf) {
#         mcsi <- which(mc_ag_cn[,tf] >= 0.4 & mcmd$cell_type == cti)
#         return(get_egc_trajectory(mc_final = mcsi, mcf = mcf, t_final = tf, t_initial = 1))
#     })
#     ctt_traj <- lapply(ctt_traj_ls, function(x) x$egc_t)
#     mcs_in <- lapply(ctt_traj_ls, function(x) lapply(x$step_m, function(y) which(rowSums(y) > 0)))
#     png(glue::glue('./figs/tf_tracebacks/{cti}_{g}.png'))
#     par(cex.main = 3)
#     ctt_all <- log2(1e-5+do.call('rbind', lapply(ctt_traj, function(x) x[g,])))
# #     print(ctt_all)
#     plot(0, xlim = c(13,18),ylim = c(min(ctt_all),max(ctt_all)), col = 'white', ylab = 'log2 RNA', xlab = 'day', main = paste(cti, g))
#     sapply(seq_along(ctt_traj), function(i) {
#         points(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
#         lines(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
#         sapply(seq_along(mcs_in[[i]]), function(j) {
#             print(mcs_in[[i]][[j]])
#             points(rep(12+j, length(mcs_in[[i]][[j]])), legc[g,mcs_in[[i]][[j]]], col = clrs[[i]])
#         })
#     })
#     legend('topright', legend = paste0("t_final = E", 16:18), col = clrs, pch = rep(1,3), lwd = rep(1,3))
#     dev.off()
# }

# g <- 'Fezf2'
# cti <- 'CPN_L2-3'

# ctt_traj <- lapply(4:6, function(tf) {
#         mcsi <- which(mc_ag_cn[,tf] >= 0.4 & mcmd$cell_type == cti)
#         get_egc_trajectory(mc_final = mcsi, mcf = mcf, t_final = tf, t_initial = 1)
#     })

# par(cex.main = 3)
# plot(0, xlim = c(13,18),ylim = c(-17,-8), col = 'white', ylab = 'log2 RNA', xlab = 'day', main = paste(cti, g))
# sapply(seq_along(ctt_traj), function(i) {
#     points(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
#     lines(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
# })
# legend('topright', legend = paste0("t_final = E", 16:18), col = clrs, pch = rep(1,3), lwd = rep(1,3))


# par(cex.main = 3)
# plot(0, xlim = c(13,18),ylim = c(-17,-8), col = 'white', ylab = 'log2 RNA', xlab = 'day', main = g)
# sapply(seq_along(ctt_traj), function(i) {
#     points(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
#     lines(13:18, log2(1e-5+ctt_traj[[i]][g,]), col = clrs[[i]])
# })
# legend('topright', legend = paste0("t_final = E", 16:18), col = clrs, pch = rep(1,3), lwd = rep(1,3))

# colSums(scpn_3_5_traj[[1]])

# g <- 'Ldb2'
# plot(0, xlim = c(13,18),ylim = c(-17,-10), col = 'white')
# sapply(seq_along(scpn_3_5_traj), function(i) {
#     points(13:18, log2(1e-5+scpn_3_5_traj[[i]][g,]), col = clrs[[i]])
#     lines(13:18, log2(1e-5+scpn_3_5_traj[[i]][g,]), col = clrs[[i]])
# })

# g <- 'Ldb2'
# plot(0, xlim = c(13,18),ylim = c(-17,-10), col = 'white')
# sapply(seq_along(scpn_4_6_traj), function(i) {
#     points(13:18, log2(1e-5+scpn_4_6_traj[[i]][g,]), col = clrs[[i]])
#     lines(13:18, log2(1e-5+scpn_4_6_traj[[i]][g,]), col = clrs[[i]])
# })

# mc_pf <- prop_mc_back_and_forth(mci = mcsi,mcf =  mcf,t_final =  6, t_initial = 1)

# mcs_in <- mcsi
# mc_pf_f <- mc_pf
# miv <- c()
# miv[[6]] <- mcs_in
# for (i in length(mc_pf_f$step_m):1) {
#     print(i)
#     print(mcs_in)
#     mc_pf_f$step_m[[i]][,which(!(1:ncol(mc_pf_f$step_m[[i]]) %in% mcs_in))] <- 0
#     if (length(mcs_in) > 1) {
#         mcs_in <- unique(unlist(apply(mc_pf_f$step_m[[i]][,mcs_in], 2, function(x) which(x > 0))))
#     } else {
#         mcs_in <- unique(which(mc_pf_f$step_m[[i]][,mcs_in] > 0))
#     }
#     miv[[i]] <- mcs_in
#     print(mcs_in)
# }

# nr <- nrow(mc_pf_f$probs)
# for (i in 1:length(miv)) {
#     mc_pf_f$probs[which(!(1:nr %in% miv[[i]])),i] <- 0
# }

# mpfpn <- t(t(mc_pf_f$probs)/colSums(mc_pf_f$probs))

# library(pheatmap)

# options(repr.plot.width = 7)
# options(repr.plot.height = 17)

# apply(mc_pf_f$probs, 2, function(x) which(x > 0))

# pheatmap(mpfpn, cluster_rows = F, cluster_cols = F)

# egc_t = mc@e_gc %*% mpfpn

# options(repr.plot.width = 7)
# options(repr.plot.height = 7)

# plot(log2(egc_t["Ldb2",] + 1e-5), ylim = c(-16,-8))

# colSums(mc_pf_f$probs)

# miv

# mcmd$cell_type[c(487, 494, 211, 578)]

# head(sort(mc_pf$probs[,6]/sum(mc_pf$probs[,6]), decreasing = T), 10)

# head(sort(rowSums(propagated_mc$step_m[[1]]), decreasing = T), 10)

# mcmd$cell_type[as.numeric(names(head(sort(rowSums(propagated_mc$step_m[[1]]), decreasing = T), 10)))]

# rs10 <- head(sort(rowSums(propagated_mc$step_m[[1]]), decreasing = T), 10)
# tapply(rs10, mcmd$cell_type[as.numeric(names(rs10))], sum)

# rs10

# cumsum(rs10)

# mcmd$cell_type[as.numeric(names(head(sort(rowSums(propagated_mc$step_m[[1]]), decreasing = T), 10)))]

# mcmd$cell_type[which(rowSums(propagated_mc$step_m[[1]]) > 0)]

# sort(colSums(propagated_mc$step_m[[2]]))

# options(repr.plot.width = 7)
# options(repr.plot.height = 7)

# plot(log2(egc_t["Sox5",] + 1e-5))

#  fig3_plot_line_graphs(plot_pdf = F)

# mct_extract_subtrajectory = function(probs_trans,types1,types2,t_list,mct,mc) {
  
#   n_mc = ncol(mc@e_gc)
#   # generate mc_forward and mc_backward matrices from probs_trans 
#   # that are used below to propagate p_mc vectors
#   mc_forward = lapply(probs_trans$step_m,function(mat_tr) {
    
#     f = rowSums(mat_tr) > 0
#     if(sum(f) > 1) {
#       mat_tr[f,] = mat_tr[f,]/rowSums(mat_tr[f,])
#     } else {
#       mat_tr[f,] = mat_tr[f,]/sum(mat_tr[f,])
#     }
    
    
#     return(mat_tr)
#   })
  
#   mc_backward = lapply(probs_trans$step_m,function(mat_tr) {
    
#     f = colSums(mat_tr) > 0 
#     if(sum(f) > 1) {
#       mat_tr[,f] = t(t(mat_tr[,f])/colSums(mat_tr[,f]))
#     } else {
#       mat_tr[,f] = mat_tr[,f]/sum(mat_tr[,f])
#     }
#     return(mat_tr)
#   })
#   mcs_type1 = which(mc@colors %in% types1)
#   mcs_type2 = which(mc@colors %in% types2)
#   v0 = rep(0,ncol(mc@e_gc))
#   probs_trans_diff = probs_trans
#   probs_trans_new = matrix(0, nrow = nrow(probs_trans$probs), ncol=ncol(probs_trans$probs))
#   step_m_new = list()
#   for (i in 1:length(probs_trans$step_m)) {
#     step_m_new[[i]] = matrix(0,nrow = nrow(probs_trans$probs),ncol = nrow(probs_trans$probs))
#   }
#   probs_new = matrix(0, nrow = nrow(probs_trans_diff$probs), ncol=ncol(probs_trans_diff$probs))
  
#   for (t1 in t_list) {
#     transition_t = probs_trans$step_m[[t1]]
#     mc_vec2 = v0
#     mc_vec2[mcs_type2] = 1
#     mc_vec1 = v0
#     mc_vec1[mcs_type1] = 1
#     f_row = (transition_t %*% mc_vec2)[,1]  > 0 
#     f_col = (mc_vec1 %*% transition_t)[1,]  > 0
#     mc_list_t1 = intersect(mcs_type1,c(1:n_mc)[f_row])
#     mc_list_t2 = intersect(mcs_type2,c(1:n_mc)[f_col])
#     max_t = ncol(probs_trans_diff$probs)
#     step_m = list()
#     probs = matrix(0, nrow = nrow(probs_trans_diff$probs), ncol=ncol(probs_trans_diff$probs))
#     p_mc_t1 = rep(0,nrow(probs_trans_diff$probs))
#     p_mc_t2 = p_mc_t1
#     mc_vec2 = v0
#     mc_vec2[mc_list_t2] = 1
#     mc_vec1 = v0
#     mc_vec1[mc_list_t1] = 1
#     p_mc_t1[mc_list_t1] = (probs_trans_diff$step_m[[t1]] %*% mc_vec2)[mc_list_t1,1]
#     p_mc_t2[mc_list_t2] = (mc_vec1 %*% probs_trans_diff$step_m[[t1]])[1,mc_list_t2]
#     probs[,t1] = p_mc_t1
#     probs[,t1 + 1] = p_mc_t2
#     probs_trans_diff$probs[,t1] = probs_trans_diff$probs[,t1] - p_mc_t1
#     probs_trans_diff$probs[,t1 + 1] = probs_trans_diff$probs[,t1 + 1] - p_mc_t2
    
#     mat_tr_t1 = Matrix(0,nrow = nrow(probs_trans$probs),ncol = nrow(probs_trans$probs),sparse = T)
#     mat_tr_t1[mc_list_t1,mc_list_t2] = probs_trans_diff$step_m[[t1]][mc_list_t1,mc_list_t2]
#     step_m[[t1]] = mat_tr_t1
#     probs_trans_diff$step_m[[t1]] = probs_trans_diff$step_m[[t1]] - mat_tr_t1
    
#     if(t1 > 1) {
#       for(i in (t1-1):1) {
#         step_m[[i]] = Matrix(t(t(as.matrix(mc_backward[[i]])) * probs[,i+1]), sparse=T)
        
#         probs_trans_diff$step_m[[i]] = probs_trans_diff$step_m[[i]] - step_m[[i]]
        
#         probs[,i] = as.matrix(mc_backward[[i]]) %*% probs[,i+1]
        
#         probs_trans_diff$probs[,i] = probs_trans_diff$probs[,i] - probs[,i]
#       }
#     }
#     if(t1 < max_t - 1) {
#       for(i in (t1+2):max_t) {
#         step_m[[i-1]] = Matrix(as.matrix(mc_forward[[i-1]]) * probs[,i-1], sparse=T)
        
#         probs_trans_diff$step_m[[i-1]] = probs_trans_diff$step_m[[i-1]] - step_m[[i-1]]
        
#         probs[,i] = t(probs[,i-1]) %*% as.matrix(mc_forward[[i-1]])
        
#         probs_trans_diff$probs[,i] = probs_trans_diff$probs[,i] - probs[,i]
#       }
#     }
#     probs_new = probs_new + probs
#     for (i in 1:length(step_m_new)) {
#       step_m_new[[i]] = step_m_new[[i]] + step_m[[i]]
#     }
#   }
#   probs_trans_new = list(probs=probs_new, step_m=step_m_new)
#   return(list(probs_trans_new = probs_trans_new,probs_trans_diff = probs_trans_diff))
# }

# plot_genes_along_trajectory3 = function(mc_prob_ls,mc,fig_dir,genes,show_sd = T,plot_pdf = T,ylim_max = NULL,ylim_min = NULL,add_gridline = F,col = "black",tag = NULL) {
  
#   legc = log2(mc@e_gc + 1e-5)
  
#   out_ls = lapply(mc_prob_ls,function(mc_prob) {
#     mc_prob = t(t(mc_prob)/colSums(mc_prob))
    
    
    
#     feat_traj = comp_e_gc_along_trajectory(mc_prob,legc)
#     print(dim(feat_traj[[1]]))
#       print(dim(feat_traj[[2]]))
#       print(feat_traj[[1]][1:4,1:4])
#       print(feat_traj[[2]][1:4,1:4])
#     e_gc_mean = feat_traj$e_gc_mean
#     e_gc_sd = feat_traj$e_gc_sd
    
    
#     df_mean = as.data.frame(as.table(e_gc_mean))
#     colnames(df_mean) = c("gene","age_group","e_gc")
#     df_sd = as.data.frame(as.table(e_gc_sd))
#     colnames(df_sd) = c("gene","age_group","sd_e_gc")
    
#     if(!identical(df_mean[,1:2],df_sd[,1:2])) {
#       stop("dataframes df_mean and df_all are not equal")
#     }
#     df_all = merge(df_mean,df_sd,by = c("gene","age_group"))
    
#     df_all$age_group = as.numeric(df_all$age_group)
#     df_all$y_min = df_all$e_gc - df_all$sd_e_gc
#     df_all$y_max = df_all$e_gc + df_all$sd_e_gc
    
#     if(!dir.exists(fig_dir)) {
#       dir.create(fig_dir)
#     }
    
#     if(is.null(ylim_max)) {
#       ylim_max = max(feat_traj$e_gc_mean)
#     }
    
#     if(is.null(ylim_min)) {
#       ylim_min = min(feat_traj$e_gc_mean)
#     }
#     return(list(df_all = df_all,ylim_max = ylim_max,ylim_min = ylim_min))
#   })
  
#   df_all_ls = lapply(c(1:length(out_ls)),function(i) {
#     a = out_ls[[i]]
#     df_all = a$df_all
    
#     traj_n = paste0(rep("I",i),collapse = "")
    
#     df_all$traj = as.character(rep(traj_n,nrow(df_all)))
    
#     return(df_all)
#   })
#   df_all = do.call(rbind,df_all_ls)
  
#   ylim_min_ls = lapply(out_ls,function(a) {
#     return(a$ylim_min)
#   })
#   ylim_min = min(unlist(ylim_min_ls))
  
#   ylim_max_ls = lapply(out_ls,function(a) {
#     return(a$ylim_max)
#   })
#   ylim_max = max(unlist(ylim_max_ls))
  
#   for (gene in genes) {
#     df_f = df_all[df_all$gene == gene,]
#       print(head(df_f))
#     gene_name = gsub(";","_",gene)
#     gene_name = gsub("/","_",gene_name)
#     if(!show_sd) {
#       out = ggplot(data = df_f,aes(x = age_group,y = e_gc,group = traj)) +
#         xlab("") +
#         ylab("") +
#         scale_linetype_manual(values = c("solid","dashed","dotted","twodash")) +
#         geom_line(size = 1,aes(linetype = traj)) +
#         geom_point(size = 2) + 
#         ggtitle(label = gene) + 
#         theme(legend.position="none") +
#         theme(panel.background = element_rect(fill = "white",color = "black")) +
#         theme(axis.text = element_text(size = 20),
#               axis.title = element_text(size = 20),
#               plot.title = element_text(size = 60,hjust = 0.5))
      
#     } else {
#       out = ggplot(data = df_f,aes(x = age_group,y = e_gc,group = traj)) +
#         xlab("") +
#         ylab("") +
#         geom_line(size = 1,aes(group = traj)) + 
#         geom_errorbar(aes(ymax = y_max,ymin = y_min),size = 1.3,width = 0.5) +
#         ggtitle(label = gene) + theme(legend.position="none") + 
#         theme(panel.background = element_rect(fill = "white",color = "black")) +
#         theme(axis.text = element_text(size = 20),
#               axis.title = element_text(size = 20),
#               plot.title = element_text(size = 100,hjust = 0.5))
      
#     }
#     if(add_gridline) {
#       out = out +
#         geom_vline(xintercept = 5, linetype = "dashed",size = 2) +
#         geom_vline(xintercept = 10, linetype = "dashed",size = 2)
#     }
#     if(is.null(tag)) {
#       fn = sprintf("%s/%s",fig_dir,gene_name)
#     } else {
#       fn = sprintf("%s/%s_%s",fig_dir,gene_name,tag)
#     }
#     if(show_sd){
#       fn =  sprintf("%s_with_sd.pdf",fn)
#     } else {
#       fn =  sprintf("%s_without_sd.pdf",fn)
#     }
#     if (plot_pdf) {
#       ggsave(plot = out,filename = fn,width = 10,height = 7)
#     } else {
#       fn = gsub(pattern = ".pdf",replacement = ".png",x = fn)
#       ggsave(plot = out,filename = fn,width = 10,height = 7,bg = "transparent")
#     } 
#   } 
# }

# comp_e_gc_along_trajectory = function(mc_prob,e_gc) {
#     mc_prob = t(t(mc_prob)/colSums(mc_prob))
#     if (is.vector(e_gc)) {
#         e_gc_mean = e_gc %*% mc_prob
#         e_gc_sd = apply(X = mc_prob,MARGIN = 2,FUN = function(p) {
#         e_gc_m = sum(e_gc * p) 
#         var_e_gc = e_gc - e_gc_m
#         var_e_gc = var_e_gc ^ 2
#         var_e_gc = sum(var_e_gc * p)
#         sd_e_gc = sqrt(var_e_gc)
#         return(sd_e_gc)
#     })
#   } else {
#         e_gc_mean = e_gc %*% mc_prob
#         e_gc_sd = apply(X = mc_prob,MARGIN = 2,FUN = function(p) {
#         e_gc_m = (e_gc %*% p)[,1]
#         var_e_gc = e_gc - e_gc_m
#         var_e_gc = var_e_gc ^ 2
#         var_e_gc = (var_e_gc %*% p)[,1]
#         sd_e_gc = sqrt(var_e_gc)
#         return(sd_e_gc)
#     })
#   }
#   return(list(e_gc_mean = e_gc_mean, e_gc_sd = e_gc_sd))
# }


# my_mctnetwork_propogate_from_t <- function (mct, mcf, t, mc_p) 
# {
#     max_t = ncol(mct@mc_t)
#     step_m = list()
#     probs = matrix(0, nrow = nrow(mct@mc_t), ncol = ncol(mct@mc_t))
#     probs[, t] = mc_p
#     if (t > 1) {
#         for (i in (t - 1):1) {
#             step_m[[i]] = Matrix(t(t(as.matrix(mcf@mc_backward[[i]])) * 
#                 probs[, i + 1]), sparse = T)
#             probs[, i] = as.matrix(mcf@mc_backward[[i]]) %*% 
#                 probs[, i + 1]
#         }
#     }
#     if (t < max_t) {
#         for (i in (t + 1):max_t) {
#             step_m[[i - 1]] = Matrix(as.matrix(mcf@mc_forward[[i - 
#                 1]]) * probs[, i - 1], sparse = T)
#             probs[, i] = t(probs[, i - 1]) %*% as.matrix(mcf@mc_forward[[i - 
#                 1]])
#         }
#     }
#     return(list(probs = probs, step_m = step_m))
# }

#  fig3_plot_line_graphs(plot_pdf = F)












