ls()

# devtools::load_all('~/src/metacell')
library(metacell)
devtools::load_all('~/src/metacell.flow')

wd = '/home/feshap/raid/proj/mmcortex'
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scdb_flow_init()
SEED = 1337
K = 16
set.seed(SEED)
scfigs_init("figs/")

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

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
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

col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')

ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))
# ann_colors[['cell_type']] = ann_colors[['cell_type']][order(match(names(ann_colors[['cell_type']]), cust_st_ord))]

# nn <- mct@network

# nn$flow <- mcf@edge_flows

# nn <- nn[nn$flow > 0 & as.numeric(nn$mc1) > 0 & as.numeric(nn$mc2) > 0 & nn$mc1 != nn$mc2,]

# mg <- mgraph@mgraph

# mg_ig <- igraph::graph_from_data_frame(mg, directed = F)

# mg_louv <- igraph::cluster_louvain(graph = mg_ig, weights = 1/mg$dist)

# mg_louv

# library(princurve)

# princurve::principal_curve

# gr_louv <- igraph::communities(mg_louv)

# nn_avg <- tidyr::drop_na(do.call('rbind', lapply(gr_louv, function(x) as.data.frame(list(x0 = mean(mc2d@mc_x[as.numeric(x)]), 
#                                            y0 = mean(mc2d@mc_y[as.numeric(x)]),
#                                             x1 = sum(mc2d@mc_x[nn$mc2[nn$mc1 %in% as.numeric(x)]]*nn$flow[nn$mc1 %in% as.numeric(x)])/sum(nn$flow[nn$mc1 %in% as.numeric(x)]),
#                                             y1 = sum(mc2d@mc_y[nn$mc2[nn$mc1 %in% as.numeric(x)]]*nn$flow[nn$mc1 %in% as.numeric(x)])/sum(nn$flow[nn$mc1 %in% as.numeric(x)]),
#                                             flow = sum(nn$flow[nn$mc1 %in% x]))))))

# dim(nn)

# nn_avg

# quantile(mc2d@mc_x[as.numeric(nn$mc1)])
# quantile(mc2d@mc_x[as.numeric(nn$mc2)])
# quantile(mc2d@mc_y[as.numeric(nn$mc1)])
# quantile(mc2d@mc_y[as.numeric(nn$mc2)])

# head(mc2d@mc_x[as.numeric(nn$mc1)])
# head(mc2d@mc_x[as.numeric(nn$mc2)])
# head(mc2d@mc_y[as.numeric(nn$mc1)])
# head(mc2d@mc_y[as.numeric(nn$mc2)])

# options(repr.plot.width = 20)
# options(repr.plot.height = 20)

# nn_avg$flow

# png('./figs/avg_flow_on_mc2d.png', h = 1200, w = 1200)
# plot(mc2d@mc_x, mc2d@mc_y, col = mcmd$color, cex = 2, pch = 16)
# x0 = mc2d@mc_x[as.numeric(nn$mc1)]
# x1 = mc2d@mc_x[as.numeric(nn$mc2)]
# y0 = mc2d@mc_y[as.numeric(nn$mc1)]
# y1 = mc2d@mc_y[as.numeric(nn$mc2)]
# ms <- (y1 - y0)/(x1 - x0)
# arrows(x0 = nn_avg$x0,y0 = nn_avg$y0, x1 = nn_avg$x1,y1 = nn_avg$y1, lwd = 15*sqrt(nn_avg$flow),cex = 1.5 )
# dev.off()# arrows(x0 = x0,y0 = y0, x1 = x1,y1 = y1, lwd = 15*sqrt(nn$flow),cex = 1.5 )

# ls()

# sapply(c('Sox2', 'Pax6'), function(g) mcell_mc2d_plot_gene('pl_cort', g))

# tgconfig::get_package_params('metacell')

# par(cex.lab = 2)
# mcell_mc_plot_gg('pl_cort', 'Bcl11b', 'Satb2', fig_fn = './figs/pl_cort_Bcl11b_vs_Satb2.png', use_egc = T, cex = 3, height = 1400, width = 1400)

legc = log2(1e-05 + mc@e_gc)

feats = scdb_gset('pl_filt_lat')
feats = names(feats@gene_set)
feats = feats[feats %in% rownames(mc@e_gc)]

legc_st = tgs_matrix_tapply(legc[feats,], mcmd$cell_type, mean)

distt = as.matrix(dist(legc_st))

options(repr.plot.width = 8)
options(repr.plot.height = 8)

pheatmap::pheatmap(distt[cust_st_ord,cust_st_ord], cluster_cols = F, cluster_rows = F)

# color_key = color_key[c('Astrocytes','early_NSC?','OPC','NSC','nNSC?','early_nNSC?','IPC','iCfuPN','iCPN_L2-3',
#               'iCPN_L5-6','CPN_L2-3','CPN_L5-6','SCPN','CthPN','Stellate_L4'),]
df = data.frame(color_key[order(match(color_key$cell_type, cust_st_ord)),])

# st_for_plot = df$st
# st_for_plot['iCPN_L2-3'] = 'immature CPN L2-3'
# st_for_plot['iCfuPN'] = 'immature CPN L2-3'


l = nrow(df)
scale_y = 2
png('./figs/legend_pl.png', width = 1200, height = 900, res=250)
plot(rep(0.93,l), scale_y*seq(l,1,-1), pch = 16, cex = 1, col = df$color, ylim = c(0.5,scale_y*l+1),
    xlab = '', 
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(rep(0.94,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 0.5, df$cell_type)
dev.off()

# top_mc_fp_gene = apply(mc@mc_fp, 1, which.max)
# top_fp_gene = apply(mc@mc_fp, 1, max)
# top_mc_fp_gene = apply(legc, 1, which.max)
# top_fp_gene = apply(legc, 1, max)
# max_genes_by_st = lapply(cust_st_ord, function(st) setNames(top_fp_gene[top_mc_fp_gene %in% mcmd$metacell[mcmd$cell_type == st]],
#                                                            names(top_fp_gene)[top_mc_fp_gene %in% mcmd$metacell[mcmd$cell_type == st]]))
fp_by_st = t(tgs_matrix_tapply(mc@mc_fp, mcmd$cell_type, mean))

top_fp_gene = apply(fp_by_st, 1, max)

head(top_fp_gene)

max_genes_by_st = lapply(cust_st_ord, function(st) setNames(top_fp_gene[top_st_fp_gene == st],
                                                           names(top_fp_gene)[top_st_fp_gene == st])
                         )
names(max_genes_by_st) = cust_st_ord

max_genes_by_st

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

mat <- scdb_mat('pl_cort')
mc2d <- scdb_mc2d('pl_cort')

colnames(mc@e_gc) <- 1:ncol(mc@e_gc)

colnames(mc@mc_fp) <- 1:ncol(mc@mc_fp)

scdb_add_mc(id = 'pl_cort', mc)

# mcell_mc2d_plot_by_factor('pl_cort', 'pl_cort', meta_field = 'day', single_plot = T)

# mcell_mc2d_plot_by_factor('pl_cort', 'pl_cort', meta_field = 'day', single_plot = F)

options(repr.plot.width = 18)
options(repr.plot.height = 18)

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
pp = pheatmap::pheatmap(pltmt, color = colorRampPalette(pltt)(l), breaks = brks, fontsize = 14,
                        cluster_cols = F, cluster_rows = F, 
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors, 
                        show_colnames = F, fontsize_row = 14)

save_pheatmap_png(pp, './figs/pl_cort_goi_hm_for_rp_new.png', h = 2200, w = 3800,r=300)

# st_fp <- as.data.frame(t(tgs_matrix_tapply(mc@mc_fp, mcmd$cell_type, mean)))

# cs <- Matrix::colSums(mat@mat)

# summary(mc_umis)
# summary(mcn_umis)

# gs <- scdb_gset('pl_filt_lat')

# c('Pcna', 'Mki67', 'Top2a') %in% names(gs@gene_set)
# gs@gene_set[c('Pcna', 'Mki67', 'Top2a')]

# dror_path <- "/net/mraid14/export/tgdata/users/drorba/ambient/out/mmcortex_new/"

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

cor_md_by_st = apply(head(legc), 1, function(x) purrr::map2(.x = x, .y = mcmd$mean_day, .f = tapply, ... = mcmd$cell_type))

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

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(legc))
# tfs

length(tfs)

tfs_hi <- tfs[matrixStats::rowMaxs(legc[tfs,]) >= -13 & 
              matrixStats::rowMaxs(legc[tfs,]) - matrixStats::rowMins(legc[tfs,]) >= 3]

length(tfs_hi)

'Zfhx4' %in% tfs_hi

tfs_pos_cor <- lapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = T), 20))
# names(tfs_pos_cor) <- c('NSC', 'IPC', 'IPC_cyc')
tfs_neg_cor <- lapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) head(sort(cor_gene_legc_w_md_in_st[tfs_hi,ct], decreasing = F), 20))
# names(tfs_neg_cor) <- c('NSC', 'IPC', 'IPC_cyc')

library(pheatmap)

mcs_ord <- unlist(sapply(c('NSC', 'IPC', 'IPC_cyc'), function(ct) {
    f <- mcmd$cell_type == ct; 
    return(which(f)[order(mcmd$mean_day[f])])}))

ppp <- pheatmap(legc[unique(c(names(unlist(tfs_pos_cor)), names(unlist(tfs_neg_cor)))),mcs_ord], 
                color = colorRampPalette(c('white', 'gold', 'red2','black'))(100),
         cluster_cols = F, cluster_rows = T,
         annotation_col = col_annot, 
         annotation_colors = ann_colors)

save_pheatmap_png(ppp, './figs/hi_var_tfs_nsc_ipc.png')

ls()

ls()

unique(names(unlist(tfs_pos_cor)))

unique(names(unlist(tfs_neg_cor)))

tfs_neg_cor

intersect(unique(names(unlist(tfs_pos_cor))),unique(names(unlist(tfs_neg_cor))))

head(cor_gene_legc_w_md_in_st[order(rowMeans(cor_gene_legc_w_md_in_st)),], 20)

head(cor_gene_legc_w_md_in_st[order(rowMeans(cor_gene_legc_w_md_in_st), decreasing = T),], 20)

cor_legc_md <- tgs_cor(t(legc), as.matrix(mcmd$mean_day), spearman = T)

abs_diff_cor_ct_cor_all <- abs(rowMeans(subset(cor_gene_legc_w_md_in_st, select = -c(Astrocytes, Oligodendrocytes))) - cor_legc_md)
abs_diff_cor_ct_cor_all <- setNames(abs_diff_cor_ct_cor_all[,1], rownames(abs_diff_cor_ct_cor_all))

head(sort(abs_diff_cor_ct_cor_all, decreasing = T))

cor_gene_legc_w_md_in_st[names(head(sort(abs_diff_cor_ct_cor_all, decreasing = T), 20)),]

cor_gene_legc_w_md_in_st[c('Ptx3', 'Cux1'),]

plot(mcmd$mean_day, legc['Ssbp3',], col = mcmd$color, pch = 16,
#      ylim = c(-13,-9), 
     cex = 3)

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

tfoi

heatmap_genes_along_trajectory = function(mc_prob,mc,fig_dir,tag,min_max_fold = 2,min_thr = -13,main = "",show_rownames = T,
                                          highlighted_genes = NULL) {
  
  
  legc = log2(mc@e_gc + 1e-7)
  
  mc_prob = t(t(mc_prob)/colSums(mc_prob))
  legc_traj = legc %*% mc_prob
  min_traj = apply(legc_traj,1,min)
  max_traj = apply(legc_traj,1,max)
  
  cond = max_traj- min_traj > min_max_fold & max_traj > min_thr & rownames(legc_traj) %in% tfoi
  
  gnms_f = rownames(legc_traj)[cond]
    print(length(gnms_f))
  legc_traj_n = (legc_traj[gnms_f,] - min_traj[gnms_f])/(max_traj[gnms_f]-min_traj[gnms_f])
  
  legc_traj_nn = legc_traj_n/rowSums(legc_traj_n)
  
  gene_mean_time = (legc_traj_nn %*% c(1:6))[,1]
  gene_mode = apply(legc_traj_nn,1,which.max)
  gene_rank = 10*gene_mode + gene_mean_time
    
  gene_ord = rev(order(gene_rank))
  
  legc_traj_n = legc_traj_n[gene_ord,]
  legc_traj_n = legc_traj_n[c(nrow(legc_traj_n):1),]
  
  nm_to_pos = seq(0,1,length.out = nrow(legc_traj_n))
  names(nm_to_pos) = rownames(legc_traj_n)
  
  
#   rownames(legc_traj_n) = substr(rownames(legc_traj_n),1,10)
  colnames(legc_traj_n) = c(1:ncol(legc_traj_n))+12
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9,name = "PuBu"))(1000)

#   pdf(file = sprintf("%s/heatmap_genes_along_traj%s.pdf",fig_dir,tag),h = 16, w = 10,useDingbats = F)
#   par(mar = c(5,12,2,1))
    p = pheatmap::pheatmap(legc_traj_n, color = shades, cluster_rows = F, cluster_cols = F, 
                           fontsize_row = 18 - nrow(legc_traj_n)/3, main = paste(tag, mcmd$cell_type[as.numeric(tag)], collapse = '--'),
#                            nrow(legc_traj_n)/30,
                           fontsize_col = 20)
    st = mcmd$cell_type[as.numeric(tag)]
    save_pheatmap_png(p, paste0('./figs/st_traj/', st, '_', tag, '.png'), h=300+50*nrow(legc_traj_n), w = 600)
#   image(x = t(legc_traj_n),col = shades,axes = F)
  
#   ind_even = seq.int(from = 2,by = 2,length.out = (ncol(legc_traj_n) %/% 2))
#   ind_odd = seq.int(from = 1,by = 2,length.out = (ncol(legc_traj_n) %/% 2 + ncol(legc_traj_n) %% 2 ))
#   axis_pos = seq(0,1,length.out = 13)
  
#   axis(side = 1,at = axis_pos[ind_even],labels = c(1:ncol(legc_traj_n))[ind_even],lwd = 0,lwd.ticks = 2,cex.axis = 2,padj = 1)
#   axis(side = 1,at = axis_pos[ind_odd],labels = c(1:ncol(legc_traj_n))[ind_odd],lwd = 0,lwd.ticks = 2,cex.axis = 2,padj = 1)
#   if(!is.null(highlighted_genes)) {
#     highlighted_genes = intersect(names(nm_to_pos),highlighted_genes)
#     axis(side = 2,at = nm_to_pos[highlighted_genes],labels = highlighted_genes,lwd = 0,lwd.ticks = 2,cex.axis = 2,las = 2)
#   }
#   dev.off()

}

plot_heatmap_trajectory = function() {
  
  mc = scdb_mc(nm)
  mat = scdb_mat(nm)
  mct_id = nm
  mct = scdb_mctnetwork(mct_id)
  mcf = scdb_mctnetflow(mct_id)
  marks = scdb_gset('pl_cort_marks_f')
  marks = names(marks@gene_set)
#   fig3_param_ls = fig3_parameters()
  
#   highlighted_genes = fig3_param_ls$highlighted_genes
  
    
#   fig_dir = "figs/paper_figs/fig3"
#   if(!dir.exists(fig_dir)) {
#     dir.create(fig_dir)
#   }
  
#   mc_ord = read.table("data/wt10.cluster_annotation/mc_order_clust.txt",sep = "\t",stringsAsFactors = F)
#   mc_ord= mc_ord$mc
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"t"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
  late_st = c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')
  st_mcs = lapply(late_st, function(x) mcmd$metacell[mcmd$cell_type %in% x])
  late_mcs = lapply(st_mcs, function(x) x[mc_ag_cn[x,ncol(mc_ag_cn)] > 0.3])
  
  df_param = data.frame(mc = setNames(sapply(late_mcs, function(x) x[[1]]), late_st),
                          edge_w_scale = c(5e-5,5e-5,5e-5,5e-5),
                          fr_scale = c(1,1,1,1),
                          max_lwd = c(20,20,20,20),
                          lfp_thr = c(2.5,2.5,2.5,3),
                          min_thr = c(-13,-13,-13,-12))
    
  
  for (i in 1:nrow(df_param)) {
    m = df_param$mc[i]
    lfp_thr = df_param$lfp_thr[i]
    min_thr = df_param$min_thr[i]
    
#     fig_dir = sprintf("figs/paper_figs/fig3/mc_%d",m)
#     if(!dir.exists(fig_dir)) {
#       dir.create(fig_dir)
#     }
    
    p_mc = rep(0,ncol(mc@e_gc))
    p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
    
#     probs = mctnetwork_propogate_from_t(mct, ncol(mc_ag_c), p_mc)
#     which(p_mc > 0)
    probs = mctnetflow_propogate_from_t(mcf, ncol(mc_ag_c), p_mc)
    mc_prob = t(t(probs$probs)/colSums(probs$probs))
    
    markers = marks[marks %in% rownames(mc@e_gc)]
#       highlighted_genes[[as.character(m)]]
    heatmap_genes_along_trajectory(mc_prob = mc_prob,mc = mc,
                                   fig_dir = fig_dir,
                                   tag = as.character(m),
                                   min_max_fold = lfp_thr,
                                   main = sprintf("Metacell %d",m),
                                   min_thr = min_thr,show_rownames = F,
                                   highlighted_genes = markers)
    
    
  }
  
}

if (!dir.exists('./figs/st_traj')) {dir.create('./figs/st_traj')}

plot_heatmap_trajectory()
















