library(metacell)
library(lpsymphony)
library(sparseMatrixStats)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

mc_rna = scdb_mc('cort6')

prom = scdb_mat('prom_cort')
gb = scdb_mat('gb_cort')
peak = scdb_mat('peak_cort')
mcmd = vroom::vroom('./BonevCollab/mcmd_cort_6.tsv')


mat = scdb_mat('cort6')

feats = scdb_gset('cort5_feats_f')

km_feats = readRDS('./data/km_feats.rds')

clust_sum = tgs_matrix_tapply(t(mc_rna@e_gc[names(feats@gene_set),]), km_feats$cluster, sum)

clust_sum_tf = readRDS('./data/clust_sum_tf.rds')

dim(clust_sum_tf)

eg_louv = readRDS('./data/flow_louv_cl.rds')

pca_feat_cl = princomp(x = t(clust_sum))

# options(repr.plot.height = 12, repr.plot.width = 12)
# png('./figs/feat_clust_pca.png', h = 1200, w = 1200, res = 150)
plot(pca_feat_cl$scores[,1], pca_feat_cl$scores[,3], col = scales::alpha(mcmd$color, 0.5), pch = 16, xlab = 'PC1 - 39.8% of variance', ylab = 'PC2 - 17.5% of variance',
     main = 'PCA of feature gene cluster values',
     cex = 0.8,
    )
# points(pca_feat_cl$scores[,1], pca_feat_cl$scores[,2], col = 'black', 
#        cex = 0.8,
#       )
# dev.off()
# options(repr.plot.height = 6, repr.plot.width = 6)

# options(repr.plot.height = 12, repr.plot.width = 12)
# png('./figs/feat_clust_pca.png', h = 1200, w = 1200, res = 150)
plot(pca_feat_cl$scores[,1], pca_feat_cl$scores[,2], col = scales::alpha(mcmd$color, 0.5), pch = 16, xlab = 'PC1 - 39.8% of variance', ylab = 'PC2 - 17.5% of variance',
     main = 'PCA of feature gene cluster values',
     cex = 0.8,
    )
# points(pca_feat_cl$scores[,1], pca_feat_cl$scores[,2], col = 'black', 
#        cex = 0.8,
#       )
# dev.off()
# options(repr.plot.height = 6, repr.plot.width = 6)

pca_feat_cl$sdev/sum(pca_feat_cl$sdev)

color_key = unique(mcmd[,c('st', 'color')])
annotation_col = data.frame(CellType = mcmd$st)
rownames(annotation_col) = mcmd$mc
ann_colors = list(CellType = setNames(color_key$color, color_key$st))

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mct = scdb_mctnetwork('cort6_comp')

flow_res = readRDS('./data/flow_res.rds')

flow_mat = flow_res$flow_mat
cl_all = flow_res$cl_all
mc_from_mcl_flow = flow_res$mc_from_mcl_flow
colnames(mc_from_mcl_flow) = 1:ncol(mc_from_mcl_flow)

genes_both = intersect(rownames(mc_from_mcl_flow), rownames(mc_rna@e_gc))
length(genes_both)

quantile(rowSums(mc_from_mcl_flow), seq(0,1,l=21))

MIN_COV = 10
min_cov_g = rowSums(mc_from_mcl_flow) > MIN_COV

quantile(colSums(mc_from_mcl_flow), seq(0,1,l=21))

cs = colSums(mc_from_mcl_flow)
gene_folds = mc_from_mcl_flow[min_cov_g,]/cs

fp_reg = 1e-05
flow_fp = (gene_folds + fp_reg)/(apply(gene_folds, 1, median) + fp_reg)

table(apply(flow_fp, 1, function(x) any(abs(log2(x)) > 2.5)))

feats = scdb_gset('cort5_feats_f')

marks = sort(rownames(flow_fp)[apply(flow_fp, 1, function(x) any(abs(log2(x)) > 2))])

marks_filt = marks[marks %in% names(feats@gene_set)]
# marks_filt

cust_st_ord = c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC','iCfuPN',
                      'iCPN/CfuPN','iCPN_L2-3','CPN_L2-3','CthPN','SCPN','CPN_L5-6')
cust_mc_ord_st = unlist(sapply(cust_st_ord, function(st) sort(mcmd$mc[mcmd$st == st])))
# cust_mc_ord_st

options(repr.plot.height = 8, repr.plot.width = 14)

plot_mat = log2(mc_rna@mc_fp[marks_filt,cust_mc_ord_st])
min_val = min(plot_mat)
max_val = max(plot_mat)
p = pheatmap::pheatmap(plot_mat, cluster_rows = T, cluster_cols = F, 
                   annotation_col = annotation_col, show_rownames = T, show_colnames = F,
                  annotation_colors = ann_colors, fontsize_row = 5,
                       breaks = c(seq(min_val, 0,l=50), seq(0.01,max_val,l=51)),
                      color = colorRampPalette(c('blue','white', 'red'))(100),
#                        color = colorRampPalette(c('white', 'yellow', 'blue', 'brown', 'black'))(1000) 
                      )
options(repr.plot.height = 6, repr.plot.width = 6)
save_pheatmap_png(p, './figs/rna_mc_hm.png')

# genes = 
genes_ord_by_rna = p$tree_row$labels[p$tree_row$order]

options(repr.plot.height = 8, repr.plot.width = 14)
plot_mat = log2(flow_fp[genes_ord_by_rna,cust_mc_ord_st])
min_val = min(plot_mat)
max_val = max(plot_mat)
p = pheatmap::pheatmap(plot_mat, cluster_rows = F, cluster_cols = F, 
                    annotation_col = annotation_col, show_rownames = T, show_colnames = F,
                    annotation_colors = ann_colors, fontsize_row = 5,
                    breaks = c(seq(min_val, 0,l=50), seq(0.01,max_val,l=51)),
                    color = colorRampPalette(c('blue','white', 'red'))(100),
#                       color = colorRampPalette(c('white', 'red', 'black'))(100),
                      )
options(repr.plot.height = 6, repr.plot.width = 6)
save_pheatmap_png(p, './figs/atac_mc_hm.png')

### Get the right featrures
gstat = scdb_gstat('test')
x = log(gstat$ds_mean)
init_filt = which(x >= -5)
y = gstat$ds_log_varmean[init_filt]
name_filt = gstat$name
x = x[init_filt]

xcut = cut(x, breaks = seq(min(x), max(x), l = 20))
top_q_inds = lapply(levels(xcut), function(l) {inds = which(xcut == l); 
                                  xfilt = y[inds]; 
                                  xtop = head(inds[order(xfilt, decreasing = T)], 30); 
                                  return(xtop)
                                 }
      )

names(top_q_inds) = levels(xcut)
feats = intersect(rownames(gstat)[init_filt[unlist(top_q_inds)]], rownames(mc_rna@mc_fp))
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

mc_day = lapply(tail(sort(unique(mat@cell_metadata$day)), -1), function(d) sort(unique(mc_rna@mc[rownames(mat@cell_metadata)[mat@cell_metadata$day == d & 
                                                                                             rownames(mat@cell_metadata) %in% names(mc_rna@mc)]])))

names(mc_day) = tail(sort(unique(mat@cell_metadata$day)), -1)

sapply(seq_along(mc_day), function(mcs, n,i) {
    mc_ord = mcs[[i]][order(mcmd$st[mcmd$mc %in% mcs[[i]]])];
    mc_ord = cust_mc_ord_st[cust_mc_ord_st %in% mc_ord]
    new_df = data.frame(setNames(annotation_col[mc_ord,], mc_ord))
    colnames(new_df) = 'CellType'
    cor_mat = tgs_cor(mc_from_mcl_flow[feats_filt,mc_ord], mc_rna@e_gc[feats_filt,mc_ord], spearman = T)
    rownames(cor_mat) = mc_ord
    colnames(cor_mat) = mc_ord
    p = pheatmap::pheatmap(cor_mat, cluster_cols=F, cluster_rows = F, main = paste0('day ', n[[i]]), 
                        show_rownames = F,
                        show_colnames = F,
                        annotation_row = new_df,
                        annotation_col = new_df, 
                        annotation_colors = list('CellType' = ann_colors$CellType[unique(new_df$CellType)]),
                        silent = T,
                        color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000),
                      );
    save_pheatmap_png(p, glue::glue('./figs/mc_cor_day_{n[[i]]}.png'))
    }, mcs = mc_day, n = names(mc_day)
      )

gb_mcl = data.frame(t(tgs_matrix_tapply(as.matrix(gb@mat), cl_all, mean)))
gb_mc = as.matrix(gb_mcl) %*% apply(flow_mat, 2, function(x) x/sum(x))

gint = intersect(rownames(gb_mc), rownames(mc_rna@e_gc))
length(gint)

cor_gb_rna = tgs_cor(gb_mc[gint,],mc_rna@e_gc[gint,], spearman = T)

cor_gb_rna_genes = tgs_cor(t(gb_mc[gint,]),t(mc_rna@e_gc[gint,]), spearman = T)

genes_to_plot = head(diag(cor_gb_rna_genes)[order(diag(cor_gb_rna_genes), decreasing = T)], 20)
genes_to_plot

options(repr.plot.height = 6, repr.plot.width = 6)
sapply(seq_along(names(genes_to_plot)), function(g,gt,i) {
#     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
#     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
#     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
#     png(glue::glue('./figs/cor_atac_rna_gene_bodies/{i}_{g[[i]]}.png'))
#     png(paste0('./figs/cor_atac_rna_genes/',i,'_',g[[i]],'.png'))
    plot(log2(gb_mc[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), ylab = 'log2 RNA e_gc', xlab = 'log2 ATAC e_gc',
#          xlim = c(min_val,max_val),
                                                                                     col = mcmd$color, pch = 16, main = paste(g[[i]],round(gt[[i]],3),sep=' - '));
                                      points(log2(gb_mc[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), col = 'black');                     
#     dev.off()
    }, g = names(genes_to_plot), gt = genes_to_plot
      )


options(repr.plot.height = 12, repr.plot.width = 14)
pheatmap::pheatmap(cor_gb_rna[cust_mc_ord_st,cust_mc_ord_st], cluster_rows=F, cluster_cols=F,annotation_colors = ann_colors,
                  annotation_col = annotation_col, annotation_row = annotation_col)
options(repr.plot.height = 6, repr.plot.width = 6)

peak_max_qs = sparseMatrixStats::rowMaxs(peak@mat)

quantile(peak_max_qs, c(seq(0,0.9,0.1),0.95,0.975,0.99,0.995,0.999,1))

peak_mat = peak@mat[peak_max_qs >= quantile(peak_max_qs, 0.9),]

peak_mcl = data.frame(t(tgs_matrix_tapply(as.matrix(peak_mat), cl_all, mean)))
peak_mc = as.matrix(peak_mcl) %*% apply(flow_mat, 2, function(x) x/sum(x))

 mc_ag = table(mc_rna@mc,mat@cell_metadata[names(mc_rna@mc),"day"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)

late_mcs_per_traj = lapply(eg_louv, function(x) x[x %in% which(mc_ag_cn[,'18'] >= 0.33)])

lapply(late_mcs_per_traj, function(x) mcmd$st[x])

get_peak_flow_across_clusters = function(cl) {
#     mcs_clust_11 = which(mc_max_clust == cl & mc_ag_cn[,tail(colnames(mc_ag_cn), 1)] >= 0.4)
    mcs_clust_11 = cl
    if (length(mcs_clust_11) > 0) {
        probs_mc_11 = mc_ag_cn[mcs_clust_11,tail(colnames(mc_ag_cn), 1)]
        probs_flow_11 = mctnetwork_propogate_from_t(mct = mct, 
                                                    t = 6, 
                                                    mc_p = ifelse(1:nrow(mc_ag_cn) %in% mcs_clust_11, probs_mc_11, 0))

        probs_flow_11_norm = probs_flow_11$probs/colSums(probs_flow_11$probs)

        peak_flow_11 = peak_mc %*% probs_flow_11_norm
        return(peak_flow_11)
    }
    else {return(NULL)}
}

flows_across_louv = lapply(late_mcs_per_traj, get_peak_flow_across_clusters)

length(unlist(tail(flows_across_louv, -1)))

flow_arr = array(tail(flows_across_louv, -1), dim = c(nrow(flows_across_louv[[2]]), 
                                                      ncol(flows_across_louv[[2]]), length(tail(flows_across_louv, -1))))

dim(flow_arr)

class(flow_arr)

ls()

max_flow_arr = apply(flow_arr, 1, max)

head(flow_arr[,,1])

lapply(tail(flows_across_louv, -1), function(x) do.call()

lapply(flows_across_louv, head)

mc_max_clust = readRDS('./data/mc_max_clust_k=24.rds')

### Reassign small mc_max_clust clusters
mc_max_clust[mc_max_clust == 1] = 5
mc_max_clust[mc_max_clust == 2] = 5
mc_max_clust[mc_max_clust == 3] = 5
mc_max_clust[mc_max_clust == 10] = 9
mc_max_clust[mc_max_clust == 13] = 6
mc_max_clust[mc_max_clust == 14] = 6
mc_max_clust[mc_max_clust == 16] = 15
mc_max_clust[mc_max_clust == 22] = 19

mcs_clust_11 = which(mc_max_clust == 11 & mc_ag_cn[,'18'] >= 0.4)
mcs_clust_11

probs_mc_11 = mc_ag_cn[mcs_clust_11,'18']
probs_mc_11

probs_flow_11 = mctnetwork_propogate_from_t(mct = mct, t = 6, mc_p = ifelse(1:nrow(mc_ag_cn) %in% mcs_clust_11, probs_mc_11, 0))

# probs_flow_11 = lapply(mcs_clust_11, function(mci) mctnetwork_propogate_from_t(mct = mct, t = 6, mc_p = probs_mc_11[[as.character(mci)]]))

probs_flow_11_norm = probs_flow_11$probs/colSums(probs_flow_11$probs)

peak_flow_11 = peak_mc %*% probs_flow_11_norm

high_var_peaks = rownames(peak_flow_11)[head(order(apply(peak_flow_11 + 1e-06, 1, function(x) 
                                                    log2(max(x)/min(x))), decreasing = T), 1000)]
head(high_var_peaks)

head(peak_flow_11[high_var_peaks,])

col_pal = sample(colors()[grep('white|gray|grey', colors(), inv=T)], 20)

hi_exp_hi_var = which(apply(peak_flow_11, 1, max) >= 0.25 & apply(peak_flow_11+1e-06, 1, function(x) log2(max(x)/min(x))) >= 2)

quantile(apply(peak_flow_11, 1, max), seq(0,1,l=21))

plot(0, xlim = c(12.5,18.5), ylim = c(-10,-1))
sapply(sample(hi_exp_hi_var, 10), function(x) {
    lines(13:18, log2(peak_flow_11[x,] + 1e-06), col = sample(colors(), 1))
})

get_peak_flow_across_clusters = function(cl) {
    mcs_clust_11 = which(mc_max_clust == cl & mc_ag_cn[,tail(colnames(mc_ag_cn), 1)] >= 0.4)
    if (length(mcs_clust_11) > 0) {
        probs_mc_11 = mc_ag_cn[mcs_clust_11,tail(colnames(mc_ag_cn), 1)]
        probs_flow_11 = mctnetwork_propogate_from_t(mct = mct, 
                                                    t = 6, 
                                                    mc_p = ifelse(1:nrow(mc_ag_cn) %in% mcs_clust_11, probs_mc_11, 0))

        probs_flow_11_norm = probs_flow_11$probs/colSums(probs_flow_11$probs)

        peak_flow_11 = peak_mc %*% probs_flow_11_norm
        return(peak_flow_11)
    }
    else {return(NULL)}
}

as.numeric(sort(unique(mc_max_clust)))

table(mc_max_clust)

peak_flow_clust_ts = lapply(as.numeric(sort(unique(mc_max_clust))), get_peak_flow_across_clusters)

peak_flow_clust_hi_exp_hi_var = lapply(peak_flow_clust_ts, function(x) {
    if (!(is.null(x))) {
        which(apply(x, 1, max) >= 0.15 & apply(x+1e-06, 1, function(y) log2(max(y)/min(y))) >= log2(10))
    }
})

names(peak_flow_clust_ts) = as.numeric(sort(unique(mc_max_clust)))

names(peak_flow_clust_hi_exp_hi_var) = names(peak_flow_clust_ts)

lapply(peak_flow_clust_hi_exp_hi_var, head)

head(peak_flow_clust_hi_exp_hi_var[[5]])

unq_dyn_enh = sort(unique(unlist(sapply(peak_flow_clust_hi_exp_hi_var, function(x) names(x)))))
put_dyn_enh_mat = matrix(0, nrow = length(unq_dyn_enh), ncol = length(peak_flow_clust_hi_exp_hi_var))
rownames(put_dyn_enh_mat) = unq_dyn_enh
colnames(put_dyn_enh_mat) = names(peak_flow_clust_hi_exp_hi_var)
for (i in 1:length(peak_flow_clust_hi_exp_hi_var)) {
    put_dyn_enh_mat[names(peak_flow_clust_hi_exp_hi_var[[i]]),i] = 1 
}

dim(put_dyn_enh_mat)

table(apply(put_dyn_enh_mat, 1, sum))

put_enh_clust_dist = as.matrix(dist(t(put_dyn_enh_mat), method = 'binary'))

put_enh_dist = as.matrix(dist(put_dyn_enh_mat, method = 'binary'))

put_enh_km = tglkmeans::TGL_kmeans(put_enh_dist, k = 32, seed = SEED)

put_enh_km$size[order(put_enh_km$size, decreasing = T)]

table(apply(put_dyn_enh_mat[put_enh_km$cluster == 28,], 1, which.max))

head(put_dyn_enh_mat[put_enh_km$cluster == 28,])

peak_flow_clust_ts[[1]]['chr1-100737142-100737642',]

col_pal = sample(grep('grey|gray|white', colors(), inv=T, v=T), length(peak_flow_clust_ts))
col_pal

enh = 'chr1-150115374-150115874'
plot(0, xlim = c(12.5,18.5), ylim = c(-20,0))
purrr::walk(seq_along(peak_flow_clust_ts), function(x, cp, i) {
    lines(13:18, log2(x[[i]][enh,] + 1e-06), col = cp[[i]]);
    legend(18.5, -4, names(peak_flow_clust_ts), fill = col_pal, xpd=T)}, 
       x = peak_flow_clust_ts, cp = col_pal)

pp = pheatmap::pheatmap(put_enh_dist[order(put_enh_km$cluster), order(put_enh_km$cluster)], 
                        color = colorRampPalette(c('black', 'white'))(100),
                        silent = T,
                   show_rownames = F, show_colnames = F,
                   cluster_cols = F, cluster_rows = F)

save_pheatmap_png(pp, './figs/put_enh_bin_dist.png')

p = pheatmap::pheatmap(1-put_enh_clust_dist, color = colorRampPalette(c('white', 'red', 'black'))(100))
save_pheatmap_png(p, './figs/fclust_bin_dist.png', w = 600, h = 600)

head(put_dyn_enh_mat)

peak_mc_fclust = t(tgs_matrix_tapply(peak_mc, mc_max_clust, mean))

peak_mc_fclust_norm = t(apply(peak_mc_fclust, 1, function(x) x/sum(x)))

min(apply(peak_mc_fclust_norm, 1, max))

head(peak_mc_fclust_norm)

head(peak_mc_fclust)

mctnetwork_propogate_from_t()

dim(peak_mc_fclust)

pmm = rowMeans(peak_mc)
pmv = rowVars(peak_mc)
pmn = pmv/(pmm**2)
plot(log2(pmm), log2(pmn))

smoothScatter(log2(pmm**2), log2(pmn))

cor_peak_tf_cl = tgs_cor(t(peak_mc), t(clust_sum_tf), spearman = T)

dim(cor_peak_tf_cl)

hist(apply(cor_peak_tf_cl, 1, function(x) max(abs(c(min(x), max(x))))))

quantile(apply(cor_peak_tf_cl, 1, function(x) max(abs(c(min(x), max(x))))), seq(0,1,l=21))

hist(apply(cor_peak_tf_cl, 1, min))

head(cor_peak_tf_cl)

head(cor_peak_tf_cl[order(cor_peak_tf_cl[,2], decreasing = T),])

head(enh_cor_tf_cl_ord_dec)

enh_cor_tf_cl_ord_dec = apply(cor_peak_tf_cl,2, function(x) order(x, decreasing = T))
enh_cor_tf_cl_ord_inc = apply(cor_peak_tf_cl,2, function(x) order(x, decreasing = F))


enh_cor_tf_cl_ord_dec[,4][[1]]

rownames(cor_peak_tf_cl)[enh_cor_tf_cl_ord_dec[,4][[1]]]

cl = 13
plot(peak_mc[rownames(cor_peak_tf_cl)[enh_cor_tf_cl_ord_dec[,cl][[3]]],], clust_sum_tf[cl,], col=mcmd$color, pch = 16)

cor_peak_feat_cl = tgs_cor(t(peak_mc), t(clust_sum), spearman = T)

dim(cor_peak_feat_cl)

hist(apply(cor_peak_feat_cl, 1, function(x) max(abs(c(min(x), max(x))))))

quantile(apply(cor_peak_feat_cl, 1, function(x) max(abs(c(min(x), max(x))))), seq(0,1,l=21))

hist(apply(cor_peak_feat_cl, 1, min))

head(cor_peak_feat_cl)



library(ggplot2)
library(dplyr)
# library(tidyr)
# library(viridis)

x = log2(peak_mc[1,] + 1e-06)
x_df = data.frame(cbind(x,mcmd$st, mcmd$color))
colnames(x_df) = c('x', 'st', 'color')
x_df[,'x'] = as.numeric(x_df[,'x'])

# x_df

plot_den_by_st = function(x) {
    x_df = data.frame(cbind(x,mcmd$st, mcmd$color))
    colnames(x_df) = c('x', 'st', 'color')
    x_df[,'x'] = as.numeric(x_df[,'x'])
#     print(setNames(unique(x_df$color), unique(x_df$st)))
    p <- x_df %>%
#         mutate(st = forcats::fct_reorder(st, x)) %>% # Reorder data
        ggplot( aes(x=st, y=x, fill=st)) +
            geom_violin(width=2.1, size=0.2, alpha = 0.5) + 
            scale_fill_manual(values = setNames(unique(x_df$color), unique(x_df$st))) + 
            scale_color_manual(values = 'black') + 

#             scale_fill_viridis(discrete=TRUE) +
#             scale_color_viridis(discrete=TRUE) +
#             theme_ipsum() +
#             theme(
#               legend.position="none"
#             ) +
            coord_flip() + # This switch X and Y axis and allows to get the horizontal version
            xlab("") +
            ylab("log2 ATAC") + 
            ylim(-10, 1)
    p
}

plot_den_by_st(log2(peak_mc[1000,] + 1e-06))

gfilt = rownames(mc_rna@e_gc)[apply(mc_rna@e_gc, 1, max) >= 5e-05]

cor_peak_rna = tgs_cor(t(peak_mc), t(mc_rna@e_gc[gfilt,]), spearman = T)

dim(cor_peak_rna)

my_tfs = c('Sox5','Nfib', 'Pou3f2')
inds_tfs = match(my_tfs, gfilt)
inds_tfs

cor_peak_ord = order(cor_peak_rna,decreasing = T)

head(cor_peak_ord)

m = nrow(cor_peak_rna)

r = ((cor_peak_ord-1) %% m) + 1
c = floor((cor_peak_ord-1) / m) + 1

head(r)

head(c)

top_peak_mat = head(cbind(r[c %in% inds_tfs],c[c %in% inds_tfs]), 250)
head(top_peak_mat)

peak_mc_st = tgs_matrix_tapply(peak_mc, mcmd$st, max)

peak_mc_st_norm = t(apply(peak_mc_st, 2, function(x) x/sum(x)))

head(peak_mc_st_norm)

head(peak_mc_st_norm[apply(peak_mc_st_norm, 1, max) >= 0.33,])

peak_ord_by_st = order(apply(peak_mc_st_norm, 1, max), decreasing = T)
peak_mc_st_norm = peak_mc_st_norm[peak_ord_by_st,]

dir.create('./figs/enh_bar')

plot_st_enh_bar = function(st) {
    ipc_peaks = which(apply(peak_mc_st_norm, 1, which.max) == which(colnames(peak_mc_st_norm) == st))
    st_name = gsub('\\/', '-', st)
#     dir.create(glue::glue('./figs/enh_bar/{st_name}'))
    sapply(head(names(ipc_peaks), 5), function(p) {
        png(glue::glue('./figs/enh_bar/{st_name}/{p}.png'), width = 1500, height = 700, res = 75)
        barplot(peak_mc[p,order(mcmd$st)], col = mcmd$color[order(mcmd$st)], main = p, xlab = 'mc', ylab = 'ATAC e_gc')
        dev.off()
    })
}

sapply(unique(mcmd$st), plot_st_enh_bar)

ipc_peaks = which(apply(peak_mc_st_norm, 1, which.max) == which(colnames(peak_mc_st_norm) == 'CthPN'))
head(ipc_peaks)

options(repr.plot.width = 15)
sapply(head(names(ipc_peaks), 20), function(p) barplot(peak_mc[p,order(mcmd$st)], col = mcmd$color[order(mcmd$st)], main = p))
# barplot(peak_mc['chr12-114267268-114267768',order(mcmd$st)], col = mcmd$color[order(mcmd$st)])
options(repr.plot.width = 6)

tbl1 = table(colnames(peak_mc_st_norm)[apply(peak_mc_st_norm[apply(peak_mc_st_norm, 1, max) >= 0.4,], 1, which.max)])
tbl1/sum(tbl1)

tbl2 = table(colnames(peak_mc_st_norm)[apply(peak_mc_st_norm, 1, which.max)])
tbl2/sum(tbl2)

hist(apply(peak_mc_st_norm, 1, max), 100)

dir.create('./figs/enh_examples')

options(repr.plot.height = 6, repr.plot.width = 6)
sapply(seq_along(1:nrow(head(top_peak_mat, 20))), function(pm,egc,tpm,cpr,i) {

    #     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
#     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
#     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
    g = rownames(egc)[tpm[i,2]]
    p = rownames(peak_mc)[tpm[i,1]]
#     png(glue::glue('./figs/cor_peak_rna/{i}_{g}_{p}.png'))
#     png(paste0('./figs/cor_peak_rna/',i,'_',g[[i]],'.png'))
    plot(log2(pm[tpm[i,1],] + 1e-06), log2(egc[tpm[i,2],] + 1e-06), ylab = 'log2 RNA e_gc', xlab = 'log2 ATAC e_gc',
#          xlim = c(min_val,max_val),
         xlim = c(-8,0),
            col = mcmd$color, pch = 16, 
            main = paste(g,p,cpr[tpm[i,1],tpm[i,2]],sep=' - '));
    points(log2(pm[tpm[i,1],] + 1e-06), log2(egc[tpm[i,2],] + 1e-06), col = 'black');                     
#     dev.off()
    }, pm = peak_mc, egc = mc_rna@e_gc[gfilt,],tpm = top_peak_mat, cpr = cor_peak_rna
      )


mc_ag = table(mc_rna@mc,mat@cell_metadata[names(mc_rna@mc),"day"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)

late_mcs = which(mc_ag_cn[,ncol(mc_ag_cn)] > 0.33)

late_sts =mcmd$st[late_mcs]
table(late_sts)

cn = as.numeric(colnames(mc_ag_cn))
t_com = apply(as.matrix(mc_ag_cn), 1, function(x) sum(x*cn)/sum(x))

calc_egc_traj = function(mc) {
    p_mc = rep(0,ncol(mc_rna@e_gc))
    p_mc[mc] = mc_ag_c[mc,ncol(mc_ag_c)]
    probs_trans = mctnetwork_propogate_from_t(mct, ncol(mc_ag_c), p_mc)
    mc_prob = t(t(probs_trans$probs)/colSums(probs_trans$probs))
    egc_traj = mc_rna@e_gc %*% mc_prob
    return(egc_traj)
}

traj_ls = lapply(late_mcs, calc_egc_traj)

traj_cor = tgs_cor(t(traj_ls[[1]][apply(traj_ls[[1]], 1, max) >= 5e-05,]), spearman = F)

traj_km = tglkmeans::TGL_kmeans(traj_cor, k = 16)

traj_km$size

length(traj_km$cluster)

km_mean = tgs_matrix_tapply(apply(traj_ls[[1]][apply(traj_ls[[1]], 1, max) >= 5e-05,], 1, function(x) x/max(x)), traj_km$cluster, mean)

pheatmap::pheatmap(km_mean, cluster_rows = F, cluster_cols = F)

traj_filt = traj_ls[[1]][apply(traj_ls[[1]], 1, max) >= 5e-05,]

length(sort(rownames(traj_filt)[traj_km$cluster == 7 & rowMaxs(traj_filt)/rowMins(traj_filt + 1e-08) >= 3]))

# goi = c('Qser1','Fzd3', 'Smarcd3', 'Syt11', 'Maz')
goi = head(sort(rownames(traj_filt)[traj_km$cluster == 7 & rowMaxs(traj_filt)/rowMins(traj_filt + 1e-08) >= 3]), 20)
plot(0,ylim = c(-20,-10), xlim = c(13,18))
# par(mar = c(5,4,4,10))
sapply(1:20, function(g, i) {lines(13:18, log2(traj_filt[g[[i]],]), col = col_pal[[i]], lwd = 2);
                                  legend(13, -9, legend = g, fill = col_pal, cex = 0.5)}, g = head(goi, 20))

head(traj_cor)

hist(log10(apply(traj_ls[[1]], 1, max) + 1e-08))

lapply(head(traj_ls), head)

names(traj_ls) = late_mcs

T_lfc = 3
traj_deg = lapply(traj_ls, function(x) rownames(x[apply(x, 1, function(y) log2((max(y) + 1e-06)/(min(y) + 1e-06))) > T_lfc,]))

lapply(head(traj_deg), head)

traj_detf = lapply(traj_deg, function(x) x[x %in% tfoi])

col_pal = sample(colors()[grep('white|gray|grey', colors(), inv=T)], 20)
col_pal

sapply(95:100, function(td, n, i) {
    tf_mat = log2(traj_ls[[i]][td[[i]],] + 1e-06);
    tf_mat = tf_mat[order(apply(tf_mat, 1, mean), decreasing = T),]
    par(mar = c(5, 4, 4, 7) + 0.1)
    plot(0,xlim = c(13,18),ylim = c(min(tf_mat), max(tf_mat)), main = paste(n[[i]], mcmd$st[as.numeric(n[[i]])]));
    sapply(1:nrow(tf_mat), function(j) lines(13:18, tf_mat[j,], col = col_pal[[j]], lwd = 1.5));
    legend(18.2, 1.05*max(tf_mat), legend = rownames(tf_mat), fill = col_pal[1:nrow(tf_mat)], xpd = T)
}, td = traj_detf, n = names(traj_detf))

length(traj_deg[[1]])

tfoi = read.delim('./data/tfoi.txt') %>% unlist %>% sort %>% as.character

sample(tfoi, 10)

sapply(sample(tfoi, 10), function(tf) plot(13:18,log2(egc_traj[tf,] + 1e-06), main = tf))

pheatmap::pheatmap(mc_prob, cluster_cols = F, cluster_rows = F)

pheatmap::pheatmap(probs_trans$probs, cluster_cols = F, cluster_rows = F)

options(repr.plot.height = 12, repr.plot.width = 8)
pheatmap::pheatmap(as.matrix(mc_ag_cn)[order(t_com),], cluster_rows = F, cluster_cols = F, 
                   annotation_row = annotation_col, annotation_colors = ann_colors,
                  color = colorRampPalette(c('white', 'red', 'black'))(100))
options(repr.plot.height = 6, repr.plot.width = 6)

tail(order(t_com))


