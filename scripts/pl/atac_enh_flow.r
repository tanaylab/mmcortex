library(metacell)
library(lpsymphony)
library(sparseMatrixStats)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

mc_rna = scdb_mc('pl_cort')

prom = scdb_mat('pl_prom_cort')
# gb = scdb_mat('pl_gb_cort')

peak = scdb_mat('pl_ig_cort')

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN/CfuPN',
                'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))

mat = scdb_mat('pl_cort')

feats = scdb_gset('pl_f')

color_key = unique(mcmd[,c('st', 'color')])
col_annot = tibble::column_to_rownames(mcmd[,c('mc', 'st')], 'mc')
# rownames(col_annot) = mcmd$mc
ann_colors = list(st = setNames(color_key$color, color_key$st))
# ann_colors[['st']] = ann_colors[['st']][order(match(names(ann_colors[['st']]), cust_st_ord))]

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mct = scdb_mctnetwork('pl_cort')

flow_res = readRDS('./data/pl_flow_res.rds')

flow_mat = flow_res$flow_mat
cl_all = flow_res$cl_all
mc_from_mcl_flow = flow_res$mc_from_mcl_flow
colnames(mc_from_mcl_flow) = 1:ncol(mc_from_mcl_flow)

genes_both = intersect(rownames(mc_from_mcl_flow), rownames(mc_rna@e_gc))
length(genes_both)

dim(peak@mat)

peak_max_qs = sparseMatrixStats::rowMaxs(peak@mat)
quantile(peak_max_qs, c(seq(0,0.9,0.1),0.95,0.975,0.99,0.995,0.999,1))
peak_mat = peak@mat[peak_max_qs >= quantile(peak_max_qs, 0.9),]

write(rownames(peak_mat), './data/peak_locs.txt')

# peak_rs = Matrix::rowSums(peak@mat)

# peak_rs2 = Matrix::rowSums(peak_mat)

# quantile(peak_rs, c(seq(0,0.9,0.1),0.95,0.975,0.99,0.995,0.999,1))

# quantile(peak_rs2, c(seq(0,0.9,0.1),0.95,0.975,0.99,0.995,0.999,1))

# cs = Matrix::colSums(peak@mat)

# cs2= Matrix::colSums(peak_mat)

# rs = Matrix::rowSums(peak@mat)

# rs2 = Matrix::rowSums(peak_mat)

# quantile(rs2, seq(0,1,0.05))

# quantile(cs2, seq(0,1,0.05))

# quantile(cs, seq(0,1,0.05))

MIN_COV = 1
min_cov_g = rowSums(mc_from_mcl_flow) > MIN_COV

quantile(colSums(mc_from_mcl_flow), seq(0,1,l=21))

cs = colSums(mc_from_mcl_flow)
gene_folds = t(t(mc_from_mcl_flow[min_cov_g,])/cs)

fp_reg = 1e-05
flow_fp = (gene_folds + fp_reg)/(apply(gene_folds, 1, median) + fp_reg)

feats = scdb_gset('pl_f')
marks = scdb_gset('pl_cort_marks_f')
marks = names(marks@gene_set)


### Get the right featrures
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

dir.create('./figs/pl_cort_mcatac_mcrna_cor_by_day')

purrr::walk(seq_along(mc_day), function(mcs, n,i) {
    mc_ord = mcs[[i]][order(mcmd$st[mcmd$mc %in% mcs[[i]]])];
    mc_ord = cust_mc_ord_st[cust_mc_ord_st %in% mc_ord]
    # new_df = data.frame(setNames(col_annot$st[mc_ord], mc_ord))
    # colnames(new_df) = 'CellType'
    cor_mat = tgs_cor(mc_from_mcl_flow[feats_filt,mc_ord], mc_rna@e_gc[feats_filt,mc_ord], spearman = T)
    rownames(cor_mat) = mc_ord
    colnames(cor_mat) = mc_ord
    p = pheatmap::pheatmap(cor_mat, cluster_cols=F, cluster_rows = F, main = paste0('day ', n[[i]]), 
                        show_rownames = F,
                        show_colnames = F,
                        annotation_row = col_annot,
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors,
#                         silent = T,
                        color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000),
                      );
    save_pheatmap_png(p, glue::glue('./figs/pl_cort_mcatac_mcrna_cor_by_day/{n[[i]]}.png'))
    }, mcs = mc_day, n = names(mc_day)
      )

# gb_mcl = data.frame(t(tgs_matrix_tapply(as.matrix(gb@mat), cl_all, mean)))
# gb_mc = as.matrix(gb_mcl) %*% apply(flow_mat, 2, function(x) x/sum(x))

# gint = intersect(rownames(gb_mc), rownames(mc_rna@e_gc))
# length(gint)

# cor_gb_rna = tgs_cor(gb_mc[gint,],mc_rna@e_gc[gint,], spearman = T)

# cor_gb_rna_genes = tgs_cor(t(gb_mc[gint,]),t(mc_rna@e_gc[gint,]), spearman = T)

# genes_to_plot = head(diag(cor_gb_rna_genes)[order(diag(cor_gb_rna_genes), decreasing = T)], 20)
# genes_to_plot

# genes_to_plot = marks_filt

# options(repr.plot.height = 6, repr.plot.width = 6)
# sapply(seq_along(names(genes_to_plot)), function(g,gt,i) {
# #     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
# #     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
# #     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
# #     png(glue::glue('./figs/cor_atac_rna_gene_bodies/{i}_{g[[i]]}.png'))
# #     png(paste0('./figs/cor_atac_rna_genes/',i,'_',g[[i]],'.png'))
#     plot(gb_mc[g[[i]],], log2(mc_rna@e_gc[g[[i]],] + 1e-06), ylab = 'log2 RNA e_gc', xlab = 'ATAC e_gc',
# #          xlim = c(min_val,max_val),
#                                                                                      col = mcmd$color, pch = 16, main = paste(g[[i]],round(gt[[i]],3),sep=' - '));
#                                       points(gb_mc[g[[i]],], log2(mc_rna@e_gc[g[[i]],] + 1e-06), col = 'black');                     
# #     dev.off()
#     }, g = names(genes_to_plot), gt = genes_to_plot
#       )

# mc_ord = mcs[[i]][order(mcmd$st[mcmd$mc %in% mcs[[i]]])];
# mc_ord = cust_mc_ord_st[cust_mc_ord_st %in% mc_ord]
# new_df = data.frame(setNames(annotation_col[mc_ord,], mc_ord))
# colnames(new_df) = 'CellType'
# cor_mat = tgs_cor(mc_from_mcl_flow[feats_filt,], log2(mc_rna@e_gc[feats_filt,] + 1e-7), spearman = T)
# rownames(cor_mat) = mc_ord
# colnames(cor_mat) = mc_ord
# p = pheatmap::pheatmap(cor_mat[cust_mc_ord_st,cust_mc_ord_st], cluster_cols=F, cluster_rows = F, 
#                     show_rownames = F,
#                     show_colnames = F,
#                     annotation_row = col_annot,
#                     annotation_col = col_annot, 
#                     annotation_colors = ann_colors,
# #                         silent = T,
#                     color = colorRampPalette(c('blue','green','yellow','red'))(1000),
#                   );
# save_pheatmap_png(p, glue::glue('./figs/pl_cort_mcatac_mcrna_cor_by_day/{n[[i]]}.png'))

peak_max_qs = sparseMatrixStats::rowMaxs(peak@mat)

# quantile(peak_max_qs, c(seq(0,0.9,0.1),0.95,0.975,0.99,0.995,0.999,1))

peak_mat = peak@mat[peak_max_qs >= quantile(peak_max_qs, 0.9),]

peak_mcl = data.frame(t(tgs_matrix_tapply(peak_mat, cl_all, mean)))

peak_mc = as.matrix(peak_mcl) %*% apply(flow_mat, 2, function(x) x/sum(x))

saveRDS(peak_mc, './data/pl_cort_peak_mc.rds')
mg = scdb_mgraph('pl_cort')
mgg = mg@mg
pmc_sm = sapply(sort(unique(mgg$mc1)), function(mci) apply(pmc[,mgg$mc2[mgg$mc1 == mci]], 1, mean))
saveRDS(pmc_sm, './data/pl_cort_peak_mc_smoothed_mg.rds')
peak_mc = readRDS('./data/pl_cort_peak_mc.rds')

# par(mfrow=c(1,2))
peak_norm = t(apply(peak_mc, 1, function(x) x/sum(x)))
peak_ent = apply(peak_norm, 1, function(x) sum(-log(x)*x))
# quantile(apply(peak_mc, 1, mean),seq(0,1,l=21))

peak_ent = peak_ent[!is.na(peak_ent)]

head(peak_norm)

top_peak_per_mc = apply(peak_norm, 2, function(x) rownames(peak_norm)[head(order(x, decreasing = T), 500)])

dim(top_peak_per_mc)

length(unique(as.character(top_peak_per_mc)))

readr::write_tsv(as.data.frame(top_peak_per_mc), './data/top_peak_per_mc.tsv')

which_max_peak_norm = apply(peak_norm, 1, which.max)

peak_st = t(tgs_matrix_tapply(peak_mc, mcmd$st, mean))

peak_st_norm = t(apply(peak_st, 1, function(x) x/sum(x)))