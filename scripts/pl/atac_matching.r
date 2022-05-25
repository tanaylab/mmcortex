library(metacell)
library(lpsymphony)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

K = 32

mc_rna = scdb_mc('pl_cort')

prom_sum = scdb_mat('pl_prom_cort')
gb_sum = scdb_mat('pl_gb_cort')

mcmd = vroom::vroom('./BonevCollab/mcmd_pl.tsv')

mcell_add_gene_stat(mat_id = 'pl_prom_cort', gstat_id = 'pl_prom_cort', force = F)
gstat = scdb_gstat('pl_prom_cort')

# mcell_add_gene_stat(mat_id = 'pl_gb_cort', gstat_id = 'pl_gb_cort', force = F)
# gstat = scdb_gstat('pl_gb_cort')

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

# png('./figs/atac_feat_select.png', h=1000,w=1000,r=150)
plot(log(gstat$ds_mean), gstat$ds_log_varmean)

points(x[sort(unlist(top_q_inds))], y[sort(unlist(top_q_inds))], col = 'red', pch = 16)

# dev.off()

plot(log(gstat$ds_mean), gstat$sz_cor)
y = gstat$sz_cor[init_filt]
points(x[sort(unlist(top_q_inds))], y[sort(unlist(top_q_inds))], col = 'red', pch = 16)


feats = intersect(rownames(gstat)[init_filt][sort(unlist(top_q_inds))], rownames(mc_rna@mc_fp))

length(feats)

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
length(feats_filt)
length(feats)

days = unique(gb_sum@cell_metadata[colnames(gb_sum@mat),'day'])
# days = unique(prom_sum@cell_metadata[colnames(prom_sum@mat),'day'])
days

day_mat = cbind(sapply(days, function(d) apply(mcmd[,grep(d, colnames(mcmd))], 1, sum)))

day_mat = t(apply(day_mat, 1, function(x) x/sum(x)))

mcmd= tibble::tibble(cbind(mcmd, day_mat))
colnames(mcmd)[colnames(mcmd) %in% 1:6] = 13:18
head(mcmd)

# make_day_mcl = function(d, mat) {
#     cond1 = mat@cell_metadata$day == d
#     cond2 = rownames(mat@cell_metadata) %in% colnames(mat@mat)
#     mat = mat@mat[,rownames(mat@cell_metadata)[cond1 & cond2]]
#     mat = as.matrix(mat)
#     aumi = scm_downsamp(mat, n = min(Matrix::colSums(mat)))
# #     cor_mat = tgs_cor(as.matrix(mat[feats_filt,]), as.matrix(mc_rna@mc_fp[feats_filt,]), spearman = T)
#     cor_mat = tgs_cor(as.matrix(aumi[feats_filt,]), log2(1e-07 + mc_rna@e_gc[feats_filt,]), spearman = T)
#     cor_km = tglkmeans::TGL_kmeans(df = cor_mat, k = K, seed = SEED)
#     # print(dim(mat))
#     # print(length(cor_km$cluster))
#     sc_mic_cl = t(tgs_matrix_tapply(mat, cor_km$cluster, mean))
# #     cor_mic_cl = tgs_cor(x = as.matrix(sc_mic_cl[feats_filt,]), y = mc_rna@mc_fp[feats_filt,], spearman = T)
#     # print(dim(sc_mic_cl))
#     # print(class(sc_mic_cl))
#     cor_mic_cl = tgs_cor(x = sc_mic_cl[feats_filt,], y = log2(1e-07 + mc_rna@e_gc[feats_filt,]), spearman = T)
#     return(list('cor_mic_cl' = cor_mic_cl, 'sc_mic_cl' = sc_mic_cl, 'cor_km' = cor_km, 'cor_mat' = cor_mat, 'mat' = mat))
# }

# day_mcls = lapply(days, function(d) make_day_mcl(d, prom_sum))
# # day_mcls = lapply(days, function(d) make_day_mcl(d, gb_sum))
# names(day_mcls) = days

# saveRDS(day_mcls, './data/pl_cort_prom_day_mcls.rds')

# saveRDS(day_mcls, './data/pl_cort_gb_day_mcls.rds')

# day_mcls = readRDS('./data/pl_cort_gb_day_mcls.rds')
day_mcls = readRDS('./data/pl_cort_prom_day_mcls.rds')

# lapply(day_mcls, function(x) sort(table(x$cor_km$cluster)))

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

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc','IPC','iCPN','iCPN_L2-3','CPN_L2-3','Stellate_L4',
                'CPN_L5_6','iCPN_L5_6','iCPN/CfuPN','iCfuPN','CthPN','SCPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))

day_mcl_cor_figs_dir = file.path(wd, 'figs', 'day_mcl_cor_mcnra')
if (!dir.exists(day_mcl_cor_figs_dir)) {dir.create(day_mcl_cor_figs_dir)}
plot_mc_pheatmap = function(cor_km, cor_mat, cor_mic_cl, day) {
                pi = pheatmap::pheatmap(cor_mic_cl[,cust_mc_ord_st], cluster_cols = F, cluster_rows = F, 
                                   main = glue::glue('day {day}'), annotation_col = annotation_col, 
                                   annotation_colors = ann_colors,
                        color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000))
                save_pheatmap_png(pi, paste0(day_mcl_cor_figs_dir, '/', day, '.png'), h=1200,w=1200)
}

lapply(seq_along(day_mcls), function(x,n,i) {
    plot_mc_pheatmap(x[[i]]$cor_km, x[[i]]$cor_mat, x[[i]]$cor_mic_cl, n[[i]])}, 
       x = day_mcls, n = names(day_mcls)
      )

# apply(sapply(day_mcls, function(x) mcmd$st[apply(x$cor_mic_cl, 1, function(x) head(order(x, decreasing = T), 5))]), 2, table)

make_mcl_mc_edges = function(d, cor_mic_cl) {
    edge_inds_mc_avail = which(mcmd[,d] > 0)
    
    edge_inds_mc = t(apply(cor_mic_cl, 2, function(x) head(order(x, decreasing = T), 10)))
#     print(dim(edge_inds_mc))
    edge_inds_mcl = t(apply(cor_mic_cl, 1, function(x) {
                                                        cor_vec = x[edge_inds_mc_avail];
                                                        inds_in_vec = head(order(cor_vec, decreasing = T), 10)
                                                        inds = edge_inds_mc_avail[inds_in_vec];
                                                        return(inds)
                                                        }
                           )
                     )

    edge_mat = matrix(0,nrow(cor_mic_cl),ncol(cor_mic_cl))

    edge_mat_mc = sapply(seq_along(1:ncol(edge_mat)), function(em, ei, i) {tmp = em[,i];
                                                                           tmp[ei[i,]] = 1;
                                                                          return(tmp)}, em = edge_mat, ei = edge_inds_mc)
    edge_mat_mc[,!(1:ncol(edge_mat_mc) %in% edge_inds_mc_avail)] = 0
    edge_mat_mcl = t(sapply(seq_along(1:nrow(edge_mat)), function(em, ei, i) {tmp = em[i,];
                                                                           tmp[ei[i,]] = 1;
                                                                          return(tmp)}, em = edge_mat, ei = edge_inds_mcl))
#     print(dim(edge_mat_mc))
#     print(dim(edge_mat_mcl))
    edge_mat_all = edge_mat_mc + edge_mat_mcl
    edge_mat_all[edge_mat_all > 1] = 1
    edge_inds = which(edge_mat_all == 1)
#     freq_mat = matrix(rep(unlist(mcmd[,d]), nrow(cor_mic_cl)), nrow(cor_mic_cl), nrow(mcmd), byrow = T)
#     edge_vals = cor_mic_cl[edge_inds] - freq_mat[edge_inds]
    edge_vals = cor_mic_cl[edge_inds]

    row_inds = matrix(rep(1:nrow(cor_mic_cl), ncol(cor_mic_cl)), nrow(cor_mic_cl), ncol(cor_mic_cl))
    col_inds = matrix(rep(1:ncol(cor_mic_cl), nrow(cor_mic_cl)), nrow(cor_mic_cl), ncol(cor_mic_cl), byrow = T)

    edge_inds_row = row_inds[edge_inds]
    edge_inds_col = col_inds[edge_inds]
#     pheatmap::pheatmap(edge_mat_mc, color = colorRampPalette(c('white', 'black'))(2),cluster_rows = F,cluster_cols = F)
#     pheatmap::pheatmap(edge_mat_mcl, color = colorRampPalette(c('white', 'black'))(2),cluster_rows = F,cluster_cols = F)
    # print(dim(edge_mat_all))
    # print(dim(annotation_col))
#     pheatmap::pheatmap(edge_mat_all, cluster_rows = F,cluster_cols = F, 
#                        annotation_col = annotation_col, annotation_colors = ann_colors,
#                        color = colorRampPalette(c('white','black'))(2))
    return(list('edge_mat_all' = edge_mat_all, 'edge_inds' = edge_inds, 'edge_vals' = edge_vals, 
                'edge_inds_row' = edge_inds_row, 'edge_inds_col' = edge_inds_col))
}

day_edges_all = lapply(seq_along(1:length(day_mcls)), function(d, dm, i) {make_mcl_mc_edges(d[[i]], dm[[i]]$cor_mic_cl)}, d = names(day_mcls), dm = day_mcls)

edge_vals = unlist(sapply(day_edges_all, function(x) x$edge_vals))

cor_mic_cl = do.call("rbind", lapply(day_mcls, function(x) x$cor_mic_cl))

edge_inds_row = do.call("c", lapply(seq_along(1:length(day_edges_all)), 
                                                        function(dea,K,i) {return(dea[[i]]$edge_inds_row + (i-1)*K)}, 
                                                                            dea = day_edges_all, 
                                                                            K = K
                                   )
                       )

edge_inds_col = do.call("c", lapply(day_edges_all, function(x) x$edge_inds_col))



num_edges = length(edge_vals)
N = nrow(cor_mic_cl)
M = ncol(cor_mic_cl)
num_edges
N
M

min_cap = diag(N + num_edges + 4*M)

source_flow = cbind(diag(N), matrix(0, nrow = N, ncol = num_edges), matrix(0, nrow = N, ncol = 4*M))  

atac_cons_mat = matrix(0, nrow = nrow(cor_mic_cl), ncol = num_edges)
atac_cons_mat = sapply(seq_along(1:nrow(atac_cons_mat)), function(acm, i) {
                                                                x = acm[i,];
                                                                x[edge_inds_row == i] = -1;
                                                                return(x)
                                                                }, acm = atac_cons_mat
                      )
atac_cons_mat = t(atac_cons_mat)

atac_cons = cbind(diag(N), atac_cons_mat, matrix(0, nrow = N, ncol = 4*M))  

mc_flow_in = matrix(0, nrow = ncol(cor_mic_cl), ncol = num_edges)
mc_flow_in = sapply(seq_along(1:nrow(mc_flow_in)), function(acm, i) {
                                                                x = acm[i,];
                                                                x[edge_inds_col == i] = 1;
                                                                return(x)
                                                                }, acm = mc_flow_in
                      )
mc_flow_in = t(mc_flow_in)

mc_cons = cbind(matrix(0, nrow = M, ncol = N),  mc_flow_in,  
    t(apply(-1 * diag(M), 2, function(x) matrix(x, nrow=4, ncol=length(x), byrow=TRUE)))
            )

source_maxflow = c(rep(1, N), rep(0, num_edges), rep(0, 4*M))  

mc_flowcaps = cbind(matrix(0, nrow=4*M, ncol = N), matrix(0, nrow = 4*M, ncol = num_edges), diag(4*M))

lhs = rbind(min_cap, source_flow, atac_cons, mc_cons, source_maxflow, mc_flowcaps)

tbl_mc_rna = table(mc_rna@mc)
p_j = tbl_mc_rna/sum(tbl_mc_rna)

p_j_vec = as.numeric(unlist(lapply(p_j, function(x) x*c(0.8, 0.2, 0.2, 1e+08))))

cl_sizes = do.call("c", lapply(day_mcls, function(x) x$cor_km$size))

p_i = cl_sizes/sum(cl_sizes)

# rhs = c(rep(0, N+num_edges+4*M), rep(1/N, N), rep(0, N), rep(0, M), 1, p_j_vec)
rhs = c(rep(0, N+num_edges+4*M), p_i, rep(0, N), rep(0, M), sum(p_i), p_j_vec)

dir = c(rep(">=", N+num_edges+4*M), rep("==", N), rep("==", N), rep("==", M), "==", rep("<=", 4*M))

types = rep('C', ncol(lhs))
length(types)
head(types)

K2 = 1e+01
K1 = 1e+00

# obj = c(rep(0, N), match(edge_vals,sort(edge_vals))/length(edge_vals), unlist(sapply(p_j, function(p) 10*p*c(-K2, -K1, K1, K2))))
# obj = c(rep(0, N), order(edge_vals)/length(edge_vals), unlist(sapply(p_j, function(p) 10*p*c(-K2, -K1, K1, K2))))
# obj = c(rep(0, N), -exp(edge_vals), rep(c(-K2, -K1, K1, K2), M))
# obj = c(rep(0, N), -edge_vals, unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))))
obj = c(rep(0, N), -edge_vals, rep(c(-K2, -K1, K1, K2), M))

# length(obj)
# length(types)
# length(rhs)
# length(dir)
# dim(lhs)

# class(obj)
# class(lhs)
# class(dir)
# class(rhs)
# class(types)

sol = lpsymphony_solve_LP(obj, lhs, dir, rhs, types = types, max = FALSE, )

sol$objval

source_sol = sol$solution[1:N]
atac_sol = tapply(sol$solution[(N+1):(N + num_edges)], edge_inds_row, list)
mc_sol = matrix(sol$solution[(N + num_edges + 1):length(sol$solution)], nrow = M, ncol = 4, byrow = T)




# quantile(as.numeric(mc_sol)*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))), seq(0.05,1,0.05))

# hist(as.numeric(mc_sol)*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))))

# sum(as.numeric(sol$solution[(N + num_edges + 1):length(sol$solution)])*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))))

# length(edge_vals)

# length(unlist(atac_sol))

# quantile(-exp(edge_vals)*as.numeric(sol$solution[(N+1):(N + num_edges)]), seq(0.05,1,0.05))

# sum(-exp(edge_vals)*as.numeric(sol$solution[(N+1):(N + num_edges)]))

# hist(-exp(edge_vals)*as.numeric(unlist(atac_sol)))

flow_mat = matrix(0, nrow(cor_mic_cl), ncol(cor_mic_cl))
for (nm in names(atac_sol)) {
    nm = as.numeric(nm)
    flow_mat[as.numeric(nm),edge_inds_col[edge_inds_row == as.numeric(nm)]] = atac_sol[[nm]]
}

colnames(flow_mat) = 1:ncol(cor_mic_cl)
rownames(flow_mat) = 1:nrow(cor_mic_cl)

tapply(apply(mc_sol, 1, sum), mcmd$st, sum)

tapply(p_j, mcmd$st, sum)

tapply(apply(mc_sol, 1, sum), mcmd$st, sum)/tapply(p_j, mcmd$st, sum)

png('./figs/pl_sum_flow_to_mc.png', res = 250, height = 1600, width = 1600)

plot(as.numeric(p_j), as.numeric(apply(flow_mat, 2, sum)), ylab = 'sum flow to metacell', xlab = 'metacell fraction',
    pch = 16, col = mcmd$color, bg='black')
points(as.numeric(p_j), as.numeric(apply(flow_mat, 2, sum)), col='black')
x = seq(0, 0.04, l=50)
lines(x, x, lty='dashed')
dev.off()

cl_all = do.call('c', lapply(seq_along(1:length(day_mcls)), function(dm, i) {dm[[i]]$cor_km$cluster + (i-1)*K}, dm = day_mcls))

sc_mic_cl_full = data.frame(t(tgs_matrix_tapply(as.matrix(prom_sum@mat), cl_all, mean)))
mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% apply(flow_mat, 2, function(x) x/sum(x))
# mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% flow_mat

mcl_flow_norm = mc_from_mcl_flow/colSums(mc_from_mcl_flow)

saveRDS(object = list('flow_mat' = flow_mat, 'cl_all' = cl_all, 'mc_from_mcl_flow' = mc_from_mcl_flow), file = './data/pl_flow_res.rds')

gb = intersect(rownames(mc_from_mcl_flow), rownames(mc_rna@e_gc))

# cor_atac_rna = tgs_cor(log2(mcl_flow_norm[gb,] + 1e-06), log2(mc_rna@e_gc[gb,] + 1e-06), spearman = T)
cor_atac_rna = tgs_cor(mc_from_mcl_flow[feats_filt,], mc_rna@mc_fp[feats_filt,], spearman = T)

rownames(cor_atac_rna) = colnames(cor_atac_rna)

options(repr.plot.height = 12, repr.plot.width = 14)
p = pheatmap::pheatmap(cor_atac_rna[cust_mc_ord_st,cust_mc_ord_st], cluster_cols = F, cluster_rows = F, 
                   annotation_row = annotation_col,
                   annotation_col = annotation_col, 
                   annotation_colors = ann_colors,
                       color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000),
                  )
save_pheatmap_png(p, './figs/pl_cor_mc_atac_rna.png')
options(repr.plot.height = 6, repr.plot.width = 6)

cor_atac_rna_genes = tgs_cor(t(mcl_flow_norm[gb,]), t(log2(mc_rna@e_gc[gb,] + 1e-06)), spearman = T)

rownames(cor_atac_rna_genes) = colnames(cor_atac_rna_genes)

flow_com = apply(flow_mat[,cust_mc_ord_st], 1, function(x) sum(x*(1:length(x)))/sum(x))

p = pheatmap::pheatmap(flow_mat[order(flow_com),cust_mc_ord_st], 
#                        border_color = 'black',
                       show_colnames = F, 
                       show_rownames = T, 
                       cluster_cols = F, 
                       cluster_rows = F,
#                        gaps_col = tail(sapply(cust_st_ord, function(u) min(which(names(cust_mc_ord_st) == u))), -1) - 1,
                    color = colorRampPalette(c('black', 'red', 'orange','yellow', 'white'))(100),
                   annotation_col = annotation_col, 
                   annotation_colors = ann_colors)

save_pheatmap_png(p,'./figs/pl_mcl_mc_flow_mat.png')

flow_mat_st = t(tgs_matrix_tapply(flow_mat[order(flow_com),cust_mc_ord_st], names(cust_mc_ord_st), sum))
options(repr.plot.height = 12)
p = pheatmap::pheatmap(flow_mat_st[,cust_st_ord], 
#                        border_color = 'black',
                       show_colnames = T, 
                       show_rownames = T, 
                       cluster_cols = F, 
                       cluster_rows = F,
                       fontsize_col = 24,
                       fontsize_row = 8,
#                        gaps_col = tail(sapply(cust_st_ord, function(u) min(which(names(cust_mc_ord_st) == u))), -1) - 1,
                    color = colorRampPalette(c('black', 'red', 'orange','yellow', 'white'))(100),
#                    annotation_col = annotation_col, 
#                    annotation_colors = ann_colors
                      )
options(repr.plot.height = 8)

save_pheatmap_png(p,'./figs/pl_flow_mat_st.png', height = 2200, width= 1200)

# genes_to_plot = head(diag(cor_atac_rna_genes)[order(diag(cor_atac_rna_genes), decreasing = T)], 50)
# genes_to_plot

# genes_to_plot2 = head(diag(cor_atac_rna_genes)[order(diag(cor_atac_rna_genes))], 20)
# genes_to_plot2

# dir.create('./figs/cor_atac_rna_genes')

# names(genes_to_plot)

# genes_to_plot

# options(repr.plot.height = 6, repr.plot.width = 6)
# sapply(seq_along(names(genes_to_plot)), function(g,gt,i) {
# #     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
# #     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
# #     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
# #     png(glue::glue('./figs/cor_atac_rna_genes/{i}_{g[[i]]}.png'))
#     png(paste0('./figs/cor_atac_rna_genes/',i,'_',g[[i]],'.png'))
#     plot(log2(mcl_flow_norm[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), 
#          ylab = 'log2 RNA e_gc', xlab = 'log2 ATAC e_gc',
#          cex.lab = 2, cex.main = 2,
# #          xlim = c(min_val,max_val),
#          col = mcmd$color, pch = 16, main = paste(g[[i]],round(gt[[i]],3),sep=' - '));
#     points(log2(mcl_flow_norm[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), col = 'black');
# #     axis(2,cex.axis = 3)
# #     axis(1,cex.axis = 3)
#     dev.off()
#     }, g = names(genes_to_plot), gt = genes_to_plot
#       )


# grep('Myt1l', rownames(mc_rna@e_gc),v=T)

# g = 'Rnd2'
# # png(paste0('./figs/cor_atac_rna_genes/', g, '.png'))
# options(repr.plot.height = 8, repr.plot.width = 8)
# plot(mcl_flow_norm[g,], log2(mc_rna@e_gc[g,] + 1e-06), 
#      ylab = 'log2 RNA e_gc', xlab = 'ATAC e_gc',
#               cex.lab = 2, cex.main = 2,
#      col = mcmd$color, pch = 16, main = paste(g,round(diag(cor_atac_rna_genes)[[g]],3),sep=' - '))
# points(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), col = 'black')

# dev.off()

# options(repr.plot.height = 6, repr.plot.width = 6)
# sapply(names(genes_to_plot2), function(g) {
# #     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
# #     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
# #     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
#     plot(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), 
# #          xlim = c(min_val,max_val),
#                                                                                      col = mcmd$color, pch = 16, main = g); 
#                                       points(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), col = 'black');
#                                     }
    #   )


