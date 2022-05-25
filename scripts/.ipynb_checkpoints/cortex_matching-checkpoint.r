library(metacell)
library(lpsymphony)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

K = 25

mc_rna = scdb_mc('cort5')

prom_sum = scdb_mat('prom_cort')

mcmd = vroom::vroom('./BonevCollab/mcmd_cort_5_iter1.tsv')

### ATAC gene correlation structure

# cor_atac_gene = tgs_cor(t(as.matrix(prom_sum@mat[gb,])), spearman = T)

# dim(cor_atac_gene)[[1]]/50



# atac_g_km = tglkmeans::TGL_kmeans(cor_atac_gene, k = 200, metric = 'pearson')

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

dim(prom_sum@mat)

# prom_ds = scm_downsamp(prom_sum@mat, min(Matrix::colSums(prom_sum@mat)))

# make_day_dist_mat = function(d) {
# #     prom_sum = scdb_mat(glue::glue('prom_cort_day_{d}'))
#     prom_mat = prom_ds[,rownames(prom_sum@cell_metadata)[prom_sum@cell_metadata$day == d]]
#     prom_mat = prom_mat[apply(prom_mat, 1, max) >= 4,]
#     return(as.matrix(tgs_dist(x = t(as.matrix(prom_mat)))))
# }

# make_day_cor_mat = function(d) {
# #     prom_sum = scdb_mat(glue::glue('prom_cort_day_{d}'))
#     prom_mat = prom_ds[,rownames(prom_sum@cell_metadata)[prom_sum@cell_metadata$day == d]]
#     prom_mat = prom_mat[apply(prom_mat, 1, max) >= 4,]
#     return(tgs_cor(x = as.matrix(prom_mat), spearman = T))
# }

# day_euc_dist_mats = lapply(unique(prom_sum@cell_metadata$day), make_day_dist_mat)

# day_cor_mats = lapply(unique(prom_sum@cell_metadata$day), make_day_cor_mat)

# day_km = lapply(day_euc_dist_mats, function(x) tglkmeans::TGL_kmeans(df = x, k = 25))

# day_km_cor = lapply(day_cor_mats, function(x) tglkmeans::TGL_kmeans(df = x, k = 25))

# days = unique(prom_sum@cell_metadata$day)
# names(day_cor_mats) = days
# names(day_euc_dist_mats) = days
# names(day_km) = days
# names(day_km_cor) = days

# sapply(seq_along(day_cor_mats), function(cm, km, d, i) {
#     plot_mat = cm[[i]][order(km[[i]]$cluster), order(km[[i]]$cluster)]
#     min_val = min(plot_mat)
#     max_val = max(plot_mat)
#     range = max_val - min_val
#     len = 1000
#     p = pheatmap(plot_mat, silent = T, 
#              show_colnames = F, 
#              show_rownames = F, 
#              cluster_cols = F, 
#              cluster_rows = F,
#             fontsize = 20,
#                  res = 150,
#             breaks = c(seq(min_val, min_val+0.1*range, l=len/5),
#                         seq(min_val+0.1*range+range/len, min_val+0.2*range, l=len/5),
#                         seq(min_val+0.2*range+range/len, min_val+0.3*range, l=len/5),
#                         seq(min_val+0.3*range+range/len, min_val+0.4*range, l=len/5),
#                         seq(min_val+0.4*range+range/len, min_val+range, l=len/5+1)
#                       ),
#              color = colorRampPalette(c('blue', 'green', 'yellow', 'red', 'black'))(len),
#              main = glue::glue('day = {d[[i]]}, Spearman dist., n_cells = {nrow(cm[[i]])}, k = 25'),
             
#             );
#     save_pheatmap_png(p, glue::glue('./figs/day_{d[[i]]}_scatac_dist_sp.png'),height = 2500, width = 2500)
# }, cm = day_cor_mats, km = day_km_cor, d = days)

# sapply(seq_along(day_cor_mats), function(cm, km, d, i) {
#     plot_mat = cm[[i]][order(km[[i]]$cluster), order(km[[i]]$cluster)]
#     min_val = min(plot_mat)
#     max_val = max(plot_mat)
#     range = max_val - min_val
#     len = 1000
#     p = pheatmap(plot_mat, silent = T,
#              show_colnames = F, 
#              show_rownames = F, 
#              cluster_cols = F, 
#              cluster_rows = F,
#                  fontsize = 20,
#                  res = 150,
#             breaks = c(seq(min_val, min_val+0.6*range, l=len/5),
#                         seq(min_val+0.6*range+range/len, min_val+0.7*range, l=len/5),
#                         seq(min_val+0.7*range+range/len, min_val+0.8*range, l=len/5),
#                         seq(min_val+0.8*range+range/len, min_val+0.9*range, l=len/5),
#                         seq(min_val+0.9*range+range/len, min_val+range, l=len/5+1)
#                       ),
#              color = colorRampPalette(c('blue', 'green', 'yellow', 'red', 'black'))(len),
#              main = glue::glue('day = {d[[i]]}, Euc. dist., n_cells = {nrow(cm[[i]])}, k = 25'),
#             );
#     save_pheatmap_png(p, glue::glue('./figs/day_{d[[i]]}_scatac_dist.png'),height = 2500, width = 2500)
# }, cm = day_euc_dist_mats, km = day_km, d = days)

nm = 'test'
mcell_add_gene_stat(mat_id = 'prom_cort', gstat_id = nm, force = F)


gstat = scdb_gstat('test')

png('./figs/atac_feat.png', h=1000,w=1000,r=150)
plot(log(gstat$ds_mean), gstat$ds_log_varmean)
dev.off()

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

# plot(x,y)
png('./figs/atac_feat_select.png', h=1000,w=1000,r=150)
plot(log(gstat$ds_mean), gstat$ds_log_varmean)

points(x[sort(unlist(top_q_inds))], y[sort(unlist(top_q_inds))], col = 'red', pch = 16)

dev.off()
# smoothScatter(log(gstat$ds_mean), gstat$ds_log_varmean)

# tapply(x, xcut, function(y) )

# scdb_ls('gset')

# mcell_gset_filter_cov(gstat_id = nm, gset_id = nm,T_tot = 100, T_top3 = 2)
# mcell_gset_filter_varmean(gstat_id = nm, gset_id = nm, )

# quantile(apply(as.matrix(prom_sum@mat), 1, max), seq(0,1,l=21))

# gb = intersect(rownames(prom_sum@mat[apply(as.matrix(prom_sum@mat), 1, max) >= 10,]), rownames(mc_rna@e_gc))

# feats = gb

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

days = unique(prom_sum@cell_metadata$day)
days

day_mat = cbind(sapply(days, function(d) apply(mcmd[,grep(d, colnames(mcmd))], 1, sum)))

day_mat = t(apply(day_mat, 1, function(x) x/sum(x)))

mcmd= tibble::tibble(cbind(mcmd, day_mat))
head(mcmd)

# sapply(unique(prom_sum@cell_metadata$day), function(d) {
#     ig_cells = rownames(prom_sum@cell_metadata)[prom_sum@cell_metadata$day == d]
#     mcell_mat_ignore_cells(new_mat_id = glue::glue('prom_cort_day_{d}'), mat_id = 'prom_cort', ig_cells = ig_cells,reverse = F)
# })

# make_day_mcl = function(d) {
# #     prom_sum = scdb_mat(glue::glue('prom_cort_day_{d}'))
#     prom_mat = prom_sum@mat[,rownames(prom_sum@cell_metadata)[prom_sum@cell_metadata$day == d]]
#     cor_mat = tgs_cor(as.matrix(prom_mat[feats_filt,]), as.matrix(mc_rna@mc_fp[feats_filt,]), spearman = T)
#     cor_km = tglkmeans::TGL_kmeans(df = cor_mat, k = K, seed = SEED)
#     sc_mic_cl = data.frame(t(tgs_matrix_tapply(as.matrix(prom_mat), cor_km$cluster, mean)))
#     cor_mic_cl = tgs_cor(x = as.matrix(sc_mic_cl[feats_filt,]), y = mc_rna@mc_fp[feats_filt,], spearman = T)
#     return(list('cor_mic_cl' = cor_mic_cl, 'sc_mic_cl' = sc_mic_cl, 'cor_km' = cor_km, 'cor_mat' = cor_mat, 'prom_mat' = prom_mat))
# }

# day_mcls = lapply(days, make_day_mcl)
# names(day_mcls) = days

# saveRDS(day_mcls, './data/cort_day_mcls.rds')

day_mcls = readRDS('./data/cort_day_mcls.rds')

color_key = unique(mcmd[,c('st', 'color')])
annotation_col = data.frame(CellType = mcmd$st)
rownames(annotation_col) = mcmd$mc
ann_colors = list(CellType = setNames(color_key$color, color_key$st))

cust_st_ord = c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC','iCfuPN',
                      'iCPN/CfuPN','iCPN_L2-3','CPN_L2-3','CthPN','SCPN','CPN_L5-6')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))
# cust_mc_ord_st

cust_mc_ord_st

plot_mc_pheatmap = function(cor_km, cor_mat, cor_mic_cl, day) {
#     print(head(rownames(cor_mat)))
#     cl_mc_df = data.frame(cbind(setNames(cor_km$cluster, rownames(cor_mat)), mc_rna@mc[rownames(cor_mat)]))
#     colnames(cl_mc_df) = c('cl', 'mc')
#     cl_mc_df$st = mcmd$st[cl_mc_df$mc]
#     tbl_mc = matrix(table(cl_mc_df[with(cl_mc_df, order(cl, mc)),c('cl', 'mc')]), 
#                     nrow = length(unique(cor_km$cluster)), 
#                     ncol = length(unique(mc_rna@mc)))
#     tbl_mc[is.na(tbl_mc)] = 0
#     rownames(tbl_mc) = 1:nrow(tbl_mc)
#     colnames(tbl_mc) = 1:ncol(tbl_mc)
#     print(head(cl_mc_df))
    p = pheatmap::pheatmap(cor_mic_cl[,cust_mc_ord_st], cluster_cols = F, cluster_rows = F, main = glue::glue('day {day}'),
                           annotation_col = annotation_col, annotation_colors = ann_colors,
#                            breaks = c(
# #                                         seq(min(cor_mic_cl), mean(min(cor_mic_cl), 0)/2,l=33),
# #                                         seq(mean(min(cor_mic_cl), 0)/2+0.01,mean(max(cor_mic_cl), 0)/2,l=33),
# #                                         seq(mean(max(cor_mic_cl), 0)/2+0.01, max(cor_mic_cl), l=34)
#                                     ),
                       color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000))
    save_pheatmap_png(p, glue::glue('./figs/mcl_mc_cor_day_{day}.png'))
#     print(head(t(apply(tbl_mc, 1, function(x) x/sum(x)))))
#     p2 = pheatmap::pheatmap(t(apply(tbl_mc, 1, function(x) x/sum(x))), cluster_cols = F, cluster_rows = F, 
#                                     annotation_col = annotation_col, annotation_colors = ann_colors,
#                    color = colorRampPalette(c('white', 'red', 'black'))(100))
    }

lapply(seq_along(day_mcls), function(x,n,i) {plot_mc_pheatmap(x[[i]]$cor_km, x[[i]]$cor_mat, x[[i]]$cor_mic_cl, n[[i]])}, x = day_mcls, n = names(day_mcls))

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
    print(dim(edge_mat_all))
    print(dim(annotation_col))
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

length(obj)
length(types)
length(rhs)
length(dir)
dim(lhs)

class(obj)
class(lhs)
class(dir)
class(rhs)
class(types)

sol = lpsymphony_solve_LP(obj, lhs, dir, rhs, types = types, max = FALSE, )

sol$objval

source_sol = sol$solution[1:N]
atac_sol = tapply(sol$solution[(N+1):(N + num_edges)], edge_inds_row, list)
mc_sol = matrix(sol$solution[(N + num_edges + 1):length(sol$solution)], nrow = M, ncol = 4, byrow = T)

# ## QC 

quantile(as.numeric(mc_sol)*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))), seq(0.05,1,0.05))

hist(as.numeric(mc_sol)*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))))

sum(as.numeric(sol$solution[(N + num_edges + 1):length(sol$solution)])*unlist(sapply(p_j, function(p) p*c(-K2, -K1, K1, K2))))

# length(edge_vals)

# length(unlist(atac_sol))

quantile(-exp(edge_vals)*as.numeric(sol$solution[(N+1):(N + num_edges)]), seq(0.05,1,0.05))

sum(-exp(edge_vals)*as.numeric(sol$solution[(N+1):(N + num_edges)]))

hist(-exp(edge_vals)*as.numeric(unlist(atac_sol)))

flow_mat = matrix(0, nrow(cor_mic_cl), ncol(cor_mic_cl))
for (nm in names(atac_sol)) {
    nm = as.numeric(nm)
    flow_mat[as.numeric(nm),edge_inds_col[edge_inds_row == as.numeric(nm)]] = atac_sol[[nm]]
}

colnames(flow_mat) = 1:ncol(cor_mic_cl)
rownames(flow_mat) = 1:nrow(cor_mic_cl)

# plot_flow = function(flow_mat) {
    x_mcl = rep(0,nrow(flow_mat))
    x_mc = rep(1,ncol(flow_mat))
    y_mcl = 1:nrow(flow_mat)
#     y_mc = seq(min(y_mcl), max(y_mcl), l=ncol(flow_mat))
    y_mc = 1:ncol(flow_mat)
    x1_edges = unlist(apply(flow_mat, 1, function(x) rep(0, length(which(x>0)))))
    y1_edges = unlist(sapply(seq_along(1:nrow(flow_mat)), function(fl, i) {rep(i, length(which(fl[i,]>0)))}, fl = flow_mat))
    x2_edges = unlist(apply(flow_mat, 1, function(x) rep(1, length(which(x>0)))))
    y2_edges = unlist(apply(flow_mat, 1, function(x) which(x>0)))
    y2_mcl_edge_mean = unlist(sapply(y_mcl, function(u) rep(mean(y2_edges[y1_edges == u]), length(which(y1_edges == u)))))
    y2_mcl_mean = unlist(sapply(y_mcl, function(u) mean(y2_edges[y1_edges == u])))
    

# png('./figs/flow_graph.png', height = 800, width = 1600, res = 100)
par(mfrow=c(1,2))
    plot(x_mc, y_mc, col = mcmd$color,pch = 16, cex = 100*p_j, xlim = c(-0.5,1.5))
    points(x_mc, y_mc, col = 'black', cex=100*p_j)
    points(x_mcl, y_mcl, cex = 100*p_i)
    segments(x1_edges, y1_edges, x2_edges, y2_edges,lwd = 100*unlist(apply(flow_mat, 1, function(x) x[x>0])))
    text(x_mcl-0.2, y_mcl, 1:length(y_mcl), cex=0.5)
#     y_mcl = unique(y1_edges[order(y2_mcl_edge_mean)])
#     y1_edges = y1_edges[order(y2_mcl_edge_mean)]
    tmp = unique(y1_edges[order(y2_mcl_edge_mean)])
    y1_edges = match(y1_edges[order(y2_mcl_edge_mean)], tmp)
    y2_edges = y2_edges[order(y2_mcl_edge_mean)]
    plot(x_mc, y_mc, col = mcmd$color,pch = 16,cex=100*p_j, xlim = c(-0.5,1.5))
    points(x_mc, y_mc, col = 'black', cex = 100*p_j)
    points(x_mcl, y_mcl, cex = 100*p_i[tmp])
    segments(x1_edges, y1_edges, x2_edges, y2_edges,lwd = 100*unlist(apply(flow_mat, 1, function(x) x[x>0]))[order(y2_mcl_edge_mean)])
    text(x_mcl-0.2, 1:length(y_mcl), tmp, cex=0.5)
# dev.off()
# }

tapply(apply(mc_sol, 1, sum), mcmd$st, sum)

tapply(p_j, mcmd$st, sum)

tapply(apply(mc_sol, 1, sum), mcmd$st, sum)/tapply(p_j, mcmd$st, sum)

# png('./figs/sum_flow_to_mc.png', res = 250, height = 1600, width = 1600)
plot(as.numeric(p_j), as.numeric(apply(flow_mat, 2, sum)), ylab = 'sum flow to metacell', xlab = 'metacell fraction',
    pch = 16, col = mcmd$color, bg='black')
points(as.numeric(p_j), as.numeric(apply(flow_mat, 2, sum)), col='black')
x = seq(0, 0.04, l=50)
lines(x, x, lty='dashed')
# dev.off()

cl_all = do.call('c', lapply(seq_along(1:length(day_mcls)), function(dm, i) {dm[[i]]$cor_km$cluster + (i-1)*K}, dm = day_mcls))

sc_mic_cl_full = data.frame(t(tgs_matrix_tapply(as.matrix(prom_sum@mat), cl_all, mean)))
mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% apply(flow_mat, 2, function(x) x/sum(x))
# mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% flow_mat

mcl_flow_norm = mc_from_mcl_flow/colSums(mc_from_mcl_flow)

saveRDS(object = list('flow_mat' = flow_mat, 'cl_all' = cl_all, 'mc_from_mcl_flow' = mc_from_mcl_flow), file = './data/flow_res.rds')

# flow_res = readRDS('./data/flow_res.rds')

# flow_mat = flow_res$flow_mat
# cl_all = flow_res$cl_all
# mc_from_mcl_flow = flow_res$mc_from_mcl_flow
# colnames(mc_from_mcl_flow) = 1:ncol(mc_from_mcl_flow)

head(mcl_flow_norm)

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
# save_pheatmap_png(p, './figs/cor_mc_atac_rna.png')
options(repr.plot.height = 6, repr.plot.width = 6)


cor_atac_rna_genes = tgs_cor(t(mcl_flow_norm[gb,]), t(log2(mc_rna@e_gc[gb,] + 1e-06)), spearman = T)

rownames(cor_atac_rna_genes) = colnames(cor_atac_rna_genes)

autocor_atac = tgs_cor(t(mcl_flow_norm[gb,]), spearman = T)

rownames(autocor_atac) = colnames(autocor_atac)

head(autocor_atac)

hist(diag(cor_atac_rna_genes))

quantile(diag(cor_atac_rna_genes), seq(0,1,l=21), na.rm = T)

hist(autocor_atac)

# row_inds_ac = matrix(rep(1:ncol(autocor_atac), nrow(autocor_atac)), nrow(autocor_atac), ncol(autocor_atac))
# col_inds_ac = matrix(rep(1:nrow(autocor_atac), ncol(autocor_atac)), nrow(autocor_atac), ncol(autocor_atac), byrow = T)

# ut_ac_a = upper.tri(autocor_atac)
# autocor_atac_tri = autocor_atac[ut_ac_a]
# autocor_atac_row = row_inds_ac[ut_ac_a]
# autocor_atac_col = col_inds_ac[ut_ac_a]

# ac_a_ord = order(autocor_atac_tri,decreasing = T)

# setNames(autocor_atac_tri[head(ac_a_ord, 50)], rownames(autocor_atac)[autocor_atac_row[head(ac_a_ord, 50)]])
# setNames(autocor_atac_tri[head(ac_a_ord, 50)], rownames(autocor_atac)[autocor_atac_col[head(ac_a_ord, 50)]])

# gene_pairs = cbind(rownames(autocor_atac)[autocor_atac_row[head(ac_a_ord, 20)]],
#                     rownames(autocor_atac)[autocor_atac_col[head(ac_a_ord, 20)]])

# gene_pairs

# sapply(1:10, function(gp, i) {
#     plot(log2(mcl_flow_norm[gp[i,1],] + 1e-06), log2(mcl_flow_norm[gp[i,2],] + 1e-06),
#                                                                                      col = mcmd$color, pch = 16, xlab = gp[i,1], ylab = gp[i,2]); 
#                                       points(log2(mcl_flow_norm[gp[i,1],] + 1e-06), log2(mcl_flow_norm[gp[i,2],] + 1e-06), col = 'black');
#                                     }, gp = gene_pairs
#       )

# options(repr.plot.height = 12, repr.plot.width = 14)
# pheatmap::pheatmap(cor_atac_rna_genes, show_colnames = F, show_rownames = F)
#                    annotation_row = annotation_col,
#                    annotation_col = annotation_col, 
#                    annotation_colors = ann_colors

flow_com = apply(flow_mat[,cust_mc_ord_st], 1, function(x) sum(x*(1:length(x)))/sum(x))

flow_com

p = pheatmap::pheatmap(flow_mat[order(flow_com),cust_mc_ord_st], 
                       border_color = 'black',
                       show_colnames = F, 
                       show_rownames = T, 
                       cluster_cols = F, 
                       cluster_rows = F,
                       gaps_col = tail(sapply(cust_st_ord, function(u) min(which(names(cust_mc_ord_st) == u))), -1) - 1,
                    color = colorRampPalette(c('black', 'red', 'orange','yellow', 'white'))(100),
                   annotation_col = annotation_col, 
                   annotation_colors = ann_colors)

save_pheatmap_png(p,'./figs/flow_mat.png')

flow_mat_st = t(tgs_matrix_tapply(flow_mat[order(flow_com),cust_mc_ord_st], names(cust_mc_ord_st), sum))

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
                                       

save_pheatmap_png(p,'./figs/flow_mat_st.png', height = 2200, width= 1200)

hist(diag(cor_atac_rna_genes), 100)

head(diag(cor_atac_rna_genes))

quantile(diag(cor_atac_rna_genes), seq(0,1,l=21), na.rm = T)

genes_to_plot = head(diag(cor_atac_rna_genes)[order(diag(cor_atac_rna_genes), decreasing = T)], 20)
genes_to_plot

genes_to_plot2 = head(diag(cor_atac_rna_genes)[order(diag(cor_atac_rna_genes))], 20)
genes_to_plot2

dir.create('./figs/cor_atac_rna_genes')

names(genes_to_plot)

genes_to_plot

options(repr.plot.height = 6, repr.plot.width = 6)
sapply(seq_along(names(genes_to_plot)), function(g,gt,i) {
#     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
#     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
#     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
#     png(glue::glue('./figs/cor_atac_rna_genes/{i}_{g[[i]]}.png'))
    png(paste0('./figs/cor_atac_rna_genes/',i,'_',g[[i]],'.png'))
    plot(log2(mcl_flow_norm[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), 
         ylab = 'log2 RNA e_gc', xlab = 'log2 ATAC e_gc',
         cex.lab = 2, cex.main = 2,
#          xlim = c(min_val,max_val),
         col = mcmd$color, pch = 16, main = paste(g[[i]],round(gt[[i]],3),sep=' - '));
    points(log2(mcl_flow_norm[g[[i]],] + 1e-06), log2(mc_rna@e_gc[g[[i]],] + 1e-06), col = 'black');
#     axis(2,cex.axis = 3)
#     axis(1,cex.axis = 3)
    dev.off()
    }, g = names(genes_to_plot), gt = genes_to_plot
      )


grep('Myt1l', rownames(mc_rna@e_gc),v=T)

g = 'Myt1l'
png(paste0('./figs/cor_atac_rna_genes/', g, '.png'))
plot(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), 
     ylab = 'log2 RNA e_gc', xlab = 'log2 ATAC e_gc',
              cex.lab = 2, cex.main = 2,
     col = mcmd$color, pch = 16, main = paste(g,round(diag(cor_atac_rna_genes)[[g]],3),sep=' - '))
points(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), col = 'black')
dev.off()

options(repr.plot.height = 6, repr.plot.width = 6)
sapply(names(genes_to_plot2), function(g) {
#     min_val_x = min(log2(mcl_flow_norm[g,] + 1e-06))
#     min_val_y = min(log2(mc_rna@e_gc[g,] + 1e-06))
#     max_val = max(c(log2(mc_rna@e_gc[g,] + 1e-06), log2(mcl_flow_norm[g,] + 1e-06)))
    plot(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), 
#          xlim = c(min_val,max_val),
                                                                                     col = mcmd$color, pch = 16, main = g); 
                                      points(log2(mcl_flow_norm[g,] + 1e-06), log2(mc_rna@e_gc[g,] + 1e-06), col = 'black');
                                    }
      )














