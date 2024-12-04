library(metacell)
library(lpsymphony)
library(pheatmap)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

K = 50
NUM_EDGE_NODE <- 10

day_mcl_path <- file.path(wd, 'output/mcatac/pl_cort_prom_day_mcls_k=50.rds')
flow_res_path <- file.path(wd, "output/mcatac/pl_cort_flow_mat.tsv")

mc_rna = scdb_mc('pl_cort')
prom_sum = scdb_mat('pl_prom_cort')
mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mcell_add_gene_stat(mat_id = 'pl_prom', gstat_id = 'pl_prom', force = F)

gstat = scdb_gstat('pl_prom')

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

png('./figs/atac_feat_select.png', h=1000,w=1000,r=150)
plot(log(gstat$ds_mean), gstat$ds_log_varmean)

points(x[sort(unlist(top_q_inds))], y[sort(unlist(top_q_inds))], col = 'red', pch = 16)

dev.off()

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

days = tail(unique(prom_sum@cell_metadata$day), -1)
days

day_mat = cbind(sapply(days, function(d) apply(mcmd[,grep(d, colnames(mcmd))], 1, sum)))

day_mat = t(apply(day_mat, 1, function(x) x/sum(x)))

mcmd= tibble::tibble(cbind(mcmd, day_mat))
colnames(mcmd)[colnames(mcmd) %in% 1:6] = 13:18
head(mcmd)
prom_mat <- prom_sum@mat
cs <- Matrix::colSums(prom_mat)

prom_mat_ds <- scm_downsamp(prom_mat, min(cs))

make_day_mcl = function(d) {
    # prom_sum = scdb_mat(glue::glue('prom_cort_day_{d}'))
    inds_in = rownames(prom_sum@cell_metadata)[prom_sum@cell_metadata$day == d]
    inds_in = inds_in[inds_in %in% colnames(prom_sum@mat)]
    prom_mat_here = prom_mat_ds[,inds_in]
    cor_mat = tgs_cor(as.matrix(prom_mat_here[feats_filt,]), as.matrix(mc_rna@mc_fp[feats_filt,]), spearman = T)
    cor_km = tglkmeans::TGL_kmeans(df = cor_mat, k = K, seed = SEED)
    sc_mic_cl = data.frame(t(tgs_matrix_tapply(as.matrix(prom_mat_here), cor_km$cluster, mean)))
    cor_mic_cl = tgs_cor(x = as.matrix(sc_mic_cl[feats_filt,]), y = mc_rna@mc_fp[feats_filt,], spearman = T)
    return(list('cor_mic_cl' = cor_mic_cl, 'sc_mic_cl' = sc_mic_cl, 'cor_km' = cor_km, 'cor_mat' = cor_mat, 'prom_mat' = prom_mat_here))
}

day_mcls = parallel::mclapply(days, make_day_mcl, mc.cores = length(days))
names(day_mcls) = days

# saveRDS(day_mcls, './data/pl_cort_prom_day_mcls.rds')
# saveRDS(day_mcls, './data/pl_cort_prom_day_mcls_k=100.rds')
saveRDS(day_mcls, day_mcl_path)
# day_mcls = readRDS('./data/pl_cort_prom_day_mcls.rds')

# day_mcls = readRDS('./data/pl_cort_prom_day_mcls_k=100.rds')
day_mcls = readRDS(day_mcl_path)

color_key = unique(mcmd[,c('cell_type', 'color')])
annotation_col = data.frame(cell_type = mcmd$cell_type)
rownames(annotation_col) = mcmd$metacell
ann_colors = list(cell_type = setNames(color_key$color, color_key$cell_type))

cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$metacell[mcmd$cell_type == st]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))
# cust_mc_ord_st

day_mcl_fig_dir = './figs/pl_cort_day_mcl'
if (!dir.exists(day_mcl_fig_dir)) {dir.create(day_mcl_fig_dir)}

plot_mc_pheatmap = function(cor_km, cor_mat, cor_mic_cl, day) {
                p = pheatmap::pheatmap(cor_mic_cl[,cust_mc_ord_st], cluster_cols = F, cluster_rows = F, silent = T, fontsize = 18, annotation_legend = F,
                                   main = glue::glue('day {day}'), annotation_col = annotation_col, 
                                   annotation_colors = ann_colors,
                        color = colorRampPalette(c('blue','green', 'white','yellow','red'))(1000))
                save_pheatmap_png(p, paste0(day_mcl_fig_dir, '/', day, '.png'), height = 2500, width = 2500, res = 150)
}

lapply(seq_along(day_mcls), function(x,n,i) {
    plot_mc_pheatmap(x[[i]]$cor_km, x[[i]]$cor_mat, x[[i]]$cor_mic_cl, n[[i]])}, 
       x = day_mcls, n = names(day_mcls)
      )

make_mcl_mc_edges = function(d, cor_mic_cl) {
    edge_inds_mc_avail = which(mcmd[,d] > 0)
    
    edge_inds_mc = t(apply(cor_mic_cl, 2, function(x) head(order(x, decreasing = T), NUM_EDGE_NODE)))
    edge_inds_mcl = t(apply(cor_mic_cl, 1, function(x) {
                                                        cor_vec = x[edge_inds_mc_avail];
                                                        inds_in_vec = head(order(cor_vec, decreasing = T), NUM_EDGE_NODE)
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
    edge_mat_all = edge_mat_mc + edge_mat_mcl
    edge_mat_all[edge_mat_all > 1] = 1
    csi <- which(colSums(edge_mat_all) >= NUM_EDGE_NODE)
    rsi <- which(rowSums(edge_mat_all) >= NUM_EDGE_NODE)
    dsc <- scm_downsamp(edge_mat_all[,csi], NUM_EDGE_NODE)
    dsr <- scm_downsamp(t(edge_mat_all[rsi,]), NUM_EDGE_NODE)
    ema_c <- edge_mat_all
    ema_c[,csi] <- as.matrix(dsc)
    ema_r <- edge_mat_all
    ema_r[rsi,] <- t(as.matrix(dsr))
    dsa <- as.matrix(ema_c) + as.matrix(ema_r)
    edge_mat_all <- dsa
    edge_mat_all[edge_mat_all > 1] = 1
    edge_inds = which(edge_mat_all == 1)
    edge_vals = cor_mic_cl[edge_inds]

    row_inds = matrix(rep(1:nrow(cor_mic_cl), ncol(cor_mic_cl)), nrow(cor_mic_cl), ncol(cor_mic_cl))
    col_inds = matrix(rep(1:ncol(cor_mic_cl), nrow(cor_mic_cl)), nrow(cor_mic_cl), ncol(cor_mic_cl), byrow = T)

    edge_inds_row = row_inds[edge_inds]
    edge_inds_col = col_inds[edge_inds]
    print(dim(edge_mat_all))
    print(dim(annotation_col))
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

p_j_vec = as.numeric(unlist(lapply(p_j, function(x) x*c(0.8, 0.2, 0.2, 2))))

cl_sizes = do.call("c", lapply(day_mcls, function(x) x$cor_km$size))

p_i = cl_sizes/sum(cl_sizes)

rhs = c(rep(0, N+num_edges+4*M), p_i, rep(0, N), rep(0, M), sum(p_i), p_j_vec)

dir = c(rep(">=", N+num_edges+4*M), rep("==", N), rep("==", N), rep("==", M), "==", rep("<=", 4*M))

types = rep('C', ncol(lhs))

K2 = 1e+04
K1 = 1e+01

obj = c(rep(0, N), -edge_vals, rep(c(-K2, -K1, K1, K2), M))

sol = lpsymphony_solve_LP(obj, lhs, dir, rhs, types = types, max = FALSE)

sol$objval

source_sol = sol$solution[1:N]
atac_sol = tapply(sol$solution[(N+1):(N + num_edges)], edge_inds_row, list)
mc_sol = matrix(sol$solution[(N + num_edges + 1):length(sol$solution)], nrow = M, ncol = 4, byrow = T)

flow_mat = matrix(0, nrow(cor_mic_cl), ncol(cor_mic_cl))
for (nm in names(atac_sol)) {
    nm = as.numeric(nm)
    flow_mat[as.numeric(nm),edge_inds_col[edge_inds_row == as.numeric(nm)]] = atac_sol[[nm]]
}

colnames(flow_mat) = 1:ncol(cor_mic_cl)
rownames(flow_mat) = 1:nrow(cor_mic_cl)

readr::write_tsv(tibble::rownames_to_column(as.data.frame(flow_mat)), flow_res_path)

# cl_all = do.call('c', lapply(seq_along(1:length(day_mcls)), function(dm, i) {dm[[i]]$cor_km$cluster + (i-1)*K}, dm = day_mcls))

# sc_mic_cl_full = data.frame(t(tgs_matrix_tapply(as.matrix(prom_sum@mat), cl_all, mean)))
# mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% apply(flow_mat, 2, function(x) x/sum(x))
# # mc_from_mcl_flow = as.matrix(sc_mic_cl_full) %*% flow_mat

# mcl_flow_norm = mc_from_mcl_flow/colSums(mc_from_mcl_flow)

# saveRDS(object = list('flow_mat' = flow_mat, 'cl_all' = cl_all, 'mc_from_mcl_flow' = mc_from_mcl_flow), file = './data/pl_flow_res.rds')

pp <- pheatmap(slanter::slanted_reorder(flow_mat), color = colorRampPalette(c('black', 'red', 'orange', 'yellow', 'white'))(100),
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         cluster_cols = F, 
         cluster_rows = F
)

save_pheatmap_png(pp, './figs/atac_cortex_matching_flow_mat.png', h = 1600, w = 1600)

metacells_from_mcl_flow <- function(flow_path, day_mcl_path) {
    day_mcls <- readRDS(day_mcl_path)
    flow_mat <- as.matrix(tibble::column_to_rownames(readr::read_tsv(flow_path), 'rowname'))
    dmcl_vec <- unlist(lapply(seq_along(day_mcls), 
                                function(x,n,i) setNames(paste0(n[[i]], '_', x[[i]]$cor_km$cluster), 
                                                            colnames(x[[i]]$prom_mat)), 
                                                            x = day_mcls, n = names(day_mcls)))
    mcls_nums <- unlist(lapply(names(day_mcls), function(d) paste0(d, "_", sort(unique(day_mcls[[1]]$cor_km$cluster)))))
    dmcl_rename <- match(dmcl_vec, mcls_nums)
    names(dmcl_rename) <- names(dmcl_vec)
    flow_mat_norm <- flow_mat/rowSums(flow_mat)
    flow_mat_cs <- apply(flow_mat_norm, 1, function(x) {
        y <- setNames(x[which(x > 0)], which(x > 0));
        ycs <- cumsum(y)
        names(ycs) <- names(y)
        return(ycs)
    })
    dmcl_rand_nums <- runif(length(dmcl_rename))
    cell_to_metacell <- as.numeric(unlist(purrr::map2(.x = dmcl_rand_nums, .y = dmcl_rename, .f = function(.x,.y) {
        names(flow_mat_cs[[.y]])[min(which(flow_mat_cs[[.y]] >= .x))]
    })))
    names(cell_to_metacell) <- names(dmcl_vec)
    return(tibble::enframe(cell_to_metacell, name = 'cell', value = 'metacell'))
}

c2mc <- metacells_from_mcl_flow(flow_path = flow_res_path, day_mcl_path = day_mcl_path)

readr::write_tsv(c2mc, "./output/mcatac/c2mc.tsv")