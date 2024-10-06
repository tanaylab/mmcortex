# devtools::load_all('~/src/metacell')
library(metacell)
devtools::load_all('~/src/metacell.flow')
library(pheatmap)
library(ComplexHeatmap)
library(Matrix)
library(matrixStats)
library(princurve)
library(umap)

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

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         
cust_st_ord2 = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
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

clrmp <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue1', 'blue4', 'purple3'))(1000)

clrmp_abs <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(1000)
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)

ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

legc = log2(1e-05 + mc@e_gc)

st_legc <- as.data.frame(t(tgs_matrix_tapply(legc, mcmd$cell_type, mean)))

feats <- scdb_gset('pl_filt_lat')
feats <- names(feats@gene_set)
feats <- feats[feats %in% rownames(mc@e_gc)]
feats_all <- scdb_gset('pl')

nsc_mcs <- which(mcmd$cell_type == 'NSC')
ipc_cyc_mcs <- which(mcmd$cell_type == 'IPC_cyc')
ipc_mcs <- which(mcmd$cell_type == 'IPC')

nsc_legc_gq <- t(apply(legc[,nsc_mcs], 1, quantile, c(0,0.01,0.05,0.1,0.9, 0.95,0.99,1)))
g_cor_cc <- lapply(c('Pcna', 'Mki67', 'Top2a', 'Mcm5'), function(g) {cor_legc_g <- tgs_cor(t(legc[,nsc_mcs]), as.matrix(legc[g,nsc_mcs]), spearman = F);
                                                        return(rownames(cor_legc_g)[cor_legc_g >= 0.5 & apply(nsc_legc_gq[,c(3,6)], 1, diff) >= 2])})

g_anti_cor_cc <- lapply(c('Pcna', 'Mki67', 'Top2a', 'Mcm5'), function(g) {cor_legc_g <- tgs_cor(t(legc[,nsc_mcs]), as.matrix(legc[g,nsc_mcs]), spearman = F);
                                                        return(rownames(cor_legc_g)[cor_legc_g <= -0.5 & apply(nsc_legc_gq[,c(3,6)], 1, diff) >= 2])})

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- intersect(c('Bcl11b', tfs$Symbol), rownames(legc))

get_genes_specific_to_mcs <- function(legc, mc_pos = NULL, mc_neg = NULL, cl_vec = NULL) {
    if (!is.null(mc_pos) && is.null(mc_neg)) {
        cl_vec <- ifelse(1:ncol(legc) %in% mc_pos, 1, 0)
    } else if (!is.null(mc_pos) && !is.null(mc_neg)) {
        # if (!(length(intersect(mc_pos, mc_neg)) == 0)) {
        #     stop('mc_pos and mc_neg intersect')
        # }
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

all_ct_genes <- get_genes_specific_to_mcs(legc, cl_vec = mcmd$cell_type)
names(all_ct_genes) <- sort(unique(mcmd$cell_type))

nsc_dyn_genes <- rownames(nsc_legc_gq)[which(apply(nsc_legc_gq[,c(3,6)], 1, diff)>= 2 & nsc_legc_gq[,6] >= -14)]
plot(nsc_legc_gq[,3], nsc_legc_gq[,6])
points(nsc_legc_gq[nsc_dyn_genes,3], nsc_legc_gq[nsc_dyn_genes,6], col ='red', pch = 16)

calc_cor_gene_legc_w_md_in_st <- function(legc, cell_type, mcmd) {
    mcs_st <- which(mcmd$cell_type == cell_type)
    cor_genes_md <- tgs_cor(t(legc[,mcs_st]), as.matrix(mcmd$mean_day[mcs_st]), spearman=T)
    names(cor_genes_md) <- rownames(legc)
    return(cor_genes_md)
}

cor_gene_legc_w_md_in_st <- do.call('cbind', lapply(unique(mcmd$cell_type), function(ct) calc_cor_gene_legc_w_md_in_st(legc, ct, mcmd)))
colnames(cor_gene_legc_w_md_in_st) <- unique(mcmd$cell_type)

cor_nsc_legc_dyn_genes <- tgs_cor(t(legc[nsc_dyn_genes,nsc_mcs]), spearman = T)
hc_cor_nsc <- hclust(dist(cor_nsc_legc_dyn_genes), method = 'ward.D2')
ct_hc_cor_nsc <- cutree(hc_cor_nsc, k = 8)

save(ct_hc_cor_nsc, file = './output/metacell_model/nsc_gene_modules.rda')

## Generate mc2d not correlated to cell cycle
## Start
feats_new <- feats[!(feats %in% names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3,6,7)]))]
gs_new <- gset_new_gset(sets = setNames(rep(1, length(feats_new)), feats_new), desc = 'pl_filt_lat and removed all cell cycle-correlated genes')

new_gs_id <- 'pl_filt_lat_cor_cc'
mc_id <- 'pl_cort'
new_mgraph_id <- 'pl_cort_not_cor_cc'

scdb_add_gset(id = new_gs_id, gset = gs_new)

gs <- scdb_gset(new_gs_id)
mc <- scdb_mc(mc_id)

feats_cc <- feats[feats %in% names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3,6,7)])]
cc_gs <- gset_new_gset(sets = setNames(rep(1, length(feats_cc)), feats_cc), desc = 'cell-cycle-correlated genes in NSC')

scdb_add_gset(id = 'pl_cort_nsc_cor_cc', gset = cc_gs)

mcell_add_cgraph_from_mat_bknn(mat_id='pl_cort',
                gset_id = 'pl_cort_nsc_cor_cc',
                graph_id='pl_cort_cor_cc',
                K=20,
                dsamp=T)


mcell_mgraph_logistic(mgraph_id=new_mgraph_id, mc_id=mc_id,feats_gset=new_gs_id)

mc = scdb_mc(mc_id)
mg = scdb_mgraph(new_mgraph_id)
mgraph = mg@mgraph

feat_genes = scdb_gset(new_gs_id)
feat_genes = names(feat_genes@gene_set)

uconf <- umap::umap.defaults

uconf$n_neighbors=7
uconf$min_dist = 0.5
uconf$spread = 1
symmetrize = F
umap_mgraph = F

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)

cgraph_id <- 'pl_cort'

xy = mc2d_comp_cell_coord(mc_id = mc_id,graph_id = cgraph_id, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_init('./scdb/', force_reinit = T)

mc2d_id <- new_mgraph_id

mc2d <- tgMC2D(mc_id, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph)

scdb_add_mc2d(mc2d_id, mc2d)
## End



load('./output/metacell_model/pcurve_nsc_vs_mat_neuro.rda')
ro <- rev(pcu$ord)

ra <- as.data.frame(tibble::column_to_rownames(tibble::enframe(ct_hc_cor_nsc, name = 'gene', value = 'cluster'), 'gene'))
ac <- list(cluster = setNames(chameleon::distinct_colors(8)$name, 1:8))
p_cor_dyn_genes_nsc <- pheatmap::pheatmap(cor_nsc_legc_dyn_genes[hc_cor_nsc$order,hc_cor_nsc$order], show_rownames = F, show_colnames = F, 
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-1,1,l=100),
                   cluster_cols = F, cluster_rows = F, 
                   annotation_col = ra, annotation_row = ra, annotation_colors = ac, silent = TRUE)
save_pheatmap_png(p_cor_dyn_genes_nsc, './output/metacell_model/figs/nsc_gene_module_analysis/cor_dyn_genes_nsc.png', h = 1600, w = 1600)

ipc_cyc_vs_ipc_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'IPC_cyc'), mc_neg = which(mcmd$cell_type == 'IPC'))

ipc_genes_vs_glia <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('IPC_cyc', 'IPC')), 
                                                mc_neg = which(mcmd$cell_type %in% c('Astrocytes', 'Oligodendrocytes')))
ipc_genes_vs_nsc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('IPC_cyc', 'IPC')), 
                                                mc_neg = which(mcmd$cell_type %in% c('NSC')))

ipc_module <- intersect(names(ipc_genes_vs_glia[ipc_genes_vs_glia >= 2]),
                    names(ipc_genes_vs_nsc[ipc_genes_vs_nsc >= 2]))


pheatmap::pheatmap(legc[ipc_module,cust_mc_ord_st2]- rowMeans(legc[ipc_module,]), 
                  # fontsize_row = 6, 
                   show_rownames = T, show_colnames = F, treeheight_col = 0, annotation_legend = F,
                   treeheight_row = 0,
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-3,3,l=100),
                   cluster_cols = F, cluster_rows = T, 
                   annotation_col = col_annot, annotation_colors = ann_colors)

astro_genes_vs_oligo <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('Oligodendrocytes')))
astro_genes_vs_nsc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('NSC')))
astro_genes_vs_ipc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc')))

astro_module <- intersect(names(astro_genes_vs_oligo[astro_genes_vs_oligo >= 2]),
                      intersect(names(astro_genes_vs_nsc[astro_genes_vs_nsc >= 2]),
                           names(astro_genes_vs_ipc[astro_genes_vs_ipc >= 2])))

write(astro_module, './output/metacell_model/nsc_gene_modules/astro_module.txt')
write(ipc_module, './output/metacell_model/nsc_gene_modules/ipc_module.txt')
write(stem_module, './output/metacell_model/nsc_gene_modules/stem_module.txt')

stem_genes_nsc_fc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in%  c('NSC')), mc_neg = which(mcmd$cell_type %in%  c('Astrocytes', 'IPC')))
stem_genes_ipc_fc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in%  c('IPC_cyc')), mc_neg = which(mcmd$cell_type %in%  c('Astrocytes', 'IPC')))

cond1 <- st_legc[,'NSC'] - st_legc[,'Astrocytes'] >= 1
cond2 <- st_legc[,'NSC'] - st_legc[,'IPC'] >= 1
cond3 <- st_legc[,'IPC_cyc'] - st_legc[,'Astrocytes'] >= 1
cond4 <- st_legc[,'IPC_cyc'] - st_legc[,'IPC'] >= 1
stem_genes <- rownames(st_legc)[(cond1 & cond2) | (cond3 & cond4)]
sort(stem_genes)

cor_stem_cc_nsc <- tgs_cor(t(legc[stem_genes,which(mcmd$cell_type %in%  c('NSC'))]), 
                       t(legc[c('Top2a', 'Mki67', 'Mcm4', 'Pcna'), which(mcmd$cell_type %in%  c('NSC'))]))
cor_stem_cc_ipc_cyc <- tgs_cor(t(legc[stem_genes,which(mcmd$cell_type %in%  c('IPC_cyc'))]), 
                       t(legc[c('Top2a', 'Mki67', 'Mcm4', 'Pcna'), which(mcmd$cell_type %in%  c('IPC_cyc'))]))

cor_stem_cc_nsc_filt <- cor_stem_cc_nsc[rowMaxs(cor_stem_cc_nsc) < .5,]
cor_stem_cc_ipc_cyc_filt <- cor_stem_cc_ipc_cyc[rowMaxs(cor_stem_cc_ipc_cyc) < .5,]

stem_module <- union(rownames(cor_stem_cc_nsc_filt), rownames(cor_stem_cc_ipc_cyc_filt))

asc_module <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 4])

cond1 <- rowMaxs(subset(as.matrix(st_legc[asc_module,]), select = -c(Astrocytes,Oligodendrocytes,NSC))) - st_legc[asc_module,'NSC'] >= 1

asc_neuro_module <- asc_module[cond1]

p_dyn_genes_nsc <- pheatmap::pheatmap(legc[names(ct_hc_cor_nsc[hc_cor_nsc$order]),nsc_mcs] - rowMeans(legc[names(ct_hc_cor_nsc[hc_cor_nsc$order]),nsc_mcs]),
                                      fontsize_row = 6, show_rownames = T, show_colnames = F, treeheight_col = 0, annotation_legend = F,
                   col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-3,3,l=100),
                   cluster_cols = T, cluster_rows = F, 
                   annotation_col = col_annot, annotation_row = ra, annotation_colors = c(ac, ann_colors))
save_pheatmap_png(p_dyn_genes_nsc, './output/metacell_model/figs/nsc_gene_module_analysis/all_nsc_dynamic_genes.png', h = 5200, w = 2600)

legc_nsc_dyn_genes_avg_cl <- tgs_matrix_tapply(t(legc[nsc_dyn_genes,nsc_mcs]), ct_hc_cor_nsc, mean)

cl1_mc_clvls <- clrmp[1+999*(legc_nsc_dyn_genes_avg_cl[1,] - min(legc_nsc_dyn_genes_avg_cl[1,]))/(max(legc_nsc_dyn_genes_avg_cl[1,]) - min(legc_nsc_dyn_genes_avg_cl[1,]))]
cl6_mc_clvls <- clrmp[1+999*(legc_nsc_dyn_genes_avg_cl[6,] - min(legc_nsc_dyn_genes_avg_cl[6,]))/(max(legc_nsc_dyn_genes_avg_cl[6,]) - min(legc_nsc_dyn_genes_avg_cl[6,]))]
cl7_mc_clvls <- clrmp[1+999*(legc_nsc_dyn_genes_avg_cl[7,] - min(legc_nsc_dyn_genes_avg_cl[7,]))/(max(legc_nsc_dyn_genes_avg_cl[7,]) - min(legc_nsc_dyn_genes_avg_cl[7,]))]

png('./output/metacell_model/figs/nsc_gene_module_analysis/cc_module_mc_scatter.png', w = 1200, h = 400)
par(mfrow = c(1,3), cex.lab = 3, mar = c(5,5,4,1), cex.axis = 2, cex.main = 3)
plot(legc_nsc_dyn_genes_avg_cl[1,], legc_nsc_dyn_genes_avg_cl[6,], xlab = 'Cluster 1 MC', ylab = 'Cluster 6 MC', main = 'Color = cluster 7', pch = 16, cex = 2.5, col = cl7_mc_clvls)
points(legc_nsc_dyn_genes_avg_cl[1,], legc_nsc_dyn_genes_avg_cl[6,], cex = 2.5)
plot(legc_nsc_dyn_genes_avg_cl[1,], legc_nsc_dyn_genes_avg_cl[7,], xlab = 'Cluster 1 MC', ylab = 'Cluster 7 MC', main = 'Color = cluster 6', pch = 16, cex = 2.5, col = cl6_mc_clvls)
points(legc_nsc_dyn_genes_avg_cl[1,], legc_nsc_dyn_genes_avg_cl[7,], cex = 2.5)
plot(legc_nsc_dyn_genes_avg_cl[6,], legc_nsc_dyn_genes_avg_cl[7,], xlab = 'Cluster 6 MC', ylab = 'Cluster 7 MC',main = 'Color = cluster 1',  pch = 16, cex = 2.5, col = cl1_mc_clvls)
points(legc_nsc_dyn_genes_avg_cl[6,], legc_nsc_dyn_genes_avg_cl[7,],cex = 2.5)
dev.off()



md_clvls <- clrmp[1+round(999*(mcmd$mean_day-13)/(18-13))]

mcell_mc2d_plot(mc2d_id = 'pl_cort_not_cor_cc', colors = md_clvls, fig_fn = './output/metacell_model/figs/pl_cort_mc2d_col_by_mean_day.png')

mc2d <- scdb_mc2d('pl_cort_not_cor_cc')

sc_t <- mat@cell_metadata[names(mc2d@sc_x), 't']

mat_ds <- scm_downsamp(mat@mat, 3000)

nsc_sc <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'NSC')]), colnames(mat_ds))

nsc_genes <- sort(rownames(all_ct_genes[['NSC']])[all_ct_genes[['NSC']][,'Astrocytes'] <= -1 & all_ct_genes[['NSC']][,'IPC_cyc'] <= -1])

cl1_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 1])
cl6_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 6])
cl7_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 7])
cl3_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 3])
cl13_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3)])
cl1_sc <- colSums(mat_ds[cl1_genes,])
cl13_sc <- colSums(mat_ds[cl13_genes,])
cl3_sc <- colSums(mat_ds[cl3_genes,])
cl6_sc <- colSums(mat_ds[cl6_genes,])
cl7_sc <- colSums(mat_ds[cl7_genes,])
cc_sc <- colSums(mat_ds[c(cl1_genes, cl3_genes, cl6_genes, cl7_genes),])
ipc_sc <- colSums(mat_ds[ipc_module,])
astro_sc <- colSums(mat_ds[astro_module,])
stem_sc <- colSums(mat_ds[stem_module,])
nsc_genes_sc <- colSums(mat_ds[nsc_genes,])


png('./output/metacell_model/nsc_gene_modules/figs/pairwise_cell_cycle_clusters_in_scs.png', h = 350, w = 1050)
par(mfrow = c(1,3), cex.lab = 2, mar = c(5,5,1,1))
plot(cl1_sc, cl6_sc, col = color_key$color[match(mcmd$cell_type[mc@mc[names(stem_sc)]], color_key$cell_type)], pch = 16, cex = 0.57, ylab = 'Cluster 6 + 3', xlab = 'Cluster 1')
plot(cl1_sc, cl7_sc, col =color_key$color[match(mcmd$cell_type[mc@mc[names(stem_sc)]], color_key$cell_type)], pch = 16, cex = 0.57, ylab = 'Cluster 7', xlab = 'Cluster 1')
plot(cl7_sc, cl6_sc, col =color_key$color[match(mcmd$cell_type[mc@mc[names(stem_sc)]], color_key$cell_type)], pch = 16, cex = 0.57, ylab = 'Cluster 6 + 3', xlab = 'Cluster 7')
dev.off()

cc_data_mat <- cbind(cl1_sc - mean(cl1_sc), cl3_sc - mean(cl3_sc), cl6_sc - mean(cl6_sc), cl7_sc - mean(cl7_sc))
pca_cc <- prcomp(x = cc_data_mat)
cc_pca_mat <- pca_cc$x
theta <- 140*pi/180
cc_mat_rot <- cc_pca_mat %*% rbind(c(cos(theta), -sin(theta), 0, 0), c(sin(theta), cos(theta), 0, 0), c(0, 0, 1,  0), c(0, 0, 0, 1))

cycling_scs <- which(cc_sc >= quantile(cc_sc, 0.8))

phase_rad <- -atan2(cc_mat_rot[,2] - mean(cc_mat_rot[names(cycling_scs),2]), cc_mat_rot[,1] - mean(cc_mat_rot[names(cycling_scs),1]))

phase_clvls <- setNames(clrmp[1+round(999*(phase_rad[rownames(cc_mat_rot)] - min(phase_rad))/(max(phase_rad) - min(phase_rad)))],
                        rownames(cc_mat_rot))

clrmp_mod <- colorRampPalette(c('white', 'orange', 'purple', 'black'))(1000)

phase_clvls_mod <- setNames(clrmp_mod[1+round(999*(phase_rad[rownames(cc_mat_rot)] - min(phase_rad))/(max(phase_rad) - min(phase_rad)))],
                        rownames(cc_mat_rot))

nsc_ipc_cyc_cells_in_cc_mat <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type %in% c('NSC', 'IPC_cyc'))]), rownames(cc_mat_rot))

phase_pcu <- principal_curve(x = cc_mat_rot[nsc_ipc_cyc_cells_in_cc_mat[order(phase_rad[match(nsc_ipc_cyc_cells_in_cc_mat, names(phase_rad))])],1:2], 
                                    smoother = "periodic_lowess", stretch = 1, approx_points = 100, thresh = 1e-5)
phase_pcu_rm <- zoo::rollmean(x = phase_pcu$s[,1:2], k = 50)
phase_pcu_rm_ss <- phase_pcu_rm[sort(sample(x = nrow(phase_pcu_rm), size = round(nrow(phase_pcu_rm)/10), replace = F)),]

dist_phase_pc_to_pcu <- as.matrix(tgstat::tgs_dist(rbind(cc_mat_rot[,1:2], phase_pcu_rm_ss[,1:2])))

rownames(dist_phase_pc_to_pcu) <- c(rownames(dist_phase_pc_to_pcu)[1:nrow(cc_mat_rot)], paste0('p', 1:nrow(phase_pcu_rm_ss)))
colnames(dist_phase_pc_to_pcu) <- c(colnames(dist_phase_pc_to_pcu)[1:nrow(cc_mat_rot)], paste0('p', 1:nrow(phase_pcu_rm_ss)))

knn_phase_pc_to_pcu <- tgstat::tgs_knn(1/dist_phase_pc_to_pcu[(nrow(cc_mat_rot)+1):ncol(dist_phase_pc_to_pcu),1:nrow(cc_mat_rot)], k = 1)
rownames(dist_phase_pc_to_pcu)[(nrow(cc_mat_rot)+1):ncol(dist_phase_pc_to_pcu)]
knn_phase_pc_to_pcu[,2] <- paste0('p', 1:nrow(phase_pcu_rm_ss))[match(knn_phase_pc_to_pcu[,2], rownames(dist_phase_pc_to_pcu))]

phase <- setNames(as.numeric(gsub('^p', '', knn_phase_pc_to_pcu[,2]))[match(rownames(cc_mat_rot), as.character(knn_phase_pc_to_pcu[,1]))],
                  rownames(cc_mat_rot))

phase_rev <- max(phase) - phase + 1

NUM_PARTITION <- 13
phase_qs <- seq(1-1e-2,max(phase),l=NUM_PARTITION)

phase_cut <- setNames(cut(phase_rev, phase_qs, labels = paste('bin', 1:(NUM_PARTITION - 1))), names(phase_rev))

rb12 <- rainbow(12)
phase_cut_color <- setNames(rainbow(12)[phase_cut], names(phase_cut))

png('./output/metacell_model/nsc_gene_modules/figs/pc2_vs_pc1_color_phase_bin.png', h = 600, w = 650, res = 100)
par(cex.axis = 1.5, cex.lab = 2, mar = c(5,5,1,1))
# plot(cc_pca_mat[,1], cc_pca_mat[,2], main = 'PCA no rotation', col = phase_clvls,  pch = 16, cex = 0.37,ylab = 'PC2', xlab = 'PC1')
plot(cc_mat_rot[,1], cc_mat_rot[,2], main = '', 
     col = phase_cut_color,
        pch = 16, cex = 0.67,ylab = 'PC2', xlab = 'PC1')
uuu <- sapply(sort(unique(phase_cut)), function(u) {
    cells <- intersect(names(phase_cut)[phase_cut == u], rownames(cc_mat_rot))
    com_u <- colMeans(cc_mat_rot[cells,1:2])
    text(com_u[[1]], com_u[[2]], labels = u, cex = 1.5)
})

lines(phase_pcu_rm_ss[,1], phase_pcu_rm_ss[,2])
legend('topleft', legend = sort(unique(phase_cut)), col = unique(phase_cut_color[order(phase_cut)]), pch = 16, cex = 1)
dev.off()

plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,1001,round(1001/6)))
{
  if (!is.null(fig_fn)) {
      png(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,1000), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 107, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 107, length(cols)+1, col=NA, border = 'black')
  text(117, show_vals_ind,cex = 1.5, labels=round(vals[show_vals_ind], 3), pos=4)
  if (!is.null(fig_fn)) {
    dev.off()
  }
}

plot_color_bar(seq(-180,180,l=1001), clrmp, fig_fn = './output/metacell_model/figs/clrmp_phase.png') 

phase <- setNames(as.numeric(gsub('^p', '', knn_phase_pc_to_pcu[,2]))[match(rownames(cc_mat_rot), as.character(knn_phase_pc_to_pcu[,1]))],
                  rownames(cc_mat_rot))

cc_genes_select <- c('Pcna', 'Mcm7', 'Mki67', 'Top2a')

mat_ds_cc_genes_select_ord_phase <- mat_ds[cc_genes_select, nsc_sc[rev(order(phase[nsc_sc]))]]

K <- 100
mat_ds_cc_genes_select_ord_phase_rm <- t(apply(mat_ds_cc_genes_select_ord_phase, 1, zoo::rollmean, k = K, na.pad = T))

pba <- setNames(as.numeric(phase_cut), names(phase_cut))

pba[pba %in% c(1,4,12)] <- '1_G1'  
pba[pba %in% 5:7] <- '2_S'  
pba[pba %in% 8:9] <- '3_G2'  
pba[pba %in% 10:11] <- '4_M'  
pba[pba %in% 2:3] <- '5_G0'

png('./output/metacell_model/nsc_gene_modules/figs/cc_genes_ord_phase_rollmean_clean_new.png', h = 1400, w = 1800, res = 100)
par(mfrow = c(2,2), mar = c(6,7,4,1), cex.lab = 3, cex.axis = 3, cex.main = 5)

bin_borders <- setNames(phase_qs[c(1,2,4,5,8,10,12, 13)], c('G1', 'G1_0', 'G1', 'S', 'G2', 'M', 'G1'))
heights_text <- c(1.5, 2.2, 2.2, 4.2)
ttt <- sapply(1:nrow(mat_ds_cc_genes_select_ord_phase), function(i) {
    plot(sort(phase[nsc_sc]), mat_ds_cc_genes_select_ord_phase[i,], 
         # col = col_key['NSC'], 
         col = 'white', 
         cex = 1, pch = 16, ylim = quantile(mat_ds_cc_genes_select_ord_phase[i,], c(0.05, 0.95)),
        xlab = '', ylab = ''
        )
    vvv <- sapply(head(seq_along(bin_borders), -1), function(j) {
        bj <- bin_borders[[j]]
        lines(rep(bj, 2), c(-1,10), lty = 2, lwd = 2)
        nmj <- names(bin_borders)[[j]]
    })
    lines(sort(phase[nsc_sc]), mat_ds_cc_genes_select_ord_phase_rm[i,], col = 'black', lwd = 3)
    if (i == 1) {legend('topleft', legend = glue::glue('Rollmean k = {K}'), lwd = 3, col = 'black', cex = 2)}
    title(xlab = 'phase', line = 4)
    title(ylab = 'Downsampled UMIs', line = 4)
})
dev.off()

ct_by_phase_bin <- as.matrix(t(table(phase_cut, mcmd$cell_type[mc@mc[names(phase_cut)]])))
nsc_sc_by_day <- lapply(13:18, function(di) nsc_sc[nsc_sc %in% rownames(mat@cell_metadata[mat@cell_metadata$day == paste0('E',di),])])
names(nsc_sc_by_day) <- 13:18

nsc_sc_by_day_vec <- setNames(unlist(sapply(names(nsc_sc_by_day), function(x) rep(x, length(nsc_sc_by_day[[x]])))), do.call('c', nsc_sc_by_day))
s_genes_sum <- colSums(mat_ds[s_genes,])
m_genes_sum <- colSums(mat_ds[m_genes,])
sm_genes_sum <- colSums(mat_ds[union(m_genes, s_genes),])

png('./output/metacell_model/nsc_gene_modules/figs/s_m_genes_per_bin_nsc.png', h = 500, w = 1500)
par(mfrow = c(1,3), cex.main = 2, cex.axis = 2, cex.lab = 2, mar = c(6,6,5,1), las = 2)
boxplot(s_genes_sum[names(nsc_sc_by_day_vec)]/length(s_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], main = 'UMIs per S gene per single NSC per bin\nn_{S genes} = 11', ylab = '', xlab = '')
title(ylab = 'UMIs per gene', line = 4)
boxplot(m_genes_sum[names(nsc_sc_by_day_vec)]/length(m_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], main = 'UMIs per M gene per single NSC per bin\nn_{M genes} = 26', ylab = 'UMIs per gene', xlab = '')

boxplot(s_genes_sum[names(nsc_sc_by_day_vec)]/length(s_genes) + m_genes_sum[names(nsc_sc_by_day_vec)]/length(m_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], 
        main = 'UMIs per S+M gene per cell per bin', ylab = 'UMIs per gene', xlab = '')
dev.off()

png('./output/metacell_model/nsc_gene_modules/figs/sum_of_s_m_genes_per_ct.png', h = 400, w = 800)
par(mfrow = c(1,2), las = 2, mar = c(10,5,2,1), cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(s_genes_sum ~ mcmd$cell_type[mc@mc[names(s_genes_sum)]], xlab = '', ylab = '', col = col_key[sort(unique(mcmd$cell_type[mc@mc[names(s_genes_sum)]]))], main = 'S genes')
title(ylab = 'Sum of downsampled UMIs', line = 3)
boxplot(m_genes_sum ~ mcmd$cell_type[mc@mc[names(m_genes_sum)]], xlab = '', ylab = '', col = col_key[sort(unique(mcmd$cell_type[mc@mc[names(s_genes_sum)]]))], main = 'M genes')
title(ylab = 'Sum of downsampled UMIs', line = 3)
dev.off()

mat_ds_avg_phase <- t(tgs_matrix_tapply(as.matrix(mat_ds[,intersect(names(phase_cut), nsc_sc)]), phase_cut[intersect(names(phase_cut), nsc_sc)], mean))
mat_ds_avg_phase_per_day <- lapply(nsc_sc_by_day, function(scdi) t(tgs_matrix_tapply(as.matrix(mat_ds[,intersect(names(phase_cut), scdi)]), 
                            phase_cut[intersect(names(phase_cut), scdi)], mean)))
gin <- which(rowMaxs(mat_ds_avg_phase) - rowMins(mat_ds_avg_phase) >= 1)

legc_avg_ct <- t(tgs_matrix_tapply(legc, mcmd$cell_type, mean))

annot_ct_exp <- as.data.frame(apply(legc_avg_ct[,c('Astrocytes', 'IPC')], 1, diff))
colnames(annot_ct_exp) <- 'IPC_vs_Astro_bias'

ipc_sc_names <- names(mc@mc)[mc@mc %in% which(mcmd$cell_type %in% c('IPC'))]
ipc_cyc_sc_names <- names(mc@mc)[mc@mc %in% which(mcmd$cell_type %in% c('IPC_cyc'))]

astro_sc_names <- names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'Astrocytes')])

sc_data_df <- dplyr::select(mat@cell_metadata[colnames(mat_ds),], day)
sc_data_df$cell_type <- mcmd$cell_type[mc@mc[rownames(sc_data_df)]]
sc_data_df$ipc <- ipc_sc[rownames(sc_data_df)]
sc_data_df$astro <- astro_sc[rownames(sc_data_df)]
sc_data_df$nsc <- nsc_genes_sc[rownames(sc_data_df)]
sc_data_df$stem <- stem_sc[rownames(sc_data_df)]
sc_data_df$cc <- cc_sc[rownames(sc_data_df)]
sc_data_df$phase <- phase[rownames(sc_data_df)]
sc_data_df$phase_rev <- phase_rev[rownames(sc_data_df)]
sc_data_df$phase_cut <- as.numeric(phase_cut[rownames(sc_data_df)])

sc_data_df <- dplyr::mutate(sc_data_df, pba = pba)

mg_bon_marks <- as.data.frame(t(sapply(apply(readr::read_csv('./BonevCollab//marker_genes.tsv'), 1, 
                                             stringr::str_split, ' '), function(x) c(x[[1]][[1]], x[[1]][[length(x[[1]])]]))))

colnames(mg_bon_marks) <- c('cell_type', 'marks')

mbm_lst <- lapply(1:nrow(mg_bon_marks), function(n) stringr::str_split(mg_bon_marks$marks[[n]], ',')[[1]])
names(mbm_lst) <- mg_bon_marks$cell_type


sc_data_df$nsc <- colSums(mat_ds[mbm_lst[['NSC']],rownames(sc_data_df)])
sc_data_df$color <- col_key[sc_data_df$cell_type]

png('./output/metacell_model/figs/ipc_module_vs_nsc_markers_in_nsc_and_ipc_cycs.png', h = 600, w = 800)
par(las = 2, mar = c(13,6,3,1), mfrow = c(1,2), cex.axis = 2, cex.main = 2, cex.lab = 2)
boxplot(nsc ~ ., data = dplyr::select(dplyr::filter(sc_data_df, cell_type %in% c('NSC', 'IPC_cyc')), nsc, cell_type, pba),  col = rep(col_key[c('IPC_cyc', 'NSC')],5), xlab = '' ,main = 'NSC marker genes', ylab = '')
title(ylab = 'UMIs', line = 4)
boxplot(ipc ~ ., data = dplyr::select(dplyr::filter(sc_data_df, cell_type %in% c('NSC', 'IPC_cyc')), ipc, cell_type, pba), col = rep(col_key[c('IPC_cyc', 'NSC')],5), xlab = '',main = 'IPC module', ylab = '')
title(ylab = 'UMIs', line = 4)
dev.off()

cor_stem_cc_nsc_filt <- cor_stem_cc_nsc[rowMaxs(cor_stem_cc_nsc) < .5,]
cor_stem_cc_ipc_cyc_filt <- cor_stem_cc_ipc_cyc[rowMaxs(cor_stem_cc_ipc_cyc) < .5,]

days <- unique(mat@cell_metadata[nsc_sc,'day'])

mari <- c(5,5,2,2)
BOXWEX <- 0.6
LINE <- 5
png('./output/metacell_model/nsc_gene_modules/figs/gene_modules_in_scs_by_day_boxplots.png', h = 1100, w = 1200)
par(mfcol = c(3,3), mar = mari, cex.main = 2, cex.lab = 2, cex.axis = 2, las = 2)
mari[[2]] <- 7
mari[[3]] <- 5
par(mar = mari)
boxplot(ipc_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX, xaxt = 'n',ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', main = 'NSCs', xlab = '', col = col_key[['NSC']])
title(ylab = "IPC module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX,xaxt = 'n',ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), ylab = '', xlab = '', col = col_key[['NSC']])
title(ylab = "Astro module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)
boxplot(stem_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX,xaxt = 'n',ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', xlab = '', col = col_key[['NSC']])
title(ylab = "Stem module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)

mari[[2]] <- 3
mari[[3]] <- 5
par(mar = mari)
boxplot(ipc_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5), ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', main = 'Astrocytes', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5),ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), ylab = '', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))
boxplot(stem_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5),ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))

mari[[3]] <- 5
mari[[2]] <- 3
par(mar = mari)
boxplot(ipc_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), main = 'IPCs',ylab = '', xlab = '', col = col_key[['IPC']])
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), xlab = '', ylab = '',col = col_key[['IPC']])
boxplot(stem_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), xlab = '',ylab = '', col = col_key[['IPC']])
dev.off()

nsc_sc_by_day <- lapply(13:18, function(di) nsc_sc[nsc_sc %in% rownames(mat@cell_metadata[mat@cell_metadata$day == paste0('E',di),])])
names(nsc_sc_by_day) <- 13:18

ipc_sc_by_day <- lapply(13:18, function(di) names(mc@mc)[mc@mc %in% which(mcmd$cell_type %in% c('IPC')) & 
                    names(mc@mc) %in% rownames(mat@cell_metadata[mat@cell_metadata$day == paste0('E',di),]) & 
                    names(mc@mc) %in% colnames(mat_ds)])
names(ipc_sc_by_day) <- 13:18

ipc_cyc_sc_by_day <- lapply(13:18, function(di) names(mc@mc)[mc@mc %in% which(mcmd$cell_type %in% c('IPC_cyc')) & 
                            names(mc@mc) %in% rownames(mat@cell_metadata[mat@cell_metadata$day == paste0('E',di),]) & 
                            names(mc@mc) %in% colnames(mat_ds)])
names(ipc_cyc_sc_by_day) <- 13:18

sc_by_day <- lapply(13:18, function(di) {y <- rownames(mat@cell_metadata[mat@cell_metadata$day == paste0('E',di),]); return(y[y %in% names(mc@mc)])})

names(sc_by_day) <- 13:18

# png('./output//metacell_model/nsc_gene_modules/figs/astro_vs_ipc_in_nsc_by_day_scatter.png', h = 400, w = 800)
mari <- c(5,5,3,1)
par(mfrow = c(2,3), mar = mari, cex.main = 3, cex.lab = 2, cex.axis = 2)
ttt <-lapply(names(nsc_sc_by_day), function(di) {
    cdi <- nsc_sc_by_day[[di]]
    plot(ipc_sc[cdi], astro_sc[cdi], 
         # col = phase_clvls[cdi],
         col = mcmd$color[mc@mc[cdi]], 
         xlab = 'IPC (ds UMIs)', ylab = "Astro (ds UMIs)", main = paste0('E',di),
         pch = 16, cex = 1.25,
         xlim = c(0, quantile(ipc_sc[nsc_sc], 1)),
         ylim  = c(0, quantile(astro_sc[nsc_sc], 1)))
         # xlim = c(0, quantile(ipc_sc[nsc_sc], 0.999)),
         # ylim  = c(0, quantile(astro_sc[nsc_sc], 0.999)))
})

mcs_here <- union(nsc_mcs, which(mcmd$cell_type %in% c('IPC', 'IPC_cyc', 'Astrocytes')))
legc_mm_nsc <- legc - rowMeans(legc[,nsc_mcs])
mat_ipc_module <- legc_mm_nsc[ipc_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]]
mat_astro_module <- legc_mm_nsc[astro_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]]
mat_stem_module <- legc_mm_nsc[stem_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]] 
hc_mat_ipc <- hclust(dist(mat_ipc_module), method = 'ward.D2')
hc_mat_astro <- hclust(dist(mat_astro_module), method = 'ward.D2')
hc_mat_stem <- hclust(dist(mat_stem_module), method = 'ward.D2')

pltmt_modules <- rbind(mat_ipc_module[hc_mat_ipc$order,], mat_astro_module[hc_mat_astro$order,], mat_stem_module[hc_mat_stem$order,])

marks <- names(scdb_gset('pl_cort_marks_f')@gene_set)
marks


# ipc_genes_not_in_goi <- setdiff(lr[lr %in% ipc_module], lr[lr %in% goi])
ipc_genes_not_in_goi <- setdiff(ipc_module, goi)
lr <- ifelse(rownames(pltmt_modules) %in% setdiff(marks, ipc_genes_not_in_goi), rownames(pltmt_modules), '')

p_pltmt_modules <-pheatmap::pheatmap(pltmt_modules, silent = T,
                                     cluster_cols = F, 
                                     fontsize_row = 14,
                                     annotation_col = col_annot, cluster_rows = F, 
                                     # labels_row = lr,
                                     # gaps_row = match(c(rownames(mat_ipc_module)[hc_mat_ipc$order][[nrow(mat_ipc_module)]],
                                     #                      rownames(mat_astro_module)[hc_mat_astro$order][[nrow(mat_astro_module)]]),
                                     #                      rownames(pltmt_modules)),
                                     gaps_col = match(unique(names(cust_mc_ord_st)[cust_mc_ord_st %in% mcs_here]), names(cust_mc_ord_st)[cust_mc_ord_st %in% mcs_here]) -1 ,
                                     annotation_legend = F, annotation_colors = ann_colors, 
                                     col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), 
                                     treeheight_row = 0, show_colnames = F)

source('./scripts/util.r')

p_pltmt_modules$gtable <- add.flag(pheatmap = p_pltmt_modules, kept.labels = lr, repel.degree = 0)

mcATAC::save_pheatmap(x = p_pltmt_modules,filename =  './output/metacell_model/nsc_gene_modules/figs/nsc_ipc_astro_stem_modules_phm_flag.png', h = 900, w = 850)

fpba <- factor(pba)


phase_cut_color2 <- factor(setNames(c(colorRampPalette(c('blueviolet', 'orangered'))(4), 'gray')[fpba], names(fpba)))

phase_color_key <- setNames(gplots::col2hex(c(colorRampPalette(c('blueviolet', 'orangered'))(4), 'gray')), levels(fpba))

png('./output/metacell_model/nsc_gene_modules/figs/pc2_vs_pc1_color_phase.png', h = 1200, w = 1350, res = 100)
par(cex.lab = 4, mar = c(8,8,1,1), cex.main = 6, cex.axis = 4)
plot(cc_mat_rot[,1], cc_mat_rot[,2], col = phase_color_key[fpba[rownames(cc_mat_rot)]], pch = 16, cex = 2,ylab = '', xlab = '', xaxt ='n')
llp <- length(phase_color_key)
legend('topleft', col = phase_color_key, legend = gsub('^\\d_', '', names(phase_color_key)), pch = rep(16,llp), cex = rep(4, llp))
axis(1, padj = 1)
title(ylab = 'PC2', line = 4)
title(xlab = 'PC1', line = 6)

dev.off()

nsc_inds <- which(sc_data_df$cell_type == 'NSC')
lupc <- length(unique(phase_cut))
bin_seq <- seq(lupc+0.5, 0.5+length(days)*lupc, 12)
png('./output/metacell_model/nsc_gene_modules/figs/modules_in_nsc_scs_by_day_and_phase.png', h = 1050, w = 1600)
par(mfrow = c(3,1), las = 2, cex.lab = 3, cex.axis = 1.5, mar = c(6,6,0.5,0.5))
boxplot(ipc ~ ., data = sc_data_df[nsc_inds,c('ipc','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$ipc[nsc_inds], 0.94, na.rm = T)),
       col = rep(rainbow(lupc), length(days)), ylab = 'IPC module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))
boxplot(astro ~ ., data = sc_data_df[nsc_inds,c('astro','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$astro[nsc_inds], 0.975, na.rm = T)),
              col = rep(rainbow(lupc), length(days)), ylab = 'Astro module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))
boxplot(stem ~ ., data = sc_data_df[nsc_inds,c('stem','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$stem[nsc_inds], 0.95, na.rm = T)),
              col = rep(rainbow(lupc), length(days)), ylab = 'Stem module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))

dev.off()

indsh <- rownames(sc_data_df)[sc_data_df$cell_type %in% c('NSC', 'IPC', 'IPC_cyc')]

pcu_nsc_ipc <- principal_curve(as.matrix(sc_data_df[indsh,c('ipc', 'nsc')]))

cells_astro_nsc <- rownames(sc_data_df)[sc_data_df$cell_type %in% c('Astrocytes', 'NSC')]
cells_astro_ss_nsc <- rownames(sc_data_df)[union(which(sc_data_df$cell_type == 'Astrocytes'), 
                                                 sample(x = which(sc_data_df$cell_type == 'NSC'), size = length(which(sc_data_df$cell_type == 'Astrocytes'))))]
cells_astro_nsc_f <- setdiff(cells_astro_nsc, rownames(sc_data_df)[sc_data_df$ipc >= 30])
indsh_f <- setdiff(indsh, rownames(sc_data_df)[sc_data_df$astro >= 30])

pcu_nsc_astro <- principal_curve(as.matrix(sc_data_df[cells_astro_ss_nsc,c('astro', 'nsc')]))

ord_pcu_nsc_astro <- order(apply(sc_data_df[cells_astro_ss_nsc,c('astro', 'nsc')], 1, diff), decreasing = T)


match_points_to_pcu <- function(x, pcu, frac_subsample = 0.1) {
    pcu <- pcu[order(pcu[,1]),]
    phase_pcu_rm_ss <- pcu[sort(sample(x = nrow(pcu), size = round(nrow(pcu)*frac_subsample), replace = F)),]
    dist_phase_pc_to_pcu <- as.matrix(tgstat::tgs_dist(rbind(x[,1:2], phase_pcu_rm_ss[,1:2])))
    rownames(dist_phase_pc_to_pcu) <- c(rownames(dist_phase_pc_to_pcu)[1:nrow(x)], paste0('p', 1:nrow(phase_pcu_rm_ss)))
    colnames(dist_phase_pc_to_pcu) <- c(colnames(dist_phase_pc_to_pcu)[1:nrow(x)], paste0('p', 1:nrow(phase_pcu_rm_ss)))
    knn_phase_pc_to_pcu <- tgstat::tgs_knn(1/dist_phase_pc_to_pcu[(nrow(x)+1):ncol(dist_phase_pc_to_pcu),1:nrow(x)], k = 1)
    rownames(dist_phase_pc_to_pcu)[(nrow(cc_mat_rot)+1):ncol(dist_phase_pc_to_pcu)]
    knn_phase_pc_to_pcu[,2] <- paste0('p', 1:nrow(phase_pcu_rm_ss))[match(knn_phase_pc_to_pcu[,2], rownames(dist_phase_pc_to_pcu))]
    matched_points <- setNames(as.numeric(gsub('^p', '', knn_phase_pc_to_pcu[,2]))[match(rownames(x), as.character(knn_phase_pc_to_pcu[,1]))],
                      rownames(x))
    return(matched_points)
}

matched_coords <- match_points_to_pcu(as.matrix(sc_data_df[indsh_f,c('ipc', 'nsc')]), pcu_nsc_ipc$s, frac_subsample = 0.1)

matched_coords_bins <- setNames(cut(matched_coords, breaks = round(seq(0, max(matched_coords), l = 20))), names(matched_coords))

matched_coords_astro <- match_points_to_pcu(as.matrix(sc_data_df[cells_astro_nsc_f,c('astro', 'nsc')]), pcu_nsc_astro$s[order(pcu_nsc_astro$s[,1]),], frac_subsample = 1)

matched_coords_bins_astro <- setNames(cut(matched_coords_astro, breaks = round(seq(0, max(matched_coords_astro), l = 20))), names(matched_coords_astro))


cti_inds <- which(names(matched_coords_bins) %in% names(mc@mc[mc@mc %in% which(mcmd$cell_type %in% c('NSC', 'IPC', 'IPC_cyc'))]))
tbl_phase_by_nsc_ipc_axis <- table(matched_coords_bins[cti_inds], sc_data_df[names(matched_coords_bins)[cti_inds],'pba'])
mat_phase_by_nsc_ipc_axis <- as.matrix(t(tbl_phase_by_nsc_ipc_axis))
mat_phase_norm <- t(t(mat_phase_by_nsc_ipc_axis)/colSums(mat_phase_by_nsc_ipc_axis))

colnames(mat_phase_norm) <- 1:ncol(mat_phase_norm)
colnames(mat_phase_by_nsc_ipc_axis) <- 1:ncol(mat_phase_by_nsc_ipc_axis)

tbl_ct_by_princurve_bin <- table(matched_coords_bins[cti_inds], sc_data_df[names(matched_coords_bins)[cti_inds],'cell_type'])

mat_ct_by_princurve_bin <- t(as.matrix(tbl_ct_by_princurve_bin))
mat_ct_by_princurve_bin_norm <- t(t(mat_ct_by_princurve_bin)/colSums(mat_ct_by_princurve_bin))
colnames(mat_ct_by_princurve_bin_norm) <- 1:ncol(mat_ct_by_princurve_bin_norm)
colnames(mat_ct_by_princurve_bin) <- 1:ncol(mat_ct_by_princurve_bin)


tbl_phase_by_nsc_astro_axis <- table(matched_coords_bins_astro, sc_data_df[names(matched_coords_bins_astro),'pba'])
mat_phase_by_nsc_astro_axis <- as.matrix(t(tbl_phase_by_nsc_astro_axis))
mat_phase_astro_norm <- t(t(mat_phase_by_nsc_astro_axis)/colSums(mat_phase_by_nsc_astro_axis))

colnames(mat_phase_astro_norm) <- 1:ncol(mat_phase_astro_norm)
colnames(mat_phase_by_nsc_astro_axis) <- 1:ncol(mat_phase_by_nsc_astro_axis)

tbl_ct_by_princurve_bin_astro <- table(matched_coords_bins_astro, sc_data_df[names(matched_coords_bins_astro),'cell_type'])

mat_ct_by_princurve_bin_astro <- t(as.matrix(tbl_ct_by_princurve_bin_astro))
mat_ct_by_princurve_bin_astro_norm <- t(t(mat_ct_by_princurve_bin_astro)/colSums(mat_ct_by_princurve_bin_astro))
colnames(mat_ct_by_princurve_bin_astro_norm) <- 1:ncol(mat_ct_by_princurve_bin_astro_norm)
colnames(mat_ct_by_princurve_bin_astro) <- 1:ncol(mat_ct_by_princurve_bin_astro)

png(glue::glue('./output/metacell_model/nsc_gene_modules/figs/phase_probs_by_astro_and_ipc_axes_princurve_bin_barplots.png'), h = 800, w = 700)
par(mar = c(4,3,3,2), cex.main = 1.5, cex.axis = 1.5)
nf <- graphics::layout(mat = matrix(c(1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8), nrow = 6, ncol = 2))

# row 1, col 1
plot(log2((1+mat_phase_by_nsc_ipc_axis['2_S',])/(1+mat_phase_by_nsc_ipc_axis['4_M',])), main = 'log2 S/M fraction ratio - IPC-NSC axis', type= 'l', ylab = '')
lines(c(-10,30), c(0,0), col = 'red', lty = 2, lwd = 2)

# row 2, col 1
barplot(mat_phase_norm[c(5,1:4),], col = phase_color_key[c(5,1:4)], ylab = '', main = 'Phase fractions by IPC-NSC princurve bin', add = F)

# row 3, col 1
barplot(mat_ct_by_princurve_bin_norm, col = col_key[rownames(mat_ct_by_princurve_bin_norm)], ylab = '', main = 'Cell type fractions by IPC-NSC princurve bin', add = F)

# row 4, col 1
yh <- log2(colSums(mat_phase_by_nsc_ipc_axis[2:4,])/colSums(mat_phase_by_nsc_ipc_axis[c(1,5),])) 
plot(yh, main = 'log2 cycling/non-cycling ratio - IPC-NSC axis', , type = 'l', xlab = '')
lines(c(-10,30), c(0,0), col = 'red', lty = 2, lwd = 2)

# row 1, col 2
par(mar = c(4,3,3,2), cex.main = 1.5, cex.axis = 1.5)
plot(log2((1+mat_phase_by_nsc_astro_axis['2_S',])/(1+mat_phase_by_nsc_astro_axis['4_M',])), main = 'log2 S/M fraction ratio - astro-NSC axis', type= 'l', ylab = '')
lines(c(-10,30), c(0,0), col = 'red', lty = 2, lwd = 2)

# row 2, col 2
barplot(mat_phase_astro_norm[c(5,1:4),], col = phase_color_key[c(5,1:4)], ylab = '', main = 'Phase fractions by Astro-NSC princurve bin', add = F)
# legend(x = 23,y = 1, legend = names(phase_color_key), col = phase_color_key, pch = 15, cex = 1.5, xpd = T)

# row 3, col 2
barplot(mat_ct_by_princurve_bin_astro_norm, col = col_key[rownames(mat_ct_by_princurve_bin_astro_norm)], ylab = '', main = 'Cell type fractions by Astro-NSC princurve bin', add = F)

# row 4, col 2
yh <- log2(colSums(mat_phase_by_nsc_astro_axis[2:4,])/colSums(mat_phase_by_nsc_astro_axis[c(1,5),])) 
plot(yh, main = 'log2 cycling/non-cycling ratio - astro-NSC axis', ylim = max(abs(yh))*c(-1,1), type = 'l', xlab = '')
lines(c(-10,30), c(0,0), col = 'red', lty = 2, lwd = 2)


dev.off()

png('./output/metacell_model/nsc_gene_modules/figs/nsc_phase_pc1_pc2_scatter_by_day.png', w = 900, h = 700)
mari <- c(5,.5,3,.5)
par(mfrow = c(2,3), cex.lab = 3, mar = mari, cex.axis = 2, cex.main = 3)
ttt <- lapply(names(nsc_sc_by_day), function(di) {
    if (di %in% c('13', '16')) {
        mari[[2]] <- 6; 
        ylabi = 'PC2'; yaxti = 's'
    } else {mari[[2]] <- 1; ylabi = ''; yaxti = 'n'}
    par(mar = mari)
    plot(cc_mat_rot[nsc_sc_by_day[[di]],1], cc_mat_rot[nsc_sc_by_day[[di]],2], pch = 16, cex = 1.3, 
         col = col_key[['NSC']],
         yaxt = yaxti,
         main = paste0('E', di), xlab = 'PC1', ylab = ylabi, 
         xlim = c(min(cc_mat_rot[,1]), max(cc_mat_rot[,1])), 
        ylim = c(min(cc_mat_rot[,2]), max(cc_mat_rot[,2])))
    sci <- nsc_sc_by_day[[di]]
    text(-10,150,adj = c(0,0.3), labels = paste0(100*round(length(which(fpba[sci] == '5_G0'))/length(sci), 2), '% G0'), cex = 3)
})
dev.off()


dir.create('./output/metacell_model/figs/nsc_gene_module_analysis')
dir.create('./output/metacell_model/figs/nsc_gene_module_analysis/pheatmaps')
dir.create('./output/metacell_model/figs/nsc_gene_module_analysis/boxplots')


vvv <- lapply(sort(unique(ct_hc_cor_nsc)), function(clj) {
    gnj <- names(ct_hc_cor_nsc)[ct_hc_cor_nsc == clj]
    # png(glue::glue('./output/metacell_model/figs/nsc_gene_module_analysis/boxplots/cluster_{clj}.png'), h = 600, w = 800)
    par(las = 2, cex.main = 2, cex.lab = 3,cex.axis = 2,mar = c(14,4,4,0.5))
    boxplot(st_legc[gnj,cust_st_ord2], col = col_key[cust_st_ord2], main = paste('Cluster', clj))
    # dev.off()
})


vvv <- lapply(sort(unique(ct_hc_cor_nsc)), function(clj) {
    gnj <- names(ct_hc_cor_nsc)[ct_hc_cor_nsc == clj]
    pltmtj <- legc[gnj,cust_mc_ord_st2] - rowMeans(legc[gnj,nsc_mcs])
    p_clj <- pheatmap::pheatmap(pltmtj, color = clrmp_rel, treeheight_row = 0, fontsize_row = 8, show_colnames = F,
                       # breaks = seq(-16.6,-8,l=100), 
                       breaks = seq(-3,3,l=100), 
                       annotation_col = col_annot, annotation_colors = ann_colors, cluster_cols = F, silent = T)
    # save_pheatmap_png(p_clj, glue::glue('./output/metacell_model/figs/nsc_gene_module_analysis/pheatmaps/cluster_{clj}.png'), h = 16*nrow(pltmtj), w = 2000)
})