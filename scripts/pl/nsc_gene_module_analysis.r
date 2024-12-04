# ls()

# devtools::load_all('~/src/metacell')
library(metacell)
# devtools::load_all('~/src/metacell.flow')
library(pheatmap)
library(ComplexHeatmap)
library(matrixStats)
library(tgstat)

wd = '/home/feshap/raid/proj/mmcortex'
setwd(wd)
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
# scdb_flow_init()
SEED = 1337
K = 16
set.seed(SEED)
scfigs_init("figs/")
doMC::registerDoMC(77)
nm = 'pl_cort'


mc = scdb_mc(nm)
mat = scdb_mat(nm)

source('./scripts/util.r')

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
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

feats = scdb_gset('pl_filt_lat')
feats = names(feats@gene_set)
feats = feats[feats %in% rownames(mc@e_gc)]

feats_all <- scdb_gset('pl')

nsc_mcs <- which(mcmd$cell_type == 'NSC')

ipc_cyc_mcs <- which(mcmd$cell_type == 'IPC_cyc')
ipc_mcs <- which(mcmd$cell_type == 'IPC')

nsc_legc_gq <- t(apply(legc[,nsc_mcs], 1, quantile, c(0,0.01,0.05,0.1,0.9, 0.95,0.99,1)))


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

cor_nsc_legc_dyn_genes <- tgs_cor(t(legc[nsc_dyn_genes,nsc_mcs]), spearman = T)

hc_cor_nsc <- hclust(dist(cor_nsc_legc_dyn_genes), method = 'ward.D2')

ct_hc_cor_nsc <- cutree(hc_cor_nsc, k = 8)

rna_mc_day <- as.matrix(mcmd[,grep('E\\d\\d', colnames(mcmd))])
rna_mc_day_norm <- rna_mc_day/rowSums(rna_mc_day)
nsc_mcs <- which(mcmd$cell_type == 'NSC')
egc_by_day <- mc@e_gc[,nsc_mcs] %*% rna_mc_day_norm[nsc_mcs,]
egc_by_day_n <- t(t(egc_by_day)/colSums(egc_by_day))
# colnames(egc_by_day_n) <- paste0('E', colnames(egc_by_day_n))
legc_by_day_n <- log2(1e-5 + egc_by_day_n)

astro_mcs <- which(mcmd$cell_type == 'Astrocytes')


load('./output/metacell_model/pcurve_nsc_vs_mat_neuro.rda')

ro <- rev(pcu$ord)

ra <- as.data.frame(tibble::column_to_rownames(tibble::enframe(ct_hc_cor_nsc, name = 'gene', value = 'cluster'), 'gene'))

ac <- list(cluster = setNames(chameleon::distinct_colors(8)$name, 1:8))

# p_cor_dyn_genes_nsc <- pheatmap::pheatmap(cor_nsc_legc_dyn_genes[hc_cor_nsc$order,hc_cor_nsc$order], show_rownames = F, show_colnames = F, annotation_legend = F,
#                    col = colorRampPalette(c('blue3', 'white', 'red3'))(100),, breaks = seq(-1,1,l=100),
#                    cluster_cols = F, cluster_rows = F, 
#                    annotation_col = ra, annotation_row = ra, annotation_colors = ac)
# save_pheatmap_png(p_cor_dyn_genes_nsc, './output/metacell_model/figs/nsc_gene_module_analysis/cor_dyn_genes_nsc.png', h = 1600, w = 1600)


ipc_cyc_vs_ipc_genes <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'IPC_cyc'), mc_neg = which(mcmd$cell_type == 'IPC'))

ipc_genes_vs_glia <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('IPC_cyc', 'IPC')), mc_neg = which(mcmd$cell_type %in% c('Astrocytes', 'OPCs')))
ipc_genes_vs_nsc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in% c('IPC_cyc', 'IPC')), mc_neg = which(mcmd$cell_type %in% c('NSC')))

ipc_module <- intersect(names(ipc_genes_vs_glia[ipc_genes_vs_glia >= 2]),
                    names(ipc_genes_vs_nsc[ipc_genes_vs_nsc >= 2]))

astro_genes_vs_oligo <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('OPCs')))
astro_genes_vs_nsc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('NSC')))
astro_genes_vs_ipc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type == 'Astrocytes'), mc_neg = which(mcmd$cell_type %in% c('IPC', 'IPC_cyc')))

astro_module <- intersect(names(astro_genes_vs_oligo[astro_genes_vs_oligo >= 2]),
                      intersect(names(astro_genes_vs_nsc[astro_genes_vs_nsc >= 2]),
                           names(astro_genes_vs_ipc[astro_genes_vs_ipc >= 2])))

stem_genes_nsc_fc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in%  c('NSC')), mc_neg = which(mcmd$cell_type %in%  c('Astrocytes', 'IPC')))
stem_genes_ipc_fc <- get_genes_specific_to_mcs(legc, mc_pos = which(mcmd$cell_type %in%  c('IPC_cyc')), mc_neg = which(mcmd$cell_type %in%  c('Astrocytes', 'IPC')))

cond1 <- st_legc[,'NSC'] - st_legc[,'Astrocytes'] >= 1
cond2 <- st_legc[,'NSC'] - st_legc[,'IPC'] >= 1

stem_genes <- rownames(st_legc)[(cond1 & cond2)]

cor_stem_cc_nsc <- tgs_cor(t(legc[stem_genes,which(mcmd$cell_type %in%  c('NSC'))]), 
                       t(legc[c('Top2a', 'Mki67', 'Mcm4', 'Pcna'), which(mcmd$cell_type %in%  c('NSC'))]))
cor_stem_cc_ipc_cyc <- tgs_cor(t(legc[stem_genes,which(mcmd$cell_type %in%  c('IPC_cyc'))]), 
                       t(legc[c('Top2a', 'Mki67', 'Mcm4', 'Pcna'), which(mcmd$cell_type %in%  c('IPC_cyc'))]))

cor_stem_cc_nsc_filt <- cor_stem_cc_nsc[matrixStats::rowMaxs(cor_stem_cc_nsc) < .5,]
cor_stem_cc_ipc_cyc_filt <- cor_stem_cc_ipc_cyc[matrixStats::rowMaxs(cor_stem_cc_ipc_cyc) < .5,]
stem_module <- union(rownames(cor_stem_cc_nsc_filt), rownames(cor_stem_cc_ipc_cyc_filt))

write(astro_module, './output/metacell_model/nsc_gene_modules/astro_module.txt')

write(ipc_module, './output/metacell_model/nsc_gene_modules/ipc_module.txt')

write(stem_module, './output/metacell_model/nsc_gene_modules/stem_module.txt')

astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')

ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')

stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')


legc_nsc_dyn_genes_avg_cl <- tgs_matrix_tapply(t(legc[nsc_dyn_genes,nsc_mcs]), ct_hc_cor_nsc, mean)

set.seed(1337)
mat_ds <- scm_downsamp(mat@mat, 3000)
save(mat_ds, file = './output/metacell_model/mat_ds.rda')


nsc_sc <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'NSC')]), colnames(mat_ds))

nsc_genes <- sort(rownames(all_ct_genes[['NSC']])[all_ct_genes[['NSC']][,'Astrocytes'] <= -1 & all_ct_genes[['NSC']][,'IPC_cyc'] <= -1])

cl1_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 1])
cl6_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 6])
cl7_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 7])
cl3_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 3])
cl13_genes <- names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3)])
cl1_sc <- Matrix::colSums(mat_ds[cl1_genes,])
cl13_sc <- Matrix::colSums(mat_ds[cl13_genes,])
cl3_sc <- Matrix::colSums(mat_ds[cl3_genes,])
# cl6_sc <- colSums(mat_ds[union(cl6_genes, cl3_genes),])
cl6_sc <- Matrix::colSums(mat_ds[cl6_genes,])
cl7_sc <- Matrix::colSums(mat_ds[cl7_genes,])
cc_sc <- Matrix::colSums(mat_ds[c(cl1_genes, cl3_genes, cl6_genes, cl7_genes),])
ipc_sc <- Matrix::colSums(mat_ds[ipc_module,])
astro_sc <- Matrix::colSums(mat_ds[astro_module,])
stem_sc <- Matrix::colSums(mat_ds[stem_module,])
nsc_genes_sc <- Matrix::colSums(mat_ds[nsc_genes,])
# oligo_sc <- colSums(mat_ds[oligo_module,])

cc_data_mat <- cbind(cl1_sc - mean(cl1_sc), cl3_sc - mean(cl3_sc), cl6_sc - mean(cl6_sc), cl7_sc - mean(cl7_sc))
pca_cc <- prcomp(x = cc_data_mat)
cc_pca_mat <- pca_cc$x
theta <- 140*pi/180
cc_mat_rot <- cc_pca_mat %*% rbind(c(cos(theta), -sin(theta), 0, 0), c(sin(theta), cos(theta), 0, 0), c(0, 0, 1,  0), c(0, 0, 0, 1))

cycling_scs <- which(cc_sc >= quantile(cc_sc, 0.8))

# phase_rad <- phase_rad_17
phase_rad <- -atan2(cc_mat_rot[,2] - mean(cc_mat_rot[names(cycling_scs),2]), cc_mat_rot[,1] - mean(cc_mat_rot[names(cycling_scs),1]))

phase_clvls <- setNames(clrmp[1+round(999*(phase_rad[rownames(cc_mat_rot)] - min(phase_rad))/(max(phase_rad) - min(phase_rad)))],
                        rownames(cc_mat_rot))

clrmp_mod <- colorRampPalette(c('white', 'orange', 'purple', 'black'))(1000)

phase_clvls_mod <- setNames(clrmp_mod[1+round(999*(phase_rad[rownames(cc_mat_rot)] - min(phase_rad))/(max(phase_rad) - min(phase_rad)))],
                        rownames(cc_mat_rot))

library(princurve)

nsc_ipc_cyc_cells_in_cc_mat <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type %in% c('NSC', 'IPC_cyc'))]), rownames(cc_mat_rot))

phase_pcu <- principal_curve(x = cc_mat_rot[nsc_ipc_cyc_cells_in_cc_mat[order(phase_rad[match(nsc_ipc_cyc_cells_in_cc_mat, names(phase_rad))])],1:2], smoother = "periodic_lowess", stretch = 1, approx_points = 100, thresh = 1e-5)

phase_pcu_rm <- zoo::rollmean(x = phase_pcu$s[,1:2], k = 50)

phase_pcu_rm_ss <- phase_pcu_rm[sort(sample(x = nrow(phase_pcu_rm), size = round(nrow(phase_pcu_rm)/10), replace = F)),]


dist_phase_pc_to_pcu <- as.matrix(tgstat::tgs_dist(rbind(cc_mat_rot[,1:2], phase_pcu_rm_ss[,1:2])))

rownames(dist_phase_pc_to_pcu) <- c(rownames(dist_phase_pc_to_pcu)[1:nrow(cc_mat_rot)], paste0('p', 1:nrow(phase_pcu_rm_ss)))
colnames(dist_phase_pc_to_pcu) <- c(colnames(dist_phase_pc_to_pcu)[1:nrow(cc_mat_rot)], paste0('p', 1:nrow(phase_pcu_rm_ss)))

knn_phase_pc_to_pcu <- tgstat::tgs_knn(1/dist_phase_pc_to_pcu[(nrow(cc_mat_rot)+1):ncol(dist_phase_pc_to_pcu),1:nrow(cc_mat_rot)], k = 1)

knn_phase_pc_to_pcu[,2] <- paste0('p', 1:nrow(phase_pcu_rm_ss))[match(knn_phase_pc_to_pcu[,2], rownames(dist_phase_pc_to_pcu))]

phase <- setNames(as.numeric(gsub('^p', '', knn_phase_pc_to_pcu[,2]))[match(rownames(cc_mat_rot), as.character(knn_phase_pc_to_pcu[,1]))],
                  rownames(cc_mat_rot))

phase_rev <- max(phase) - phase + 1

NUM_PARTITION <- 13
phase_qs <- seq(1-1e-2,max(phase),l=NUM_PARTITION)

phase_cut <- setNames(cut(phase_rev, phase_qs, labels = paste('bin', 1:(NUM_PARTITION - 1))), names(phase_rev))

save(phase, phase_rev, phase_cut, file = './output/metacell_model/nsc_gene_modules/phase_info.rda')

load('./output/metacell_model/nsc_gene_modules/phase_info.rda')

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

fpba <- factor(pba)

fpba_by_day_mat <- sapply(tail(unique(mat@cell_metadata$day),-1), function(di) {cells_di <- intersect(colnames(mat@mat), rownames(mat@cell_metadata[mat@cell_metadata$day == di,]));
                                                             cells_di <- intersect(cells_di, names(mc@mc[mc@mc %in% nsc_mcs]));
                                                   table(fpba[cells_di])})
fpba_by_day_mat_norm <- t(t(fpba_by_day_mat)/colSums(fpba_by_day_mat))



sc_data_df <- dplyr::select(mat@cell_metadata[colnames(mat_ds),], day)
sc_data_df$cell_type <- mcmd$cell_type[mc@mc[rownames(sc_data_df)]]
sc_data_df$ipc <- ipc_sc[rownames(sc_data_df)]
sc_data_df$astro <- astro_sc[rownames(sc_data_df)]
sc_data_df$nsc <- nsc_genes_sc[rownames(sc_data_df)]
sc_data_df$stem <- stem_sc[rownames(sc_data_df)]
sc_data_df$cc <- cc_sc[rownames(sc_data_df)]
sc_data_df$phase <- phase[rownames(sc_data_df)]
sc_data_df$phase_rev <- phase_rev[rownames(sc_data_df)]
# sc_data_df$oligo <- oligo_sc[rownames(sc_data_df)]
sc_data_df$phase_cut <- as.numeric(phase_cut[rownames(sc_data_df)])

sc_data_df <- dplyr::mutate(sc_data_df, pba = pba)

tbl_pba_by_ct <- t(table(sc_data_df$cell_type, sc_data_df$pba))
tbl_pba_by_ct_norm <- t(t(tbl_pba_by_ct)/colSums(tbl_pba_by_ct))
tbl_pba_by_ct_norm

rownames(tbl_pba_by_ct_norm) <- gsub('\\d_', '', rownames(tbl_pba_by_ct_norm))

mg_bon_marks <- as.data.frame(t(sapply(apply(readr::read_csv(file.path(wd, 'BonevCollab//marker_genes.tsv')), 1, 
                                             stringr::str_split, ' '), function(x) c(x[[1]][[1]], x[[1]][[length(x[[1]])]]))))

colnames(mg_bon_marks) <- c('cell_type', 'marks')

mbm_lst <- lapply(1:nrow(mg_bon_marks), function(n) stringr::str_split(mg_bon_marks$marks[[n]], ',')[[1]])
names(mbm_lst) <- mg_bon_marks$cell_type

sc_data_df$nsc <- Matrix::colSums(mat_ds[mbm_lst[['NSC']],rownames(sc_data_df)])
sc_data_df$color <- col_key[sc_data_df$cell_type]

indsh <- rownames(sc_data_df)[sc_data_df$cell_type %in% c('NSC', 'IPC', 'IPC_cyc')]

pcu_nsc_ipc <- principal_curve(as.matrix(sc_data_df[indsh,c('ipc', 'nsc')]))

cells_astro_nsc <- rownames(sc_data_df)[sc_data_df$cell_type %in% c('Astrocytes', 'NSC')]
cells_astro_ss_nsc <- rownames(sc_data_df)[union(which(sc_data_df$cell_type == 'Astrocytes'), 
                                                 sample(x = which(sc_data_df$cell_type == 'NSC'), size = length(which(sc_data_df$cell_type == 'Astrocytes'))))]
cells_astro_nsc_f <- setdiff(cells_astro_nsc, rownames(sc_data_df)[sc_data_df$ipc >= 30])
indsh_f <- setdiff(indsh, rownames(sc_data_df)[sc_data_df$astro >= 30])

pcu_nsc_astro <- principal_curve(as.matrix(sc_data_df[cells_astro_ss_nsc,c('astro', 'nsc')]))


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


# cti_inds <- which(names(matched_coords_bins_astro) %in% names(mc@mc[mc@mc %in% which(mcmd$cell_type %in% c('NSC', 'A'))]))
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



save(mat_ds,
    pba, 
    cc_mat_rot, 
    ct_hc_cor_nsc,
    mat_phase_by_nsc_ipc_axis, 
    mat_phase_norm,
    mat_ct_by_princurve_bin_norm,
    mat_phase_by_nsc_astro_axis,
    mat_phase_astro_norm,
    mat_ct_by_princurve_bin_astro_norm,
file = './output/metacell_model/nsc_gene_modules/fig2_data.rda')



save(nsc_dyn_genes,
    cor_nsc_legc_dyn_genes,
    mat_ds,
    hc_cor_nsc,
    phase,
    phase_cut,
    phase_pcu_rm_ss,
    sc_data_df,
    mat_ds_cc_genes_select_ord_phase,
    mat_ds_cc_genes_select_ord_phase_rm,
    cc_mat_rot, 
    ct_hc_cor_nsc,
    matched_coords_bins,
    matched_coords_bins_astro,
    legc_by_day_n,
file = './output/metacell_model/nsc_gene_modules/figs2_data.rda')

