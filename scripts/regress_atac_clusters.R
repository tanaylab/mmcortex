### Script to compare clusters obtained from scATAC-derived microclusters 
## and microclusters derived from scATAC-RNA metacell correlations
wd <- "/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/"
setwd(wd)
library(matrixStats)
library(pheatmap)
library(metacell)
scdb_init('scdb', f=T)
mat_prom <- scdb_mat('pl_prom_cort')
options(gmax.data.size = 1e+9)
devtools::load_all("~/src/mcATAC/")
library(prego)
gset_genome('mm10')
doMC::registerDoMC(cores = 70)
mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))                         
col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))


aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
cl_raw <- aaa$km_a_legc
cl_rna <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
mcl_all <- aaa$mcl_all


mat_tbl <- as.matrix(table(cl_raw$cluster, cl_rna$cluster))
colnames(mat_tbl) <- paste0('raw_', colnames(mat_tbl))
rownames(mat_tbl) <- paste0('rna_', rownames(mat_tbl))
mat_tbl_norm_col <- t(t(mat_tbl)/colSums(mat_tbl))


ord <- order(apply(mat_tbl, 1, sum))
ord_col <- order(colSums(mat_tbl))
ord_row_norm <- order(rowSums(mat_tbl_norm_col))
ord_col_norm <- order(colSums(mat_tbl_norm_col))
ord_row_max <- order(apply(mat_tbl, 1, max))
ord_row_max_norm <- order(apply(mat_tbl_norm_col, 1, max))
ord_col_max <- order(apply(mat_tbl, 2, max))

sc_cor_kms <- readRDS('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/microcluster_assignment.RDS')

day_mcl_prom <- t(do.call('rbind', 
                    lapply(1:length(sc_cor_kms), 
                            function(i) tgs_matrix_tapply(as.matrix(mat_prom@mat[,colnames(aaa$sc_res[[i]]$sc_cor)]), sc_cor_kms[[i]]$cluster, sum)
                            )))

nm <- "pl_prom_cort"

gstat = scdb_gstat(nm)

x = log(gstat$ds_mean)
init_filt = which(x >= -4)
y = gstat$ds_log_varmean[init_filt]
name_filt = gstat$name
x = x[init_filt]
xcut = cut(x, breaks = seq(min(x), max(x), l = 50))
top_q_inds = lapply(levels(xcut), function(l) {inds = which(xcut == l); 
                                  xfilt = y[inds]; 
                                  xtop = head(inds[order(xfilt, decreasing = T)], 30); 
                                  return(xtop)
                                 }
      )

names(top_q_inds) = levels(xcut)
mc_rna <- scdb_mc('pl_cort')
feats = intersect(rownames(gstat)[init_filt[unlist(top_q_inds)]], rownames(mc_rna@mc_fp))
cor_mcl_mc <- tgs_cor(day_mcl_prom[feats,], mc_rna@mc_fp[feats,], spearman=T)
rownames(cor_mcl_mc) <- sapply(13:18, function(x) paste0(x, '_', 1:30))

order_columns_hc_com <- function(mat, hc_cut_deg) {
    hcm <- hclust(dist(t(mat)), method = 'ward.D2')
    hc_ct <- cutree(hcm, hc_cut_deg)
    com_hc <- sapply(unique(hc_ct), function(ci) {
                        vec <- rowMeans(mat[,hc_ct == ci]); 
                        return(sum(vec*1:length(vec))/sum(vec))
                        })
    ord <- as.numeric(do.call('c', lapply(unique(hc_ct)[order(com_hc)], function(ci) which(hc_ct == ci))))
    return(ord)
}

amc_st_assn <- mcmd$cell_type[apply(cor_mcl_mc, 1, which.max)]
amc_cat_assn <- ifelse(amc_st_assn %in% c('NSC', 'IPC', 'IPC_cyc', 'Astrocytes', 'Oligodendrocytes'), amc_st_assn, 'Neuron')
amc_cat_assn <- ifelse(amc_st_assn %in% c('Astrocytes', 'Oligodendrocytes'), 'NSC', amc_cat_assn)
amc_ord <- do.call('c', lapply(c('NSC' ,'IPC', 'Neuron'), function(ct) {
    vec <- which(amc_cat_assn == ct);
    com_ord <- order_columns_hc_com(t(cor_mcl_mc[vec,]), 5);
    return(setNames(vec[com_ord], rep(ct, length(vec))))
    }))

row_annot <- as.data.frame(amc_cat_assn)
rownames(row_annot) <- colnames(mcl_all)
colnames(row_annot) <- 'cell_type'
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))
ann_colors[['cell_type']] <- c(ann_colors[['cell_type']], setNames('yellow', 'Neuron'))
pltmt <- cor_mcl_mc
rownames(pltmt) <- colnames(mcl_all)

mcl_all_norm <- t(t(mcl_all)/Matrix::colSums(mcl_all))
mcl_all_sum_atac <- t(log2(1e-5 + tgs_matrix_tapply(t(mcl_all_norm), cl_raw$cluster, mean)))
var_thresh_up <- -15
var_thresh_dn <- -16
varp_atac <- which(colMaxs(mcl_all_sum_atac) > var_thresh_up & 
                    colMins(mcl_all_sum_atac) < var_thresh_dn)

pltmt_atac <- mcl_all_sum_atac[amc_ord,varp_atac]
ord_col_atac <- order_columns_hc_com(2**pltmt_atac, 4)

sp_varp_atac <- do.call('c', lapply(varp_atac[ord_col_atac], function(pcl) {vec <- rownames(mcl_all)[which(cl_raw$cluster == pcl)]; 
                                    setNames(vec, rep(pcl, length(vec)))}))

### prego regression 
# get sequences
seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(mcl_all))
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# Identify promoter-proximal peaks from each group
tss <- gintervals.load('intervs.global.tss')
nei_sp_atac_prom <- gintervals.neighbors(seq_coords[match(sp_varp_atac, seq_coords$peak_name),], tss, maxneighbors = 1, mindist = -1e3, maxdist = 1e3)
sp_varp_atac_prom <- sp_varp_atac[sp_varp_atac %in% nei_sp_atac_prom$peak_name]
sp_varp_atac_dist <- sp_varp_atac[!(sp_varp_atac %in% nei_sp_atac_prom$peak_name)]

rsp_mat <- sapply(unique(names(sp_varp_atac_dist)), function(u) ifelse(names(sp_varp_atac_dist) == u, 1, 0))

seq_select <- seqs_all[match(sp_varp_atac_dist, rownames(mcl_all))]
rm(list = ls()[!(ls() %in% c(grep('G', ls(), v=T), 'wd', 'seq_select', 'rsp_mat'))])
print(head(seq_select))
print(length(seq_select))
print(head(rsp_mat))
print(dim(rsp_mat))
# regress on ATAC peak clusters
res <- plyr::llply(1:ncol(rsp_mat), function(i) prego::regress_pwm.cv(sequences=seq_select, 
                                   response = rsp_mat[,i], 
                                   nfolds = 5,
                                        min_kmer_cor = 0.01,
                                        unif_prior = 0.1,
                                        score_metric = "r2",
                                        final_metric = "ks",
                                        multi_kmers = F,
                                        match_with_db=T, 
                                   use_sample = T, 
                                        parallel=T
                                        ), .parallel = T)

save(res, file = './output/sequence_modeling/reg_atac_clusters_score_metric_r2_nfold_5.rda')


### Make ATAC cluster fastas

library(seqinr)
dir.create('./output/sequence_modeling/atac_cluster_fasta')
write_cluster_fasta <- function(cluster,name = cluster,sequences, sequence_clusters) {
    print(cluster)
    fg_names <- as.character(sequence_clusters[names(sequence_clusters) %in% cluster])
    print(head(fg_names))
    fg_seqs <- as.list(as.character(sequences[fg_names]))
    print(head(fg_seqs))
    # fg_fo <- paste0('./output/sequence_modeling/atac_cluster_fasta/', cluster, '_fg.fa')
    fg_fo <- paste0('./output/sequence_modeling/atac_cluster_fasta/', name, '_fg.fa')
    bg_names <- as.character(sequence_clusters[!(names(sequence_clusters) %in% cluster)])
    bg_seqs <- as.list(as.character(sequences[bg_names]))
    # print(head(bg_seqs))
    # bg_fo <- paste0('./output/sequence_modeling/atac_cluster_fasta/', cluster, '_bg.fa')
    bg_fo <- paste0('./output/sequence_modeling/atac_cluster_fasta/', name, '_bg.fa')
    write.fasta(sequences = fg_seqs, names = fg_names, file.out = fg_fo, open = "w", nbchar = 60, as.string = T)
    write.fasta(sequences = bg_seqs, names = bg_names, file.out = bg_fo, open = "w", nbchar = 60, as.string = TRUE)
}
names(seqs_all) <- rownames(mcl_all)
write_cluster_fasta(c), name = 'IPC', seqs_all, sp_varp_atac)

un <- unique(names(sp_varp_atac))
plyr::laply(un[!(un %in% c("1", "7"))], function(u) {
    dirn <- paste0("output/sequence_modeling/motifResults/", u)
    if (!(dir.exists(dirn))) { dir.create(dirn)}
    cmd <- paste0("~/src/homer/bin/findMotifs.pl output/sequence_modeling/atac_cluster_fasta/", u, "_fg.fa fasta output/sequence_modeling/motifResults/", u, "/ -fasta output/sequence_modeling/atac_cluster_fasta/", u, "_bg.fa -noknown")
    # print(cmd)
    system(cmd)
}, .parallel = T)

