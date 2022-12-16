library(dplyr)
library(pheatmap)
library(metacell)
library(matrixStats)
scdb_init('scdb', f=T)
devtools::load_all("~/src/mcATAC/")
doMC::registerDoMC(cores = 70)
misha.ext::gset_genome('mm10')
load("~/raid/proj/mmcortex/output/mcatac/peak_indices_var_atac_clust.rda")
intervs_energy <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')
ie_mat <- subset(intervs_energy, select = -c(chrom, start, end, peak_name))
rownames(ie_mat) <- peak_names(intervs_energy[,c('chrom', 'start', 'end')], tad_based =F)
aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
scatac_km <- aaa$km_a_legc

seq_coords <- intervs_energy[,c('chrom', 'start', 'end')]
seq_coords$peak_name <- peak_names(seq_coords, tad_based = F)
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100) %>% 
                select(chrom, start, end) 
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# Identify promoter-proximal peaks from each group
tss <- gintervals.load('intervs.global.tss')
nei_sp_atac_prom <- gintervals.neighbors(seq_coords[match(sp_varp_atac, seq_coords$peak_name),], tss, maxneighbors = 1, mindist = -1e3, maxdist = 1e3)
sp_varp_atac_prom <- sp_varp_atac[sp_varp_atac %in% nei_sp_atac_prom$peak_name]
sp_varp_atac_dist <- sp_varp_atac[!(sp_varp_atac %in% nei_sp_atac_prom$peak_name)]

mc_rna <- scdb_mc('pl_cort')
tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- c('Bcl11b', tfs$Symbol)
tfs <- tfs[tfs %in% rownames(mc_rna@e_gc)]
fp_tfs <- tfs[tfs %in% rownames(mc_rna@mc_fp)[rowMaxs(mc_rna@mc_fp) > 3]]
fp_tf_motif <- unlist(plyr::llply(fp_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(ie_mat), v=T, ign=T), .parallel = T))

calculate_d_in_cluster <- function(vec, cl_vec,
                 parallel = getOption("mcatac.parallel"), 
                 alternative = "less", nc = getOption("mcatac.parallel.nc")) {
    # print(head(vec))
    # print(head(cl_vec))
    return(plyr::laply(unique(cl_vec), function(ci) ks.test(vec[cl_vec == ci], vec[cl_vec != ci], alternative = alternative)$statistic, .parallel = parallel))
}

iemf <- ie_mat[sp_varp_atac_dist,fp_tf_motif]
iemfp <- ie_mat[sp_varp_atac_prom,fp_tf_motif]
res <- plyr::laply(fp_tf_motif, function(x) calculate_d_in_cluster(iemf[,x], names(sp_varp_atac_dist)), .parallel = TRUE)
rownames(res) <- fp_tf_motif
colnames(res) <- unique(names(sp_varp_atac_dist))
resp <- plyr::laply(fp_tf_motif, function(x) calculate_d_in_cluster(iemfp[,x], names(sp_varp_atac_prom)), .parallel = TRUE)
rownames(resp) <- fp_tf_motif
colnames(resp) <- unique(names(sp_varp_atac_prom))


motif_maxs <- setNames(rownames(res)[apply(res, 2, which.max)], colnames(res))
amd <- all_motif_datasets()
res_reg <- plyr::llply(names(motif_maxs), function(i) {
    prego::regress_pwm(sequences=seqs_all[match(sp_varp_atac_dist, rownames(mcl_all))], 
                            response = ifelse(names(sp_varp_atac_dist) == i, 1, 0), 
                            motif = dplyr::filter(amd, motif == motif_maxs[[i]]),
                                        unif_prior = 0.1,
                                        score_metric = "ks",
                                        final_metric = "ks",
                                        use_sge = F, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
}, .parallel = T)

res_reg_r2 <- plyr::llply(names(motif_maxs), function(i) {
    prego::regress_pwm(sequences=seqs_all[match(sp_varp_atac_dist, rownames(mcl_all))], 
                            response = ifelse(names(sp_varp_atac_dist) == i, 1, 0), 
                            motif = dplyr::filter(amd, motif == motif_maxs[[i]]),
                                        unif_prior = 0.1,
                                        score_metric = "r2",
                                        final_metric = "ks",
                                        use_sge = F, 
                                        match_with_db=T, use_sample = T, 
                                        parallel=T)
}, .parallel = T)

# res <- plyr::laply(fp_tf_motif, function(x) mean(ie_mat[,x]), .parallel = TRUE)
inds <- which(rowMaxs(res) >= 0.2)
ri <- res[inds,]
pd <- pheatmap(ri[order(apply(ri, 1, which.max), decreasing=F),], 
            fontsize_row = 6, cluster_cols = F, cluster_rows = F, silent = F)
indsp <- which(rowMaxs(res) >= 0.2)
rp <- resp[indsp,]
pp <- pheatmap(rp[order(apply(rp, 1, which.max), decreasing=T),], 
            fontsize_row = 10, cluster_cols = F, cluster_rows = F, silent = F)
ctd <- cutree(pd$tree_row, k = 12)
ctp <- cutree(pp$tree_row, k = 12)

dmot <- sapply(unique(ctd), function(x) fp_tf_motif[which(ctd == x)[which.max(rowMaxs(res)[ctd == x])]])
iemfn <- iemf - rowMaxs(as.matrix(iemf))
iemfn[iemfn < -10] <- -10

iemf2 <- ie_mat[sp_varp_atac_dist,]
res_all <- plyr::laply(1:ncol(iemf2), function(x) calculate_d_in_cluster(iemf2[,x], names(sp_varp_atac_dist)), .parallel = TRUE)
raq98 <- apply(iemf, 2, quantile, probs = 0.98)
ra_98_bin <- sapply(1:length(raq98), function(i) as.numeric(iemf[,i] >= raq98[[i]]))
colnames(ra_98_bin) <- colnames(iemf)
ra_98_sum_clust <- tgs_matrix_tapply(t(ra_98_bin), names(sp_varp_atac_dist), sum)
ra_98_lfc <- log2(0.1+ra_98_sum_clust/(0.02*as.numeric(table(names(sp_varp_atac_dist)))))
p_r9l <- pheatmap(ra_98_lfc, color = colorRampPalette(c('blue4', 'white', 'red4'))(100), show_colnames = F,
                breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)))
save_pheatmap(p_r9l, 'output/sequence_modeling/figs/all_peaks_ra_98_lfc.png', width = 2500, height = 2500, res = 200)
inds <- which(colMaxs(ra_98_lfc) - colMins(ra_98_lfc) >= 2)
cor_mot <- tgs_cor(ra_98_lfc[,inds], spearman = T)
p_cor_mot <- pheatmap(cor_mot, breaks = seq(-1,1,l=101), colorRampPalette(c('blue4', 'white', 'red4'))(100), show_colnames = F, show_rownames = F)
save_pheatmap(p_cor_mot, 'output/sequence_modeling/figs/ra_98_lfc_cor_motifs.png', width = 1800, height = 1800, res = 100)

hccm <- hclust(tgs_dist(cor_mot))
ctmot <- cutree(hccm, k = 20)
mot_max_ks <- colMaxs(ra_98_lfc) - colMins(ra_98_lfc)
names(mot_max_ks) <- rownames(res_all)
motifs_to_take <- sapply(sort(unique(ctmot)), function(k) {
                k_inds <- which(ctmot == k); 
                ctmk <- ctmot[k_inds]; 
                return(ctmk[which.max(mot_max_ks[names(ctmk)])])
                })
pltmt <- t(ra_98_lfc[order(match(rownames(ra_98_lfc), names(sp_varp_atac_dist))),names(motifs_to_take)])
p_plt <- pheatmap(pltmt, cluster_rows = T,cluster_cols = F, silent = F)
svat <- table(names(sp_varp_atac_dist))
pltmt_rep <- do.call('cbind', lapply(colnames(pltmt), function(i) matrix(rep(pltmt[,i], svat[[i]]), nrow = nrow(pltmt))))
rownames(pltmt_rep) <- rownames(pltmt)
colnames(pltmt_rep) <- sp_varp_atac_dist
pprep <- pheatmap(pltmt_rep[p_plt$tree_row$order,], annotation_col = pk_col_annot_atac, annotation_colors = ann_colors,
                color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                breaks = c(seq(-3,0,l=50), seq(0.01,3,l=51)), 
                show_colnames = F,
                cluster_rows = F,
                cluster_cols = F,
                fontsize_row = 8)
save_pheatmap(pprep, 'output/sequence_modeling/figs/inter_cluster_ra_98_lfc.png', width = 3400, height = 800, res = 200)

motif_plots <- lapply(rownames(pltmt_rep)[p_plt$tree_row$order], function(nm) prego::plot_pssm_logo(filter(amd, motif == nm), title = nm))
names(motif_plots) <- rownames(pltmt_rep)[p_plt$tree_row$order]
gall <- marrangeGrob(grobs = motif_plots, ncol =1, nrow = length(motif_plots))
ggsave('./output/sequence_modeling/figs/gall.png', gall, height = 40, width = 8, units = 'in')



piemf2 <- pheatmap(res_all[rowMaxs(res_all) >= 0.27,unique(names(sp_varp_atac_dist))], cluster_cols = F)
save_pheatmap(piemf2, 'output/sequence_modeling/figs/inter_cluster_ks_all_motifs_colmax_over_0_27.png', width = 1800, height = 800, res = 200)



int_select <- coords_exp[seq_coords$peak_name %in% sp_varp_atac_dist,]
rownames(int_select) <- seq_coords$peak_name[seq_coords$peak_name %in% sp_varp_atac_dist]
int_select2 <- int_select[order(match(rownames(int_select), sp_varp_atac_dist)),]
mex_ag <- rownames(res_all)[rowMaxs(res_all) >= 0.25]
# vals <- gextract_pwm(intervals = int_select2,motifs = mex_ag, prior = 0.01)
vals <- gextract_pwm(intervals = int_select2,motifs = names(motifs_to_take), prior = 0.01)
vm <- t(as.matrix(subset(vals, select = -c(chrom, start, end))))
vmn <- vm - rowMaxs(vm)
vmn[vmn < -10] <- -10

cor_vm <- tgs_cor(t(vm), spearman=T)
pcor_vm <- pheatmap(cor_vm,color = colorRampPalette(c('blue4', 'white', 'red4'))(100), fontsize_col = 14, fontsize_row =14,
                breaks = seq(-1,1,l=101))
save_pheatmap(pcor_vm, 'output/sequence_modeling/figs/cor_energies_motifs_to_take.png', width = 2800, height = 2800, res = 150)

gall2 <- marrangeGrob(grobs = motif_plots[order(match(names(motif_plots), rownames(cor_vm)[pcor_vm$tree_row$order]))], ncol =1, nrow = length(motif_plots))
ggsave('./output/sequence_modeling/figs/gall2.png', gall, height = 40, width = 8, units = 'in')


energy_cluster <- tgs_matrix_tapply(vm, names(sp_varp_atac_dist), mean)
accessibility_cluster <- tgs_matrix_tapply(t(mcl_all_norm[sp_varp_atac_dist,]), names(sp_varp_atac_dist), mean)
cor_energy_acc <- tgs_cor(energy_cluster,accessibility_cluster)
pltmt_cor_en_acc <- t(cor_energy_acc[,amc_ord])
p_cor_en_acc <- pheatmap(pltmt_cor_en_acc[,order(apply(1+pltmt_cor_en_acc, 2, function(x) sum(x*1:length(x))/sum(x)))], 
                show_rownames = F, 
                cluster_cols = F,
                cluster_rows = F, 
                color = colorRampPalette(c('blue4', 'white', 'red4'))(100), 
                fontsize_col = 10, 
                fontsize_row =10, 
                annotation_row = row_annot, 
                annotation_colors = ann_colors,
                breaks = seq(-1,1,l=101))
save_pheatmap(p_cor_en_acc, 'output/sequence_modeling/figs/cor_energies_accessibility.png', width = 2800, height = 2000, res = 200)

gall3 <- marrangeGrob(grobs = motif_plots[rev(order(match(names(motif_plots), colnames(pltmt_cor_en_acc)[order(apply(1+pltmt_cor_en_acc, 2, function(x) sum(x*1:length(x))/sum(x)))])))], ncol =1, nrow = length(motif_plots))
ggsave('./output/sequence_modeling/figs/gall3.png', gall3, height = 40, width = 8, units = 'in')


load('./output/mcatac/ann_colors.rda')
load('./output/mcatac/pk_col_annot_atac.rda')

pvm <- pheatmap(vmn, cluster_rows = T, cluster_cols = F, show_colnames=F, fontsize_row = 6, 
            annotation_col = pk_col_annot_atac, annotation_colors = ann_colors)
save_pheatmap(pvm, 'output/sequence_modeling/figs/select_motifs_energy_values_norm_max_0.png', width = 2800, height = 1200, res = 150)

mttt <- apply(vmn, 1, zoo::rollmean, k = 50)
pvmn2 <- pheatmap(t(mttt), cluster_cols = F, cluster_rows = T, annotation_col = pk_col_annot_atac, annotation_colors = ann_colors, show_colnames = F,fontsize_row = 10)
save_pheatmap(pvmn2, 'output/sequence_modeling/figs/select_motifs_energy_values_norm_max_0_rollmean_50.png', width = 2800, height = 1000, res = 150)


pheatmap(t(apply(t(iemfn[sp_varp_atac_dist,dmot]), 1, zoo::rollmean, k = 10)), cluster_cols = F, show_colnames = F)

ks_or_prego <- apply(cbind(colMaxs(ri), sapply(res_reg_r2, function(x) x$ks$stat)), 1, which.max)

best_pssms <- lapply(1:length(ks_or_prego), function(i) ifelse(ks_or_prego[[i]] == 1, 
                return(filter(amd, motif == rownames(ri)[which.max(ri[,i])])), 
                return(res_reg_r2[[i]]$pssm)))
names(best_pssms) <- colnames(ri)
best_pssms <- lapply(seq_along(best_pssms), function(i) {
                x <- best_pssms[[i]]; 
                x$motif <- names(best_pssms)[[i]]; 
                select(relocate(x, motif, .before = pos), c('motif', 'pos', 'A', 'C', 'G', 'T'))
                })

pssm_df <- do.call('rbind', best_pssms)
int_select <- coords_exp[seq_coords$peak_name %in% sp_varp_atac_dist,]
rownames(int_select) <- seq_coords$peak_name[seq_coords$peak_name %in% sp_varp_atac_dist]
int_select2 <- int_select[order(match(rownames(int_select), sp_varp_atac_dist)),]
rex_pssm <- gextract_pwm(intervals = int_select2,dataset = pssm_df, prior = 0.01)
rpmat <- as.matrix(subset(rex_pssm, select = -c(chrom, start, end)))

df1 <- filter(amd, motif %in% motif_maxs)
pc1 <- gextract_pwm(intervals = int_select2,dataset = df1, prior = 0.01)
df2 <- do.call('rbind', lapply(seq_along(res_reg_r2), function(i) {y <- res_reg_r2[[i]]$pssm; y$motif <- colnames(ri)[[i]]; relocate(y, motif, .before = pos)}))
pc2 <- gextract_pwm(intervals = int_select2,dataset = df2, prior = 0.01)


pc1mat <- as.matrix(subset(pc1, select = -c(chrom, start, end)))
pc1mat <- t(pc1mat) - colMaxs(pc1mat)
pc1mat[pc1mat < -10]  <- -10
quantile(rowMaxs(pc1mat))

pc2mat <- as.matrix(subset(pc2, select = -c(chrom, start, end)))
pc2mat <- t(pc2mat) - colMaxs(pc2mat)
pc2mat[pc2mat < -10]  <- -10
dev.new();pheatmap(pc1mat[,order(match(colnames(pc1mat), colnames(ri)))], cluster_cols = F, cluster_rows = F, show_colnames=F)
dev.new();pheatmap(pc2mat[,order(match(colnames(pc2mat), colnames(ri)))], cluster_cols = F, cluster_rows = F, show_colnames=F)
