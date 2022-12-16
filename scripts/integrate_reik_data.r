### Integrate reik data
devtools::load_all("~/src/mcATAC/")
library(pheatmap)
gset_genome('mm10')
load( "/net/mraid14/export/tgdata/users/atanay/proj/enhflow/scdb/rna_md.Rmd")
ct_neuro <- c('Epiblast', 'Rostral neural plate','Neural crest', 'Definitive ectoderm', 'Forebrain/Midbrain/Hindbrain')
ct_neuro_mc <- lapply(ct_neuro, function(ct) {y <- md$metacell[md$cell_type == ct]; as.numeric(levels(y)[y])})
names(ct_neuro_mc) <- ct_neuro
rk_trk <- grep('marginal|DELETE|smoothed|R1|R2|grouped', gtrack.ls('wt_reik'), inv=T, v=T)
rk_trk <- rk_trk[order(as.numeric(gsub('wt_reik\\.mc', '', rk_trk)))]
aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
peaks <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(aaa$mcl_all))
reik_all_peaks <- gextract(rk_trk[unlist(ct_neuro_mc)], peaks, iterator = peaks)

reik_neuro_peaks <- subset(reik_all_peaks, select = -c(chrom, start, end, intervalID))
rownames(reik_neuro_peaks) <- misha.ext::convert_misha_intervals_to_10x_peak_names(reik_all_peaks[,1:3])
reik_neuro_peaks[is.na(reik_neuro_peaks)] <- 0

load('./output/mcatac/peak_indices_var_atac_clust.rda')
rnn <- t(t(reik_neuro_peaks)/colSums(reik_neuro_peaks))
rnnl <- log2(1e-7 + rnn)
colnames(rnnl) <- as.numeric(gsub('wt_reik\\.mc', '', colnames(rnnl)))
color_key <- unique(md[,c('cell_type', 'color')])
row_annot_reik <- tibble::column_to_rownames(md[,c('metacell','cell_type')], 'metacell')
row_annot_reik <- as.data.frame(as.character(row_annot_reik[order(as.numeric(rownames(row_annot_reik))),]))
colnames(row_annot_reik)  <- 'cell_type'
rownames(row_annot_reik)  <- as.numeric(rownames(row_annot_reik))


# pheatmap(t(rnnl[sp_varp_atac,]), 
#                 cluster_rows = F, cluster_cols = F, 
#                 show_colnames=F, show_rownames=F, 
#                 annotation_row = row_annot_reik, annotation_colors=ann_colors)

rnnl_sum_ct <- tgs_matrix_tapply(t(rnnl[sp_varp_atac,]), names(sp_varp_atac), mean)
mca_svad <- readRDS('./output/mcatac/mcl_all_norm_sp_varp_atac_dist.rds')
load('./output/mcatac/atac_metacell_order_annotation.rda')
ann_colors[['cell_type']] <- c(ann_colors[['cell_type']], tibble::deframe(color_key[color_key$cell_type %in% ct_neuro,]))
reik_mmc_merge <- t(as.matrix(bind_cols(rnnl[sp_varp_atac_dist,], as.matrix(log2(1e-5+t(mca_svad[amc_ord,]))))))
colnames(reik_mmc_merge) <- rownames(aaa$mcl_all[sp_varp_atac_dist,])
ra_all <- bind_rows(row_annot_reik, row_annot)
prmm <- pheatmap(t(reik_mmc_merge[sp_varp_atac_dist,]), annotation_row = ra_all, annotation_colors = ac_filt,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames=F, show_colnames=F,
                 annotation_col = pk_col_annot_atac)
save_pheatmap(prmm, './output/mcatac/figs/var_atac_peaks_mmcortex_and_reik.png', height = 2000, width = 3500, res = 150)

reik_mmc_merge_sum_cl <- tgs_matrix_tapply(reik_mmc_merge, names(sp_varp_atac_dist), mean)
p_reik_avg_cl <- pheatmap(t(reik_mmc_merge_sum_cl), cluster_rows = F, cluster_cols =F, annotation_row = ra_all, annotation_colors = ann_colors)
save_pheatmap(p_reik_avg_cl, './output/mcatac/figs/var_atac_peaks_avg_clust_mmcortex_and_reik.png', height = 1200, width = 2000, res = 150)

## ENCODE early mouse brain ATAC data
misha::gtrack.import('encode_e11_mouse_brain_atac', 
                    description='ENCFF326ULQ.bigWig, ENCODE mouse brain E11.5 ATAC-seq', 
                    file = './encode_forebrain_neural_tube_e10_e14_data/E11/ENCFF326ULQ.bigWig',
                    binsize = 1)

misha::gtrack.import('encode_e12_mouse_brain_atac', 
                    description='ENCFF541FKK.bigWig, ENCODE mouse brain E12.5 ATAC-seq', 
                    file = './encode_forebrain_neural_tube_e10_e14_data/E12/ENCFF541FKK.bigWig',
                    binsize = 1)

misha::gtrack.import('encode_e13_mouse_brain_atac', 
                    description='ENCFF633ETU.bigWig, ENCODE mouse brain E13.5 ATAC-seq', 
                    file = './encode_forebrain_neural_tube_e10_e14_data/E13/ENCFF633ETU.bigWig',
                    binsize = 1)


enc_peaks_list <- list(enc11 = enc11, enc12 = enc12, enc13 = enc13)
save(enc_peaks_list, file = './output/mcatac/enc_peaks_list.rda')

reik_mmc_enc_merge <- t(as.matrix(bind_cols(rnnl[sp_varp_atac_dist,], atac_enc_mmc_peaks_ln[sp_varp_atac_dist,], as.matrix(log2(1e-5+t(mca_svad[amc_ord,]))))))
colnames(reik_mmc_enc_merge) <- rownames(aaa$mcl_all[sp_varp_atac_dist,])
enc_row_annot <- as.data.frame(list('cell_type' = c('ENCODE_E11', 'ENCODE_E12','ENCODE_E13')))
rarf <- dplyr::filter(row_annot_reik, cell_type %in% ct_neuro) %>% mutate(dummy <- 1:nrow(.))
rarf <- rarf[order(match(rownames(rarf), colnames(rnnl))),]
hcmcs <- hclust(tgs_dist(as.matrix(log2(1e-5+mca_svad)))
ct_hc <- cutree(hcmcs, k = 8)
ra_all <- bind_rows(dplyr::select(rarf, cell_type), enc_row_annot, tibble::column_to_rownames(tibble::enframe(as.character(ct_hc), name = 'amc', value = 'cell_type'), 'amc'))
all_data_avg_ct <- tgs_matrix_tapply(t(reik_mmc_enc_merge), ra_all$cell_type, mean)
all_data_avg_ct <- all_data_avg_ct[c(ct_neuro, enc_row_annot$cell_type, 1:8),]
p_adac <- pheatmap(all_data_avg_ct, cluster_cols = F, show_colnames = F, cluster_rows = F)
save_pheatmap(p_adac, './output/mcatac/figs/var_atac_peaks_avg_clust_reik_encode_mmcortex.png', height = 1500, width = 2500, res = 150)

seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(aaa$mcl_all))
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))


## GC content
nuc_content <- gextract('(seqG + seqC)/2', intervals = coords_exp, iterator = coords_exp, colnames = 'gc_content')
nuc_content$peak_name <- peak_names(nuc_content, tad_based=F)
nuc_content$cluster <- aaa$km_a_legc$cluster
cg_check <- gscreen("seq.CG == 1", intervals = coords_exp, iterator = 1)
nei_coords_exp_cg_check <- gintervals.neighbors(cg_check, coords_exp, maxdist = 0, mindist = 0)
peak_cpg <- table(nei_coords_exp_cg_check$peak_name)
coords_exp$peak_cpg <- peak_cpg[coords_exp$peak_name]
coords_exp$peak_cpg[is.na(coords_exp$peak_cpg)] <- 0
coords_exp$cpg_over_median <- coords_exp$peak_cpg/median(coords_exp$peak_cpg)

## Early-late
library(misha)
library(misha.ext)
gset_genome("mm10")
gvtrack.create("tor", "Encode.esd3.replichip.rep1", "avg")
gvtrack.iterator("tor", sshift = -15000, eshift = 15000)
late_regs <- gscreen("tor <= 0")
early_regs <- gscreen("tor > 0")
early_regs <- dplyr::mutate(early_regs, early_name = peak_names(early_regs, tad_based=F))
late_regs <- dplyr::mutate(late_regs, late_name = peak_names(late_regs, tad_based=F))
nei_early <- gintervals.neighbors(coords_exp, early_regs, mindist=0,maxdist=0)
nei_late <- gintervals.neighbors(coords_exp, late_regs, mindist=0,maxdist=0)
coords_exp$is_early <- ifelse(coords_exp$peak_name %in% nei_early$peak_name, 1, 0)
ce_acoords_exp$is_late <- ifelse(coords_exp$peak_name %in% nei_late$peak_name, 1, 0)

coords_exp$gc_content <- nuc_content$gc_content[match(misha.ext::convert_misha_intervals_to_10x_peak_names(coords_exp), 
                                            nuc_content$peak_name)]
coords_exp$gc_over_median <- coords_exp$gc_content/median(coords_exp$gc_content)

ce_avg_clust <- as.data.frame(tgs_matrix_tapply(t(coords_exp[which(rownames(coords_exp) %in% sp_varp_atac_dist),
                        c('gc_content', 'peak_cpg', 'cpg_over_median', 'gc_over_median', 'is_early', 'is_late')]), 
                        aaa$km_a_legc$cluster[which(rownames(coords_exp) %in% sp_varp_atac_dist)], mean))
ce_avg_clust$early_over_median <- ce_avg_clust$is_early/median(ce_avg_clust$is_early)
ce_avg_clust$late_over_median <- ce_avg_clust$is_late/median(ce_avg_clust$is_late)

mca_svad <- readRDS('./output/mcatac/mcl_all_norm_sp_varp_atac_dist.rds')
sp_varp_atac_dist <- colnames(mca_svad)
names(sp_varp_atac_dist) <- names(sp_varp_atac)[match(sp_varp_atac_dist, sp_varp_atac)]


hcce <- hclust(tgs_dist(ce_avg_clust[as.numeric(unique(names(sp_varp_atac))),]), 
                method = 'ward.D2')
png('./output/mcatac/figs/peak_cluster_cg_cpg_e_l_log2_over_median.png', height = 1500, w = 600, res=  150)
par(mfrow = c(4,1), las = 2, mar = c(6,5,4,1), cex.lab = 1.5)
purrr::walk(c('cpg_over_median', 'gc_over_median', 'early_over_median', 'late_over_median'), 
                function(x) barplot(log2(ce_avg_clust[,x]), 
                                names.arg = rownames(ce_avg_clust), ylab = paste0('log2 ', x)))
dev.off()

png('./output/mcatac/figs/peak_cluster_cg_cpg_e_l.png', height = 1500, w = 600, res=  150)
par(mfrow = c(4,1), las = 2, mar = c(6,5,4,1), cex.lab = 1.5)
purrr::walk(c('peak_cpg', 'gc_content', 'is_early', 'is_late'), 
                function(x) barplot(ce_avg_clust[,x], 
                                names.arg = rownames(ce_avg_clust), ylab = x))
dev.off()



