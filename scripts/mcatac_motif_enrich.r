devtools::load_all("~/src/mcATAC/")
library(prego)
library(misha.ext)
library(matrixStats)
gset_genome('mm10')
prego::set_parallel(40)
options(gmax.data.size = 1e+9)
setwd('~/raid/proj/mmcortex/')
doMC::registerDoMC(70)
### Get peaks
# mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks.rds')

# ### Get cluster assignments
mca_rna_km <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peak_km_a_legc.rds')
aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
scatac_km <- aaa$km_a_legc

# ### Get sequences
# seq_coords <- mca@peaks[,c('chrom', 'start', 'end', 'peak_name')]
# # coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100)
# # seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# intervs_energy <- gextract_pwm(intervals = seq_coords, dataset = all_motif_datasets(), prior = 0.01)
# saveRDS(intervs_energy, './output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')
intervs_energy <- readRDS('./output/sequence_modeling/mmcortex_feat_peak_motif_energy.rds')


gen_random_genome_intervals <- function(num_peaks = 1e+5,
                                                peak_width = 2e+2,
                                                bp_from_chrom_edge_to_avoid = 3e+6,
                                                datasets_of_interest = NULL,
                                                motif_regex = NULL,
                                                motif_tracks = NULL,
                                                parallel = TRUE) {
    ALLGENOME[[1]] <- ALLGENOME[[1]][!grepl("_", ALLGENOME[[1]]$chrom), ]
    chrom_lens <- apply(ALLGENOME[[1]][, 2:3], 1, diff)
    chrom_fracs <- setNames(chrom_lens / sum(chrom_lens), ALLGENOME[[1]][, 1])
    sample_seqs <- mapply(chrom_fracs, names(chrom_fracs), chrom_lens, FUN = function(x, y, z) {
        chrom <- rep(y, round(x * num_peaks))
        start <- sample.int(n = z, size = length(chrom))
        end <- start + peak_width
        return(as.data.frame(rbind(chrom, start, end)))
    })
    sample_seqs <- as.data.frame(do.call("rbind", lapply(sample_seqs, t)))
    sample_seqs[, 2:3] <- apply(sample_seqs[, 2:3], 2, as.numeric)
    sample_seqs <- sample_seqs[with(sample_seqs, order(chrom, start, end)), ]
    end_shift <- ALLGENOME[[1]][match(sample_seqs$chrom, ALLGENOME[[1]][, 1]), 3] - bp_from_chrom_edge_to_avoid
    sample_seqs <- sample_seqs[sample_seqs$start >= bp_from_chrom_edge_to_avoid & sample_seqs$end <= end_shift, ]
    sample_seqs <- PeakIntervals(sample_seqs)
    return(sample_seqs)
}

# rg_ints <- gen_random_genome_intervals(num_peaks = 2.5e+4, peak_width = 300)
rg_ints_tmp <- gen_random_genome_intervals(num_peaks = 1e+4, peak_width = 300)

# rg_energy <- gextract_pwm(intervals = rg_ints, dataset = all_motif_datasets(), prior = 0.01)
# saveRDS(rg_energy, './output/sequence_modeling/random_genome_25k_motif_energy.rds')

rg_energy <- readRDS('./output/sequence_modeling/random_genome_25k_motif_energy.rds')

#' Calculate Kolmogorov-Smirnov D statistics between two interval sets with motif energies
#'
#' This function does a one-sided KS test between a foreground set of peaks (\code{pssm_fg}) and a background set \code{pssm_bg}.
#' The option \code{alternative == "less"}, checks the null hypothesis that the foreground distribution is not less than the
#' background distribution (applicable when looking for motif enrichment; for anti-enrichment, \code{alternative == 'greater'},
#' see ks.test documentation for further details)
#' @param pssm_fg motif energies calculated for a certain set of motifs on a PeakIntervals/ScATAC/McATAC object
#' @param pssm_bg a background set of intervals (e.g. random genome, all ENCODE enhancers etc.) that include all/subset of the motifs (columns) in pssm_fg
#' @param fg_clustering a vector of cluster assignments for the foreground peaks (e.g. from \code{gen_atac_peak_clust})
#' @param parallel (optional) - whether to use parallelize computations
#' @param nc (optional) - number of cores for parallel computations
#' @inheritParams stats::ks.test
#' @return if \code{fg_clustering == TRUE}, returns a matrix of clusters x motifs (rows x columns) with the D-statistic for each combination
#' @examples
#' \dontrun{
#' pssm_fg <- generate_motif_pssm_matrix(my_atac_mc, datasets_of_interest = "jaspar")
#' pssm_bg <- gen_random_genome_peak_motif_matrix(num_peaks = nrow(my_atac_mc@peaks), datasets_of_interest = "jaspar")
#' d_vs_rg <- calculate_d_stats(pssm_fg, pssm_bg)
#' peak_clust <- gen_atac_peak_clust(my_atac_mc, k = 12)
#' d_vs_rg_cl <- calculate_d_stats(pssm_fg, pssm_bg, fg_clustering = peak_clust)
#' }
#' @export
calculate_d_stats <- function(pssm_fg, pssm_bg, fg_clustering = NULL, parallel = getOption("mcatac.parallel"), alternative = "greater", nc = getOption("mcatac.parallel.nc")) {
    cols_fg <- grep("chrom|start|end$|interval|peak_name", colnames(pssm_fg), ignore.case = T, invert = T, value = T)
    cols_bg <- grep("chrom|start|end$|interval|peak_name", colnames(pssm_bg), ignore.case = T, invert = T, value = T)
    cols_both <- intersect(cols_fg, cols_bg)
    if (!parallel) {
        nc <- pmax(2, round(0.1 * nc))
    }
    if (!is.null(fg_clustering)) {
        ks_test_results <- parallel::mclapply(cols_both, FUN = function(x, i) {
            return(tapply(x[, i], fg_clustering, function(y) suppressWarnings(ks.test(y, pssm_bg[, i], alternative = alternative))))
        }, mc.cores = nc, x = pssm_fg)
        ks_d <- sapply(ks_test_results, function(x) sapply(x, function(y) y$statistic))
        colnames(ks_d) <- cols_both
    } else {
        ks_test_results <- parallel::mclapply(cols_both, function(x) {
            ks.test(x = pssm_fg[, x], y = pssm_bg[, x], alternative = alternative)
        },
        mc.cores = nc
        )
        ks_d <- setNames(sapply(ks_test_results, function(x) x$statistic), cols_both)
    }
    return(ks_d)
}

d_atac <- calculate_d_stats(pssm_fg = intervs_energy, pssm_bg = rg_energy, fg_clustering = scatac_km$cluster, alternative = 'greater')
d_rna <- calculate_d_stats(pssm_fg = intervs_energy, pssm_bg = rg_energy, fg_clustering = mca_rna_km$cluster, alternative = 'greater')

png('./output/sequence_modeling/figs/d_stats_atac.png', width = 14, height = 6, units = 'in', res = 300)
pheatmap(d_atac, fontsize_col = 2)
dev.off()

png('./output/sequence_modeling/figs/d_stats_rna.png', width = 14, height = 6, units = 'in', res = 300)
pheatmap(d_rna, fontsize_col = 2)
dev.off()

png('./output/sequence_modeling/figs/qq_plot_atac_rna_cluster_Ds.png', width = 6, height = 6, units = 'in', res = 300)
plot(quantile(d_atac, (0:100)/100), quantile(d_rna, (0:100)/100), xlab = 'D-statistic quantiles in ATAC', ylab = 'D-statistic quantiles in  RNA')
abline(0,1,lty=2,col='red')
dev.off()


png('./output/sequence_modeling/figs/atac_vs_rna_d_maxs.png', width = 7, height = 7, units = 'in', res = 300)
plot(colMaxs(d_atac), colMaxs(d_rna), main = 'maximal D stats across clusters', ylab = 'Max D-statistic (RNA)', xlab = 'Max D-statistic (ATAC)', pch = 16, cex = 0.5)
abline(0, 1, col = 'red',lty=2)
dev.off()

png('./output/sequence_modeling/figs/atac_vs_rna_d_means.png', width = 7, height = 7, units = 'in', res = 300)
plot(colMeans(d_atac), colMeans(d_rna), main = 'mean D stats across clusters', ylab = 'mean D-statistic (RNA)', xlab = 'mean D-statistic (ATAC)', pch = 16, cex = 0.5)
abline(0, 1, col = 'red',lty=2)
dev.off()

library(matrixStats)

library(metacell)
scdb_init(base_dir = './scdb', force_reinit = T)
mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
mc <- scdb_mc('pl_cort')
nsc_mc <- which(mcmd$cell_type == 'NSC')
legc <- log2(1e-5+ mc@e_gc)
inds_d <- which(colMaxs(d_atac) >= 0.4)

tfs = readr::read_tsv('~/raid/Mus_musculus_TF.txt')
tfs <- c('Bcl11b', tfs$Symbol)
tfs <- tfs[tfs %in% rownames(mc@e_gc)]
nsc_tfs <- tfs[tfs %in% rownames(mc@mc_fp)[rowMaxs(mc@mc_fp[,nsc_mc]) > 3]]
fp_tfs <- tfs[tfs %in% rownames(mc@mc_fp)[rowMaxs(mc@mc_fp) >= 3]]
# nsc_tf_motif <- unlist(plyr::llply(nsc_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(d_atac)[inds_d], v=T, ign=T), .parallel = T))
# nsc_tf_motif <- c(nsc_tf_motif, unlist(plyr::llply(c('emx', 'pax', 'dmrt'), function(g) grep(g, colnames(d_atac), v=T, ign=T), .parallel = T)))
nsc_tf_motif <- unlist(plyr::llply(fp_tfs, function(g) grep(paste0(g, '$|', g, '_'), colnames(d_atac)[inds_d], v=T, ign=T), .parallel = T))

load('~/raid/proj/mmcortex/output/mcatac/peak_indices_var_atac_clust.rda')
nsc_motif_pwm <- t(as.matrix(intervs_energy[sp_varp_atac,nsc_tf_motif]))
colnames(nsc_motif_pwm) <- misha.ext::convert_misha_intervals_to_10x_peak_names(intervs_energy[,1:3])[sp_varp_atac]
count_quantile_hits_in_clusters <- function(mat, cl, q_thresh = 0.98) {
    mat_q <- plyr::aaply(mat, 1, function(x) ecdf(x)(x), .parallel = T)
    mat_count_q_thresh <- tgs_matrix_tapply(mat_q, cl, function(x) sum(x >= q_thresh)/(length(x)*(1-q_thresh)))    
    return(mat_count_q_thresh)
}
nsc_motif_count_q_98 <- count_quantile_hits_in_clusters(nsc_motif_pwm, names(sp_varp_atac))
## Plot stacked barplot matrix of motif enrichment
sva_ht <- as.matrix(table(names(sp_varp_atac))/sum(table(names(sp_varp_atac))))
sva_mat <- matrix(rep(sva_ht, ncol(nsc_motif_count_q_98)), 
                nrow = ncol(nsc_motif_count_q_98), 
                ncol = length(sva_ht),
                byrow=T)
rownames(sva_mat) <- colnames(nsc_motif_count_q_98)
colnames(sva_mat) <- unique(names(sp_varp_atac))
cut_vec <- cut(nsc_motif_count_q_98, breaks = brks)
nsc_cut <- matrix(match(cut_vec, levels(cut_vec)),nrow = nrow(nsc_motif_count_q_98), ncol = ncol(nsc_motif_count_q_98))

nsc_motif_pwm_norm <- nsc_motif_pwm
nsc_motif_pwm_norm <- nsc_motif_pwm_norm - rowMaxs(nsc_motif_pwm_norm)
nsc_motif_pwm_norm[nsc_motif_pwm_norm < -15] <- -15
nsc_motif_lq <- plyr::aaply(nsc_motif_pwm, 1, function(x) {
    vec <- -log2(1 - ecdf(x)(x)); 
    vec[is.infinite(vec)] <- max(vec[!is.infinite(vec)]); 
    return(vec)}, 
    .parallel=T)
colnames(nsc_motif_lq) <- colnames(nsc_motif_pwm)




mcmd <- readr::read_tsv('~/raid/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv')
ann_colors = list('cell_type' = setNames(unique(mcmd$color), unique(mcmd$cell_type)))

pk_col_annot_atac <- as.data.frame(tibble::enframe(sp_varp_atac))
rownames(pk_col_annot_atac) <- colnames(nsc_motif_pwm)
pk_col_annot_atac <- dplyr::select(pk_col_annot_atac, name)
colnames(pk_col_annot_atac) <- 'cluster_atac'
ann_colors[['cluster_atac']] <- setNames(gplots::col2hex(sample(grep('white|gray|grey', colors(), inv=T, v=T), 
                                    length(unique(pk_col_annot_atac$cluster_atac)))), 
                                    unique(pk_col_annot_atac$cluster_atac))

p_pwm_norm <- pheatmap::pheatmap(nsc_motif_pwm_norm, cluster_cols = F, show_colnames=F,fontsize_row = 6,
                annotation_col = pk_col_annot_atac, annotation_colors=ann_colors)


p_pwm_lq <- pheatmap::pheatmap(nsc_motif_lq, cluster_cols = F, show_colnames=F,
                    main = '-log2(1-quantile(PWM))',
                     annotation_col = pk_col_annot_atac, annotation_colors=ann_colors)



tf_motif <- plyr::llply(tfs, function(g) grep(g, colnames(d_atac)[inds_d]), .parallel = T)
atac_d_filt <- d_atac[,colnames(d_atac)[inds_d][unique(unlist(tf_motif))]]
cor_mc_motif_2 <- tgs_cor(egc_cl, atac_d_filt, spearman =T)

png('./output/sequence_modeling/figs/correlation_of_d_with_mcatac.png', width = 10, height = 10, units = 'in', res = 300)
pheatmap(cor_mc_motif_2[cust_mc_ord_st,], cluster_rows = F, 
            annotation_row = col_annot, 
            annotation_colors = ann_colors, 
            color = rev(colorRampPalette(c('red', 'yellow', 'black'))(1000)), breaks = seq(0.5,1,l=1001))
dev.off()

