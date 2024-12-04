SEED = 1337
set.seed(SEED)
wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
library(metacell)
scdb_init('scdb', f=T)
devtools::load_all("~/src/mcATAC/")
gset_genome("mm10")
gdb.reload()

write_sc_counts_from_fragments(fragments_file = './scatac_data/fragments_filtered.tsv.gz', 
            out_dir = './data/frag_reads', 
            cell_names = colnames(mat@mat), 
            genome = "mm10")

# metacells_from_mcl_flow <- function(flow_path, day_mcl_path) {
#     flow_file <- readRDS(flow_path)
#     day_mcls <- readRDS(day_mcl_path)
#     flow_mat <- flow_file$flow_mat
#     dmcl_vec <- unlist(lapply(seq_along(day_mcls), function(x,n,i) setNames(paste0(n[[i]], '_', x[[i]]$cor_km$cluster), colnames(x[[i]]$prom_mat)), x = day_mcls, n = names(day_mcls)))
#     mcls_nums <- unlist(lapply(names(day_mcls), function(d) paste0(d, "_", sort(unique(day_mcls[[1]]$cor_km$cluster)))))
#     dmcl_rename <- match(dmcl_vec, mcls_nums)
#     names(dmcl_rename) <- names(dmcl_vec)
#     flow_mat_norm <- flow_mat/rowSums(flow_mat)
#     flow_mat_cs <- apply(flow_mat_norm, 1, function(x) {
#         y <- setNames(x[which(x > 0)], which(x > 0));
#         ycs <- cumsum(y)
#         names(ycs) <- names(y)
#         return(ycs)
#     })
#     dmcl_rand_nums <- runif(length(dmcl_rename))
#     cell_to_metacell <- as.numeric(unlist(purrr::map2(.x = dmcl_rand_nums, .y = dmcl_rename, .f = function(.x,.y) {
#         names(flow_mat_cs[[.y]])[min(which(flow_mat_cs[[.y]] >= .x))]
#     })))
#     names(cell_to_metacell) <- names(dmcl_vec)
#     return(tibble::enframe(cell_to_metacell, name = 'cell', value = 'metacell'))
# }

# c2mc <- metacells_from_mcl_flow(flow_path = './data/pl_flow_res.rds', day_mcl_path = './data/pl_cort_prom_day_mcls.rds')
c2mc <- readr::read_tsv("./output/mcatac/c2mc.tsv")

flow_res_path <- file.path(wd, "output/mcatac/pl_cort_flow_mat.tsv")
flow_res_results <- readr::read_tsv(flow_res_path)
frr_mat <- as.matrix(dplyr::select(flow_res_results, -rowname))
rownames(frr_mat) <- flow_res_results$rowname
glia_ct <- c('Astrocytes', 'OPCs')
mc_glia <- which(mcmd$cell_type %in% glia_ct)
nnz_frcs <- apply(frr_mat[c('16_16', '16_20'),], 1, function(x) {sum_glia <- sum(x[mc_glia]); x[mc_glia] <- 0; nnz_mc <- which(x > 0); x_nnz_frc <- x[nnz_mc]/sum(x[nnz_mc]); return(setNames(x_nnz_frc, nnz_mc))})
nnz_cumsum <- lapply(nnz_frcs, function(x) cumsum(sort(x)))
flow_by_ct <- t(tgs_matrix_tapply(frr_mat, mcmd$cell_type, sum))
flow_by_ct_norm <- flow_by_ct/rowSums(flow_by_ct)

famc_new <- do.call('rbind', lapply(c(602, 603), function(x) {
    scah <- c2mc$cell[c2mc$metacell == x]; 
    amch <- rownames(frr_mat)[which(frr_mat[,x] > 0)];
    sca_rn <- runif(n = length(scah))
    famc_new <- as.numeric(sapply(seq_along(scah), function(i) mc_new <- names(nnz_cumsum[[amch]][nnz_cumsum[[amch]] >= sca_rn[[i]]])[[1]]))
    # print(cbind(sca_rn, famc_new))
    return(tibble::enframe(setNames(famc_new, scah), name = 'cell', value = 'metacell'))
}))

c2mc_new <- c2mc
c2mc_new$metacell[match(famc_new$cell, c2mc_new$cell)] <- famc_new$metacell


scc <- scc_read(path = './data/frag_reads/')
# mcc <- scc_to_mcc(sc_counts = scc, cell_to_metacell = dplyr::rename(c2mc, cell_id = cell))
mcc <- scc_to_mcc(sc_counts = scc, cell_to_metacell = dplyr::filter(dplyr::select(dplyr::rename(c2mc_new, cell_id = cell), cell_id, metacell), c2mc_new$cell_id %in% scc@cell_names))
mcc@metadata <- mcmd
mcc@rna_egc <- mc@e_gc
mcc_write(object=mcc, out_dir='./output/mcatac/mmcortex_mcc', overwrite = T)
mcc <- mcc_read('./output/mcatac/mmcortex_mcc/')
mcc_to_tracks(mc_counts = mcc, track_prefix = "mmcortex", overwrite = T, create_marginal_track = T)

# pks_95 <- call_peaks(marginal_track = 'mmcortex.marginal', 
#                     quantile_thresh = 0.95, 
#                     min_umis = 100, 
#                     split_peaks = T, 
#                     target_size = 300, 
#                     max_peak_size = 1e+3, 
#                     very_long = 2.4e+3, 
#                     min_peak_size = 100,  
#                     seed = 1337)
# tracks <- gtrack.ls('mmcortex\\.mc')

# tracks_prc <- gsub("mmcortex", "mmcortex_prc", tracks)
# mapply(tracks_prc, tracks, FUN = function(x,y) {
#     misha.ext::gtrack.create_dirs(x)
#     gvtrack.create(x, y, "global.percentile")
# })

# prc <- 0.975
# screen_peaks <- parallel::mclapply(tracks_prc, function(x) {
#     df <- gscreen(expr = glue::glue("{x} > {prc}"), intervals = gintervals.all(), iterator = 10)
#     df$mc <- gsub("mmcortex_prc.mc", "", x)
#     return(df)
# }, mc.cores = 60)

# sp_centers <- lapply(screen_peaks, function(x) {y <- misha.ext::gintervals.centers(x); y$value <- 1; y$mc <- x$mc; return(y)})

# sp_center_all <- do.call('rbind', sp_centers)

# sp_center_all$peak_name <- peak_names(sp_center_all, tad_based = FALSE)
# sp_center_all$peak_name <- gsub(":|-", "_", sp_center_all$peak_name)

# peak_sums <- tapply(sp_center_all$value, sp_center_all$peak_name, sum)

# sp_center_unique <- tibble::enframe(peak_sums, name = "peak_name", value = "sum")

# sp_center_intervs <- misha.ext::mat_to_intervs(tibble::column_to_rownames(sp_center_unique, "peak_name"))

# gdb.reload()
# gtrack.rm('mmcortex_peak_center_marginal', f=T)
# gdb.reload()

# gtrack.create_sparse(track = "mmcortex_peak_center_marginal", 
#                     description = "Count of peak centers across metacells from mmcortex scATAC data", 
#                     intervals = sp_center_intervs, 
#                     values = sp_center_intervs$sum)

# gdb.reload()
# gtrack.rm('mmcortex_peak_center_marginal_dense', f=T)
# gdb.reload()

# gtrack.create(track = "mmcortex_peak_center_marginal_dense", 
#                 expr = "mmcortex_peak_center_marginal", 
#                 description = "Dense count of peak centers across metacells from mmcortex scATAC data", 
#                 iterator = 1)

# gdb.reload()
# gtrack.rm('mmcortex_peak_center_marginal_smooth', f=T)
# gdb.reload()

# gtrack.modify(track = "mmcortex_peak_center_marginal_dense", 
#                 expr = "ifelse(is.na(mmcortex_peak_center_marginal_dense), 0, mmcortex_peak_center_marginal_dense)", 
#                 intervals = ALLGENOME[[1]])

# gtrack.smooth(track = "mmcortex_peak_center_marginal_smooth", 
#                 expr = "mmcortex_peak_center_marginal_dense", 
#                 winsize = 100, 
#                 iterator = 10,
#                 description = "Smoothed count of peak centers across metacells from mmcortex 10X scATAC data")
# sp <- gscreen("-log2(1 - mmcortex_peak_center_marginal_smooth) >= 0.1", iterator = 10)
# pl_ig_cort <- scdb_mat('pl_ig_cort')
# rn <- rownames(pl_ig_cort@mat)
# pks_ig_cort <- misha.ext::convert_10x_peak_names_to_misha_intervals(rn)
# nei_sp_pks_ig <- gintervals.neighbors(sp, pks_ig_cort, 
#                                         mindist = 0, 
#                                         maxdist = 0, 
#                                         maxneighbors = 10)
# sp$len <- sp$end - sp$start
# sp_filt <- dplyr::anti_join(sp, unique(nei_sp_pks_ig[,1:3]), by = c('chrom', 'start', 'end')) %>% 
#                 misha.ext::gintervals.centers(.) %>% 
#                 dplyr::mutate(start = start - 150, end = end + 149, len = end - start)
# spf_can <- gintervals.canonic(sp_filt) %>% mutate(len = end - start)
# np <- split_peaks_arbitrarily(dplyr::filter(spf_can, len >= 1e+3), max_peak_size = 300)
# pks_to_add <- dplyr::bind_rows(dplyr::filter(spf_can, len < 1e+3), np) %>% 
#                     arrange(chrom, start, end) %>% 
#                     mutate(len = end - start) %>% 
#                     PeakIntervals
# pks_all <- dplyr::bind_rows(dplyr::select(pks_95, chrom, start, end), 
#                             dplyr::select(pks_to_add, chrom, start, end)) %>% 
#                             dplyr::mutate(len = end - start) %>% 
#                             dplyr::arrange(chrom, start, end)
# pks_all <- unique(pks_all[,c("chrom", "start", "end")])
# pks_all$peak_name <- peak_names(pks_all)

# saveRDS(pks_all, './data/all_peaks_for_mmcortex_atac_mc.rds')

pks_all <- readRDS('./data/all_peaks_for_mmcortex_atac_mc.rds')

atac_mc <- mcc_to_mcatac(mc_counts = mcc, 
                            peaks = pks_all, 
                            metadata = dplyr::rename(mcmd, metacell = mc, cell_type = st))

mc_rna <- scdb_mc('pl_cort')

atac_mc <- add_mc_rna(atac_mc, mc_rna)

saveRDS(atac_mc, './data/mmcortex_atac_mc.rds')