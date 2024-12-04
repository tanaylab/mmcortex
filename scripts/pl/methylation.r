wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)

devtools::load_all("~/src/mcATAC/")
my_genome <- "mm10"
gset_genome(my_genome)

library(metacell)
scdb_init('scdb', f=T)
devtools::load_all("~/src/mcATAC/")
library(matrixStats)
library(pheatmap)

options(gmax.data.size = 1e+9)

mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
mcmd <- mcmd[!(mcmd$metacell %in% 602:603),]
cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]], 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))

cts <- c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3')

color_key <- unique(mcmd[,c('cell_type', 'color')])
color_key <- color_key[match(cust_st_ord, color_key$cell_type),]
col_key <- tibble::deframe(color_key)
col_annot <- as.data.frame(mcmd[,c('cell_type', 'mean_day')])
ann_colors <- list(cell_type = tibble::deframe(color_key), mean_day = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green3', 'blue', 'purple'))(100),
                                                                              seq(13,18,100)))

source('./scripts/util.r')

cpg_folder <- normalizePath('./cpg_methylation')

cpg_files <- grep('rep|mpra|dev|enh', list.files(cpg_folder), inv=T, v=T)
# cpg_files <- grep('rep', list.files(cpg_folder), inv=F, v=T)
cpg_files

cpg_filepaths <- paste0(cpg_folder, '/', cpg_files)
cpg_filepaths

# track_names <- gsub('\\.bedGraph.gz', '', basename(cpg_filepaths))
track_names <- gsub('\\.bedGraph', '', basename(cpg_filepaths))
track_names

# create_meth_tracks <- function(file_path, track_name) {
#     df <- vroom::vroom(file = file_path, delim = '\t', skip = 1, col_names = c('chrom', 'start', 'end', 'avg1', 'meth', 'unmeth'))
#     mdf <- dplyr::mutate(df, cov = meth + unmeth, avg = round(meth/(meth + unmeth), 3))
#     dir.create(paste0('/home/aviezerl/mm10/tracks/mmcortex/', track_name))
#     gdb.reload()
#     gtrack.create_sparse(track = paste0('mmcortex.', track_name, '.cov'), description = 'Neural stem cells CpG methylation data from Bonev lab', intervals =  dplyr::select(mdf, chrom, start, end), values = dplyr::pull(mdf, cov))
#     gtrack.create_sparse(track = paste0('mmcortex.', track_name, '.avg'), description = 'Neural stem cells CpG methylation data from Bonev lab', intervals =  dplyr::select(mdf, chrom, start, end), values = dplyr::pull(mdf, avg))
#     gtrack.create_sparse(track = paste0('mmcortex.', track_name, '.meth'), description = 'Neural stem cells CpG methylation data from Bonev lab', intervals =  dplyr::select(mdf, chrom, start, end), values = dplyr::pull(mdf, meth))
#     gtrack.create_sparse(track = paste0('mmcortex.', track_name, '.unmeth'), description = 'Neural stem cells CpG methylation data from Bonev lab', intervals =  dplyr::select(mdf, chrom, start, end), values = dplyr::pull(mdf, unmeth))
# }


# ttt <- sapply(tail(seq_along(track_names), -1), function(i) create_meth_tracks(cpg_filepaths[[i]], track_names[[i]]))

meth_tracks <- grep('tracks', gtrack.ls('NSC_meth_CpG.avg'), inv= T, v=T)

meth_cov_tracks <- grep('tracks', gtrack.ls('NSC_meth_CpG.cov'), inv= T, v=T)

mca <- readRDS(file.path(wd, 'output/mcatac/mmcortex_mcatac_feat_peaks_no_602_603.rds'))

a_legc <- log2(1e-5 + t(t(mca@egc)/colSums(mca@egc)))

nsc_peaks_vs_ipc <- get_genes_specific_to_mcs(a_legc, mc_pos = as.character(mcmd$metacell[which(mcmd$cell_type == 'NSC')]), mc_neg = as.character(mcmd$metacell[which(mcmd$cell_type %in% c('IPC', 'IPC_cyc'))]))
nsc_peaks_vs_astro <- get_genes_specific_to_mcs(a_legc, mc_pos = as.character(mcmd$metacell[which(mcmd$cell_type == 'NSC')]), mc_neg = as.character(mcmd$metacell[which(mcmd$cell_type %in% c('Astrocytes'))]))

ipc_peaks_vs_neuro <- get_genes_specific_to_mcs(a_legc, mc_pos = as.character(mcmd$metacell[which(mcmd$cell_type == 'IPC')]), 
                                                mc_neg = as.character(mcmd$metacell[which(mcmd$cell_type %in% c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3'))]))

nsc_peaks <- union(names(nsc_peaks_vs_ipc[nsc_peaks_vs_ipc >= 0.75]), names(nsc_peaks_vs_astro[nsc_peaks_vs_astro >= 0.75]))
length(nsc_peaks)

ipc_peaks <- union(names(nsc_peaks_vs_ipc[nsc_peaks_vs_ipc <= -0.75]), names(ipc_peaks_vs_neuro[ipc_peaks_vs_neuro >= 0.75]))
length(ipc_peaks)

neuro_peaks <- names(ipc_peaks_vs_neuro[ipc_peaks_vs_neuro <= -0.75])
length(neuro_peaks)

astro_peaks_vs_oligo <- get_genes_specific_to_mcs(a_legc, mc_pos = as.character(mcmd$metacell[which(mcmd$cell_type == 'Astrocytes')]), mc_neg = as.character(mcmd$metacell[which(mcmd$cell_type %in% c('OPCs'))]))
astro_peaks_vs_neuro <- get_genes_specific_to_mcs(a_legc, mc_pos = as.character(mcmd$metacell[which(mcmd$cell_type == 'Astrocytes')]), mc_neg = as.character(mcmd$metacell[which(mcmd$cell_type %in% c('CthPN', 'SCPN', 'CPN_L5_6', 'CPN_L2-3'))]))
astro_peaks <- union(names(nsc_peaks_vs_astro[nsc_peaks_vs_astro <= -0.75]), names(astro_peaks_vs_oligo[astro_peaks_vs_oligo >= 0.75]))
                       

nsc_peaks <- setdiff(nsc_peaks, multunion(ipc_peaks, astro_peaks, neuro_peaks))
ipc_peaks <- setdiff(ipc_peaks, multunion(nsc_peaks, astro_peaks, neuro_peaks))
neuro_peaks <- setdiff(neuro_peaks, multunion(nsc_peaks, astro_peaks, ipc_peaks))
astro_peaks <- setdiff(astro_peaks, multunion(ipc_peaks, neuro_peaks, nsc_peaks))

avg_meth_all <- gextract(meth_tracks, intervals = mca@peaks, iterator = mca@peaks)
cov_meth_all <- gextract(meth_cov_tracks, intervals = mca@peaks, iterator = 1)

cov_meth_all_peaks <- tgs_matrix_tapply(t(cov_meth_all[,grep('meth', colnames(cov_meth_all))]), cov_meth_all$intervalID, sum, na.rm = T)

avg_meth_all <- avg_meth_all[which(rowMins(cov_meth_all_peaks) >= 20),]

avg_meth_all$peak_name <- mca@peaks$peak_name[avg_meth_all$intervalID]
rownames(avg_meth_all) <- avg_meth_all$peak_name
# avg_meth_ct_peaks[,'ct'] <- unlist(purrr::map(stringr::str_split(ct_peaks_df$peak_name1, '--'), 1))
y <- colnames(avg_meth_all)[grep('E\\d\\d', colnames(avg_meth_all))]
colnames(avg_meth_all)[grep('E\\d\\d', colnames(avg_meth_all))] <- unlist(purrr::map(stringr::str_split(purrr::map(stringr::str_split(y, '\\.'), 2), '_'), 1))

save(avg_meth_all, file = './output/methylation/avg_meth_all.rda')

load(file = './output/methylation/avg_meth_all.rda')

load('./output/mcatac/var_peaks_after_enh_prom_separation.rda')

# load('./output//mcatac/mmcortex_feat_peak_variable_peak_clusters.rda')

tss <- gintervals.load('intervs.global.tss')

mc_rna <- scdb_mc('pl_cort')

tss2 <- tss[!duplicated(tss$geneSymbol),]

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

prom_peaks <- gintervals.neighbors(mcp, tss2, mindist = -1e+3, maxdist = 1e+3, maxneighbors = 1)

dist_peaks <- mcp[!(mcp$peak_name %in% prom_peaks$peak_name),]

prom_peaks <- prom_peaks$peak_name
dist_peaks <- dist_peaks$peak_name
yv <- intersect(intersect(dist_peaks, var_peaks), rownames(avg_meth_all))
yv2 <- intersect(setdiff(dist_peaks, var_peaks), rownames(avg_meth_all))

nei_peaks_tss <- tidyr::drop_na(gintervals.neighbors(mcp, tss2, 
                                                     mindist = -1e+6, maxdist = 1e+6, maxneighbors = 1))

mpra_lib <- readr::read_tsv('./cpg_methylation//st_and_temporal_and_e14_shadow_enh_9-1-22.tsv')
mpra_lib <- as.data.frame(dplyr::relocate(mpra_lib, rowname, .after = end))

ct_seq <- unlist(purrr::map(stringr::str_split(mpra_lib$rowname, '\\.'), 1))
ct_seq <- ifelse(ct_seq  %in%  c('IPC', 'NSC', 'Mature'), ct_seq, 'shadow')

dir_seq <- as.character(purrr::map(stringr::str_split(mpra_lib$rowname, '\\.'), 2))

mpra_lib$ct <- ct_seq
mpra_lib$dir <- dir_seq

mcp <- as.data.frame(dplyr::select(mca@peaks, chrom, start, end, peak_name))

nei_mpra_mcp <- tidyr::drop_na(gintervals.neighbors1(mpra_lib, mcp, maxdist = 0, mindist = 0, maxneighbors = 3))

meth_qs <- c(-1e-2,1e-2,5e-2,seq(0.1,1,0.1))

agg_id <- readr::read_csv('./scatac_data//aggregation_id.csv')
mca@cell_to_metacell$agg_id <- as.numeric(unlist(purrr::map(stringr::str_split(mca@cell_to_metacell$cell_id, '-'), 2)))
agg_id_day <- data.frame(cbind(1:14, rep(12:18, each = 2)))
colnames(agg_id_day) <- c('agg_id', 'day')
mca@cell_to_metacell$day <- agg_id_day$day[mca@cell_to_metacell$agg_id]

atac_mc_day <- table(mca@cell_to_metacell$metacell, mca@cell_to_metacell$day)
atac_mc_day_norm <- atac_mc_day/rowSums(atac_mc_day)
nsc_mcs <- which(mcmd$cell_type == 'NSC')
egc_by_day <- mca@egc[,nsc_mcs] %*% atac_mc_day_norm[nsc_mcs,]
egc_by_day_n <- t(t(egc_by_day)/colSums(egc_by_day))
colnames(egc_by_day_n) <- paste0('E', colnames(egc_by_day_n))
a_legc_by_day_n <- log2(1e-5 + egc_by_day_n)

nsc_inc_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] >= 1)], 
        rownames(avg_meth_all))
nsc_inc_atac <- a_legc_by_day_n[nsc_inc_peaks,]
nsc_dec_peaks <- intersect(rownames(a_legc_by_day_n)[which(a_legc_by_day_n[,'E17'] - a_legc_by_day_n[,'E13'] <= -1)], 
        rownames(avg_meth_all))
nsc_dec_atac <- a_legc_by_day_n[nsc_dec_peaks,]

nsc_asc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_asc')
nsc_asc_atac <- a_legc_by_day_n[nsc_asc_nei$peak_name,]

nsc_desc_nei <- dplyr::filter(nei_mpra_mcp, ct == 'NSC' & dir == 'seqs_desc')
nsc_desc_atac <- a_legc_by_day_n[nsc_desc_nei$peak_name,]

eb <- intersect(rownames(avg_meth_all), rownames(a_legc_by_day_n))
cn <- paste0('E', 13:17)
mx <- as.matrix(a_legc_by_day_n[eb,cn])
my <- as.matrix(avg_meth_all[eb,cn])

asc_in <- intersect(nsc_asc_nei$peak_name, rownames(mx))
desc_in <- intersect(nsc_desc_nei$peak_name, rownames(mx))

vasc <- colMeans(a_legc[nsc_asc_nei$peak_name,])
vdesc <- colMeans(a_legc[nsc_desc_nei$peak_name,])

legc <- log2(1e-5 + mc_rna@e_gc)

pltmt <- rbind(vasc,
               vdesc,
               legc[sort(grep('os|Sox10|Sox13|Sox7|Sox17|Sox15|Sox18|Sox30',
                         # c(
                             grep('Tet|Dnmt',rownames(mc_rna@e_gc), v=T) 
                           # sort(grep('Eomes|Pou3f2|Neurog2|Sox|Hox', rownames(mc_rna@e_gc), v=T)))
                        , inv=T, v=T)),mcmd$metacell]
               )

rownames(pltmt) <- c('inc_peaks', 'dec_peaks', tail(rownames(pltmt), -2))

nsc_accessible_peaks <- rownames(a_legc_by_day_n)[apply(a_legc_by_day_n, 1, max) >= -17.5]

dist_peaks_acc <- multintersect(dist_peaks, nsc_accessible_peaks, rownames(avg_meth_all))
prom_peaks_acc <- multintersect(prom_peaks, nsc_accessible_peaks, rownames(avg_meth_all))

save(mpra_lib,
        pltmt,
        prom_peaks,
        dist_peaks,
        a_legc_by_day_n,
        avg_meth_all,
        yv,
        yv2,
        dist_peaks_acc,
        prom_peaks_acc,
        nsc_peaks,
        ipc_peaks,
        astro_peaks,
        neuro_peaks,
        meth_qs,
        file = './output/methylation/fig4_meth_data.rda')