library(misha)
library(dplyr)
gsetroot('/home/aviezerl/mm9')
options(gmax.data.size = 1e9)

wd = '/net/mraid14//export/tgdata/users/yonshap/proj/mmcortex/'
output_dir = file.path(wd, 'data')
if (!dir.exists(output_dir)) {dir.create(output_dir)}
tracks_dir = file.path(output_dir, 'track_data')
if (!dir.exists(tracks_dir)) {dir.create(tracks_dir)}
SEED = 1337
set.seed(SEED)
setwd(wd)

prcs = c(0.95, 0.99, 0.995, 0.999)

motif_tracks = gtrack.ls('motifs_10bp')

jolma_tracks = gtrack.ls('jolma_10bp\\.')

tracks_all = c(motif_tracks, jolma_tracks)

tfoi = read.delim('./data/tfoi.txt') %>% unlist %>% sort %>% as.character
tfoi_match = lapply(tfoi, function(tf) grep(tf, tracks_all,v=T, ignore.case = T))
names(tfoi_match) = tfoi
# tfoi_match

# peaks_all = read.delim('./data/peak_locs.txt', header = F)
# peaks_all = as.data.frame(do.call(rbind, sapply(apply(peaks_all, 1, stringr::str_split, '-'), identity)))
# colnames(peaks_all) = c('chrom', 'start', 'end')
# peaks_all[,c('start', 'end')] = apply(peaks_all[,c('start', 'end')], 2, as.numeric)

# chain_path = '/home/aviezerl/proj/ebdnmt/rawdata/import/Weber_Nature_Communication_2020/mm10ToMm9.over.chain.fixed1'

# peaks_all_mm9 = gintervals.liftover(peaks_all, chain = chain_path)

# tracks_dir

tracks_folders = list.files(tracks_dir)
head(tracks_folders)

# list.files(file.path(tracks_dir, tracks_folders[[1]]))

# file_fold = file.path(tracks_dir, tracks_folders[[1]])
# interv_max = readr::read_tsv(file.path(file_fold, list.files(file_fold)[[1]]), show_col_types = FALSE)
# interv_max = unique(interv_max[,1:4])

top_peaks_per_mc = readr::read_tsv('./data/top_peak_per_mc.tsv', show_col_types = F)
# head(top_peaks_per_mc)

top_peak_dfs_list = lapply(1:ncol(top_peaks_per_mc), function(i) {
    interv_mat = data.frame(do.call('rbind', stringr::str_split(unlist(top_peaks_per_mc[,i]), '-')))
    colnames(interv_mat) = c('chrom', 'start', 'end')
    interv_mat[,c('start', 'end')] = apply(interv_mat[,c('start', 'end')], 2, as.numeric)
    interv_mat = interv_mat[with(interv_mat, order(chrom, start, end)),]
})

purrr::walk(unlist(tfoi_match), function(f) gvtrack.create(vtrack = paste0(f, '_max'), src = f, func = 'max'))

count = 1
mc_peak_tf_energy = lapply(top_peak_dfs_list, function(x, count) {
    print(count)
    sapply(unlist(tfoi_match), function(trk) {
        return(gextract(paste0(trk, '_max'), intervals = x))
    })
    count = count + 1;
    }, count
)

saveRDS(mc_peak_tf_energy, './data/mc_peak_tf_energy.rds')

# quantiles_all = sapply(tracks_folders, function(f) {
#     file_fold = file.path(tracks_dir, f)
#     interv_max = readr::read_tsv(file.path(file_fold, list.files(file_fold)[[1]]), show_col_types = FALSE)
#     qs = quantile(interv_max$val, c(seq(0,0.95,0.05), 0.99,0.995,0.999,1))
#     return(qs)
# })

# head(quantiles_all)

# as.numeric(gsub('%', '', rownames(quantiles_all)))

# qa_km=  tglkmeans::TGL_kmeans(t(quantiles_all[(n-6):n,]), k = 4)

# col_key = data.frame(cbind(sort(unique(qa_km$cluster)), sample(grep('white|gray|grey', colors(), v=T, inv=T), max(unique(qa_km$cluster)))))
# colnames(col_key) = c('cl', 'color')
# col_key

# qa_col = setNames(col_key$color[match(qa_km$cluster, col_key$cl)], qa_km$cluster)
# qa_col

# plot(0, xlim = c(0, 100), ylim = c(-50,0))
# n = nrow(quantiles_all)
# purrr::walk(1:ncol(quantiles_all), function(i) {
#     x = as.numeric(gsub('%', '', rownames(quantiles_all)));
# #     x = tail(x, 6);
#     y = quantiles_all[,i];
# #     y = tail(y, 6)
# #     y = (y - min(y))/(max(y) - min(y));
#     lines(x,y, cex = 1.5, lwd = 1,
#                                 col = qa_col[[i]])
# #     points(x,y, cex = 1.5, col = 'black', pch =21)

# }
#        )

# motif_enrich = readr::read_csv('./data/mmcortex_motif_enrich.csv', show_col_types = F)

# tf_q = unique(motif_enrich$gene_egc[motif_enrich$qval <= 1])
# tf_q = sort(tf_q[!is.na(tf_q)])
# # tf_q = tf_q[tf_q %in% tfoi]

# mc_motif = tidyr::pivot_wider(select(filter(motif_enrich, gene_egc %in% tf_q), gene_egc, mc, rel_enrich), names_from = mc, values_from = rel_enrich
# #                               names_from = mc, values_from = rel_enrich
#                              )

# motif_max = tibble::column_to_rownames(
#                     tidyr::drop_na(
#                         as.data.frame(t(apply(mc_motif, 1, function(x) sapply(x, function(y) max(y)))))
#                                     ),
#                                               'gene_egc')

# motif_max = motif_max[,order(as.numeric(colnames(motif_max)))]
# motif_max = motif_max[order(rownames(motif_max)),]
# rn = rownames(motif_max)
# motif_max = apply(motif_max, 2, as.numeric)
# rownames(motif_max) = rn

# head(motif_max)

# dim(quantiles_all)

# head(interv_max)

# quantile(interv_max$val, seq(0,1,0.05))

# dim(unique(interv_max[,1:4]))

# head(interv_max)


