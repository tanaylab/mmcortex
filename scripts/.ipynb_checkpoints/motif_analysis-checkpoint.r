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

peaks_all = read.delim('./data/peak_locs.txt', header = F)
peaks_all = as.data.frame(do.call(rbind, sapply(apply(peaks_all, 1, stringr::str_split, '-'), identity)))
colnames(peaks_all) = c('chrom', 'start', 'end')
peaks_all[,c('start', 'end')] = apply(peaks_all[,c('start', 'end')], 2, as.numeric)

top_peak_per_mc = readr::read_tsv('./data/top_peak_per_mc.tsv')

top_peak_long = tidyr::pivot_longer(top_peak_per_mc, cols = everything(), names_to = 'mc') %>% arrange(mc)
top_peak_long[,c('chrom', 'start', 'end')] = do.call(rbind, 
               stringr::str_split(t(top_peak_long[,'value']), '-'))
top_peak_long = top_peak_long[,c('chrom', 'start', 'end', 'mc')]
top_peak_long[,c('start', 'end', 'mc')] = apply(top_peak_long[,c('start', 'end', 'mc')], 2, as.numeric)

chain_path = '/home/aviezerl/proj/ebdnmt/rawdata/import/Weber_Nature_Communication_2020/mm10ToMm9.over.chain.fixed1'

peaks_all_mm9 = gintervals.liftover(peaks_all, chain = chain_path)

ints_per_mc_mm9 = # data.frame(cbind(
    gintervals.liftover(intervals = top_peak_long, chain_path)
# , top_peak_long$mc))

ints_per_mc_mm9$mc = top_peak_long$mc[ints_per_mc_mm9$intervalID]

tracks_to_test = #sample(
                        c(unlist(tfoi_match), sample(tracks_all[!(tracks_all %in% unlist(tfoi_match))], length(tfoi_match)))
                    #    10)

# mm9_dELS = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/data_other/gastru_enhancer_yang_cell_res_2019/output/mm9_dELS_intervals_300.csv')
mm9_dELS = peaks_all_mm9

my_peaks_file = ints_per_mc_mm9

intervs = rbind(select(mutate(mm9_dELS, set = 'bg', mc = -1), c('chrom', 'start', 'end', 'mc', 'set')), 
                select(mutate(my_peaks_file, set = 'fg'), c('chrom', 'start', 'end', 'mc', 'set'))) %>% 
                arrange(chrom, start, end)

# print(head(intervs))

# glob_val = gquantiles(paste0(track, '_max'), percentiles = prcs, intervals = intervs, iterator = intervs)

calc_dfc_by_cl = function(intervs, track, cl, glob_val) {
#     intervs_cl = filter(intervs, (cluster == cl) | (cluster == -1))
    # print(cl)
    intervs_cl = filter(intervs, (mc == cl) | (mc == -1))
    dfc_raw = intervs_cl %>% group_by(set) %>%
        summarise(n = n(), n_ok = sum(val >= glob_val), `.groups` = )
    dfc = tidyr::pivot_longer(dfc_raw, c('n', 'n_ok')) %>% 
            mutate(var_name = paste(set, name, sep='_')) %>% 
            select(var_name, value) %>% 
            tibble::column_to_rownames('var_name') %>% t
    dfc = data.frame(dfc)
    rownames(dfc) = paste(track, 'mc', cl, sep = '_')
    dfc$mc = cl
    return(dfc)
}

count_motif_thresh = function(track, intervs, prc) {
    print(track)
    max_track = gvtrack.create(paste0(track, '_max'), track, 'max')
    interv_max = gextract(paste0(track, '_max'), intervals = intervs, iterator = intervs) %>% 
        rename(val = 4) %>%
        arrange(chrom, start, end)
    intervs$val = interv_max$val
    glob_val = gquantiles(paste0(track, '_max'), percentiles = prc, intervals = intervs, iterator = intervs)
    track_dir = file.path(tracks_dir,track)
    if (!dir.exists(track_dir)) {dir.create(track_dir)}
    readr::write_tsv(interv_max, file.path(track_dir, 'interv_max.tsv'))
    clusts = unique(intervs$mc)
    clusts = sort(clusts[clusts > 0])
    cl_enrich = purrr::map_dfr(clusts, function(x) calc_dfc_by_cl(intervs, track, x, glob_val))
    # print(head(cl_enrich))
    return(cl_enrich)
}

# count_motif_thresh = function(track, intervs, prc) {
#     print(track)
#     max_track = gvtrack.create(paste0(track, '_max'), track, 'max')
#     interv_max = gextract(paste0(track, '_max'), intervals = intervs, iterator = intervs) %>% rename(val = 4)
#     intervs = left_join(intervs, interv_max, by=c('chrom', 'start', 'end'))
#     track_dir = file.path(tracks_dir,track)
#     if (!dir.exists(track_dir)) {dir.create(track_dir)}
#     readr::write_tsv(interv_max, file.path(track_dir, 'interv_max.tsv'))
#     glob_val = gquantiles(paste0(track, '_max'), percentiles = prc, intervals = intervs, iterator = intervs)
#     dfc_raw = intervs %>% group_by(set) %>%
#         summarise(n = n(), n_ok = sum(val >= glob_val))
#     dfc = tidyr::pivot_longer(dfc_raw, c('n', 'n_ok')) %>% 
#             mutate(var_name = paste(set, name, sep='_')) %>% 
#             select(var_name, value) %>% 
#             tibble::column_to_rownames('var_name') %>% t
#     rownames(dfc) = track
#     return(data.frame(dfc))
# }


motif_c = purrr::map_dfr(tracks_to_test, function(x) count_motif_thresh(x, intervs, prcs[[2]]))
readr::write_excel_csv(motif_c, file.path(output_dir, 'motif_vals.csv'))

motif_c$pv = phyper(motif_c[,'fg_n_ok'], 
       motif_c[,'fg_n_ok'] + motif_c[,'bg_n_ok'], 
       motif_c[,'fg_n'] + motif_c[,'bg_n'] - motif_c[,'fg_n_ok'] - motif_c[,'bg_n_ok'], 
       motif_c[,'fg_n'], 
       lower.tail = FALSE)
motif_c$qval = p.adjust(motif_c$pv)

motif_c$rel_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (motif_c$bg_n_ok / motif_c$bg_n)
motif_c$abs_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (1 - prcs[[2]])

motif_c = tibble::rownames_to_column(motif_c)

readr::write_excel_csv(motif_c, file.path(output_dir, 'motif_enrich_all.csv'))
# readr::write_excel_csv(motif_c, file.path(output_dir, 'motif_enrich_test.csv'))


