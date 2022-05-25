library(misha)
library(dplyr)
gsetroot('/home/aviezerl/mm10')
options(gmax.data.size = 1e+9)
#gdb.reload()

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

jaspar_tracks = gtrack.ls('jaspar_10bp')

cisbp_tracks = gtrack.ls('cis_bp_10bp')

tracks_all = c(motif_tracks, jolma_tracks, jaspar_tracks, cisbp_tracks)

tracks_all_tfs = stringr::str_split(tracks_all, '\\.') %>% purrr::map(2) %>% stringr::str_split('_') %>% purrr::map(1) %>% unlist
#tracks_all = cisbp_tracks

tfoi = sort(c('Bcl11b', read.delim('./data/tfoi.txt') %>% unlist %>% sort %>% as.character))
#tfoi_match = lapply(tfoi, function(tf) grep(tf, tracks_all,v=T, ignore.case = T))
tfoi_match = lapply(tfoi, function(tf) tracks_all[grep(tf, tracks_all_tfs, ignore.case = T)])
names(tfoi_match) = tfoi
# tfoi_match

peaks_all = read.delim('./data/peak_locs.txt', header = F)
peaks_all = as.data.frame(do.call(rbind, sapply(apply(peaks_all, 1, stringr::str_split, '-'), identity)))
colnames(peaks_all) = c('chrom', 'start', 'end')
peaks_all[,c('start', 'end')] = apply(peaks_all[,c('start', 'end')], 2, as.numeric)

# for (i in c(100,500,2500)) {
for (i in c(10000)) {
    n_peaks = i
    print(as.character(c('Running analysis for ', i, ' peaks')))
    top_peak_per_mc = readr::read_tsv(glue::glue('./data/top_{n_peaks}_peak_per_mc.tsv'))

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
    
    tracks_to_test = c(
				unlist(lapply(tfoi_match, function(x) {if (length(x) > 2) {sample(x,2)} else {x}})),
				 sample(tracks_all[!(tracks_all %in% unlist(tfoi_match))], length(tfoi_match))
			)
#unlist(tfoi_match)
    print(tracks_to_test)
#sample(
                        #    10)

    # mm9_dELS = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/data_other/gastru_enhancer_yang_cell_res_2019/output/mm9_dELS_intervals_300.csv')
    mm9_dELS = peaks_all_mm9

    my_peaks_file = ints_per_mc_mm9

    intervs = rbind(select(mutate(mm9_dELS, set = 'bg', mc = -1), c('chrom', 'start', 'end', 'mc', 'set')), 
                    select(mutate(my_peaks_file, set = 'fg'), c('chrom', 'start', 'end', 'mc', 'set'))) %>% 
                    arrange(chrom, start, end)

    intervs_dedup = intervs[!duplicated(intervs[,1:3]),]
    print('Creating max-tracks')
    sapply(tracks_to_test, function(trk) gvtrack.create(paste0(trk, '_max'), trk, 'max'))
    num_groups = ceiling(length(tracks_to_test)/50)
    trk_groups = sample.int(4, length(tracks_to_test), replace = T)
    print('Calculating motif energies per intervals')
    intervs_max = do.call('cbind', lapply(1:num_groups, function(grp) {
			print(c('group ', grp)); 
			trk_all_raw = gextract(tracks_to_test[trk_groups == grp], intervals = intervs_dedup, iterator = 10)
            trk_all_raw = trk_all_raw[with(trk_all_raw, order(chrom, start, end)),]
            trk_all_avg = tgs_matrix_tapply(t(subset(trk_all_raw, select = -c(chrom, start, end, intervalID))), trk_all_raw$intervalID, function(x) log(sum(exp(x))))
            return(trk_all_avg)
            # gextract(paste0(tracks_to_test[trk_groups == grp], '_max'), intervals = intervs_dedup, iterator = intervs_dedup)
									}
					)
			)
    intervs_max = cbind(intervs_dedup, intervs_max)
    colnames(intervs_max) = c('chrom', 'start', 'end', tracks_to_test[order(trk_groups)])
#     intervs_max = intervs_max[,c(1:3,which(!colnames(intervs_max) %in% c('chrom', 'start', 'end', 'intervalID')))]
#     intervs_max = intervs_max[,c(1:3, 3+order(sapply(colnames(intervs_max[,-c(1:3)]), function(cn) grep(cn, paste0(tracks_to_test, '_max')))))]
#     print('Getting global values/quantiles')
#     prc = prcs[[2]]
#     glob_vals = apply(intervs_max[,-c(1:3)], 2, function(x) quantile(x, prc))

#     intervs = left_join(intervs, intervs_max, by=c('chrom', 'start', 'end'))
#     colnames(intervs)[-c(1:5)] = tracks_to_test
#     # print(head(intervs))

#     calc_dfc_by_cl = function(intervs, track, cl, glob_val) {
#     #     intervs_cl = filter(intervs, (cluster == cl) | (cluster == -1))
#     # print(cl)
 
#         intervs_cl = filter(intervs, (mc == cl) | (mc == -1))
#         dfc_raw = intervs_cl %>% select(track, set) %>% group_by(set) %>% rename(val = track) %>% summarise(n = n(), n_ok = sum(val >= glob_val))
#         dfc = tidyr::pivot_longer(dfc_raw, c('n', 'n_ok')) %>%
#             mutate(var_name = paste(set, name, sep='_')) %>%
#             select(var_name, value) %>%
#             tibble::column_to_rownames('var_name') %>% t
#         dfc = data.frame(dfc)
#         rownames(dfc) = paste(track, 'mc', cl, sep = '_')
#         dfc$mc = cl
#         return(dfc)
#     }        
    
#     clusts = unique(intervs$mc)
#     clusts = sort(clusts[clusts > 0])
#     print('Calculating metacell enrichment')
#     n_ok_mat = sapply(seq_along(intervs[,-c(1:5)]), function(x, gv, i) {
#     #     print(length(x[,i]))
#     #     print(length(x$mc))
#         tapply(x[,i], intervs$mc, function(y) length(which(y > gv[[i]])))
#                }, x = intervs[,-c(1:5)], gv = glob_vals)

#     # head(n_ok_mat)

#     colnames(n_ok_mat) = colnames(intervs[-c(1:5)])

#     # nrow(peaks_all_mm9)

#     bg_rat = n_ok_mat[1,]/nrow(peaks_all_mm9)

#     ints_per_mc_tbl = as.numeric(table(ints_per_mc_mm9$mc))

#     fg_rat = n_ok_mat[-c(1),]/ints_per_mc_tbl

#     p_tbl = phyper(n_ok_mat[-c(1),], 
#            n_ok_mat[-c(1),] + n_ok_mat[1,], 
#            ints_per_mc_tbl + nrow(peaks_all_mm9) - n_ok_mat[-c(1),] - n_ok_mat[1,],
#            ints_per_mc_tbl
#           )

#     q_tbl = matrix(p.adjust(p_tbl), nrow(p_tbl), ncol(p_tbl))

#     # grep('eomes', 
#          colnames(fg_rat)[apply(q_tbl, 2, function(x) any(x< 0.2))]
#     #                             , ign =T)

#     re_tbl = t(apply(fg_rat, 1, function(x) x/bg_rat))
#     abs_e_tbl = fg_rat/(1-prc)

#     res = list('intervs' = intervs, 'glob_vals' = glob_vals, 'bg_n_ok' = bg_rat, 'fg_n_ok' = fg_rat, 're' = re_tbl, 'p_tbl' = p_tbl, 'q_tbl' = q_tbl)
#     saveRDS(res, './data/motif_enrich_10000.rds')
#     #motif_c = purrr::map_dfr(clusts, function(mc) purrr::map2_dfr(tracks_to_test, glob_vals, function(.x, .y) calc_dfc_by_cl(intervs, .x, mc, .y)))
    
#     #motif_c = purrr::map_dfr(tracks_to_test, function(x) count_motif_thresh(x, intervs, prcs[[2]]))
#     #readr::write_excel_csv(motif_c, file.path(output_dir, glue::glue('motif_vals_{n_peaks}_test_new.csv')))

#     #motif_c$pv = phyper(motif_c[,'fg_n_ok'], 
#     #    motif_c[,'fg_n_ok'] + motif_c[,'bg_n_ok'], 
#     #    motif_c[,'fg_n'] + motif_c[,'bg_n'] - motif_c[,'fg_n_ok'] - motif_c[,'bg_n_ok'], 
#     #    motif_c[,'fg_n'], 
#     #    lower.tail = FALSE)
#     #motif_c$qval = p.adjust(motif_c$pv)

#     #motif_c$rel_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (motif_c$bg_n_ok / motif_c$bg_n)
#     #motif_c$abs_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (1 - prcs[[2]])

#     #motif_c = tibble::rownames_to_column(motif_c)

#     #readr::write_excel_csv(motif_c, file.path(output_dir, glue::glue('motif_enrich_all_{n_peaks}_test_new.csv')))
# # readr::write_excel_csv(motif_c, file.path(output_dir, 'motif_enrich_test.csv'))
# }

