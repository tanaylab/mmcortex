library(gpwm)
library(misha)
library(dplyr)
gsetroot('/home/aviezerl/mm9')
options(gmax.data.size = 1e9)
options(dplyr.summarise.inform = FALSE)
wd = '/net//mraid14/export//tgdata//users//yonshap//proj//mmcortex/'

setwd(wd)

res = vroom::vroom(file.path(wd, 'km_clust_motifs_k_32_q_0.9.csv'))

tfoi = read.delim('./data/tfoi.txt', sep= '\n', header=F) %>% unlist %>% as.character %>% sort

all_tracks = gtrack.ls('motifs_10bp')
jol_tracks = gtrack.ls('jolma_10bp')

all_tracks_in = sort(all_tracks[lapply(toupper(tfoi), function(x) grep(x, toupper(all_tracks))) %>% unlist %>% unique])
jol_tracks_in = sort(jol_tracks[lapply(toupper(tfoi), function(x) grep(x, toupper(jol_tracks))) %>% unlist %>% unique])

tf_tracks = c(all_tracks_in, jol_tracks_in)
head(tf_tracks)
tail(tf_tracks)

mm9_dELS = vroom::vroom('~/raid/data_other//gastru_enhancer_yang_cell_res_2019/output//mm9_dELS_intervals_300.csv')

mmcortex_peaks = readRDS('./data/TF_motif.RDS')

tracks_all = tf_tracks

intervs = rbind(select(mutate(mm9_dELS, set = 'bg', cluster = -1), c('chrom', 'start', 'end', 'set', 'cluster')), 
                select(mutate(mmcortex_peaks, set = 'fg'), c('chrom', 'start', 'end', 'set', 'cluster'))) #%>% 
#                 rename(start = start_new, end = end_new))
# intervs$cluster = replace_na(intervs$cluster, -1)
prcs = c(0.95, 0.99, 0.995, 0.999)
# glob_val = gquantiles(paste0(track, '_max'), percentiles = prcs, intervals = intervs, iterator = intervs)

calc_dfc_by_cl = function(intervs, track, cl, glob_val) {
    intervs_cl = filter(intervs, (cluster == cl) | (cluster == -1))
    dfc_raw = intervs_cl %>% group_by(set) %>%
        summarise(n = n(), n_ok = length(which(val >= glob_val)))
#     dfc_raw[is.na(dfc_raw)] = 0
    dfc = tidyr::pivot_longer(dfc_raw, c('n', 'n_ok')) %>% 
            mutate(var_name = paste(set, name, sep='_')) %>% 
            select(var_name, value) %>% 
            tibble::column_to_rownames('var_name') %>% t
    dfc = data.frame(dfc)
    rownames(dfc) = paste(track, 'cl', cl, sep = '_')
    dfc$cluster = cl
    return(dfc)
}

count_motif_thresh = function(track, intervs, prc) {
    print(track)
    max_track = gvtrack.create(paste0(track, '_max'), track, 'max')
    interv_max = gextract(paste0(track, '_max'), intervals = intervs, iterator = intervs) %>% rename(val = 4)
    intervs = left_join(intervs, interv_max, by=c('chrom', 'start', 'end'))
    glob_val = gquantiles(paste0(track, '_max'), percentiles = prc, intervals = intervs, iterator = intervs)
    clusts = unique(intervs$cluster)
    clusts = clusts[clusts > 0]
    cl_enrich = purrr::map_dfr(clusts, function(x) calc_dfc_by_cl(intervs, track, x, glob_val))
    return(cl_enrich)
}

motif_c = purrr::map_dfr(tf_tracks, function(x) count_motif_thresh(x, intervs, prcs[[2]]))

motif_c$pv = phyper(motif_c[,'fg_n_ok'], 
       motif_c[,'fg_n_ok'] + motif_c[,'bg_n_ok'], 
       motif_c[,'fg_n'] + motif_c[,'bg_n'] - motif_c[,'fg_n_ok'] - motif_c[,'bg_n_ok'], 
       motif_c[,'fg_n'], 
       lower.tail = FALSE)
motif_c$qval = p.adjust(motif_c$pv)

motif_c$rel_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (motif_c$bg_n_ok / motif_c$bg_n)
motif_c$abs_enrich = (motif_c$fg_n_ok / motif_c$fg_n) / (1 - prcs[[2]])

motif_c = tibble::rownames_to_column(motif_c)

# motif_c$motif = stringr::str_split(motif_c$rowname, '\\.|_') %>% lapply(function(x) paste0(x[3:(length(x)-2)], collapse='_'))
# motif_c



readr::write_tsv(motif_c, file.path(wd, 'data','motif_enrich_k=32_q_0.9.tsv'))