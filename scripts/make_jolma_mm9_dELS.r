library(tgutil)
library(misha)
library(glue)
library(dplyr)
library(gpwm)
doMC::registerDoMC(20)
wd = '/home/feshap/raid/proj/mmcortex'
gsetroot('/net/mraid14/export/tgdata/db/tgdb/mm9/trackdb/')
window_len = 300
jol_tracks = gtrack.ls('jolma_10bp')

clip_vals = function(df) {    
    for (chr in unique(ALLGENOME[[1]]$chrom)) { df[(df$chrom == chr) & (df$end > 
        ALLGENOME[[1]][ALLGENOME[[1]]$chrom == chr,'end']),'end'] = ALLGENOME[[1]][ALLGENOME[[1]]$chrom == chr,'end'] }
    for (chr in unique(ALLGENOME[[1]]$chrom)) { df[(df$chrom == chr) & (df$start < 0),'start'] = 0 }
    return(df)
}

mm9_dELS_intervals = gintervals.load('ENCODE_SCREEN_ccREs_dELS')
uniform_half_len = window_len/2
mm9_dELS = mm9_dELS_intervals %>% mutate(mid = round((end+start)/2)) %>% 
                        mutate(start_new = mid - uniform_half_len, end_new = mid + uniform_half_len) %>%
                        select(chrom, start_new, end_new) %>% 
                        rename(start = start_new, end = end_new)


mm9_dELS = clip_vals(mm9_dELS) %cache_df% file.path(wd, glue('mm9_dELS_intervals_{window_len}.csv'))



res = gpwm.extract(jol_tracks, mm9_dELS, gtrack.ls("jolma_10bp"), parallel = T) %cache_df% 
        file.path(wd, 'jolma_mm9_dELS.csv'))