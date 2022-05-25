library(misha)
gsetroot('/home/aviezerl/mm10')
mmc_path = '/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/tracks/mmcortex//tracks'
setwd(mmc_path)
options(gmultitasking=T, gmax.data.size=1e+9)
doMC::registerDoMC(cores = 32)

# day_dirs = dir(mmc_path)
# day_dirs

ex_tracks = lapply(2:7, function(n) gtrack.ls(glue::glue('(mmcortex).*(E1{n})')))

names(ex_tracks) = paste0('E', 12:17)

purrr::walk(ex_tracks, function(tr) {
    prefix = stringr::str_extract(tr[[1]], '(.*NSC)')
    contact_tracks = paste0(prefix, c('.rep1', '.rep2'))
    gtrack.create(paste0(prefix, '_sum_reps_1e4'), "sum of reps1+2 at 10kbp resolution", 
              iterator = rep(1e+4,2),expr = glue::glue('{contact_tracks[[1]]} + {contact_tracks[[2]]}'))
})