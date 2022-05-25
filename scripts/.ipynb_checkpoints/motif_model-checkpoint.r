
prc = 0.999
EXT = 250
n_peaks = 500
tfoi = sort(c('Bcl11b', unlist(read.table('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex/data/tfoi.txt'))))
tfs_extract = c('Eomes', 'Pou3f2')
tfoi


library(dplyr)
library(misha)
gsetroot('/home/aviezerl/mm9')
library(metacell)
library(glmnet)

scdb_init('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex/scdb')
chain_path = '/home/aviezerl/proj/ebdnmt/rawdata/import/Weber_Nature_Communication_2020/mm10ToMm9.over.chain.fixed1'
mc = scdb_mc('pl_cort')
mcmd = readr::read_tsv('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex//BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN/CfuPN',
                'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))


cat('Importing top peaks per mc...\n')
top_peak_per_mc = readr::read_tsv(glue::glue('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex//data//top_{n_peaks}_peak_per_mc.tsv'))

top_peak_long = tidyr::pivot_longer(top_peak_per_mc, cols = everything(), names_to = 'mc') %>% arrange(mc)
top_peak_long[,c('chrom', 'start', 'end')] = do.call(rbind, 
               stringr::str_split(t(top_peak_long[,'value']), '-'))
top_peak_long = top_peak_long[,c('chrom', 'start', 'end', 'mc')]
top_peak_long[,c('start', 'end', 'mc')] = apply(top_peak_long[,c('start', 'end', 'mc')], 2, as.numeric)


ints_per_mc_mm9 = gintervals.liftover(intervals = top_peak_long, chain_path)

ints_per_mc_mm9$mc = top_peak_long$mc[ints_per_mc_mm9$intervalID]

tfoi_trk = tfoi
trks_motifs = unlist(sapply(c('motifs_10bp', 'jaspar', 'jolma', 'cis_bp_10bp'), gtrack.ls))
tf_tracks = unlist(purrr::map(sapply(tfoi_trk, function(x) grep('oct4', grep(x, trks_motifs, v=T, ign = T), ign=T, v=T, inv=T)), 1))

cat('Creating max tracks of TF binding energy...\n')
purrr::walk(tf_tracks, function(trk) gvtrack.create(paste0(trk, '_max'), src = trk, func = 'max'))

intervs = ints_per_mc_mm9

intervs[,'start'] = intervs[,'start'] - EXT
intervs[,'end'] = intervs[,'end'] + EXT

get_top_ints = function(track) {
    intervs_dedup = intervs[!duplicated(intervs[,1:3]),]
    max_track = gvtrack.create(paste0(track, '_max'), track, 'max')
    interv_max = gextract(paste0(track, '_max'), intervals = intervs_dedup, 
                          iterator = intervs_dedup) %>% rename(val = 4) %>% 
                arrange(chrom, start, end)
    intervs = left_join(intervs, interv_max, by=c('chrom', 'start', 'end'))

    glob_val = gquantiles(paste0(track, '_max'), percentiles = prc, intervals = intervs_dedup, iterator = intervs_dedup)

    intervs_f = intervs[intervs$val >= glob_val,]
    intervs_f = intervs_f[!duplicated(intervs_f[,1:3]),]
#     intervs_f$tf = purrr::map(stringr::str_split(track, '\\.'), function(x) unlist(x[length(x)]))
    return(intervs_f)
}

cat('Importing intergenic peak accessibility matrix...\n')
pmc_mm9 = readRDS('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex//data/pl_cort_peak_mc_smoothed_mg_mm9.rds')

pmc_mm9$start = pmc_mm9$start - EXT
pmc_mm9$end = pmc_mm9$end + EXT

tracks_to_extract = tf_tracks[unlist(sapply(tfs_extract, function(tf) grep(tf, tf_tracks, ign=T)))]
tf_track_ints = lapply(tracks_to_extract, get_top_ints)

all_ints = do.call('rbind', tf_track_ints)
all_ints = all_ints[!duplicated(all_ints[,1:3]),]


cat('Extracting TF binding energy for intergenic peaks...\n')
all_hits_all_tfs = gextract(paste0('exp(', tf_tracks, '_max)'), intervals = all_ints[,1:3], iterator = all_ints[,1:3])

tfoi_trk_in = tfoi_trk[sapply(tfoi_trk, function(x) length(grep(x, colnames(all_hits_all_tfs), ign=T)) > 0)]

cat('Creating matrix of RNA-TF energy geometric means per peak per metacell...\n')
rna_xx = apply(t(mc@e_gc[tfoi_trk_in,]), 1, function(x) apply(all_hits_all_tfs[,4:(ncol(all_hits_all_tfs)-1)], 1, 
                                                                                    function(y) apply(rbind(y, x), 2, function(yx) exp(mean(log(yx))))))

rna_xxx_r = matrix(rna_xx, ncol = length(tfoi_trk_in), nrow = ncol(mc@e_gc)*nrow(all_hits_all_tfs), byrow = T)


cat('Organizing train/test data...\n')
jn = dplyr::right_join(pmc_mm9, all_hits_all_tfs[,1:3], by=c('chrom', 'start', 'end'))

y = as.numeric(unlist(jn[,5:ncol(jn)]))

sample_inds = sort(sample(1:nrow(rna_xxx_r), size = 0.05*nrow(rna_xxx_r)))

x = as.matrix(rna_xxx_r[!(1:nrow(rna_xxx_r) %in% sample_inds),])
y_train = y[!(1:nrow(rna_xxx_r) %in% sample_inds)]

x_test = as.matrix(rna_xxx_r[sample_inds,])
y_test = y[sample_inds]

cat('Running model...\n')
lambdas <- 10**seq(6, -6, by = -.1)

cv_ridge <- cv.glmnet(x, y_train, alpha = 1, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}

ridge_reg_opt = glmnet(x, y_train, lambda = optimal_lambda, alpha = 1, family = 'gaussian')

predictions_train <- predict(ridge_reg_opt, s = optimal_lambda, newx = x)
predictions_test <- predict(ridge_reg_opt, s = optimal_lambda, newx = x_test)

print('Results (on train set):')
print(eval_results(y_train, predictions_train, x))
print('Results (on test set):')
print(eval_results(y_test, predictions_test, x_test))

# smoothScatter(log2(y_train + min(y_train[y_train > 0])), log2(predictions_train + min(predictions_train[predictions_train > 0])))