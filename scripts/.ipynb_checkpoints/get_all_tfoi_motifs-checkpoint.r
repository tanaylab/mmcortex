library(misha)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')

EXT = 133

tfoi = sort(c('Bcl11b', unlist(read.delim('./data/tfoi.txt', header = F))))
#tfoi
tfoi_trk = tfoi

dels = gintervals.load('ENCODE_SCREEN_ccREs_dELS')

diffs = dels$end - dels$start
mids = round(0.5*(dels$start + dels$end))
dels$start = mids - EXT
dels$end = mids + EXT
dels = dels[with(dels, order(chrom, start, end)),]

for (chri in sort(unique(ALLGENOME[[1]][,'chrom']))) {
	chri_inds = which(dels$chr == chri)
	inds_clip = which(dels$end[chri_inds] > ALLGENOME[[1]][ALLGENOME[[1]]$chrom == chri,'end'])
	dels$end[chri_inds[inds_clip]] = ALLGENOME[[1]][ALLGENOME[[1]]$chrom == chri,'end']
}


trks_motifs = unlist(sapply(c('motifs_10bp', 'jaspar_10bp', 'jolma_10bp', 'cis_bp_10bp'), gtrack.ls))
tf_tracks = unlist(sapply(tfoi_trk, function(x) grep('oct4', grep(x, trks_motifs, v=T, ign = T), ign=T, v=T, inv=T)))
tf_tracks

trk_all_raw = gextract(tf_tracks, intervals = dels, iterator = 10)
trk_all_raw = trk_all_raw[with(trk_all_raw, order(chrom, start, end)),]

trk_all_avg = tgs_matrix_tapply(t(subset(trk_all_raw, select = -c(chrom, start, end, intervalID))), trk_all_raw$intervalID, function(x) log(sum(exp(x))))

# purrr::walk(tf_tracks, function(trk) gvtrack.create(paste0(trk, '_max'), trk, 'max'))

# trk_maxs = gextract(paste0(tf_tracks, '_max'), intervals = dels, iterator = dels, file = './data/motif_max_ENCODE_SCREEN_dELS_266bp.tsv')
trk_log_sum_exps = cbind(dels, trk_all_avg)

# saveRDS(trk_maxs, './data/motif_max_ENCODE_SCREEN_dELS_266bp.rds')

saveRDS(trk_log_sum_exps, './data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds')
