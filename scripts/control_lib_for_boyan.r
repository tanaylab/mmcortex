ls()

wd = '/net/mraid14//export/tgdata/users//yonshap/proj//mmcortex'

setwd(wd)

library(misha)
gsetroot('/home/aviezerl/mm10')

library(dplyr)

library("Biostrings")

s = readDNAStringSet("./data/scrContr_preFilt.fasta")

homer_key = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/homer.key',col_names = c('key', 'track', 'dummy'))
homer_data = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/homer.data', col_names = c('key', 'pos', 'A', 'C', 'G', 'T'))
cis_bp_key = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.key',col_names = c('key', 'track', 'dummy'))
cis_bp_data = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.data', col_names = c('key', 'pos', 'A', 'C', 'G', 'T'))
jolma_key = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/jolma.key',col_names = c('key', 'track', 'dummy'))
jolma_data = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/jolma.data', col_names = c('key', 'pos', 'A', 'C', 'G', 'T'))
cis_key = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.key',col_names = c('key', 'track', 'dummy'))
cis_data = readr::read_tsv('/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.data', col_names = c('key', 'pos', 'A', 'C', 'G', 'T'))
jaspar_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.key', header = F, col.names = c('key', 'track', 'dummy'))
jaspar_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))

tracks = as.character(unlist(read.delim('./data/motif_tracks_selected_k=64_n=64.txt')))
tracks

get_pwm = function(trk) {
    db = stringr::str_split(trk, '\\.') %>% purrr::map(1) %>% unlist
    db = gsub('_10bp' ,'', stringr::str_split(trk, '\\.') %>% purrr::map(1))
    tf = stringr::str_split(trk, '\\.') %>% purrr::map(2) %>% unlist
    print(c(trk, db, tf))
    if (db == 'jaspar') {
        pwm = t(jaspar_data[jaspar_data$key[jaspar_key$track == tf],3:6])}
    else if (db == 'motifs') {
        pwm = t(homer_data[homer_data$key[homer_key$track == tf],3:6])}
    else if (db == 'jolma') {
        pwm = t(jolma_data[jolma_data$key[jolma_key$track == tf],3:6])}
    else {
        pwm = t(cis_data[cis_data$key[cis_key$track == tf],3:6])}
    return(t(pwm))
}

motif_mats = lapply(tracks, function(trk) get_pwm(trk))

names(motif_mats) = tracks

motif_mats_log = lapply(motif_mats, function(mat) {mat_n = apply(mat + 0.01, 1, function(x) x/sum(x)); return(apply(mat_n + 0.01,1, log))})

log_sum_log = function(l1,l2) {
	if(l1 > l2) {
		if(!is.infinite(l2)) {
			l1 = l1 + log(1+ exp(l2-l1))
		}
	} else {
		if(is.infinite(l1)) {
			l1 = l2
		} else {
			l1 = l2 + log(1 + exp(l1-l2));
		}
	}
}

run_motif_on_window_rc = function(mm, seq) {
                                        seq_rc = as.character(reverseComplement(DNAString(seq)));
                                        seq_s = unlist(strsplit(seq, ''));
                                         seq_rc_s = unlist(strsplit(as.character(seq_rc), ''));
                                        vals_to_sum = sapply(seq_along(seq_s), function(i) mm[i,seq_s[[i]]])
                                        sum_sense = sum(vals_to_sum)
                                        sum_rc = sum(sapply(seq_along(seq_rc_s), function(i) mm[i,seq_rc_s[[i]]]))
                                        return(log_sum_log(sum_sense, sum_rc))
                                                           
                                                           }

sc = lapply(s, as.character)

names(sc) = 1:length(sc)

calc_seq_track_val = function(trk, seq) {
    motif_vals_rc = sapply(1:(nchar(seq) - nrow(motif_mats_log[[trk]])), 
       function(i) run_motif_on_window_rc(motif_mats_log[[trk]], 
                                        substr(seq, i, i+nrow(t(motif_mats_log[[trk]])-1))
                                      )
                    )
    l = length(motif_vals_rc)
    grp_vec = c(sort(rep(1:floor(l/10), 10)), rep(ceiling(l/10), l %% 10))
    return(log(sum(exp(tapply(motif_vals_rc, grp_vec, mean)))))
}

vals_all = mclapply(X = 1:length(sc), FUN = function(seq, i) {cat('i\n'); sapply(tracks, function(trk) calc_seq_track_val(trk, seq[[i]]))}, mc.cores = 52, seq = sc)

vals_all_df = do.call('rbind', vals_all)

readr::write_tsv(as.data.frame(vals_all_df), './data/scrContr_preFilt_motifs.tsv')

# rownames(vals_all_df) = 

# vals_all_df

# length(vals_all)

# head(vals_all)

# vals_all[[1]]

# motif_vals_rc = sapply(1:(nchar(sc[[1]]) - nrow(motif_mats_log[[1]])), 
# # motif_vals = sapply(1:10, 
#        function(i) run_motif_on_window_rc(motif_mats_log[[1]], 
#                                         substr(sc[[1]], i, i+nrow(t(motif_mats_log[[1]]))-1)
#                                       )
#                     )

# l = length(motif_vals_rc)
# grp_vec = c(sort(rep(1:floor(l/10), 10)), rep(ceiling(l/10), l %% 10))

# grp_vec

# length(grp_vec)

# length(motif_vals_rc)

# mv_mean = tapply(motif_vals_rc, grp_vec, mean)


# mv_lse = log(sum(exp(tapply(motif_vals_rc, grp_vec, mean))))

# mv_lse

# length(sort(rep(1:ceiling(nchar(sc[[1]])/10), 10)))

# a



