wd = '/net//mraid14//export//tgdata/users/yonshap/data_other/TF_motif/cis-BP/'
setwd(wd)

pwm_files = list.files('./pwms_all_motifs')

tf_info = readr::read_tsv('./TF_Information.txt')

#tf_info_plus = readr::read_tsv('./TF_Information_all_motifs.txt')

# tf_info_plus = tf_info

# tfoi = c('Bcl11b', sort(unlist(read.table('/net//mraid14//export//tgdata/users/yonshap/proj/mmcortex/data/tfoi.txt'))))
# #tfoi
# tf_info_plus$uniqueID = apply(tf_info_plus[,c('TF_Name', 'DBID...14')], 1, stringr::str_c, collapse = '_')

# tf_info_plus = tf_info_plus[tf_info_plus$TF_Name %in% tfoi & tf_info_plus$Motif_ID != '.',]

# tf_info_plus = tf_info_plus[!duplicated(tf_info_plus$uniqueID),]

# tf_info = tf_info_plus
#tf_info$uniqueID = apply(tf_info[,c('TF_Name', 'DBID...14')], 1, stringr::str_c, collapse = '_')

#tf_info = tf_info[tf_info$TF_Name %in% tfoi & tf_info$Motif_ID != '.',]


out_key = '/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.key'
out_data = '/net/mraid14/export/tgdata/db/tgdb/mm10/trackdb/pssms/cis-bp.data'

#readr::write_tsv(data.frame(do.call('rbind', lapply(seq_along(tf_info$uniqueID), function(uid, i) {
#	c(i, gsub('\\.|-|\\$', '_', uid[[i]]), 1)}, uid = tf_info$uniqueID))), file = out_key, col_names = F)

cis_bp_data = do.call('rbind', lapply(seq_along(tf_info$uniqueID), function(uid, i) {
    mid = tf_info$Motif_ID[tf_info$uniqueID == uid[[i]]]
#    print(mid)
    mot = read.table(file = paste0('./pwms_all_motifs/', mid, '.txt'), header = T)
#    print(head(mot))
    
#     print(head(mot))
    if (nrow(mot) > 0) {
        mot = apply(mot, 2, as.numeric)
        mot = round(mot, 3)
        mot[,1] = mot[,1]-1
        mot = data.frame(cbind(rep(i, nrow(mot)), mot))
    #     print(mot)
    #     write(x = paste(i, uid[[i]], 1, collapse = '\t'), file = out_key, append = T)
    #     readr::write_tsv(x = mot, file= out_data, append = T)
        return(mot)
    }
}, uid = tf_info$uniqueID
                            )
                      )
head(cis_bp_data)
readr::write_tsv(x = cis_bp_data, file= out_data, append = F, col_names = F)

cis_bp_key = data.frame(do.call('rbind', lapply(seq_along(tf_info$uniqueID), function(uid, i) {c(i, gsub('\\.|-', '_', uid[[i]]), 1)}, uid = tf_info$uniqueID)))
head(cis_bp_key)
cis_bp_key = cis_bp_key[as.numeric(cis_bp_key[,1]) %in% cis_bp_data[,1],]

readr::write_tsv(cis_bp_key, file = out_key, col_names = F)
key_file = readr::read_tsv(out_key,col_names = c('key', 'track', 'dummy'))
print(head(key_file))

data_file = readr::read_tsv(out_data, col_names = c('key', 'pos', 'A', 'C', 'G', 'T'))
print(head(data_file))

library(misha)
gsetroot('/home/aviezerl/mm10')
gdb.reload()
options(gmax.data.size = 1e+09, gmultitasking=T)
options(try.outFile = stdout()) 

tracks_dir = file.path(GROOT, 'tracks/cis_bp')
if (!dir.exists(tracks_dir)) {dir.create(tracks_dir)}

tracks_dir_10bp = file.path(GROOT, 'tracks/cis_bp_10bp')
if (!dir.exists(tracks_dir_10bp)) {dir.create(tracks_dir_10bp)}


#print(grep('eomes|fezf2', key_file$track, v=T, ign=T))
#print(key_file$track[grep('eomes|fezf2', key_file$track, ign=T)])
#purrr::walk(key_file$key[grep('eomes|fezf2', key_file$track, ign=T)], function(kf, df, i) {
purrr::walk(key_file$key, function(kf, df, i) {
    ind = which(key_file$key == i)
    if (!gtrack.exists(track = paste0('cis_bp.', kf[ind,2]))) {
       print(paste0('creating track ', kf[ind,2]))
       try(gtrack.create_pwm_energy(track = paste0('cis_bp.', kf[ind,2]), description = paste0('CIS-BP motif data for (TF_Name)_(Motif_ID): ', kf[ind,2]), prior = 0.01, iterator = 1,
                                    pssmset = 'cis-bp', pssmid = i))
	}
    if (!gtrack.exists(track = paste0('cis_bp_10bp.', kf[ind,2]))) {
        print(paste0('creating track ', kf[ind,2], ' -- 10bp'))
        try(gtrack.create_pwm_energy(track = paste0('cis_bp_10bp.', kf[ind,2]), description = paste0('CIS-BP motif data for (TF_Name)_(Motif_ID): ', kf[ind,2]), iterator = 10,prior = 0.01,
                                    pssmset = 'cis-bp', pssmid = i))
        }
#     )
}, kf = key_file, df = data_file)}
gdb.reload()
