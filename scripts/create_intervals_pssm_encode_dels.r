library(dplyr)
library(misha)
# library(parallel)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')
SEED = 1337
set.seed(SEED)

jolma_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jolma.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
motifs_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/homer.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
jaspar_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
cis_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/cis-bp.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
jolma_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jolma.key', header = F, col.names = c('key', 'track', 'dummy'))
motifs_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/homer.key', header = F, col.names = c('key', 'track', 'dummy'))
jaspar_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.key', header = F, col.names = c('key', 'track', 'dummy'))
cis_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/cis-bp.key', header = F, col.names = c('key', 'track', 'dummy'))
kmers_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/kmers.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
kmers_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/kmers.key', header = F, col.names = c('key', 'track', 'dummy'))

pmc_inp = gintervals.load('ENCODE_SCREEN_ccREs_dELS')
mids = round(apply(pmc_inp[,c('start', 'end')],1,mean))
pmc_inp$start = mids-250
pmc_inp$end = mids+250
# pmc = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')
# pmc_rn = as.data.frame(do.call('rbind', sapply(rownames(pmc), stringr::str_split, '-')))
# pmc_rn[,2:3] = apply(pmc_rn[,2:3], 2, as.numeric)
# colnames(pmc_rn) = c('chrom', 'start','end')
# encode_dels$enh_id = 1:nrow(encode_dels)
# dels_nei_pmc = gintervals.neighbors(encode_dels, pmc_rn, maxdist = 0, maxneighbors=5)
# pmc_inp = encode_dels[!(1:nrow(encode_dels) %in% dels_nei_pmc$enh_id),1:3]

print(table(abs(apply(pmc_inp[,c('start', 'end')], 1, diff))))

chroms_missing = ALLGENOME[[1]]$chrom[!(ALLGENOME[[1]]$chrom %in% unique(pmc_inp$chrom))]
fake_seqs = do.call('rbind', lapply(chroms_missing, function(chrom) {
    coord = round(ALLGENOME[[1]]$end[match(chrom, ALLGENOME[[1]]$chrom)]/2)
    return(as.data.frame(list('chrom' = chrom, 'start' = coord - 5e+2, 'end' = coord + 5e+2)))
}))
if (nrow(fake_seqs) > 0) {
    pmc_inp = rbind(pmc_inp, fake_seqs)
}

print(head(fake_seqs))

# trks_motifs = mapply(list(cis_key, jaspar_key, jolma_key, motifs_key), c('jolma', 'jaspar', 'motifs', 'cis_bp'), FUN = function(x,y) {
#     return(sapply(x$track, function(z) paste0(y, '.', z)))
# }) %>% unlist

trks_motifs = mapply(list(kmers_key), c('kmers'), FUN = function(x,y) {
    return(sapply(x$track, function(z) paste0(y, '.', z)))
})
# print(head(trks_motifs))
# print(grep('-|\\$', trks_motifs, v=T))
# trks_motifs = gsub('-|\\$', '_', trks_motifs)
# print(head(trks_motifs))
gdb.reload()

# print(head(pmc_inp))
# print(tail(pmc_inp))
create_intervals_pssm = function(trks) {
    ## Create tracks
    sapply(seq_along(trks), function(trks, i) {
        if (i%%10 == 0) {cat(paste0(i/length(trks), '\n'))}
        pssmset = stringr::str_split(trks[[i]], '\\.') %>% purrr::map(1) %>% unlist
        trk_clean = stringr::str_split(trks[[i]], '\\.') %>% purrr::map(2) %>% unlist
        # print(pssmset)
        # print(pssmset)
        
        if (pssmset == 'jolma') {pssmid = jolma_key$key[jolma_key$track == trk_clean]; pssmset_fix = pssmset}
        else if (pssmset == 'motifs') {pssmid = motifs_key$key[motifs_key$track == trk_clean]; pssmset_fix = 'homer'}
        else if (pssmset == 'jaspar') {pssmid = jaspar_key$key[jaspar_key$track == trk_clean]; pssmset_fix = pssmset}
        else if (pssmset == 'cis_bp') {pssmid = cis_key$key[cis_key$track == trk_clean]; pssmset_fix = 'cis-bp'}
        else if (pssmset == 'kmers') {pssmid = kmers_key$key[kmers_key$track == trk_clean]; pssmset_fix = pssmset}
        # print(pssmid)
        trknm = paste0(pssmset, '.', trk_clean, '_500bp')
        trknm = gsub('-|\\$', '_', trknm)
        # print(pssmset_fix)
        # print(trknm)
        # if(!is.numeric(pssmid)) {
        
            # }
        if (!gtrack.exists(trknm) && is.numeric(pssmid)) {
            # print(trks[[i]])
            # print(pssmset)
            # print(pssmid)
            gtrack.create_pwm_energy(track = trknm, description = 'tmp', pssmset = pssmset_fix, pssmid = pssmid, prior = 0.01, iterator = pmc_inp)
        }
    }, trks = trks)
    ## Extract tracks
    # a = gextract(paste0(trks, '_tmp'), intervals=gintervals.all(), iterator=pmc_inp)
#     a = gextract(paste0(trks, '_tmp'), intervals=pmc_inp, iterator=1)
    # sapply(paste0(trks, '_tmp'), gtrack.rm, f=T)
    # return(a)
    gdb.reload()
}

aa = create_intervals_pssm(trks_motifs)

trks_motifs = gtrack.ls('kmers.*_500bp')
trk_inds = sort(1+(1:length(trks_motifs))%%4)
aa_lst = lapply(sort(unique(trk_inds)), function(u) gextract(trks_motifs[trk_inds == u], intervals=pmc_inp, iterator=pmc_inp))
aa_all = plyr::join_all(aa_lst, by=c('chrom', 'start', 'end'), type='left')
# saveRDS(aa_all, './data/encode_dels_all_motifs_500bp.rds')
saveRDS(aa_all, './data/encode_dels_kmers_motifs_500bp.rds')
# if (length(gtrack.ls('_500bp')) > 0) {gugu = sapply(gtrack.ls('_500bp'), gtrack.rm, f=T)}

