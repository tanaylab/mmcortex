library(dplyr)
library(misha)
# library(parallel)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
gsetroot('/home/aviezerl/mm10')
# options(gmultitasking = FALSE)

SUFFIX = '_tmp_200bp'

# trks_motifs = grep('10bp', unlist(sapply(c('motifs', 'jaspar', 'jolma', 'cis_bp'), gtrack.ls)), inv=T, v=T)

jolma_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jolma.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
motifs_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/homer.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
jaspar_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))
cis_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/cis-bp.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))


kmers_data = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/kmers.data', header = F, col.names = c('key', 'pos', 'A', 'C', 'G', 'T'))

jolma_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jolma.key', header = F, col.names = c('key', 'track', 'dummy'))
motifs_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/homer.key', header = F, col.names = c('key', 'track', 'dummy'))
jaspar_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/jaspar.key', header = F, col.names = c('key', 'track', 'dummy'))
cis_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/cis-bp.key', header = F, col.names = c('key', 'track', 'dummy'))

kmers_key = read.delim('/net//mraid14//export//tgdata//db//tgdb//mm10/trackdb/pssms/kmers.key', header = F, col.names = c('key', 'track', 'dummy'))


pmc = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')

pmc_rn = as.data.frame(do.call('rbind', sapply(rownames(pmc), stringr::str_split, '-')))

colnames(pmc_rn) = c('chrom', 'start', 'end')
pmc_rn[,2:3] = apply(pmc_rn[,2:3], 2, as.numeric)

pmc_int = gintervals(pmc_rn[,1], pmc_rn[,2], pmc_rn[,3])
chroms_missing = ALLGENOME[[1]]$chrom[!(ALLGENOME[[1]]$chrom %in% unique(pmc_int$chrom))]
fake_seqs = do.call('rbind', lapply(chroms_missing, function(chrom) {
    coord = round(ALLGENOME[[1]]$end[match(chrom, ALLGENOME[[1]]$chrom)]/2)
    return(as.data.frame(list('chrom' = chrom, 'start' = coord - 5e+2, 'end' = coord + 5e+2)))
}))
print(fake_seqs)
pmc_inp = rbind(pmc_int, fake_seqs)
mids = round((pmc_inp$start + pmc_inp$end)/2)
pmc_inp$start = mids - 100
pmc_inp$end = mids + 100

trks_motifs = mapply(list(cis_key, jaspar_key, jolma_key, motifs_key), c('cis_bp', 'jaspar', 'jolma', 'motifs'), FUN = function(x,y) {
    return(sapply(x$track, function(z) paste0(y, '.', z)))
})

# trks_motifs = mapply(list(kmers_key), c('kmers'), FUN = function(x,y) {
#     return(sapply(x$track, function(z) paste0(y, '.', z)))
# })
gdb.reload()

print('start creating tracks')
create_intervals_pssm = function(trks) {
    ## Create tracks
    sapply(seq_along(trks), function(trks, i) {
        if (i%%100 == 0) {message(paste0(i/length(trks), '\n'))}
        pssmset = stringr::str_split(trks[[i]], '\\.') %>% purrr::map(1) %>% unlist
        trk_clean = stringr::str_split(trks[[i]], '\\.') %>% purrr::map(2) %>% unlist
        # print(pssmset)
        if (pssmset == 'jolma') {pssmid = jolma_key$key[jolma_key$track == trk_clean]; pssmset_fix = pssmset}
        else if (pssmset == 'motifs') {pssmid = motifs_key$key[motifs_key$track == trk_clean]; pssmset_fix = 'homer'}
        else if (pssmset == 'jaspar') {pssmid = jaspar_key$key[jaspar_key$track == trk_clean]; pssmset_fix = pssmset}
        else if (pssmset == 'cis_bp') {pssmid = cis_key$key[cis_key$track == trk_clean]; pssmset_fix = 'cis-bp'}
        else if (pssmset == 'kmers') {pssmid = kmers_key$key[kmers_key$track == trk_clean]; pssmset_fix = pssmset}
        # print(pssmid)
        trknm = paste0(pssmset, '.', trk_clean, SUFFIX)
        trknm = gsub('-|\\$', '_', trknm)
        if (!gtrack.exists(trknm)) {
            gtrack.create_pwm_energy(track = trknm, description = 'tmp', pssmset = pssmset_fix, pssmid = pssmid, prior = 0.01, iterator = pmc_inp)
        }
    }, trks = trks)
    ## Extract tracks
    gdb.reload()
    
#     a = gextract(paste0(trks, '_tmp'), intervals=pmc_inp, iterator=1)
    # sapply(paste0(trks, '_tmp'), gtrack.rm, f=T)
    # return(a)
}

aa = create_intervals_pssm(unlist(trks_motifs))
# trks_motifs = gtrack.ls('kmers.*_tmp')
# trk_inds = sort(1+(1:length(trks_motifs))%%4)
trk_inds = unlist(do.call('c', lapply(1:4, function(i) 
                    if (i!=4) {
                        rep(i, ceiling(length(trks_motifs)/4))
                    } 
                    else {
                        rep(i, floor(length(trks_motifs)/4))
                    }
)))
trks_motifs = gtrack.ls('SUFFIX')
aa_lst = lapply(sort(unique(trk_inds)), function(u) gextract(trks_motifs[trk_inds == u], intervals=pmc_inp, iterator=pmc_inp))
aa_all = plyr::join_all(aa_lst, by=c('chrom', 'start', 'end'), type='left')
# a = gextract(gtrack.ls('kmers.*_tmp'), intervals=pmc_inp, iterator=pmc_inp)

# saveRDS(a, './data/peak_mc_all_motifs.rds')
# saveRDS(a, './data/peak_mc_kmers_motifs.rds')
saveRDS(aa_all, './data/peak_mc_motifs_200bp.rds')