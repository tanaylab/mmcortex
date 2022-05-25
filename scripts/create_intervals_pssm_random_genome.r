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

N_CONTROL = 1e+5

chrom_lens = apply(ALLGENOME[[1]][,2:3], 1, diff)
chrom_fracs = setNames(chrom_lens/sum(chrom_lens), ALLGENOME[[1]][,1])
# chrom_fracs[order(as.numeric(gsub('chr', '', names(chrom_fracs))))]
# chrom_fracs

sample_seqs = mapply(chrom_fracs, names(chrom_fracs), chrom_lens, FUN = function(x, y, z) {
    chrom = rep(y, round(x*N_CONTROL)); start = sample.int(n = z, size =length(chrom)); end = start + 500;
    return(as.data.frame(rbind(chrom, start, end)))
})

sample_seqs = as.data.frame(do.call('rbind', lapply(sample_seqs, t)))
sample_seqs[,2:3] = apply(sample_seqs[,2:3], 2, as.numeric)
sample_seqs = sample_seqs[with(sample_seqs, order(chrom, start, end)),]
sample_seqs = sample_seqs[sample_seqs$start >= 3e+6,]
inds_remove = which(sample_seqs$end > ALLGENOME[[1]]$end[match(sample_seqs$chrom, ALLGENOME[[1]]$chrom)])
sample_seqs = sample_seqs[!(1:nrow(sample_seqs) %in% inds_remove),]

pmc_inp = sample_seqs

chroms_missing = ALLGENOME[[1]]$chrom[!(ALLGENOME[[1]]$chrom %in% unique(pmc_inp$chrom))]
fake_seqs = do.call('rbind', lapply(chroms_missing, function(chrom) {
    coord = round(ALLGENOME[[1]]$end[match(chrom, ALLGENOME[[1]]$chrom)]/2)
    return(as.data.frame(list('chrom' = chrom, 'start' = coord - 5e+2, 'end' = coord + 5e+2)))
}))
if (nrow(fake_seqs) > 0) {
    pmc_inp = rbind(pmc_inp, fake_seqs)
}
mids = round((pmc_inp$start + pmc_inp$end)/2)
pmc_inp$start = mids - 100
pmc_inp$end = mids + 100
print(head(fake_seqs))

trks_motifs = mapply(list(cis_key, jaspar_key, jolma_key, motifs_key), c('cis_bp', 'jaspar', 'jolma', 'motifs'), FUN = function(x,y) {
    return(sapply(x$track, function(z) paste0(y, '.', z)))
}) %>% unlist
print(head(trks_motifs))

# trks_motifs = mapply(list(kmers_key), c('kmers'), FUN = function(x,y) {
#     return(sapply(x$track, function(z) paste0(y, '.', z)))
# })
# print(grep('-|\\$', trks_motifs, v=T))
# trks_motifs = gsub('-|\\$', '_', trks_motifs)
# # print(head(trks_motifs))
gdb.reload()

# print(head(pmc_inp))
# print(tail(pmc_inp))
create_intervals_pssm = function(trks) {
    ## Create tracks
    sapply(seq_along(trks), function(trks, i) {
        if (i%%100 == 0) {cat(paste0(i/length(trks), '\n'))}
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
        trknm = paste0(pssmset, '.', trk_clean, '_rg_200bp')
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

# trks_motifs = gtrack.ls('kmers.*_rand_gen')
# trk_inds = sort(1+(1:length(trks_motifs))%%4)
trk_inds = unlist(do.call('c', lapply(1:4, function(i) 
                    if (i!=4) {
                        rep(i, ceiling(length(trks_motifs)/4))
                    } 
                    else {
                        rep(i, floor(length(trks_motifs)/4))
                    }
)))
trks_motifs = gtrack.ls('_rg_200bp')
aa_lst = lapply(sort(unique(trk_inds)), function(u) gextract(trks_motifs[trk_inds == u], intervals=pmc_inp, iterator=pmc_inp))
aa_all = plyr::join_all(aa_lst, by=c('chrom', 'start', 'end'), type='left')
# aaa = gextract(trks_motifs, intervals=gintervals.all(), iterator=pmc_inp)
# saveRDS(aa_all, './data/random_genome_100k_seqs_all_motifs.rds')
# saveRDS(aa_all, './data/random_genome_100k_seqs_kmers_motifs.rds')
saveRDS(aa_all, './data/random_genome_100k_seqs_motifs_200bp.rds')
# if (length(gtrack.ls('tmp')) > 0) {gugu = sapply(gtrack.ls('tmp'), gtrack.rm, f=T)}