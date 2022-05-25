ls()

library(metacell)
library(dbscan)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

MAX_DIST = 2.5e+06
NUM_TARGET = 2e+4
NUM_FINAL = 12000
EXT = 133

mc_rna = scdb_mc('pl_cort')

chain_path = '/home/aviezerl/proj/ebdnmt/rawdata/import/Weber_Nature_Communication_2020/mm10ToMm9.over.chain.fixed1'

mcmd = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN/CfuPN',
                'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))
cust_mc_ord_st_ord_md = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))[order(mcmd$mean_day[mcmd$st == st])]
                                      )
                                      )

# mat = scdb_mat('pl_cort')


library(pheatmap)

color_key = unique(mcmd[,c('st', 'color')])
md_clrmp = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100)
# col_annot = data.frame(st = mcmd$st, md = round(mcmd$mean_day))
col_annot = data.frame(st = mcmd$st, md = mcmd$mean_day)
rownames(col_annot) = mcmd$mc
# ann_colors = list(st = setNames(color_key$color, color_key$st), md = setNames(md_clrmp, 13:18))
ann_colors = list(st = setNames(color_key$color, color_key$st), md = md_clrmp)
ann_colors[['st']] = ann_colors[['st']][order(match(names(ann_colors[['st']]), cust_st_ord))]

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pmc = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')

pmc_rn = as.data.frame(do.call('rbind', stringr::str_split(rownames(pmc), '-')))
colnames(pmc_rn) = c('chrom', 'start', 'end')
pmc_rn[,c('start', 'end')] = apply(pmc_rn[,c('start', 'end')], 2, as.numeric)
head(pmc_rn)

prox_enh = unlist(read.delim('./data/mmcortex_proximal_enhancers_mm10.txt', h = F))
dist_enh = unlist(read.delim('./data/mmcortex_distal_enhancers_mm10.txt', h = F))

pmc_rn = pmc_rn[rownames(pmc) %in% sort(c(prox_enh, dist_enh)),]
pmc = pmc[rownames(pmc) %in% sort(c(prox_enh, dist_enh)),]


nsc_mc = mcmd$mc[mcmd$st == 'NSC']
ipc_mc = mcmd$mc[grep('IPC', mcmd$st)]
mat_mc = mcmd$mc[mcmd$st %in% c('CPN_L2-3','SCPN','CthPN','CPN_L5_6')]
unique(mcmd$st[mat_mc])

bulk_st = c(setNames(rep('NSC', length(nsc_mc)), nsc_mc),
            setNames(rep('IPC', length(ipc_mc)), ipc_mc),
            setNames(rep('Mature', length(mat_mc)), mat_mc))

bulk_st[as.character((1:nrow(mcmd))[!(1:nrow(mcmd)) %in% as.numeric(names(bulk_st))])] = 'Other'

bulk_st = bulk_st[order(as.numeric(names(bulk_st)))]

head(bulk_st)

pmc_norm = t(apply(pmc, 1, function(x) x/sum(x)))

pmc_norm_st = t(tgs_matrix_tapply(pmc_norm, bulk_st, mean))

pmc_norm_st_diff12 = apply(apply(pmc_norm_st, 1, function(x) x[head(order(x, decreasing = T), 2)]), 2, function(y) diff(rev(y))/mean(y))

norm_st_maxs = colnames(pmc_norm_st)[apply(pmc_norm_st[pmc_norm_st_diff12 >= 0.33,], 1, which.max)]
enh_nsc = which(norm_st_maxs == 'NSC')
enh_ipc = which(norm_st_maxs == 'IPC')
enh_mat = which(norm_st_maxs == 'Mature')

table(colnames(pmc_norm_st)[apply(pmc_norm_st[pmc_norm_st_diff12 >= 0.33,], 1, which.max)])

length(enh_nsc)
length(enh_ipc)
length(enh_mat)

lst_enh = list('NSC' = enh_nsc, 'IPC' = enh_ipc, 'Mature neurons' = enh_mat)
lst_mc = list('NSC' = nsc_mc, 'IPC' = ipc_mc, 'Mature neurons' = mat_mc)

eps = 1e-08

enh_cor_md = lapply(1:3, function(le, lm, i) {
    md = mcmd$mean_day[lm[[i]]]
    cor_md_enh = tgs_cor(t(pmc_norm[le[[i]],lm[[i]]]), as.matrix(md), spearman = T)
    ord_cor = order(cor_md_enh**2, decreasing = T)
    min_max_rat = log2((apply(pmc_norm[le[[i]],lm[[i]]], 1, max) + eps)/(apply(pmc_norm[le[[i]],lm[[i]]], 1, min) + eps))
    ret_mat = data.frame(do.call('cbind', list('enh' = rownames(pmc)[le[[i]][ord_cor]], 
                                               'seq_name' = paste0(names(le)[[i]], '_', 1:length(ord_cor)),
                                               'cor_md' = cor_md_enh[ord_cor], 
                                               'min_max_rat' = min_max_rat[ord_cor])))
    ret_mat[,3:ncol(ret_mat)] = apply(ret_mat[,3:ncol(ret_mat)], 2, as.numeric)
#     v = ret_mat[,'enh']
#     ret_mat = apply(subset(ret_mat, select = -enh), 2, as.numeric)
#     rownames(ret_mat) = v
    return(ret_mat)
},  le = lst_enh, lm = lst_mc)

COR_THRESH = 0.6
diff_enh_per_st = do.call('cbind', lapply(enh_cor_md, function(x) list('asc' = length(which(x[,'cor_md'] >= COR_THRESH)),
                                    'desc' = length(which(x[,'cor_md'] <= -COR_THRESH)))))
colnames(diff_enh_per_st) = c('NSC', 'IPC', 'Mature')

diff_enh_per_st

# rat = NUM_TARGET/sum(as.numeric(diff_enh_per_st))

# num_seqs_to_take = round(apply(diff_enh_per_st, 2, as.numeric)*rat)
# num_seqs_to_take = rbind(num_seqs_to_take, rep(50, ncol(num_seqs_to_take)))
# # num_seqs_to_take['unc',] = 50
# rownames(num_seqs_to_take) = c('asc', 'desc', 'unc')

# num_seqs_to_take = apply(num_seqs_to_take, 2, lapply, identity)
# num_seqs_to_take

# num_seqs_to_take = list('NSC' = list('asc' = 600, 'desc' = 600, 'unc' = 50),
#                        'IPC' = list('asc' = 600, 'desc' = 600, 'unc' = 50),
#                        'Mature' = list('asc' = 600, 'desc' = 600, 'unc' = 50))

num_seqs_to_take = apply(diff_enh_per_st, 2, function(x) {x[x > 900] = 900; return(as.list(x))})
num_seqs_to_take = lapply(num_seqs_to_take, function(x) {x$unc = 3500 - x$asc - x$desc; return(x)})

num_seqs_to_take

seqs_to_take = lapply(seq_along(num_seqs_to_take), function(x,n,i) {
    seqs_asc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'cor_md'], decreasing = T)], x[[i]]$asc), 
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = T),'seq_name'], x[[i]]$asc)) 
    seqs_desc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'cor_md'], decreasing = F)], x[[i]]$desc),
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = F),'seq_name'], x[[i]]$desc)) 
    seqs_unc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'min_max_rat'], decreasing = F)], x[[i]]$unc),
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = F),'seq_name'], x[[i]]$unc))
    seqs_unc = seqs_unc[!(seqs_unc %in% c(seqs_desc, seqs_asc))]
#     names(seqs_asc) = 1:length(seqs_asc)
#     names(seqs_desc) = 1:length(seqs_desc)
#     names(seqs_unc) = 1:length(seqs_unc)
    return(list('seqs_asc' = seqs_asc, 'seqs_desc' = seqs_desc, 'seqs_unc' = seqs_unc))
}, x = num_seqs_to_take,n = names(num_seqs_to_take))

names(seqs_to_take) = names(num_seqs_to_take)

seqs_unlist = do.call('c', lapply(seqs_to_take, unlist))

seqs_ints = data.frame(do.call('rbind', stringr::str_split(seqs_unlist, '-')))

seqs_ints$seq_name = names(seqs_unlist)

seqs_ints = seqs_ints[!duplicated(seqs_ints[,1:3]),]

colnames(seqs_ints) = c('chrom', 'start', 'end', 'seq_name')
seqs_ints[,2:3] = apply(seqs_ints[,2:3], 2, as.numeric)
head(seqs_ints)

legc = log2(mc_rna@e_gc + 1e-07)

hi_exp_genes = sapply(list(nsc_mc, ipc_mc, mat_mc), function(mci) apply(legc[,mci], 1, median))

head(hi_exp_genes[order(apply(hi_exp_genes, 1, max) - apply(hi_exp_genes, 1, min), decreasing = T),])

colnames(hi_exp_genes) = c('NSC', 'IPC', 'Mature')

hi_genes_st = lapply(seq_along(1:ncol(hi_exp_genes)), function(i) rownames(hi_exp_genes)[which(hi_exp_genes[,i] >= -14)])

intersect_all = sort(unique(unlist(hi_genes_st)))

gsetroot('/home/aviezerl/mm10')

# gsetroot('/home/aviezerl/mm9')

# seqs_ints_mm9 = gintervals.liftover(intervals = seqs_ints, chain = chain_path)

# seqs_ints_mm9$seq_name = seqs_ints[seqs_ints_mm9$intervalID, 'seq_name']
# seqs_ints_mm9 = seqs_ints_mm9[!duplicated(seqs_ints_mm9$seq_name),]

# seqs_ints_mm9_mids = round(0.5*(seqs_ints_mm9$start + seqs_ints_mm9$end))
# seqs_ints_mm9$start = seqs_ints_mm9_mids - 133
# seqs_ints_mm9$end = seqs_ints_mm9_mids + 133

# seqs_ints_mm9 = seqs_ints_mm9[with(seqs_ints_mm9, order(chrom, start, end)),]

seqs_ints_mids = round(0.5*(seqs_ints$start + seqs_ints$end))
seqs_ints$start = seqs_ints_mids - 133
seqs_ints$end = seqs_ints_mids + 133


# tss = gintervals.load('intervs.global.tss')

tss = gintervals.load('tss')

tss = tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]

hi_exp_tss = tss[tss$geneSymbol %in% intersect_all,]

trk_s = unlist(read.delim('./data/motif_tracks_selected_k=64_n=64.txt'))

stringr::str_split(trk_s, '\\.') %>% purrr::map(2) %>% stringr::str_split('_') %>% purrr::map(1) %>% unlist %>% sort

# purrr::walk(trk_s, function(trk) gvtrack.create(paste0(trk, '_max'), trk, 'max'))

trk_all_raw = gextract(trk_s, intervals = seqs_ints, iterator = 10)
trk_all_raw = trk_all_raw[with(trk_all_raw, order(chrom, start, end)),]

cn = c('chrom', 'start', 'end', 'intervalID')

trk_all_avg = tgs_matrix_tapply(t(subset(trk_all_raw, select = -c(chrom, start, end, intervalID))), trk_all_raw$intervalID, function(x) log(sum(exp(x))))

seqs_motifs = cbind(seqs_ints, trk_all_avg)

# seqs_motifs = gextract(paste0(trk_s, '_max'), intervals = seqs_ints_mm9, iterator = seqs_ints_mm9)

seqs_motifs = seqs_motifs[,!duplicated(colnames(seqs_motifs))]

# seqs_motifs = dplyr::left_join(seqs_motifs, subset(seqs_ints_mm9, select = -intervalID), by=c('chrom', 'start', 'end'))

head(seqs_motifs)

nrow(seqs_motifs)

# saveRDS(seqs_motifs, glue::glue('./data/seqs_ints_motif_trk_s_{NUM_TARGET}.rds'))
# saveRDS(seqs_motifs, glue::glue('./data/seqs_ints_motif_trk_s_{NUM_TARGET}_mm10.rds'))
saveRDS(seqs_motifs, glue::glue('./data/seqs_ints_trk_s_900_per_cat_mm10.rds'))

# seqs_motifs = readRDS(glue::glue('./data/seqs_ints_motif_trk_s_{NUM_TARGET}.rds'))

# seqs_motifs = readRDS(glue::glue('./data/seqs_ints_motif_trk_s_{NUM_TARGET}_mm10.rds'))
seqs_motifs = readRDS(glue::glue('./data/seqs_ints_trk_s_900_per_cat_mm10.rds'))

nrow(seqs_motifs)

nrow(seqs_motifs)

# dels_motifs = readRDS('./data/motif_max_ENCODE_SCREEN_dELS.rds')

# dels_motifs = vroom::vroom('./data/motif_max_ENCODE_SCREEN_dELS_266bp.tsv')

dels_motifs = readRDS('./data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds')

dels_motifs$seq_name = 1:nrow(dels_motifs)

dels_motifs_nei_seqs = gintervals.neighbors(dels_motifs[,c('chrom', 'start', 'end', 'seq_name')], pmc_rn, mindist = -1e+03, maxdist = 1e+3, maxneighbors = 10)

# dels_motifs_nei_seqs = gintervals.neighbors(dels_motifs[,c('chrom', 'start', 'end', 'seq_name')], seqs_motifs[,1:3], mindist = -1e+03, maxdist = 1e+3, maxneighbors = 10)

dels_motifs = dplyr::anti_join(dels_motifs, dels_motifs_nei_seqs[,1:3])

dels_motifs = as.data.frame(dels_motifs)
# dels_motifs$intervalID = 1:nrow(dels_motifs)

rna_nei = gintervals.neighbors(dels_motifs, hi_exp_tss, maxneighbors = 1, mindist = 0, maxdist = MAX_DIST)

dels_motifs = dplyr::left_join(dels_motifs, subset(rna_nei, select = c(chrom, start, end, dist)), by = c('chrom', 'start', 'end'))

if (length(which(is.na(dels_motifs$dist))) > 0) {dels_motifs$dist[is.na(dels_motifs$dist)] = MAX_DIST}

# dels_nei_ints = gintervals.neighbors(dels_motifs[,1:3], seqs_motifs[,1:3], maxdist = 1e+5, maxneighbors = 2)

# ints_nei_dels = gintervals.neighbors(seqs_motifs[,1:3], dels_motifs[,1:3], maxdist = 5e+2, maxneighbors = 2)

# motifs_df = seqs_motifs
# cn = c('chrom', 'start', 'end', 'seq_name')
# seqs_nei_5k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 5e+3, maxneighbors = 50)
# seqs_nei_25k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 2.5e+4, maxneighbors = 50)
# seqs_nei_250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+4 + 1, maxdist = 2.5e+5, maxneighbors = 50)
# seqs_nei_1250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+5 + 1, maxdist = 1.5e+6, maxneighbors = 50)

# tbl_nei_5k = table(unique(seqs_nei_5k[,c('seq_name', 'geneSymbol')])$seq_name)
# tbl_nei_25k = table(unique(seqs_nei_25k[,c('seq_name', 'geneSymbol')])$seq_name)
# tbl_nei_250k = table(unique(seqs_nei_250k[,c('seq_name', 'geneSymbol')])$seq_name)
# tbl_nei_1250k = table(unique(seqs_nei_1250k[,c('seq_name', 'geneSymbol')])$seq_name)

# #     nei_genes_5k = 
# #     print(head(seqs_nei_25k))
# tbl_nei_5k = table(unique(seqs_nei_5k[,c('intervalID', 'geneSymbol')])$intervalID)
# tbl_nei_25k = table(unique(seqs_nei_25k[,c('intervalID', 'geneSymbol')])$intervalID)
# tbl_nei_250k = table(unique(seqs_nei_250k[,c('intervalID', 'geneSymbol')])$intervalID)
# tbl_nei_1250k = table(unique(seqs_nei_1250k[,c('intervalID', 'geneSymbol')])$intervalID)

# motifs_df$nei_5k = as.numeric(unlist(tbl_nei_5k[motifs_df$seq_name]))

# motifs_df$nei_5k[is.na(motifs_df$nei_5k)] = 0

# motifs_df$nei_25k = as.numeric(unlist(tbl_nei_25k[motifs_df$seq_name]))
# motifs_df$nei_25k[is.na(motifs_df$nei_25k)] = 0

# motifs_df$nei_250k = as.numeric(unlist(tbl_nei_250k[motifs_df$seq_name]))
# motifs_df$nei_250k[is.na(motifs_df$nei_250k)] = 0

# motifs_df$nei_1250k = as.numeric(unlist(tbl_nei_1250k[motifs_df$seq_name]))
# motifs_df$nei_1250k[is.na(motifs_df$nei_1250k)] = 0

add_tbl_nei = function(motifs_df) {
    cn = c('chrom', 'start', 'end', 'seq_name')
    seqs_nei_5k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 5e+3, maxneighbors = 50)
    seqs_nei_25k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 2.5e+4, maxneighbors = 50)
    seqs_nei_250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+4 + 1, maxdist = 2.5e+5, maxneighbors = 50)
    seqs_nei_1250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+5 + 1, maxdist = 1.5e+6, maxneighbors = 50)
    
#     nei_genes_5k = 
#     print(head(seqs_nei_25k))
    tbl_nei_5k = table(unique(seqs_nei_5k[,c('seq_name', 'geneSymbol')])$seq_name)
    tbl_nei_25k = table(unique(seqs_nei_25k[,c('seq_name', 'geneSymbol')])$seq_name)
    tbl_nei_250k = table(unique(seqs_nei_250k[,c('seq_name', 'geneSymbol')])$seq_name)
    tbl_nei_1250k = table(unique(seqs_nei_1250k[,c('seq_name', 'geneSymbol')])$seq_name)

    motifs_df$nei_5k = as.numeric(unlist(tbl_nei_5k[motifs_df$seq_name]))
    motifs_df$nei_5k[is.na(motifs_df$nei_5k)] = 0

    motifs_df$nei_25k = as.numeric(unlist(tbl_nei_25k[motifs_df$seq_name]))
    motifs_df$nei_25k[is.na(motifs_df$nei_25k)] = 0

    motifs_df$nei_250k = as.numeric(unlist(tbl_nei_250k[motifs_df$seq_name]))
    motifs_df$nei_250k[is.na(motifs_df$nei_250k)] = 0

    motifs_df$nei_1250k = as.numeric(unlist(tbl_nei_1250k[motifs_df$seq_name]))
    motifs_df$nei_1250k[is.na(motifs_df$nei_1250k)] = 0
    
    return(motifs_df)
}

seqs_motifs = add_tbl_nei(seqs_motifs)

dels_motifs = add_tbl_nei(dels_motifs)

# dels_inds = unique(unlist(sapply(colnames(seqs_motifs)[colnames(seqs_motifs) != 'seq_name'], function(x) grep(x, colnames(dels_motifs), ign=T, v=F)[[1]])))
# dels_inds = unique(unlist(sapply(colnames(seqs_motifs)[colnames(seqs_motifs) != 'seq_name'], function(x) grep(x, colnames(dels_motifs), ign=T, v=F)[[1]])))
dels_inds = unique(unlist(sapply(colnames(seqs_motifs), function(x) grep(x, colnames(dels_motifs), ign=T, v=F)[[1]])))

# dels_inds = unique(unlist(sapply(gsub('_max', '', colnames(seqs_motifs)[colnames(seqs_motifs) != 'seq_name']), function(x) grep(x, gsub('_max', '', colnames(dels_motifs)), ign=T, v=F)[[1]])))

dels_motifs = dels_motifs[,c(dels_inds, grep('dist', colnames(dels_motifs)))]

motif_ecdfs = apply(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist)), 2, function(x) ecdf(x))

seqs_q = sapply(1:ncol(subset(seqs_motifs, select = -c(chrom, start, end, seq_name))), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = subset(seqs_motifs, select = -c(chrom, start, end, seq_name)))

rownames(seqs_q) = gsub(' ', '', apply(seqs_motifs[,1:3], 1, paste0, collapse = '-'))

dels_q = sapply(1:ncol(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist))), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = subset(dels_motifs, select = -c(chrom, start, end, seq_name)))

rownames(dels_q) = gsub(' ', '', apply(dels_motifs[,1:3], 1, paste0, collapse = '-'))

seqs_m_z = apply(seqs_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

dels_m_z = apply(dels_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

head(seqs_m_z)

head(dels_m_z)

# seqs_m_z = apply(subset(seqs_motifs, select = -c(chrom, start, end, seq_name)), 2, function(x) {v = -log2(1-ecdf(x)(x)); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

# dels_m_z = apply(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist)), 2, function(x) {v = -log2(1-ecdf(x)(x)); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

# ncol(seqs_m_z)

# ncol(dels_m_z)

# colnames(seqs_m_z) == colnames(dels_m_z)

# colnames(seqs_m_z)[colnames(seqs_m_z) != colnames(dels_m_z)]
# colnames(dels_m_z)[colnames(seqs_m_z) != colnames(dels_m_z)]

knn_new = tgs_cor_knn(t(seqs_m_z), t(dels_m_z), knn = 100, spearman = F)

# knn_new$max_diff = apply(knn_new, 1, function(x) max(seqs_m_z[x[['col1']],] - dels_m_z[x[['col2']],]))

knn_new$max_diff = mapply(knn_new$col1, knn_new$col2, FUN = function(x, y) {
    vec1 = seqs_m_z[x,]; vec2 = dels_m_z[y,]; 
#     ord = order(vec1, decreasing = T);
    max(abs(vec1 - vec2))
    }
                          )

# max_ind = apply(knn_new, 1, function(x) {
#     ord = order(seqs_m_z[x[['col1']],], decreasing = T); 
#     return(which.max(seqs_m_z[x[['col1']],ord] - dels_m_z[x[['col2']],ord]))
# }
#                )

# table(max_ind)

# head(knn_new)

# hist(knn_new$max_diff[knn_new$rank == 1])

# hist(knn_new$max_diff)

# hist(knn_new$cor[knn_new$rank == 1])

# library(dbscan)

# knn_z = kNN(dels_m_z, k = 4, query = seqs_m_z)

# nrow(knn_z$id)

saveRDS(knn_new, './data/enh_dELS_log_sum_exp_cor_knn_log_q_{NUM_TARGET}_mm10.rds')

knn_new = readRDS('./data/enh_dELS_log_sum_exp_cor_knn_log_q_{NUM_TARGET}_mm10.rds')

# knn_new$max_diff = mapply(knn_new$col1, knn_new$col2, FUN = function(x, y) {
#     vec1 = seqs_m_z[x,]; vec2 = dels_m_z[y,]; 
# #     ord = order(vec1, decreasing = T);
#     max(abs(vec1 - vec2))
#     }
#                           )

# saveRDS(knn_z, './data/enh_dELS_log_sum_exp_knn_log_q_{NUM_TARGET}_mm10.rds')

# knn_z = readRDS('./data/enh_dELS_log_sum_exp_knn_log_q_{NUM_TARGET}_mm10.rds')

# knn_z = readRDS('./data/enh_dELS_log_sum_exp_cor_knn_log_q_{NUM_TARGET}_mm10.rds')

mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')

tads = readRDS('./data/npc_tads_egc.rds')

cn_tads = readRDS('./data/cn_tads_egc.rds')

(sum(apply(ALLGENOME[[1]][,2:3], 1, diff)) - sum(apply(tads[,2:3], 1, diff)))/sum(apply(ALLGENOME[[1]][,2:3], 1, diff))
sum(apply(ALLGENOME[[1]][,2:3], 1, diff)) - sum(apply(tads[,2:3], 1, diff))
(sum(apply(ALLGENOME[[1]][,2:3], 1, diff)) - sum(apply(cn_tads[,2:3], 1, diff)))/sum(apply(ALLGENOME[[1]][,2:3], 1, diff))
sum(apply(ALLGENOME[[1]][,2:3], 1, diff)) - sum(apply(cn_tads[,2:3], 1, diff))

# cn_tads = readRDS('./data/cn_tads_egc.rds')

FOLD_THRESH = 4

st_vec = setNames(mcmd$mc, mcmd$st)
names(st_vec)[st_vec %in% nsc_mc] = 'NSC'
names(st_vec)[st_vec %in% ipc_mc] = 'IPC'
names(st_vec)[st_vec %in% mat_mc] = 'Mature'
names(st_vec)[!(st_vec %in% c(nsc_mc, ipc_mc, mat_mc))] = 'NA'

tad_st = t(tgs_matrix_tapply(as.matrix(tads[,as.character(1:nrow(mcmd))]), names(st_vec), mean))

FOLD_THRESH = 4

tad_cf = apply(tad_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) {
    y = order(x, decreasing = T);
    ifelse(x[[y[[1]]]]/x[[y[[2]]]] >= FOLD_THRESH, c('NSC', 'IPC', 'Mature')[which.max(x)], NA)
})

tad_neg_cf = apply(tad_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) {
    y = order(x, decreasing = F);
    ifelse(x[[y[[1]]]]/x[[y[[2]]]] <= 1/FOLD_THRESH, c('NSC', 'IPC', 'Mature')[which.min(x)], NA)
})

diff_tads = tads[which(!is.na(tad_neg_cf)), c('chrom', 'start', 'end', 'tadID')]

diff_tads$type_silenced = tad_neg_cf[which(!is.na(tad_neg_cf))]

diff_tad_nei_dels = gintervals.neighbors(diff_tads, dels_motifs, maxneighbors = 10000, maxdist = 0)

pmc_rn$peak = rownames(pmc)

# sum_mc_atac_across_tads = function(tads, tss, pmc) {
#     tads_genes = gintervals.neighbors(tads, dplyr::select(tss, chrom, start, end, peak), maxneighbors = 2000, maxdist = 0)
#     tads_genes_vec = tapply(tads_genes$peak, tads_genes$tadID, c)
#     tads$peaks = tads_genes_vec[match(tads$tadID, as.numeric(names(tads_genes_vec)))]    
#     genes_tads_egc = t(sapply(tads$peaks, function(x) if(length(x) > 1) {return(apply(pmc[x,], 2, sum))} else if (length(x) == 1) {return(pmc[x,])} else {return(rep(NA, ncol(pmc)))}))
#     return(cbind(tads, genes_tads_egc))
# }

# tads_atac = sum_mc_atac_across_tads(dplyr::select(tads, chrom, start, end, tadID), pmc_rn, pmc)

# pmc_prom = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')

# pmc_prom = pmc_prom[!(rownames(pmc_prom) %in% sort(c(prox_enh, dist_enh))),]

# pmc_prom_rn = as.data.frame(do.call('rbind', stringr::str_split(rownames(pmc_prom), '-')))

# pmc_prom_rn[,2:3] = apply(pmc_prom_rn[,2:3], 2, as.numeric)

# pmc_prom_rn$peak = rownames(pmc_prom)

# colnames(pmc_prom_rn) = c('chrom', 'start', 'end', 'peak')

# tads_prom = sum_mc_atac_across_tads(dplyr::select(tads, chrom, start, end, tadID), pmc_prom_rn, pmc_prom)

# # cn_tads_atac = sum_mc_atac_across_tads(dplyr::select(cn_tads, chrom, start, end, tadID), pmc_rn, pmc)

# # saveRDS(tads_atac, './data/npc_tads_atac.rds')
# # saveRDS(cn_tads_atac, './data/cn_tads_atac.rds')

# tad_atac_cf = apply(tad_atac_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) {
# #     y = order(x, decreasing = T);
#     ifelse(max(x)/mean(x) >= FOLD_THRESH/2, c('NSC', 'IPC', 'Mature')[which.max(x)], NA)
# })

# tad_atac_neg_cf = apply(tad_atac_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) {
# #     y = order(x, decreasing = F);
#     ifelse(min(x)/mean(x) <= 2/FOLD_THRESH, c('NSC', 'IPC', 'Mature')[which.min(x)], NA)
# })

# tads_sig_atac = which(apply(tad_atac_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) max(x + 1e-07)/min(x + 1e-07)) >= 2)

# tads_sig_rna = which(apply(tad_st[,c('NSC', 'IPC', 'Mature')], 1, function(x) max(x + 1e-09)/min(x + 1e-09)) >= 4)

# tads_drop_na = tidyr::drop_na(tads)

# kw_res = lapply(1:nrow(tads_drop_na), function(i) kruskal.test(unlist(tads_drop_na[i,as.character(1:nrow(mcmd))]), g = names(st_vec)))

tads_order_rna = apply(tad_st[,c('NSC', 'IPC', 'Mature')], 1, order)

# tads_order_atac = apply(tad_atac_st[,c('NSC', 'IPC', 'Mature')], 1, order)

# min_max_tbl = rbind(table(colnames(tad_st)[unlist(apply(tad_st, 1, which.max))]),
#     table(colnames(tad_st)[unlist(apply(tad_st, 1, which.min))]),
#     table(colnames(tad_st)[unlist(apply(tad_atac_st, 1, which.max))]),
#     table(colnames(tad_st)[unlist(apply(tad_atac_st, 1, which.min))]))

# rownames(min_max_tbl) = c('max_in_RNA', 'min_in_RNA', 'max_in_ATAC', 'min_in_ATAC')

# min_max_tbl

# pheatmap(min_max_tbl, cluster_cols = F, cluster_rows = F, fontsize = 20)

# options(repr.plot.width = 10, repr.plot.height = 10)

tads = dplyr::mutate(tads, len = end - start) %>% dplyr::relocate(len, .before = Compartment)
# tad_atac_st = t(tgs_matrix_tapply(as.matrix(tads[,as.character(1:nrow(mcmd))]), names(st_vec), mean))


# tad_all = as.data.frame(cbind(tad_st[,c('NSC', 'IPC', 'Mature')], tad_atac_st[,c('NSC', 'IPC', 'Mature')]))

# tad_all = median(tads$len)*tad_all/tads$len
# na_inds = which(apply(tad_all, 1, function(x) any(is.na(x))))
# tad_all = t(tidyr::drop_na(tad_all))
# # tad_all[1:3,] = log2(tad_all[1:3,] + 1e-010)
                      
# rownames(tad_all)[1:3] = paste0(rownames(tad_all)[1:3], '_RNA')
# rownames(tad_all)[4:6] = paste0(rownames(tad_all)[4:6], '_ATAC')
# print(dim(tad_all))

# tad_all = t(apply(tad_all, 1, function(x) x/sum(x)))
# # tad_all = t(apply(tad_all, 1, function(x) ecdf(x)(x)))
# tad_all = t(apply(tad_all, 1, function(x) log10(x+1e-010)))
# # tad_all = t(apply(tad_all, 1, function(x) (x - mean(x))/sd(x)))
# # tad_all = t(apply(tad_all, 1, function(x) ifelse(x > 5, 5, x)))
# # tad_all = t(apply(tad_all, 1, function(x) ifelse(x < -5, -5, x)))

# # tad_all = tad_all[c('NSC_ATAC', 'NSC_RNA', 'IPC_ATAC', 'IPC_RNA', 'Mature_ATAC', 'Mature_RNA'),]

# options(repr.plot.width = 25, repr.plot.height = 10)

# p = pheatmap(tad_all, cluster_rows = F, show_colnames = F, fontsize = 20)
                  
# # save_pheatmap_png(p, './figs/tad_length_and_frac_norm.png', h = 800, w = 3200)

# # p = pheatmap(tad_all, cluster_rows = F, show_colnames = F, fontsize = 20)
# # save_pheatmap_png(p, './figs/tad_st_log_rna_log_atac_no_z.png', h = 800, w = 6400)

# # p = pheatmap(t(tad_all), show_rownames = F, fontsize = 20)

#                   # save_pheatmap_png(p, './figs/tad_st_rna_atc_hm_clipped_5_vert.png', h = 6400, w = 800)
                  
# options(repr.plot.width = 8, repr.plot.height = 8)


tads_na_inds = which(apply(tads[,as.character(1:nrow(mcmd))], 1, function(x) any(is.na(x))))

tads = tads[!apply(tads[,as.character(1:nrow(mcmd))], 1, function(x) any(is.na(x))),]

tad_max = log2(1e-7 + apply(tads[,as.character(1:nrow(mcmd))], 1, max))
tad_min = log2(1e-7 + apply(tads[,as.character(1:nrow(mcmd))], 1, min))

names(tad_min) = tads$tadID
names(tad_max) = tads$tadID

xcut = cut(tad_min, breaks = seq(min(tad_min)-1, max(tad_min), l=21))

hi_var_tad = unlist(lapply(levels(xcut), function(q) names(tad_min)[which(xcut == q & tad_max >= quantile(tad_max[xcut == q], 0.75))]))

const_lo_tad = unlist(lapply(levels(xcut)[1:11], function(q) names(tad_min)[which(xcut == q & tad_max <= quantile(tad_max[xcut == q], 0.15))]))
const_hi_tad = unlist(lapply(levels(xcut)[14:20], function(q) names(tad_min)[which(xcut == q & tad_max <= quantile(tad_max[xcut == q], 0.25))]))

tad_annot = c(setNames(rep('hi_var', length(hi_var_tad)), hi_var_tad),
             setNames(rep('const_lo', length(const_lo_tad)), const_lo_tad),
             setNames(rep('const_hi', length(const_hi_tad)), const_hi_tad))


tad_annot[as.character(names(tad_min)[!(names(tad_min) %in% as.numeric(names(tad_annot)))])] = NA
tad_annot = tad_annot[order(as.numeric(names(tad_annot)))]
tad_annot

marks = scdb_gset('pl_cort_marks_f')

marks = names(marks@gene_set)

tads_genes_vec = unlist(mapply(tads$genes, names(tads$genes), FUN = function(x, y) setNames(x, rep(y, length(x)))))

names(tads_genes_vec) = as.numeric(purrr::map(stringr::str_split(names(tads_genes_vec), '\\.'), 2))

tad_marks_tbl = table(names(tads_genes_vec)[tads_genes_vec %in% marks])

dir.create('./figs/marker_tads_by_st')
dir.create('./figs/marker_tads_by_time')

options(repr.plot.width = 18)
eps = 1e-07
sapply(as.numeric(names(tad_marks_tbl)), function(n) {x = tads[match(n, tads$tadID),-c(1:7)]; 
#                                                       x = log2((x + eps)/(min(x) + eps));
#                                                       ord = order(mcmd$mean_day)
                                                      coords = paste0(tads[match(n, tads$tadID),1:3], collapse = '-');
#                                                       ord = order(mcmd$color);
                                                      ord = cust_mc_ord_st_ord_md
                                                      x = as.numeric(x[ord]); 
                                                      gns = unlist(tads$genes[[match(n, tads$tadID)]]);
                                                      gns_mrks = gns[gns %in% marks];
#                                                       print(gns_mrks)
                                                      png(paste0('./figs/marker_tads_by_st/tad_', n, '_', paste0(gns_mrks, collapse='_'), '.png'), w = 1600, h = 800)
#                                                       png(paste0('./figs/marker_tads_by_time/tad_', n, '_', paste0(gns_mrks, collapse='_'), '.png'), w = 1600, h = 800)
                                                      par(mar = c(5,4,6,2))
                                                      barplot(x, col = mcmd$color[ord]);
#                                                       print(coords)
                                                      title(paste0(c(coords, paste0(c('All genes: ', paste0(gns[!(gns %in% marks)], collapse = ', '), 'Marker(s): ', paste0(gns_mrks, collapse = ', '))))));
                                                      dev.off()})
options(repr.plot.width = 8)

plot(tad_min, tad_max)
points(tad_min[hi_var_tad], tad_max[hi_var_tad], col = 'red')
points(tad_min[const_lo_tad], tad_max[const_lo_tad], col = 'blue')
points(tad_min[const_hi_tad], tad_max[const_hi_tad], col = 'green')
legend('bottomright', legend = c('Hi-var', 'Const. lo', 'Const. hi'), col = c('red', 'blue', 'green'), pch = c(1,1,1))

tad_st_all_st = t(tgs_matrix_tapply(as.matrix(tads[,as.character(1:nrow(mcmd))]), mcmd$st, mean))

pheatmap(log2(1e-09 + tad_st_all_st[which(tad_annot == 'hi_var'),]), breaks = seq(-20,-6,l=100))
pheatmap(log2(1e-09 + tad_st_all_st[which(tad_annot == 'const_lo'),]), breaks = seq(-20,-6,l=100))
pheatmap(log2(1e-09 + tad_st_all_st[which(tad_annot == 'const_hi'),]), breaks = seq(-20,-6,l=100))

tads_nei_dels = gintervals.neighbors(dplyr::mutate(dplyr::select(dels_motifs, chrom, start, end, seq_name), rowname = 1:nrow(dels_motifs)), dplyr::select(tads, chrom, start, end, tadID), maxdist=0)

# diff_tads_atac = tads_atac[which(!is.na(tad_atac_neg_cf)), c('chrom', 'start', 'end', 'tadID')]

# diff_tads_atac$type_silenced = tad_atac_neg_cf[which(!is.na(tad_atac_neg_cf))]

# diff_tad_atac_nei_dels = gintervals.neighbors(diff_tads_atac, dels_motifs, maxneighbors = 10000, maxdist = 0)

# length(diff_tad_atac_nei_dels$seq_name)

# length(unique(diff_tad_atac_nei_dels$seq_name))

nrow(knn_new)

sapply(3:8, function(j) length(which(knn_new$max_diff <= j)))

knn_new = knn_new[knn_new$max_diff <= 4,]

nrow(knn_new)

# knn_new = knn_new[!duplicated(knn_new$col2),]

knn_new$native_seq = seqs_motifs$seq_name[knn_new$col1]

dels_rn = gsub(' ', '', apply(dels_motifs[,1:3], 1, paste0, collapse = '-'))

knn_new$peak = dels_rn[knn_new$col2]

knn_new$tadID = tads_nei_dels$tadID[match(knn_new$col2, tads_nei_dels$rowname)]

knn_new$tad_annot = tad_annot[as.character(knn_new$tadID)]

head(knn_new)

chosen_ind_list = as.list(setNames(rep(NA, length(sort(unique(knn_new$col1)))), sort(unique(knn_new$col1))))

for (i in sort(unique(knn_new$col1))) {
    inds = which(knn_new$col1 == i & !is.na(knn_new$tad_annot))
    blklst = which(knn_new$col2[inds] %in% unlist(chosen_ind_list))
    inds = inds[!(1:length(inds) %in% blklst)]
    cands = tapply(knn_new$col2[inds], knn_new$tad_annot[inds], function(x) as.numeric(x[[1]]))
    if (length(cands) >= 2 & all(!(unlist(cands) %in% unlist(chosen_ind_list)))) {
        chosen_ind_list[[as.character(i)]] = cands
    }
#     print(cands)
}

saveRDS(chosen_ind_list, './data/mpra_chosen_ind_list.rds')

table(sapply(chosen_ind_list, length))

sum(table(sapply(chosen_ind_list, length)))

length(which(sapply(chosen_ind_list, length) > 1))

seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) paste0(x[1:2], collapse = '_')) %>% table %>% unlist
        

seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) x[[1]]) %>% table
        

seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) x[[1]]) %>% table %>% sum
        

# chosen_inds = sapply(sort(unique(knn_new$col1)), function(i) {
# #     inds = which(knn_new$col1 == i & !is.na(knn_new$tad_order_rna))
# #     return(tapply(knn_new$col2[inds], knn_new$tad_order_rna[inds], function(x) x[[1]]))
#     inds = which(knn_new$col1 == i & !is.na(knn_new$tad_annot))
#     return(tapply(knn_new$col2[inds], knn_new$tad_annot[inds], function(x) as.numeric(x[1:2])))
# }
#                   )

# names(chosen_inds) = seqs_motifs$seq_name[sort(unique(knn_new$col1))]

# ci_len = lapply(chosen_inds, length)

# ci2_tads = lapply(chosen_inds[ci_len == 2], function(x) paste0(names(x), collapse = '_')) %>% unlist %>% table

# length(hi_var_tad)
# length(const_hi_tad)
# length(const_lo_tad)

# (length(hi_var_tad)+length(const_hi_tad))/nrow(tads)
# (length(const_hi_tad)+length(const_lo_tad))/nrow(tads)
# (length(hi_var_tad)+length(const_lo_tad))/nrow(tads)

# ci2_tads/sum(ci2_tads)

# ci3_tbl = names(chosen_inds[ci_len == 3]) %>% stringr::str_split('\\.') %>% sapply(function(x) x[[1]]) %>% table
# ci3_tbl = setNames(as.numeric(ci3_tbl), names(ci3_tbl))
# ci3_tbl

# take_ci2 = 1000 - ci3_tbl
# take_ci2

# head(names(chosen_inds[ci_len == 2]))

# ci_st = stringr::str_split(names(chosen_inds), '\\.') %>% sapply(function(x) x[[1]])
# head(ci_st)

# ci_pat = stringr::str_split(names(chosen_inds), '\\.') %>% sapply(function(x) x[[2]])
# head(ci_pat)

# chosen_ci2 = names(chosen_inds)[c(head(which(ci_len == 2 & ci_st == 'NSC' & ci_pat == 'seqs_desc'), take_ci2[['NSC']]),
#                                   head(which(ci_len == 2 & ci_st == 'IPC'), take_ci2[['IPC']]),
#                                   head(which(ci_len == 2 & ci_st == 'Mature'), take_ci2[['Mature']]))]

# chosen_ci3 = names(chosen_inds)[ci_len == 3]

# chosen_ci2 %>% stringr::str_split('\\.') %>% sapply(function(x) paste0(x[1:2], collapse = '_')) %>% table

# chosen_all = c(chosen_ci3, chosen_ci2)

# head(chosen_all)

# chosen_all %>% stringr::str_split('\\.') %>% sapply(function(x) paste0(x[1:2], collapse = '_')) %>% table

# length(chosen_all)

# control_seqs = lapply(chosen_inds[chosen_all], function(x) lapply(x, function(y) y[[1]]))
# # control_seqs = control_seqs[!is.na(control_seqs)]

names(chosen_ind_list) = seqs_motifs$seq_name[as.numeric(names(chosen_ind_list))]

tad3_inds = which(sapply(chosen_ind_list, length) == 3)
lib_nat = seqs_motifs[match(names(chosen_ind_list)[tad3_inds], seqs_motifs$seq_name),c('chrom', 'start', 'end', 'seq_name')]
lib_nat = dplyr::mutate(lib_nat, is_native = TRUE)
lib_cont = as.data.frame(cbind(dels_motifs[unlist(chosen_ind_list[tad3_inds]), c('chrom', 'start', 'end')], 
                               names(unlist(chosen_ind_list[tad3_inds]))))
colnames(lib_cont) = c('chrom', 'start', 'end', 'seq_name')
lib_cont = dplyr::mutate(lib_cont, is_native = FALSE)

head(lib_nat)
head(lib_cont)

final_library = as.data.frame(rbind(lib_nat, lib_cont))
rownames(final_library) = 1:nrow(final_library)

final_library$type = stringr::str_split(final_library$seq_name, '\\.') %>% purrr::map(1) %>% unlist
final_library$pattern = stringr::str_split(final_library$seq_name, '\\.') %>% purrr::map(2) %>% unlist
final_library$source = NA
final_library$source[!final_library$is_native] = gsub('.const_hi|.const_lo|.hi_var', '', final_library$seq_name[!final_library$is_native])
final_library$coords = gsub(' ', '', apply(final_library[,1:3], 1, paste0, collapse = '-'))

head(final_library)

tail(final_library)

fl_coords_mod = gsub(' ', '', apply(dplyr::mutate(final_library[,1:3], start = start - 117, end = end+117), 1, paste0, collapse = '-'))

head(fl_coords_mod)

sum(sapply(enh_cor_md, function(x) length(which(x$enh %in% fl_coords_mod))))

head(unlist(seqs_to_take))

sum(unlist(lapply(seqs_to_take, lapply, length)))

e14_data = readr::read_tsv('./data/E14_MPRA.tsv')

e14_data$seq_name = gsub(' ', '', apply(e14_data[,1:3], 1, paste0, collapse = '-'))

e14_data = dplyr::relocate(.data = e14_data, seq_name, .after = end)

colnames(dels_motifs)

e14_nei_dels = gintervals.neighbors(dels_motifs[,1:4], e14_data[,1:4], maxdist = 2e+3)

head(e14_nei_dels)

e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')] = apply(subset(e14_data, select = c(NSC_pvalue, IPC_pvalue, PN_pvalue)), 2, p.adjust, method = 'BH')

e14_enh_sig = e14_data$seq_name[apply(e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')], 1, function(x) any(x < 0.01))]


length(which(apply(e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')], 1, function(x) any(x < 0.01))))

apply(e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')], 2, function(x) length(which(x < 0.01)))

apply(e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')], 2, hist)

e14_motifs = gextract(trk_s, intervals = e14_data[,1:3], iterator = 10)

e14_motifs = e14_motifs[with(e14_motifs, order(chrom, start, end)),]

e14_motifs_avg = tgs_matrix_tapply(t(subset(e14_motifs, select = -c(chrom, start, end, intervalID))), e14_motifs$intervalID, function(x) log(sum(exp(x))))

e14_motifs = cbind(e14_data[,1:4], e14_motifs_avg)

saveRDS(e14_motifs, './data/E14_motifs.rds')

e14_motifs = readRDS('./data/E14_motifs.rds')

e14_motifs = add_tbl_nei(motifs_df = e14_motifs)

head(e14_motifs)

e14_q = sapply(1:ncol(e14_motifs[,-c(1:4)]), function(x,i) motif_ecdfs[[i]](x[,i]), x = e14_motifs[,-c(1:4)])

head(e14_q)

e14_m_z = apply(e14_q, 2, function(x) -log2(1-x))

e14_knn = tgs_cor_knn(t(e14_m_z), t(dels_m_z), knn = 100)

e14_knn$max_diff = mapply(e14_knn$col1, e14_knn$col2, FUN = function(x, y) {
    vec1 = e14_m_z[x,]; vec2 = dels_m_z[y,]; 
#     ord = order(vec1, decreasing = T);
    max(abs(vec1 - vec2))
    }
                          )

e14_knn = e14_knn[e14_knn$max_diff <= 4,]

e14_knn = subset(e14_knn, subset = !(e14_knn$col2 %in% match(e14_nei_dels$seq_name, dels_motifs$seq_name)))

e14_knn$seq_name = e14_motifs$seq_name[e14_knn$col1]

head(dels_motifs)

e14_knn$twin_name = gsub(' ', '', apply(dels_motifs[e14_knn$col2,1:3], 1, paste0, collapse = '-'))

head(e14_knn)

e14_knn = e14_knn[e14_knn$seq_name %in% e14_enh_sig,]

e14_seqs_inds = tapply(e14_knn$max_diff, e14_knn$col1, which.min)

head(e14_seqs_inds)

e14_seqs_choose = as.data.frame(cbind(e14_motifs$seq_name[as.numeric(names(e14_seqs_inds))], 
                                      mapply(as.numeric(names(e14_seqs_inds)), e14_seqs_inds, FUN = function(x,y) return(e14_knn$twin_name[which(e14_knn$col1 == x)[[y]]])),
                                      mapply(as.numeric(names(e14_seqs_inds)), e14_seqs_inds, FUN = function(x,y) return(e14_knn$cor[which(e14_knn$col1 == x)[[y]]])),  
                                      mapply(as.numeric(names(e14_seqs_inds)), e14_seqs_inds, FUN = function(x,y) return(e14_knn$max_diff[which(e14_knn$col1 == x)[[y]]]))  
                                             ))

colnames(e14_seqs_choose) = c('E14_enhancer', 'ENCODE_twin', 'cor', 'max_diff')

e14_seqs_choose_sort = 

head(dplyr::arrange(e14_seqs_choose, max_diff))

head(e14_motifs)

head(dels_motifs)

length(e14_seqs_choose)





dels_motifs[46154,]

cor(abs(e14_knn$cor), 1/(abs(e14_knn$max_diff) + 0.1))

head(e14_knn[order(e14_knn$max_diff),], 20)

inds_core = sample(which(final_library$is_native == T), 2500)

core_source = final_library$seq_name[inds_core]

lib_core = final_library[final_library$seq_name %in% core_source | final_library$source %in% core_source,]

lib_rest = final_library[!(final_library$seq_name %in% core_source | final_library$source %in% core_source),]

nrow(lib_core)

nrow(lib_rest)/4

sum(as.numeric(table(lib_core[,c('type', 'pattern')])))

table(lib_core[,c('type', 'pattern')])/4

head(lib_core)

head(lib_rest)

tail(lib_rest)

readr::write_tsv(lib_core, './data/mpra_core_library_take_2.tsv', )
readr::write_tsv(lib_rest, './data/mpra_spare_library_take_2.tsv')

lctest = readr::read_tsv('./data/mpra_core_library_take_2.tsv')

head(lctest)

tail(lctest)





# nrow(final_library)

# final_library$seq_name %>% stringr::str_split('\\.') %>% sapply(function(x) paste0(x[1:3], collapse = '.')) %>% table %>% sort(decreasing = T) %>% head(50)

# final_library$seq_name %>% stringr::str_split('\\.') %>% sapply(function(x) paste0(x[1:3], collapse = '.')) %>% table %>% table

# length(which(duplicated(final_library[,1:3])))



# ## RNA expression in vicinity of putative enhancers

# cn = c('chrom', 'start', 'end', 'seq_name')

# pmc_st_rm = data.frame(sapply(list(nsc_mc, ipc_mc, mat_mc), function(mci) rowMeans(pmc[,mci])))
# rownames(pmc_st_rm) = gsub(' ', '', apply(pmc_rn, 1, paste0, collapse='-'))
# colnames(pmc_st_rm) = c('NSC', 'IPC', 'Mature')

# get_nei_peak_list = function(nei_df, motifs_df) {
#     v = tapply(nei_df[,8], nei_df[,4], function(x) c(x))
#     v = v[motifs_df$seq_name]
# #     v = ifelse(as.character(sort(unique(motifs_df$intervalID))) %in% names(v), v[as.character(sort(unique(motifs_df$intervalID)))], NA)
# #                v[as.character(which(!(sort(unique(motifs_df$intervalID)) %in% as.numeric(names(v)))))] = NA
# #     v = v[order(as.numeric(names(v)))]
#     return(v)
# }

# add_atac_umis = function(motifs_df, pmc_orig, pmc_st_rm) {
#     seq_st = unlist(purrr::map(stringr::str_split(motifs_df$seq_name, '\\.'), 1))
#     seqs_nei_5k = gintervals.neighbors(motifs_df[,cn], pmc_orig, mindist = 0, maxdist = 5e+03, maxneighbors = 50)
#     seqs_nei_25k = gintervals.neighbors(motifs_df[,cn], pmc_orig, mindist = 5e+03, maxdist = 25e+03, maxneighbors = 50)
#     seqs_nei_250k = gintervals.neighbors(motifs_df[,cn], pmc_orig, mindist = 25e+03, maxdist = 2.5e+5, maxneighbors = 50)
#     seqs_nei_1250k = gintervals.neighbors(motifs_df[,cn], pmc_orig, mindist = 2.5e+5, maxdist = 1.25e+6, maxneighbors = 50)
#     motifs_df$peaks_5k = get_nei_peak_list(seqs_nei_5k, motifs_df)
#     motifs_df$peaks_25k = get_nei_peak_list(seqs_nei_25k, motifs_df)
#     motifs_df$peaks_250k = get_nei_peak_list(seqs_nei_250k, motifs_df)
#     motifs_df$peaks_1250k = get_nei_peak_list(seqs_nei_1250k, motifs_df)
#     motifs_df$atac_5k = unlist(purrr::map2(motifs_df$peaks_5k, seq_st, function(.x, .y) sum(pmc_st_rm[.x,.y])))
#     motifs_df$atac_25k = unlist(purrr::map2(motifs_df$peaks_25k, seq_st, function(.x, .y) sum(pmc_st_rm[.x,.y])))
#     motifs_df$atac_250k = unlist(purrr::map2(motifs_df$peaks_250k, seq_st, function(.x, .y) sum(pmc_st_rm[.x,.y])))
#     motifs_df$atac_1250k = unlist(purrr::map2(motifs_df$peaks_1250k, seq_st, function(.x, .y) sum(pmc_st_rm[.x,.y])))
#     motifs_df[,grep('atac', colnames(motifs_df))] = apply(motifs_df[,grep('atac', colnames(motifs_df))], 2, function(x) ifelse(is.na(x), 0, x))
#     return(motifs_df)
# }


# get_nei_gene_list = function(nei_df, motifs_df, genes_egc) {
#     nei_gene_list = as.list(tapply(nei_df[,'geneSymbol'],nei_df[,'seq_name'], function(x) {y = unique(c(x)); return(y[y %in% genes_egc])}))
#     nei_gene_list = nei_gene_list[motifs_df$seq_name]
#     names(nei_gene_list) = 1:length(nei_gene_list)
# #     nei_gene_list[nei_gene_list == NULL] = c()
#     return(nei_gene_list)
# }

# add_genes_umis = function(motifs_df, genes_egc, egc_st_rm) {
# #     genes_mat = rownames(mat@mat)
# #     mat_rm = Matrix::rowMeans(mat@mat)
#     cn = c('chrom', 'start', 'end', 'seq_name')
#     seq_st = unlist(purrr::map(stringr::str_split(motifs_df$seq_name, '\\.'), 1))
#     seqs_nei_5k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 5e+3, maxneighbors = 50)
#     seqs_nei_25k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 2.5e+4, maxneighbors = 50)
#     seqs_nei_250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+4 + 1, maxdist = 2.5e+5, maxneighbors = 50)
#     seqs_nei_1250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+5 + 1, maxdist = 1.5e+6, maxneighbors = 50)
#     motifs_df$genes_5k = get_nei_gene_list(seqs_nei_5k, motifs_df, genes_egc)
#     motifs_df$genes_25k = get_nei_gene_list(seqs_nei_25k, motifs_df, genes_egc)
#     motifs_df$genes_250k = get_nei_gene_list(seqs_nei_250k, motifs_df, genes_egc)
#     motifs_df$genes_1250k = get_nei_gene_list(seqs_nei_1250k, motifs_df, genes_egc)
#     motifs_df$umis_5k = unlist(purrr::map2(motifs_df$genes_5k, seq_st, function(.x, .y) ifelse(length(.x) > 0, sum(egc_st_rm[.x,.y]), 0)))
#     motifs_df$umis_25k = unlist(purrr::map2(motifs_df$genes_25k, seq_st, function(.x, .y) ifelse(length(.x) > 0, sum(egc_st_rm[.x,.y]), 0)))
#     motifs_df$umis_250k = unlist(purrr::map2(motifs_df$genes_250k, seq_st, function(.x, .y) ifelse(length(.x) > 0, sum(egc_st_rm[.x,.y]), 0)))
#     motifs_df$umis_1250k = unlist(purrr::map2(motifs_df$genes_1250k, seq_st, function(.x, .y) ifelse(length(.x) > 0, sum(egc_st_rm[.x,.y]), 0)))
#     return(motifs_df)
# }

# rm(list = ls()[grep('seqs_nei|motifs_df|seq_st|nei_df', ls())])

# genes_egc = rownames(mc_rna@e_gc)

# egc_st_rm = sapply(list(nsc_mc, ipc_mc, mat_mc), function(mci) {
#     rowMeans(mc_rna@e_gc[, mci])
# })

# egc_st_rm = egc_st_rm*1e+05

# colnames(egc_st_rm) = c('NSC','IPC','Mature')

# seqs_motifs_add = add_genes_umis(seqs_motifs, genes_egc, egc_st_rm)

# seqs_motifs_add = add_atac_umis(seqs_motifs_add, as.data.frame(cbind(pmc_rn, 1:nrow(pmc))), pmc_st_rm)

# dels_motifs_f = dels_motifs[del_ids,]
# dels_motifs_f$seq_name = seqs_motifs$seq_name
# # dels_motifs_f = dels_motifs_f[with(dels_motifs_f, order(chrom, start, end)),]
# # rownames(dels_motifs_f) = 1:nrow(dels_motifs_f)

# dels_motifs_add = add_genes_umis(dels_motifs_f, genes_egc, egc_st_rm)

# dels_motifs_add = add_atac_umis(dels_motifs_add, as.data.frame(cbind(pmc_rn, 1:nrow(pmc))), pmc_st_rm)

# options(repr.plot.width = 8, repr.plot.height = 8)

# pltmt = tgs_cor(as.matrix(seqs_motifs_add[,grep('nei|umis|atac|intervalID', colnames(seqs_motifs_add))]), 
#        as.matrix(dels_motifs_add[,grep('nei|umis|atac|dist', colnames(dels_motifs_add))]), spearman=T)
# pltmt
# pheatmap(pltmt, cluster_cols = F, cluster_rows = F, fontsize = 14)

# add_df_all = cbind(seqs_motifs_add, dels_motifs_add)

# add_df_ints = add_df_all[,grep('chrom|start|end', colnames(add_df_all))]

# add_df_ints[,4] = as.character(add_df_ints[,4])

# colnames(add_df_ints) = c('chrom', 'start', 'end', 'chrom', 'start', 'end')

# head(add_df_ints)

# nrow(add_df_ints)

# nrow(gintervals.neighbors(add_df_ints[,1:3], add_df_ints[,4:6], maxdist=100))

# motif_add_nei = gintervals.neighbors(dels_motifs_add[,1:3], pmc_rn[,1:3], maxdist = 1e+3)

# nrow(motif_add_nei)

# eps= 1e-04

# # dist_inds = which(dels_motifs_add$dist >= MAX_DIST)
# dist_inds = which(dels_motifs_add$dist >= 1.5e+6)

# length(dist_inds)

# trk_cols = intersect(grep('jaspar|cis|jolma|motifs', colnames(seqs_motifs_add), v=T),
#                           grep('jaspar|cis|jolma|motifs', colnames(dels_motifs_add), v=T))
# # trk_cols

# ind = 17
# plot(as.numeric(dels_motifs_add[ind,trk_cols]), as.numeric(seqs_motifs_add[ind,trk_cols]))

# lines(c(min(as.numeric(dels_motifs_add[ind,trk_cols])), max(as.numeric(dels_motifs_add[ind,trk_cols]))),
#      c(min(as.numeric(dels_motifs_add[ind,trk_cols])), max(as.numeric(dels_motifs_add[ind,trk_cols]))),
#      lty = 'dashed', lwd = 1.5)
# text(min(as.numeric(dels_motifs_add[ind,trk_cols]))+5, max(as.numeric(seqs_motifs_add[ind,trk_cols])), 
#      labels = cor(as.numeric(dels_motifs_add[ind,trk_cols]), as.numeric(seqs_motifs_add[ind,trk_cols])))

# res = sapply(1e+5*seq(1,25,1), function(d) {
#     dist_inds = which(dels_motifs_add$dist >= d)
#     print(length(dist_inds))
#     num_lst = list(
#         '5k' = c(mean(log2(eps+seqs_motifs_add$umis_5k[dist_inds])), mean(log2(eps+dels_motifs_add$umis_5k[dist_inds]))),
#         '25k' = c(mean(log2(eps+seqs_motifs_add$umis_25k[dist_inds])), mean(log2(eps+dels_motifs_add$umis_25k[dist_inds]))),
#         '250k' = c(mean(log2(eps+seqs_motifs_add$umis_250k[dist_inds])), mean(log2(eps+dels_motifs_add$umis_250k[dist_inds]))),
#         '1250k' = c(mean(log2(eps+seqs_motifs_add$umis_1250k[dist_inds])), mean(log2(eps+dels_motifs_add$umis_1250k[dist_inds])))
#         )
#     num_df = do.call('cbind', num_lst)
#     apply(num_df, 2, diff)
#                 }
#        )

# # apply(res, 1, order)

# res

# plot(0, xlim = c(0.5e+5, 27e+5), ylim = c(-4,2))
# apply(res, 1, function(xi) {
#     lines(seq(1,25)*1e+5, xi, col = sample(grep('white|grey|gray', colors(), inv=T, v=T), 1))
# })

# umi_df = add_df_all[,c(1:3,grep('umis', colnames(add_df_all)))]
# umi_df[,-c(1:3)] = umi_df[,-c(1:3)] + 1e-02
# mean_lfc = apply(umi_df[,-c(1:3)], 1, function(x) mean(log2(x[1:4]/x[5:8])))

# atac_df = add_df_all[,c(1:3,grep('atac', colnames(add_df_all)))]
# atac_df[,-c(1:3)] = atac_df[,-c(1:3)] + 1e-02
# mean_lfc_atac = apply(atac_df[,-c(1:3)], 1, function(x) mean(log2(x[1:4]/x[5:8])))

# eps = 1e-02

# all_umi_df = add_df_all[,c(1:3,grep('umis|atac', colnames(add_df_all)))]
# all_umi_df[,-c(1:3)] = all_umi_df[,-c(1:3)] + 1e-02
# mean_lfc_all = apply(all_umi_df[,-c(1:3)], 1, function(x) mean(c(log2(x[1:4]/x[9:12]),
#                                                                      log2(x[5:8]/x[13:16]))))

# length(which(mean_lfc >= log2(10) | mean_lfc_atac >= log2(10)))

# length(which(mean_lfc >= log2(10)))
             
# length(which(mean_lfc_atac >= log2(10)))

# lfc_inds_to_select = which(mean_lfc >= log2(10) | mean_lfc_atac >= log2(10))

# length(lfc_inds_to_select)

# length(unique(lfc_inds_to_select))

# hist(as.numeric(rownames(add_df_all[lfc_inds_to_select,grep('chrom|start|end', colnames(add_df_all))])))

# grep('chrom|start|end|seq_name$', colnames(add_df_all))

# lfc_df = add_df_all[lfc_inds_to_select,grep('chrom|start|end|seq_name$', colnames(add_df_all))]

# check_unique_frac(lfc_df)

# head(lfc_df)

# head(rbind(as.matrix(lfc_df[,1:3]), as.matrix(lfc_df[,5:7])))

# check_unique_frac(rbind(as.matrix(lfc_df[,1:3]), as.matrix(lfc_df[,5:7])))



# lfc_df = as.data.frame(rbind(cbind(as.matrix(lfc_df[,1:4]), rep('native', nrow(lfc_df))), cbind(as.matrix(lfc_df[,5:8]), rep('control', nrow(lfc_df)))))

# colnames(lfc_df) = c('chrom', 'start', 'end', 'seq_name', 'source')

# lfc_df[,c('start', 'end')] = apply(lfc_df[,c('start', 'end')], 2, as.numeric)

# head(lfc_df)

# # gsetroot('/home/aviezerl/mm10')

# # chain_path_mm10 = '/home/feshap/raid/data_other//mm9ToMm10.over.chain'
# # # gintervals.load_chain(chain_path_mm10)

# # lfc_df_mm10 = gintervals.liftover(lfc_df, chain = chain_path_mm10)

# # lfc_df$source = lfc_df$source[lfc_df_mm10$intervalID]

# # lfc_df$seq_name = lfc_df$seq_name[lfc_df_mm10$intervalID]

# # bad_len_names = lfc_df_mm10$seq_name[lfc_df_mm10$end - lfc_df_mm10$start != 266]
# # bad_len_names

# # no_twin_names = names(table(lfc_df_mm10$seq_name)[table(lfc_df_mm10$seq_name) < 2])
# # no_twin_names

# # lfc_df_mm10_filt = lfc_df_mm10[!(lfc_df_mm10$seq_name %in% c(no_twin_names, bad_len_names)),]

# num_seqs_vec = setNames(unlist(num_seqs_to_take), gsub('\\.', '_seqs_', names(unlist(num_seqs_to_take))))

# num_seqs_vec

# num_taken_by_lfc = as.matrix(table(sapply(stringr::str_split(lfc_df$seq_name, '\\.'), function(x) paste0(x[1:2], collapse='_'))))/2

# num_left_tot = NUM_FINAL - nrow(lfc_df)
# num_left_tot

# num_desired = round(num_seqs_vec*num_left_tot/sum(num_seqs_vec))
# num_desired

# rest_of_seqs = sapply(seq_along(num_desired), function(x, n, i) {
#     spl = unlist(stringr::str_split(n[[i]], '_'))
#     st = spl[[1]]
#     pt = paste0(spl[2:3], collapse='_')
#     v = seqs_to_take[[st]][[pt]]
#     nv = paste0(st, '.', pt, '.', names(v))
#     v_ret = head(v[!(nv %in% lfc_df$seq_name)], x[[i]])
#     nv_ret = head(nv[!(nv %in% lfc_df$seq_name)], x[[i]])
#     return(setNames(v_ret, nv_ret))
# }, x = num_desired, n = names(num_desired))

# rest_df = do.call('rbind', lapply(rest_of_seqs, function(x) {
#     y = as.data.frame(do.call('rbind', stringr::str_split(x, '-')))
#     rownames(y) = names(x)
#     return(y)
# }
#                            )      
#                  )

# colnames(rest_df) = c('chrom', 'start', 'end')
# rest_df[,c('start', 'end')] = apply(rest_df[,c('start', 'end')], 2, as.numeric)

# rest_mids = round(0.5*(rest_df$start + rest_df$end))
# rest_df$start = rest_mids - EXT
# rest_df$end = rest_mids + EXT

# rest_df$seq_name = rownames(rest_df)
# rest_df$source = 'rest'

# final_library = rbind(lfc_df, rest_df)
# final_library = cbind(final_library, do.call('rbind', stringr::str_split(final_library$seq_name, '\\.')))
# final_library = final_library[with(final_library, order(chrom, start, end)),]
# colnames(final_library) = c('chrom', 'start', 'end', 'seq_name', 'source', 'st', 'enh_type', 'enh_num')
# rownames(final_library) = 1:nrow(final_library)
# nrow(final_library)

# length(unique(spare_inds))

# nrow(unique(final_library[,c('chrom', 'start', 'end')]))

# nrow(final_library[,c('chrom', 'start', 'end')])

# check_unique_frac = function(df) {
#     n1 = nrow(df[,c('chrom', 'start', 'end')])
#     n2 = nrow(unique(df[,c('chrom', 'start', 'end')]))
#     print(c(n1, n2, n2/n1))
# }

# spare_inds = sample(which(final_library$source == 'rest'), 1800)
# spare_library = final_library[spare_inds,]
# core_library = final_library[!(1:nrow(final_library) %in% spare_inds),]

# nrow(seqs_motifs_add)
# nrow(unique(seqs_motifs_add[,c('chrom', 'start', 'end')]))

# head(add_df_all)

# head(core_library)

# nrow(core_library)

# nrow(unique(core_library[,c('chrom', 'start', 'end')]))

# cols_to_take = grep('chrom|start|end|seq_name|jolma|jaspar|cis|motifs', colnames(seqs_motifs_add), v=T)
# motifs_df_all = as.data.frame(rbind(seqs_motifs_add[,cols_to_take],
#                      dels_motifs_add[,cols_to_take]))
# library_nat_cnt_inds = which(core_library$source %in% c('native', 'control'))
# motifs_df_library = dplyr::left_join(core_library[library_nat_cnt_inds,])

# dim(motifs_df_all)

# cbind(colnames(seqs_motifs_add), colnames(dels_motifs_add))

# all_motifs_

# plot_final_library = function(sm, dm) {
    
# }







# head(final_library)

# readr::write_tsv(core_library, './data/enh_mpra_core_library.tsv')

# readr::write_tsv(spare_library, './data/enh_mpra_spare_library.tsv')



# pmc = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')
# prox_enh = unlist(read.delim('./data/mmcortex_proximal_enhancers_mm10.txt', h = F))
# dist_enh = unlist(read.delim('./data/mmcortex_distal_enhancers_mm10.txt', h = F))

# pmc = pmc[rownames(pmc) %in% sort(c(prox_enh, dist_enh)),]

# nrow(final_library[final_library$source == 'control',])

# nrow(final_library[final_library$source == 'native',])

# nei_control = gintervals.neighbors(lfc_df[lfc_df$source == 'control',], pmc_rn[,1:3], maxdist = 1e+3)

# gsetroot('/home/aviezerl/mm10')

# head(pmc_rn)

# nei_control = gintervals.neighbors(final_library[final_library$source == 'control',], pmc_rn, maxdist = 1e+3)

# nrow(nei_control)

# head(prox_enh)

# seqs_plot = seqs_unlist[unique(final_library$seq_name[with(final_library, order(source, st, enh_type))])]
# # seqs_plot = seqs_unlist[]

# ann_colors$source = setNames(gplots::col2hex(c('red', 'yellow', 'blue')), unique(final_library$source))

# ann_colors$enh_type = setNames(gplots::col2hex(c('purple', 'green', 'orange')), unique(final_library$enh_type))

# row_annot = final_library[match(unique(final_library$seq_name[with(final_library, order(source, st, enh_type))]), final_library$seq_name),c('source', 'st', 'enh_type')]

# row_annot$st[row_annot$st == 'Mature'] = 'CPN_L2-3'









# rownames(row_annot) = seqs_plot

# head(row_annot)

# pltmt = pmc_norm[seqs_plot,c(nsc_mc[order(mcmd$mean_day[nsc_mc])], ipc_mc[order(mcmd$mean_day[ipc_mc])], mat_mc[order(mcmd$mean_day[mat_mc])])]

# colnames(pltmt) = c(nsc_mc[order(mcmd$mean_day[nsc_mc])], ipc_mc[order(mcmd$mean_day[ipc_mc])], mat_mc[order(mcmd$mean_day[mat_mc])])

# options(repr.plot.width = 18, repr.plot.height = 18)

# eps = 1e-04

# dim(apply(pltmt, 1, function(x) log2((x+eps)/(median(x) + eps))))

# pltmt = t(apply(pltmt, 1, function(x) log2((x+eps)/(median(x) + eps))))
# l = 100

# brks = c(seq(min(pltmt), 0, length.out = l/2),
#         seq(0+(max(pltmt) - min(pltmt))/l, max(pltmt), l=l/2+1))
# clrmp = colorRampPalette(c('blue4', 'white', 'red4'))(l)

# pp = pheatmap(pltmt,
#          cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_colors = ann_colors,
#          annotation_col = col_annot, 
#          annotation_row = row_annot,
#         color = clrmp,
#          breaks = brks
# #          color = colorRampPalette(c('white', 'yellow', 'red', 'brown', 'black'))(1000)
#         )

# save_pheatmap_png(pp, './figs/enh_mpra_library_hm.png', h = 3000, w = 2500)

# misha.ext::fwrite_ucsc(intervals = final_library[,1:3], file = './data/mmcortex_enh_library.bed', name = 'mmcortex_enh_MPRA_library', type = 'bed', append = F, span = 10)

# mean(log2(eps+seqs_motifs_add$umis_5k[dist_inds])) -mean(log2(eps+dels_motifs_add$umis_5k[dist_inds]))

# mean(log2(eps+seqs_motifs_add$umis_25k[dist_inds])) - mean(log2(eps+dels_motifs_add$umis_25k[dist_inds]))

# mean(log2(eps+seqs_motifs_add$umis_250k[dist_inds])) - mean(log2(eps+dels_motifs_add$umis_250k[dist_inds]))

# mean(log2(eps+seqs_motifs_add$umis_1250k[dist_inds])) - mean(log2(eps+dels_motifs_add$umis_1250k[dist_inds]))

# options(repr.plot.width = 4, repr.plot.height = 4)

# plot(density(log2(eps+seqs_motifs_add$umis_5k[dist_inds])))
# lines(density(log2(eps+dels_motifs_add$umis_5k[dist_inds])), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_25k[dist_inds])))
# lines(density(log2(eps+dels_motifs_add$umis_25k[dist_inds])), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_250k[dist_inds])))
# lines(density(log2(eps+dels_motifs_add$umis_250k[dist_inds])), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_1250k[dist_inds])))
# lines(density(log2(eps+dels_motifs_add$umis_1250k[dist_inds])), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_5k)))
# lines(density(log2(eps+dels_motifs_add$umis_5k)), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_25k)))
# lines(density(log2(eps+dels_motifs_add$umis_25k)), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_250k)))
# lines(density(log2(eps+dels_motifs_add$umis_250k)), col = 'red')

# plot(density(log2(eps+seqs_motifs_add$umis_1250k)))
# lines(density(log2(eps+dels_motifs_add$umis_1250k)), col = 'red')




