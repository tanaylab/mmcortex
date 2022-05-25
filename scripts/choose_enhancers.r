gsetroot('/home/aviezerl/mm10')

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

bulk_st = c(setNames(rep('NSC', length(nsc_mc)), nsc_mc),
            setNames(rep('IPC', length(ipc_mc)), ipc_mc),
            setNames(rep('Mature', length(mat_mc)), mat_mc))

bulk_st[as.character((1:nrow(mcmd))[!(1:nrow(mcmd)) %in% as.numeric(names(bulk_st))])] = 'Other'

bulk_st = bulk_st[order(as.numeric(names(bulk_st)))]


pmc_m = apply(pmc, 1, mean)
pmc_sd = apply(pmc, 1, sd)

inv_q_pmc_m = ecdf(pmc_m)

pmc_norm = t(apply(pmc, 1, function(x) x/sum(x)))

pmc_norm_mc_max = apply(pmc_norm, 1, which.max)

pmc_norm_st = t(tgs_matrix_tapply(pmc_norm, bulk_st, mean))


pmc_norm_st_diff12 = apply(apply(pmc_norm_st, 1, function(x) x[head(order(x, decreasing = T), 2)]), 2, function(y) diff(rev(y))/mean(y))

norm_st_maxs = colnames(pmc_norm_st)[apply(pmc_norm_st[pmc_norm_st_diff12 >= 0.33,], 1, which.max)]
enh_nsc = which(norm_st_maxs == 'NSC')
enh_ipc = which(norm_st_maxs == 'IPC')
enh_mat = which(norm_st_maxs == 'Mature')

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

num_seqs_to_take = apply(diff_enh_per_st, 2, function(x) {x[x > 475] = 475; return(as.list(x))})
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

length(which(duplicated(seqs_ints[,1:3])))

legc = log2(mc_rna@e_gc + 1e-07)

hi_exp_genes = sapply(list(nsc_mc, ipc_mc, mat_mc), function(mci) apply(legc[,mci], 1, median))

colnames(hi_exp_genes) = c('NSC', 'IPC', 'Mature')

hi_genes_st = lapply(seq_along(1:ncol(hi_exp_genes)), function(i) rownames(hi_exp_genes)[which(hi_exp_genes[,i] >= -14)])

intersect_all = sort(unique(unlist(hi_genes_st)))

seqs_ints_mids = round(0.5*(seqs_ints$start + seqs_ints$end))
seqs_ints$start = seqs_ints_mids - 133
seqs_ints$end = seqs_ints_mids + 133

tss = gintervals.load('tss')

tss = tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]

hi_exp_tss = tss[tss$geneSymbol %in% intersect_all,]

trk_s = unlist(read.delim('./data/motif_tracks_selected_k=64_n=64.txt'))

# trk_all_raw = gextract(trk_s, intervals = seqs_ints, iterator = 10)
# trk_all_raw = trk_all_raw[with(trk_all_raw, order(chrom, start, end)),]

# cn = c('chrom', 'start', 'end', 'intervalID')

# trk_all_avg = tgs_matrix_tapply(t(subset(trk_all_raw, select = -c(chrom, start, end, intervalID))), trk_all_raw$intervalID, function(x) log(sum(exp(x))))

# seqs_motifs = cbind(seqs_ints, trk_all_avg)

# seqs_motifs = seqs_motifs[,!duplicated(colnames(seqs_motifs))]

# saveRDS(seqs_motifs, glue::glue('./data/seqs_ints_trk_s_600_per_cat_mm10.rds'))

# nrow(seqs_motifs)

seqs_motifs = readRDS(glue::glue('./data/seqs_ints_trk_s_600_per_cat_mm10.rds'))

dels_motifs = readRDS('./data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds')

dels_motifs$seq_name = 1:nrow(dels_motifs)

dels_motifs_nei_seqs = gintervals.neighbors(dels_motifs[,c('chrom', 'start', 'end', 'seq_name')], pmc_rn, mindist = -1e+03, maxdist = 1e+3, maxneighbors = 10)

dels_motifs = dplyr::anti_join(dels_motifs, dels_motifs_nei_seqs[,1:3])

dels_motifs = as.data.frame(dels_motifs)

rna_nei = gintervals.neighbors(dels_motifs, hi_exp_tss, maxneighbors = 1, mindist = 0, maxdist = MAX_DIST)

dels_motifs = dplyr::left_join(dels_motifs, subset(rna_nei, select = c(chrom, start, end, dist)), by = c('chrom', 'start', 'end'))

if (length(which(is.na(dels_motifs$dist))) > 0) {dels_motifs$dist[is.na(dels_motifs$dist)] = MAX_DIST}

add_tbl_nei = function(motifs_df) {
    cn = c('chrom', 'start', 'end', 'seq_name')
    seqs_nei_5k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 5e+3, maxneighbors = 50)
    seqs_nei_25k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 0, maxdist = 2.5e+4, maxneighbors = 50)
    seqs_nei_250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+4 + 1, maxdist = 2.5e+5, maxneighbors = 50)
    seqs_nei_1250k = gintervals.neighbors(motifs_df[,cn], tss, mindist = 2.5e+5 + 1, maxdist = 1.5e+6, maxneighbors = 50)
    
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

dels_inds = unique(unlist(sapply(colnames(seqs_motifs), function(x) grep(x, colnames(dels_motifs), ign=T, v=F)[[1]])))

dels_motifs = dels_motifs[,c(dels_inds, grep('dist', colnames(dels_motifs)))]

seqs_m_z = apply(subset(seqs_motifs, select = -c(chrom, start, end, seq_name)), 2, function(x) {v = -log2(1-ecdf(x)(x)); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

dels_m_z = apply(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist)), 2, function(x) {v = -log2(1-ecdf(x)(x)); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

knn_new = tgs_cor_knn(t(seqs_m_z), t(dels_m_z), knn = 60, spearman = F)

knn_new$max_diff = apply(knn_new, 1, function(x) max(seqs_m_z[x[['col1']],] - dels_m_z[x[['col2']],]))

saveRDS(knn_new, './data/enh_dELS_log_sum_exp_cor_knn_log_q_{NUM_TARGET}_mm10.rds')

knn_new = readRDS('./data/enh_dELS_log_sum_exp_cor_knn_log_q_{NUM_TARGET}_mm10.rds')

tads = readRDS('./data/npc_tads_egc.rds')

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

tads_order_rna = apply(tad_st[,c('NSC', 'IPC', 'Mature')], 1, order)


tads = dplyr::mutate(tads, len = end - start) %>% dplyr::relocate(len, .before = Compartment)

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
# tad_annot

plot(tad_min, tad_max)
points(tad_min[hi_var_tad], tad_max[hi_var_tad], col = 'red')
points(tad_min[const_lo_tad], tad_max[const_lo_tad], col = 'blue')
points(tad_min[const_hi_tad], tad_max[const_hi_tad], col = 'green')
legend('bottomright', legend = c('Hi-var', 'Const. lo', 'Const. hi'), col = c('red', 'blue', 'green'), pch = c(1,1,1))

tad_st_all_st = t(tgs_matrix_tapply(as.matrix(tads[,as.character(1:nrow(mcmd))]), mcmd$st, mean))

# tgs_cor(t(tad_all), spearman = T)

tads_nei_dels = gintervals.neighbors(dplyr::mutate(dplyr::select(dels_motifs, chrom, start, end, seq_name), rowname = 1:nrow(dels_motifs)), dplyr::select(tads, chrom, start, end, tadID), maxdist=0)

knn_new = knn_new[knn_new$max_diff <= 3,]

knn_new$native_seq = seqs_motifs$seq_name[knn_new$col1]

dels_rn = gsub(' ', '', apply(dels_motifs[,1:3], 1, paste0, collapse = '-'))

knn_new$peak = dels_rn[match(knn_new$col2, tads_nei_dels$rowname)]

knn_new$tadID = tads_nei_dels$tadID[match(knn_new$col2, tads_nei_dels$rowname)]

knn_new$type_silenced = tad_neg_cf[tads_nei_dels$tadID[match(knn_new$col2, tads_nei_dels$rowname)]]

knn_new$tad_order_rna = apply(tads_order_rna, 2, paste0, collapse='-')[knn_new$tadID]

knn_new$tad_annot = tad_annot[as.character(knn_new$tadID)]

chosen_inds = sapply(sort(unique(knn_new$col1)), function(i) {
    inds = which(knn_new$col1 == i & !is.na(knn_new$tad_annot))
    return(tapply(knn_new$col2[inds], knn_new$tad_annot[inds], function(x) as.numeric(x[1:2])))
}
                  )

length(sort(unique(knn_new$col1)))

length(chosen_inds)

names(chosen_inds) = seqs_motifs$seq_name[sort(unique(knn_new$col1))]

length(unlist(chosen_inds))

length(which(is.na(unlist(chosen_inds))))

head(chosen_inds)

ci_len = lapply(chosen_inds, length)

head(ci_len)

table(unlist(ci_len))

ci3_tbl = names(chosen_inds[ci_len == 3]) %>% stringr::str_split('\\.') %>% sapply(function(x) x[[1]]) %>% table
ci3_tbl = setNames(as.numeric(ci3_tbl), names(ci3_tbl))

take_ci2 = 1000 - ci3_tbl

ci_st = stringr::str_split(names(chosen_inds), '\\.') %>% sapply(function(x) x[[1]])
ci_pat = stringr::str_split(names(chosen_inds), '\\.') %>% sapply(function(x) x[[2]])

chosen_ci2 = names(chosen_inds)[c(head(which(ci_len == 2 & ci_st == 'NSC' & ci_pat == 'seqs_desc'), take_ci2[['NSC']]),
                                  head(which(ci_len == 2 & ci_st == 'IPC'), take_ci2[['IPC']]),
                                  head(which(ci_len == 2 & ci_st == 'Mature'), take_ci2[['Mature']]))]
chosen_ci3 = names(chosen_inds)[ci_len == 3]

chosen_all = c(chosen_ci3, chosen_ci2)

control_seqs = lapply(chosen_inds[chosen_all], function(x) lapply(x, function(y) y[[1]]))

lib_nat = seqs_motifs[match(chosen_all, seqs_motifs$seq_name),c('chrom', 'start', 'end', 'seq_name')]
lib_cont = as.data.frame(cbind(dels_motifs[unlist(control_seqs), c('chrom', 'start', 'end')], names(unlist(control_seqs))))
colnames(lib_cont) = c('chrom', 'start', 'end', 'seq_name')

final_library = as.data.frame(rbind(lib_nat, lib_cont))

length(which(duplicated(final_library[,1:3])))
# chosen_ind_mat = matrix(NA, nrow = length(chosen_inds), ncol = 6)

# colnames(chosen_ind_mat) = unique(unlist(lapply(chosen_inds, function(x) names(x))))

# rownames(chosen_ind_mat) = gsub(' ', '', apply(seqs_motifs[,1:3], 1, paste0, collapse = '-'))

# for (i in 1:length(chosen_inds)) {chosen_ind_mat[i,names(chosen_inds[[i]])] = dels_rn[chosen_inds[[i]]]}

# head(seqs_motifs)

# head(chosen_ind_mat)

# head(knn_new)

# tail(knn_new)







# length(chosen_ind_mat[!is.na(chosen_ind_mat)])

# nrow(pmc_rn)

# length(dels_rn)

# max(chosen_ind_mat[!is.na(chosen_ind_mat)])

# head(do.call('rbind', lapply(chosen_inds, function(x) ))

# hist(sapply(chosen_inds, length))

# length(chosen_inds)

# table(apply(tads_order_rna, 2, paste0, collapse='-'))

# table(apply(tads_order_atac, 2, paste0, collapse='-'))

# tad_order_stat_df = as.data.frame(cbind(table(apply(tads_order_rna, 2, paste0, collapse='-')), table(apply(tads_order_atac, 2, paste0, collapse='-'))))
# colnames(tad_order_stat_df) = c('RNA', 'ATAC')

# tad_order_sig_stat_df = as.data.frame(cbind(table(apply(tads_order_rna[,tads_sig_rna], 2, paste0, collapse='-')), table(apply(tads_order_atac[,tads_sig_atac], 2, paste0, collapse='-'))))
# colnames(tad_order_sig_stat_df) = c('RNA', 'ATAC')

# knitr::kable(tad_order_sig_stat_df, 'simple')

# knitr::kable(tad_order_stat_df, 'simple')

# head(order(knn_new$cor, decreasing = T))

# # n = 100
# # ind = (n*30)+1
# ind = 14431
# plot(seqs_m_z[knn_new$col1[ind],],dels_m_z[knn_new$col2[ind],])
# lines(c(min(seqs_m_z[knn_new$col1[ind],]), max(seqs_m_z[knn_new$col1[ind],])),
#        c(min(seqs_m_z[knn_new$col1[ind],]), max(seqs_m_z[knn_new$col1[ind],])),
#        lty = 'dashed')

# table(tapply(knn_new$tad_order_atac, knn_new$col1, function(x) length(unique(x))))

# table(tapply(knn_new$tad_order_rna, knn_new$col1, function(x) length(unique(x))))

# nrow(tidyr::drop_na(knn_new))

# nrow(knn_new)

# hist(tapply(knn_new$cor, knn_new$col1, max))

# hist(tapply(knn_new$cor, knn_new$col1, min))

# head(tapply(knn_new$tad_order_atac, knn_new$col1, function(x) unique(x)))

# head(tapply(knn_new$tad_order_rna, knn_new$col1, function(x) unique(x)))

# table(knn_new$type_silenced)

# knn_filt = knn_new[knn_new$col2 %in% union(diff_tad_nei_dels$seq_name, diff_tad_atac_nei_dels$seq_name),]

# head(seqs_motifs)

# knn_filt$native_seq = seqs_motifs$seq_name[knn_filt$col1]

# length(which(is.na(knn_filt$col2)))

# max(knn_filt$col2)

# max(tads_nei_dels$rowname)

# nrow(dels_motifs)

# length(which(is.na(tads_nei_dels$tadID[match(knn_filt$col2, tads_nei_dels$rowname)])))

# length(which(!is.na(tads_nei_dels$tadID[match(knn_filt$col2, tads_nei_dels$rowname)])))

# hist(tads_nei_dels$tadID[match(knn_filt$col2, tads_nei_dels$rowname)])

# length(tads_nei_dels$tadID[match(knn_filt$col2, tads_nei_dels$rowname)])

# head(tad_atac_neg_cf)

# knn_filt$type_silenced = tad_atac_neg_cf[tads_nei_dels$tadID[match(knn_filt$col2, tads_nei_dels$rowname)]]

# nrow(knn_filt)

# nrow(tidyr::drop_na(knn_filt))

# head(knn_filt)

# table(knn_filt$type_silenced)

# hist(tapply(knn_filt$cor, knn_filt$col1, max))

# length(unique(knn_filt$col1))

# nrow(seqs_motifs)

# tail(rownames(dels_motifs))

# max(dels_motifs$seq_name)





# seqs_dels_bind = do.call('cbind', list(seqs_motifs[,1:3], 
#                                        dels_motifs[knn_z$id[,1],1:3], 
#                                        dels_motifs[knn_z$id[,2],1:3], 
#                                        dels_motifs[knn_z$id[,3],1:3], 
#                                        dels_motifs[knn_z$id[,4],1:3])
# #                                        dels_motifs[knn_z$id[,5],1:3])
#                         )

# eq1 = seqs_dels_bind[,1] == seqs_dels_bind[,4] & abs(seqs_dels_bind[,2] - seqs_dels_bind[,5]) <= MAX_DIST
# eq2 = seqs_dels_bind[,1] == seqs_dels_bind[,7] & abs(seqs_dels_bind[,2] - seqs_dels_bind[,8]) <= MAX_DIST
# eq3 = seqs_dels_bind[,1] == seqs_dels_bind[,10] & abs(seqs_dels_bind[,2] - seqs_dels_bind[,11]) <= MAX_DIST
# eq4 = seqs_dels_bind[,1] == seqs_dels_bind[,13] & abs(seqs_dels_bind[,2] - seqs_dels_bind[,14]) <= MAX_DIST
# # eq5 = seqs_dels_bind[,1] == seqs_dels_bind[,16] & abs(seqs_dels_bind[,2] - seqs_dels_bind[,17]) <= MAX_DIST

# eq_df = cbind(eq1, eq2, eq3, eq4)

# table(apply(eq_df, 1, function(x) which(!x)[[1]]))

# del_ids_inds = apply(eq_df, 1, function(x) which(!x)[[1]])

# del_ids = sapply(seq_along(del_ids_inds), function(x,kid,i) kid[i,x[[i]]], x = del_ids_inds, kid = knn_z$id)

# # del_ids = ifelse(eq1, knn_z$id[,2], knn_z$id[,1])
# match_dels = dels_motifs[del_ids,1:3]

# # del_ids = ifelse(eq1, knn_q$id[,2], knn_q$id[,1])
# # match_dels = dels_motifs[del_ids,1:3]

# match_df = cbind(seqs_motifs[,1:3], match_dels)

# match_df_all = cbind(seqs_motifs, dels_motifs[del_ids,])

# hist(match_df_all$dist)

# tgs_cor(as.matrix(match_df_all[,grep('nei', colnames(match_df_all))[1:4]]), as.matrix(match_df_all[,grep('nei', colnames(match_df_all))[5:8]]), spearman = T)

# # cor_motifs_dels = tgs_cor(t(seqs_m_z), t(dels_m_z[del_ids,]), spearman = T)

# # cor_motifs_dels = tgs_cor(t(seqs_m_z), t(dels_m_z[del_ids,]), spearman = F)

# # cor_motifs_dels = tgs_cor(t(seqs_m_q), t(dels_m_q[del_ids,]), spearman = T)

# options(repr.plot.width = 8)

# # ind = 1
# sapply(sample(1:5000, 10), function(ind) {
#     x = seqs_m_z[ind,]
#     y = dels_m_z[del_ids[[ind]],]
#     plot(x, y, main = ind)
#     lines(c(min(x), max(x)), c(min(x), max(x)), lty = 'dashed', lwd = 1.5)
# })

# length(table(colnames(seqs_m_z)[apply(seqs_m_z, 1, which.max)]))

# length(table(colnames(dels_m_z)[apply(dels_m_z[del_ids,], 1, which.max)]))

# # norm_df = apply(rbind(  table(colnames(seqs_m_z)[apply(seqs_m_z, 1, which.max)]),
# #         table(colnames(dels_m_z)[apply(dels_m_z[del_ids,], 1, which.max)])), 1, function(x) x/sum(x))
# norm_df = t(rbind(  table(colnames(seqs_m_z)[apply(seqs_m_z[,1:(ncol(seqs_m_z)-4)], 1, which.max)]),
#             table(colnames(dels_m_z)[apply(dels_m_z[del_ids,1:(ncol(dels_m_z)-4)], 1, which.max)])))

# plot(log2(norm_df[,1]), log2(norm_df[,2]))
# lines(c(min(log2(norm_df[,1])), max(log2(norm_df[,1]))), 
#      c(min(log2(norm_df[,1])), max(log2(norm_df[,1]))), lty = 'dashed', lwd = 1.5)

# cor(log2(norm_df[,1]), log2(norm_df[,2]), method = 'spearman')

# quantile(diag(cor_motifs_dels) - apply(cor_motifs_dels, 1, mean), (0:20)/20)

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


