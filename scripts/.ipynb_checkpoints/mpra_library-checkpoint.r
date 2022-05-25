ls()

gsetroot('/home/aviezerl/mm10')
library(metacell)
library(dbscan)
library(pheatmap)
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')

MAX_DIST = 2.5e+06
NUM_SHADOW = 1000

KNN_PATH = './data/enh_dELS_lse_cor_knn_fixed.rds'
TADS_PATH = './data/npc_tads_egc.rds'
MCMD_PATH = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv'
PEAK_MC_PATH = './data/pl_cort_peak_mc_smoothed_mg.rds'
PROX_ENH_PATH = './data/mmcortex_proximal_enhancers_mm10.txt'
DIST_ENH_PATH = './data/mmcortex_distal_enhancers_mm10.txt'
TRACKS_PATH = './data/motif_tracks_selected_k=64_n=64.txt'
SEQS_MOTIFS_PATH = './data//seqs_motifs_fixed.rds'
DELS_MOTIFS_PATH = './data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds'
CHOSEN_INDS_PATH = './data/mpra_chosen_ind_list_fixed.rds'
TYPE_AND_TEMPORAL_AND_E14_SHADOW_LIB_PATH = './data/st_and_temporal_and_e14_shadow_enh.tsv'
FINAL_LIBRARY_PATH = './data/mpra_library_take_3.tsv'


### Import stuff

mc_rna = scdb_mc('pl_cort') # metacell object

# metacell metadata
mcmd = vroom::vroom(MCMD_PATH)

# Peak accessibility by metacell matrix (smoothed over mgraph)
pmc = readRDS(PEAK_MC_PATH)

# Convert rownames of peak-mc matrix to misha intervals
pmc_rn = as.data.frame(do.call('rbind', stringr::str_split(rownames(pmc), '-')))
colnames(pmc_rn) = c('chrom', 'start', 'end')
pmc_rn[,c('start', 'end')] = apply(pmc_rn[,c('start', 'end')], 2, as.numeric)
head(pmc_rn)

prox_enh = unlist(read.delim(PROX_ENH_PATH, h = F))
dist_enh = unlist(read.delim(DIST_ENH_PATH, h = F))

# Remove promoter-proximal peaks (this doesn't exclude exonic peaks)
pmc_rn = pmc_rn[rownames(pmc) %in% sort(c(prox_enh, dist_enh)),]
pmc = pmc[rownames(pmc) %in% sort(c(prox_enh, dist_enh)),]
pmc_norm = t(apply(pmc, 1, function(x) x/sum(x)))

## Define some variables for filtering/ordering metacells
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN/CfuPN',
                'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))
cust_mc_ord_st_ord_md = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))[order(mcmd$mean_day[mcmd$st == st])]))

nsc_mc = mcmd$mc[mcmd$st == 'NSC']
ipc_mc = mcmd$mc[grep('IPC', mcmd$st)]
mat_mc = mcmd$mc[mcmd$st %in% c('CPN_L2-3','SCPN','CthPN','CPN_L5_6')]
# ast_mc = mcmd$mc[mcmd$st == 'Astrocytes']
# ol_mc = mcmd$mc[mcmd$st == 'Oligodendrocytes']

# Create vector of "bulk" cell type identities
bulk_st = c(setNames(rep('NSC', length(nsc_mc)), nsc_mc),
            setNames(rep('IPC', length(ipc_mc)), ipc_mc),
            setNames(rep('Mature', length(mat_mc)), mat_mc)
            # setNames(rep('Astrocytes', length(ast_mc)), ast_mc),
            # setNames(rep('Oligodendrocytes', length(ol_mc)), ol_mc)
            )

bulk_st[as.character((1:nrow(mcmd))[!(1:nrow(mcmd)) %in% as.numeric(names(bulk_st))])] = 'Other'
bulk_st = bulk_st[order(as.numeric(names(bulk_st)))]

# Average normalized accessibility over bulk cell types
pmc_norm_st = t(tgs_matrix_tapply(pmc_norm, bulk_st, mean))

# This vector is a measure of specificity of each peak to bulk cell types
# It does so by taking the top 2 bulk cell types per peak and calculating the ratio between their difference and their mean value
pmc_norm_st_diff12 = apply(apply(pmc_norm_st, 1, function(x) x[head(order(x, decreasing = T), 2)]), 2, function(y) diff(rev(y))/mean(y))

# Find the top cell type per peak and assign peaks to cell types subject to specificity threshold
norm_st_maxs = colnames(pmc_norm_st)[apply(pmc_norm_st, 1, which.max)]
enh_nsc = which(norm_st_maxs == 'NSC' & pmc_norm_st_diff12 >= 0.3)
enh_ipc = which(norm_st_maxs == 'IPC' & pmc_norm_st_diff12 >= 0.3)
enh_mat = which(norm_st_maxs == 'Mature' & pmc_norm_st_diff12 >= 0.3)

# lst_enh - list of type-specific peaks
lst_enh = list('NSC' = enh_nsc, 'IPC' = enh_ipc, 'Mature_neurons' = enh_mat)
lst_mc = list('NSC' = nsc_mc, 'IPC' = ipc_mc, 'Mature_neurons' = mat_mc)

eps = 1e-08

## This list of dataframes stored some stats about each type-specific peak
# cor_md_enh - how much accessibility per metacell (in that cell type) is correlated with metacell mean day
# min_max_rat - a measure of variance within that cell type (min-max ratio); those with low ratio are considered "unchanging"/constant
enh_cor_md = lapply(1:3, function(le, lm, i) {
    md = mcmd$mean_day[lm[[i]]]
    cor_md_enh = tgs_cor(t(pmc_norm[le[[i]],lm[[i]]]), as.matrix(md), spearman = T)
    ord_cor = order(cor_md_enh**2, decreasing = T)
    min_max_rat = log2((apply(pmc_norm[le[[i]],lm[[i]]], 1, max) + eps)/(apply(pmc_norm[le[[i]],lm[[i]]], 1, min) + eps))
    ret_mat = data.frame(do.call('cbind', list('enh' = rownames(pmc)[le[[i]][ord_cor]], 
                                               'seq_name' = paste0(names(le)[[i]], '_', 1:length(ord_cor)),
                                               'pmc_norm_st_diff12' = pmc_norm_st_diff12[le[[i]][ord_cor]],
                                               'cor_md' = cor_md_enh[ord_cor], 
                                               'min_max_rat' = min_max_rat[ord_cor])))
    ret_mat[,3:ncol(ret_mat)] = apply(ret_mat[,3:ncol(ret_mat)], 2, as.numeric)
    return(ret_mat)
},  le = lst_enh, lm = lst_mc)


COR_THRESH = 0.5
diff_enh_per_st = do.call('cbind', lapply(enh_cor_md, function(x) list('asc' = length(which(x[,'cor_md'] >= COR_THRESH)),
                                    'desc' = length(which(x[,'cor_md'] <= -COR_THRESH)))))
colnames(diff_enh_per_st) = c('NSC', 'IPC', 'Mature')

diff_enh_per_st

num_seqs_to_take = apply(diff_enh_per_st, 2, function(x) {x[x > 1500] = 1500; return(as.list(x))})
num_seqs_to_take = lapply(num_seqs_to_take, function(x) {x$unc = 4500 - x$asc - x$desc; return(x)})

num_seqs_to_take

seqs_to_take = lapply(seq_along(num_seqs_to_take), function(x,n,i) {
    seqs_asc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'cor_md'], decreasing = T)], x[[i]]$asc), 
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = T),'seq_name'], x[[i]]$asc)) 
    seqs_desc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'cor_md'], decreasing = F)], x[[i]]$desc),
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = F),'seq_name'], x[[i]]$desc)) 
    seqs_unc = setNames(head(rownames(enh_cor_md[[i]])[order(enh_cor_md[[i]][,'min_max_rat'], decreasing = F)], x[[i]]$unc),
                        head(enh_cor_md[[i]][order(enh_cor_md[[i]][,'cor_md'], decreasing = F),'seq_name'], x[[i]]$unc))
    seqs_unc = seqs_unc[!(seqs_unc %in% c(seqs_desc, seqs_asc))]
    return(list('seqs_asc' = seqs_asc, 'seqs_desc' = seqs_desc, 'seqs_unc' = seqs_unc))
}, x = num_seqs_to_take,n = names(num_seqs_to_take))

## Rearrange sequences for downstream analyses

names(seqs_to_take) = names(num_seqs_to_take)
seqs_unlist = do.call('c', lapply(seqs_to_take, unlist))
seqs_ints = data.frame(do.call('rbind', stringr::str_split(seqs_unlist, '-')))
seqs_ints[,2:3] = apply(seqs_ints[,2:3], 2, as.numeric)
colnames(seqs_ints) = c('chrom', 'start', 'end')
seqs_ints$seq_name = names(seqs_unlist)
seqs_ints = seqs_ints[!duplicated(seqs_ints[,1:3]),]
colnames(seqs_ints) = c('chrom', 'start', 'end', 'seq_name')
seqs_ints[,2:3] = apply(seqs_ints[,2:3], 2, as.numeric)
seqs_ints_mids = round(0.5*(seqs_ints$start + seqs_ints$end))
seqs_ints$start = seqs_ints_mids - 133
seqs_ints$end = seqs_ints_mids + 133

## Import motif tracks
trk_s = unlist(read.delim(TRACKS_PATH))

# Extract motif energy for all native sequences in 10bp resolution
trk_all_raw = gextract(trk_s, intervals = seqs_ints, iterator = 10)
trk_all_raw = trk_all_raw[with(trk_all_raw, order(chrom, start, end)),]

#cn = c('chrom', 'start', 'end', 'intervalID')

# Do log-sum-exp transformation for each enhancer (represent each track in each sequence by one number instead of 26-27)
trk_all_avg = tgs_matrix_tapply(t(subset(trk_all_raw, select = -c(chrom, start, end, intervalID))), trk_all_raw$intervalID, function(x) log(sum(exp(x))))

# Attach mean motif energies to intervals
seqs_motifs = cbind(seqs_ints, trk_all_avg)

# Remove duplicate columns/tracks (may be remnant of previous iterations)
seqs_motifs = seqs_motifs[,!duplicated(colnames(seqs_motifs))]
saveRDS(seqs_motifs, glue::glue('./data/seqs_motifs_fixed.rds'))


seqs_motifs = readRDS(SEQS_MOTIFS_PATH)

## Import pre-calculated motif energies for ENCODE enhancers
dels_motifs = readRDS(DELS_MOTIFS_PATH)
dels_motifs$seq_name = 1:nrow(dels_motifs)

## Remove ENCODE enhancers that are too close to native enhancers
dels_motifs_nei_seqs = gintervals.neighbors(dels_motifs[,c('chrom', 'start', 'end', 'seq_name')], pmc_rn, mindist = -1e+03, maxdist = 1e+3, maxneighbors = 10)
dels_motifs = dplyr::anti_join(dels_motifs, dels_motifs_nei_seqs[,1:3])
dels_motifs = as.data.frame(dels_motifs)

### Add distance to highly-expressed genes to distal ENCODE ennhancers
legc = log2(mc_rna@e_gc + 1e-07)
hi_exp_genes = sapply(list(nsc_mc, ipc_mc, mat_mc), function(mci) apply(legc[,mci], 1, median))
colnames(hi_exp_genes) = c('NSC', 'IPC', 'Mature')
hi_genes_st = lapply(seq_along(1:ncol(hi_exp_genes)), function(i) rownames(hi_exp_genes)[which(hi_exp_genes[,i] >= -14)])
intersect_all = sort(unique(unlist(hi_genes_st)))
tss = gintervals.load('tss')
tss = tss[tss$geneSymbol %in% rownames(mc_rna@e_gc),]
hi_exp_tss = tss[tss$geneSymbol %in% intersect_all,]
rna_nei = gintervals.neighbors(dels_motifs, hi_exp_tss, maxneighbors = 1, mindist = -MAX_DIST, maxdist = MAX_DIST)
dels_motifs = dplyr::left_join(dels_motifs, subset(rna_nei, select = c(chrom, start, end, dist)), by = c('chrom', 'start', 'end'))
if (length(which(is.na(dels_motifs$dist))) > 0) {dels_motifs$dist[is.na(dels_motifs$dist)] = MAX_DIST}

### Function to count TSSs in windows surrounding enhancers of interest 
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
    motifs_df$nei_25k = as.numeric(unlist(tbl_nei_25k[motifs_df$seq_name]))
    motifs_df$nei_250k = as.numeric(unlist(tbl_nei_250k[motifs_df$seq_name]))
    motifs_df$nei_1250k = as.numeric(unlist(tbl_nei_1250k[motifs_df$seq_name]))
    motifs_df$nei_5k[is.na(motifs_df$nei_5k)] = 0
    motifs_df$nei_25k[is.na(motifs_df$nei_25k)] = 0
    motifs_df$nei_250k[is.na(motifs_df$nei_250k)] = 0
    motifs_df$nei_1250k[is.na(motifs_df$nei_1250k)] = 0
    
    return(motifs_df)
}

seqs_motifs = add_tbl_nei(seqs_motifs)
dels_motifs = add_tbl_nei(dels_motifs)

## Match column names in native and ENCODE enhancer motif energy matrices
dels_inds = unique(unlist(sapply(colnames(seqs_motifs), function(x) grep(x, colnames(dels_motifs), ign=T, v=F)[[1]])))
dels_motifs = dels_motifs[,c(dels_inds, grep('dist', colnames(dels_motifs)))]

## Generate empirical CDFs for motif energies in ENCODE enhancers
motif_ecdfs = apply(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist)), 2, function(x) ecdf(x))

## Find the quantiles of each motif in each native/ENCODE enhancer in the ENCODE ECDFs
seqs_q = sapply(1:ncol(subset(seqs_motifs, select = -c(chrom, start, end, seq_name))), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = subset(seqs_motifs, select = -c(chrom, start, end, seq_name)))
rownames(seqs_q) = gsub(' ', '', apply(seqs_motifs[,1:3], 1, paste0, collapse = '-'))
dels_q = sapply(1:ncol(subset(dels_motifs, select = -c(chrom, start, end, seq_name, dist))), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = subset(dels_motifs, select = -c(chrom, start, end, seq_name)))
rownames(dels_q) = gsub(' ', '', apply(dels_motifs[,1:3], 1, paste0, collapse = '-'))

## Transform to -log2(1-quantile)
seqs_m_z = apply(seqs_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})
dels_m_z = apply(dels_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})

## Find KNN for each native enhancer in encode enhancers
knn_new = tgs_cor_knn(t(seqs_m_z), t(dels_m_z), knn = 200, spearman = F)

## Find the largest difference in each native-ENCODE enhancer pair from KNN
knn_new$max_diff = mapply(knn_new$col1, knn_new$col2, FUN = function(x, y) {
    vec1 = seqs_m_z[x,]; 
    vec2 = dels_m_z[y,]; 
    return(max(abs(vec1 - vec2)))
})

saveRDS(knn_new, KNN_PATH)

knn_new = readRDS(KNN_PATH)

## Annotate TADs
tads = readRDS(TADS_PATH)
tads = dplyr::mutate(tads, len = end - start) %>% dplyr::relocate(len, .before = Compartment)
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
# table(tad_annot)
tad_st_all_st = t(tgs_matrix_tapply(as.matrix(tads[,as.character(1:nrow(mcmd))]), mcmd$st, mean))
tads_nei_dels = gintervals.neighbors(dplyr::mutate(dplyr::select(dels_motifs, chrom, start, end, seq_name), rowname = 1:nrow(dels_motifs)), dplyr::select(tads, chrom, start, end, tadID), maxdist=0)
rownames(tads_nei_dels) = gsub(' ', '', apply(tads_nei_dels[,1:3], 1, paste0, collapse = '-'))

## Filter knn dataframe and add TAD information for candidate shadow enhancers
knn_new = knn_new[knn_new$max_diff <= 4,]
knn_new$native_seq = seqs_motifs$seq_name[knn_new$col1]
dels_rn = gsub(' ', '', apply(dels_motifs[,1:3], 1, paste0, collapse = '-'))
knn_new$peak = dels_rn[knn_new$col2]
knn_new$tadID = tads_nei_dels$tadID[match(knn_new$col2, rownames(tads_nei_dels))]
knn_new$tad_annot = tad_annot[as.character(knn_new$tadID)]

## Find best shadow enhancer for each native enhancer
chosen_ind_list = as.list(setNames(rep(NA, length(sort(unique(knn_new$col1)))), sort(unique(knn_new$col1))))
for (i in sort(unique(knn_new$col1))) {
    inds = which(knn_new$col1 == i & !is.na(knn_new$tad_annot))
    blklst = which(knn_new$col2[inds] %in% unlist(chosen_ind_list))
    inds = inds[!(1:length(inds) %in% blklst)]
    cands = tapply(knn_new$col2[inds], knn_new$tad_annot[inds], function(x) as.numeric(x[[1]]))
    if (length(cands) >= 2 & all(!(unlist(cands) %in% unlist(chosen_ind_list)))) {
        chosen_ind_list[[as.character(i)]] = cands
    }
}
seqs_motifs$rn = gsub(' ', '', apply(seqs_motifs[,1:3], 1, paste0, collapse = '-'))
names(chosen_ind_list) = seqs_motifs$seq_name[match(sort(unique(knn_new$col1)), seqs_motifs$rn)]


saveRDS(chosen_ind_list, CHOSEN_INDS_PATH)
chosen_ind_list = readRDS(CHOSEN_INDS_PATH)

## Some QC printouts
table(sapply(chosen_ind_list, length))
sum(table(sapply(chosen_ind_list, length)))
length(which(sapply(chosen_ind_list, length) > 1))

seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) paste0(x[1:2], collapse = '_')) %>% table %>% unlist
seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) x[[1]]) %>% table
seqs_motifs$seq_name[which(sapply(chosen_ind_list, length) == 3)] %>% stringr::str_split('\\.') %>% 
            sapply(function(x) x[[1]]) %>% table %>% sum

## Select final library
tad3_inds = which(sapply(chosen_ind_list, length) == 3)
lib_nat = seqs_motifs[match(names(chosen_ind_list)[tad3_inds], seqs_motifs$seq_name),c('chrom', 'start', 'end', 'seq_name')]
lib_nat = dplyr::mutate(lib_nat, is_native = TRUE)
lib_cont = as.data.frame(cbind(dels_motifs[unlist(chosen_ind_list[tad3_inds]), c('chrom', 'start', 'end')], 
                               names(unlist(chosen_ind_list[tad3_inds]))))
colnames(lib_cont) = c('chrom', 'start', 'end', 'seq_name')
lib_cont = dplyr::mutate(lib_cont, is_native = FALSE)

final_library = as.data.frame(rbind(lib_nat, lib_cont))
rownames(final_library) = 1:nrow(final_library)

final_library$type = stringr::str_split(final_library$seq_name, '\\.') %>% purrr::map(1) %>% unlist
final_library$pattern = stringr::str_split(final_library$seq_name, '\\.') %>% purrr::map(2) %>% unlist
final_library$source = NA
final_library$source[!final_library$is_native] = gsub('.const_hi|.const_lo|.hi_var', '', final_library$seq_name[!final_library$is_native])
final_library$coords = gsub(' ', '', apply(final_library[,1:3], 1, paste0, collapse = '-'))

head(final_library)
tail(final_library)

## Check that all enhnacer coordinates in final library are in original type-specific enhancer lists
fl_coords_mod = gsub(' ', '', apply(dplyr::mutate(final_library[,1:3], start = start - 117, end = end + 117), 1, paste0, collapse = '-'))
sum(sapply(enh_cor_md, function(x) length(which(x$enh %in% fl_coords_mod))))


## Import and analyze E14 MPRA data
e14_data = readr::read_tsv('./data/E14_MPRA.tsv')
e14_data$seq_name = gsub(' ', '', apply(e14_data[,1:3], 1, paste0, collapse = '-'))
e14_data = dplyr::relocate(.data = e14_data, seq_name, .after = end)
e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')] = apply(subset(e14_data, select = c(NSC_pvalue, IPC_pvalue, PN_pvalue)), 2, p.adjust, method = 'BH')
e14_enh_sig = e14_data$seq_name[apply(e14_data[,c('NSC_qv', 'IPC_qv', 'PN_qv')], 1, function(x) any(x < 0.01))]


## Find neighbors in ENCODE dELS in order to remove them later
e14_nei_dels = gintervals.neighbors(dels_motifs[,1:4], e14_data[,1:4], maxdist = 2e+3)

## Extract motif energies and transform to -log(1-q) where q is quantile in ENCODE ECDF
e14_motifs = gextract(trk_s, intervals = e14_data[,1:3], iterator = 10)
e14_motifs = e14_motifs[with(e14_motifs, order(chrom, start, end)),]
e14_motifs_avg = tgs_matrix_tapply(t(subset(e14_motifs, select = -c(chrom, start, end, intervalID))), e14_motifs$intervalID, 
                                   function(x) log(sum(exp(x))))
e14_motifs = cbind(e14_data[,1:4], e14_motifs_avg)
saveRDS(e14_motifs, './data/E14_motifs.rds')
e14_motifs = readRDS('./data/E14_motifs.rds')
e14_motifs = add_tbl_nei(motifs_df = e14_motifs)
e14_q = sapply(1:ncol(e14_motifs[,-c(1:4)]), function(x,i) motif_ecdfs[[i]](x[,i]), x = e14_motifs[,-c(1:4)])
rownames(e14_q) = e14_motifs$seq_name
e14_m_z = apply(e14_q, 2, function(x) -log2(1-x))

## Find KNN and annotate
e14_knn = tgs_cor_knn(t(e14_m_z), t(dels_m_z), knn = 100)
e14_knn$max_diff = mapply(e14_knn$col1, e14_knn$col2, FUN = function(x, y) {
    vec1 = e14_m_z[x,]; vec2 = dels_m_z[y,]; 
    return(max(abs(vec1 - vec2)))
})
e14_knn = e14_knn[e14_knn$max_diff <= 4,]
e14_knn = subset(e14_knn, subset = !(e14_knn$col2 %in% match(e14_nei_dels$seq_name, dels_motifs$seq_name)))
e14_knn$seq_name = e14_motifs$seq_name[e14_knn$col1]
e14_knn$twin_name = gsub(' ', '', apply(dels_motifs[e14_knn$col2,1:3], 1, paste0, collapse = '-'))
e14_knn = e14_knn[e14_knn$seq_name %in% e14_enh_sig,]

e14_seqs_inds = tapply(e14_knn$max_diff, as.character(e14_knn$col1), which.min)
e14_seqs_choose = as.data.frame(list(names(e14_seqs_inds), 
                                      mapply(names(e14_seqs_inds), e14_seqs_inds, FUN = function(x,y) return(e14_knn$twin_name[which(e14_knn$col1 == x)[[y]]])),
                                      mapply(names(e14_seqs_inds), e14_seqs_inds, FUN = function(x,y) return(e14_knn$cor[which(e14_knn$col1 == x)[[y]]])),  
                                      mapply(names(e14_seqs_inds), e14_seqs_inds, FUN = function(x,y) return(e14_knn$max_diff[which(e14_knn$col1 == x)[[y]]]))  
                                             ))
colnames(e14_seqs_choose) = c('E14_enhancer', 'ENCODE_twin', 'cor', 'max_diff')
e14_seqs_choose_sort = e14_seqs_choose[order(e14_seqs_choose$cor, decreasing = T),]
e14_seqs_choose_sort = e14_seqs_choose_sort[!duplicated(e14_seqs_choose_sort$ENCODE_twin),]
e14_lib = as.data.frame(do.call('rbind', sapply(e14_seqs_choose_sort$ENCODE_twin, stringr::str_split, '-')))
rownames(e14_lib) = e14_seqs_choose_sort$E14_enhancer
colnames(e14_lib) = c('chrom', 'start', 'end')
e14_lib[,2:3] = apply(e14_lib[,2:3], 2, as.numeric)
stt_lib = as.data.frame(do.call('rbind', sapply(unlist(seqs_to_take), stringr::str_split, '-')))
colnames(stt_lib) = c('chrom', 'start', 'end')
stt_lib[,2:3] = apply(stt_lib[,2:3], 2, as.numeric)
stt_lib$start = stt_lib$start + 117
stt_lib$end = stt_lib$end - 117
stt_e14_final_lib = tibble::rownames_to_column(as.data.frame(rbind(stt_lib, head(e14_lib, NUM_SHADOW))))

readr::write_tsv(stt_e14_final_lib, TYPE_AND_TEMPORAL_AND_E14_SHADOW_LIB_PATH)



inds_core = which(final_library$is_native == T)
core_source = final_library$seq_name[inds_core]
lib_core = final_library[final_library$seq_name %in% core_source | final_library$source %in% core_source,]

readr::write_tsv(lib_core, FINAL_LIBRARY_PATH)