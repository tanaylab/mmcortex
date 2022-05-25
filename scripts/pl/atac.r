library(metacell)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
library(misha)

gsetroot('/home/aviezerl/mm10')
scdb_init(file.path(wd, 'scdb'), force_reinit = T)
scfigs_init(file.path(wd, 'figs'))

atac_files = list.files(file.path(wd, 'scatac_data'))
# atac_files

agg_id = vroom::vroom(file.path(wd, 'scatac_data', atac_files[[2]]))

mc = scdb_mc("pl")
lfp = log2(mc@mc_fp)
mc2d = scdb_mc2d("pl")

promo_peaks = readRDS(file.path(wd, 'scatac_data/prom_counts.RDS'))
gb_peaks = readRDS(file.path(wd, 'scatac_data/genebody_counts.RDS'))
ig_peaks = readRDS(file.path(wd, 'scatac_data/peak_counts.RDS'))

make_mat = function(mat, nm) {
    rep_nums = sapply(stringr::str_split(colnames(mat), '-'),function(x) x[[2]])
    reps = agg_id$sample_name[match(rep_nums, agg_id$library_id)]
    cmd = tibble::column_to_rownames(data.frame(cbind(colnames(mat), rep_nums, reps)), 'V1')
    cmd$day = as.numeric(gsub('E', '', sapply(stringr::str_split(cmd$reps, '_'), function(x) x[[1]])))
    mat_raw = scm_new_matrix(mat, cmd)
    scdb_add_mat(nm, mat_raw)
}

make_mat(promo_peaks, 'pl_prom_raw')
make_mat(gb_peaks, 'pl_gb_raw')
make_mat(ig_peaks, 'pl_ig_raw')

csi = Matrix::colSums(ig_peaks)
csp = Matrix::colSums(promo_peaks)
csg = Matrix::colSums(gb_peaks)
q_max = 0.98
q_min = 0.1

cells_ig = colnames(ig_peaks)[csi >= quantile(csi, q_min) & csi <= quantile(csi, q_max)]
cells_gb = colnames(gb_peaks)[csg >= quantile(csg, q_min) & csg <= quantile(csg, q_max)]
cells_prom = colnames(promo_peaks)[csp >= quantile(csp, q_min) & csp <= quantile(csp, q_max)]

good_cells = intersect(cells_gb, intersect(cells_ig, cells_prom))


# gb_md = strsplit(names(rownames(gb_peaks)), '-') %>% data.frame %>% t %>% data.frame
# gb_md = gb_md[!duplicated(rownames(gb_peaks)),]
# colnames(gb_md) = c('chrom', 'start', 'end')
# gb_md[,c('start', 'end')] = apply(gb_md[,c('start', 'end')], 2, as.numeric)
# gb_md = mutate(gb_md, len = abs(end-start))
# rownames(gb_md) = rownames(gb_peaks)[!duplicated(rownames(gb_peaks))]
# gb_norm = gb_peaks[!duplicated(rownames(gb_peaks)),]/gb_md$len


mcell_add_gene_stat(mat_id = 'prom_raw', gstat_id = 'prom_raw', force = T)

gstat = scdb_gstat('test_prom')
plot(log(gstat$ds_mean), gstat$ds_log_varmean)
x = log(gstat$ds_mean)
init_filt = which(x >= -5)
y = gstat$ds_log_varmean[init_filt]
name_filt = gstat$name
x = x[init_filt]
xcut = cut(x, breaks = seq(min(x), max(x), l = 30))
top_q_inds = lapply(levels(xcut), function(l) {inds = which(xcut == l);
                                  xfilt = y[inds];
                                  xtop = head(inds[order(xfilt, decreasing = T)], 40);
                                  return(xtop)
                                 }
      )
names(top_q_inds) = levels(xcut)

var_genes_atac = rownames(gstat)[init_filt[unlist(top_q_inds)]]

feats = scdb_gset('pl_f')

gb = union(var_genes_atac, feats)



nms = intersect(rownames(lfp), rownames(promo_peaks))
high_var = intersect(rownames(lfp)[apply(lfp, 1, function(x) any(x > 0.8))], rownames(promo_peaks))
feats = scdb_gset('all_feats_f')
feats_int = intersect(names(feats@gene_set), rownames(promo_peaks))

aumi = scm_downsamp(promo_peaks[,good_cells], quantile(csi, q_min))

aumi_n = aumi - Matrix::rowMeans(aumi)
# asd = apply(as.matrix(aumi), 1, sd)
# aumi_n = (aumi - Matrix::rowMeans(aumi))/asd

bumi = scm_downsamp(umis = gb_peaks[,good_cells], n = quantile(csg, q_min))

bumi_n = bumi - Matrix::rowMeans(bumi)

# bsd = apply(as.matrix(bumi), 1, sd)
# bumi_n = (bumi - Matrix::rowMeans(bumi))/bsd
aumi_mat = as.matrix(aumi_n)
bumi_mat = as.matrix(bumi_n)

cell_names = intersect(promo_peaks@Dimnames[[2]], good_cells)
atac_rep = lapply(cell_names, function(x) sapply(stringr::str_split(x, pattern = '-'), "[[", 2)) %>% unlist %>% as.numeric
atac_metadata = data.frame(cbind(cell_names, atac_rep))
atac_metadata[,'batch_set_id'] = agg_id[atac_rep,'sample_name']
atac_metadata = tibble::column_to_rownames(atac_metadata, 'cell_names')

atac_metadata$day = gsub('E', '', gsub('_rep\\d', '', atac_metadata$batch_set_id)) %>% as.numeric

mat = scdb_mat('all_rec_bon_1')
mc_day = mat@cell_metadata[names(mc@mc), 'day']
names(mc_day) = mc@mc

day_mc = split(mc_day, names(mc_day))

day_mc_stat = purrr::map_dfc(day_mc, function(x) c(mean(x), median(x), sd(x))) %>% t %>% data.frame

day_mc_stat = day_mc_stat[order(as.numeric(rownames(day_mc_stat))),]
colnames(day_mc_stat) = c('mean', 'median','sd')

mc_md = vroom::vroom('./BonevCollab/mc_metadata_new.tsv')

mcs_with_day = lapply(sort(unique(mc_day)), function(u) sort(unique(mc@mc[mc_day == u])))

names(mcs_with_day) = sort(unique(mc_day))

dms_cut = cut(day_mc_stat$mean, breaks = seq(11,18,0.5))

lev_13_18 = levels(dms_cut)[4:length(levels(dms_cut))]

mc_lev = lapply(lev_13_18, function(l) which(dms_cut == l))
names(mc_lev) = lev_13_18
# mc_lev

mc_lev_merge = list('13' = unlist(mc_lev[1:2]), 
                    '14' = unlist(mc_lev[3:4]), 
                    '15' = unlist(mc_lev[5:6]), 
                    '16' = unlist(mc_lev[7:8]), 
                    '17' = unlist(mc_lev[9:10]),
                    '18' = unlist(mc_lev[11]))

mat_lev = split(rownames(atac_metadata)[rownames(atac_metadata) %in% colnames(aumi_n)], 
                atac_metadata[rownames(atac_metadata) %in% colnames(aumi_n),'day'])

# cross_strat = lapply(names(mc_lev_merge), function(l) tgs_cor(mc@mc_fp[nms,mc_lev_merge[[l]]], as.matrix(aumi_n[nms,mat_lev[[l]]])))

# cross_strat_gb = lapply(names(mc_lev_merge), function(l) tgs_cor(mc@mc_fp[nms,mc_lev_merge[[l]]], as.matrix(bumi_n[nms,mat_lev[[l]]])))

# cross_strat = lapply(names(mc_lev_merge), function(l) tgs_cor(mc@mc_fp[nms,mc_lev_merge[[l]]], exp(aumi_mat[nms,mat_lev[[l]]])))

# cross_strat_gb = lapply(names(mc_lev_merge), function(l) tgs_cor(mc@mc_fp[nms,mc_lev_merge[[l]]], exp(bumi_mat[nms,mat_lev[[l]]])))

# cross_strat = lapply(seq_along(mc_lev_merge), function(l,n,m,a,i) {
#                     if (i >= 2 & i < length(l)) {inds = (i-1):(i+1)} 
#                     else if (i == 1) {inds = i:(i+1)}
#                     else {inds = (i-1):i}
#                     inds_mc = unlist(l[inds])
# #                     inds_sc = unlist(m[inds])
#                     tgs_cor(mc@mc_fp[nms,inds_mc], exp(a[nms,m[[n[[i]]]]]), pairwise.complete.obs = T)
#                     }, l = mc_lev_merge, n = names(mc_lev_merge), m = mat_lev, a = aumi_mat)

# head(exp(bumi_mat[nms,head(mat_lev[[names(mc_lev_merge)[[seq_along(mc_lev_merge)[[1]]
#                                                         ]]
#                                    ]]
#                            )
#                  ]
#         )
#     )

# cross_strat_gb = lapply(seq_along(mc_lev_merge), function(l,n,m,b,i) {
#                     if (i >= 2 & i < length(l)) {inds = (i-1):(i+1)} 
#                     else if (i == 1) {inds = i:(i+1)}
#                     else {inds = (i-1):i}
#                     print(i)
#                     print(n[[i]])
#                     print(head(m[[n[[i]]]]))
#                     print(length(which(!(m[[n[[i]]]] %in% colnames(b)))))
#                     inds_mc = unlist(l[inds])
# #                     inds_sc = unlist(m[inds])
#                     tgs_cor(mc@mc_fp[nms,inds_mc], exp(b[nms,m[[n[[i]]]]]))
#                     }, l = mc_lev_merge, n = names(mc_lev_merge), m = mat_lev, b = bumi_mat)

cross_strat = lapply(seq_along(mcs_with_day), function(l,n,m,a,i) {
                    inds_mc = unlist(l[[i]])
                    tgs_cor(mc@mc_fp[nms,inds_mc], exp(a[nms,m[[n[[i]]]]]))
                    }, l = mcs_with_day, n = names(mcs_with_day), m = mat_lev, a = aumi_mat)

cross_strat_gb = lapply(seq_along(mcs_with_day), function(l,n,m,b,i) {
                    inds_mc = unlist(l[[i]])
                    tgs_cor(mc@mc_fp[nms,inds_mc], exp(b[nms,m[[n[[i]]]]]))
                    }, l = mcs_with_day, n = names(mcs_with_day), m = mat_lev, b = bumi_mat)

hit_a_strat = lapply(cross_strat, function(x) setNames(as.numeric(apply(x, 2, max)), colnames(x))) %>% unlist

hit_b_strat = lapply(cross_strat_gb, function(x) setNames(as.numeric(apply(x, 2, max)), colnames(x))) %>% unlist

hit_a_strat_i = lapply(cross_strat, function(x) setNames(as.numeric(rownames(x)[apply(x, 2, which.max)]), colnames(x))) %>% unlist

hit_b_strat_i = lapply(cross_strat_gb, function(x) setNames(as.numeric(rownames(x)[apply(x, 2, function(y) which.max(y)[1])]), colnames(x))) %>% unlist

hit_ab_i = ifelse(hit_a_strat > hit_b_strat, hit_a_strat_i, hit_b_strat_i)
hit_ab_i[is.na(hit_b_strat)] = hit_a_strat_i[is.na(hit_b_strat)]
hit_ab_i[is.na(hit_a_strat)] = hit_b_strat_i[is.na(hit_a_strat)]

hit_ab = ifelse(hit_a_strat > hit_b_strat, hit_a_strat, hit_b_strat)
hit_ab[is.na(hit_b_strat)] = hit_a_strat[is.na(hit_b_strat)]
hit_ab[is.na(hit_a_strat)] = hit_b_strat[is.na(hit_a_strat)]

hist(hit_ab)

st_a = mc_md$Bonev_annotation[hit_a_strat_i]
st_b = mc_md$Bonev_annotation[hit_b_strat_i]

mismatch_cells = names(hit_a_strat_i)[st_a != st_b]
table(ifelse(hit_a_strat[mismatch_cells] > hit_b_strat[mismatch_cells], T, F))

mismatch_cells = names(hit_a_strat_i)[st_a != st_b]
table(ifelse(hit_a_strat > hit_b_strat, T, F))

hist(hit_a_strat)

hist(hit_b_strat)

dim(promo_peaks)

prob = seq(0.05, 1, 0.05)
qh_df = data.frame(cbind(
    quantile(hit_a_strat, probs = prob, na.rm = T),
    quantile(hit_b_strat, probs = prob, na.rm = T),
    quantile(hit_ab, probs = prob, na.rm = T)
                        )
                  )

qh_df

# atac_max_mc = lapply(cross_strat, function(x) rownames(x)[apply(x, 2, which.max)]) %>% unlist

# atac_max_cor = lapply(cross_strat, function(x) apply(x, 2, max)) %>% unlist

# cross_a = tgs_cor(log2(mc@e_gc[nms,] + 1e-6), as.matrix(aumi_n)[nms,])
# hit_a = apply(cross_a,2,max)
# hit_a_i = apply(cross_a,2,which.max)

# cross_a = tgs_cor(mc@mc_fp[nms,], as.matrix(aumi_n)[nms,])
# hit_a = apply(cross_a,2,max)
# hit_a_i = apply(cross_a,2,which.max)

promo_in_aumi = promo_peaks[,colnames(promo_peaks) %in% colnames(aumi_n)]
gb_in_aumi = gb_peaks[,colnames(gb_peaks) %in% colnames(bumi_n)]

prom_mat = scm_new_matrix(mat = promo_in_aumi, 
                          cell_metadata = atac_metadata[rownames(atac_metadata) %in% colnames(aumi_n),])
gb_mat = scm_new_matrix(mat = gb_in_aumi, 
                          cell_metadata = atac_metadata[rownames(atac_metadata) %in% colnames(bumi_n),])

prom_mat@cell_metadata[names(hit_ab_i),'st'] = mc_md$Bonev_annotation[hit_ab_i]
prom_mat@cell_metadata[names(hit_ab_i),'day'] = gsub('E', '', gsub('_rep\\d', '', prom_mat@cell_metadata[names(hit_ab_i),'day']))
prom_mat@cell_metadata[names(hit_ab_i),'mc'] = hit_ab_i

prom_mat@cell_metadata[names(hit_ab_i),'mc2'] = hit_ab_i
prom_mat@cell_metadata[names(hit_ab_i),'st2'] = mc_md$Bonev_annotation[hit_ab_i]

gb_mat@cell_metadata[names(hit_ab_i),'st'] = mc_md$Bonev_annotation[hit_ab_i]
gb_mat@cell_metadata[names(hit_ab_i),'day'] = gsub('E', '', gsub('_rep\\d', '', gb_mat@cell_metadata[names(hit_ab_i),'batch_set_id']))
gb_mat@cell_metadata[names(hit_ab_i),'mc'] = hit_ab_i

prom_mat@cell_metadata[names(hit_ab_i),'cortical'] = ifelse(mc_md$Comment[hit_ab_i] == 'cortical', T, F)
gb_mat@cell_metadata[names(hit_ab_i),'cortical'] = ifelse(mc_md$Comment[hit_ab_i] == 'cortical', T, F)


prom_mat@cell_metadata$hit_ab = hit_ab[rownames(prom_mat@cell_metadata)]
# prom_mat@cell_metadata$hit_ab = hit_ab_i[rownames(prom_mat@cell_metadata)]

dim(prom_mat@mat)

options(repr.plot.width = 15)
boxplot(hit_ab ~ st, prom_mat@cell_metadata, cex.axis = 0.5)

cort_atac_sc = rownames(prom_mat@cell_metadata)[prom_mat@cell_metadata$cortical == T & prom_mat@cell_metadata$day != 12]
prom_cort = scm_new_matrix(mat = prom_mat@mat[,cort_atac_sc],
                          cell_metadata = prom_mat@cell_metadata[cort_atac_sc,])
gb_cort = scm_new_matrix(mat = gb_mat@mat[,cort_atac_sc], 
                          cell_metadata = gb_mat@cell_metadata[cort_atac_sc,])

dim(prom_cort@mat)

scdb_add_mat('prom_cort', prom_cort)

scdb_add_mat('gb_cort', gb_cort)

tbl_mc = table(hit_ab_i)
tbl_mc = tbl_mc[order(tbl_mc, decreasing = T)]
hist(tbl_mc, 100)

tbl_st = table(prom_mat@cell_metadata$st[prom_mat@cell_metadata$day != 12])
tbl_st = tbl_st[order(tbl_st, decreasing = T)]/sum(tbl_st)

tbl_st

tbl_st_rna = table(mat@cell_metadata$st_bon[mat@cell_metadata$day != 12])
tbl_st_rna = tbl_st_rna[names(tbl_st)]/sum(tbl_st_rna)

tbl_st_rna

plot(as.numeric(tbl_st), as.numeric(tbl_st_rna), 
     col = mcmd$color[match(names(tbl_st), mcmd$st)], 
     pch = 16, cex = 2)
lines(seq(0, 0.2, 0.1), seq(0, 0.2, 0.1))

cor(as.numeric(tbl_st), as.numeric(tbl_st_rna), method = 'spearman')

dim(filter(prom_mat@cell_metadata, cortical == F & day != 12))

cort_atac_sc = rownames(prom_mat@cell_metadata)[prom_mat@cell_metadata$cortical == T & prom_mat@cell_metadata$day != 12]
length(cort_atac_sc)

saveRDS(cort_atac_sc, file = './data/cort_atac_sc.rds')

## Map sc ATAC to all MCs
lapply(tail(seq_along(mat_lev), -1), function(n, l, mc_md, i) {png(paste0('./figs/atac_umap_day_all_mc', n[[i]], '.png')); 
                                        sc_inds = l[[i]];
                                        mcs_in = hit_a_strat_i[names(hit_a_strat_i) %in% sc_inds]
                                        tblt = table(mcs_in)
                                        plot(mc2d@sc_x, mc2d@sc_y, col='gray', cex = 0.3, cex.main = 1,
                                             main = paste0('day ', n[[i]], '\n n_sc = ', length(sc_inds),
                                                          ', n_mc = ', length(tblt)
#                                                            , ', mapped_mc = ', length(which(tblt > 0)), 
#                                                            ', non_mapped_mc = ', length(which(tblt == 0))
                                                          ));
                                        points(mc2d@mc_x[unique(mcs_in)], mc2d@mc_y[unique(mcs_in)], pch=21, 
                                             bg=mc@colors[unique(mcs_in)],
                                              cex=pmin(tblt/20, 2));
#                                         sc_inds = names(mc2d@sc_x)[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% l[[i]]])]
#                                         points(mc2d@sc_x[sc_inds], mc2d@sc_y[sc_inds], col='gray', cex = 0.3)
                                       dev.off()}, n = names(mat_lev), l = mat_lev, mc_md = mc_md)#, 
#                                       cex=pmin(tabulate(hit_a_strat_i[hit_a_strat_i %in% l],nbins=max(mc@mc)),2)))

## Map sc ATAC to MCs by rolling day
# seq_along(mc_lev_merge), function(l,n,m,a,i) {
#                     if (i >= 2 & i < length(l)) {inds = (i-1):(i+1)} 
#                     else if (i == 1) {inds = i:(i+1)}
#                     else {inds = (i=1):i}
#                     inds_mc = unlist(l[inds])
lapply(seq_along(mc_lev_merge), function(l, m, n, i) {png(paste0('./figs/atac_umap_roll_day_', n[[i]], '.png')); 
                                        if (i >= 2 & i < length(l)) {inds = (i-1):(i+1)} 
                                        else if (i == 1) {inds = i:(i+1)}
                                        else {inds = (i-1):i}
                                        l_inds = sort(unlist(l[inds]))
#                                         l_inds = sort(l[[i]]);
                                        sc_inds = m[[n[[i]]]]
                                        tblt = table(hit_ab_i[sc_inds])
                                        plot(mc2d@sc_x, mc2d@sc_y, col='gray', cex = 0.3, cex.main = 1,
                                             main = paste0('day ', n[[i]], '\n n_sc = ', length(sc_inds),
                                                          ', n_mc = ', length(l_inds), ', mapped_mc = ', length(tblt), 
                                                           ', non_mapped_mc = ', length(l_inds[!(l_inds %in% as.numeric(names(tblt)))])
                                                          ));
                                        points(mc2d@mc_x[as.numeric(names(tblt))], mc2d@mc_y[as.numeric(names(tblt))], pch=21, 
                                             bg=mc@colors[as.numeric(names(tblt))],
                                              cex=pmin(tblt/20, 2));
#                                         sc_inds = names(mc2d@sc_x)[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% l[[i]]])]
#                                         points(mc2d@sc_x[sc_inds], mc2d@sc_y[sc_inds], col='gray', cex = 0.3)
                                       dev.off()}, n = names(mc_lev_merge), l = mc_lev_merge, m = mat_lev)#, 
#                                       cex=pmin(tabulate(hit_a_strat_i[hit_a_strat_i %in% l],nbins=max(mc@mc)),2)))

i = 13
l = mcs_with_day
sc_inds = rownames(prom_mat@cell_metadata)[prom_mat@cell_metadata$day == i]
mc_mapped = unique(prom_mat@cell_metadata$mc[prom_mat@cell_metadata$day == i]);
tblt = table(hit_ab_i[hit_ab_i %in% mc_mapped])
l_inds = l[[as.character(i)]]

## Map sc ATAC to MCs by day
lapply(seq(13,18), function(l, n, i) {
    png(paste0('./figs/atac_umap_day_', i, '.png')); 
                                        sc_inds = rownames(prom_mat@cell_metadata)[prom_mat@cell_metadata$day == i]
                                        mc_mapped = unique(prom_mat@cell_metadata$mc[prom_mat@cell_metadata$day == i]);
                                        tblt = table(hit_ab_i[sc_inds])
                                        l_inds = l[[as.character(i)]]
                                        plot(mc2d@sc_x, mc2d@sc_y, col='gray', cex = 0.3, cex.main = 1,
                                             main = paste0('day ', i, '\n n_sc = ', sum(tblt),
                                                          ', n_mc = ', length(l_inds), ', mapped_mc = ', 
                                                           length(mc_mapped), 
                                                           ', non_mapped_mc = ', length(l_inds[!(l_inds %in% mc_mapped)])
                                                          ));
                                        points(mc2d@mc_x[l_inds], mc2d@mc_y[l_inds], pch=21, 
                                             bg=mc@colors[l_inds],
                                              cex=pmin(tblt/20, 2));
#                                         sc_inds = names(mc2d@sc_x)[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% l[[i]]])]
#                                         points(mc2d@sc_x[sc_inds], mc2d@sc_y[sc_inds], col='gray', cex = 0.3)
                                       dev.off()
                                                  }, n = names(mcs_with_day), l = mcs_with_day)#, 
#                                       cex=pmin(tabulate(hit_a_strat_i[hit_a_strat_i %in% l],nbins=max(mc@mc)),2)))

lapply(seq_along(mc_lev_merge), function(l, n, i) {png(paste0('./figs/atac_umap_full_day_', n[[i]], '.png')); 
                                        plot(mc2d@sc_x, mc2d@sc_y, col='gray', cex = 0.3,  main = paste0('day ', n[[i]], ' full'));
                                        points(mc2d@mc_x[l[[i]]], mc2d@mc_y[l[[i]]], pch=21, 
                                             bg=mc@colors[l[[i]]], main = paste0('day ', n[[i]]));
#                                               cex=pmin(tabulate(hit_a_strat_i[hit_a_strat_i %in% l[[i]]],nbins=max(mc@mc))/20, 2));
#                                         sc_inds = names(mc2d@sc_x)[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% l[[i]]])]
#                                         points(mc2d@sc_x[sc_inds], mc2d@sc_y[sc_inds], col='gray', cex = 0.3)
                                       dev.off()}, n = names(mc_lev_merge), l = mc_lev_merge)#, 
#                                       cex=pmin(tabulate(hit_a_strat_i[hit_a_strat_i %in% l],nbins=max(mc@mc)),2)))

tbl_mc_atac = table(prom_mat@cell_metadata$mc[!is.na(prom_mat@cell_metadata$mc)])
tbl_mc_rna = table(mc@mc[mc@mc %in% names(tbl_mc_atac)])

# mc = scdb_mc('all_rec_bon_1')

# mc2d = scdb_mc2d('all_umap')

# mcell_mc2d_plot_by_factor(mc2d_id = 'all_umap', mat_id = 'prom', meta_field = 'day')

lapply(seq_along(mc_lev_merge), function(i) head(hit_a_strat_i[hit_a_strat_i %in% mc_lev_merge[[i]]]))

options(repr.plot.width = 10, repr.plot.height = 6)
xx = tbl_mc_atac/sum(tbl_mc_atac)
yy = tbl_mc_rna/sum(tbl_mc_rna)
zz = seq(0,0.01, 0.0005)
cor(as.vector(xx),as.vector(yy))
plot(as.vector(xx),as.vector(yy))
lines(zz, zz)

tbl_prom = prom_mat@cell_metadata %>% dplyr::group_by(day) %>% dplyr::select(st) %>% table

tbl_prom2 = prom_mat@cell_metadata %>% dplyr::group_by(day) %>% dplyr::select(st2) %>% table

mat = scdb_mat('all_rec_bon_1')

tbl_mat = mat@cell_metadata %>% dplyr::group_by(day) %>% dplyr::select(st_bon) %>% table 

mat_day_st = tidyr::pivot_wider(data.frame(tbl_mat), names_from = 'st_bon', values_from = 'Freq')
prom_day_st = tidyr::pivot_wider(data.frame(tbl_prom), names_from = 'st', values_from = 'Freq')
prom_day_st2 = tidyr::pivot_wider(data.frame(tbl_prom2), names_from = 'st2', values_from = 'Freq')
cols_both = colnames(mat_day_st)[colnames(mat_day_st) %in% colnames(prom_day_st)]

mat_day_st[,cols_both]
prom_day_st[,cols_both]
prom_day_st2[,cols_both]

tbl_at = table(prom_mat@cell_metadata$st)/sum(table(prom_mat@cell_metadata$st))

tbl_r = table(mat@cell_metadata$st_bon)/sum(table(mat@cell_metadata$st_bon))

names(tbl_r)[!(names(tbl_r) %in% st_both)]

st_both = intersect(names(tbl_at), names(tbl_r))

cor(as.numeric(tbl_at[st_both]/sum(tbl_at[st_both])), as.numeric(tbl_r[st_both]/sum(tbl_r[st_both])), method='spearman')

png('./figs/st_atac_strat_vs_rna.png')
options(repr.plot.width = 6, repr.plot.height = 6)
plot(as.numeric(tbl_at[st_both]), as.numeric(tbl_r[st_both]), xlab = 'atac_frac', ylab = 'rna_frac')
zzz = seq(0, 0.2, 0.01)
lines(zzz, zzz)
dev.off()

calc_hit_i = function(gset) {
    print(head(gset))
    cross_feats = tgs_cor(lfp[gset,], as.matrix(aumi_n)[gset,])
    hit = apply(cross_feats,2,max)
    hit_i = apply(cross_feats,2,which.max)
    return(hit_i)
}

hit_i_df = purrr::map_dfr(list(nms, high_var, feats_int), function(x) calc_hit_i(x))

hit_tab = data.frame(apply(as.matrix(hit_i_df), 1, tabulate))
head(hit_tab)

tgs_cor(hit_tab, spearman = T)

thresh = 20
hit_filt = dplyr::filter(hit_tab, ((X1 <= thresh) | (X2 <= thresh) | (X3 <= thresh)))

thresh = seq(0,50)
thresh_df = data.frame(lapply(thresh, function(y) apply(hit_tab, 2, function(x) length(which(x <= y)))))

png(filename = './figs/less_than_x.png', width = 500, height = 500)
plot(0:50, thresh_df[1,], type = 'line', col = 'red', xlab = 'x', ylab = '# of mc having less than x cells')
lines(0:50, thresh_df[2,], col = 'green')
lines(0:50, thresh_df[3,], col = 'blue')
legend('topleft', legend = c('all', 'high_var', 'feature'), 
       col = c('red', 'green', 'blue'), lwd = 2,
      )
dev.off()

plot(hit_tab[,1], hit_tab[,2], pch = 16, col = mc@colors)
plot(hit_tab[,1], hit_tab[,3], pch = 16, col = mc@colors)
plot(hit_tab[,3], hit_tab[,2], pch = 16, col = mc@colors)

hit_diff = apply(hit_tab, 1, diff)

head(hit_tab)

png('./figs/atac_on_rna_proj.png', height = 1200, width = 1200)
plot(mc2d@mc_x, mc2d@mc_y, pch=21, bg=mc@colors, cex=pmin(tabulate(hit_i,nbins=max(mc@mc))/5,2))
dev.off()

promo_filt = promo_peaks[,dimnames(promo_peaks)[[2]] %in% dimnames(aumi)[[2]]]

atac_mat = scdb_mat('prom')

mcell_new_mc(mc_id = 'prom', mc = hit_i, outliers = list(), scmat = atac_mat)

atac_on_mcs = t(tgs_matrix_tapply(promo_filt, hit_i, sum))

head(atac_on_mcs)

q_aom = data.frame(quantile(colSums(atac_on_mcs), (0:20)/20))
head(q_aom)
