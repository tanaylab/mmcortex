library(metacell)
library(purrr)
library(dplyr)
library(umap)
library(glue)
library(pheatmap)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'all'
SEED = 1337
set.seed(SEED)


db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")




mc = scdb_mc('all_rec_bon_1')
mat = scdb_mat('all_rec_bon_1')
# marks = scdb_gset('all_markers_f')
# feats = scdb_gset('all_feats_f')
# feats = feats@gene_set

mcmd = vroom::vroom(file.path(wd, 'BonevCollab', 'mc_metadata_new.tsv'))

color_key = c(
    
    'nNSC?' = gplots::col2hex('turquoise1'),
    'NSC' = gplots::col2hex('turquoise3') ,
    'early_nNSC?' = gplots::col2hex('turquoise4') ,
    'early_NSC?' = gplots::col2hex('springgreen3') ,
    'IPC' = gplots::col2hex('mediumblue') ,
    
    'CPN_L2-3' = gplots::col2hex('thistle2'),
    'CPN_L5-6' = gplots::col2hex('slategray3'),
    'CthPN' = gplots::col2hex('sienna4') ,
    'iCfuPN' = gplots::col2hex('violetred1') ,
    'iCPN_L2-3' = gplots::col2hex('darkseagreen2'),
    'iCPN_L5-6' = gplots::col2hex('salmon'),
    
    'SCPN' = gplots::col2hex('orange1') ,
    'Stellate_L4' = gplots::col2hex('firebrick3'),

    'CR' = gplots::col2hex('bisque4'),
    'MSN' = gplots::col2hex('olivedrab2'),
    'CGE_IPC' = gplots::col2hex('olivedrab3'),
    'CGE_N' = gplots::col2hex('khaki'),    
    'CGE_NSC' = gplots::col2hex('orange3'),    
    'IN_CGE' = gplots::col2hex('orange4'),
    'IN_MGE' = gplots::col2hex('peachpuff2'),
    'NPC_GE' = gplots::col2hex('olivedrab4'),

    'Astrocytes' = filter(mcmd, inferred_subtype == 'Astrocyte') %>% select(color) %>% unique %>% as.character ,
    'OPC' = filter(mcmd, inferred_subtype == 'OPCs') %>% select(color) %>% unique %>% as.character ,
    
    'Microglia' = filter(mcmd, inferred_subtype == 'Microglia') %>% select(color) %>% unique %>% as.character ,
    'Endothelial' = filter(mcmd, inferred_subtype == 'Vascular Endothelial') %>% select(color) %>% unique %>% as.character ,
    'Pericyte' = filter(mcmd, inferred_subtype == 'Pericyte') %>% select(color) %>% unique %>% as.character 
)

# color_key

mcmd$color3 = color_key[mcmd$Bonev_annotation]
mcmd = relocate(mcmd, color3, .after = color_bon)

# mc_cort = mc
# mc_cort@colors = mcmd$color3
# mc_cort@annots = mcmd$Bonev_annotation


c2r_doub = read.delim(file = './data/doublet_names_cort_4.txt', sep = '\n', header = F) %>% t %>% as.vector
c2r_non_cort = rownames(mat@cell_metadata)[mat@cell_metadata$cortical == FALSE]
c2r_non_cort = c2r_non_cort[!is.na(c2r_non_cort)]
c2r_CR = names(mc@mc)[mc@mc %in% mcmd$mc[mcmd$Bonev_annotation == 'CR']]
c2r_12 = rownames(mat@cell_metadata)[mat@cell_metadata$day == 12]
c2r_hi_umi = colnames(mat@mat[,Matrix::colSums(mat@mat) >= 3.5*1e+04])
c2r_new = unique(c(c2r_non_cort, c2r_CR, c2r_12, c2r_hi_umi, c2r_doub))

# ## c2r - cells to remove/ignore from mc object
# mc_ignore_cells = function(mc, c2r) {
#     frac_in_mc_thresh = 0.001
#     mc_filt = mc
# #     mc_filt@outliers = c(mc_filt@outliers, c2r)
#     tbl_12 = table(mc_filt@mc[names(mc_filt@mc) %in% c2r])
#     tbl_all_filt = table(mc_filt@mc)[names(table(mc_filt@mc)) %in% names(tbl_12)]
#     tbl_frac = tbl_12/tbl_all_filt
#     mcs_to_remove = names(tbl_frac)[tbl_frac >= frac_in_mc_thresh] %>% as.numeric
#     mc_filt@mc = mc_filt@mc[!(names(mc_filt@mc) %in% c2r)]
#     mc_filt@mc = mc_filt@mc[!(mc_filt@mc %in% mcs_to_remove)]
#     mc_filt@mc_fp = mc_filt@mc_fp[,-mcs_to_remove]
#     mc_filt@e_gc = mc_filt@e_gc[,-mcs_to_remove]
#     mc_filt@cov_gc = mc_filt@cov_gc[,-mcs_to_remove]
#     mc_filt@annots = mc_filt@annots[-mcs_to_remove]
#     mc_filt@colors = mc_filt@colors[-mcs_to_remove]
#     mc_filt@n_bc = mc_filt@n_bc[,-mcs_to_remove]
#     return(mc_filt)
# }

# mc_cort_f = mc_ignore_cells(mc_cort, c2r_new)

# c2r_mc = names(mc_cort@mc)[!(names(mc_cort@mc) %in% names(mc_cort_f@mc))]

# mcmd_n = filter(mcmd, mc %in% as.numeric(unique(mc_cort_f@mc)))

# scdb_add_mc('cort', mc_cort_f)
# readr::write_tsv(x = mcmd_n, file = './BonevCollab/mcmd_3.tsv', quote_escape = F)

# c2r_mat = union(mat@ignore_cells, union(c2r_new, c2r_mc))
c2r_mat = union(mat@ignore_cells, c2r_new)

length(c2r_mat)

length(unique(c2r_mat))

mcell_mat_ignore_cells(new_mat_id = 'cort5', mat_id = 'all_rec_bon_1', ig_cells = c2r_mat, reverse = F)

nm = 'cort5'
library(glue)

mcell_add_gene_stat(mat_id = nm, gstat_id = nm, force=T)

mcell_gset_filter_varmean(gset_id=glue("{nm}_feats"), gstat_id=nm, T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id=glue("{nm}_feats"), gstat_id=nm, T_tot=100, T_top3=4)
mcell_plot_gstats(gstat_id=nm, gset_id=glue("{nm}_feats"))

mcell_gset_split_by_dsmat(gset_id=glue("{nm}_feats"), mat_id=nm, K = 96, force = T)

mcell_plot_gset_cor_mats(gset_id=glue("{nm}_feats"), scmat_id=nm)

mat = scdb_mat(nm)
mat = as.matrix(mat@mat)
rn = rownames(mat)
cort_feats = scdb_gset(paste0(nm, '_feats'))

genes = c('Top2a', 'Mki67', 'Ube2c', 'Pcna')
# genes = c('Xist')
get_top_cor_mods = function(g, mat, cort_feats) {
    top_what = 100
    cor_g = tgs_cor(as.matrix(mat[g,]), t(mat), spearman = T)
    rn_ord = rn[order(cor_g, decreasing = T)]
    print(glue('How many genes from each cluster in feature genes are in top {top_what} genes correlated with {g}'))
    print(sort(table(cort_feats@gene_set[names(cort_feats@gene_set) %in% head(rn_ord, top_what)]), decreasing = T))
}
lapply(genes, get_top_cor_mods, mat, cort_feats)


ifn1_genes = c('Isg15', 'Wars', 'Ifit1')
cell_cyc = c('Mki67', 'Pcna', 'Hist1h', 'Smc4', 'Mcm3', 'Top2a')
stress = c('Fos', 'Hsp90ab1', 'Hspa1a', 'Hif1a')
misc = c('Xist', 'Tsix', 'Rps', 'Rpl')
star_genes = c(ifn1_genes, cell_cyc, stress, misc)

# feats = scdb_gset(glue("{nm}_feats"))

genes_in = map(star_genes, function(x) grep(x, names(cort_feats@gene_set), ignore.case = F, v = T)) %>% unlist
uu = unique(cort_feats@gene_set[names(cort_feats@gene_set) %in% genes_in])

star_in = lapply(star_genes, function(x) grep(x, names(cort_feats@gene_set), v=T)) %>% unlist
print(paste0('star_in = ', star_in))
clusts_to_remove = cort_feats@gene_set[star_in] %>% unique
sort(clusts_to_remove)

# xist_cor = tgs_cor(as.matrix(mat['Xist',]), t(mat), spearman = T)

# names(xist_cor) = rownames(mat)

# head(xist_cor[order(xist_cor, decreasing = T)], 20)

# sapply(c(23, 38,43, 62), function(n) sort(names(cort_feats@gene_set)[cort_feats@gene_set %in% n]))

# genes = c('Top2a', 'Mki67', 'Ube2c', 'Pcna')
genes = c('Xist')
get_top_cor_mods = function(g, mat, cort_feats) {
    top_what = 100
    cor_g = tgs_cor(as.matrix(mat[g,]), t(mat), spearman = T)
    rn_ord = rn[order(cor_g, decreasing = T)]
    print(glue('How many genes from each cluster in feature genes are in top {top_what} genes correlated with {g}'))
    print(sort(table(cort_feats@gene_set[names(cort_feats@gene_set) %in% head(rn_ord, top_what)]), decreasing = T))
}
lapply(genes, get_top_cor_mods, mat, cort_feats)


# mat = scdb_mat(nm)
# mat = as.matrix(mat@mat)
# rn = rownames(mat)
# cort_feats = scdb_gset('cort_feats')

# genes = c('Top2a', 'Mki67', 'Ube2c', 'Pcna')
# get_top_cor_mods = function(g, mat, tel_feats) {
#     top_what = 250
#     cor_g = tgs_cor(as.matrix(mat[g,]), t(mat), spearman = T)
#     rn_ord = rn[order(cor_g, decreasing = T)]
#     print(glue('How many genes from each cluster in feature genes are in top {top_what} genes correlated with {g}'))
#     print(sort(table(cort_feats@gene_set[names(cort_feats@gene_set) %in% head(rn_ord, top_what)]), decreasing = T))
# }
# lapply(genes, get_top_cor_mods, mat, tel_feats)


### Not sure about:
## 94 - Fos/Fosb/Egr1/2
## 83 - Has two Hist1h but also Calb2 (didn't remove)
## 16 - Somewhat correlated with Top2a, Mki67, Ube2c, doesn't have CC genes in it but has many marker genes (didn't remove)
# clusts_final = c(7, 20, 21, 24, 26, 37, 63, 72, 79, 94, 95)

### cort_4
# clusts_final = c(7,17,23, 25,38,43, 62, 65, 81, 92, 95)?

### cort5
clusts_final = c(7,22,36,58,62,71,75,81,91,95)

mcell_gset_remove_clusts(gset_id = paste0(nm, '_feats'), filt_clusts = clusts_final, 
                         new_id = paste0(nm, '_feats_lateral'), reverse=T)
mcell_gset_remove_clusts(gset_id = paste0(nm, '_feats'), filt_clusts = clusts_final, 
                         new_id = paste0(nm, '_feats_f'), reverse=F)

mcell_add_cgraph_from_mat_bknn(mat_id=nm,
                gset_id = paste0(nm, '_feats_f'),
                graph_id=nm,
                K=80,
                dsamp=T)

mcell_coclust_from_graph_resamp(
                coc_id=paste0(nm, '_coc500'),
                graph_id=nm,
                min_mc_size=15,
                p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(
                coc_id=paste0(nm, '_coc500'),
                mat_id=nm,
                mc_id=nm,
                K=30, min_mc_size=20, alpha=2)

# mcell_mc_split_filt(new_mc_id='cort_f',
#             mc_id='cort',
#             mat_id='cort',
#             T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(gset_id=paste0(nm, '_markers'), mc_id=nm)

marks = scdb_gset(paste0(nm, '_markers'))
marks_f = marks

lat = scdb_gset(paste0(nm, '_feats_lateral'))

marks_f@gene_set = marks_f@gene_set[!(names(marks_f@gene_set) %in% names(lat@gene_set))]
scdb_add_gset(paste0(nm, '_markers_f'), marks_f)

mcell_mc_plot_marks(mc_id=nm, gset_id=paste0(nm, '_markers_f'), 
                    mat_id=nm, 
#                     lateral_gset=paste0(nm, '_feats_lateral')
                   )

# mcell_mc2d_force_knn(mc2d_id='test_all',mc_id='all', graph_id='all')

# tgconfig::set_param("mcell_mgraph_max_confu_deg",10,"metacell")
# mcell_mgraph_logistic(mgraph_id = 'all_k_10_log', mc_id = 'all_f', 
# 			feats_gset = 'all_feats_f')

mat_old = scdb_mat('all_rec_bon_1')
mat_new = scdb_mat(nm)
mc_old = scdb_mc('all_rec_bon_1')
mc_new = scdb_mc(nm)

feats_old = scdb_gset('all_feats_f')
feats_new = scdb_gset(paste0(nm, '_feats_f'))
genes_both = intersect(names(feats_old@gene_set), names(feats_new@gene_set))
cor_mc = tgs_cor(mc_new@e_gc[genes_both,], mc_old@e_gc[genes_both,], spearman = T)

length(genes_both)

table(mcmd$Bonev_annotation)

mcmd = vroom::vroom('./BonevCollab//mc_metadata_new.tsv')

##Control projection - on mc object with all days and non-cortical
max_st = mcmd[apply(cor_mc, 1, which.max), 'Bonev_annotation'] %>% unlist
table(max_st)[order(table(max_st), decreasing = T)]

# st_tbl_list = lapply(unique(mc_new@mc), function(u) table(mat_new@cell_metadata[names(mc_new@mc)[mc_new@mc == u], 'st_bon']))

# st_tbl_list_n = lapply(st_tbl_list, function(x) x/sum(x))

# hist(lapply(st_tbl_list_n, max) %>% unlist %>% as.numeric, main = 'Fraction of biggest subtype in metacells', xlab = 'Fraction')

mcmd_3 = vroom::vroom('./BonevCollab//mcmd_cort_070221.tsv')
mc_old = scdb_mc('cort_new_2')
mat_old = scdb_mc('cort_new')
feats_old = scdb_gset('cort_feats_f')

cor_mc = tgs_cor(mc_new@e_gc[genes_both,], mc_old@e_gc[genes_both,], spearman = T)
max_st = mcmd_3[apply(cor_mc, 1, which.max), 'st'] %>% unlist



table(max_st)[order(table(max_st), decreasing = T)]

color_key = unique(dplyr::select(mcmd_3, st, color))
color_key

mc_new@colors = color_key$color[match(max_st, color_key$st)]
names(mc_new@colors) = color_key$st[match(max_st, color_key$st)]

head(mc_new@colors)

scdb_add_mc(id = nm, mc_new)

hist(apply(cor_mc, 1, max))

mcell_mc2d_plot_by_factor(mc2d_id = nm, mat_id = nm, meta_field = 'day', single_plot = T)

nm = 'cort5'

library(metacell)
scdb_init('scdb', force_reinit = T)

library(umap)

scfigs_init('figs')

tgconfig::set_param("mcell_mgraph_max_confu_deg",5,"metacell")
mcell_mgraph_logistic(mgraph_id = nm, mc_id = nm, 
			feats_gset = paste0(nm, '_feats_f'))


mc = scdb_mc(nm)
gset = scdb_gset(paste0(nm, '_feats_f'))
feat_genes = names(gset@gene_set)

graph_id = nm
mc2d_id = nm
symmetrize = F
umap_mgraph = F

mgraph = scdb_mgraph(nm)

mgraph = mgraph@mgraph


# mcmd = vroom::vroom('./BonevCollab//mcmd_3_leak.tsv')
# mcmd = dplyr::mutate(mcmd, rn = 1:nrow(mcmd)) %>% dplyr::relocate(rn, .after = 'mc')
# mgraph_new = mgraph
# mgraph_new@mgraph$mc1 = as.numeric(mcmd$mc[mgraph_new@mgraph$mc1])
# mgraph_new@mgraph$mc2 = as.numeric(mcmd$mc[mgraph_new@mgraph$mc2])
# mgraph = mgraph_new@mgraph

uconf = umap.defaults
uconf$n_neighbors=10
uconf$min_dist = 0.5
uconf$spread = 1

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)

xy = mc2d_comp_cell_coord(mc_id = nm,graph_id = nm, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)

scdb_add_mc2d(mc2d_id, tgMC2D(nm, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph))

# mc = scdb_mc('cort_new')
tgconfig::set_param('mcell_mc2d_width', 2000, 'metacell')
tgconfig::set_param('mcell_mc2d_height', 2000, 'metacell')
tgconfig::set_param('mcell_mc2d_cex', 2, 'metacell')

mcell_mc2d_plot(mc2d_id = nm, edge_w = 1, colors = mc@colors, legend_pos = 'topleft', fn_suf = 'large')

mcell_mc2d_plot_by_factor(mc2d_id = nm, mat_id = nm, meta_field = 'day', colors = mc_new@colors, single_plot = F)


st_color_vec = color_key$color[match(max_st, color_key$st)]

mc = scdb_mc(nm)

# mc_new@colors = color_key3$color3[match(max_st, color_key3$Bonev_annotation)]
# names(mc_new@colors) = color_key3$Bonev_annotation[match(max_st, color_key3$Bonev_annotation)]
mat = scdb_mat(nm)
uu = sort(unique(mc@mc))

mc_mean_day = lapply(uu, function(u) mean(mat@cell_metadata[names(mc@mc)[mc@mc == u],'day'])) %>% unlist %>% as.numeric %>% round(digits = 3)

mc_metadata = rbind(as.vector(max_st), mc_mean_day, st_color_vec) %>% t %>% as_tibble
colnames(mc_metadata) = c('st', 'mean_day', 'color')
# rownames(mc_metadata) = 1:dim(mc_metadata)[[1]]
mc_metadata = tibble::rownames_to_column(mc_metadata)
colnames(mc_metadata) = c('mc', tail(colnames(mc_metadata), -1))
mc_metadata = cbind(mc_metadata, t(mc@n_bc))
md_cut = cut(mc_mean_day, breaks = seq(floor(min(mc_mean_day)), ceiling(max(mc_mean_day)), length.out = 13))
mc_metadata = mc_metadata %>% mutate(time_bin = match(md_cut, levels(md_cut))) %>% relocate(time_bin, .after = mean_day)
head(mc_metadata)

readr::write_tsv(mc_metadata, file = './BonevCollab/mcmd_cort_5.tsv', quote_escape = FALSE)

mcmd = vroom::vroom('./BonevCollab/mcmd_cort.tsv')

head(mcmd)

library(metacell)
scdb_init('./scdb', force_reinit = T)
nm = 'cort6'  ## function from Markus
  mat_id = nm
  mc_id = nm

  tag = nm
  m_0 = 0.006
  s_0 = 0.002
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
mc2d = scdb_mc2d(tag)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = Matrix::colSums(m@mat)
  s_tot = Matrix::colSums(m@mat[s_genes,])
  m_tot = Matrix::colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
# mc_cc
    p_coldens =densCols(x = s_score,y = m_score,colramp = colorRampPalette(c("lightgray","blue3", "red", "yellow")))
  #p_coldens =densCols(x = s_score,y = m_score,colramp = colorRampPalette(c("white",rev(inferno(5)))))
  #p_coldens =densCols(x = s_score,y = m_score)
  	
    dir.create(paste0('./figs/cc_', tag))
  png(sprintf("figs/cc_%s/cc_coldens_scores.png", tag), w=600, h=600)
  plot(s_score, m_score, pch=19, main = "S phase vs M phase UMIs",cex=0.4,
       xlab = "S phase score",ylab = "M phase score",
       col = p_coldens)
  dev.off()
  
  png(sprintf("figs/cc_%s/cc_scores.png", tag), w=600, h=600)
  plot(s_score, m_score, pch=19, main = "S phase vs M phase UMIs",cex=0.1,
       xlab = "S phase score",ylab = "M phase score")
  points(s_score[f], m_score[f], pch=19, cex=0.1, col="darkred")
  dev.off()
  shades = colorRampPalette(c("white","lightgray","lightblue","blue", "red", "yellow"))
  png(sprintf("figs/cc_%s/cc_sm_scores.png", tag), w=600, h=600)
  
  smoothScatter(s_score, m_score, colramp=shades, pch=19, cex=0.1,
                main = "S phase vs M phase UMIs",xlab = "S phase score",ylab = "M phase score")
  abline(a = m_0,b = - m_0/s_0)
  dev.off()
  
  
  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  png(sprintf("figs/cc_%s/2d_cc.png", tag), w=800, h=800)
  plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.4, col=ifelse(f[names(mc2d@sc_x)], "lightgray","black"))
  points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=2.5, bg=shades[101 - mc_cc])
  dev.off()
  
  png(sprintf("figs/cc_%s/bars_cc.png", tag), w=3600, h=1500)
  barplot(mc_cc, col=mc@colors, las=2, cex.names=0.7)
  dev.off()

mc_md_new_filt = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')

#   mc_md_new_filt = mc_md_new[mc_md_new$Comment == 'cortical',]

# mc_cc = get_mc_cc()
mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

colnames(mc_cc) = c('mc', 'cc_score')

mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))

md_cut = cut(mc_md_new_filt$mean_day, breaks = 13:18)
print(levels(md_cut))
print(table(md_cut))
prol_rate = 2
mult_by = 2 ** (prol_rate * (100 - mc_cc$cc_score)/100) 

frac_old = tapply(rep(1,nrow(mc_md_new_filt)), md_cut, function(x) x/sum(x)) %>% unsplit(md_cut)
frac_new = tapply(frac_old*mult_by, md_cut, function(x) x/sum(x), simplify = F) %>% unsplit(md_cut)
# frac_new = frac_old*mult_by
print(head(frac_old))
print(head(frac_new))
mc_cc$leak = pmax(0, (frac_old - frac_new)/frac_old)

mc_md_new_filt = dplyr::left_join(mc_md_new_filt, mc_cc, by = 'mc')

mc_cc_n = (100 - mc_cc$cc_score)/100

mc_cc_n

mc_size = table(mc@mc)/sum(table(mc@mc))

tapply(mc_cc_n, mc_md_new_filt$st, mean)

nsc_cc = tapply(mc_cc_n[mc_md_new_filt$st == 'NSC'], mc_md_new_filt$mean_day[mc_md_new_filt$st == 'NSC'], mean)

lss = lowess(as.numeric(names(nsc_cc)), nsc_cc, 4)

plot(as.numeric(names(nsc_cc)), nsc_cc)
lines(lss$x, lss$y, col = 'red')

hist(mc_cc$leak)

options(repr.plot.width = 18)
barplot(height = mc_cc$leak, col = mc_md_new_filt$color)
options(repr.plot.width = 6)

# mg = scdb_mgraph('cort')
# mg@mc_id = 'cort_new'
# scdb_add_mgraph('cort_new', mg)

mat = scdb_mat('cort5')

scdb_add_mat('cort6', mat)

library(lpsymphony)
library(metacell)
library(tidyverse)
library("Matrix")
scdb_init("scdb/", force_reinit=T)
nm = 'cort6'

get_mc_cc = function() {
  ## function from Markus
  mat_id = nm
  mc_id = nm

  tag = nm
  m_0 = 0.01
  s_0 = 0.005
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  return(mc_cc)
}

add_leak_to_md = function(){
  mc_md_new_filt = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')

#   mc_md_new_filt = mc_md_new[mc_md_new$Comment == 'cortical',]

  mc_cc = get_mc_cc()
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))

  md_cut = cut(mc_md_new_filt$mean_day, breaks = 13:18)
  print(levels(md_cut))
  print(table(md_cut))
  prol_rate = 2
  mult_by = 2 ** (prol_rate * (100 - mc_cc$cc_score)/100) 
  frac_old = tapply(rep(1,nrow(mc_md_new_filt)), md_cut, function(x) x/sum(x)) %>% unsplit(md_cut)
  frac_new = tapply(mult_by, md_cut, function(x) x/sum(x), simplify = F) %>% unsplit(md_cut)
  print(head(frac_old))
  print(head(frac_new))
  mc_cc$leak = pmax(0, (frac_old - frac_new)/frac_old)
  mc_md_new_filt = left_join(mc_md_new_filt, mc_cc, by = 'mc')

  write_tsv(x = mc_md_new_filt, file = './BonevCollab/mcmd_cort_6_leak.tsv', quote_escape = F)
}

build_singemb_net = function(mat_id,mc_id,mgraph_id,net_id,fig_dir,
                             age_field = "day",
                             mc_leak = NULL,
                             capacity_var_factor = NULL, 
                             t_exp = 1,T_cost = 1e+5,
                             flow_tolerance = 0.01,
                             network_color_ord = NULL,
                             mc_ord = NULL,
                             off_capacity_cost1 = 1,
                             off_capacity_cost2 = 1000,
                             k_norm_ext_cost = 1,
                             k_ext_norm_cost = 1,
                             k_ext_ext_cost = 1) {
  options(error = utils::recover)
  mat = scdb_mat(mat_id)
  mc  = scdb_mc(mc_id)
  mgraph = scdb_mgraph(mgraph_id)
  md = mat@cell_metadata
  cell_time = mat@cell_metadata$day %>% as.vector
  names(cell_time) = rownames(mat@cell_metadata)

  if(is.null(mc_leak)) {
    leak = rep(0,max(mc@mc))
  }
  if(is.null(capacity_var_factor)) {
    capacity_var_factor = rep(0.25,max(mc@mc))
  }
  
  if(is.null(mc_leak)) {
    leak = rep(0, max(mc@mc))
  } else {
    leak = mc_leak
  }

  
  mcell_new_mctnetwork(net_id = net_id,
                       mc_id = mc_id,
                       mgraph_id = mgraph_id,
                       cell_time = cell_time)
  mct = scdb_mctnetwork(net_id)

  #computing manifold costs (based on mgraph distances)
  mct = mctnetwork_comp_manifold_costs(mct,t_exp=t_exp, T_cost=T_cost)
  message("computed manifold costs")
  
  #generating network structure	
  mct = mctnetwork_gen_network(mct, mc_leak = leak,capacity_var_factor = capacity_var_factor,
                               k_norm_ext_cost = k_norm_ext_cost,k_ext_norm_cost = k_ext_norm_cost,k_ext_ext_cost = k_ext_ext_cost,
                               off_capacity_cost1 = off_capacity_cost1,off_capacity_cost2 = off_capacity_cost2)	
  message("generated network")
  
  #solving the flow problem
  mct = mctnetwork_gen_mincost_flows(mct, flow_tolerance = flow_tolerance)
  message("solved network flow problem")
  
  #compute propagatation forward and background
  mct = mctnetwork_comp_propagation(mct)
  
  #adding back the object with the network and flows
  scdb_add_mctnetwork(paste0(net_id, '_comp'), mct)
  
  # mct = scdb_mctnetwork(paste0(net_id, '_comp'))
  
  source('~/raid/proj/mmcortex/scripts/mctnetwork_plot_net.r')

  mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mcmd_cort_6_leak.tsv')
    
  color_key = unique(mc_md[,c('st', 'color')]) 
  color_key = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('st')
  color_key = color_key[c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC','iCfuPN',
                      'iCPN/CfuPN','iCPN_L2-3','CPN_L2-3','CPN_L5-6','CthPN','SCPN'),]

  mctnetwork_plot_net(mct_id = paste0(net_id, '_comp'), colors_ordered = color_key$color,
                      fn=sprintf("%s/%s_net.png",fig_dir,net_id), 
                      # mc_ord = ord_vec,
                      h = 2500,w = 2000)
    
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  message("plotted the network")
}

main_build_network = function() {
  options(error = utils::recover)
  age_field = "day"
  mat_id = nm
  mc_id = nm
  mgraph_id = nm
  net_id = nm
  fig_dir = paste0("./figs/", nm, ".net")

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  mc = scdb_mc(mc_id)
  capacity_var_factor = rep(0.25,ncol(mc@e_gc))
  
  # next define the mc_leak parameter
  
  add_leak_to_md()

  mc_md = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_cort_6_leak.tsv')
  mc_leak = mc_md$leak

    
  build_singemb_net(mat_id = mat_id,
                    mc_id = mc_id,
                    mgraph_id = mgraph_id,
                    net_id = net_id,
                    fig_dir = fig_dir,
                    age_field = age_field,
                    mc_leak = mc_leak,
                    # mc_ord = mc_ord,
                    capacity_var_factor = capacity_var_factor,
                    k_norm_ext_cost = 2,
                    k_ext_norm_cost = 2,
                    k_ext_ext_cost = 100,
                    flow_tolerance = 0.01,
                    )
  
  
}



main_build_network()

# ## Plot custom network

# library(dplyr)
# library(metacell)

# net_id = 'cort_new_2'
# source('~/raid/proj/mmcortex/scripts/mctnetwork_plot_net.r')
#  mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mcmd_cort_070221_leak.tsv')
#     fig_dir = './figs/cort_new_2.net'
#   color_key = unique(mc_md[,c('st', 'color')]) 
#   color_key = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('st')
# #   color_key = color_key[c('Astrocytes','early_NSC?','OPC','NSC','nNSC?','early_nNSC?','IPC_cyc', 'IPC','iCfuPN','iCPN_L2-3',
# #                   'iCPN_L5-6','CPN_L2-3','CPN_L5-6','SCPN','CthPN','Stellate_L4'),]
# color_key = color_key[c('OPC','Astrocytes','NSC','nNSC?','IPC_cyc', 'IPC','iCfuPN',
#                       'iCPN/CfuPN','iCPN_L2-3','CPN_L2-3','CthPN','SCPN','CPN_L5-6'),]

#   mctnetwork_plot_net(mct_id = paste0(net_id, '_comp'), colors_ordered = color_key$color,
#                       fn=sprintf("%s/%s_net_large.png",fig_dir,net_id), 
#                       # mc_ord = ord_vec,
#                       h = 5000,w = 4000)

library(lpsymphony)
library(metacell)
library(tidyverse)
library("Matrix")
scdb_init("scdb/", force_reinit=T)
nm = 'cort6'

get_mc_cc = function() {
  ## function from Markus
  mat_id = nm
  mc_id = nm

  tag = nm
  m_0 = 0.01
  s_0 = 0.005
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  return(mc_cc)
}

add_leak_to_md = function() {
  mc_md_new_filt = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')

#   mc_md_new_filt = mc_md_new[mc_md_new$Comment == 'cortical',]

  mc_cc = get_mc_cc()
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))

#   md_cut = cut(mc_md_new_filt$mean_day, breaks = 13:18)
#   print(levels(md_cut))
#   print(table(md_cut))
#   prol_rate = 2
#   mult_by = 2 ** (prol_rate * (100 - mc_cc$cc_score)/100) 
#   frac_old = tapply(rep(1,nrow(mc_md_new_filt)), md_cut, function(x) x/sum(x)) %>% unsplit(md_cut)
#   frac_new = tapply(mult_by, md_cut, function(x) x/sum(x), simplify = F) %>% unsplit(md_cut)
#   print(head(frac_old))
#   print(head(frac_new))
#   mc_cc$leak = pmax(0, (frac_old - frac_new)/frac_old)
  mc_cc$leak = rep(0, nrow(mc_md_new_filt))
  mc_md_new_filt = left_join(mc_md_new_filt, mc_cc, by = 'mc')

  write_tsv(x = mc_md_new_filt, file = './BonevCollab/mcmd_cort_6_leak_exp.tsv', quote_escape = F)
}

build_singemb_net = function(mat_id,mc_id,mgraph_id,net_id,fig_dir,
                             age_field = "day",
                             mc_leak = NULL,
                             capacity_var_factor = NULL, 
                             t_exp = 1,T_cost = 1e+5,
                             flow_tolerance = 0.01,
                             network_color_ord = NULL,
                             mc_ord = NULL,
                             off_capacity_cost1 = 10,
                             off_capacity_cost2 = 1000,
                             k_norm_ext_cost = 1,
                             k_ext_norm_cost = 1,
                             k_ext_ext_cost = 1) {
  options(error = utils::recover)
  mat = scdb_mat(mat_id)
  mc  = scdb_mc(mc_id)
  mgraph = scdb_mgraph(mgraph_id)
  md = mat@cell_metadata
  cell_time = mat@cell_metadata$day %>% as.vector
  names(cell_time) = rownames(mat@cell_metadata)

  if(is.null(mc_leak)) {
    leak = rep(0,max(mc@mc))
  }
  if(is.null(capacity_var_factor)) {
    capacity_var_factor = rep(0.25,max(mc@mc))
  }
  
  if(is.null(mc_leak)) {
    leak = rep(0, max(mc@mc))
  } else {
    leak = mc_leak
  }

  
  mcell_new_mctnetwork(net_id = net_id,
                       mc_id = mc_id,
                       mgraph_id = mgraph_id,
                       cell_time = cell_time)
  mct = scdb_mctnetwork(net_id)

  #computing manifold costs (based on mgraph distances)
  mct = mctnetwork_comp_manifold_costs(mct,t_exp=t_exp, T_cost=T_cost)
  message("computed manifold costs")
  
  #generating network structure	
  mct = mctnetwork_gen_network(mct, mc_leak = leak,capacity_var_factor = capacity_var_factor,
                               k_norm_ext_cost = k_norm_ext_cost,k_ext_norm_cost = k_ext_norm_cost,k_ext_ext_cost = k_ext_ext_cost,
                               off_capacity_cost1 = off_capacity_cost1,off_capacity_cost2 = off_capacity_cost2)	
  message("generated network")
  
  #solving the flow problem
  mct = mctnetwork_gen_mincost_flows(mct, flow_tolerance = flow_tolerance)
  message("solved network flow problem")
  
  #compute propagatation forward and background
  mct = mctnetwork_comp_propagation(mct)
  
  #adding back the object with the network and flows
  scdb_add_mctnetwork(paste0(net_id, '_exp_comp'), mct)
  
  # mct = scdb_mctnetwork(paste0(net_id, '_comp'))
  
  source('~/raid/proj/mmcortex/scripts/mctnetwork_plot_net.r')

  mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mcmd_cort_6_leak_exp.tsv')
    
  color_key = unique(mc_md[,c('st', 'color')]) 
  color_key = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('st')
  color_key = color_key[c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC',
                      'iCPN/CfuPN','iCPN','iCPN_L2-3','CPN_L2-3','CPN_L5-6','iCfuPN','CthPN','SCPN'),]

  mctnetwork_plot_net(mct_id = paste0(net_id, '_exp_comp'), colors_ordered = color_key$color,
                      fn=sprintf("%s/%s_net_exp_leak.png",fig_dir,net_id), 
                      # mc_ord = ord_vec,
                      h = 2500,w = 2000)
    
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  message("plotted the network")
}

main_build_network = function() {
  options(error = utils::recover)
  age_field = "day"
  mat_id = nm
  mc_id = nm
#   mgraph_id = nm
  mgraph_id = 'cort6'
  net_id = nm
  fig_dir = paste0("./figs/")

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  mc = scdb_mc(mc_id)
  capacity_var_factor = rep(0.05,ncol(mc@e_gc))
  
  # next define the mc_leak parameter
  
  add_leak_to_md()

  mc_md = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_cort_6_leak_exp.tsv')
  mc_leak = mc_md$leak

    
  build_singemb_net(mat_id = mat_id,
                    mc_id = mc_id,
                    mgraph_id = mgraph_id,
                    net_id = net_id,
                    fig_dir = fig_dir,
                    age_field = age_field,
                    mc_leak = mc_leak,
                    # mc_ord = mc_ord,
                    capacity_var_factor = capacity_var_factor,
                    k_norm_ext_cost = 1,
                    k_ext_norm_cost = 1,
                    k_ext_ext_cost = 10,
                    flow_tolerance = 0.05,
                    )
  
  
}



main_build_network()

# mc = scdb_mc('cort6')

# mba = as.matrix(Matrix::readMM('./data//MBA_clust_data.mtx'))

# mba_col = tibble::column_to_rownames(read.delim('./data/MBA_clust_col_att.tsv', sep='\t'), 'X')

# mba_genes = read.delim('./data/MBA_clust_genes.tsv', sep='\t', header = F)
# mba_genes = gsub('\'', '', unlist(mba_genes))
# names(mba_genes) = NA

# gb = intersect(mba_genes, rownames(mc@e_gc))

# feats = scdb_gset('cort6_feats_f')

# feats = names(feats@gene_set)

# gb_feats = intersect(feats, gb)

# rownames(mba) = mba_genes

# head(mba)

# mba_mc_cor = tgs_cor(mba[gb_feats,], mc@e_gc[gb_feats,], spearman = T)

# max_ind = apply(mba_mc_cor,2, which.max)


# max_cl = mba_col['Subclass',max_ind]

# head(max_cl)

# mcmd = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')

# sc_st_df = data.frame(t(rbind(max_cl, mcmd$st)))
# colnames(sc_st_df) = c('sc', 'st')

# tapply(sc_st_df$sc, sc_st_df$st, table)

# abi = read.delim('./data/allen_icx_hcf_cluster_medians.csv', sep = ',', header = 1)

# abi = as.matrix(tibble::column_to_rownames(abi, 'feature'))

# colnames(abi) = purrr::map(stringr::str_split(colnames(abi), '_'), 2)

# sort(colnames(abi))

# head(abi)

# gb = intersect(rownames(abi), rownames(mc@e_gc))
# gb_feats = intersect(feats, gb)

# abi_mc_cor = tgs_cor(abi[gb_feats,], mc@e_gc[gb_feats,], spearman = T)

# max_ind = apply(abi_mc_cor, 2, which.max)

# max_cl = colnames(abi)[max_ind]

# sc_st_df = data.frame(t(rbind(max_cl, mcmd$st)))

# colnames(sc_st_df) = c('sc', 'st')

# tapply(sc_st_df$sc, sc_st_df$st, table)

# hist(apply(abi_mc_cor, 2, max))

feats = scdb_gset('cort5_feats_f')
mc = scdb_mc('cort6')
mc2d = scdb_mc2d('cort6')
mcmd = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')
# mgi = vroom::vroom('~/raid/gene_association.mgi', skip = 3)

# colnames(mgi) = c('Database','Designation','MGI Marker Accession ID','Mouse Marker Symbol','NOT Designation',
#                     'GO Term ID', 'MGI Reference Accession ID', 'GO Evidence Code', 'Inferred From','Ontology',
#                   'Mouse Marker Name','Mouse Marker Synonyms','Mouse Marker Type','Taxon','Modification Date',
#                   'Assigned By','Annotation Extension','Gene Product')

# mig_genes = unique(unlist(mgi[grep('0001764', unlist(mgi[,'NOT Designation'])),'MGI Marker Accession ID']))



mat = scdb_mat('cort6')

cust_st_ord = c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN',
                      'iCPN_L2-3','CPN_L2-3','CPN_L5-6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))

feats_cor = tgs_cor(t(mc@e_gc[names(feats@gene_set),]), spearman = T)

km_feats = tglkmeans::TGL_kmeans(feats_cor, k = 32, seed = SEED)

saveRDS(km_feats,'./data/km_feats.rds')

p = pheatmap::pheatmap(feats_cor[order(km_feats$cluster),order(km_feats$cluster)], 
                       fontsize_col = 4, fontsize_row = 4,
                        cluster_cols = F, 
                        cluster_rows = F,silent = T
                      )


save_pheatmap_png(p, './figs/test_feats_cor.png', 2500, 2500, 200)


clust_sum = tgs_matrix_tapply(t(mc@e_gc[names(feats@gene_set),]), km_feats$cluster, mean)

pca_feat_cl = princomp(x = t(clust_sum))

# options(repr.plot.height = 12, repr.plot.width = 12)
png('./figs/feat_clust_pca.png', h = 1200, w = 1200, res = 150)
plot(pca_feat_cl$scores[,1], pca_feat_cl$scores[,2], col = mcmd$color, pch = 16, xlab = 'PC1 - 39.8% of variance', ylab = 'PC2 - 17.5% of variance',
     main = 'PCA of feature gene cluster values',
     cex = 0.8,
    )
points(pca_feat_cl$scores[,1], pca_feat_cl$scores[,2], col = 'black', 
       cex = 0.8,
      )
dev.off()
# options(repr.plot.height = 6, repr.plot.width = 6)

pca_feat_cl$sdev/sum(pca_feat_cl$sdev)

head(clust_sum)

legc = log2(mc@e_gc + 1e-06)
legc_norm = t(apply(legc, 1, function(x) x - median(x)))

clust_sum = tgs_matrix_tapply(t(legc_norm[names(feats@gene_set),]), km_feats$cluster, mean)

clust_sum_st = t(tgs_matrix_tapply(clust_sum, mcmd$st, mean))

plot_mat = clust_sum_st
min_val = min(plot_mat)
max_val = max(plot_mat)
# print(min_val)
# print(max_val)
# 0+(max_val-min_val)/51
p = pheatmap(plot_mat[,cust_st_ord], main = 'Mean feature module enrichment by type', cluster_cols = F,
                     breaks = c(seq(0, max_val/2, l=50),
                               seq(max_val/2+(max_val-max_val/2)/51, max_val, l=51)),
                 color = colorRampPalette(c('white', 'yellow', 'red'))(100))
save_pheatmap_png(p, './figs/feat_module_resid.png')

dir.create('./figs/feats_mc2d')

clrmp = colorRampPalette(c('white', 'yellow', 'red', 'brown', 'black'))(100)

mc_mean_day = tapply(mat@cell_metadata[names(mc@mc),'day'], mc@mc, mean)

cor_clust_mc_day = tgs_cor(t(clust_sum), matrix(mc_mean_day), spearman = T)
setNames(cor_clust_mc_day[order(abs(cor_clust_mc_day), decreasing = T)], rownames(cor_clust_mc_day)[order(abs(cor_clust_mc_day), decreasing = T)])

cor_clust_st_ord = tgs_cor(t(clust_sum), matrix(cust_mc_ord_st), spearman = T)
setNames(cor_clust_st_ord[order(abs(cor_clust_st_ord), decreasing = T)], 
         rownames(cor_clust_st_ord)[order(abs(cor_clust_st_ord), decreasing = T)])

sort(names(feats@gene_set)[km_feats$cluster %in% c(30,16,10)])

dir.create('./figs/cor_feat_clust_mc_day/')

sapply(1:nrow(clust_sum), function(x) {
    png(glue::glue('./figs/cor_feat_clust_mc_day/{x}.png'), width = 1200, h = 1200, res = 150)
    plot(clust_sum[x,], mc_mean_day, pch = 16, col = mcmd$color, main = glue::glue('clust {x} - cor = {round(cor_clust_mc_day[[x]], digits = 3)}'))
    points(clust_sum[x,], mc_mean_day, col = 'black')
    dev.off()
})

plot_gene_clust = function(cl) {

    clr_vals = clrmp[round(100*((clust_sum[cl,] - min(clust_sum[cl,]))/(max(clust_sum[cl,]) - min(clust_sum[cl,]))))]
    # clr_vals = clrmp[round(100*(mc_mean_day - min(mc_mean_day, na.rm = T))/(max(mc_mean_day, na.rm = T) - min(mc_mean_day, na.rm = T)))]
    png(glue::glue('./figs/feats_mc2d/{cl}.png'), width = 1200, height = 1200, res = 150)

    plot(mc2d@sc_x, mc2d@sc_y, col = scales::alpha(mcmd$color[mc@mc[names(mc2d@sc_x)]], 0.2), 
         main = glue::glue("cluster {cl}"),
         pch = 16, cex = 0.3,
        )
    points(mc2d@mc_x, mc2d@mc_y, col = clr_vals, pch = 16)
        points(mc2d@mc_x, mc2d@mc_y, col = 'black', lwd = 0.5,
         )
    dev.off()
}

sapply(1:nrow(clust_sum), plot_gene_clust)

plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6))
{
  if (!is.null(fig_fn)) {
#     .plot_start(fig_fn, 400, 400)
      png(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

  if (is.null(show_vals_ind)) {
    show_vals_ind = rep(T, length(cols))
  }
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=round(vals[show_vals_ind], 3), pos=4)
#   text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
}

sapply(1:nrow(clust_sum), function(x) plot_color_bar(seq(min(clust_sum[x,]),max(clust_sum[x,]),l=101), clrmp, fig_fn = glue::glue('./figs/feats_mc2d/colorbar_{x}.png')))

cloi = feats@gene_set[km_feats$cluster == 13]
cloi[order(names(cloi))]

nm = 'cort6'

mcmd = vroom::vroom('./BonevCollab//mcmd_cort_5.tsv')
# readr::write_tsv(mcmd, './BonevCollab//mcmd_cort_6.tsv')

# mcmd = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')
color_key = unique(dplyr::select(mcmd, st, color))

color_key

mc = scdb_mc(nm)
mct = scdb_mctnetwork(paste0(nm, '_comp'))

net = mct@network[mct@network$flow > 0,]

net = net[net$mc1 > 0 & 
          net$mc2 > 0 & 
          net$time1 > 0 & 
          net$time2 < max(net$time2) &
          net$time1 != net$time2 &
          net$mc1 != net$mc2,]

net[,c('mc1','mc2')] = apply(net[,c('mc1','mc2')], 2, as.numeric)

bad_icpn = unique(net$mc1[net$mc1 %in% mcmd$mc[mcmd$st == 'iCPN/CfuPN'] & net$mc2 %in% mcmd$mc[mcmd$st == 'IPC']])
outflow_bad_icpn = unique(net$mc2[net$mc1 %in% bad_icpn])
unique(mcmd$st[outflow_bad_icpn])

bad_ipc = unique(net$mc1[net$mc1 %in% mcmd$mc[mcmd$st == 'IPC'] & net$mc2 %in% mcmd$mc[mcmd$st == 'IPC_cyc']])
bad_ipc
outflow_bad_ipc = unique(net$mc2[net$mc1 %in% bad_ipc])
unique(mcmd$st[outflow_bad_ipc])

# all_to_cfupn = setNames(unique(net[net$mc2 %in% mcmd$mc[mcmd$st %in% c('CthPN', 'SCPN')],'mc1']), 
#                         mcmd$st[unique(net[net$mc2 %in% mcmd$mc[mcmd$st %in% c('CthPN', 'SCPN')],'mc1'])])

# all_to_cfupn[!(names(all_to_cfupn) %in% c('SCPN', 'CthPN'))]

# icfupn_target_st = sapply(all_to_cfupn, function(mc) mcmd$st[net$mc2[net$mc1 == mc]])
# names(icfupn_target_st) = all_to_cfupn
# good_icfupn = sort(as.numeric(names(icfupn_target_st)[sapply(icfupn_target_st, function(x) !any(grepl('^CPN', x)))]))
# good_icfupn = good_icfupn[good_icfupn %in% mcmd$mc[mcmd$st == 'iCfuPN']]

set1 = mcmd$mc[mcmd$st == 'iCfuPN']
set2 = unique(net$mc1[net$mc2 %in% mcmd$mc[mcmd$st == 'SCPN']])
set3 = unique(net$mc1[net$mc2 %in% mcmd$mc[mcmd$st == 'CthPN']])
set4 = unique(net$mc1[net$mc2 %in% mcmd$mc[mcmd$st == 'iCPN/CfuPN']])
first_order = intersect(set1, intersect(set2[!(set2 %in% set4)], set3[!(set3 %in% set4)]))
first_order = c(first_order, c(161, 186))
# second_order = good_icfupn[!(good_icfupn %in% set4)]
# good_icfupn

bad_icpnl23 = unique(net$mc1[net$mc1 %in% mcmd$mc[mcmd$st == 'iCPN_L2-3'] & net$mc2 %in% mcmd$mc[mcmd$st == 'CPN_L2-3']])
length(bad_icpnl23)
length(mcmd$mc[mcmd$st == 'iCPN_L2-3'])

mcmd_new = mcmd

color_key = unique(mcmd[,c('st', 'color')])

mcmd_new$st[mcmd$st == 'iCPN/CfuPN'] = 'iCPN'
mcmd_new$color[mcmd_new$st == 'iCPN'] = color_key$color[color_key$st == 'iCPN/CfuPN']

mcmd_new$st[mcmd$st == 'iCfuPN' & !(mcmd$mc %in% first_order)] = 'iCPN/CfuPN'

mcmd_new$color[mcmd$mc %in% first_order] = gplots::col2hex('lightskyblue')

mcmd_new$st[mcmd$mc %in% bad_icpn] = 'IPC'
mcmd_new$color[mcmd$mc %in% bad_icpn] = color_key$color[color_key$st == 'IPC']

# mcmd_new$st[mcmd$mc %in% bad_ipc] = 'IPC_cyc'
# mcmd_new$color[mcmd$mc %in% bad_ipc] = color_key$color[color_key$st == 'IPC_cyc']

mc@colors = mcmd_new$color
scdb_add_mc('cort6', mc)

source('./scripts/mctnetwork_plot_net.r')

color_key = unique(mcmd_new[,c('st', 'color')])

color_key

color_key_new = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('st')
  color_key_new = color_key_new[c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC',
                      'iCPN/CfuPN','iCPN','iCPN_L2-3','CPN_L2-3','CPN_L5-6','iCfuPN','SCPN','CthPN'),]

color_key_new

readr::write_tsv(x = mcmd_new, file = './BonevCollab//mcmd_cort_6.tsv')

mctnetwork_plot_net(mct_id = 'cort6_comp', 
#                     mc_ord = new_mc_ord,
                    colors_ordered = color_key_new$color,
                      fn=sprintf("%s/%s_net.png",fig_dir='./figs',net_id = 'cort6_comp'), 
                      # mc_ord = ord_vec,
                      h = 2500,w = 2000)

mcell_mc2d_plot(mc2d_id = 'cort6', colors = mcmd_new$color)

mcell_mc2d_plot(mc2d_id = 'cort6', filt_mc = 
                ifelse(1:nrow(mcmd) %in% good_icfupn[good_icfupn %in% mcmd$mc[mcmd$st == 'iCfuPN']], TRUE, FALSE),
               fn_suf = 'test_good_icfupn')



ipc_raw = mcmd$mc[mcmd$st == 'IPC' | mcmd$st == 'IPC_cyc']

ipc_cyc_mcs = ipc_raw[mc_cc[ipc_raw] < 80]
ipc_cyc_mcs
mcmd$st[ipc_cyc_mcs] = 'IPC_cyc'

df_ann = data.frame(data.table::rbindlist(
                                            lapply(unique(mg$mc1), 
                                                   function(m) {st_old = mcmd$st[m];
                                                        tbl = table(mcmd$st[mg$mc2[mg$mc1 == m & mg$mc2 != m]]);
                                                        st_new = names(tbl[order(tbl, decreasing = T)])[[1]]
                                                        return(list('st_old' = st_old, 'st_new' = st_new))
                                                                }
                                                   )
                                            )
                      )


head(df_ann)

df_ann$st_new[df_ann$st_old == 'OPC'] = 'OPC'
df_ann$st_new[df_ann$st_old == 'Astrocytes'] = 'Astrocytes'
df_ann$st_new[df_ann$st_old == 'IPC_cyc'] = 'IPC_cyc'
df_ann$st_new[is.na(df_ann$st_new)] = df_ann$st_old[is.na(df_ann$st_new)]
df_ann$st_new[df_ann$st_new == 'IPC_cyc'] = df_ann$st_old[df_ann$st_new == 'IPC_cyc']

length(which(df_ann[,1] != df_ann[,2]))
head(df_ann[which(df_ann[,1] != df_ann[,2]),])

mcell_mc2d_plot(mc2d_id = nm)

mcell_mc2d_plot(mc2d_id = nm, colors = color_key$color[match(df_ann$st_old, color_key$st)], fn_suf = 'no_reannot')

mcmd$st = df_ann$st_new
mcmd$color = color_key$color[match(mcmd$st, color_key$st)]

readr::write_tsv(mcmd, './BonevCollab//mcmd_cort_5_iter1.tsv')

mc = scdb_mc(nm)

mc@colors = mcmd$color

scdb_add_mc(nm, mc)

# mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mcmd_cort_5_iter1.tsv')
mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab/mcmd_cort_6.tsv')

head(mc_md)

    
color_key = unique(mc_md[,c('st', 'color')]) 
color_key = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('st')
color_key = color_key[c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC',
                      'iCPN/CfuPN','iCPN','iCPN_L2-3','CPN_L2-3','CPN_L5-6','iCfuPN','SCPN','CthPN'),]
color_key$ord = 1:nrow(color_key)
color_key = tibble::rownames_to_column(color_key)
colnames(color_key) = c('cell_type', 'color', 'i', 'ord')

color_key

readr::write_tsv(color_key, '/net//mraid14//export//tgdata//users/aviezerl/proj/mmcortex/cell_type_annot_YS.tsv', quote_escape = F)

ys = vroom::vroom('/net//mraid14//export//tgdata//users/aviezerl/proj/mmcortex/cell_type_annot_YS.tsv')

avz = vroom::vroom('/net//mraid14//export//tgdata//users/aviezerl/proj/mmcortex/cell_type_annot.tsv')

ys

mcmd = vroom::vroom('./BonevCollab//mcmd_cort_5.tsv')
color_key = unique(mcmd[,c('st', 'color')])



# color_key = color_key[c('Astrocytes','early_NSC?','OPC','NSC','nNSC?','early_nNSC?','IPC','iCfuPN','iCPN_L2-3',
#               'iCPN_L5-6','CPN_L2-3','CPN_L5-6','SCPN','CthPN','Stellate_L4'),]
df = data.frame(color_key)

# st_for_plot = df$st
# st_for_plot['iCPN_L2-3'] = 'immature CPN L2-3'
# st_for_plot['iCfuPN'] = 'immature CPN L2-3'


l = nrow(df)
scale_y = 2
png('./figs/legend_cort5_new.png', width = 1000, height = 750, res=250)
plot(rep(0.93,l), scale_y*seq(l,1,-1), pch = 16, cex = 1, col = df$color, ylim = c(0.5,scale_y*l+1),
    xlab = '', 
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(rep(0.94,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 0.5, df$st)
dev.off()

mcell_mc2d_plot_by_factor(mc2d_id = 'cort', mat_id = 'cort', meta_field = 'day', 
                          single_plot = F)

mcmd3 = vroom::vroom('./BonevCollab//mcmd_cort_new.tsv')

st_inds = lapply(sort(unique(mcmd3$st)), function(x) which(mcmd3$st == x))
names(st_inds) = sort(unique(mcmd3$st))
st_inds = st_inds[order(unlist(lapply(st_inds, length)), decreasing = T)]

save_pheatmap_png <- function(x, filename, width=1000, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

nm = 'cort5'
mc = scdb_mc(nm)

mcmd = vroom::vroom('./BonevCollab//mcmd_cort_5_iter1.tsv')

cust_st_ord = c('OPC','Astrocytes','NSC','IPC_cyc', 'IPC','iCfuPN',
                      'iCPN/CfuPN','iCPN_L2-3','CPN_L2-3','CthPN','SCPN','CPN_L5-6')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))

mcell_mc_plot_marks(mc_id = nm, gset_id = glue::glue('{nm}_markers_f'), mat_id = nm, fig_fn = glue::glue('./figs/{nm}_marks_hm.png'), mc_ord = cust_mc_ord_st)

marks = scdb_gset(glue::glue('{nm}_markers_f'))
marks = names(marks@gene_set)

sapply(marks, function(g) mcell_mc2d_plot_gene(mc2d_id = nm, gene = g, show_legend = F))

tfoi = unlist(read.delim(file = './data/tfoi.txt'))

# tfoi
# marks = names(marks@gene_set)
sapply(tfoi[!(tfoi %in% marks)], function(g) mcell_mc2d_plot_gene(mc2d_id = nm, gene = g, show_legend = F))

mcell_mc2d_plot_gene(mc2d_id = nm, gene = 'Nhlh1', )

max_st = tgs_matrix_tapply(mc@mc_fp, mcmd$st, max)

max_st_norm = apply(max_st, 1, function(x) x/max(x))

st_genes_inds = apply(max_st_norm, 2, function(x) rownames(max_st_norm)[head(order(x, decreasing = T), 5)])

st_genes_ord = apply(max_st_norm, 2, function(x) rownames(max_st_norm)[order(x, decreasing = T)])

st_tfs = apply(st_genes_ord, 2, function(x) x[x %in% tfoi])

head(st_tfs)

st_genes_inds

sapply(colnames(st_genes_inds), function(col) sapply(st_genes_inds[,col], function(x) x))

sort(unique(unlist(data.frame(st_genes_inds))))

head(max_st_norm)

dir.create('./figs/cort_substruct')

marks = scdb_gset('cort_markers_f')

plot_substruct = function(st, n, i, mc, mat, marks) {
    fig_dir = file.path(wd, 'figs', 'cort_substruct', paste(i, n, sep = '_'))
    if (!dir.exists(fig_dir)) {dir.create(fig_dir)}
    sc_inds = names(mc@mc[mc@mc %in% st])
    ext_meta = as.character(mat@cell_metadata[sc_inds,'day'])
    names(ext_meta) = sc_inds
    if (length(st) >= 4) {
        cor_fp = tgs_cor(mc@mc_fp[,st], spearman = T)
        p = pheatmap(cor_fp, silent = T, title = n)
        save_pheatmap_png(x = p, 
                          filename = file.path(fig_dir, 'mc_cor_hm.png'), 
                          height = 1200, width = 1200, res=75)
        mcell_mc_plot_marks(mc_id = 'cort_new', 
                            gset_id = 'cort_markers_f', 
                            mat_id = 'cort_new', 
                            fig_fn = file.path(fig_dir, 'marks_hm.png'), 
                            add_metadata = 'day',
#                             ext_metadata = ext_meta, 
                            focus_mcs = st)
    } else {
        sc_inds = names(mc@mc[mc@mc %in% st])
        cor_sc = tgs_cor(as.matrix(mat@mat[,sc_inds]), spearman = T)
        p = pheatmap(cor_sc, silent = T)
        save_pheatmap_png(x = p, 
                          filename = file.path(fig_dir, 'sc_cor_hm.png'), 
                          height = 1200, width = 1200, res=75)
        if (length(st) > 1) {
            mcell_mc_plot_marks(mc_id = 'cort_new', 
                                gset_id = 'cort_markers_f', 
                                mat_id = 'cort_new', 
                                add_metadata = 'day',
#                                 ext_metadata = ext_meta, 
                                fig_fn = file.path(fig_dir, 'marks_hm.png'), 
                                focus_mcs = as.vector(st))
        } else {
            mat_log = log2(mat@mat[names(marks@gene_set),sc_inds]+0.1) %>% as.matrix
            mat_norm = mat_log - rowMeans(mat_log)
            p2 = pheatmap(mat_norm, silent = T)
            save_pheatmap_png(x = p2, 
                          filename = file.path(fig_dir, 'marks_hm.png'), 
                          height = 1200, width = 1200, res=75)
        }
    }
}

mc = scdb_mc('cort_new')
mat = scdb_mat('cort_new')

lapply(seq_along(st_inds), function(x, n, i) {plot_substruct(x[[i]], n[[i]], i, mc, mat, marks)}, x=st_inds, n=names(st_inds))

mcell_mc2d_plot_by_factor(mc2d_id = 'cort', mat_id = 'cort_new', meta_field = 'st_cort', single_plot = T, colors = mcmd3$color)
mcell_mc2d_plot_by_factor(mc2d_id = 'cort', mat_id = 'cort_new', meta_field = 'st_cort', single_plot = F, colors = mcmd3$color)

mat@cell_metadata$st_cort = mcmd3$st[mc@mc[rownames(mat@cell_metadata)]]

scdb_add_mat('cort', mat)

length(mat@cell_metadata$st_cort[!is.na(mat@cell_metadata$st_cort)])

tail(mat@cell_metadata, 30)

library(officer)
ppt = read_pptx()

print_slide = function(n, i, ppt) {
    path_hm = file.path(wd, 'figs', 'cort_substruct', paste(i, n, sep='_'), 'marks_hm.png')
    path_2d = file.path(wd, 'figs', 'cort.by_st_cort', paste0(paste('cort.2d_proj', n, sep='_'), '.png'))
#     print(path_2d)
    img_hm = external_img(src = path_hm, height = 7, width = 5)
    img_2d = external_img(src = path_2d, height = 4.6, width = 4.6)
#     ppt = add_slide(ppt)
    ppt = add_slide(ppt, layout = "Two Content", master = "Office Theme")
    ppt = ph_with(x = ppt, value = img_hm, location = ph_location(left = 4.75, 
                                    top = 0.5, width = 5, height = 7))
    ppt = ph_with(x = ppt, value = img_2d, location = ph_location(left = 0.25, 
                            top = 2, width = 4.2, height = 4.2))
}

lapply(seq_along(st_inds), function(n, i) {print_slide(n[[i]], i, ppt)}, n=names(st_inds))

print(ppt, target = file.path(wd, 'figs', 'Reannot_cort.pptx'))

mc_subset_plot_gene = function(mc, mc2d, mc_md, type, gene) {
    st_fold = file.path('figs', 'st_marks_cort', type)
    if (!dir.exists(st_fold)) {dir.create(st_fold)}
    png(paste0(file.path(st_fold, gene), '.png'))
    mcs = select(filter(mc_md, st == type), 'mc') %>% unlist
    mc_x = mc2d@mc_x[names(mc2d@mc_x) %in% mcs]
    mc_y = mc2d@mc_y[names(mc2d@mc_y) %in% mcs]
    sc_x = mc2d@sc_x[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% mcs])]
    sc_y = mc2d@sc_y[names(mc2d@sc_y) %in% names(mc@mc[mc@mc %in% mcs])]
    plot(sc_x, sc_y, pch = 16, col = 'gray', main = gene)
    col_map = colorRampPalette(colors = c('white', 'red'))(100)
    vals = mc@mc_fp[gene,mcs]
    vals_t = 100*(vals - min(vals))/(max(vals) - min(vals))
    points(mc_x, mc_y, pch = 16, cex = 3, col = col_map[vals_t])
    text(mc_x, mc_y, mcs, cex = 0.7, col = 'black')
    dev.off()
}

mc = scdb_mc('cort')
mc2d = scdb_mc2d('cort')
mcmd = vroom::vroom('./BonevCollab//mcmd_cort.tsv')
marks = scdb_gset('cort_markers_f')
feats = scdb_gset('cort_feats_f')


# head(mcfp, 10)

color_key = unique(select(mcmd, st, color))
color_key = rbind(color_key, list(st = 'IPC_cyc', color = gplots::col2hex('olivedrab1')))
color_key

ipc_cc = which((mc_cc <= 80) & (mcmd$st == 'IPC'))
icfupn = c(109, 447, 448, 456)
icpn_l23 = c(175, 181, 343)
cpn_l23 = c(306, 308, 462)
scpn = c(382, 442)
icpn_l56 = c(351, 352)
ipc = c(171)



changes_list = list(
    c(which((mc_cc <= 80) & (mcmd$st == 'IPC'))),
    c(109, 447, 448, 456),
    c(175, 181, 343),
    c(306, 308, 462),
    c(382, 442),
    c(351, 352),
    c(171)
)

names(changes_list) = c('IPC_cyc', 'iCfuPN', 'iCPN_L2-3', 'CPN_L2-3', 'SCPN', 'iCPN_L5-6', 'IPC')

changes_list[[1]] = as.vector(changes_list[[1]])

lapply(changes_list, paste, collapse = ', ')

change_st_col = function(inds, type, color_key, mcmd) {
    mcmd$st[inds] = type
    mcmd$color[inds] = color_key$color[color_key$st == type]
#     names(mc_new@colors[ipc_cc]) = 'IPC_cyc'
#     mc_new@colors[ipc_cc] = gplots::col2hex('olivedrab1')
    return(mcmd)
}

# rapply(seq_along(changes_list), function(x, n, i) {change_st_col(x[[i]], n[[i]], color_key, mcmd)}, 
#        x=changes_list, n=names(changes_list))

for (i in seq_along(changes_list)) {
    mcmd = change_st_col(changes_list[[i]], names(changes_list)[[i]], color_key, mcmd)
}

for (i in seq_along(changes_list)) {
    print(mcmd[changes_list[[i]],])
}

df_mgraph = data.frame(mgraph)

# mcs_check = c(447, 448, 456)
mcs_check = 231

df_m_filt = filter(df_mgraph, mc1 %in% mcs_check | mc2 %in% mcs_check )

df_m_filt

mcmd[unique(c(df_m_filt$mc1, df_m_filt$mc2)),] %>% arrange(mc)

mc = scdb_mc('cort')
mc_new = mc
mc_new@colors = mcmd$color
scdb_add_mc('cort_new', mc_new)

mcell_mc2d_plot(mc2d_id = 'cort', colors = mc_new@colors)

readr::write_tsv(mcmd, './BonevCollab/mcmd_cort_new.tsv')

mat = scdb_mat('cort')

mat@cell_metadata$st_cort = mcmd$st[mc@mc[rownames(mat@cell_metadata)]]

scdb_add_mat('cort_new', mat)

cell_cyc = c('Mki67', 'Pcna', 'Hist1h', 'Smc4', 'Mcm3', 'Top2a')
lapply(cell_cyc, function(x) {y = grep(x, rownames(mc@mc_fp), v=T); lapply(y, function(i) mc_subset_plot_gene(mc, mc2d, mcmd, 'IPC', i))})

# options(repr.plot.height = 8, repr.plot.width = 8)

plot_diff_genes_type = function(type, mc, mc_md, mc2d){
    mcfp = mc@mc_fp[names(marks@gene_set), select(filter(mc_md, st == type), 'mc') %>% unlist]

    diff = apply(mcfp,1, max) - apply(mcfp,1, min)
    gvar = apply(mcfp,1, var)
    mcfp = mcfp[order(gvar, decreasing = T),]
    lapply(rownames(head(mcfp, 20)), function(x) mc_subset_plot_gene(mc, mc2d, mc_md, type, x))
}

dir.create('./figs/st_marks_cort')

mcmd_new = vroom::vroom('./BonevCollab//mcmd_cort_070221.tsv')
mcmd_old = vroom::vroom('./BonevCollab//mcmd_cort_new.tsv')
mc = scdb_mc('cort_new')
mc2d = scdb_mc2d('cort')
mc@colors = mcmd_new$color
scdb_add_mc('cort_new_2', mc)
mc2d@mc_id = 'cort_new_2'
scdb_add_mc2d('cort_new_2', mc2d)
mg = scdb_mgraph('cort_new')
mg@mc_id = 'cort_new_2'
scdb_add_mgraph('cort_new_2', mg)
mcell_mc2d_plot('cort_new_2')



mat_all = scdb_mat('all_rec_bon_1')
mc_all = scdb_mc('all_rec_bon_1')

dim(mc_all@e_gc)

mcmd = vroom::vroom('./BonevCollab//mc_metadata_new.tsv')

cr = mcmd$mc[mcmd$Bonev_annotation == 'CR']
cr

nm = 'cort5'

mat = scdb_mat(nm)

mc = scdb_mc(nm)

mcmd = vroom::vroom('./BonevCollab//mcmd_cort_5_iter1.tsv')

ipc_cyc = which(mcmd$st == 'IPC_cyc')
ipc = which(mcmd$st == 'IPC')
nsc = which(mcmd$st == 'NSC')
ipc_cyc
ipc
nsc

egc_cyc = setNames(apply(mc@e_gc[,ipc_cyc], 1, mean), rownames(mc@e_gc))
egc_ipc = setNames(apply(mc@e_gc[,ipc], 1, mean), rownames(mc@e_gc))
egc_nsc = setNames(apply(mc@e_gc[,nsc], 1, mean), rownames(mc@e_gc))

diff1 = egc_ipc - egc_cyc
diff2 = egc_nsc - egc_cyc

diffn1 = (egc_ipc - egc_cyc)/(egc_cyc + 1e-6)
diffn2 = (egc_nsc - egc_cyc)/(egc_cyc + 1e-6)

diffn12 = (egc_cyc - egc_ipc)/(egc_ipc + 1e-6)
diffn22 = (egc_cyc - egc_nsc)/(egc_nsc + 1e-6)

diffnn = log2((egc_cyc - apply(cbind(egc_nsc, egc_ipc), 1, max))/(egc_nsc + 1e-6))

ratio1 = log2((egc_ipc + 1e-06)/(egc_cyc + 1e-06))
ratio2 = log2((egc_nsc + 1e-06)/(egc_cyc + 1e-06))

fold_df = data.frame(t(data.frame(lapply(1:5, function(x) c(length(ratio1[ratio1 > x]),length(ratio2[ratio2 > x]))))))
rownames(fold_df) = 2**(1:5)
fold_df = tibble::rownames_to_column(fold_df)
colnames(fold_df) = c('fold', '# IPC>IPC_cyc', '# NSC>IPC_cyc')
fold_df

egc_cyc['Hpca']
egc_ipc['Hpca']
egc_nsc['Hpca']

head(ratio1[order(ratio1, decreasing = T)])
head(ratio2[order(ratio2, decreasing = T)])

head(diffnn[order(diffnn, decreasing = T)])

log2(tail(ratio1[order(ratio1)], 10))
log2(tail(ratio2[order(ratio2)], 10))

log2(tail(diffn1[order(diffn1)], 10))
log2(tail(diffn2[order(diffn2)], 10))

length(which(log2(diffn1[order(diffn1)]) > 2))
length(which(log2(diffn2[order(diffn2)]) > 3))

log2(head(diffn12[order(diffn12, decreasing = T)], 10))
log2(head(diffn22[order(diffn22, decreasing = T)], 10))

ranks_12 = setNames(match(names(diffn12), names(diffn12[order(diffn12, decreasing = T)])), names(diffn12))
ranks_22 = setNames(match(names(diffn22), names(diffn22[order(diffn22, decreasing = T)])), names(diffn22))

head(match(names(diffn12), names(diffn12[order(diffn12)])))

match(names(diffn12), names(diffn12[order(diffn12, decreasing = T)]))[which(names(diffn12) == 'Top2a')]

dim(cbind(ranks_12, ranks_22))

head(ranks_12)
head(ranks_22)

ranks_both = apply(cbind(ranks_12, ranks_22), 1, mean)

head(ranks_both)

head(ranks_both[order(ranks_both)], 20)

mcell_mc_plot_gg(mc_id = nm, g1 = 'Mt3', g2 = 'Shisa9')

sapply(c(nnsc,bad_mcs), function(x) hist(Matrix::colSums(mat@mat[,names(mc@mc)[mc@mc == x]]), 50, main = x))

cs_all = Matrix::colSums(mat@mat)
cs_bad = Matrix::colSums(mat@mat[,names(mc@mc)[mc@mc %in% bad_mcs]])
cs_nn = Matrix::colSums(mat@mat[,names(mc@mc)[mc@mc %in% nnsc]])

apply()

qa = quantile(cs_all, probs = seq(0.05, 1, 0.05))
qb = quantile(cs_bad, probs = seq(0.05, 1, 0.05))
qn = quantile(cs_nn, probs = seq(0.05, 1, 0.05))
qa
qb
qn

png('./figs/hist_all_sc.png')
h1 = hist(cs_all, 100, main = 'all sc', xlab = 'UMIs/cell')
dev.off()
png('./figs/hist_bad_mc.png')
h2 = hist(cs_bad, 50, main = paste('mcs', paste0(bad_mcs, collapse = ', ')), xlab = 'UMIs/cell')
dev.off()
png('./figs/hist_nnsc.png')
h3 = hist(cs_nn, 50, main = 'nNSC?', xlab = 'UMIs/cell')
dev.off()

data.frame(cbind(qa, qb, qn))

bad_mcs

png('./figs/bad_mc_cdf.png')
plot(ecdf(cs_all), col = 'red', ylab = 'cdf', xlab = 'UMIs/cell', main = '', cex.lab = 1.5)
lines(ecdf(cs_bad) ,col = 'blue')
lines(ecdf(cs_nn) ,col = 'green',)
legend(x = 6*1e+04, y = 0.9, legend = c('all_sc', paste('mcs', paste0(bad_mcs, collapse = ', ')), 'nNSC?'), 
      fill = c('red', 'blue', 'green'))
dev.off()
