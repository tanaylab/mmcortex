library(metacell)
library(glue)
library(tidyverse)
library(gplots)
library(Matrix)


wd = '/home/feshap/raid/proj/mmcortex'
nm = 'merge'
nm_f = paste0(nm, '_f')
set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")

# mat = scdb_mat(nm)
# tot_umis = Matrix::colSums(mat@mat)
# mito_umis = Matrix::colSums(mat@mat[grep('^mt-', rownames(mat@mat), value = T),])
# ery_genes = grep('Hba|Hbb', rownames(mat@mat), v=T)
# ery_umi = apply(mat@mat[ery_genes,], 2, max)
# cells_ignore = colnames(mat@mat)[(mito_umis/tot_umis > 0.15) | (ery_umi > 10)]
# mcell_mat_ignore_cells(nm_f, nm, cells_ignore)

# nms = rownames(mat@mat)
# ig_genes = c(grep("^IGJ", nms, ignore.case = T, v=T), grep("^IGH",nms, ignore.case = T, v=T), 
#                  grep("^IGK", nms, ignore.case = T,  v=T), grep("^IGL", nms,  ignore.case = T, v=T))
# bad_genes = unique(c(grep('mt-', nms, value = T), grep("^MT-", nms, v=T), 
#                      grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),'Malat1', "Neat1","Tmsb4nm", "Tmsb10", ig_genes))

# mcell_mat_ignore_genes(nm_f, nm_f, ig_genes = bad_genes)
# mcell_mat_ignore_small_cells(nm_f, nm_f, 4000)

merge = scdb_mat(paste0(nm_f))
scdb_add_mat('all', merge)

# nm = 'all'

# mcell_add_gene_stat(nm, nm, force=T)

# mcell_gset_filter_varmean(gset_id=glue("{nm}_feats"), gstat_id=nm, T_vm=0.08, force_new=T)
# mcell_gset_filter_cov(gset_id=glue("{nm}_feats"), gstat_id=nm, T_tot=100, T_top3=2)
# mcell_plot_gstats(gstat_id=nm, gset_id=glue("{nm}_feats"))

# mcell_gset_split_by_dsmat(gset_id=glue("{nm}_feats"), mat_id=nm, K = 96, force = T)

# mcell_plot_gset_cor_mats(gset_id=glue("{nm}_feats"), scmat_id=nm)

# ifn1_genes = c('Isg15', 'Wars', 'Ifit1')
# cell_cyc = c('Mki67', 'Hist1h', 'Pcna', 'Smc4', 'Mcm3', 'Top2a')
# stress = c('Fos', 'Hsp90ab1', 'Hspa1a', 'Hif1a')
# misc = c('Xist', 'Tsix')
# star_genes = c(ifn1_genes, cell_cyc, stress, misc)

# feats = scdb_gset(glue("{nm}_feats"))

# genes_in = map(star_genes, function(x) grep(x, names(feats@gene_set), ignore.case = F, v = T)) %>% unlist

# uu = unique(feats@gene_set[names(feats@gene_set) %in% genes_in])
# sort(uu)

# clusts_output = '9 14 16 19 36 41 56 69 78 94'

# clusts_to_remove = c(14, 16, 36, 41, 69, 94)

# gset_nm = glue('{nm}_feats')
# mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = clusts_to_remove, 
#                          new_id = glue('{nm}_lateral'), reverse=T)
# mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = clusts_to_remove, 
#                          new_id = gset_nm, reverse=F)

# mcell_add_cgraph_from_mat_bknn(mat_id=nm,
#                 gset_id = glue("{nm}_feats"),
#                 graph_id=nm,
#                 K=100,
#                 dsamp=T)

# mcell_coclust_from_graph_resamp(
#                 coc_id=glue("{nm}_coc500"),
#                 graph_id=nm,
#                 min_mc_size=15,
#                 p_resamp=0.75, n_resamp=500)

# mcell_mc_from_coclust_balanced(
#                 coc_id=glue("{nm}_coc500"),
#                 mat_id=nm,
#                 mc_id=nm,
#                 K=30, min_mc_size=20, alpha=2)

# mcell_mc_split_filt(new_mc_id=paste0(nm, '_f'),
#             mc_id=nm,
#             mat_id=nm,
#             T_lfc=3, plot_mats=F)

# mcell_gset_from_mc_markers(gset_id=glue("{nm}_markers"), mc_id=nm)

# mcell_mc_plot_marks(mc_id=paste0(nm, '_f'), gset_id=paste0(nm, '_markers'), 
#                     mat_id=nm, lateral_gset=paste0(nm, '_lateral'))

# mcell_mc2d_force_knn(mc2d_id='test_all',mc_id=paste0(nm, '_f'), graph_id=nm)


# feats = scdb_gset(glue('{nm}_feats'))

# mc = scdb_mc(nm)

# mmc_feats_mat = mc@e_gc[names(feats@gene_set),]
# mmc_feats_mat = mmc_feats_mat[order(rownames(mmc_feats_mat)),]

# mba_clust_mm = Matrix::readMM('./data/MBA_clust_data.mtx')

# mba_genes = read.delim('./data/MBA_clust_genes.tsv', header = FALSE, quote = "") %>% 
#                 apply(1, function(x) gsub("'", "", x))

# mba_col_att = read_tsv('./data/MBA_clust_col_att.tsv')
# mba_col_att = column_to_rownames(mba_col_att, var = 'X1')
# mba_cell_types = mba_col_att['Class',]
# mba_cell_subtypes = mba_col_att['Subclass',]
# mba_markers = mba_col_att['MarkerGenes',]

# all_mba_markers = mba_markers %>% unlist %>% str_split(pattern = ' ') %>% unlist %>% unique
# markers_filt = scdb_gset('all_markers')
# markers_filt = names(markers_filt@gene_set)
# mm_markers_in_mba_inds = map(markers_filt, function(x) grep(x, mba_markers))
# mm_cell_types = map(mm_markers_in_mba_inds, function(x) unlist(mba_cell_types[x]))
# mm_cell_subtypes = map(mm_markers_in_mba_inds, function(x) unlist(mba_cell_subtypes[x]))
# mba_norm = mba_clust_mm/Matrix::colSums(mba_clust_mm)
# mba_genes_in_feats = mba_genes %in% names(feats@gene_set)
# dups_to_remove = which(mba_genes %in% mba_genes[mba_genes_in_feats][duplicated(mba_genes[mba_genes_in_feats])])[c(2,4)]
# mba_genes_in_feats[dups_to_remove] = FALSE
# mba_feats_mat = mba_norm[mba_genes_in_feats,]
# rownames(mba_feats_mat) = mba_genes[mba_genes_in_feats]
# mba_feats_mat = mba_feats_mat[order(rownames(mba_feats_mat)),]
# mats_bind = cbind(as.matrix(mba_feats_mat), mmc_feats_mat)
# mats_cor_sp = tgs_cor(as.matrix(mba_feats_mat), mmc_feats_mat, spearman = T)

# max_clust = apply(mats_cor_sp, 2, which.max)
# max_subtypes = mba_cell_subtypes[max_clust]
# max_types = mba_cell_types[max_clust]

# color_mat = expand.grid(1:26,1:26)
# color_mat = color_mat[1:length(colors()),]
# color_mat['color'] = colors()
# colnames(color_mat) = c('x', 'y', 'color')

# base_color_df = data.frame(rbind(
#         c('Neuron', filter(color_mat, ((x == 1) & (y == 20))) %>% select(color)),
#         c('Neuroblast', filter(color_mat, ((x == 13) & (y == 6))) %>% select(color)),
#         c('Radial glia', filter(color_mat, ((x == 16) & (y == 21))) %>% select(color)),
#         c('Oligo', filter(color_mat, ((x == 13) & (y == 17))) %>% select(color)),
#         c('Immune', filter(color_mat, ((x == 7) & (y == 16))) %>% select(color)),
#         c('Glioblast', filter(color_mat, ((x == 24) & (y == 23))) %>% select(color)),
#         c('Blood', filter(color_mat, ((x == 5) & (y == 6))) %>% select(color)),
#         c('Vascular', filter(color_mat, ((x == 1) & (y == 3))) %>% select(color))
#         )
# )
# base_color_df[,'rgb'] = col2hex(base_color_df[,'color'])
# colnames(base_color_df)[[1]] = 'type'

# neuro_to_blast_colors = colorRampPalette(colors = c('olivedrab2', 'magenta', 'gold1'), bias = 0.5)(22)
# subtype_color_df = data.frame(rbind(
#         cbind(c(c('Early NIPCs'),
#         c('Cycling Forebrain Inhibitory Neuroblast'),
#         c('Upper layer cortex neuron'),
#         c('Forebrain Neurons'),        
#         c('MSN Neuro'),    
#         c('Forebrain Inhibitory Neuro'),
#         c('Cycling Forebrain Neuroblast'),
#         c('Radialglia-like Progenitor'),    
#         c('NIPCs'),
#         c('Early Telencephalic Excitatory Neurons'),
#         c('Forebrain Neuroblasts'),
#         c('Cortex or Hyppocampal Neuroblast'),
#         c('CGE Neuro'),
#         c('MGE Neuroblast'),
#         c('MSNs'),
#         c('Hypothalamic Peptidergic'),
#         c('Thalamic Inhibitory Neuroblasts'),
#         c('Ventral Telencephalic Neuroblasts'),
#         c('Ventral Thelencephalic NIPCS'),
        
#         c('Cortex or Hippocampal Neuroblast'),    
#         c('Cortex or Hyppocampal Neuroblast'),
#         c('NIPCs Forebrain')    ), neuro_to_blast_colors),  

#         c('Middle dorsal forebrain radial glia', filter(color_mat, ((x == 16) & (y == 21))) %>% select(color)),
#         c('Cycling Forebrain Radialglia', filter(color_mat, ((x == 19) & (y == 21))) %>% select(color)),
#         c('Middle forebrain radial glia', filter(color_mat, ((x == 19) & (y == 21))) %>% select(color)),
#         c('Heterogeneus Radialglia', filter(color_mat, ((x == 20) & (y == 21))) %>% select(color)),
        
#         c('Astrocyte', filter(color_mat, ((x == 1) & (y == 6))) %>% select(color)),
        
#         c('PreOPC', filter(color_mat, ((x == 13) & (y == 17))) %>% select(color)),
#         c('OPCs', filter(color_mat, ((x == 12) & (y == 16))) %>% select(color)),
    
#         c('Microglia', filter(color_mat, ((x == 7) & (y == 16))) %>% select(color)),
        
#         c('AstroEpendymal', filter(color_mat, ((x == 24) & (y == 23))) %>% select(color)),
        
#         c('Erythrocyte 1', filter(color_mat, ((x == 4) & (y == 6))) %>% select(color)),
#         c('Erythrocyte 2', filter(color_mat, ((x == 6) & (y == 6))) %>% select(color)),
        
#         c('Cycling Pericyte', filter(color_mat, ((x == 3) & (y == 3))) %>% select(color)),
#         c('Vascular Endothelial', filter(color_mat, ((x == 1) & (y == 3))) %>% select(color))
#         )
# )

# subtype_color_df = cbind(subtype_color_df, as.character(matrix(ifelse(grepl('#', subtype_color_df[,'neuro_to_blast_colors']), 
#                                             subtype_color_df[,'neuro_to_blast_colors'], col2hex(subtype_color_df[,'neuro_to_blast_colors'])))))
# colnames(subtype_color_df) = c('subtype', 'color', 'rgb')

# types_in_data = base_color_df[,'type']
# subtypes_in_data = subtype_color_df[,'subtype']
# type_inds = match(max_types, types_in_data)
# subtype_inds = match(max_subtypes, subtypes_in_data)