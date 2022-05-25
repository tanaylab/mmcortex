library(metacell)
library(glue)
library(tidyverse)
library(Matrix)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'all'
nm_f = paste0(nm, '_f')
set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")

mat = scdb_mat('all')

mc = scdb_mc('all')

feats = scdb_gset('all_feats_f')

markers = scdb_gset('all_markers')

lateral = scdb_gset('all_lateral')

mat@cell_metadata[,'day'] = gsub('E', '', gsub('_rep\\d_L\\d', '', mat@cell_metadata$batch_set_id)) %>% as.numeric

mmc_feats_mat = mc@e_gc[names(feats@gene_set),]
mmc_feats_mat = mmc_feats_mat[order(rownames(mmc_feats_mat)),]

mba_clust_mm = Matrix::readMM('./data/MBA_clust_data.mtx')

mba_genes = read.delim('./data/MBA_clust_genes.tsv', header = FALSE, quote = "") %>% 
                apply(1, function(x) gsub("'", "", x))

mba_col_att = read_tsv('./data/MBA_clust_col_att.tsv')
mba_col_att = column_to_rownames(mba_col_att, var = 'X1')
mba_cell_types = mba_col_att['Class',]
mba_cell_subtypes = mba_col_att['Subclass',]
mba_markers = mba_col_att['MarkerGenes',]

all_mba_markers = mba_markers %>% unlist %>% str_split(pattern = ' ') %>% unlist %>% unique
markers_filt = scdb_gset('all_markers')
markers_filt = names(markers_filt@gene_set)
mm_markers_in_mba_inds = map(markers_filt, function(x) grep(x, mba_markers))
mm_cell_types = map(mm_markers_in_mba_inds, function(x) unlist(mba_cell_types[x]))
mm_cell_subtypes = map(mm_markers_in_mba_inds, function(x) unlist(mba_cell_subtypes[x]))
mba_norm = mba_clust_mm/Matrix::colSums(mba_clust_mm)
mba_genes_in_feats = mba_genes %in% names(feats@gene_set)
dups_to_remove = which(mba_genes %in% mba_genes[mba_genes_in_feats][duplicated(mba_genes[mba_genes_in_feats])])[c(2,4)]
mba_genes_in_feats[dups_to_remove] = FALSE
mba_feats_mat = mba_norm[mba_genes_in_feats,]
rownames(mba_feats_mat) = mba_genes[mba_genes_in_feats]
mba_feats_mat = mba_feats_mat[order(rownames(mba_feats_mat)),]
mats_bind = cbind(as.matrix(mba_feats_mat), mmc_feats_mat)
mats_cor_sp = tgs_cor(as.matrix(mba_feats_mat), mmc_feats_mat, spearman = T)

max_clust = apply(mats_cor_sp, 2, which.max)
max_subtypes = mba_cell_subtypes[max_clust]
max_types = mba_cell_types[max_clust]

subtype_by_type = lapply(unique(unlist(max_types)), function(x) unique(unlist(max_subtypes[max_types == x])))
names(subtype_by_type) = unique(unlist(max_types))
subtype_by_type = stack(subtype_by_type) %>% rename(type = ind, subtype = values) %>% relocate(type)
col = subtype_by_type[,'type']
subtype_by_type[,'new_type'] = ifelse((col == 'Neuron') | (col == 'Neuroblast'), 'Neuron_Neuroblast', levels(col)[col]) 
subtype_by_type = relocate(subtype_by_type, new_type, .after = type)
color_key = c('Neuron_Neuroblast' = list(c('olivedrab2', 'blue3', 'red2', 'gold1')),
                'Radial glia' = list(c('pink', 'deeppink4')),
                'Oligo' = list(c('lightseagreen', 'darkgreen')),
                'Immune' = list(c('lemonchiffon', 'lemonchiffon3')),
                'Glioblast' = list(c('slateblue2', 'blue')),
                'Blood' = list(c('firebrick1', 'firebrick4')),
                'Vascular' = list(c('chocolate1', 'chocolate4'))
             )

subtype_by_type[,'subtype_rgb'] = map(unique(subtype_by_type[, 'new_type']), function(x) colorRampPalette(colors = 
                                            color_key[[x]])(dim(filter(subtype_by_type, new_type ==x))[[1]])) %>% unlist
subtype_by_type[,'type_rgb'] = map(unique(subtype_by_type[, 'new_type']), function(x) colorRampPalette(colors = 
                                            color_key[[x]][[1]])(dim(filter(subtype_by_type, new_type ==x))[[1]])) %>% unlist
subtype_by_type

st_color_vec = subtype_by_type[match(as.matrix(max_subtypes), subtype_by_type[,'subtype']), 'subtype_rgb']

t_color_vec = subtype_by_type[match(as.matrix(max_types), subtype_by_type[,'type']), 'type_rgb']

uu = sort(unique(mc@mc))




mc_mean_day = lapply(uu, function(u) mean(mat@cell_metadata[names(mc@mc)[mc@mc == u],'day'])) 
pal = colorRampPalette(colors = c('purple2', 'orange3'))(100)
day_quantiles = 100*round(rank(as.numeric(mc_mean_day))/length(mc_mean_day), digits = 2)
batch_color_vec = pal[day_quantiles]

write(x = batch_color_vec, file = './batch_mc_color_vec.tsv')
write_tsv(x = subtype_by_type, file = './st_color_df.tsv')
write(x = st_color_vec, file = './st_mc_color_vec.tsv')

cols = read_tsv('./st_mc_color_vec.tsv') 
cols = c(colnames(cols), as.matrix(cols))

mc@colors = cols

scdb_add_mc('all', mc)

feat_exp = rownames_to_column(data.frame(mc@e_gc[rownames(mc@e_gc) %in% names(feats@gene_set),]))
feat_exp = feat_exp[order(feat_exp[,'rowname']),]

write_excel_csv(feat_exp, col_names = FALSE, file = file.path(wd, 'BonevCollab/mc_feature_gene_mat.xls'), quote_escape = FALSE)

mc_metadata = rbind(max_subtypes, mc_mean_day, st_color_vec) %>% t %>% as_tibble
colnames(mc_metadata) = c('inferred_subtype', 'mean_day', 'color')
# rownames(mc_metadata) = 1:dim(mc_metadata)[[1]]
mc_metadata = rownames_to_column(mc_metadata)
colnames(mc_metadata) = c('mc', tail(colnames(mc_metadata), -1))
mc_metadata = cbind(mc_metadata, t(mc@n_bc))
write_tsv(mc_metadata, file = file.path(wd, 'BonevCollab/mc_metadata.tsv'), quote_escape = FALSE)
write_excel_csv(mc_metadata, file = file.path(wd, 'BonevCollab/mc_metadata.xls'), quote_escape = FALSE)

mc_of_c = mc_metadata$inferred_subtype[mc@mc]
names(mc_of_c) = names(mc@mc)
mc_of_c = mc_of_c %>% data.frame
mat@cell_metadata = left_join(rownames_to_column(mat@cell_metadata), 
                                rownames_to_column(mc_of_c), 
                                by='rowname') %>% 
                                rename('st' = '.') 
mat@cell_metadata = column_to_rownames(mat@cell_metadata, 'rowname')

scdb_add_mat('all', mat)
