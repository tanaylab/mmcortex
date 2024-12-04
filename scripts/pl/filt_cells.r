library(metacell)

wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)
nm = 'pl'

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))


mat = scdb_mat('merge')

mito_genes = grep('^mt-', rownames(mat@mat),v=T)
hb_genes = grep('^Hb[ab]', rownames(mat@mat), v=T)


mito = Matrix::colSums(mat@mat[mito_genes,])
hb = Matrix::colSums(mat@mat[hb_genes,])
tot = Matrix::colSums(mat@mat)

mito_ig = colnames(mat@mat)[mito/tot >= 0.15]

tot_ig = colnames(mat@mat)[tot <= 1000 | tot >= 30000]

hb_ig = colnames(mat@mat)[hb/tot >= 0.005]

ig_12 = rownames(mat@cell_metadata)[grep('E12', mat@cell_metadata$amp_batch_id)]

df_pred = readRDS('./data/df_pred_pl.rds')
doub_ig = df_pred$doublet_names[df_pred$doublet_preds == 'Doublet']

cells_filt = unique(c(mito_ig, tot_ig, hb_ig, doub_ig, ig_12))

mcell_mat_ignore_cells(nm, nm, cells_filt, reverse=F)
mcell_mat_ignore_genes(nm, nm, mito_genes, reverse=F)

scdb_init(db_path, force_reinit = T)
mat = scdb_mat(nm)
cmd = mat@cell_metadata
cmd$day = as.character(purrr::map(stringr::str_split(cmd$amp_batch_id, '_'), 1))
cmd$t = as.numeric(gsub('^E', '', cmd$day))

mat@cell_metadata = cmd
scdb_add_mat(nm, mat)