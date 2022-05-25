library(metacell)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'

scdb_init(file.path(wd, 'scdb'), force_reinit = T)
scfigs_init(file.path(wd, 'figs'))

atac_files = list.files(file.path(wd, 'scatac_data'))
# atac_files

agg_id = vroom::vroom(file.path(wd, 'scatac_data/aggregation_id.csv'))

print('importing mats')
promo_peaks = readRDS(file.path(wd, 'scatac_data/prom_counts.RDS'))
gb_peaks = readRDS(file.path(wd, 'scatac_data/genebody_counts.RDS'))
ig_peaks = readRDS(file.path(wd, 'scatac_data/peak_counts.RDS'))
print('done importing')

print('creating metacells mats')
make_mat = function(mat, nm) {
    rep_nums = sapply(stringr::str_split(colnames(mat), '-'),function(x) x[[2]])
    reps = agg_id$sample_name[match(rep_nums, agg_id$library_id)]
    cmd = tibble::column_to_rownames(data.frame(cbind(colnames(mat), rep_nums, reps)), 'V1')
    cmd$day = as.numeric(gsub('E', '', sapply(stringr::str_split(cmd$reps, '_'), function(x) x[[1]])))
    print(dim(mat))
    if (length(which(duplicated(rownames(mat)))) > 0) {
            mat = tgs_matrix_tapply(t(as.matrix(mat)), rownames(mat), sum)
    }
    print(dim(mat))
    mat_raw = scm_new_matrix(mat, cmd)
    scdb_add_mat(nm, mat_raw)
}

make_mat(promo_peaks, 'pl_prom_raw')
make_mat(gb_peaks, 'pl_gb_raw')
make_mat(ig_peaks, 'pl_ig_raw')

print('done creating metacells mats')

csi = Matrix::colSums(ig_peaks)
csp = Matrix::colSums(promo_peaks)
csg = Matrix::colSums(gb_peaks)
q_max = 0.98
q_min = 0.1

mat = scdb_mat('pl_prom_raw')
cmd = mat@cell_metadata

cells_ig = colnames(ig_peaks)[csi >= quantile(csi, q_min) & csi <= quantile(csi, q_max) & cmd$day != 12]
cells_gb = colnames(gb_peaks)[csg >= quantile(csg, q_min) & csg <= quantile(csg, q_max) & cmd$day != 12]
cells_prom = colnames(promo_peaks)[csp >= quantile(csp, q_min) & csp <= quantile(csp, q_max) & cmd$day != 12]

good_cells = intersect(cells_gb, intersect(cells_ig, cells_prom))

table(cmd[good_cells,'day'])

print('filtering cells')
mcell_mat_ignore_cells('pl_prom','pl_prom_raw',good_cells,T)
mcell_mat_ignore_cells('pl_gb','pl_gb_raw',good_cells,T)
mcell_mat_ignore_cells('pl_ig','pl_ig_raw',good_cells,T)
print('done filtering cells')



# gb_md = strsplit(names(rownames(gb_peaks)), '-') %>% data.frame %>% t %>% data.frame
# gb_md = gb_md[!duplicated(rownames(gb_peaks)),]
# colnames(gb_md) = c('chrom', 'start', 'end')
# gb_md[,c('start', 'end')] = apply(gb_md[,c('start', 'end')], 2, as.numeric)
# gb_md = mutate(gb_md, len = abs(end-start))
# rownames(gb_md) = rownames(gb_peaks)[!duplicated(rownames(gb_peaks))]
# gb_norm = gb_peaks[!duplicated(rownames(gb_peaks)),]/gb_md$len

