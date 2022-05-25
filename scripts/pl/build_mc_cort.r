library(metacell)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'pl'
nm_new = paste0(nm, '_cort')

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

mat = scdb_mat(nm_new)
feats = scdb_gset(paste0(nm, '_f'))

mcell_add_cgraph_from_mat_bknn(mat_id=nm_new,
                gset_id = paste0(nm, '_f'),
                graph_id=nm_new,
                K=100,
                dsamp=T)

mcell_coclust_from_graph_resamp(
                coc_id=nm_new,
                graph_id=nm_new,
                min_mc_size=10,
                p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(mc_id=nm_new, coc_id=nm_new, mat_id=nm_new, K=50, min_mc_size=15)

mcell_gset_from_mc_markers(gset_id = paste0(nm_new, '_marks'), mc_id = nm_new)

mcell_gset_from_mc_markers(gset_id = paste0(nm_new, '_marks_f'), mc_id = nm_new, 
        blacklist_gset_id = paste0(nm, '_lateral'))

mc = scdb_mc(nm_new)

feats = scdb_gset(paste0(nm, '_f'))

feats_f = feats
feats_f@gene_set = feats_f@gene_set[names(feats_f@gene_set) %in% rownames(mc@e_gc)]

scdb_add_gset(nm_new, feats_f)