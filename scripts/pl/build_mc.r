library(metacell)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'pl'

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

mat = scdb_mat(nm)
feats = scdb_gset(paste0(nm, '_f'))

mcell_add_cgraph_from_mat_bknn(mat_id=nm,
                gset_id = paste0(nm),
                graph_id=paste0(nm),
                K=100,
                dsamp=T)

mcell_coclust_from_graph_resamp(
                coc_id=paste0(nm),
                graph_id=nm,
                min_mc_size=10,
                p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(mc_id=nm, coc_id=nm, mat_id=nm, K=50, min_mc_size=15)
                