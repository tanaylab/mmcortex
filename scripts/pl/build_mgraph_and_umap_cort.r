library(metacell)
library(umap)
uconf = umap.defaults
umap_mgraph = F

wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm_new = 'pl_cort'
# nm_new = paste0(nm, '_cort')

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

tgconfig::set_param('mcell_mgraph_max_confu_deg', 7, 'metacell')

mcell_mgraph_logistic(mgraph_id=nm_new, mc_id=nm_new,feats_gset=nm_new)

mc = scdb_mc(nm_new)
mg = scdb_mgraph(nm_new)
mgraph = mg@mgraph

feat_genes = scdb_gset(nm_new)
feat_genes = names(feat_genes@gene_set)

uconf$n_neighbors=10
uconf$min_dist = 0.5
uconf$spread = 1
symmetrize = F
umap_mgraph = F

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)
xy = mc2d_comp_cell_coord(mc_id = nm_new,graph_id = nm_new, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_init('./scdb/', force_reinit = T)
scdb_add_mc2d(nm_new, tgMC2D(nm_new, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph))
mcell_mc2d_plot(nm_new)
