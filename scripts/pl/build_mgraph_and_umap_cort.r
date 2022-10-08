library(metacell)
library(umap)
uconf = umap.defaults
umap_mgraph = F

wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
# nm_new = 'pl_cort_dns'
nm_new = 'pl_cort'
# nm_new = paste0(nm, '_cort')
mc_id <- nm_new
cgraph_id <- 'pl_cort'
mgraph_id <- nm_new
mc2d_id <- nm_new
feats_gset_id <- nm_new

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

tgconfig::set_param('mcell_mgraph_max_confu_deg', 7, 'metacell')

## Fix gene set
gs <- scdb_gset(feats_gset_id)
mc <- scdb_mc(mc_id)
gs@gene_set <- gs@gene_set[names(gs@gene_set) %in% rownames(mc@e_gc)]
scdb_add_gset(feats_gset_id, gs)

mcell_mgraph_logistic(mgraph_id=mgraph_id, mc_id=mc_id,feats_gset=feats_gset_id)

mc = scdb_mc(mc_id)
mg = scdb_mgraph(mgraph_id)
mgraph = mg@mgraph

feat_genes = scdb_gset(feats_gset_id)
feat_genes = names(feat_genes@gene_set)

uconf$n_neighbors=7
uconf$min_dist = 0.5
uconf$spread = 1
symmetrize = F
umap_mgraph = F

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)
xy = mc2d_comp_cell_coord(mc_id = mc_id,graph_id = cgraph_id, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_init('./scdb/', force_reinit = T)
scdb_add_mc2d(mc2d_id, tgMC2D(mc2d_id, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph))
mcell_mc2d_plot(mc2d_id)
