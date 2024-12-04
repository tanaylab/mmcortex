library(metacell)
library(umap)
uconf = umap.defaults
umap_mgraph = F

wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'pl'

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

tgconfig::set_param('mcell_mgraph_max_confu_deg', 3, 'metacell')


mc = scdb_mc(nm)
feats = scdb_gset(nm)
feats_f = feats
feats_f@gene_set = feats_f@gene_set[names(feats_f@gene_set) %in% rownames(mc@e_gc)]
scdb_add_gset(paste0(nm, '_in_mc'), feats_f)

mcell_mgraph_logistic(mgraph_id=nm, mc_id=nm,feats_gset=paste0(nm,'_f'))

mc = scdb_mc(nm)
mg = scdb_mgraph(nm)
mgraph = mg@mgraph

feat_genes = scdb_gset(paste0(nm, '_f'))
feat_genes = names(feat_genes@gene_set)

uconf$n_neighbors=10
uconf$min_dist = 0.5
uconf$spread = 1
symmetrize = F
umap_mgraph = F

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)
xy = mc2d_comp_cell_coord(mc_id = nm,graph_id = nm, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_init('./scdb/', force_reinit = T)
scdb_add_mc2d(nm, tgMC2D(nm, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph))
mcell_mc2d_plot(nm)
