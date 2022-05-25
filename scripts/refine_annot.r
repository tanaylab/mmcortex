library(metacell)

wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'all'
set.seed(1337)


db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")

mc_hc = mcell_mc_hclust_confu(mc_id = nm, graph_id = nm)
mc_sup = mcell_mc_hierarchy(mc_id = nm, mc_hc = mc_hc, T_gap = 0.04)

mcell_mc_plot_hierarchy(mc_id = nm, graph_id = nm,
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width = 3000, height = 10000,
                        min_nmc=2, show_mc_ids = T)


