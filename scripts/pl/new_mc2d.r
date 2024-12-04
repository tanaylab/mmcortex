
library(metacell)
scdb_init('./scdb')
load('./output/metacell_model/nsc_gene_modules.rda')


feats = scdb_gset('pl_filt_lat')
feats = names(feats@gene_set)
feats = feats[feats %in% rownames(mc@e_gc)]


feats_new <- feats[!(feats %in% names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3,6,7)]))]
gs_new <- gset_new_gset(sets = setNames(rep(1, length(feats_new)), feats_new), desc = 'pl_filt_lat and removed all cell cycle-correlated genes')

new_gs_id <- 'pl_filt_lat_cor_cc'
mc_id <- 'pl_cort'
new_mgraph_id <- 'pl_cort_not_cor_cc'

scdb_add_gset(id = new_gs_id, gset = gs_new)

gs <- scdb_gset(new_gs_id)
mc <- scdb_mc(mc_id)

feats_cc <- feats[feats %in% names(ct_hc_cor_nsc[ct_hc_cor_nsc %in% c(1,3,6,7)])]
cc_gs <- gset_new_gset(sets = setNames(rep(1, length(feats_cc)), feats_cc), desc = 'cell-cycle-correlated genes in NSC')

scdb_add_gset(id = 'pl_cort_nsc_cor_cc', gset = cc_gs)

mcell_add_cgraph_from_mat_bknn(mat_id='pl_cort',
                gset_id = 'pl_cort_nsc_cor_cc',
                graph_id='pl_cort_cor_cc',
                K=20,
                dsamp=T)


mcell_mgraph_logistic(mgraph_id=new_mgraph_id, mc_id=mc_id,feats_gset=new_gs_id)

mc = scdb_mc(mc_id)
mg = scdb_mgraph(new_mgraph_id)
mgraph = mg@mgraph

feat_genes = scdb_gset(new_gs_id)
feat_genes = names(feat_genes@gene_set)

uconf <- umap::umap.defaults

uconf$n_neighbors=7
uconf$min_dist = 0.5
uconf$spread = 1
symmetrize = F
umap_mgraph = F

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)

cgraph_id <- 'pl_cort'

xy = mc2d_comp_cell_coord(mc_id = mc_id,graph_id = cgraph_id, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_init('./scdb/', force_reinit = T)

mc2d_id <- new_mgraph_id

mc2d <- tgMC2D(mc_id, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph)

scdb_add_mc2d(mc2d_id, mc2d)

mcell_mc2d_plot(mc2d_id, colors = mcmd$color)

set_param(param = 'mc_plot_device', value = 'pdf', package = 'metacell')
mcell_mc2d_plot(mc2d_id, colors = mcmd$color, fig_fn = './output/metacell_model/figs/pl_cort_not_cor_cc_mc2d.pdf')

md_clvls <- clrmp[1+round(999*(mcmd$mean_day-13)/(18-13))]
set_param(param = 'mc_plot_device', value = 'pdf', package = 'metacell')
mcell_mc2d_plot(mc2d_id = 'pl_cort_not_cor_cc', colors = md_clvls, fig_fn = './output/metacell_model/figs/pl_cort_mc2d_col_by_mean_day.pdf')
