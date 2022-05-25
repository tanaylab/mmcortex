library(metacell)
library(glue)
library(tidyverse)

wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'all'
set.seed(1337)

setwd(wd)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")

mcell_add_cgraph_from_mat_bknn(mat_id=nm,
                gset_id = glue("{nm}_feats_f"),
                graph_id=nm,
                K=100,
                dsamp=T)

mcell_coclust_from_graph_resamp(
                coc_id=glue("{nm}_coc500"),
                graph_id=nm,
                min_mc_size=15,
                p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(
                coc_id=glue("{nm}_coc500"),
                mat_id=nm,
                mc_id=nm,
                K=30, min_mc_size=20, alpha=2)

mcell_mc_split_filt(new_mc_id=paste0(nm, '_f'),
            mc_id=nm,
            mat_id=nm,
            T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(gset_id=glue("{nm}_markers"), mc_id=paste0(nm, '_f'))

mcell_mc_plot_marks(mc_id=paste0(nm, '_f'), gset_id=paste0(nm, '_markers'), 
                    mat_id=nm, lateral_gset=paste0(nm, '_lateral'))

mcell_mc2d_force_knn(mc2d_id='test_all',mc_id='all', graph_id='all')

tgconfig::set_param("mcell_mgraph_max_confu_deg",10,"metacell")
mcell_mgraph_logistic(mgraph_id = 'all_k_10_log', mc_id = 'all_f', 
			feats_gset = 'all_feats_f')