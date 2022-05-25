library("devtools")
load_all("metacell")
source("scripts/generic_mc.r")
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")

scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")


bad_genes = read.table("data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)
bad_genes = bad_genes[,1]

mat_nm = "sing_emb_wt10_1_2"
mc_id = paste(mat_nm,"_bs500f",sep="")

#generate_mc(mat_nm, color_key=NA,add_bad_genes = bad_genes)
generate_mc(mat_nm, color_key=NA,add_bad_genes = bad_genes,recompute = T)

gotg_atlas = mcell_gen_atlas(mat_id = "emb_gotg", 
									mc_id = "emb_gotg_bs500f", 
									gset_id  = "emb_gotg", 
									mc2d_id= "emb_gotg_bs500f")

wt_atlas = mcell_gen_atlas(mat_id = "sing_emb_wt10", 
                           mc_id = "sing_emb_wt10_recolored2", 
                           gset_id  = "sing_emb_wt10", 
                           mc2d_id= "sing_emb_wt10_recolored2")

cmp = mcell_proj_on_atlas(mat_id = mat_nm, mc_id = mc_id, atlas = wt_atlas, fig_cmp_dir = paste("figs/wt_atlas_1_2",mat_nm,sep=""), recolor_mc_id = mc_id)
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")
mcell_mc2d_force_knn(mc2d_id = mc_id,mc_id = mc_id,graph_id = mat_nm)
mcell_mc2d_plot(mc2d_id = mc_id)

mcell_mc2d_plot_by_factor(mc_id, mat_nm, "embryo", single_plot = F)
mcell_mc2d_plot_by_factor(mc_id, mat_nm, "Sort.Date", single_plot = T)
