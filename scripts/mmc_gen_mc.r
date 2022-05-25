# devtools::load_all("/net/mraid14/export/tgdata/users/atanay/proj/embdyn/metacell")
library('metacell')
scdb_init("scdb")

scfigs_init("figs/")
tgconfig::override_params("config/mmc.yaml",package="metacell")

source("scripts/generic_mc.r")

if(1) {
bad_genes = c(as.character(read.table("data/filtered_genes_mmc.txt", sep="\t")$x))
generate_mc("all_cells", color_key="config/base_color.txt",
			recompute=T, 
			Knn = 80, Knn_core=15,
#			Knn = 100, Knn_core=18,
			min_mc_sz=15,
			add_bad_genes = bad_genes)
}
