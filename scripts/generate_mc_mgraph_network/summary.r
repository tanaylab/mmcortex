library(devtools)
library(lpsymphony)
load_all("/net/mraid14/export/tgdata/users/atanay/proj/embdyn/metacell/")
scdb_init("scdb/", force_reinit=T)
scfigs_init("figs")

# short summary of generating a network

# first generate an mgraph
# in the paper we used "mcell_mgraph_max_confu_deg" = 4 but you can set it much higher, e.g. 20 or 30
# you need to set the mgraph_id, gset_id and mat_id in the script
source("gastru_scripts/generate_mc_mgraph_network/gen_mgraph.r")

# next visualize the mgraph - we tried out our usual mcell_mc2d_force_knn layout as well as the umap package

# you should check in the 2D visualization that there are no disconnected components unless you really think that they are separated lineage-wise
source("gastru_scripts/generate_mc_mgraph_network/gen_mgraph_umap.r")

# next generate the network
# adjust all the mc_id, net_id etc.
# you need to adjust the leak parameter in the function main_build_network
# specify the age_field. This the metadata entry in the mat object that assigns to each cell an age group. I think it should be integers 1,2,3,4...
# no need to change the capacity_var_factor at this stage.

main_build_network()


