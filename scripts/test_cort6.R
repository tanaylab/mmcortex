library(devtools)
load_all("metacell/metacell/")
load_all("metacell.flow/")
#scdb_init("scrna_db/",force_reinit = T)
scfigs_init("figs")
scdb_flow_init()



mat = scdb_mat("cort6")
mc_id = "cort6"
mgraph_id = "cort6"

mc = scdb_mc(mc_id)
cell_time = mat@cell_metadata[names(mc@mc),"day"]
cell_time = cell_time - 12
names(cell_time) = names(mc@mc)


mgraph = scdb_mgraph(mgraph_id)

net_id = "test_cortex"

time_age_groups = c(13:18)

mcell_mctnet_from_mgraph(net_id = net_id,mgraph_id = mgraph_id,cell_time = cell_time,mc_proliferation_rate_fn = "data/mmcortex_prol_rates.tsv",time_age_groups = time_age_groups)


flow_id = "test_cortex"
fig_dir = "figs/"

flow_tolerance = 0.05

message("generate flows")

mcell_new_mctnetflow(flow_id, net_id, 
                     init_mincost = T, flow_tolerance=flow_tolerance)

message("solved network flow problem")

mcf = scdb_mctnetflow(flow_id)

#compute propagatation forward and background
mcf = mctnetflow_comp_propagation(mcf)

#adding back the object with the network and flows
scdb_add_mctnetflow(flow_id, mcf)

colors_ordered = unique(mc@colors)
colors_ordered = c("#FA8072","#9FB6CD","#FF3E96","#FFA500","#8B4726","#87CEFA","#C0FF3E","#0000CD","#B4EEB4","#EED2EE","#00C5CD","#C4658D","#006400")

mm_mctnetwork_plot_net(mct_id = "test_cortex",flow_id = "test_cortex",fn = "test_cortex25.png",
  show_over_under_flow = T,
  #     show_axes=F,
  #     mc_ord = rank,
  colors_ordered = colors_ordered
  #     plot_mc_ids = T,
  #     h = 8000,
)



