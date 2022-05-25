wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library(devtools)
# library("metacell")
library(Matrix)
load_all("~/src/metacell/")
load_all("~/src/metacell.flow/")
scdb_init("scdb/",force_reinit = T)
scfigs_init("figs")
scdb_flow_init()

SEED = 1337
set.seed(SEED)

nm = 'pl_cort'

mat = scdb_mat(nm)

mc = scdb_mc(nm)



mgraph_id = nm

cell_time = mat@cell_metadata[names(mc@mc),"t"]
cell_time = cell_time - 12
names(cell_time) = names(mc@mc)

net_id = paste0(nm, '_test')

time_age_groups = c(13:18)



get_mc_cc = function() {
  ## function from Markus
  mat_id = nm
  mc_id = nm

  tag = nm
  m_0 = 0.006
  s_0 = 0.002
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))

  return(mc_cc)
}

mc_cc = get_mc_cc()


cc_score = mc_cc$cc_score

max_div_per_day = 2
prol_df = dplyr::mutate(mc_cc, proliferation_rate = max_div_per_day*(100-mc_cc$cc_score)/100) %>% 
            select(1,3) %>% 
            tibble::rownames_to_column(.)
# prol_df

readr::write_tsv(prol_df, file = './data/pl_cort_prol_rates.tsv')


# occ2 = 1e+02
# capacity_var_factor = 0.5

# occ2 = 10**seq(1,7,1)
occ2 = 10**seq(8,10,1)
flow_tol = seq(0,0.02,0.01)
cap_var = c(0.05, seq(0.1,0.5,0.1))
param_gs = expand.grid(occ2, flow_tol, cap_var)
colnames(param_gs) = c('occ2', 'flow_tol', 'cap_var')
# print(nrow(param_gs))
# print(round(0.02*nrow(param_gs)))
param_gs = dplyr::sample_n(param_gs, round(0.5*nrow(param_gs)))

flow_id = net_id
fig_dir = "figs/"
flow_fig_dir = './figs/flow_grid_search'
gs_data_dir = './data/flow_grid_search'
gs_mcf_dir = './data/flow_grid_search_mcf'
if (!dir.exists(flow_fig_dir)) {dir.create(flow_fig_dir)}
if (!dir.exists(gs_data_dir)) {dir.create(gs_data_dir)}
if (!dir.exists(gs_mcf_dir)) {dir.create(gs_mcf_dir)}

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late', 'iCPN_early','iCPN_late',
                          'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
color_key = unique(mcmd[,c('st', 'color')])
color_ord = color_key$color[match(cust_st_ord, color_key$st)]

scpn_mc = which(mcmd$st == 'SCPN')
cthpn_mc = which(mcmd$st == 'CthPN')
cpn_l56_mc = which(mcmd$st == 'CPN_L5_6')
cpn_l23_mc = which(mcmd$st == 'CPN_L2-3')
mature_mc = list(scpn_mc, cthpn_mc, cpn_l56_mc, cpn_l23_mc)
mat_mc = unlist(mature_mc)

run_gs = function(occ2, flow_tol, cap_var) {
    scdb_init('scdb', force_reinit = T)
    mcell_mctnet_from_mgraph(net_id = net_id,mgraph_id = mgraph_id,cell_time = cell_time,
                            capacity_var_factor = cap_var,
                            mc_proliferation_rate_fn = "data/pl_cort_prol_rates.tsv",
                            time_age_groups = time_age_groups,
                             off_capacity_cost2 = occ2
                            )




    flow_tolerance = flow_tol
    max_flow_tolerance = flow_tol
                                  
    message("generate flows")

    mct = scdb_mctnetwork(net_id)

    mcell_new_mctnetflow(flow_id, net_id, 
                        init_mincost = T, flow_tolerance=flow_tolerance, max_flow_tolerance = max_flow_tolerance)

    message("solved network flow problem")
    # scdb_init('scdb', force_reinit = T)
    mcf = scdb_mctnetflow(flow_id)

    #compute propagatation forward and background
    mcf = mctnetflow_comp_propagation(mcf)

    #adding back the object with the network and flows
    scdb_add_mctnetflow(flow_id, mcf)
    mcf = scdb_mctnetflow(flow_id)

    st_outflow = lapply(mature_mc, function(x) {
          lapply(mcf@mc_forward, function(y) apply(y[x,mat_mc[!(mat_mc %in% x)]], 1, sum))
          }
          )
    save(st_outflow, file=paste0(gs_data_dir, glue::glue('/occ2={occ2}_cap_var={cap_var}_flow_tol={flow_tol}.Rda')))
    save(mcf, file=paste0(gs_mcf_dir, glue::glue('/mcf_occ2={occ2}_cap_var={cap_var}_flow_tol={flow_tol}.Rda')))
    flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0])
    mctnetwork_plot_net(net_id, net_id, flow_thresh = flow_thresh,fn = paste0(flow_fig_dir, 
          glue::glue('/occ2={occ2}_cap_var={cap_var}_flow_tol={flow_tol}.png')), colors_ordered=color_ord)

}

for (i in 1:nrow(param_gs)) {
    run_gs(occ2 = param_gs[i,'occ2'], flow_tol = param_gs[i,'flow_tol'], cap_var = param_gs[i,'cap_var'])
}