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
mc_mean_day = tapply(cell_time, mc@mc, mean) + 12

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

eomes_egc = mc@e_gc['Eomes',]
eomes_egc_lin = (eomes_egc - min(eomes_egc))/(max(eomes_egc) - min(eomes_egc))

max_div_per_day = 2
# prol_df = dplyr::mutate(mc_cc, proliferation_rate = max_div_per_day*(100-mc_cc$cc_score)/100) %>% 
#             select(1,3) %>% 
#             tibble::rownames_to_column(.)
prol_df = dplyr::mutate(mc_cc, proliferation_rate = as.numeric(max_div_per_day*(100-mc_cc$cc_score)*(1 - eomes_egc_lin)*(1 - 0.4*plogis(mc_mean_day,mean(c(13,18)),.15))/100)) %>% 
            select(1,3) %>% 
            tibble::rownames_to_column(.)

readr::write_tsv(x=prol_df, file = './data/pl_cort_prol_rates.tsv')


occ2 = 1e+04
# capacity_var_factor = 0.5

mcell_mctnet_from_mgraph(net_id = net_id,mgraph_id = mgraph_id,cell_time = cell_time,
                         mc_proliferation_rate_fn = "data/pl_cort_prol_rates.tsv",
                         time_age_groups = time_age_groups,
                         off_capacity_cost2 = occ2
                         )


flow_id = net_id
fig_dir = "figs/"

# flow_tolerance = 0.01
# max_flow_tolerance = 0.05
flow_tolerance = 0
max_flow_tolerance = 0
    
message("generate flows")

mct = scdb_mctnetwork(net_id)

# mcf = scdb_mctnetflow(net_id)

# differentiation_penalty_factor = 1e+5

# terminal_genes = c('Mapt','Mef2c','Meg3', 'Stmn2')

# mat_mod = mc@mc_fp[terminal_genes,]
# mat_z = exp(t(apply(mat_mod, 1, function(x) (x - mean(x))/sd(x)))**2)

# mc_penalty  = differentiation_penalty_factor * apply(mat_z, 2, sum)
# # mc_penalty = mc_penalty/min(mc_penalty)
# mcmfc_new = mct@mc_manif_cost
# mcmfc_new$cost = ifelse(mct@mc_manif_cost$mc1 != mct@mc_manif_cost$mc2, mct@mc_manif_cost$cost*mc_penalty[mcmfc_new$mc1], mct@mc_manif_cost$cost)
                    
# mct@mc_manif_cost = mcmfc_new
# scdb_add_mctnetwork(net_id, mct)

# mct@mc_manif_cost = mcmfc_new
# scdb_add_mctnetwork(net_id, mct)

mcell_new_mctnetflow(flow_id, net_id, 
                     init_mincost = T, flow_tolerance=flow_tolerance, max_flow_tolerance = max_flow_tolerance)

message("solved network flow problem")

mcf = scdb_mctnetflow(flow_id)

#compute propagatation forward and background
mcf = mctnetflow_comp_propagation(mcf)

#adding back the object with the network and flows
scdb_add_mctnetflow(flow_id, mcf)

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
color_key = unique(mcmd[,c('st', 'color')])
color_ord = color_key$color[match(cust_st_ord, color_key$st)]
mctnetwork_plot_net(net_id, net_id, fn = './figs/pl_cort_test_plot_net_scale_0.15_coef_0.4.png', colors_ordered=color_ord)