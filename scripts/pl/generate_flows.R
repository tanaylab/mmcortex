wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
library("metacell")
library(Matrix)
devtools::load_all("~/src/metacell.flow/")
scdb_init("scdb/",force_reinit = T)
scfigs_init("figs")
scdb_flow_init()

SEED = 1337
set.seed(SEED)

nm = 'pl_cort'
mat_id <- 'pl_cort'
mc_id <- nm
mgraph_id = nm
# net_id = paste0(nm, '_test')
net_id = nm
flow_id = net_id
fig_dir = "figs/"
occ2 = 1e+04
max_div_per_day = 2

mat = scdb_mat(mat_id)
mc = scdb_mc(mc_id)

cell_time = mat@cell_metadata[names(mc@mc),"t"]
cell_time = cell_time - 12
names(cell_time) = names(mc@mc)
mc_mean_day = tapply(cell_time, mc@mc, mean) + 12

time_age_groups = c(13:18)

get_mc_cc = function() {
  ## function from Markus
  tag = nm
  m_0 = 0.0025
  s_0 = 0.001
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
# print(head(cc_score))
save(cc_score, file = glue::glue(file.path(wd, 'data/{nm}_cc_score.rda')))
# load(glue::glue(file.path(wd, 'data/{nm}_cc_score.rda')))
# print(head(cc_score))
eomes_egc <- mc@e_gc['Eomes',]
eomes_egc_lin <- (eomes_egc - min(eomes_egc))/(max(eomes_egc) - min(eomes_egc))
eomes_factor <- (1-0.5*plogis(eomes_egc_lin, location = 0.1, scale = 0.01))
mature_genes <- c('Mef2c', 'Mapt')
mature_genes_egc <- colSums(mc@e_gc[mature_genes,])
mature_genes_egc_lin <- (mature_genes_egc - min(mature_genes_egc))/(max(mature_genes_egc) - min(mature_genes_egc))
# mature_genes_factor <- (1-0.5*plogis(mature_genes_egc_lin, location = 0.25, scale = .025))
mature_genes_factor <- 1

# prol_df = dplyr::mutate(mc_cc, proliferation_rate = max_div_per_day*(100-mc_cc$cc_score)/100) %>% 
#             select(1,3) %>% 
#             tibble::rownames_to_column(.)

prol_df = dplyr::mutate(mc_cc, proliferation_rate = (2**as.numeric(max_div_per_day*((101-mc_cc$cc_score)/100)*eomes_factor))*mature_genes_factor) %>% 
            select(1,3) %>% 
            tibble::rownames_to_column(.)
# prol_df = dplyr::mutate(mc_cc, proliferation_rate = as.numeric(max_div_per_day*(100-mc_cc$cc_score)*(1 - eomes_egc_lin)*(1 - 0.4*plogis(mc_mean_day,mean(c(13,18)),.15))/100)) %>% 
#             select(1,3) %>% 
#             tibble::rownames_to_column(.)

# readr::write_tsv(x=prol_df, file = glue::glue('./data/{nm}_prol_rates.tsv'))

capacity_var_factor <- 0.2

mcell_mctnet_from_mgraph(net_id = net_id,mgraph_id = mgraph_id,cell_time = cell_time,
                         mc_proliferation_rate_fn = glue::glue('./data/{nm}_prol_rates.tsv'),
                         time_age_groups = time_age_groups,
                         capacity_var_factor = capacity_var_factor,
                         off_capacity_cost2 = occ2
                         )


flow_tolerance = 0.04
max_flow_tolerance = 0.04
# flow_tolerance = 0
# max_flow_tolerance = 0
    
message("generate flows")

mct = scdb_mctnetwork(net_id)

mcell_new_mctnetflow(flow_id, net_id, 
                     init_mincost = T, flow_tolerance=flow_tolerance, max_flow_tolerance = max_flow_tolerance)

message("solved network flow problem")

mcf = scdb_mctnetflow(flow_id)

#compute propagatation forward and background
mcf = mctnetflow_comp_propagation(mcf)

#adding back the object with the network and flows
scdb_add_mctnetflow(flow_id, mcf)

mcmd = readr::read_tsv('~/raid/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv')
# color_key = unique(mcmd[,c('cell_type', 'color')])
cust_st_ord = c('OPCsC_cyc', 'IPC','iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
# mcmd = readr::read_tsv('./BonevCollab/mcmd_pl_cort_freeze_1_9_22.tsv')
cust_st_ord_mc = unlist(lapply(cust_st_ord, function(u) 
    setNames(which(mcmd$cell_type == u)[order(mcmd$mean_day[which(mcmd$cell_type == u)])], rep(u, length(which(mcmd$cell_type == u))))))
# 
# color_ord = color_key$color[match(cust_st_ord, color_key$cell_type)]

source('./scripts/util.r')

mctnetwork_plot_net_YSh(net_id, net_id, fn = glue::glue('./output/metacell_flow/figs/{flow_id}.png'), mc_ord = cust_st_ord_mc, flow_thresh =0)

set_param(param = 'mc_plot_device', value = 'pdf', package = 'metacell')
mctnetwork_plot_net_YSh(net_id, net_id, plot_pdf = TRUE, fn = glue::glue('./output/metacell_flow/figs/{flow_id}.pdf'), mc_ord = cust_st_ord_mc, flow_thresh =0)
