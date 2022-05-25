
# library(devtools)
library(lpsymphony)
# load_all("/net/mraid14/export/tgdata/users/atanay/proj/embdyn/metacell/", export_all=T)
library(metacell)
library(tidyverse)
library("Matrix")
scdb_init("scdb/", force_reinit=T)
# tgconfig::override_params("./config/mmc.yaml","metacell")

get_mc_cc = function() {
  ## function from Markus
  mat_id = 'cort'
  mc_id = 'cort'

  tag = 'cort'
  m_0 = 0.01
  s_0 = 0.005
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
  return(mc_cc)
}

add_leak_to_md = function(){
  mc_md_new = vroom::vroom('./BonevCollab//mcmd_3.tsv')

  mc_md_new_filt = mc_md_new[mc_md_new$Comment == 'cortical',]

  mc_cc = get_mc_cc()
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))

  md_cut = cut(mc_md_new_filt$mean_day, breaks = 13:18)
  print(levels(md_cut))
  print(table(md_cut))
  prol_rate = 2
  mult_by = 2 ** (prol_rate * (100 - mc_cc$cc_score)/100) 
  frac_old = tapply(rep(1,nrow(mc_md_new_filt)), md_cut, function(x) x/sum(x)) %>% unsplit(md_cut)
  frac_new = tapply(mult_by, md_cut, function(x) x/sum(x), simplify = F) %>% unsplit(md_cut)
  print(head(frac_old))
  print(head(frac_new))
  mc_cc$leak = pmax(0, (frac_old - frac_new)/frac_old)
  mc_md_new_filt = left_join(mc_md_new_filt, mc_cc, by = 'mc')

  write_tsv(x = mc_md_new_filt, file = './BonevCollab/mcmd_3_leak.tsv', quote_escape = F)
}

build_singemb_net = function(mat_id,mc_id,mgraph_id,net_id,fig_dir,
                             age_field = "day",
                             mc_leak = NULL,
                             capacity_var_factor = NULL, 
                             t_exp = 1,T_cost = 1e+5,
                             flow_tolerance = 0.01,
                             network_color_ord = NULL,
                             mc_ord = NULL,
                             off_capacity_cost1 = 1,
                             off_capacity_cost2 = 1000,
                             k_norm_ext_cost = 1,
                             k_ext_norm_cost = 1,
                             k_ext_ext_cost = 1) {
  options(error = utils::recover)
  mat = scdb_mat(mat_id)
  mc  = scdb_mc(mc_id)
  mgraph = scdb_mgraph(mgraph_id)
  md = mat@cell_metadata
  cell_time = mat@cell_metadata$day %>% as.vector
  names(cell_time) = rownames(mat@cell_metadata)

  if(is.null(mc_leak)) {
    leak = rep(0,max(mc@mc))
  }
  if(is.null(capacity_var_factor)) {
    capacity_var_factor = rep(0.25,max(mc@mc))
  }
  
  if(is.null(mc_leak)) {
    leak = rep(0, max(mc@mc))
  } else {
    leak = mc_leak
  }

  
  mcell_new_mctnetwork(net_id = net_id,
                       mc_id = mc_id,
                       mgraph_id = mgraph_id,
                       cell_time = cell_time)
  mct = scdb_mctnetwork(net_id)

  #computing manifold costs (based on mgraph distances)
  mct = mctnetwork_comp_manifold_costs(mct,t_exp=t_exp, T_cost=T_cost)
  message("computed manifold costs")
  
  #generating network structure	
  mct = mctnetwork_gen_network(mct, mc_leak = leak,capacity_var_factor = capacity_var_factor,
                               k_norm_ext_cost = k_norm_ext_cost,k_ext_norm_cost = k_ext_norm_cost,k_ext_ext_cost = k_ext_ext_cost,
                               off_capacity_cost1 = off_capacity_cost1,off_capacity_cost2 = off_capacity_cost2)	
  message("generated network")
  
  #solving the flow problem
  mct = mctnetwork_gen_mincost_flows(mct, flow_tolerance = flow_tolerance)
  message("solved network flow problem")
  
  #compute propagatation forward and background
  mct = mctnetwork_comp_propagation(mct)
  
  #adding back the object with the network and flows
  scdb_add_mctnetwork(paste0(net_id, '_comp'), mct)
  
  # mct = scdb_mctnetwork(paste0(net_id, '_comp'))
  
  source('~/raid/proj/mmcortex/scripts/mctnetwork_plot_net.r')

  #to plot the "big" network diagram
  # if(is.null(network_color_ord)) {
  #   mctnetwork_plot_net(mct_id = net_id, 
  #                       fn=sprintf("%s/%s_net.png",fig_dir,net_id), 
  #                       colors_ordered=network_color_ord,h = 4000,w = 2000)
  # } else {
  #   if(!is.null(mc_ord)) {
  #     mctnetwork_plot_net(mct_id = net_id, 
  #                         fn=sprintf("%s/%s_net.png",fig_dir,net_id), 
  #                         mc_ord = mc_ord,h = 4000,w = 2000)
  #   } else {
  mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mcmd_3_leak.tsv')
# mc_md = vroom::vroom('~/raid/proj/mmcortex/BonevCollab//mc_metadata.tsv')

  # ord_vec = mc_md %>% mutate(nrow = 1:nrow(mc_md)) %>% arrange(color_bon) %>% select(nrow) %>% unlist
  # ord_vec = mc_md %>% mutate(nrow = 1:nrow(mc_md)) %>% select(nrow) %>% unlist
  color_key = unique(mc_md[,c('Bonev_annotation', 'color3')]) 
  color_key = color_key %>% mutate(i = 1:nrow(color_key)) %>% tibble::column_to_rownames('Bonev_annotation')
  color_key = color_key[c('Astrocytes','early_NSC?','OPC','NSC','nNSC?','early_nNSC?','IPC','iCfuPN','iCPN_L2-3',
                  'iCPN_L5-6','CPN_L2-3','CPN_L5-6','SCPN','CthPN','Stellate_L4'),]

  mctnetwork_plot_net(mct_id = paste0(net_id, '_comp'), colors_ordered = color_key$color3,
                      fn=sprintf("%s/%s_net.png",fig_dir,net_id), 
                      # mc_ord = ord_vec,
                      h = 2500,w = 2000)
    # }

    
  # }
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  




  message("plotted the network")
}

main_build_network = function() {
  options(error = utils::recover)
  age_field = "day"
  mat_id = "cort"
  mc_id = "cort"
  mgraph_id = "cort"
  net_id = "cort"
  fig_dir = "./figs/cort.net"

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  mc = scdb_mc(mc_id)
  capacity_var_factor = rep(0.25,ncol(mc@e_gc))
  
  # next define the mc_leak parameter
  
  add_leak_to_md()

  mc_md = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_3_leak.tsv')
  mc_leak = mc_md$leak
  # mc_leak = rep(0, nrow(mc_md))
  # mc_ordd = mc_md$mc



  # mc_leak = rep(0,ncol(mc@e_gc))
  # f_extra = mc@colors == "#F6BFCB" | mc@colors == "#7F6874"
  # f_endo = mc@colors == "#0F4A9C" | mc@colors == "#EF5A9D" | mc@colors == "#F397C0" | mc@colors == "#c19f70"
  # mc_leak[f_extra] = 0.17
  # mc_leak[f_endo] = 0.12
  
  
  # you can also set network_color_ord or mc_ord for the visualisation of the network
  # mc_md = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mc_metadata.tsv')
  # mc_ordd = order(mc_md$mean_day)
  # mc_ordd = 1:nrow(mc_md)
  build_singemb_net(mat_id = mat_id,
                    mc_id = mc_id,
                    mgraph_id = mgraph_id,
                    net_id = net_id,
                    fig_dir = fig_dir,
                    age_field = age_field,
                    mc_leak = mc_leak,
                    # mc_ord = mc_ord,
                    capacity_var_factor = capacity_var_factor,
                    k_norm_ext_cost = 2,
                    k_ext_norm_cost = 2,
                    k_ext_ext_cost = 100,
                    flow_tolerance = 0.01,
                    )
  
  
}



main_build_network()