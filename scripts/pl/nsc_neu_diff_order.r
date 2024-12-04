wd = '/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)

library(metacell)
library(tgstat)
devtools::load_all("~/src/metacell.flow")

scdb_init('scdb', f = T)

source('scripts/util.r')

mc <- scdb_mc('pl_cort')
mct <- scdb_mctnetwork('pl_cort')
mcf <- scdb_mctnetflow('pl_cort')
mc2d <- scdb_mc2d('pl_cort_not_cor_cc')

legc <- log2(1e-05 + mc@e_gc)


mcmd <- readr::read_tsv(file.path(wd, 'BonevCollab/mcmd_pl_cort.tsv'))
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))


nsc_mcs <- which(mcmd$cell_type == 'NSC')
ipc_cyc_mcs <- which(mcmd$cell_type == 'IPC_cyc')
ipc_mcs <- which(mcmd$cell_type %in% c('IPC'))
neuron_mcs <- which(mcmd$cell_type %in% grep('i', grep('PN', unique(mcmd$cell_type), v=T), inv = T, v=T))
cfupn_mcs <- which(mcmd$cell_type %in% c('SCPN', 'CthPN'))
cpn_mcs <- which(mcmd$cell_type %in% c('CPN_L2-3','CPN_L5_6'))

glia_mcs <- which(mcmd$cell_type %in% grep('cyte|OPC', unique(mcmd$cell_type), v=T))

nsc_genes_vs_ipc <- get_genes_specific_to_mcs(legc = legc, mc_pos = nsc_mcs, mc_neg = ipc_cyc_mcs)

nsc_genes_vs_glia <- get_genes_specific_to_mcs(legc = legc, mc_pos = nsc_mcs, mc_neg = glia_mcs)

cfupn_vs_ipc <- get_genes_specific_to_mcs(legc = legc, mc_pos = cfupn_mcs, mc_neg = ipc_mcs)
cpn_vs_ipc <- get_genes_specific_to_mcs(legc = legc, mc_pos = cpn_mcs, mc_neg = ipc_mcs)

nsc_module <- union(names(nsc_genes_vs_ipc[nsc_genes_vs_ipc > 2]), names(nsc_genes_vs_glia[nsc_genes_vs_glia > 2]))

cfupn_module <- names(cfupn_vs_ipc[cfupn_vs_ipc >1])
cpn_module <- names(cpn_vs_ipc[cpn_vs_ipc >1])

neuron_genes <- get_genes_specific_to_mcs(legc = legc, mc_pos = neuron_mcs, mc_neg = ipc_mcs)

neuron_module <- names(neuron_genes[neuron_genes > 2])

neu_score <- colSums(legc[neuron_module,])
nsc_score <- colSums(legc[nsc_module,])
vnn <- neu_score - nsc_score
bins_nsc_neu <- setNames(droplevels(cut(vnn, breaks = seq(min(vnn), max(vnn), l = 21))), names(nsc_score))

levels(bins_nsc_neu) <- 1:length(levels(bins_nsc_neu))

tbl_vnn_ct <- table(bins_nsc_neu, mcmd$cell_type)
tbl_vnn_ct <- tbl_vnn_ct[,!(colnames(tbl_vnn_ct) %in% c('Astrocytes', 'OPCs'))]
tbl_vnn_ct_norm <- tbl_vnn_ct/rowSums(tbl_vnn_ct)
barplot(t(tbl_vnn_ct_norm), col = col_key[colnames(tbl_vnn_ct_norm)])

tbl_vnn_ct_cfupn <- tbl_vnn_ct[,!(colnames(tbl_vnn_ct) %in% c('iCPN_early', 'iCPN_late', 'CPN_L2-3', 'CPN_L5_6'))]

tbl_vnn_ct_cfupn_norm <- tbl_vnn_ct_cfupn/rowSums(tbl_vnn_ct_cfupn)

tbl_vnn_ct_cpn <- tbl_vnn_ct[,!(colnames(tbl_vnn_ct) %in% c('iCPN/CfuPN', 'iCfuPN', 'SCPN', 'CthPN'))]

tbl_vnn_ct_cpn_norm <- tbl_vnn_ct_cpn/rowSums(tbl_vnn_ct_cpn)




mc_pf <- 0.1
ct_prop_ls <- lapply(list(c('CthPN', 'SCPN'), c('CPN_L5_6', 'CPN_L2-3')), function(cti) {
    lapply(2:5, function(tt) {
        mc_p = rep(0,ncol(mc@e_gc))
        mc_p[mcmd$cell_type %in% cti] = mct@mc_t[mcmd$cell_type %in% cti,tt]
        mc_p = mc_pf*mc_p/sum(mc_p)
        card_prop = mctnetflow_propogate_from_t(mcf = mcf, t = tt, mc_p = mc_p)
        return(card_prop)
    })
})

cfupn_traj_mcs <-  which(rowSums(do.call('cbind', lapply(ct_prop_ls[[1]], function(x) x$probs))) > 1e-4)

cpn_traj_mcs <-  as.numeric(mcmd$metacell[mcmd$cell_type %in% c('NSC', 'IPC', 'IPC_cyc', 'iCPN_early', 'iCPN_late', 'CPN_L5_6', 'CPN_L2-3') & 
                                          mcmd$mean_day >= 14 & mcmd$mean_day <= 18])

cfupn_traj_mcs_ord <- cfupn_traj_mcs[order(bins_nsc_neu[cfupn_traj_mcs])]
pcu_cfupn <- princurve::principal_curve(x = cbind(mc2d@mc_x[cfupn_traj_mcs_ord], 
                                              mc2d@mc_y[cfupn_traj_mcs_ord]))

cpn_traj_mcs_ord <- cpn_traj_mcs[order(bins_nsc_neu[cpn_traj_mcs])]
pcu_cpn <- princurve::principal_curve(cbind(mc2d@mc_x[cpn_traj_mcs_ord], 
                                            mc2d@mc_y[cpn_traj_mcs_ord]))

save(vnn, bins_nsc_neu, 
     tbl_vnn_ct_norm, 
     tbl_vnn_ct_cfupn_norm, 
     tbl_vnn_ct_cpn_norm, 
     pcu_cfupn,
     pcu_cpn,
     file = './output/metacell_model/diff_order_data.rda')
