library('metacell')
library(dplyr)
library(pheatmap)
wd = '/home/feshap/raid/proj/mmcortex'
setwd(wd)
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
devtools::load_all('~/src//metacell.flow')
scdb_flow_init()

SEED = 1337
set.seed(SEED)
scfigs_init("figs/")

nm = 'pl_cort'
# nm_test <- 'pl_cort_test_colors_reannot'
nm_test <- 'pl_cort'
mc_id <- nm_test
mat_id <- nm
mct_id <- nm_test
flow_id <- nm_test


fig_width <- 2000
fig_height <- 2000

tgconfig::set_param('mcp_heatmap_height', fig_width, 'metacell')
tgconfig::set_param('mcp_heatmap_width', fig_height, 'metacell')

mc2d_path <- glue::glue('./figs/{nm_test}_mc2d.png')
flow_plot_path <- glue::glue('./figs/{nm_test}_flow.png')
flow_plot_w_id_path <- glue::glue('./figs/{nm_test}_flow_mc_id.png')

mc = scdb_mc(mc_id)
mat = scdb_mat(mat_id)
legc <- log2(1e-5 + mc@e_gc)
mct = scdb_mctnetwork(mct_id)
mcf = scdb_mctnetflow(flow_id)

mg_bon = read.delim('./BonevCollab//marker_genes.tsv', sep='\t') %>% 
            apply(1, function(x) stringr::str_split(x, pattern = ' ')) %>%
            purrr::map(1) %>% purrr::map(function(x) c(x[[1]], stringr::str_split(x[[length(x)]], ',')))

neuron_st = c('CPN_L2-3','CPN_L5_6','SCPN','CthPN')
cpn = c('CPN_L5_6', 'CPN_L2-3')
cfupn = c('CthPN', 'SCPN')

### Import marker gene set ###
mg_bon = unique(data.frame(cbind(purrr::map(mg_bon, 1), purrr::map(mg_bon, 2))))
colnames(mg_bon) = c('st', 'marks')
mg_bon = mg_bon[mg_bon$st %in% c('NSC', 'IPC', 'CthPN', 'SCPN', 'CPN_L2-3', 'Stellate_L4', 'CPN_L5_6', 'OPCs', 'Astrocytes'),]

mg_bon$marks[mg_bon$st == 'Astrocytes'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'Astrocytes'])))
mg_bon$marks[mg_bon$st == 'OPCs'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'OPCs'])))
mg_bon$marks[mg_bon$st == 'SCPN'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'SCPN'])))
mg_bon$marks[mg_bon$st == 'CthPN'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CthPN'], 'Foxp2', 'Hs3st4', 'Prickle1')))
mg_bon$marks[mg_bon$st == 'CPN_L2-3'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L2-3'], c('Eif1b', 'Gpm6a'))))
mg_bon$marks[mg_bon$st == 'CPN_L5_6'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L5_6'], 'Gucy1a1', 'Gucy1b1', 'Syt4')))

x = unlist(mg_bon$marks[mg_bon$st == 'CthPN'])
x = list(x[x != 'Sox5'])
mg_bon$marks[mg_bon$st == 'CthPN'] = x

### Define additional properties of metacells ###
add_df = data.frame(rbind(
    list(st = 'Mature', marks = c('Mef2c', 'Mapt', 'Stmn2')),
    list(st = 'Immature', marks = c('Neurod1', 'Sstr2', 'Rnd2')),
    list(st = 'Cycling', marks = c('Pcna', 'Top2a', 'Mki67', 'Ube2c'))
))

mg_bon = rbind(mg_bon, add_df)
mg_bon_out = as.data.frame(apply(mg_bon, 2, function(x) sapply(x, function(y) paste0(y, collapse = ', '))))
readr::write_csv(mg_bon_out, file = "./data/pl_mg_bon.csv")

fp_marks <- mc@mc_fp[sort(unique(unlist(mg_bon$marks))),]

### Calculate additional properties of metacells ###
st_mark_mc = sapply(mg_bon$marks, function(x) apply(fp_marks[as.character(x),], 2, mean))
st_mark_mc = data.frame(st_mark_mc)
colnames(st_mark_mc) = unlist(mg_bon$st)
st_mark_mc$CPN = apply(st_mark_mc[,cpn], 1, max)
st_mark_mc$CfuPN = apply(st_mark_mc[,cfupn], 1, max)
cat_mark_mc = st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')]
cat_mark_mc$cpn_cfu_rat = log2(cat_mark_mc$CPN/cat_mark_mc$CfuPN)
cat_mark_mc$mat_imm_rat = log2(cat_mark_mc$Mature/cat_mark_mc$Immature)
st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')] = c()


### Tabulate metacells by single cells from each day ###
day_mat = matrix(0, nrow=length(unique(mc@mc)), ncol = length(unique(mat@cell_metadata[names(mc@mc),'day'])))
colnames(day_mat) = sort(unique(mat@cell_metadata[names(mc@mc),'day']))
rownames(day_mat) = 1:nrow(day_mat)
for (i in 1:nrow(day_mat)) {
    tbl = table(mat@cell_metadata[names(mc@mc)[mc@mc == i], 'day'])
    day_mat[i,names(tbl)] = tbl
}

max_day = apply(day_mat, 1, which.max)

mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort_freeze_1_9_22.tsv')

color_key <- unique(mcmd[,c('st', 'color')])
colnames(color_key) <- c('cell_type', 'color')

options(repr.plot.width = 16, repr.plot.height = 16)

plot(cat_mark_mc$mat_imm_rat, cat_mark_mc$cpn_cfu_rat, col = mc@colors, pch = 16, cex = 2)
options(repr.plot.width = 8, repr.plot.height = 8)

### Define thresholds for assigning cell types based on markers and derived quantities ###
cond1 = apply(st_mark_mc, 1, max) > 2                                   ## filter for high st specificity
cond2 = cat_mark_mc$mat_imm_rat > 2 & 
                    legc['Mef2c',] >= -10.5 & 
                    legc['Mapt',] >= -10.5 
                    ## filter for mature mcs
cond3 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] %in% 
                c('NSC', 'IPC', 'Astrocytes', 'OPCs')       ## Check neural or glial
cond7 = cat_mark_mc$cpn_cfu_rat >= log2(2)                              ## filter for CPN-tending
cond8 = cat_mark_mc$cpn_cfu_rat <= -log2(2)                             ## filter for CfuPN-tending
cond9 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'NSC'
cond10 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'IPC'
cond11 = cat_mark_mc$Cycling >= 2
cond12 = cond2 & (cond7 | cond8)

st = ifelse(cond1 & cond3, colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)], NA)
st[is.na(st) & cond9] = 'NSC'
st[cond10 & cond11] = 'IPC_cyc'
st[st == 'Astrocytes' & st_mark_mc$Astrocytes < 4 & rowSums(day_mat[,c('E13','E14','E15', 'E16')]) >= 4] <- 'NSC'

mature_inds = as.numeric(which((cond1 & cond2) & is.na(st)))
st[mature_inds] = colnames(st_mark_mc)[apply(st_mark_mc[mature_inds,], 1, which.max)]

net = mct@network
net$flow = mcf@edge_flows
net = net[net$flow > 0,]
net = net[net$type1 == 'norm_f',]
net = net[net$type2 != 'sink',]

mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"day"])
mc_ag_n = mc_ag/rowSums(mc_ag)
mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
ag_cn_com = apply(mc_ag_cn, 1, function(x) sum(x*(13:18))/sum(x))
nc = ncol(mc_ag_cn)

#####################################################
### Don't remove - how gene modules were calculated

# mcmd_new <- readr::read_tsv('./BonevCollab/mcmd_pl_cort_freeze_1_9_22.tsv')
# fp_st <- t(tgs_matrix_tapply(mc@mc_fp, mcmd_new$st, mean))
# gene_modules <- setNames(lapply(1:ncol(fp_st), function(i) names(head(fp_st[order(fp_st[,i], decreasing = T),i], 10))), colnames(fp_st))
# gene_modules <- gene_modules[names(gene_modules) != 'IPC_late']
# save(gene_modules, file = './data/gene_modules_mcmd_pl_cort.Rda')

### Don't remove - how gene modules were calculated
#####################################################


## Load gene modules and assign remaining metacells based on module scores

load('./data/gene_modules_mcmd_pl_cort.Rda')

mc_module_score <- as.data.frame(sapply(gene_modules, function(mi) colSums(mc@mc_fp[mi,])))

st[which(is.na(st))] <- colnames(mc_module_score)[apply(mc_module_score[which(is.na(st)),], 1, which.max)]

pairs <- as.data.frame(list(ct1 = c('SCPN', 'iCfuPN', 'iCPN/CfuPN', 'iCPN/CfuPN',  'iCfuPN','CPN_L2-3', 'iCfuPN'), 
                            ct2 = c('CPN_L5_6', 'iCPN/CfuPN', 'iCPN_early', 'iCPN_late', 'iCPN_late', 'CPN_L5_6','CPN_L5_6')))

st_fp <- t(tgs_matrix_tapply(mc@mc_fp, st, mean))

c2p <- which(st %in% c('SCPN', 'CthPN','iCfuPN', 'iCPN/CfuPN', 'iCPN_late', 'iCPN_early', 'CPN_L5_6', 'CPN_L2-3'))
for (i in 1:nrow(pairs)) {
    g1m2 <- head(rownames(st_fp[order(st_fp[,pairs$ct1[[i]]] - st_fp[,pairs$ct2[[i]]], decreasing = T),]), 10)
    g2m1 <- head(rownames(st_fp[order(st_fp[,pairs$ct2[[i]]] - st_fp[,pairs$ct1[[i]]], decreasing = T),]), 10)
    x <- colSums(legc[g1m2,])
    y <- colSums(legc[g2m1,])
    st[which(x > y & st %in% pairs[i,])] <- pairs[i,'ct1']
    st[which(y > x & st %in% pairs[i,])] <- pairs[i,'ct2']
    cur_cols <- color_key$color[match(st, color_key$cell_type)]   
    png(glue::glue('./figs/cell_type_signature_pairs/{gsub("/", "-", pairs$ct1[[i]])}_vs_{gsub("/", "-", pairs$ct2[[i]])}.png'), h = 1600, w = 1600)
    plot(colSums(legc[g1m2,c2p]), colSums(legc[g2m1,c2p]), col = cur_cols[c2p], pch = 16, cex = 2, xlab = pairs$ct1[[i]], ylab = pairs$ct2[[i]])
    text(colSums(legc[g1m2,c2p]), -.5+colSums(legc[g2m1,c2p]), labels = c2p, cex = 1)
    abline(0,1,lty=2,col='red')
    dev.off()
}

tot_outflow_thresh <- 0.25
st_bad <- c('CPN_L5_6', 'CPN_L2-3', 'iCPN_late')
mcs_that_flow_to_CPN <- which(colSums(do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(x[,st %in% st_bad])))) >= tot_outflow_thresh)

st_bad <- c('SCPN', 'CfuPN')
mcs_that_flow_to_CfuPN <- which(colSums(do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(x[,st %in% st_bad])))) >= tot_outflow_thresh)
imm_inds <- which(st %in% c('iCPN/CfuPN', 'iCfuPN', 'iCPN_early', 'iCPN_late'))
icpn_cfupn_inds <- intersect(imm_inds, intersect(mcs_that_flow_to_CPN, mcs_that_flow_to_CfuPN))

icpn_inds <- which(st %in% c('iCPN_early', 'iCPN_late'))
bad_icpn_inds <- intersect(icpn_inds, mcs_that_flow_to_CfuPN)
icpn_late_inds <- setdiff(intersect(imm_inds,mcs_that_flow_to_CPN), union(which(st %in% cpn), icpn_cfupn_inds))
bad_icfupn_inds <- which(st == 'iCfuPN')[which(st == 'iCfuPN') %in% union(icpn_cfupn_inds, icpn_late_inds)]
st[union(bad_icfupn_inds,bad_icpn_inds)] <- 'iCPN/CfuPN'

st_new <- st


### Manual changes following MCView analysis
st_new[c(190,191,193,194)] <- 'iCPN_late'
st_new[c(444,445,448,858)] <- 'iCPN/CfuPN'
st_new[c(234,560,561)] <- 'iCfuPN'
st_new[c(621)] <- 'NSC'
# st_new[c()] <- 'CPN_L5_6'
# st_new[c()] <- 'CPN_L2-3'
# st_new[c()] <- 'CthPN'
# st_new[c()] <- 'iCPN_early'
st <- st_new


### From here on it's just saving data and plotting

color_key = color_key %>% tibble::column_to_rownames('cell_type')

color_key$ord = 1:nrow(color_key)
color_key = tibble::rownames_to_column(color_key)
color_key = color_key %>% mutate(i = 1:nrow(color_key)) 
colnames(color_key) = c('cell_type', 'color', 'i', 'ord')
types_df = data.frame(cbind(1:length(st), st))
colnames(types_df) = c('metacell', 'cell_type')
readr::write_tsv(color_key, '/net/mraid20/export/tgdata/users/aviezerl/proj/mmcortex/cell_type_annot_YS.tsv')
types_df = data.frame(cbind(1:length(st), st))
colnames(types_df) = c('metacell', 'cell_type')
types_df$mc_age <- 12+apply(day_mat, 1, function(x) sum(x*(1:ncol(day_mat)))/sum(x))
readr::write_csv(types_df, '/net/mraid20/export/tgdata/users/aviezerl/proj/mmcortex/metacell_types.csv')

mc@colors = color_key$color[match(st, color_key[,'cell_type'])]
scdb_add_mc(nm_test, mc)

cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc','IPC','iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st_new = unlist(lapply(cust_st_ord, function(s) setNames(which(st == s), 
                                                                  rep(s, length(which(st == s))))))


mcell_mc2d_plot(nm_test, colors = mc@colors, fig_fn = mc2d_path)
mct@mc_id = nm_test
mcf@net_id = nm_test

scdb_add_mctnetwork(nm_test, mct)
scdb_add_mctnetflow(nm_test, mcf)


mctnetwork_plot_net(mct_id = mct_id, 
                    flow_id = flow_id, 
                    fn = flow_plot_path, 
                    # flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0]),  
                    colors_ordered = color_key$color[match(cust_st_ord, color_key[,'cell_type'])], 
                    plot_mc_ids = F,
                    h = fig_height, 
                    w = fig_width)
mctnetwork_plot_net(mct_id = mct_id, 
                    flow_id = flow_id, 
                    fn = flow_plot_w_id_path, 
                    # flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0]),
                    colors_ordered = color_key$color[match(cust_st_ord, color_key[,'cell_type'])], 
                    plot_mc_ids = T,
                    h = 2*fig_height, 
                    w = 2*fig_width)

mcmd_new = data.frame(cbind(as.numeric(1:length(st)), st, color_key$color[match(st, color_key[,'cell_type'])]))
colnames(mcmd_new) = c('metacell', 'cell_type', 'color')
mcmd_new$color[is.na(mcmd_new$color)] = gplots::col2hex('seagreen3')
load(file.path(wd, 'data', glue::glue('{nm}_cc_score.rda')))
mcmd_new$cc_score <- cc_score
rownames(mcmd_new) = 1:nrow(mcmd_new)
mcmd_new = dplyr::left_join(tibble::rownames_to_column(mcmd_new), 
            tibble::rownames_to_column(data.frame(day_mat)), by='rowname') %>% 
            tibble::column_to_rownames('rowname')

mcmd_new$mean_day <- types_df$mc_age

tb = cut(mcmd_new$mean_day, breaks = seq(13,18,l=14))
mcmd_new$time_bin = match(tb, levels(tb))
readr::write_tsv(x=mcmd_new, file=glue::glue('./BonevCollab/mcmd_{nm_test}.tsv'))