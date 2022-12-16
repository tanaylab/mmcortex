library('metacell')
library(dplyr)

wd = '/home/feshap/raid/proj/mmcortex'
setwd(wd)
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)

SEED = 1337
K = 16
set.seed(SEED)
scfigs_init("figs/")



nm = 'pl_cort'

tgconfig::set_param('mcp_heatmap_height', 2500, 'metacell')
tgconfig::set_param('mcp_heatmap_width', 4500, 'metacell')

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mc = scdb_mc(nm)

mat = scdb_mat(nm)

mg_bon = read.delim('./BonevCollab//marker_genes.tsv', sep='\t') %>% 
            apply(1, function(x) stringr::str_split(x, pattern = ' ')) %>%
            purrr::map(1) %>% purrr::map(function(x) c(x[[1]], stringr::str_split(x[[length(x)]], ',')))

mg_bon = unique(data.frame(cbind(purrr::map(mg_bon, 1), purrr::map(mg_bon, 2))))
colnames(mg_bon) = c('st', 'marks')
mg_bon = mg_bon[mg_bon$st %in% c('NSC', 'IPC', 'CthPN', 'SCPN', 'Stellate_L4', 'CPN_L2-3', 'CPN_L5_6', 'Oligodendrocytes', 'Astrocytes'),]

# mg_bon$marks[mg_bon$st == 'Astrocytes'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'Astrocytes'], 'Slc1a3', 'Cst3')))
# mg_bon$marks[mg_bon$st == 'Oligodendrocytes'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'Oligodendrocytes'], 'Cst3')))
mg_bon$marks[mg_bon$st == 'Astrocytes'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'Astrocytes'], 'Fabp7', 'Slc1a3')))
mg_bon$marks[mg_bon$st == 'Oligodendrocytes'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'Oligodendrocytes'], "Egr1")))
mg_bon$marks[mg_bon$st == 'SCPN'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'SCPN'])))
mg_bon$marks[mg_bon$st == 'CthPN'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CthPN'], 'Foxp2', 'Hs3st4', 'Prickle1')))
mg_bon$marks[mg_bon$st == 'CPN_L2-3'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L2-3'], c('Eif1b', 'Gpm6a'))))
mg_bon$marks[mg_bon$st == 'CPN_L5_6'] = list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L5_6'], 'Gucy1a1', 'Gucy1b1', 'Syt4')))

x = unlist(mg_bon$marks[mg_bon$st == 'CthPN'])
x = list(x[x != 'Sox5'])
mg_bon$marks[mg_bon$st == 'CthPN'] = x

                                         
add_df = data.frame(rbind(
    list(st = 'Mature', marks = c('Mef2c', 'Mapt', 'Stmn2')),
    list(st = 'Immature', marks = c('Neurod1', 'Sstr2', 'Rnd2')),
    list(st = 'Cycling', marks = c('Pcna', 'Top2a', 'Mki67', 'Ube2c'))
))

mg_bon = rbind(mg_bon, add_df)
mg_bon

mg_bon_out = as.data.frame(apply(mg_bon, 2, function(x) sapply(x, function(y) paste0(y, collapse = ', '))))

readr::write_csv(mg_bon_out, file = "./data/pl_mg_bon.csv")

legcz_marks = mc@mc_fp[sort(unique(unlist(mg_bon$marks))),]

st_mark_mc = sapply(mg_bon$marks, function(x) apply(legcz_marks[as.character(x),], 2, mean))
st_mark_mc = data.frame(st_mark_mc)
colnames(st_mark_mc) = unlist(mg_bon$st)
st_mark_mc$CPN = apply(st_mark_mc[,c('CPN_L2-3', 'CPN_L5_6')], 1, max)
st_mark_mc$CfuPN = apply(st_mark_mc[,c('SCPN', 'CthPN')], 1, max)
cat_mark_mc = st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')]
cat_mark_mc$cpn_cfu_rat = log2(cat_mark_mc$CPN/cat_mark_mc$CfuPN)
cat_mark_mc$mat_imm_rat = log2(cat_mark_mc$Mature/cat_mark_mc$Immature)
st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')] = c()

day_mat = matrix(0, nrow=length(unique(mc@mc)), ncol = length(unique(mat@cell_metadata[names(mc@mc),'day'])))
colnames(day_mat) = sort(unique(mat@cell_metadata[names(mc@mc),'day']))
rownames(day_mat) = 1:nrow(day_mat)
for (i in 1:nrow(day_mat)) {
    tbl = table(mat@cell_metadata[names(mc@mc)[mc@mc == i], 'day'])
    day_mat[i,names(tbl)] = tbl
}


max_day = apply(day_mat, 1, which.max)

cond1 = apply(st_mark_mc, 1, max) > 2                                   ## filter for high st specificity
cond2 = cat_mark_mc$mat_imm_rat > 2                               ## filter for mature mcs
cond3 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] %in% 
                c('NSC', 'IPC', 'Astrocytes', 'Oligodendrocytes')       ## Check neural or glial
# cond5 = cat_mark_mc$mat_imm_rat < -1                                    ## filter for immature mcs
# cond4 = abs(cat_mark_mc$mat_imm_rat) <= 1                               ## filter for intermediate maturity
# cond6 = abs(cat_mark_mc$cpn_cfu_rat) < log2(2)                          ## filter for indeterminate CfuPN/CPN state
cond7 = cat_mark_mc$cpn_cfu_rat >= log2(2)                              ## filter for CPN-tending
cond8 = cat_mark_mc$cpn_cfu_rat <= -log2(2)                             ## filter for CfuPN-tending
cond9 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'NSC'
cond10 = colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'IPC'
cond11 = cat_mark_mc$Cycling >= 2
cond12 = cond2 & (cond7 | cond8)

st = ifelse(cond1 & cond3, colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)], NA)

st[is.na(st) & cond9] = 'NSC'
st[cond10 & cond11] = 'IPC_cyc'

mature_inds = as.numeric(which((cond1 & cond2) & is.na(st)))
st[mature_inds] = colnames(st_mark_mc)[apply(st_mark_mc[mature_inds,], 1, which.max)]

neuron_st = c('CPN_L2-3','Stellate_L4','CPN_L5_6','SCPN','CthPN')
cpn = c('CPN_L5_6', 'CPN_L2-3','Stellate_L4')
cfupn = c('CthPN', 'SCPN')

devtools::load_all('~/src//metacell.flow')

scdb_flow_init()

mat = scdb_mat(nm)

mct = scdb_mctnetwork(nm)
mcf = scdb_mctnetflow(nm)

net = mct@network
net$flow = mcf@edge_flows

net = net[net$flow > 0,]
net = net[net$type1 == 'norm_f',]
net = net[net$type2 != 'sink',]

mcmd = vroom::vroom('./BonevCollab//mcmd_pl.tsv')

mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"day"])
mc_ag_n = mc_ag/rowSums(mc_ag)
mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)

ag_cn_com = apply(mc_ag_cn, 1, function(x) sum(x*(13:18))/sum(x))

nc = ncol(mc_ag_cn)

mcs_by_st = lapply(unique(st), function(u) setNames(mcmd$mean_day[which(st == u)], which(st == u)))
names(mcs_by_st) = unique(st)
# mcs_by_st

st_mc_t = lapply(mcs_by_st[neuron_st[neuron_st %in% names(mcs_by_st)]], function(x) round(x)-12)
names(st_mc_t) = names(mcs_by_st[neuron_st[neuron_st %in% names(mcs_by_st)]])
# st_mc_t

# mcs_by_st

get_st_flows = function(i, md) {
#     print(c(i, md[[i]], ncol(imm_flows[[i]]$probs)))
    sm = apply(imm_flows[[i]]$probs[,md[[i]]:ncol(imm_flows[[i]]$probs)], 1, sum);
    setNames(sapply(seq_along(st_mc_t), function(x,i) {
    sum(sm[as.numeric(unlist(names(x[[i]])))])
        }, x = st_mc_t), names(st_mc_t))
}

get_imm_flows = function(st) {
    y =  which(is.na(st))
    mdy = max_day[y]
    y =  y[mdy<max(mdy)]
    mdy = mdy[mdy<max(mdy)]
    imm_flows = lapply(seq_along(y), function(y,md,i) {
        p = mc_ag_cn[,md[[i]]];
        p = ifelse(1:length(p) == y[[i]], p, 0); 
        return(mctnetflow_propogate_from_t(mcf, md[[i]], mc_p = p))
    }, md = mdy, y = y)

    names(imm_flows) = y
    return(imm_flows)

#     sum_imm_flows = lapply(seq_along(y), function(j) get_st_flows(j, mdy))
#     names(sum_imm_flows) = y
#     return(sum_imm_flows)
}

# assign_mcs = function(st) {
    
y =  which(is.na(st))
mdy = max_day[y]
y =  y[mdy<max(mdy)]
mdy = mdy[mdy<max(mdy)]
imm_flows = lapply(seq_along(y), function(y,md,i) {
    p = mc_ag_cn[,md[[i]]];
    p = ifelse(1:length(p) == y[[i]], p, 0); 
    return(mctnetflow_propogate_from_t(mcf, md[[i]], mc_p = p))
}, md = mdy, y = y)

names(imm_flows) = y

sum_imm_flows = lapply(seq_along(y), function(j) get_st_flows(j, mdy))
names(sum_imm_flows) = y

sum_flow_mat = data.frame(t(sapply(sum_imm_flows, function(x) x)))
colnames(sum_flow_mat) = neuron_st[neuron_st %in% names(mcs_by_st)]
rownames(sum_flow_mat) = as.numeric(rownames(sum_flow_mat))

# st_1 = apply(sum_flow_mat, 1, function(x) if (length(which(x > 0)) == 1) {return(which.max(x))}
#     else {return(NA)})
# st_1 = st_1[!unlist(lapply(st_1, is.na))]
# st_1 = setNames(colnames(sum_flow_mat)[st_1], names(st_1))
# st_1

# tbl = table(unlist(apply(sum_flow_mat, 1, function(x) if (length(which(x > 0)) == 1) {return(which.max(x))})))
# names(tbl) = neuron_st[as.numeric(names(tbl))]
# tbl

cpn_tot = apply(sum_flow_mat[,colnames(sum_flow_mat) %in% cpn], 1, sum)
cfupn_tot = apply(sum_flow_mat[,colnames(sum_flow_mat) %in% cfupn], 1, sum)

# icpn_l23_inds = as.numeric(names(st_1)[st_1 %in% c('Stellate_L4','CPN_L2-3')])
# icpn_l56_inds = as.numeric(names(st_1)[st_1 == 'CPN_L5_6'])
icpn_inds = rownames(sum_flow_mat)[cpn_tot > 0 & cfupn_tot == 0 
#                                    & !(rownames(sum_flow_mat) %in% union(icpn_l23_inds, icpn_l56_inds))
                                  ] %>% as.numeric %>% unique
icfupn_inds = c(rownames(sum_flow_mat)[cpn_tot == 0 & cfupn_tot > 0]
#                 , as.numeric(names(st_1)[st_1 %in% c('SCPN','CthPN')])
               ) %>% as.numeric %>% unique
icpn_cfupn_inds = rownames(sum_flow_mat)[cpn_tot > 0 & cfupn_tot > 0 
#                                          & !(rownames(sum_flow_mat) %in% union(icpn_l23_inds, icpn_l56_inds))
                                        ] %>% as.numeric %>% unique



# length(icpn_l23_inds)
# length(icpn_l56_inds)
# length(icpn_cfupn_inds)
# length(icpn_inds)
# length(icfupn_inds)
# length(intersect(icpn_cfupn_inds, icfupn_inds))
# length(intersect(icpn_inds, icfupn_inds))
# length(intersect(icpn_cfupn_inds, icpn_inds))
# length(intersect(icpn_inds, icpn_l23_inds))
# length(intersect(icpn_inds, icpn_l56_inds))
# length(intersect(icpn_cfupn_inds, icpn_l23_inds))
# length(intersect(icpn_cfupn_inds, icpn_l56_inds))
# length(intersect(icpn_l23_inds, icpn_l56_inds))
# length(intersect(icpn_l23_inds, icfupn_inds))
# length(intersect(icfupn_inds, icpn_l56_inds))
# length(c(icpn_inds, icpn_cfupn_inds, icpn_l23_inds, icpn_l56_inds, icfupn_inds))
# length(union(icpn_inds, union(icpn_cfupn_inds, union(icpn_l23_inds, union(icpn_l56_inds, icfupn_inds)))))

# # st[icpn_l23_inds] = 'iCPN_L2-3'
# # st[icpn_l56_inds] = 'iCPN_L5_6'
st[icpn_inds] = 'iCPN_late'
st[icfupn_inds] = 'iCfuPN'
st[icpn_cfupn_inds] = 'iCPN/CfuPN'
# st[is.na(st)] = 'iCPN'

feats = scdb_gset('pl_f')
feats = names(feats@gene_set)
marks = scdb_gset('pl_cort_marks_f')
marks = names(marks@gene_set)
legc = log2(1e-07 + mc@e_gc)
legc_marks = legc[marks,]
SEED = 1337
K = 16
# mark_fp_km = tglkmeans::TGL_kmeans(t(legc[feats[feats %in% rownames(mc@e_gc)],]), seed = SEED, k = K)
mark_fp_km = tglkmeans::TGL_kmeans(t(legc[union(marks, unlist(mg_bon$marks)),]), seed = SEED, k = K)

mean_over_km = t(tgs_matrix_tapply(legc, mark_fp_km$cluster, mean))

mok_hc = hclust(dist(t(mean_over_km)))

mok_hc_ord = unlist(lapply(mok_hc$order, function(u) mcmd$mc[mark_fp_km$cluster == u]))

st[mark_fp_km$cluster == 1 & (is.na(st) | !(st %in% c('SCPN', 'CPN_L5_6')))] = 'CPN_L5_6'
st[mark_fp_km$cluster == 1 & log2(mc@e_gc['Satb2',]) > -15 & log2(mc@e_gc['Bcl11b',]) > -14  & log2(mc@e_gc['Bcl11b',]) < -11] = 'iCPN_late'
st[mark_fp_km$cluster == 2 & (is.na(st) | st %in% c('iCPN_L5_6','iCPN/CfuPN'))] = 'iCfuPN'
st[mark_fp_km$cluster == 3] = 'iCfuPN'
st[mark_fp_km$cluster == 4 & (is.na(st) | !(st %in% cpn))] = 'iCPN_late'
st[mark_fp_km$cluster == 5] = 'iCPN/CfuPN'
st[mark_fp_km$cluster == 7] = 'iCPN_early'
st[mark_fp_km$cluster == 6 & is.na(st)] = 'iCPN_early'
st[mark_fp_km$cluster == 15] = 'IPC_late'


cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
# cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
#                                                                   rep(st, length(which(mcmd$st == st))))))
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(st == s), 
                                                                  rep(s, length(which(st == s)))
                                                                )
                                )
                        )                         

color_key = vroom::vroom('./data/pl_cort_color_key.tsv')
color_key = color_key %>% tibble::column_to_rownames('st')
color_key$ord = 1:nrow(color_key)
color_key = tibble::rownames_to_column(color_key)
color_key = color_key %>% mutate(i = 1:nrow(color_key)) 
colnames(color_key) = c('cell_type', 'color', 'i', 'ord')
# print(color_key)
readr::write_tsv(color_key, '/net/mraid14/export/tgdata/users/aviezerl/proj/mmcortex/cell_type_annot_YS.tsv')
types_df = data.frame(cbind(1:length(st), st))
colnames(types_df) = c('metacell', 'cell_type')
readr::write_tsv(types_df, '/net/mraid14/export/tgdata/users/aviezerl/proj/mmcortex/types_df.tsv')
mc@colors = color_key$color[match(st, color_key[,'cell_type'])]


print(head(mc@colors))
scdb_add_mc('pl_cort', mc)

mc = scdb_mc('pl_cort')

mct@mc_id = 'pl_cort'
scdb_add_mctnetwork('pl_cort', mct)

mcf@net_id = 'pl_cort'

scdb_add_mctnetflow('pl_cort', mcf)

scdb_init(db_path, force_reinit = T)      

mctnetwork_plot_net(mct_id = 'pl_cort', flow_id = 'pl_cort', fn = './figs/pl_cort_flow.png', flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0]),  
    colors_ordered = color_key$color[match(cust_st_ord, color_key[,'cell_type'])], plot_mc_ids = F,h = 4000, w = 4000)
mctnetwork_plot_net(mct_id = 'pl_cort', flow_id = 'pl_cort', fn = './figs/pl_cort_flow_mc_id.png', flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0]),  
    colors_ordered = color_key$color[match(cust_st_ord, color_key[,'cell_type'])], plot_mc_ids = T,h = 4000, w = 4000)

genes_var = rownames(mc@mc_fp)[apply(mc@mc_fp, 1, function(x) length(which(x >= 3)) >= 4 & 
                                     !(st[which.max(x)] %in% c('Astrocytes', 'Oligodendrocytes')))]
# genes_var

lat = scdb_gset('pl_lateral')

length(genes_var[!(genes_var %in% names(lat@gene_set))])
genes_var = genes_var[!(genes_var %in% names(lat@gene_set))]

cor_all = tgs_cor(legc[feats[feats %in% rownames(legc)],], spearman = T)
sl_ord = slanter::slanted_orders(cor_all)
tfs = vroom::vroom('~/raid/Mus_musculus_TF.txt')

tfs = tfs$Symbol
tfs = tfs[tfs %in% rownames(mc@e_gc)]
tfs = tfs[which(apply(mc@mc_fp[tfs,], 1, function(x) max(x) > 3))]

mcell_mc2d_plot('pl_cort')

mcell_mc_plot_marks(mc_id = nm, mc_ord = cust_mc_ord_st, gset_id = 'pl_cort_marks_f', mat_id = nm, 
                    fig_fn = './figs/pl_cort_hm_ord_st_w_km.png', 
                    ext_metadata = setNames(mark_fp_km$cluster[mc@mc], names(mc@mc)))

mcell_mc_plot_marks('pl_cort', gset_id = 'pl_cort_marks_f', mat_id = nm, mc_ord = cust_mc_ord_st, 
                                fig_fn = './figs/pl_cort_hm_by_day.png', add_metadata = 'day')

mcell_mc_plot_marks('pl_cort', gene_list = tfs[order(apply(legc[tfs,sl_ord$rows], 1, which.max))], mat_id = nm, mc_ord = sl_ord$rows, 
                                fig_fn = './figs/pl_cort_hm_tfoi_sl_ord.png',
                   ext_metadata = setNames(mark_fp_km$cluster[mc@mc], names(mc@mc)))

mcell_mc_plot_marks('pl_cort', gene_list = genes_var[order(apply(mc@mc_fp[genes_var,cust_mc_ord_st], 1, function(x) sum(x*1:length(x))/sum(x)))], mat_id = nm, mc_ord = sl_ord$rows, 
                                fig_fn = './figs/pl_cort_hm_genes_var_sl_ord.png',
                   ext_metadata = setNames(mark_fp_km$cluster[mc@mc], names(mc@mc)))


# st_flows = lapply(st_mc_t, function(x) {lapply(seq_along(x), function(mci,n,i) {
#                                p = mc_ag_cn[,n[[i]]];
#                                p = ifelse(1:length(p) == mci[[i]], p, 0); 
#                                mctnetflow_propogate_from_t(mcf = mcf, t = n[[i]], mc_p = p)}, n = x, mci = as.numeric(names(x))
#                                               )
#                                         }
#                   )

# st_flow_mc2 = lapply(st_flows, function(x) {sort(unique(unlist(sapply(x, function(y) apply(y$probs, 2, function(z) which(z > 0))))))})

# st_flow_mc_mat = matrix(0, nrow=length(unique(mc@mc)), ncol = length(neuron_st))
# colnames(st_flow_mc_mat) = neuron_st

# st_flow_mc3 = lapply(seq_along(st_flows), function(x,y,i) {
#                     sort(unique(unlist(sapply(seq_along(x[[i]]), function(j) which(x[[i]][[j]]$probs[,y[[i]][[j]]-1] > 0)))))
#                                                                    }, x=st_flows, y=st_mc_t) 
# names(st_flow_mc3) = neuron_st

# st_flow_mc_mat3 = matrix(0, nrow=length(unique(mc@mc)), ncol = length(neuron_st))
# colnames(st_flow_mc_mat3) = neuron_st

# for (i in 1:length(st_flow_mc2)) {
#     st_flow_mc_mat[st_flow_mc2[[i]],i] = 1
# }

# for (i in 1:length(st_flow_mc3)) {
#     st_flow_mc_mat3[st_flow_mc3[[i]],i] = 1
# }


# icpn_l23_inds = which(st_flow_mc_mat3[,'CPN_L2-3'] > 0 & 
#                       apply(st_flow_mc_mat3[,neuron_st[!(neuron_st == 'CPN_L2-3')]], 1, sum) == 0 & is.na(st))
# st[icpn_l23_inds] = 'iCPN_L2-3'
# icpn_l56_inds = which(st_flow_mc_mat3[,'CPN_L5_6'] > 0 & 
#                       apply(st_flow_mc_mat3[,neuron_st[!(neuron_st == 'CPN_L5_6')]], 1, sum) == 0 & is.na(st))
# st[icpn_l56_inds] = 'iCPN_L5_6'

# icfupn_inds = which(apply(st_flow_mc_mat[,cpn], 1, sum) == 0 & apply(st_flow_mc_mat[,cfupn], 1, sum) > 0 & is.na(st))
# st[icfupn_inds] = 'iCfuPN'

# icpn_inds = which(apply(st_flow_mc_mat[,cpn], 1, sum) > 0 & apply(st_flow_mc_mat[,cfupn], 1, sum) == 0 & is.na(st))
# icpn_inds = icpn_inds[!(icpn_inds %in% union(icpn_l23_inds, icpn_l56_inds))]
# st[icpn_inds] = 'iCPN'


# icpn_cfupn_inds = which(apply(st_flow_mc_mat[,cpn], 1, sum) > 0 & apply(st_flow_mc_mat[,cfupn], 1, sum) > 0 & is.na(st))
# st[icpn_cfupn_inds] = 'iCPN/CfuPN'
                             
# st[is.na(st)] = 'iCPN'


# mcmd_old = vroom::vroom('./BonevCollab//mcmd_cort_6.tsv')
# color_key = unique(mcmd_old[,c('st','color')])
mcmd_new = data.frame(cbind(as.numeric(1:length(st)), st, color_key$color[match(st, color_key[,'cell_type'])]))
colnames(mcmd_new) = c('mc', 'st', 'color')
mcmd_new$color[is.na(mcmd_new$color)] = gplots::col2hex('yellow')

mcmd_new = dplyr::left_join(tibble::rownames_to_column(mcmd_new), 
            tibble::rownames_to_column(data.frame(day_mat)), by='rowname') %>% 
            tibble::column_to_rownames('rowname')

mcmd_new$mean_day = mcmd$mean_day

tb = cut(mcmd_new$mean_day, breaks = seq(13,18,l=14))
mcmd_new$time_bin = match(tb, levels(tb))

head(mcmd_new)
mc@colors = mcmd_new$color
scdb_add_mc('pl_cort', mc)
mcell_mc2d_plot('pl_cort', colors = mcmd_new$color)
mcell_mc2d_plot_by_factor('pl_cort', 'pl_cort', 'day', colors = mcmd_new$color, single_plot = T)

readr::write_tsv(x=mcmd_new, file='./BonevCollab/mcmd_pl_cort.tsv')