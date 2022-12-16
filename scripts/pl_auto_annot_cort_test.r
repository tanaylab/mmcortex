library('metacell')
library(dplyr)
library(pheatmap)
wd <- '/home/feshap/raid/proj/mmcortex'
setwd(wd)
db_path <- file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)

SEED <- 1337
K <- 20
set.seed(SEED)
scfigs_init("figs/")



nm <- 'pl_cort'

tgconfig::set_param('mcp_heatmap_height', 2500, 'metacell')
tgconfig::set_param('mcp_heatmap_width', 4500, 'metacell')

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mc <- scdb_mc(nm)

mat <- scdb_mat(nm)

mg_bon <- read.delim('./BonevCollab//marker_genes.tsv', sep='\t') %>% 
            apply(1, function(x) stringr::str_split(x, pattern = ' ')) %>%
            purrr::map(1) %>% purrr::map(function(x) c(x[[1]], stringr::str_split(x[[length(x)]], ',')))

mg_bon <- unique(data.frame(cbind(purrr::map(mg_bon, 1), purrr::map(mg_bon, 2))))
colnames(mg_bon) <- c('st', 'marks')
mg_bon <- mg_bon[mg_bon$st %in% c('NSC', 'IPC', 'CthPN', 'SCPN', 'CPN_L2-3', 'CPN_L5_6', 'Oligodendrocytes', 'Astrocytes'),]

# mg_bon$marks[mg_bon$st == 'Astrocytes'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'Astrocytes'], 'Slc1a3', 'Cst3')))
# mg_bon$marks[mg_bon$st == 'Oligodendrocytes'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'Oligodendrocytes'], 'Cst3')))
mg_bon$marks[mg_bon$st == 'Astrocytes'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'Astrocytes'], 'Fabp7', 'Slc1a3', 'Cst3')))
mg_bon$marks[mg_bon$st == 'Oligodendrocytes'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'Oligodendrocytes'], "Egr1")))
mg_bon$marks[mg_bon$st == 'SCPN'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'SCPN'])))
mg_bon$marks[mg_bon$st == 'CthPN'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'CthPN'], 'Foxp2', 'Hs3st4', 'Prickle1')))
mg_bon$marks[mg_bon$st == 'CPN_L2-3'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L2-3'], c('Eif1b', 'Gpm6a'))))
mg_bon$marks[mg_bon$st == 'CPN_L5_6'] <- list(unlist(c(mg_bon$marks[mg_bon$st == 'CPN_L5_6'], 'Gucy1a1', 'Gucy1b1', 'Syt4', "Cux1", "Cux2")))

tmp <- unlist(mg_bon$marks[mg_bon$st == 'CthPN'])
tmp <- list(tmp[tmp != 'Sox5'])
mg_bon$marks[mg_bon$st == 'CthPN'] <- tmp

                                         
add_df <- data.frame(rbind(
    list(st = 'Mature', marks = c('Mef2c', 'Mapt', 'Stmn2')),
    list(st = 'Immature', marks = c('Neurod1', 'Sstr2', 'Rnd2')),
    list(st = 'Cycling', marks = c('Pcna', 'Top2a', 'Mki67', 'Ube2c'))
))

mg_bon <- rbind(mg_bon, add_df)

mg_bon_out <- as.data.frame(apply(mg_bon, 2, function(x) sapply(x, function(y) paste0(y, collapse = ', '))))

readr::write_csv(mg_bon_out, file = "./data/pl_mg_bon.csv")

legc <- log2(1e-5 + mc@e_gc)
fp_marks <- mc@mc_fp[sort(unique(unlist(mg_bon$marks))),]

st_mark_mc <- sapply(mg_bon$marks, function(x) apply(fp_marks[as.character(x),], 2, mean))
st_mark_mc <- data.frame(st_mark_mc)
colnames(st_mark_mc) <- unlist(mg_bon$st)
st_mark_mc$CPN <- apply(st_mark_mc[,c('CPN_L2-3', 'CPN_L5_6')], 1, max)
st_mark_mc$CfuPN <- apply(st_mark_mc[,c('SCPN', 'CthPN')], 1, max)
cat_mark_mc <- st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')]
cat_mark_mc$cpn_cfu_rat <- log2(cat_mark_mc$CPN/cat_mark_mc$CfuPN)
cat_mark_mc$mat_imm_rat <- log2(cat_mark_mc$Mature/cat_mark_mc$Immature)
st_mark_mc[,c('CPN', 'CfuPN', 'Mature', 'Immature', 'Cycling')] <- c()

day_mat <- matrix(0, nrow=length(unique(mc@mc)), ncol = length(unique(mat@cell_metadata[names(mc@mc),'day'])))
colnames(day_mat) <- sort(unique(mat@cell_metadata[names(mc@mc),'day']))
rownames(day_mat) <- 1:nrow(day_mat)
for (i in 1:nrow(day_mat)) {
    tbl <- table(mat@cell_metadata[names(mc@mc)[mc@mc == i], 'day'])
    day_mat[i,names(tbl)] <- tbl
}


max_day <- apply(day_mat, 1, which.max)

cond1 <- apply(st_mark_mc, 1, max) > 2                                   ## filter for high st specificity
cond2 <- cat_mark_mc$mat_imm_rat > 2 & 
                    legc['Mef2c',] >= -10.5 & 
                    legc['Mapt',] >= -10.5                                ## filter for mature mcs
cond3 <- colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] %in% 
                c('NSC', 'IPC', 'Astrocytes', 'Oligodendrocytes')       ## Check neural or glial
cond7 <- cat_mark_mc$cpn_cfu_rat >= log2(2)                              ## filter for CPN-tending
cond8 <- cat_mark_mc$cpn_cfu_rat <= -log2(2)                             ## filter for CfuPN-tending
cond9 <- colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'NSC'
cond10 <- colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)] == 'IPC'
cond11 <- cat_mark_mc$Cycling >= 2
cond12 <- cond2 & (cond7 | cond8)

st <- ifelse(cond1 & cond3, colnames(st_mark_mc)[apply(st_mark_mc, 1, which.max)], NA)

st[is.na(st) & cond9] <- 'NSC'
st[cond10 & cond11] <- 'IPC_cyc'
st[st == 'Astrocytes' & st_mark_mc$Astrocytes < 15] <- 'NSC'

mature_inds <- as.numeric(which((cond1 & cond2) & is.na(st)))
st[mature_inds] <- colnames(st_mark_mc)[apply(st_mark_mc[mature_inds,], 1, which.max)]

neuron_st <- c('CPN_L2-3','CPN_L5_6','SCPN','CthPN')
cpn <- c('CPN_L5_6', 'CPN_L2-3')
cfupn <- c('CthPN', 'SCPN')

devtools::load_all('~/src//metacell.flow')
scdb_flow_init()
mat <- scdb_mat(nm)
mct <- scdb_mctnetwork(nm)
mcf <- scdb_mctnetflow(nm)
net <- mct@network
net$flow <- mcf@edge_flows
net <- net[net$flow > 0,]
net <- net[net$type1 == 'norm_f',]
net <- net[net$type2 != 'sink',]

mcmd <- vroom::vroom('./BonevCollab//mcmd_pl.tsv')
mc_ag <- table(mc@mc,mat@cell_metadata[names(mc@mc),"day"])
mc_ag_n <- mc_ag/rowSums(mc_ag)
mc_ag_c <- t(t(mc_ag)/colSums(mc_ag))
mc_ag_cn <- mc_ag_c/rowSums(mc_ag_c)

ag_cn_com <- apply(mc_ag_cn, 1, function(x) sum(x*(13:18))/sum(x))
nc <- ncol(mc_ag_cn)

mcs_by_st <- lapply(unique(st), function(u) setNames(mcmd$mean_day[which(st == u)], which(st == u)))
names(mcs_by_st) <- unique(st)

st_mc_t <- lapply(mcs_by_st[neuron_st[neuron_st %in% names(mcs_by_st)]], function(x) round(x)-12)
names(st_mc_t) <- names(mcs_by_st[neuron_st[neuron_st %in% names(mcs_by_st)]])



get_st_flows = function(i, md) {
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
}

    
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

cpn_tot = apply(sum_flow_mat[,colnames(sum_flow_mat) %in% cpn], 1, sum)
cfupn_tot = apply(sum_flow_mat[,colnames(sum_flow_mat) %in% cfupn], 1, sum)

icpn_inds = rownames(sum_flow_mat)[cpn_tot > 0 & cfupn_tot == 0 
                                  ] %>% as.numeric %>% unique
icfupn_inds = c(rownames(sum_flow_mat)[cpn_tot == 0 & cfupn_tot > 0]
               ) %>% as.numeric %>% unique
icpn_cfupn_inds = rownames(sum_flow_mat)[cpn_tot > 0 & cfupn_tot > 0 
                                        ] %>% as.numeric %>% unique


st[icpn_inds] = 'iCPN_late'
st[icfupn_inds] = 'iCfuPN'
st[icpn_cfupn_inds] = 'iCPN/CfuPN'

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

st_bad <- c('CPN_L5_6', 'CPN_L2-3')
mcs_that_flow_to_CPN <- which(colSums(do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(x[,st %in% st_bad])))) >= 0.2)
st_bad <- c('SCPN', 'CfuPN')
mcs_that_flow_to_CfuPN <- which(colSums(do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(x[,st %in% st_bad])))) >= 0.2)

st[intersect(mcs_that_flow_to_CPN, mcs_that_flow_to_CfuPN)] <- 'iCPN/CfuPN'


mcmd_new <- vroom::vroom('./BonevCollab//mcmd_pl_cort.tsv')
color_key <- unique(mcmd_new[,c('st', 'color')])
col_annot <- as.data.frame(as.matrix(mcmd_new$st))
colnames(col_annot) <- "st"
col_annot$cl <- mark_fp_km$cluster
ann_colors <- list(st = tibble::deframe(color_key), 
                    cl = setNames(sample(grep('white|gray|grey', colors(), inv=T, v=T), 
                    length(unique(mark_fp_km$cluster))), 1:(length(unique(mark_fp_km$cluster)))))
col_annot$st_new <- ifelse(is.na(st), "None", st)
ann_colors[["st_new"]] <- c(ann_colors[['st']], setNames('seagreen3', "None"))

                           

ppp <- pheatmap(legc[marks,mok_hc_ord], cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors)
save_pheatmap_png(ppp, './figs/new_pl_cort_reannot_flows.png', h = 2000, w = 2000)
mc2 <- mc
mc2@colors <- color_key$color[match(st, color_key$st)]
scdb_add_mc(id = 'pl_cort_test_colors_reannot', mc2)

mct@mc_id <- "pl_cort_test_colors_reannot"

scdb_add_mctnetwork(id = "pl_cort_test_colors_reannot", mct)

mcf@net_id <- "pl_cort_test_colors_reannot"
scdb_add_mctnetflow("pl_cort_test_colors_reannot", mcf)

cust_st_ord <- c('Oligodendrocytes', 'Astrocytes', 'NSC_cyc', 'NSC', 'IPC_cyc', 'IPC', 'IPC_late', 'iCPN_early', 'iCPN_late', 
                'CPN_L2-3', 'CPN_L5_6','iCPN/CfuPN', 'iCfuPN','SCPN', 'CthPN')
cust_st_ord

cust_st_ord_mc <- unlist(lapply(cust_st_ord, function(u) 
    setNames(which(st == u), rep(u, length(which(st == u))))))

# cust_st_ord_mc <- unlist(lapply(cust_st_ord, function(u) setNames(which(mcmd$cell_type == u), rep(u, length(which(mcmd$cell_type == u))))))

mctnetwork_plot_net(mct_id = "pl_cort_test_colors_reannot", mc_ord = cust_st_ord_mc,
                    flow_id = "pl_cort_test_colors_reannot", 
                    fn = "./figs/pl_cort_mct_plot_test_colors_reannot_2.png")

pheatmap(legc[marks,mok_hc_ord], cluster_cols = F, annotation_col = col_annot, annotation_colors = ann_colors)

mcell_mc2d_plot(mc2d_id = 'pl_cort', colors = color_key$color[match(st, color_key$st)], fig_fn = './figs/pl_cort_mc2d_test_reannot.png')



# get_st_flow <- function(sti, st, day_mat) {
#     st_mcs <- which(st == sti)
#     if (length(st_mcs) > 0) {
#         st_p <- day_mat[st_mcs,]
#         days <- which(apply(st_p, 2, sum, na.rm=T)/colSums(day_mat) >= 0.01)
#         days <- days[days != "1"]
#         st_prop_flow_by_day <- lapply(days, function(di) {
#             p <- mc_ag_cn[,di];
#             p <- ifelse(1:length(p) %in% st_mcs, p, 0)
#             p <- p/sum(p)
#             return(mctnetflow_propogate_from_t(mcf = mcf, t =  di, mc_p = p))
#                                                          }
#                                       )
#         names(st_prop_flow_by_day) <- as.numeric(days)
#         return(st_prop_flow_by_day)
#     }
# }

# mcmd_new <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')

# get_mcs_flowing_to_cell_types <- function(cell_types, day_ma, mcmd) {
#     st_flows <- lapply(cell_types, get_st_flow, st, day_mat)
#     names(st_flows) <- cell_types
#     mcs_flowing_to_st_by_day <- lapply(st_flows, function(z) {
#         yi <- z;
#         nyi <- names(z);
#         mapply(FUN = function(yi, nyi) unique(yi[[2]][[as.numeric(nyi)-1]]@i), yi, nyi)
#     })
#     names(mcs_flowing_to_st_by_day) <- cell_types
#     abcd <- lapply(mcs_flowing_to_st_by_day, function(x) unique(unlist(x)))
#     sum_flow_mat <- as.data.frame(sapply(abcd, function(x) {
#         y <- rep(0, nrow(mcmd)); 
#         y[as.numeric(names(table(x)))] <- as.numeric(table(x)); 
#         return(y)
#     }))
#     colnames(sum_flow_mat) <- cell_types
#     return(sum_flow_mat)
# }


# sum_flow_mat_neurons <- get_mcs_flowing_to_cell_types(neuron_st, day_mat, mcmd_new)

# cpn_tot <- apply(sum_flow_mat_neurons[,colnames(sum_flow_mat_neurons) %in% cpn], 1, sum)
# cfupn_tot <- apply(sum_flow_mat_neurons[,colnames(sum_flow_mat_neurons) %in% cfupn], 1, sum)

# icpn_inds <- rownames(sum_flow_mat_neurons)[cpn_tot > 0 & cfupn_tot == 0 & is.na(st)
#                                   ] %>% as.numeric %>% unique
# icfupn_inds <- c(rownames(sum_flow_mat_neurons)[cpn_tot == 0 & cfupn_tot > 0 & is.na(st)]
#                ) %>% as.numeric %>% unique
# icpn_cfupn_inds <- rownames(sum_flow_mat_neurons)[cpn_tot > 0 & cfupn_tot > 0 & is.na(st)
#                                         ] %>% as.numeric %>% unique

# st[icpn_inds] <- 'iCPN_late'
# st[icfupn_inds] <- 'iCfuPN'
# st[icpn_cfupn_inds] <- 'iCPN/CfuPN'
# imm_st <- c('iCPN/CfuPN', 'iCPN_late', 'iCfuPN')

# sum_flow_mat_imm <- get_mcs_flowing_to_cell_types(imm_st, day_mat, mcmd_new)
# sum_flow_mat_all <- cbind(sum_flow_mat_neurons, sum_flow_mat_imm)

# mature_types <- c('NSC', 'IPC', 'IPC_cyc', 'CthPN', 'SCPN', 'CPN_L2-3', 'CPN_L5_6', 'Oligodendrocytes', 'Astrocytes')

# st[st %in% imm_st] <- NA
# icpn_early_inds <- which(!(st %in% mature_types) &
#                     rowSums(sum_flow_mat_all[,cpn]) == 0 & 
#                     rowSums(sum_flow_mat_all[,cfupn]) == 0 & 
#                     sum_flow_mat_all$iCPN_late >0 & 
#                     sum_flow_mat_all$iCfuPN == 0 & 
#                     sum_flow_mat_all$`iCPN/CfuPN` == 0)
# icpn_late_inds <- which(!(st %in% mature_types) &
#                     rowSums(sum_flow_mat_all[,cfupn]) == 0 & 
#                     (rowSums(sum_flow_mat_all[,cpn]) > 0 |
#                     sum_flow_mat_all$iCPN_late > 0) & 
#                     sum_flow_mat_all$iCfuPN == 0 & 
#                     sum_flow_mat_all$`iCPN/CfuPN` == 0)
# icfupn_inds <- which(!(st %in% mature_types) &
#                     rowSums(sum_flow_mat_all[,cpn]) == 0 & 
#                     sum_flow_mat_all$`iCPN/CfuPN` == 0 &
#                     sum_flow_mat_all$iCPN_late  == 0 & 
#                     (rowSums(sum_flow_mat_all[,cfupn]) > 0 |
#                     sum_flow_mat_all$iCfuPN > 0))
# icpn_cfupn_inds <- which(!(st %in% mature_types) &
#                     rowSums(sum_flow_mat_all[,cpn]) == 0 & 
#                     rowSums(sum_flow_mat_all[,cfupn]) == 0 & 
#                     (sum_flow_mat_all$`iCPN/CfuPN` > 0 |
#                     sum_flow_mat_all$iCPN_late > 0 |
#                     sum_flow_mat_all$iCfuPN > 0))
# st[union(icpn_early_inds, which(is.na(st)))] <- 'iCPN_early'
# st[union(icpn_late_inds, which(is.na(st)))] <- 'iCPN_late'
# st[union(icfupn_inds, which(is.na(st)))] <- 'iCfuPN'
# st[union(icpn_cfupn_inds, which(is.na(st)))] <- 'iCPN/CfuPN'

# st[which(st == 'iCPN_early' & rowSums(day_mat[,c('E13','E14','E15')]) > 10)] <- 'iCfuPN'
# st[which(st == 'iCPN_early' & day_mat[,'E15'] > 10)] <- 'iCPN/CfuPN'
# st[which(st == 'iCPN_late' & rowSums(day_mat[,c('E13','E14')]) > 10)] <- 'iCPN/CfuPN'
