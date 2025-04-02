## Boilerplate
## Start
library(metacell)
# devtools::load_all("~/src/metacell.flow")
library(metacell.flow)
library(princurve)
library(vioplot)
# wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex'
wd <- '.'
# setwd(wd)
set.seed(1337)

db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))
source(file.path(wd, 'scripts/util.r'))

nm = 'pl'

mat = scdb_mat(nm)
mc = scdb_mc(nm)

feats = scdb_gset(paste0(nm, '_f'))
feats = names(feats@gene_set)

mg_bon_marks <- as.data.frame(t(sapply(apply(readr::read_csv(file.path(wd, 'input/marker_genes.tsv')), 1, 
                                             stringr::str_split, ' '), function(x) c(x[[1]][[1]], x[[1]][[length(x[[1]])]]))))

colnames(mg_bon_marks) <- c('cell_type', 'marks')

mbm_lst <- lapply(1:nrow(mg_bon_marks), function(n) stringr::str_split(mg_bon_marks$marks[[n]], ',')[[1]])
names(mbm_lst) <- mg_bon_marks$cell_type

legc <- log2(1e-5 + mc@e_gc)

mcmd = vroom::vroom(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )))
col_annot = mcmd[,c('metacell', 'cell_type', 'mean_day')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')

ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

device <- 'pdf'
dir.create('./output/paper_figs/FigS1/')
fig_s1ab_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1AB.{device}'))
fig_s1c_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1C.{device}'))
fig_s1d_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1D.{device}'))
fig_s1e_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1E.{device}'))
fig_s1f_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1F.{device}'))
fig_s1g_path <- file.path(wd, glue::glue('output/paper_figs/FigS1/FigS1G.{device}'))


## Fig S1A
mat_cs <- setNames(Matrix::colSums(mat@mat), colnames(mat@mat))

mat_unique_genes <- setNames(diff(mat@mat@p), colnames(mat@mat))

pdf(fig_s1ab_path, h = 300/71, w = 600/71)
par(las = 2, mar = c(8, 5, 2, 1), mfrow = c(1,2), cex.lab = 1, cex.axis = 1)
boxplot(log2(mat_cs) ~ mat@cell_metadata[names(mat_cs),'batch_set_id'], ylab = '', xlab = '')
title(ylab = 'log2 UMIs', line= 3)
title(xlab = 'Batch', line = 7)

## Fig S1B
boxplot(mat_unique_genes ~ mat@cell_metadata[names(mat_unique_genes),'batch_set_id'], ylab = '', xlab = '')
title(ylab = 'Number of unique genes', line= 4)
title(xlab = 'Batch', line = 7)
dev.off()

## Fig S1C
non_cort_genes <- c('Reln', 'Lhx5','Gad2', 'Sst' , 'Lhx6', 'Nrxn3','Gsx2','Dlx2','Aif1', 'C1qb', 'Hexb', 'Igfbp7')
non_cort_mcs <- sort(unique(unlist(sapply(non_cort_genes, function(g) which(mc@mc_fp[g,] >= 1.5)))))

non_cort_cells <- names(mc@mc)[mc@mc %in% non_cort_mcs]

cells_ignore_new <- union(non_cort_cells, mat@ignore_cells)

ca <- as.data.frame(as.numeric(1:ncol(legc) %in% non_cort_mcs))
colnames(ca) <- 'is_non_cort'
rownames(ca) <- 1:nrow(ca)
ac <- list(is_non_cort = setNames(c('black', 'red'), c(0,1)))

p_filt_non_cort_rna <- pheatmap::pheatmap(legc[c('Fabp7',unique(unlist(mbm_lst))),], annotation_col = ca, annotation_colors = ac, clustering_method = 'ward.D2', fontsize = 14,
                   show_colnames = F, treeheight_col = 0, treeheight_row = 0, col = colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(100))

save_pheatmap_pdf(p_filt_non_cort_rna, fig_s1c_path, h = 600/71, w = 1000/71)


## Fig S1D
mc_cc = get_mc_cc(mat_id = 'pl_cort', mc_id = 'pl_cort', mc2d_id = 'pl_cort_not_cor_cc', plot_mc2d = FALSE)

pdf(fig_s1d_path, h = 350/71, w = 700/71)
par(las = 2, mar = c(12,7,1,1), cex.axis= 1.5, cex.lab = 2)
boxplot(100 - mc_cc$cc_score ~ factor(mcmd$cell_type, levels = cust_st_ord), col = col_key[cust_st_ord], ylab = '', xlab = '')
title(xlab= 'Cell type', line = 9)
title(ylab= 'Cell cycle score', line = 4)
dev.off()

## Fig S1E
ct_day_mat <- tgstat::tgs_matrix_tapply(t(mcmd[,grep('E\\d\\d', colnames(mcmd))]), mcmd$cell_type, sum)

ct_day_mat_norm <- t(t(ct_day_mat)/colSums(ct_day_mat))

pdf(fig_s1e_path, h = 1000/71, w = 600/71)
par(cex.axis = 2.5, cex.main = 3, mar = c(5,5,4,1))
barplot(ct_day_mat_norm[cust_st_ord,], col = color_key$color[match(cust_st_ord, color_key$cell_type)], yaxt = 'n', main = 'Cell type fraction by time point')
dev.off()


## Fig S1F
mcf <- scdb_mctnetflow('pl_cort')
ct_flow_out_ls <- lapply(cust_st_ord, function(cti) do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgstat::tgs_matrix_tapply(x[which(mcmd$cell_type == cti),], mcmd$cell_type, sum, na.rm = T), na.rm = T)/sum(colSums(x, na.rm = T), na.rm = T))))
ct_flow_out_ls_n <- lapply(ct_flow_out_ls, function(x) {y <- x/rowSums(x, na.rm = T); y[is.na(y)] <- 0; return(y)})
names(ct_flow_out_ls_n) <- cust_st_ord

ct_flow_in_ls <- lapply(cust_st_ord, function(cti) do.call('rbind', lapply(mcf@mc_backward, function(x) rowSums(tgstat::tgs_matrix_tapply(t(x[,which(mcmd$cell_type == cti)]), mcmd$cell_type, sum, na.rm = T), na.rm = T)/sum(rowSums(x, na.rm = T), na.rm = T))))
ct_flow_in_ls_n <- lapply(ct_flow_in_ls, function(x) {y <- x/rowSums(x, na.rm = T); y[is.na(y)] <- 0; return(y)})
names(ct_flow_in_ls_n) <- cust_st_ord

aa <- paste0('E', 13:17)
bb <- paste0('E', 14:18)
xlabs <- apply(cbind(aa, rep('->', length(aa)), bb), 1, paste, collapse = ' ')

pdf(fig_s1f_path, h = 2600/71, w = 700/71)
par(mfrow = c(length(cust_st_ord),2), las = 2, mar = c(1,6,2,1), cex.main = 2, cex.axis = 2, cex.lab = 2)
for (ct in cust_st_ord) {
    i <- match(ct, names(ct_flow_out_ls_n))
    xi <- ct_flow_in_ls_n[[i]]
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,6,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, xlab = '', ylab = '',
         main = glue::glue('{ct} incoming flows'),  xaxt = 'n')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Relative flow', line = 4, cex.lab = 2)
    dummy <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
    xi <- ct_flow_out_ls_n[[i]]
    if (ct == tail(cust_st_ord, 1)) {par(mar = c(10,6,2,1))}
    plot(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], 
         type = 'l', ylim = c(0,max(xi)), lwd = 3, 
         main = glue::glue('{ct} outgoing flows'), ylab = '', xaxt = 'n', xlab = '')
    if (ct == tail(cust_st_ord, 1)) {axis(1, at = 1:5, labels = xlabs)}
    title(ylab = 'Relative flow', line = 4, cex.lab = 2)
    dummy <- sapply(colnames(xi), function(ct) {
            lines(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
            points(1:5, xi[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
    })
}
dev.off()
