## Boilerplate
## Start
library(metacell)
devtools::load_all('~/src/metacell.flow')
library(ComplexHeatmap)

wd = '/home/feshap/raid/proj/mmcortex'
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scdb_flow_init()
SEED = 1337
K = 16
set.seed(SEED)
scfigs_init("figs/")
doMC::registerDoMC(60)
nm = 'pl_cort'
mc2d_id = 'pl_cort_not_cor_cc'
source('./scripts/util.r')


mc = scdb_mc(nm)
mat = scdb_mat(nm)
mct = scdb_mctnetwork(nm)
mcf = scdb_mctnetflow(nm)
mgraph <- scdb_mgraph(nm)
mc2d <- scdb_mc2d(mc2d_id)

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )))

cust_st_ord2 = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st2 = unlist(lapply(cust_st_ord2, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         
goi = c('Pou3f1', 'Pou3f2', 'Cux1', 'Cux2', 'Neurod1', 'Neurog2', 'Id4',
         'Eomes', 'Hes1', 'Apoe', 'Sox5', 'Tbr1', 'Foxp2', 'Foxp1', 'Nfia', 'Islr2', 
         'Zbtb20', 'Bcl11b', 'Fezf2', 'Satb2', 'Mef2c', 'Nhlh1', 'Tle4',
        'Rnd2',  'Runx1t1', 'Mapt', 'Mki67', 'Pcna',
        'Fabp7', 'Olig1', 'Ldb2', 'Gadd45g', 'Syt4')
marks_filt = goi

m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")
cc_genes <- union(m_genes, s_genes)

col_annot = mcmd[,c('metacell', 'cell_type', 'mean_day')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')

clrmp <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue1', 'blue4', 'purple3'))(1000)

clrmp_abs <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(1000)
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)

ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

legc = log2(1e-05 + mc@e_gc)
## End

fig_2a_path <- './output/paper_figs/Fig1/Fig1A.png'
fig_2a_legend_path <- './output/paper_figs/Fig1/Fig1A_legend.png'
fig_2b_path <- './output/paper_figs/Fig1/Fig1B.png'
fig_2c_path <- './output/paper_figs/Fig1/Fig1C.png'
fig_2c_legend_path <- './output/paper_figs/Fig1/Fig1C_legend.png'
fig_2d_path <- './output/paper_figs/Fig1/Fig1D.png'
fig_2d_color_bar_path <- './output/paper_figs/Fig1/Fig1D_color_bar.png'
fig_2e_path <- './output/paper_figs/Fig1/Fig1E.png'
fig_2f_path <- './output/paper_figs/Fig1/Fig1F.png'


## Fig 2A
mcs_here <- union(nsc_mcs, which(mcmd$cell_type %in% c('IPC', 'IPC_cyc', 'Astrocytes')))
legc_mm_nsc <- legc - rowMeans(legc[,nsc_mcs])
mat_ipc_module <- legc_mm_nsc[ipc_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]]
mat_astro_module <- legc_mm_nsc[astro_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]]
mat_stem_module <- legc_mm_nsc[stem_module,cust_mc_ord_st[cust_mc_ord_st %in% mcs_here]] 
hc_mat_ipc <- hclust(dist(mat_ipc_module), method = 'ward.D2')
hc_mat_astro <- hclust(dist(mat_astro_module), method = 'ward.D2')
hc_mat_stem <- hclust(dist(mat_stem_module), method = 'ward.D2')

pltmt_modules <- rbind(mat_ipc_module[hc_mat_ipc$order,], mat_astro_module[hc_mat_astro$order,], mat_stem_module[hc_mat_stem$order,])

marks <- names(scdb_gset('pl_cort_marks_f')@gene_set)

ipc_genes_not_in_goi <- setdiff(ipc_module, goi)
lr <- ifelse(rownames(pltmt_modules) %in% setdiff(marks, ipc_genes_not_in_goi), rownames(pltmt_modules), '')

p_pltmt_modules <-pheatmap::pheatmap(pltmt_modules, silent = T,
                                     cluster_cols = F, 
                                     fontsize_row = 14,
                                     annotation_col = col_annot, cluster_rows = F, 

                                     gaps_col = match(unique(names(cust_mc_ord_st)[cust_mc_ord_st %in% mcs_here]), names(cust_mc_ord_st)[cust_mc_ord_st %in% mcs_here]) -1 ,
                                     annotation_legend = F, annotation_colors = ann_colors, 
                                     col = colorRampPalette(c('blue3', 'white', 'red3'))(100), breaks = seq(-3,3,l=100), 
                                     treeheight_row = 0, show_colnames = F)



p_pltmt_modules$gtable <- add.flag(pheatmap = p_pltmt_modules, kept.labels = lr, repel.degree = 0)

mcATAC::save_pheatmap(x = p_pltmt_modules,filename =  './output/metacell_model/nsc_gene_modules/figs/nsc_ipc_astro_stem_modules_phm_flag.png', h = 900, w = 850)


## Fig 2B
astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')
ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')
stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')

set.seed(1337)
mat_ds <- scm_downsamp(mat@mat, 3000)
nsc_sc <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'NSC')]), colnames(mat_ds))
ipc_sc <- colSums(mat_ds[ipc_module,])
astro_sc <- colSums(mat_ds[astro_module,])
stem_sc <- colSums(mat_ds[stem_module,])
ipc_sc_names <- names(mc@mc)[mc@mc %in% which(mcmd$cell_type %in% c('IPC'))]
astro_sc_names <- names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'Astrocytes')])
days <- paste0('E', 13:18)


mari <- c(5,5,2,2)
BOXWEX <- 0.6
LINE <- 5
png(fig_2b_path, h = 1100, w = 1200)
par(mfcol = c(3,3), mar = mari, cex.main = 2, cex.lab = 2, cex.axis = 2, las = 2)
mari[[2]] <- 7
mari[[3]] <- 5
par(mar = mari)
boxplot(ipc_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX, xaxt = 'n',ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', main = 'NSCs', xlab = '', col = col_key[['NSC']])
title(ylab = "IPC module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX,xaxt = 'n',ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), ylab = '', xlab = '', col = col_key[['NSC']])
title(ylab = "Astro module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)
boxplot(stem_sc[nsc_sc] ~ setNames(mat@cell_metadata[nsc_sc,'day'], nsc_sc), boxwex = BOXWEX,xaxt = 'n',ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', xlab = '', col = col_key[['NSC']])
title(ylab = "Stem module (ds UMIs)", line = LINE)
axis(1, at = 1:6, labels = days)

mari[[2]] <- 3
mari[[3]] <- 5
par(mar = mari)
boxplot(ipc_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5), ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', main = 'Astrocytes', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5),ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), ylab = '', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))
boxplot(stem_sc[astro_sc_names] ~ setNames(mat@cell_metadata[astro_sc_names,'day'], astro_sc_names), boxwex = BOXWEX, xlim = c(-2.5,3.5),ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), ylab = '', xlab = '', col = col_key[['Astrocytes']])
axis(1, at = seq(-2,3,1), labels = paste0('E', 13:18))

mari[[3]] <- 5
mari[[2]] <- 3
par(mar = mari)
boxplot(ipc_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(ipc_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), main = 'IPCs',ylab = '', xlab = '', col = col_key[['IPC']])
mari[[3]] <- 2
par(mar = mari)
boxplot(astro_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(astro_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], .999, na.rm = T)), xlab = '', ylab = '',col = col_key[['IPC']])
boxplot(stem_sc[ipc_sc_names] ~ setNames(mat@cell_metadata[ipc_sc_names,'day'], ipc_sc_names), boxwex = BOXWEX, ylim = c(0,quantile(stem_sc[c(ipc_sc_names, nsc_sc, astro_sc_names)], 0.99, na.rm = T)), xlab = '',ylab = '', col = col_key[['IPC']])

dev.off()


## Fig 2C



## Fig 2D

## Fig 2E

## Fig 2F