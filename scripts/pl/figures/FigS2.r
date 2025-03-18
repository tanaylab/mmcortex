library(metacell)
devtools::load_all('~/src/metacell.flow')
library(ComplexHeatmap)
library(matrixStats)
wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
db_path <- file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scdb_flow_init()
SEED <- 1337
K <- 16
set.seed(SEED)
scfigs_init("figs/")
doMC::registerDoMC(60)
nm <- 'pl_cort'
source(file.path(wd,'scripts/util.r'))

mc <- scdb_mc(nm)
mat <- scdb_mat(nm)

mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
color_key <- unique(mcmd[,c('cell_type', 'color')])

cust_st_ord <- c('OPCs', 'Astrocytes', 'NSC', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )))

cust_st_ord2 <- c('OPCs', 'Astrocytes', 'NSC', 'IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st2 <- unlist(lapply(cust_st_ord2, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         
goi <- c('Pou3f1', 'Pou3f2', 'Cux1', 'Cux2', 'Neurod1', 'Neurog2', 'Id4',
         'Eomes', 'Hes1', 'Apoe', 'Sox5', 'Tbr1', 'Foxp2', 'Foxp1', 'Nfia', 'Islr2', 
         'Zbtb20', 'Bcl11b', 'Fezf2', 'Satb2', 'Mef2c', 'Nhlh1', 'Tle4',
        'Rnd2',  'Runx1t1', 'Mapt', 'Mki67', 'Pcna',
        'Fabp7', 'Olig1', 'Ldb2', 'Gadd45g', 'Syt4')
marks_filt <- goi

m_genes <- c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes <- c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")
cc_genes <- union(m_genes, s_genes)

col_annot <- mcmd[,c('metacell', 'cell_type', 'mean_day')]
col_annot <- tibble::column_to_rownames(col_annot, 'metacell')

clrmp <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue1', 'blue4', 'purple3'))(1000)

clrmp_abs <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(1000)
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)

ann_colors <- list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

legc <- log2(1e-05 + mc@e_gc)


dir.create('./output/paper_figs/FigS2/')
device <- 'pdf'
fig_s2a_path <- glue::glue('./output/paper_figs/FigS2/FigS2A.{device}')
fig_s2b_path <- glue::glue('./output/paper_figs/FigS2/FigS2B.{device}')
fig_s2b_legend_path <- glue::glue('./output/paper_figs/FigS2/FigS2B_legend.{device}')
fig_s2c_path <- glue::glue('./output/paper_figs/FigS2/FigS2C.{device}')
fig_s2d_path <- glue::glue('./output/paper_figs/FigS2/FigS2D.{device}')
fig_s2e_path <- glue::glue('./output/paper_figs/FigS2/FigS2E.{device}')
fig_s2f_path <- glue::glue('./output/paper_figs/FigS2/FigS2F.{device}')
fig_s2g_path <- glue::glue('./output/paper_figs/FigS2/FigS2G.{device}')
fig_s2h_path <- glue::glue('./output/paper_figs/FigS2/FigS2H.{device}')
fig_s2i_path <- glue::glue('./output/paper_figs/FigS2/FigS2I.{device}')
fig_s2j_path <- glue::glue('./output/paper_figs/FigS2/FigS2J.{device}')

astro_module <- readLines('./output/metacell_model/nsc_gene_modules/astro_module.txt')
ipc_module <- readLines('./output/metacell_model/nsc_gene_modules/ipc_module.txt')
stem_module <- readLines('./output/metacell_model/nsc_gene_modules/stem_module.txt')

load('./output/metacell_model/nsc_gene_modules/figs2_data.rda')
load('./output/metacell_model/nsc_gene_modules/phase_info.rda')

## Fig S2A
st_legc <- as.data.frame(t(tgs_matrix_tapply(legc, mcmd$cell_type, mean)))


cluster_names <- setNames(c('Cell cycle 1', 
                            'Temp. decreasing 1', 
                            'Cell cycle 2',
                            'Temp. increasing 1',
                            'Temp. increasing 2',
                            'Cell cycle 3',
                            'Cell cycle 4',
                            'Temp. decreasing 2'), sort(unique(ct_hc_cor_nsc)))

## Gene module table for MCV and supp table 1
all_genes_in_modules <- multunion(names(ct_hc_cor_nsc), astro_module, ipc_module, stem_module)
gene_module_table <- tibble::enframe(ct_hc_cor_nsc, name = 'gene', value = 'nsc_gene_module') 
gene_module_table$nsc_gene_module_name <- cluster_names[gene_module_table$nsc_gene_module]
gene_module_table[,c('IPC', 'astro', 'stem')] <- cbind(ifelse(gene_module_table$gene %in% ipc_module, TRUE, FALSE), 
                                                        ifelse(gene_module_table$gene %in% astro_module, TRUE, FALSE), 
                                                        ifelse(gene_module_table$gene %in% stem_module, TRUE, FALSE))
gene_module_table[,colnames(legc_by_day_n)] <- legc_by_day_n[gene_module_table$gene,]
gene_module_table <- gene_module_table %>% dplyr::arrange(nsc_gene_module_name, gene)

readr::write_tsv(gene_module_table, './output/metacell_model/nsc_gene_modules/supp_table_1_nsc_gene_modules.tsv')


pdf(fig_s2a_path, h = 500/71, w = 1000/71)


EXPAND_FACTOR <- 3
RATIO <- 1.3
layout_mat = matrix(c(rep(1:4, EXPAND_FACTOR), rep(5:8, round(EXPAND_FACTOR*RATIO))),
                nrow = EXPAND_FACTOR + round(EXPAND_FACTOR*RATIO),
                ncol = 4, 
                byrow = T)
layout(layout_mat)
mari <- c(9,6,3,0.5)
par(las = 2, cex.main = 2, cex.lab = 2, cex.axis = 1.52, mar = mari)
vvv <- lapply(sort(cluster_names), function(cnj) {
    clj <- as.numeric(names(cluster_names)[cluster_names == cnj])
    gnj <- names(ct_hc_cor_nsc)[ct_hc_cor_nsc == clj]
    if (grepl('cell', cnj, ign = T)) {
        xaxti <- 'n'
        mari[[1]] <- 0.5
    } else {
        xaxti <- 's'
        mari[[1]] <- 10
    }
    par(mar = mari)
    boxplot(st_legc[gnj,cust_st_ord2], col = col_key[cust_st_ord2], 
            main = cluster_names[[clj]], 
            ylab = '',
            xaxt = xaxti,
            ylim = quantile(unlist(st_legc[gnj,cust_st_ord2]), c(0.1,0.96)))
    title(ylab = 'Mean RNA', line = 4)
})
dev.off()

cluster_names <- setNames(gsub(' ', '\n', cluster_names), names(cluster_names))


## Fig S2B
ct_hc_cor_nsc_h <- setNames(cluster_names[ct_hc_cor_nsc], names(ct_hc_cor_nsc))

ca <- columnAnnotation(df = tibble::column_to_rownames(tibble::enframe(ct_hc_cor_nsc_h[hc_cor_nsc$order], name = 'gene', value = 'cluster'), 'gene'), 
                       show_legend = c('cluster' = F),
                       col =  list(cluster = setNames(chameleon::distinct_colors(8)$name, cluster_names[1:8])))

ra <- rowAnnotation(df = tibble::column_to_rownames(tibble::enframe(ct_hc_cor_nsc_h[hc_cor_nsc$order], name = 'gene', value = 'cluster'), 'gene'),
                                           show_legend = c('cluster' = F),
                   col =  list(cluster = setNames(chameleon::distinct_colors(8)$name, cluster_names[1:8])))

ac <- list(cluster = setNames(chameleon::distinct_colors(8)$name, 1:8))

ch_cor_dyn_genes_nsc <- ComplexHeatmap::Heatmap(cor_nsc_legc_dyn_genes[hc_cor_nsc$order,hc_cor_nsc$order], name = ' ',
                                                column_split = ct_hc_cor_nsc_h[hc_cor_nsc$order],
                                                row_split = ct_hc_cor_nsc_h[hc_cor_nsc$order],
                                                top_annotation = ca, left_annotation = ra,
                                                show_row_names = F, show_column_names = F, 
                                                col = circlize::colorRamp2(colors = c('blue3', 'white', 'red3'), breaks = c(-1,0,1)),
                   cluster_columns = F, cluster_rows = F)

pdf(fig_s2b_path, h = 10, w = 10)
draw(ch_cor_dyn_genes_nsc)
dev.off()


## Fig S2C

tbl_pba_by_ct <- t(table(sc_data_df$cell_type, sc_data_df$pba))
tbl_pba_by_ct_norm <- t(t(tbl_pba_by_ct)/colSums(tbl_pba_by_ct))

rownames(tbl_pba_by_ct_norm) <- gsub('\\d_', '', rownames(tbl_pba_by_ct_norm))

p_pba_by_ct <- pheatmap::pheatmap(tbl_pba_by_ct_norm[,cust_st_ord], cluster_cols = F, 
                    annotation_legend = F, col = clrmp_abs, fontsize = 12,
                   treeheight_row = 0, treeheight_col = 0)
save_pheatmap_pdf(p_pba_by_ct, fig_s2c_path, h = 500/71, w = 800/71)

## Fig S2D
nsc_sc <- intersect(names(mc@mc[mc@mc %in% which(mcmd$cell_type == 'NSC')]), colnames(mat_ds))

nsc_sc_by_day <- lapply(tail(sort(unique(mat@cell_metadata$day)),-1), function(di) intersect(nsc_sc, rownames(mat@cell_metadata)[mat@cell_metadata$day == di]))
names(nsc_sc_by_day) <- tail(sort(unique(mat@cell_metadata$day)),-1)

nsc_sc_by_day_vec <- setNames(unlist(sapply(names(nsc_sc_by_day), function(x) rep(x, length(nsc_sc_by_day[[x]])))), do.call('c', nsc_sc_by_day))


pdf(fig_s2d_path, h = 7, w = 28)
par(mfrow = c(1,4), mar = c(6,7,4,1), cex.lab = 3, cex.axis = 3, cex.main = 5)

NUM_PARTITION <- 13
phase_qs <- seq(1-1e-2,max(phase),l=NUM_PARTITION)
bin_borders <- setNames(phase_qs[c(1,2,4,5,8,10,12, 13)], c('G1', 'G1_0', 'G1', 'S', 'G2', 'M', 'G1'))
heights_text <- c(1.5, 2.2, 2.2, 4.2)
ttt <- sapply(1:nrow(mat_ds_cc_genes_select_ord_phase), function(i) {
    plot(sort(phase[nsc_sc]), mat_ds_cc_genes_select_ord_phase[i,], 
         col = 'white', 
         cex = 1, pch = 16, ylim = quantile(mat_ds_cc_genes_select_ord_phase[i,], c(0.05, 0.95)),
        xlab = '', ylab = ''
        )
    vvv <- sapply(head(seq_along(bin_borders), -1), function(j) {
        bj <- bin_borders[[j]]
        lines(rep(bj, 2), c(-1,10), lty = 2, lwd = 2)
        nmj <- names(bin_borders)[[j]]
    })
    lines(sort(phase[nsc_sc]), mat_ds_cc_genes_select_ord_phase_rm[i,], col = 'black', lwd = 3)
    if (i == 1) {legend('topleft', legend = glue::glue('Rollmean k = {K}'), lwd = 3, col = 'black', cex = 2, bg = 'white')}
    title(main = rownames(mat_ds_cc_genes_select_ord_phase)[[i]])
    title(xlab = 'phase', line = 4)
    title(ylab = 'Downsampled UMIs', line = 4)
})
dev.off()


## Fig S2E
s_genes_sum <- Matrix::colSums(mat_ds[s_genes,])

m_genes_sum <- Matrix::colSums(mat_ds[m_genes,])


pdf(fig_s2e_path, h = 500/71, w = 1500/71)
par(mfrow = c(1,3), cex.main = 2, cex.axis = 2, cex.lab = 2, mar = c(6,6,5,1), las = 2)
boxplot(s_genes_sum[names(nsc_sc_by_day_vec)]/length(s_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], 
            main = 'UMIs per S gene per single NSC per bin\nn_{S genes} = 11', ylab = '', xlab = '')
title(ylab = 'UMIs per gene', line = 4)
boxplot(m_genes_sum[names(nsc_sc_by_day_vec)]/length(m_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], 
            main = 'UMIs per M gene per single NSC per bin\nn_{M genes} = 26', ylab = 'UMIs per gene', xlab = '')
boxplot(s_genes_sum[names(nsc_sc_by_day_vec)]/length(s_genes) + m_genes_sum[names(nsc_sc_by_day_vec)]/length(m_genes) ~ phase_cut[names(nsc_sc_by_day_vec)], 
        main = 'UMIs per S+M gene per cell per bin', ylab = 'UMIs per gene', xlab = '')
dev.off()

## Fig S2F
days <- unique(mat@cell_metadata[nsc_sc,'day'])

nsc_inds <- which(sc_data_df$cell_type == 'NSC')
lupc <- length(unique(phase_cut))
bin_seq <- seq(lupc+0.5, 0.5+length(days)*lupc, 12)
pdf(fig_s2f_path, h = 1050/71, w = 1600/71)
par(mfrow = c(3,1), las = 2, cex.lab = 3, cex.axis = 1.5, mar = c(6,6,0.5,0.5))
boxplot(ipc ~ ., data = sc_data_df[nsc_inds,c('ipc','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$ipc[nsc_inds], 0.94, na.rm = T)),
       col = rep(rainbow(lupc), length(days)), ylab = 'IPC module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))
boxplot(astro ~ ., data = sc_data_df[nsc_inds,c('astro','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$astro[nsc_inds], 0.975, na.rm = T)),
              col = rep(rainbow(lupc), length(days)), ylab = 'Astro module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))
boxplot(stem ~ ., data = sc_data_df[nsc_inds,c('stem','phase_cut', 'day')], ylim = c(0,quantile(sc_data_df$stem[nsc_inds], 0.95, na.rm = T)),
              col = rep(rainbow(lupc), length(days)), ylab = 'Stem module ds UMIs', xlab = '')
vvv <- sapply(bin_seq, function(i) lines(rep(i,2), c(0,1000), lwd = 3))
dev.off()



## Fig S2G
nsc_late <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day > 16.5])
nsc_early <- as.character(mcmd$metacell[mcmd$cell_type == 'NSC' & mcmd$mean_day < 14.5])
nsc_late_rna <- rowMeans(legc[,as.numeric(nsc_late)])
nsc_early_rna <- rowMeans(legc[,as.numeric(nsc_early)])
astro_rna <- rowMeans(legc[,mcmd$metacell[mcmd$cell_type == 'Astrocytes']])


pdf(fig_s2g_path, h = 500/71, w = 1000/71)
par(mfrow = c(1,2), cex.lab = 1.52, cex.main = 1.5)
plot(nsc_late_rna, astro_rna, pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) RNA', ylab = 'Astrocytes RNA', main = 'Astro vs late NSC - RNA')
abline(a =-2,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+2,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-1,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '1 LFC','2 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(nsc_late_rna, astro_rna, method = 'pearson')
text(-13.5,-6, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
plot(nsc_late_rna, nsc_early_rna, pch = 16, cex = .15, xlab = 'NSC (mean day > 16.5) RNA', ylab = 'NSC (mean day < 14.5) RNA', main = 'Early vs late NSC - RNA')
abline(a =-2,b = 1,col='red', lty=  2, lwd=  1)
abline(a =+2,b = 1,col='red', lty=  2, lwd=  1)
abline(a =-1,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+1,b = 1,col='green', lty=  2, lwd=  1)
abline(a =+0,b = 1,col='blue', lty=  2, lwd=  1)
legend('bottomright', legend = c('0 LFC', '1 LFC','2 LFC'), col = c('blue', 'green', 'red'), lty = 2, lwd =1)
corh <- cor(nsc_late_rna, nsc_early_rna, method = 'pearson')
text(-13.5,-9, labels = paste0('R^2 = ', signif(corh**2, 2)), cex = 1.5)
dev.off()


## Fig S2H
### Plot IPC vs stem and astro vs stem signatures

pdf(fig_s2h_path, w = 400/71, h = 400/71)
par(mfrow = c(1,1), mar = c(5,5,1,1), cex.lab = 1, cex.axis = 1)
plot(sc_data_df[names(matched_coords_bins_astro),'astro'], sc_data_df[names(matched_coords_bins_astro),'stem'], col = mcmd$color[mc@mc[names(matched_coords_bins_astro)]], ylab = 'Stem module UMIs', xlab = 'Astro module UMIs', pch =16, cex = 0.5)
dev.off()


## Fig S2I
astro_mcs <- which(mcmd$cell_type == 'Astrocytes')

oligo_mcs <- which(mcmd$cell_type == 'OPCs')

cpnl23_mcs <- which(mcmd$cell_type == 'CPN_L2-3')

gh <- names(ct_hc_cor_nsc[ct_hc_cor_nsc == 4])
pltmt <- cbind(legc_by_day_n[gh,], rowMeans(legc[gh,astro_mcs]), 
            rowMeans(legc[gh,oligo_mcs]), rowMeans(legc[gh,cpnl23_mcs]))
colnames(pltmt)[(ncol(pltmt)-2):ncol(pltmt)] <- c('Astrocytes', 'OPCs', 'CPN_L2-3')
genes_hi <- rownames(pltmt)[which(rowMaxs(subset(pltmt, select = -c(Astrocytes, OPCs, `CPN_L2-3`))) > pltmt[,'Astrocytes'])]

genes_norm <- setdiff(rownames(pltmt), genes_hi)

pltmt2 <- rbind(pltmt[genes_hi[hclust(dist(pltmt[genes_hi,] - pltmt[genes_hi,'Astrocytes']), method = 'ward.D2')$order],],
                pltmt[genes_norm[hclust(dist(pltmt[genes_norm,] - pltmt[genes_norm,'Astrocytes']), method = 'ward.D2')$order],])

p_nsc_cl4 <- pheatmap::pheatmap(pltmt2 - pltmt2[,'Astrocytes'], 
                            gaps_row = length(genes_hi),
                                cluster_rows = F,
                                cluster_cols = F, 
                                col = clrmp_rel, breaks = seq(-2,2,l=1000), 
                                clustering_method = 'ward.D2', treeheight_row = 0, fontsize_row = 4)



save_pheatmap_pdf(p_nsc_cl4, fig_s2i_path, h= 1200/71, w = 600/71)
