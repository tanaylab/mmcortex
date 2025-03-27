# wd <- '/net//mraid20//export/tgdata/users/yonshap/proj/mmcortex/'
wd <- '.'
setwd(wd)
library(shaman)
library(viridis)
library(cowplot)
library(stringr)
library(patchwork)
library(ggrastr)
library(ggplot2)
library(ggplotify)
library(rtracklayer)
library(Gviz)
genome <- 'mm10'
trackdb <- genome
misha::gdb.create_genome(genome)
misha::gsetroot(genome)
source(file.path(wd, 'scripts/util.r'))
theme_set(theme_cowplot())


load(file = file.path(wd, 'output/methylation/avg_meth_all.rda'))

load(file = file.path(wd, 'output/hic/fig5_data.rda'))

mcmd <- readr::read_tsv(file.path(wd, 'output/metacell_model/mcmd_pl_cort.tsv'))

color_key <- unique(mcmd[,c('cell_type', 'color')])
col_key <- tibble::deframe(color_key)

cust_st_ord <- c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st <- unlist(lapply(cust_st_ord, function(st) setNames(as.character(mcmd$metacell[which(mcmd$cell_type == st)[order(mcmd$mean_day[which(mcmd$cell_type == st)])]]), 
                                                                  rep(st, length(which(mcmd$cell_type == st))))))



palette.breaks = function(n, colors, breaks){
  colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
  
  for(i in 2:(length(colors)) ){
    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
  }
  colspec = c( colspec,
               colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
  )
  colspec
}

col.pbreaks <<- c(20,35,50,65,75,85,95)        #Original
col.pos <<- palette.breaks(200 , c("lightgrey","lavenderblush2",
                                "#f8bfda","lightcoral",
                                "red","orange","yellow"), col.pbreaks)
col.nbreaks <<- c(20,35,50,65,75,85,95)
col.neg <<- rev(palette.breaks(100 , c("lightgrey", 
                                    "powderblue", "cornflowerblue", 
                                    "blue","blueviolet", "#8A2BE2", 
                                    "#4B0082"), col.nbreaks ))
col.scores <<- c(col.neg, "lightgrey", col.pos)



dir.create(file.path(wd, 'output/paper_figs/Fig5/'))
device <- 'pdf'
fig_5a_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5A.{device}'))
fig_5b_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5B.{device}'))
fig_5c_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5C.{device}'))
fig_5d_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5D.{device}'))
fig_5e_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5E.{device}'))
fig_5f_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5F.{device}'))
fig_5g_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5G.{device}'))
fig_5h_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5H.{device}'))
fig_5i_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/Fig5I.{device}'))
shaman_color_bar_path <- file.path(wd, glue::glue('output/paper_figs/Fig5/shaman_color_bar.{device}'))

## Fig 5A

p_tad_borders_all <- pheatmap::pheatmap(pltmt - rowMeans(pltmt, na.rm = T), gaps_row = nrow(hvmm_inc),
                                        cluster_rows = F, cluster_cols = F, 
                                        col = colorRampPalette(c('blue3', 'white', 'red3'))(100), 
                                        show_rownames = F, fontsize_col =  20,
                   breaks = seq(-0.35,0.35,l=100))
save_pheatmap_pdf(p_tad_borders_all, fig_5a_path, h = 800/100, w = 450/100)

clrmp_trks <- colorRampPalette(c('yellow', 'blue'))(5)

plot_insulation_and_shaman_series <- function(i, hvmi, file_path) {
    inti_nm <- paste0(hvm[hvmi,1:3], collapse = '_')
    
    # tadi_fld <- fig_5c_path
    # if (!dir.exists(tadi_fld)) {dir.create(tadi_fld)}
    SHIFT <- 50e+4
    gints1d_p <- dplyr::mutate(hvm[hvmi,], start = start -SHIFT, end = end + SHIFT)
    gints2d_p <- dplyr::mutate(hvm_2d[hvmi,], start1 = start1 -SHIFT, end1 = end1 + SHIFT,
                              start2 = start2 -SHIFT, end2 = end2 + SHIFT)

    nei_gints1d_ins_prc_all <- misha::gintervals.neighbors(gints1d_p, ins_prc_all, mindist = 0, maxdist = 0, maxneighbors = 1e+7)
    inds <- nei_gints1d_ins_prc_all$intervalID
    
    # pdf(glue::glue('{tadi_fld}/Fig5C{i}a.pdf'), h = 750/71, w = 2000/71)
    pdf(file_path, h = 750/71, w = 2000/71)
    par(mar = c(7,10,1,1), cex.lab = 5, cex.main = 5, cex.axis = 5)
    yall <- -1 * ins_prc_all[inds,grep('E\\d\\d', colnames(ins_prc_all))]
    plot(ins_prc_all[inds,2], -1 * ins_prc_all[inds,grep('E13', colnames(ins_prc_all))], ylim = c(min(yall, na.rm = T), max(yall, na.rm = T)), bg = 'white',
         type = 'l', 
         lwd = 6,
         xaxt = 'n',
         col = clrmp_trks[[1]], pch = 16, cex = 1.5, 
         xlab = '',
         main = '',
         ylab = ''
         )
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E14', colnames(ins_prc_all))], col = clrmp_trks[[2]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E15', colnames(ins_prc_all))], col = clrmp_trks[[3]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E16', colnames(ins_prc_all))], col = clrmp_trks[[4]], pch = 16, lwd = 6)
    lines(ins_prc_all[inds,2], -1*ins_prc_all[inds,grep('E17', colnames(ins_prc_all))], col = clrmp_trks[[5]], pch = 16, lwd = 6)
    title(xlab = unique(as.character(nei_gints1d_ins_prc_all$chrom)), line = 6)
    title(ylab = 'Insulation', line = 5)
    axis(1, padj = 0.5, hadj = -0.25)
    legend('topleft', cex = 4, col = clrmp_trks, legend = paste0('E', 13:17), lwd = 6)

    dev.off()

}

hvm_ls <- c('86', '214')

## Fig 5B -- Hapln1
plot_insulation_and_shaman_series(1, hvm_ls[[1]], file_path = fig_5b_path)

## Fig5C - from Boyan
p_hapln1 <- plotRegion(interv = interv_hapln1,genome = genome, hic_data = hic_data_hapln1, 
            plotHicTracks, plotGeneTrack, trackdb, scoreTracksToExtract, point_cex = 0.75)
pdf(fig_5c_path, h = 600/71, w = 510/71)
print(p_hapln1)
dev.off()

## Shaman legend
shaman_colors <- shaman_score_pal()

plot_color_bar(vals = seq(-100,100,l=length(shaman_colors)), 
            cols = shaman_colors, height = 500/71, width = 250/71, device = pdf,
            show_vals_ind = seq(1,length(shaman_colors),l=11), fig_fn = shaman_color_bar_path)


## Fig 5D
nsc_glia_mcs <- mcmd$metacell[mcmd$cell_type %in% c('NSC', 'Astrocytes', 'OPCs')]

g <- 'Hapln1'
MAR <- c(6,6,1,1)
CEX <- 2
pdf(fig_5d_path, h = 350/71, w = 550/71)
par(cex.lab = 2, cex.axis = 1.5, mar = MAR)
plot(mcmd$mean_day[nsc_glia_mcs], legc[g,nsc_glia_mcs], col = mcmd$color[nsc_glia_mcs], cex = CEX,
            pch = 16, xlab = 'Metacell mean day', ylab =paste0(g, ' expression'))
dev.off()



## Fig 5E -- Tnc
plot_insulation_and_shaman_series(2, hvm_ls[[2]], file_path = fig_5e_path)


## Fig5F
p_tnc <- plotRegion(interv = interv_tnc, genome = genome, hic_data = hic_data_tnc, 
            plotHicTracks, plotGeneTrack, trackdb, scoreTracksToExtract, point_cex = 0.75)
pdf(fig_5f_path, h = 600/71, w = 510/71)
print(p_tnc)
dev.off()


## Fig 5G
pdf(fig_5g_path, h = 350/71, w = 550/71)
par(cex.lab = 2, cex.axis = 1.5, mar = MAR)
g <- 'Tnc'
plot(mcmd$mean_day[nsc_glia_mcs], legc[g,nsc_glia_mcs], col = mcmd$color[nsc_glia_mcs], cex = CEX,
            pch = 16, xlab = 'Metacell mean day', ylab = paste0(g, ' expression'))
dev.off()


## Fig 5H

pdf(fig_5h_path, h = 600/71, w = 600/71)

par(mar = c(5,5,3,1), cex.lab = 2, cex.axis = 1.5, cex.main = 2)
LWD <- 3
plot(head(qs, -1), tbl1_rm, type = 'l', col = 'purple', 
     ylim = c(0,2e-2), lwd = 3, 
     ylab = 'Probability of neighbors per region type', 
     xlab = 'Distance [bp]') 
lines(head(qs, -1), tbl2_rm, col = 'orange', type= 'l', lwd = 3)
legend('bottomright', legend = c('Deinsulating', 'Insulating'), col = c('orange', 'purple'), 
               lty = c(1,1), lwd = rep(LWD,2), cex = 1.5)

dev.off()


# Fig 5I

pdf(fig_5i_path, h = 450/71, w = 950/71)

mari <- c(5,5,1,1)
par(mfcol = c(1,3), mar = mari, cex.lab = 2, cex.axis = 2)
Y_DFL <- 7
DELTA_Y <- 2
ncl <- setNames(unique(pval_sig_df$cre_type), names(score_lst))
gtf <- setNames(c("NSC TSS", 'Astro TSS', 'IPC TSS'), c('nsc_gene', 'astro_gene', 'ipc_gene'))

sss <- lapply(seq_along(score_lst), function(i) {
    y <- score_lst[[i]]
    nm <- names(score_lst)[[i]]
    lapply(c('nsc_gene', 'astro_gene', 'ipc_gene'), function(gt) {

        f <- y$gene_type == gt
        pltmt <- y[f,grep('score', colnames(y))] - rowMeans(as.matrix(y[f,grep('score', colnames(y))]))
        pv_df_h <- pval_sig_df[pval_sig_df$cre_type == ncl[names(score_lst)[[i]]] & pval_sig_df$tss_type == gt,]
        if (gt == 'nsc_gene')  {
            mari[[2]] <- 5
        } else {
            mari[[2]] <- 2
        }
        
        par(mar = mari)
        boxplot(pltmt, ylim = c(-10, 23), 
                xaxt = 'n',
                ylab = '',
                xlab = '', main = '')
        axis(1, at = 1:5, labels = paste0('E', 13:17))
        if (i == 1) {
          title(ylab = 'delta SHAMAN D from mean', line = 3)
        }
        
          title(xlab = 'Time point', line = 3)
        lines(c(-5,10), rep(0,2), col = 'red')
        text(4, 22, gtf[[gt]], cex = 2)
        if (nrow(pv_df_h) >= 1) {
            npi <- get_num_asterisks(as.numeric(pv_df_h[,'pval']))
            
            for (i in 1:nrow(pv_df_h)) {
                xh <- as.numeric(pv_df_h[i,c('x1','x2')])

                lines(xh, rep(Y_DFL+DELTA_Y*i - 1,2), lwd = 3.5, col = 'darkgray')

                plot_asterisks(x = mean(xh), y = Y_DFL +DELTA_Y*i,  spacing_factor = 0.075, npj = npi[[i]])
                
            }
        }
    })
})

dev.off()