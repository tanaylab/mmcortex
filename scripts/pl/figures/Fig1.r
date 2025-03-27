library(metacell)
# devtools::load_all('~/src/metacell.flow')
library(metacell.flow)
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

mc = scdb_mc(nm)
mat = scdb_mat(nm)
mcf = scdb_mctnetflow(nm)
mc2d <- scdb_mc2d(mc2d_id)
load(file.path(wd, 'output/metacell_model/fig_1_data.rda'))

source('./scripts/util.r')

mcmd = vroom::vroom('./output/metacell_model/mcmd_pl_cort.tsv')
col_key <- tibble::deframe(unique(mcmd[,c('cell_type', 'color')]))
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC','iCPN/CfuPN', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s)[order(mcmd$mean_day[which(mcmd$cell_type == s)])], 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )))

cust_st_ord2 = c('OPCs','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
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


clrs <- c('red3', 'orange3', 'yellow3', 'green4', 'blue3', 'purple2')
clrmp <- colorRampPalette(clrs)(1000)

clrmp_abs <- colorRampPalette(c('white', 'orange', 'red', 'purple', 'black'))(1000)
brks_abs <- seq(-16.6,-10, l=1000)

clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)
brks_rel <- seq(-3,3, l=1000)

ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])),
                 'mean_day' = setNames(colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))(100),
                                      seq(13,18,l=100)))

legc = log2(1e-05 + mc@e_gc)
## End

device <- 'pdf'
fig_1a_path <- glue::glue('./output/paper_figs/Fig1/Fig1A.{device}')
fig_1a_legend_path <- glue::glue('./output/paper_figs/Fig1/Fig1A_legend.{device}')
fig_1b_path <- glue::glue('./output/paper_figs/Fig1/Fig1B.{device}')
fig_1c_path <- glue::glue('./output/paper_figs/Fig1/Fig1C.{device}')
fig_1c_legend_path <- glue::glue('./output/paper_figs/Fig1/Fig1C_legend.{device}')
fig_1d_path <- glue::glue('./output/paper_figs/Fig1/Fig1D.{device}')
fig_1d_color_bar_path <- glue::glue('./output/paper_figs/Fig1/Fig1D_color_bar.{device}')
fig_1e_path <- glue::glue('./output/paper_figs/Fig1/Fig1E.{device}')
fig_1f_path <- glue::glue('./output/paper_figs/Fig1/Fig1F.{device}')


## Fig 1A

#mc2d
mc2d_id <- 'pl_cort_not_cor_cc'

min_edge_l <-0
edge_w <- .31
short_edge_w <-0

mcp_2d_cex <- tgconfig::get_param(param = "mcell_mc2d_cex", package = 'metacell')
sc_cex <- tgconfig::get_param(param = "sc_cex", package = 'metacell')
pdf(fig_1a_path, h = 1500/71, w = 1500/71)
plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=mcmd$color[mc@mc[names(mc2d@sc_x)]], 
        cex=sc_cex, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
fr <- mc2d@graph$mc1
to <-mc2d@graph$mc2
dx <-mc2d@mc_x[fr]-mc2d@mc_x[to]
dy <-mc2d@mc_y[fr]-mc2d@mc_y[to]
f <-sqrt(dx*dx+dy*dy) > 0
segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
                        lwd=ifelse(f, edge_w, short_edge_w))
points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=mcmd$color)
text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)
dev.off()


#legend
df = data.frame(color_key[order(match(color_key$cell_type, cust_st_ord), decreasing = F),])
l = nrow(df)
scale_y = 1
pdf(fig_1a_legend_path, width = 600/71, height = 800/71)
par(mar = c(4,1,4,0), bty = 'n')
plot(rep(0.9,l), scale_y*seq(l,1,-1), pch = 16, cex = 5, col = df$color, ylim = c(0,scale_y*l+1),
     xlim = c(-1,60),xlab = '', ylab = '',xaxt = 'n',yaxt = 'n')
text(rep(6,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 2.81, df$cell_type)
dev.off()


## Fig 1B
flow_thresh = min(mcf@edge_flows[mcf@edge_flows > 0])
tgconfig::set_param(param = 'mc_plot_device', value = 'pdf', package = 'metacell')
mctnetwork_plot_net_YSh(nm, nm, plot_pdf = TRUE, h = 32, w = 24, fn = fig_1b_path, 
        mc_ord = cust_mc_ord_st, flow_thresh = flow_thresh)

## Fig 1C


md_clvls <- clrmp[1+round(999*(mcmd$mean_day-13)/(18-13))]
days <- paste0('E', 13:18)
pdf(fig_1c_path, h = 1500/71, w = 1500/71)
par(cex.main = 6, cex.lab = 2, mar = c(5,5,4,3))
plot(mc2d@sc_x, mc2d@sc_y, pch = 16, col = md_clvls[mc@mc[names(mc2d@sc_x)]], 
        main = '', xaxt = 'n', yaxt  = 'n', bty = 'n', xlab = '', ylab = '')
legend('topright', legend = days, col = clrs, pch = rep(16, length(clrs)), cex = 5, bty = 'n')
dev.off()



## Fig 1D

pdf(fig_1d_path, w=8, h=8)
shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)

plot(mc2d@sc_x, mc2d@sc_y, pch=16, cex=0.8, col=ifelse(f[names(mc2d@sc_x)], "lightgray","black"),
        bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=1.8, bg=shades[101 - mc_cc$cc_score])
dev.off()

## Fig 1D legend
min_val = min(mc_cc$cc_score)
max_val = max(mc_cc$cc_score)
plot_color_bar(seq(min_val-1, max_val,l=101), shades, 
                        height = 400/71, width = 400/71, device = pdf, 
                        fig_fn = fig_1d_color_bar_path)


## Fig 1E
pltmt = mc@mc_fp[goi,cust_mc_ord_st]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')

top_ha <- anno_simple(x = mcmd$cell_type[cust_mc_ord_st], col = ann_colors[['cell_type']])
md_col_fun <- circlize::colorRamp2(breaks = 13:18, 
                                    colors = rev(c('purple', 'blue', 'green', 'yellow', 'orange', 'red')))
bottom_ha = anno_lines(x = col_annot$mean_day[cust_mc_ord_st], 
                       axis_param = list(gp = gpar(fontsize = 10)),
                                                    height =unit(1, 'cm'))

ch <- Heatmap(pltmt, name = "marker\nheatmap",
              top_annotation = HeatmapAnnotation(`cell type` = top_ha), 
              bottom_annotation = HeatmapAnnotation(`mean day` = bottom_ha, show_legend = F), 
              col = circlize::colorRamp2(breaks = c(minv, 0, maxv/2, maxv), colors = pltt), show_column_names = F, 
              row_names_gp = gpar(fontsize = 14),
              column_split = factor(names(cust_mc_ord_st), levels = cust_st_ord), column_gap = unit(.02, 'mm'),
                              column_title_gp = gpar(fontsize = 0),
                        cluster_columns = F, cluster_rows = F
                        )

pdf(fig_1e_path, w = 1000/71, h = 600/71)
draw(ch)
dev.off()

## Fig 1F

nsc_flow_out <- do.call('rbind', lapply(mcf@mc_forward, function(x) rowSums(tgstat::tgs_matrix_tapply(x[which(mcmd$cell_type == 'NSC'),], mcmd$cell_type, sum))/sum(colSums(x))))
nsc_flow_out_n <- nsc_flow_out/rowSums(nsc_flow_out)
aa <- paste0('E', 13:17)
bb <- paste0('E', 14:18)
xlabs <- apply(cbind(aa, rep('->', length(aa)), bb), 1, paste, collapse = ' ')

pdf(fig_1f_path, h = 550/71, w = 550/71)
par(las = 2, mar = c(12,8,1,1), cex.main = 2, cex.axis = 2, cex.lab = 1.5)
ct <- 'NSC'
plot(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], 
     type = 'l', ylim = c(0,max(nsc_flow_out_n)), lwd = 3, 
      ylab = '', xaxt = 'n', xlab = '')
axis(1, at = 1:5, labels = xlabs)
title(ylab = 'Relative flow', line = 4, cex.lab = 2)
dummy <- sapply(colnames(nsc_flow_out_n), function(ct) {
        lines(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], lwd = 3)
        points(1:5, nsc_flow_out_n[,ct], col = color_key$color[color_key$cell_type == ct], pch = 16, cex = 2)
})
dev.off()