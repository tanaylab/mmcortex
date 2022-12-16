library(metacell)
scdb_init('scdb')
source('./scripts/paper_scripts/util.r')

library(Matrix)
nm = 'pl_cort'
mat_id = nm
mc_id = nm

mc <- scdb_mc(mc_id)
mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
color_key = unique(mcmd[,c('cell_type', 'color')])

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s)))
                                                                )
                                )
                        )                         

goi = c('Pou3f2', 'Cux1', 'Cux2', 'Neurod1', 'Neurog2', 'Id4',
         'Eomes', 'Hes1', 'Apoe', 'Sox5', 'Tbr1', 'Foxp1', 'Nfia', 'Islr2', 
         'Zbtb20', 'Bcl11b', 'Fezf2', 'Satb2', 'Mef2c', 'Nhlh1', 'Tle4',
        'Rnd2',  'Runx1t1', 'Mapt', 'Mki67', 'Pcna',
        'Fabp7', 'Olig1', 'Ldb2', 'Gadd45g', 'Syt4')
marks_filt = goi

col_annot = mcmd[,c('metacell', 'cell_type')]
col_annot = tibble::column_to_rownames(col_annot, 'metacell')
ann_colors = list('cell_type' = tibble::deframe(unique(mcmd[,c('cell_type', 'color')])))

pdf('./paper_figs/fig1/fig1.pdf', width = 8.3, height = 8.3)

mc2d_bord_r <- 0.6
mc2d_bord_l <- 0.02
mc2d_bord_b <- 0.35
mc2d_bord_t <- 0.98
legend_bord_l <- mc2d_bord_r*0.8
legend_bord_b <- 0.7


par(fig = c(mc2d_bord_l,mc2d_bord_r, mc2d_bord_b,mc2d_bord_t), mar = rep(0.2,4))
bb <- my_mcell_mc2d_plot(nm,
                            show_mcid = F,
                            edge_w = 0.05,
                            short_edge_w = 0,
                            min_edge_l = 5,
                            sc_cex = 0.2
)
par(fig = c(legend_bord_l,mc2d_bord_r, legend_bord_b,mc2d_bord_t), new = T)
df = data.frame(color_key[order(match(color_key$cell_type, cust_st_ord)),])
l = nrow(df)
scale_y = 2
plot(rep(0.93,l), scale_y*seq(l,1,-1), pch = 16, cex = 1.5, col = df$color, xlim = c(0.92, 1), ylim = c(0.5,scale_y*l+1),
    xlab = '', 
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(rep(0.94,l), scale_y*seq(l,1,-1), adj = c(0, 0.5), cex = 0.75, df$cell_type)

tag = nm
m_0 = 0.0025
s_0 = 0.001
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
mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

colnames(mc_cc) = c('mc', 'cc_score')

mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))
mc2d = scdb_mc2d('pl_cort')

shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)

cc_l <- mc2d_bord_r + 0.02
cc_r <- 0.98
cc_t <- mc2d_bord_t
cc_b <- 0.7
par(fig = c(cc_l, cc_r, cc_b, cc_t), mar = rep(0.2,4), new = T)
plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.1, col=ifelse(f[names(mc2d@sc_x)], "lightgray","black"), 
            main = 'Cell-cycle score',xlab = '', 
            ylab = '',
            xaxt = 'n',
            yaxt = 'n'
            )
points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=0.75, bg=shades[101 - mc_cc$cc_score])

cc_l_l <- 0.9*cc_r
cc_l_r <- cc_r
cc_l_t <- cc_t
cc_l_b <- 0.9*cc_t
par(fig = c(cc_l_l, cc_l_r, cc_l_b, cc_t), mar = rep(0.2,4), new = T)
plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6)) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
    rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
    rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
    text(19, show_vals_ind,cex = 0.2, labels=round(vals[show_vals_ind], 3), pos=4)
}
shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
min_val = min(mc_cc$cc_score)
max_val = max(mc_cc$cc_score)
plot_color_bar(seq(min_val-1, max_val,l=101), shades)

gg_l <- cc_l
gg_r <- cc_r
gg_t <- cc_b - 0.02
gg_b <- mc2d_bord_b
par(fig = c(gg_l, gg_r, gg_b, gg_t), mgp = c(0.75,0.15,0), cex.axis = 0.5, cex.lab = 1, mar = c(2,2,0.02,1), new = T)
mcell_mc_plot_gg('pl_cort', 'Bcl11b', 'Satb2', text_cex = 0, cex = 1, use_egc = T, main = 'log2 mRNA fraction per metacell')

pltmt = mc@mc_fp[goi,cust_mc_ord_st]
pltmt = pltmt[order(apply(pltmt, 1, which.max)),]
pltmt = log2(pltmt)
minv = min(pltmt)
maxv = max(pltmt)
l = 100
pltt = c('blue4', 'white', 'red4', 'black')
brks = c(seq(minv,0,l=round(l/(length(pltt)-1))),
         seq(0+1/l,maxv/2,l=round(l/(length(pltt)-1))+1),
        seq(maxv/2+1/l,maxv,l=round(l/(length(pltt)-1))+1))
pp = pheatmap::pheatmap(pltmt, color = colorRampPalette(pltt)(l), breaks = brks, fontsize = 14, 
                        annotation_legend = F,
                        cluster_cols = F, cluster_rows = F, 
                        annotation_col = col_annot, 
                        annotation_colors = ann_colors, 
                        show_colnames = F, fontsize_row = 14, silent = T)

hm_l <- mc2d_bord_l
hm_r <- cc_r
hm_t <- mc2d_bord_b - 0.02
hm_b <- 0.02
par(fig = c(hm_l, hm_r, hm_b, hm_t), new = T)
grid::grid.draw(pp$gtable)

dev.off()