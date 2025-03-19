library(beeswarm)

# wd <- '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/'
wd <- '.'
setwd(wd)
library(misha)
gsetroot('/home/aviezerl/mm10')
source(file.path(wd, 'scripts/util.r'))

load(file = './output/sequence_modeling/fig_6_data.rda')

dir.create('./output/paper_figs/Fig6/')
device <- 'pdf'
fig_6a_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6A.{device}'))
fig_6b_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6B.{device}'))
shap_color_bar_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/shap_color_bar.{device}'))
fig_6c_path <- file.path(wd, glue::glue('./output/paper_figs/Fig6/Fig6C.{device}'))

## Fig 6A

pdf(fig_6a_path, h = 500/71, w = 600/71)
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5, mar = c(5,5,2,1))
plot(delta_ipc_nsc[names(pred_meth_all)], pred_meth_all, pch = 16, cex = 0.25,
             ylab = 'Predicted delta ATAC', xlab = 'Observed delta ATAC', 
                main = '')
abline(0,1,col = 'red', lwd = 2)
text(-1.2, 1.2, labels = paste0('R^2 = ', round(mean(ev_all$Rsquare), 2)), cex = 2)
dev.off()


## Fig 6B

pdf(fig_6b_path, h = 1050/71, w = 1500/71)
par(las = 2, mar = c(14,12,2,1), cex.axis = 3, cex.lab = 6)
beeswarm(value ~ lvl, 
         method = 'compactswarm',
         spacing = 0.2,
         corral = 'wrap',
         data = samp_df,
         pwcol = samp_df$color, 
         ylim = c(-0.9,0.55), ylab = '', xlab = '')
title( ylab = 'SHAP value', line = 7)
grid(lwd = 5)
v <- norm_integral_over_feature_contrib
text(match(names(v), levels(samp_df$lvl)), -0.45, -0.25, col = 'black', labels = round(v,2), cex = 3, srt = 90)

dev.off()

## SHAP value legend
clrmp_rel <- colorRampPalette(c('blue3', 'white','red3'))(1000)

plot_color_bar(vals = seq(0,1,l=length(clrmp_rel)), 
            cols = clrmp_rel, height = 500/71, width = 250/71, device = pdf,
            show_vals_ind = seq(1,length(clrmp_rel),l=11), fig_fn = shap_color_bar_path)


## Fig 6C

pdf(fig_6c_path, w=900/71,h=900/71)

par(las = 2, mar = c(10,8,2,1), mfrow = c(3,1), cex.lab = 2, cex.axis = 1.5)
plot(-1,-1, xlim=c(1.1,6.8), ylim=c(-16.6,-15.0), ylab="", xlab='', xaxt = 'n', yaxt = 'n')

## NSC/IPC ATAC vs methylation and E-box energy
for(m in seq(-0.05,0.55,0.1)) {
        f = !is.na(cres$methylation) & cres$methylation > m & cres$methylation < m+0.1
        boxplot(split(cres$NSC[f], e_box_rng[f]), add=T, boxwex = 0.05, 
                                        at=m*1.4+(1:5)*1.2, xaxt='n', col=col_key[['NSC']], cex=0.1)
        boxplot(split(cres$IPC[f], e_box_rng[f]), add=T, boxwex = 0.05, yaxt = 'n',
                                        at=m*1.4+0.07+(1:5)*1.2, xaxt='n', col=col_key[['IPC']], cex=0.1 )
        axis(1, at = m*1.4+0.035+(1:5)*1.2, labels = rep(glue::glue("({m},{m+0.1})"), 5))
}
text(1:5*1.2+0.4, rep(-14.85,5), labels = glue::glue('E-box energy q=[{seq(0,0.8,0.2)},{seq(0.2,1,0.2)}]'), xpd = T, cex=1.5)
title(xlab = 'methylation', line = 8)
title(ylab = 'Normalized accessibility', line = 5)
legend(y = -15.15, x = 1.8, legend = c('NSC', 'IPC'), fill = col_key[c('NSC', 'IPC')], cex = 1, xpd = T, bg = 'white')
abline(h=-16)
abline(h=-16.25)
# dev.off()

par(mar = c(6,8,2,1))
## NSC/IPC ATAC vs # of proximal CREs with high E-box/T-box energy and E-box energy
plot(-1,-1, xlim=c(1,5.8), ylim=c(-16.6,-15.0), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
for(prox in 0:7) {
        f = (reg_both_a7) == prox
        boxplot(split(cres$NSC[f], e_box_rng[f]), add=T, boxwex = 0.04, 
                                                at=prox*0.12+(1:5), xaxt='n', col=col_key[['NSC']], cex=0.1)
        boxplot(split(cres$IPC[f], e_box_rng[f]), add=T, boxwex = 0.04,
                                                at=prox*0.12+0.05+(1:5), xaxt='n', col=col_key[['IPC']], cex=0.1 )
        if (prox %in% c(0,7)) {
            axis(1, at = prox*0.12+(1:5), labels =  rep(prox, 5))
        }
}
arrows(x0 = (1:5)+0.1, y0 = -16.75, x1 = (1:5) + 0.75, y1 = -16.75, xpd = T, lwd = 1.5, angle = 20)
text((1:5)+ 0.4, rep(-14.85,5), labels = glue::glue('E-box energy q=[{seq(0,0.8,0.2)},{seq(0.2,1,0.2)}]'), cex=1.5, xpd = T)
text((1:5) + 0.375, rep(-16.83,5), labels = '# of neighbors', cex=1.5, xpd = T)
title(ylab = 'Normalized accessibility', line =5)
title(xlab = '# of high-affinity E-box or T-box proximal CREs', line = 4)
abline(h=-16)
abline(h=-16.25)
# dev.off()


## NSC/IPC ATAC vs # of proximal activating elements and E-box energy
par(mar = c(9,8,2,1))
plot(-1,-1, xlim=c(1,5.8), ylim=c(-16.6,-15.0), xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
for(prox in 0:6) {
        f = (reg_d_a6) == prox
        boxplot(split(cres$NSC[f], e_box_rng[f]), add=T, boxwex = 0.04, yaxt = 's',
                                                at=prox*0.12+1:5, xaxt='n', 
                                                col=col_key[['NSC']], cex=0.1)
        boxplot(split(cres$IPC[f], e_box_rng[f]), add=T, boxwex = 0.04, yaxt = 'n',
                                                at=prox*0.12+0.05+1:5, xaxt='n', 
                                                col=col_key[['IPC']], cex=0.1 )
        if (prox %in% c(0,6)) {
            axis(1, at = prox*0.12+1:5, labels = rep(prox, 5))
        }
}
arrows(x0 = (1:5) + 0.1, y0 = rep(-16.75, 5), x1 = (1:5) + 0.675, y1 = rep(-16.75, 5), xpd = T, lwd = 1.5, angle = 20)
arrows(x0 = 1.1, y0 = -17.2, x1 = 5.9, y1 = -17.2, xpd = T, lwd = 2, angle = 20)
text(3.5, -17.35, labels = 'E-box energy range', cex=2, xpd = T)
text((1:5) + 0.35, rep(-16.85,5), labels = '# of neighbors', cex=1.5, xpd = T)
text((1:5)+ 0.4, rep(-14.85,5), labels = glue::glue('E-box energy q=[{seq(0,0.8,0.2)},{seq(0.2,1,0.2)}]'), xpd = T, cex=1.5)
title(ylab = 'Normalized accessibility', line =5)
title(xlab = '# of activating proximal CREs', line = 4)
abline(h=-16)
abline(h=-16.25)

dev.off()