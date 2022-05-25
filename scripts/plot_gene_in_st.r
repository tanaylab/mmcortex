library(metacell)
library(dplyr)
set.seed(1337)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
db_path = file.path(wd, 'scdb')
scdb_init(db_path, force_reinit = T)
scfigs_init("figs/")



mc_subset_plot_gene = function(mc, mc2d, mc_md, type, gene) {
    # st_fold = file.path('figs', 'st_marks', type)
    # if (!dir.exists(st_fold)) {dir.create(st_fold)}
    # png(paste0(file.path(st_fold, gene), '.png'))
    mcs = select(filter(mc_md, Bonev_annotation == type), 'mc') %>% unlist
    mc_x = mc2d@mc_x[names(mc2d@mc_x) %in% mcs]
    mc_y = mc2d@mc_y[names(mc2d@mc_y) %in% mcs]
    sc_x = mc2d@sc_x[names(mc2d@sc_x) %in% names(mc@mc[mc@mc %in% mcs])]
    sc_y = mc2d@sc_y[names(mc2d@sc_y) %in% names(mc@mc[mc@mc %in% mcs])]
    plot(sc_x, sc_y, pch = 16, col = 'gray', main = gene)
    col_map = colorRampPalette(colors = c('white', 'red'))(100)
    vals = mc@mc_fp[gene,mcs]
    vals_t = 100*(vals - min(vals))/(max(vals) - min(vals))
    points(mc_x, mc_y, pch = 16, cex = 3, col = col_map[vals_t])
    text(mc_x, mc_y, mcs, cex = 0.7, col = 'black')
    # dev.off()
}

mc = scdb_mc('all_rec_bon_1')
mc2d = scdb_mc2d('all_umap')
mc_md = vroom::vroom('./BonevCollab//mc_metadata_new.tsv')
marks = scdb_gset('all_markers_f')
feats = scdb_gset('all_feats_f')
type = 'IPC'
gene = 'Hes5'

mc_subset_plot_gene(mc, mc2d, mc_md, type, gene)


# mcfp = mc@mc_fp[names(marks@gene_set), select(filter(mc_md, Bonev_annotation == type), 'mc') %>% unlist]

# diff = apply(mcfp,1, max) - apply(mcfp,1, min)
# gvar = apply(mcfp,1, var)
# mcfp = mcfp[order(gvar, decreasing = T),]
# lapply(rownames(head(mcfp, 100)), function(x) mc_subset_plot_gene(mc, mc2d, mc_md, type, x))