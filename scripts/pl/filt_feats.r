library(metacell)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'pl'

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))



mcell_add_gene_stat(mat_id = nm, gstat_id = nm, force=T)

mcell_gset_filter_varmean(gset_id=nm, gstat_id=nm, T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id=nm, gstat_id=nm, T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id=nm, gset_id=nm)

mcell_gset_split_by_dsmat(gset_id=nm, mat_id=nm, K = 96, force = T)

mcell_plot_gset_cor_mats(gset_id=nm, scmat_id=nm)

ifn1_genes = c('Isg15', 'Wars', 'Ifit1')
cell_cyc = c('Mki67', 'Pcna', 'Hist1h', 'Smc4', 'Mcm3', 'Top2a')
stress = c('Fos', 'Hsp90ab1', 'Hspa1a', 'Hif1a')
misc = c('Xist', 'Tsix')
star_genes = c(ifn1_genes, cell_cyc, stress, misc)

rb = grep('Rps|Rpl', rownames(mat@mat), v=T)
star_rb = c(star_genes, rb)

feats = scdb_gset(nm)
mat = scdb_mat(nm)


feat_cors = tgs_cor(t(as.matrix(mat@mat[names(feats@gene_set),])), spearman = T)
star_gene_cor_genes = sapply(star_rb, function(g) {
    if (g %in% rownames(feat_cors)) {
        tbl = sort(table(feats@gene_set[head(names(feat_cors[g,order(feat_cors[g,], decreasing = T)]), 100)]), decreasing = T)
                                     }
    
})

num_clust_hits = sapply(star_gene_cor_genes[!sapply(star_gene_cor_genes,is.null)], 
                function(x) names(x)[x>= 10]) %>% unlist %>% table

cust_exclude = c()

clusts_to_remove = as.numeric(names(num_clust_hits)[num_clust_hits >= 2 & !(names(num_clust_hits) %in% cust_exclude)])
clusts_to_remove = c(clusts_to_remove, feats@gene_set[names(feats@gene_set) == 'Xist'])

feats_removed = sort(names(feats@gene_set)[feats@gene_set %in% clusts_to_remove])
print('Genes removed from feature genes:')
feats_removed

mcell_gset_remove_clusts(gset_id = nm, filt_clusts = clusts_to_remove, 
                         new_id = paste0(nm, '_lateral'), reverse=T)

lat = scdb_gset(paste0(nm, '_lateral'))
lat@gene_set = c(lat@gene_set, setNames(rep(NA, length(star_rb)), star_rb))
scdb_add_gset(paste0(nm, '_lateral'), lat)

mcell_gset_remove_clusts(gset_id = nm, filt_clusts = clusts_to_remove, 
                         new_id = paste0(nm, '_f'), reverse=F)