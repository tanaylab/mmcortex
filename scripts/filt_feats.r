library(metacell)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
nm = 'pl'

set.seed(1337)


db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))



mcell_add_gene_stat(mat_id = nm, gstat_id = nm, force=T)

mcell_gset_filter_varmean(gset_id=glue("{nm}_feats"), gstat_id=nm, T_vm=0.1, force_new=T)
mcell_gset_filter_cov(gset_id=glue("{nm}_feats"), gstat_id=nm, T_tot=100, T_top3=4)
mcell_plot_gstats(gstat_id=nm, gset_id=glue("{nm}_feats"))

mcell_gset_split_by_dsmat(gset_id=glue("{nm}_feats"), mat_id=nm, K = 96, force = T)

mcell_plot_gset_cor_mats(gset_id=glue("{nm}_feats"), scmat_id=nm)

ifn1_genes = c('Isg15', 'Wars', 'Ifit1')
cell_cyc = c('Mki67', 'Pcna', 'Hist1h', 'Smc4', 'Mcm3', 'Top2a')
stress = c('Fos', 'Hsp90ab1', 'Hspa1a', 'Hif1a')
misc = c('Xist', 'Tsix', 'Rps', 'Rpl')
star_genes = c(ifn1_genes, cell_cyc, stress, misc)

feats = scdb_gset(glue("{nm}_feats"))

genes_in = map(star_genes, function(x) grep(x, names(feats@gene_set), ignore.case = F, v = T)) %>% unlist
uu = unique(feats@gene_set[names(feats@gene_set) %in% genes_in])

star_in = lapply(star_genes, function(x) grep(x, names(feats@gene_set), v=T)) %>% unlist
print(paste0('star_in = ', star_in))
clusts_to_remove = feats@gene_set[star_in] %>% unique
clusts_to_remove = c(70, 83, clusts_to_remove)   ## clusts added from CC analysis after constructing metacells
print(paste0('clusts_to_remove = ', as.character(clusts_to_remove)))
gset_nm = glue('{nm}_feats')
print(head(table(feats@gene_set)))
mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = clusts_to_remove, 
                         new_id = glue('{nm}_lateral'), reverse=T)
mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = clusts_to_remove, 
                         new_id = glue('{gset_nm}_f'), reverse=F)