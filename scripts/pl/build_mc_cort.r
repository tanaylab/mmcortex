library(metacell)


wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)
db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))
nm = 'pl'
nm_new = paste0(nm, '_cort')
set.seed(1337)

mat_id <- nm_new
feats_gset_id <- paste0(nm, '_filt_lat')
marks_gset_id <- paste0(nm_new, '_marks')
coc_id <- nm_new
mc_id <- nm_new
cgraph_id <- nm_new

# mcell_add_cgraph_from_mat_bknn(mat_id=nm_new,
# mcell_add_cgraph_from_mat_bknn(mat_id=mat_id,
#                 gset_id = feats_gset_id,
#                 graph_id=nm_new,
#                 K=50,
#                 dsamp=T)

# mcell_coclust_from_graph_resamp(
#                 coc_id=nm_new,
#                 graph_id=nm_new,
#                 min_mc_size=10,
#                 p_resamp=0.75, n_resamp=500)

# mcell_mc_from_coclust_balanced(mc_id=nm_new, coc_id=nm_new, mat_id=nm_new, K=20, min_mc_size=11)
mcell_mc_from_coclust_balanced(mc_id=mc_id, coc_id=coc_id, mat_id=mat_id, K=20, min_mc_size=10)

mcell_gset_from_mc_markers(gset_id = marks_gset_id, mc_id = mc_id)

mcell_gset_from_mc_markers(gset_id = paste0(marks_gset_id, '_f'), mc_id = mc_id, 
        blacklist_gset_id = paste0(nm, '_lateral'))

mc = scdb_mc(mc_id)

feats = scdb_gset(feats_gset_id)

feats_f = feats
feats_f@gene_set = feats_f@gene_set[names(feats_f@gene_set) %in% rownames(mc@e_gc)]

scdb_add_gset(nm_new, feats_f)

### Send MC to object to Dror Bar for denoising ###
### Substitute denoised UMI matrix (divided by MC sizes) instead of mc@e_gc
### Calculate @mc_fp by hand

my_mc_compute_fp = function(mc_mat, norm_by_mc_size=T, min_total_umi=10) {
	f_g_cov = rowSums(mc_mat) > min_total_umi
	clust_geomean = mc_mat[f_g_cov,]
	rownames(clust_geomean) = rownames(mc_mat)[f_g_cov]
	if (norm_by_mc_size) {
		mc_meansize = colSums(mc_mat)
		ideal_cell_size = pmin(1000, median(mc_meansize))
		g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
	} else {
		g_fp = clust_geomean
	}
	#normalize each gene
	fp_reg = 0.1
	#0.1 is defined here because 0.1*mean_num_of_cells_in_cluster
	#is epxected to be 3-7, which means that we regulairze
	#umicount in the cluster by 3-7.
	g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)
	return(g_fp_n)
}

my_mc_compute_e_gc= function(mc_mat, norm_by_mc_meansize=T) {
	f_g_cov = rowSums(mc_mat) > 10
	e_gc <- mc_mat[f_g_cov,]
	rownames(e_gc) = rownames(mc_mat)[f_g_cov]
	if (norm_by_mc_meansize) {
		mc_meansize = colSums(e_gc)
		e_gc = t(t(e_gc)/as.vector(mc_meansize))
	}
	return(e_gc)
}

dnmc <- anndata::read_h5ad(file.path(wd, 'output/metacell_model', 'denoise_mc.h5ad'))
dnmc_egc <- my_mc_compute_e_gc(t(dnmc$X))
colnames(dnmc_egc) <- 1:ncol(dnmc_egc)
dnmc_fp <- my_mc_compute_fp(t(dnmc$X))
colnames(dnmc_fp) <- 1:ncol(dnmc_fp)

### Problematic MCs were found based on co-expression of callosal and
### corticofugal markers (e.g. Satb2 vs. Bcl11b) while having high neuron maturity markers
### see "mg_bon" dataframe for markers

# problematic_mcs <- sort(unique(c(221, 222, 217, 492, 216, 212, 488, 489, 218, 234, 491, 481, 485, 220, 485, 213, 214, 215,205, 486)))

# mc@e_gc <- dnmc_egc[, -problematic_mcs]
# colnames(mc@e_gc) <- 1:ncol(mc@e_gc)
# mc@mc_fp <- dnmc_fp[, -problematic_mcs]
# colnames(mc@mc_fp) <- 1:ncol(mc@mc_fp)

# print('Dimensions before removing outliers')
# print(dim(mc@e_gc))
# mc <- mc_set_outlier_mc(mc, problematic_mcs)
# print('Dimensions after removing outliers')
# print(dim(mc@e_gc))
scdb_add_mc(mc_id, mc)
scdb_init(db_path, force_reinit = T)