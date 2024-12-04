library(metacell)

wd = '/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex'
set.seed(1337)

db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
scfigs_init(file.path(wd, "figs"))

nm = 'pl'

mat = scdb_mat(nm)

mc = scdb_mc(nm)

# mg_bon = read.delim('./BonevCollab//marker_genes.tsv', sep='\t') %>% 
#             apply(1, function(x) stringr::str_split(x, pattern = ' ')) %>%
#             purrr::map(1) %>% purrr::map(function(x) c(x[[1]], stringr::str_split(x[[length(x)]], ',')))

# mg_bon = unique(data.frame(cbind(purrr::map(mg_bon, 1), purrr::map(mg_bon, 2))))
# colnames(mg_bon) = c('st', 'marks')
# mg_bon


# feats = scdb_gset(paste0(nm, '_f'))
# feats = names(feats@gene_set)

non_cort_genes = c('Reln', 'Lhx5','Gad2', 'Sst' , 'Lhx6', 'Nrxn3','Gsx2','Dlx2','Aif1', 'C1qb', 'Hexb', 'Igfbp7')
non_cort_mcs = sort(unique(unlist(sapply(non_cort_genes, function(g) which(mc@mc_fp[g,] >= 1.5)))))
non_cort_mcs
non_cort_cells = names(mc@mc)[mc@mc %in% non_cort_mcs]

cells_ignore_new = union(non_cort_cells, mat@ignore_cells)

nm_new = paste0(nm, '_cort')

mcell_mat_ignore_cells(new_mat_id = nm_new, mat_id = nm, ig_cells = cells_ignore_new, reverse = F)




# marks_all = sort(unique(unlist(mg_bon$marks)))
# marks_all = marks_all[!(marks_all %in% c('Apoe','Sox2', 'Vim'))]

# legcz_marks = mc@mc_fp[marks_all,]
# mark_scores = sapply(mg_bon$marks, function(x) {
#     gns = as.character(stringr::str_split(x, ', '));
#     gns = gns[!(gns %in% c('Apoe','Sox2', 'Vim'))]
#     print(gns)
#     return(apply(legcz_marks[gns,], 2, mean))
# })

# colnames(mark_scores) = mg_bon$st

# top_3_st = apply(mark_scores, 1, function(x) {ord = order(x, decreasing = T);
#                                   return(head(colnames(mark_scores)[ord], 3))})

# npc_mc = which(top_3_st[1,] %in% c('NPC_GE','NPC_CGE'))
# npc_ge_marks_mat = t(mc@mc_fp[c('Pou3f2', 'Rnd2', 'Sstr2', 'Nrp1', 'Dlx1','Dlx2', 'Gsx2', 'Lhx9'),npc_mc])

# mcs_not_npc = as.numeric(head(rownames(slanter::slanted_reorder(npc_ge_marks_mat)), 7))
# # mcs_not_npc

# non_cort_st = c('IN_MGE', 'PC', 'CR', 'MSN', 'NPC_GE', 'NPC_CGE', 'IN_CGE', 'Microglia')

# mcs_non_cort = which(top_3_st[1,] %in% non_cort_st)
# mcs_non_cort = mcs_non_cort[!(mcs_non_cort %in% mcs_not_npc)]
# # mcs_non_cort

# ge_genes = c('Dlx1', 'Dlx2', 'Dlx5','Lhx6', 'Gad2',  'Gad1', 'Nrxn3')
# ge_marks_in_cort = sort(apply(t(mc@mc_fp[ge_genes,!(1:ncol(mc@mc_fp) %in% mcs_non_cort)]), 1, max))         
# bad_mcs = as.numeric(names(ge_marks_in_cort[ge_marks_in_cort >= 1.25]))

# mcs_non_cort = union(mcs_non_cort, bad_mcs)

# t(top_3_st[,mcs_non_cort])


# non_cort_cells = names(mc@mc)[mc@mc %in% mcs_non_cort]

