### Filter scATAC profiles and build McCounts object

library(metacell)
devtools::load_all("~/src/mcATAC/")
gset_genome('mm10')
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
nm <- "pl_cort_feat_peaks"

mcmd = vroom::vroom('./BonevCollab/mcmd_pl_cort.tsv')
cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC', 'iCPN_early','iCPN_late',
                      'CPN_L2-3','CPN_L5_6','iCPN/CfuPN','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(s) setNames(which(mcmd$cell_type == s), 
                                                                  rep(s, length(which(mcmd$cell_type == s))))))                         


load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/multi_mmcortex.Rda",v=T)
# load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_mmcortex.Rda",v=T)
# load("/net/mraid14/export/tgdata/users/atanay/proj/mmcortex/work0922/data/atac_clsts_annot_mmcortex.Rda",v=T)

#all intervals
head(multi_model$atac_intervs)
#those in clusters annotated as variable
# feat_peak = multi_model$atac_intervs[acn$f_var_peak,]
feat_peak = multi_model$atac_intervs


promo_peaks = readRDS(file.path(wd, 'scatac_data/prom_counts.RDS'))
batch_ids <- setNames(unlist(purrr::map(stringr::str_split(promo_peaks@Dimnames[[2]], '-'), 2)), 
                promo_peaks@Dimnames[[2]])
batch_names <- setNames(agg_id$sample_name[match(batch_ids, agg_id$library_id)], 
                promo_peaks@Dimnames[[2]])


write_sc_counts_from_fragments(fragments_file = './scatac_data/fragments_filtered.tsv.gz', overwrite =T, 
            out_dir = './data/frag_reads_28122022', 
            cell_names = promo_peaks@Dimnames[[2]][!(batch_ids %in% 1:2)], 
            genome = "mm10")

scc <- scc_read('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads_28122022/')

scmat <- scc_extract(scc = scc, intervals = feat_peak)

cs <- Matrix::colSums(scmat)

mt <- seq(4e+3, 8e+3, 100)
tlov <- sapply(mt, function(tlo) length(which(cs >= tlo & cs <= 10*tlo)))
mtm <- mt[[which.max(tlov)]]
good_cells <- which(cs >= mtm & cs <= 10*mtm)

## get multiplet ids
wd <- '~/raid/proj/mmcortex/output/AMULET/'
setwd(wd)
rep_folders <- file.path(wd, grep('rep', list.files('.'), v=T))
mbc_list <- lapply(tail(rep_folders, -2), function(x) unlist(read.delim(file.path(x, 'MultipletBarcodes_01.txt'), header = F)))
mbc_list_fix <- lapply(seq_along(mbc_list), function(i) gsub('-1$', paste0('-', i+2), mbc_list[[i]]))
mbc_vec <- unique(unlist(mbc_list_fix))
good_cells <- good_cells[!(good_cells %in% mbc_vec)]

scmat_f <- scmat[, good_cells]
scmd <- tibble::enframe(batch_names[good_cells], name = 'cell_id', value = 'batch') %>%
        tibble::column_to_rownames('cell_id')
scmd$cell_id <- rownames(scmd)
# promo_peaks@cell_metadata[good_cells,]

mat_new <- scm_new_matrix(mat = scmat_f[,cb], cell_metadata=scmd[cb,], stat_type='umi')
nm_f <- paste0(nm, '_f')
scdb_add_mat(nm_f, mat_new)