library(metacell)
library(glue)
library(tidyverse)
library(pheatmap)
library(vroom)
library(gplots)

wd = '/home/feshap/raid/proj/mmcortex'
data_path = file.path(wd, 'scrna_data')
dirs = list.dirs(data_path)
db_path = file.path(wd, 'scdb')
if (!dir.exists(db_path)) {dir.create(db_path)}
scdb_init(db_path, force_reinit = T)
set.seed(1337)
scfigs_init("figs/")


# walk(tail(dirs, -1), function(x)
#      mcell_import_scmat_10x(mat_nm = last(str_split(x, '\\/')[[1]]), base_dir = x))

samples = gsub('mat\\.', '', scdb_ls('mat'))

create_merged_batches = function(samples) {
    samples_batch = list(list(samples[1:2]), list(samples[3:4]), list(samples[5:6]), list(samples[7:8]),
                    list(samples[9:10]), list(samples[11:12]), list(samples[13:14]))
    mats = map(samples, function(x) scdb_mat(x))
    names(mats) = samples
    mats_batch = list(list(mats[1:2]), list(mats[3:4]), list(mats[5:6]), list(mats[7:8]),
                    list(mats[9:10]), list(mats[11:12]), list(mats[13:14]))
    mats_merge = map(mats_batch, function(x) scm_merge_mats(x[[1]][[1]], x[[1]][[2]]))
    names(mats_merge) = map(samples_batch, function(x) gsub('_rep\\d', '', x[[1]][[1]]))
    walk(names(mats_merge), function(x) scdb_add_mat(x, mats_merge[x][[1]]))
    return(mats_merge)
}

mats_merge = create_merged_batches(samples)

merge_all = scm_merge_mats(mats_merge)
scdb_add_mat(id='merge', merge_all)