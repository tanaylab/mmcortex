library(metacell)
library(dplyr)
set.seed(1337)
doMC::registerDoMC(20)
wd = '/net/mraid14//export/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)

scdb_init('scdb', force_reinit=T)
scfigs_init('figs')

# mat = scdb_mat('cort5')
# mc = scdb_mc('cort5')

mat = scdb_mat('peak_cort')
mc = scdb_mc('peak_cort')

mct = scdb_mctnetwork('cort5_comp')
mcmd = vroom::vroom('./BonevCollab//mcmd_cort_5_iter1.tsv')

flows = mct@network[mct@network$flow > 0 & mct@network$time1 > 0 & 
        mct@network$time2 < 7 & mct@network$time1 != mct@network$time2,]

mc_parents = function(f, mc, t, mc_p) {
    if (t >= min(f$time1))
        new_mcs = f$mc1[f$time2 == t & f$mc2 %in% mc & f$mc1 != -2]
        if (length(new_mcs > 0)) {
            df2 = t(data.frame(setNames(rbind(new_mcs, rep(as.numeric(t-1), length(new_mcs))), c('mc', 't'))))
            mc_p = rbind(mc_p, mc_parents(f, new_mcs, t-1, df2))
        }
    return(mc_p)
}

get_mc_tree = function(mc, t, flows, umc_list) {
    df1 = t(data.frame(setNames(c(mc,t), c('mc', 't'))))
    res = unique(data.frame(apply(mc_parents(flows, mc, t, df1), 2, as.numeric)))
    res[res$mc %in% umc_list,]
    return(res)
}
mc_day = purrr::map_dfr(unique(mat@cell_metadata$day), function(u) 
        tapply(rownames(mat@cell_metadata)[mat@cell_metadata$day == u], 
           mc@mc[rownames(mat@cell_metadata)[mat@cell_metadata$day == u]], length)
                        ) %>% data.frame
mc_day[is.na(mc_day)] = 0
rownames(mc_day) = unique(mat@cell_metadata$day)
colnames(mc_day) = gsub('X', '', colnames(mc_day))
cnn = as.numeric(colnames(mc_day))
cn_max = 1:max(cnn)
df_new = data.frame(matrix(0, ncol = length(cn_max[!(cn_max %in% cnn)]), nrow = nrow(mc_day)))
colnames(df_new) = as.character(cn_max[!(cn_max %in% cnn)])
df_new = cbind(mc_day, df_new)
mc_day = df_new[,order(as.numeric(colnames(df_new)))]

umc_list = unique(mc@mc)
mc_trees = lapply(sort(unique(flows$mc2[flows$time2 == max(flows$time2)])), 
                function (mci) get_mc_tree(mci, max(flows$time2), flows, umc_list))
names(mc_trees) = sort(unique(flows$mc2[flows$time2 == max(flows$time2)]))
mc_trees = lapply(mc_trees, function(tree) {return(tree[apply(tree, 1, function(x) mc_day[as.character(12 + x[[2]]), as.character(x[[1]])] > 0),])})
# mc_trees = readRDS('./data/mc_trees.RDS')

print('done getting trees')
mc_sc_list = lapply(mc_trees, function(tree) setNames(apply(tree, 1, function(x) 
                    rownames(mat@cell_metadata)[mat@cell_metadata$day == 12+x[[2]] & 
                                            rownames(mat@cell_metadata) %in% names(mc@mc)[mc@mc == x[[1]]]
                                            ]), 
                                            tree$mc)
                    )
names(mc_sc_list) = names(mc_trees)
mc_sc_list = lapply(mc_sc_list, function(x) x[sapply(x, function(y) length(y) > 0)])
# mc_trees = purrr::map2(mc_trees, mc_sc_list, function(tree,sc_list) return(x[sapply(y, function(yi) length(yi) > 0),]))
saveRDS(mc_trees, './data/mc_trees_peak.RDS')

sc_named_list = lapply(mc_sc_list, function(x) lapply(seq_along(x), function(n,x,i) {
                                                                    setNames(rep(n[[i]], length(x[[i]])), x[[i]])
                                                                                    }, x=x, n = 1:length(x)) %>% unlist
                        )
names(sc_named_list) = names(mc_trees)
# mc_sc_list = readRDS('./data/mc_sc_list.RDS')
print('done getting sc_lists')

saveRDS(mc_sc_list, './data/mc_sc_list_peak.RDS')
saveRDS(sc_named_list, './data/sc_named_list_peak.RDS')

# tfoi = read.delim('./data/tfoi.txt') %>% unlist %>% sort %>% as.character
# tfoi = c('4930506M07Rik', tfoi)
# tfoi = grep('Shtn1', tfoi, invert=T, v=T)

rs = Matrix::rowMeans(mat@mat)
top_enh = head(rownames(mat@mat)[order(rs, decreasing = T)], 10000)

mc_trees = purrr::map2(mc_trees, sc_named_list, function(tree, sc_list) {
                        print(head(tree, 2))
                        cols_k = colnames(mat@mat) %in% names(unlist(sc_list))
                        mat_filt = mat@mat[top_enh,cols_k];
                        cs = Matrix::colSums(mat@mat[,cols_k])
                        mat_filt = mat_filt/cs
                        # quantile(mat_filt, seq(0, 1, 0.2))
                        return(cbind(tree, tgs_matrix_tapply(mat_filt[,names(sc_list)], sc_list, mean)))
                        }
)
print('done computing tf flows')

saveRDS(mc_trees, './data/mc_trees_peak.RDS')
