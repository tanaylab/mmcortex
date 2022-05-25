library(misha)
library(pheatmap)
gsetroot('/home/aviezerl/mm10')
setwd('/home/feshap/raid/proj/mmcortex/')
mm_esd = readRDS('./data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds')
mm_esd = tidyr::drop_na(as.data.frame(mm_esd))
K = 64

save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#dir.create('./figs/TF_motif_track_selection')
#dir.create('./figs/TF_motif_track_selection/hists')

#sapply(colnames(mm_esd)[4:(ncol(mm_esd)-1)], function(nm) {cni = purrr::map(stringr::str_split(purrr::map(stringr::str_split(nm, '\\.'), 2), '_'), 1); png(paste0('./figs/TF_motif_track_selection/hists/', cni, '_', nm, '.png')); hist(mm_esd[,nm], 50, main = nm); dev.off()})

vec_cuts = sapply(colnames(mm_esd)[4:(ncol(mm_esd)-1)], function(nm) {vec = mm_esd[,nm]; brks = seq(min(vec), max(vec), l=51); vec_cut = cut(vec, breaks = brks); return(vec_cut)})
vec_cuts = tidyr::drop_na(data.frame(vec_cuts))

ent_trk = apply(vec_cuts, 2, function(x) {sum_cut = as.numeric(table(x)); sum_norm = sum_cut/sum(sum_cut); ent = sum(-log(sum_norm)*sum_norm); return(ent)})
# print(names(ent_trk))
# ent_trk_nms = sapply(sapply(names(ent_trk), function(x) gsub('_max\\.?\\d?|\\.{3,4}\\d?', '', x)), function(y) gtrack.ls(y)[[1]])
ent_trk_nms = sapply(names(ent_trk), function(y) {v = gtrack.ls(y); return(ifelse(length(v) > 0, v[[1]], NA))})

names(ent_trk) = ent_trk_nms 
readr::write_tsv(data.frame(cbind(names(ent_trk), ent_trk)), './data/track_entropy_mm10.txt')
cor_motif = tgs_cor(as.matrix(mm_esd[,grep('jolma_10bp|motifs_10bp|cis_bp_10bp|jaspar_10bp', colnames(mm_esd))]), spearman=T)
print(colnames(cor_motif)[grep('jolma_10bp|motifs_10bp|cis_bp_10bp|jaspar_10bp', colnames(cor_motif), inv=T, v=T)])

motif_km  = tglkmeans::TGL_kmeans(cor_motif, k=K)
col_annot = data.frame(motif_km$cluster)
rownames(col_annot) = rownames(cor_motif)
colnames(col_annot) = 'cl'
km_col_key = data.frame(cbind(1:K, sample(grep('white|grey|gray', colors(), inv=T, v=T), K)))
colnames(km_col_key) = c('cl', 'color')

ann_colors = list('cl' = setNames(gplots::col2hex(km_col_key$color), km_col_key$cl))
ord = order(motif_km$cluster)
p =  pheatmap(cor_motif[ord,ord], fontsize = 8, annotation_colors = ann_colors, annotation_row = col_annot, annotation_col = col_annot, cluster_rows=F, cluster_cols=F)
save_pheatmap_png(p, './figs/TF_motif_track_selection/cor_motif_mm10.png', h = 4000, w = 4000)

trk_select = sapply(sort(unique(motif_km$cluster)), function(i) rownames(cor_motif)[which(motif_km$cluster == i)[[1]]])
print(trk_select)
# nms_match = sapply(sapply(trk_select, function(x) gsub('_max|\\.{3,4}\\d?', '', x)), function(y) gtrack.ls(y)[[1]])
nms_match = sapply(trk_select, function(y) {v = gtrack.ls(y); return(ifelse(length(v) > 0, v[[1]], NA))})
print(length(nms_match))
nms_match = nms_match[!duplicated(nms_match)]
print(length(nms_match))
write(nms_match[!duplicated(nms_match)], glue::glue('./data/motif_tracks_selected_k={K}_n={length(nms_match)}.txt'))

print(table(unlist(purrr::map(stringr::str_split(purrr::map(stringr::str_split(trk_select, '\\.'), 2), '_'), 1))))
