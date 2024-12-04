library(metacell)
library(lpsymphony)
library(pheatmap)
devtools::load_all("~/src/mcATAC")
SEED = 1337
set.seed(SEED)
wd = '/net/mraid20/export/tgdata/users/atanay/proj/mmcortex/work0922/'
setwd(wd)
scdb_init(base_dir = './scdb', force_reinit = T)
scfigs_init('./figs')
save_pheatmap_png <- function(x, filename, width=2500, height=2500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

K = 50
NUM_EDGE_NODE <- 10

load("data/multi_mmcortex.Rda",v=T)
load("data/atac_clsts_mmcortex.Rda",v=T)
load("data/atac_clsts_annot_mmcortex.Rda",v=T)

#all intervals
head(multi_model$atac_intervs)
#those in clusters annotated as variable
feat_peak = multi_model$atac_intervs[acn$f_var_peak,]

scc <- scc_read("/net/mraid20/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads/")

scmat <- scc_extract(scc =scc, intervals = feat_peak)
