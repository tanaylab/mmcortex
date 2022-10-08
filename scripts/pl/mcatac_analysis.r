devtools::load_all("~/src/mcATAC/")
library(pheatmap)
setwd('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex')
my_genome <- "mm10"
gset_genome(my_genome)
options(gmax.data.size = 1e+9)
options(mcatac.parallel = 10)

mcmd <- readr::read_tsv('./BonevCollab/mcmd_pl_cort.tsv')
marg_track <- gtrack.ls('mmcortex.marg')
qt <- 0.875
max_peak_size <- 500
min_umis <- 8


peaks_split <- call_peaks(marg_track, 
                            split_peaks = TRUE, 
                            quantile_thresh = qt, 
                            min_umis = min_umis, 
                            max_peak_size = max_peak_size, 
                            genome = my_genome)
peaks_split$peak_name <- peak_names(peaks_split)

mcc <- mcc_read("./output/mcatac/mmcortex_mcc/")

mca <- mcc_to_mcatac(mc_counts=mcc, peaks=peaks_split, metadata = mcmd)

mx <- apply(mca@egc, 1, max)
mn <- apply(mca@egc, 1, min)

mcaf <- subset_peaks(mca, mca@peaks[which(log2(mx - mn)) >= 4,])

mcaf_egc_z <- (mcaf@egc - rowMeans(mcaf@egc))/matrixStats::rowSds(mcaf@egc)
mcaf_z_km <- tglkmeans::TGL_kmeans(mcaf_egc_z, k=48, seed = 1337)

color_key <- unique(mcmd[,c('cell_type', 'color')])
col_annot <- as.data.frame(mcmd[,'cell_type'])
ann_colors <- list(cell_type = tibble::deframe(color_key))

diversity <- function(x,q) {
    dm <- dim(x)
    if (length(dm) <= 2) {
        if (min(dim(x)) == 1) {
            if (sum(x) != 1) {
                cli::cli_warn("Vector should be normalized to 1. Normalizing")
                x <- x/sum(x)
                div <- sum(x**q)
            }
        } else {
            if (!all(colSums(x) == 1)) {
                cli::cli_warn("Matrix columns should be normalized to 1. Normalizing")
                x <- t(t(x)/colSums(x))
            }
            div <- rowSums(x**q)
        }
        div <- div**(1/(1-q))
        return(div)
    } else {
        stop("Input must be a vector or matrix")
    }
}