### Check runtime of scc_extract
devtools::load_all("~/src/mcATAC/")

gset_genome('mm10')
options(gmax.data.size = 1e+9)

aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
mmc_peaks <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(aaa$mcl_all))
mmc_peaks$peak_name <- peak_names(mmc_peaks, tad_based = F)

print('here1')

reik_peaks <- readr::read_csv("/net/mraid14/export/tgdata/users/aviezerl/proj/motif_reg/output/all_peaks_tor.tsv")
reik_peaks$peak_name <- peak_names(reik_peaks, tad_based = F)
nei_reik_mcc_peaks <- gintervals.neighbors(as.data.frame(reik_peaks), mmc_peaks, maxdist = 0, mindist = 0, maxneighbors = 1)
nei_mmc_reik_peaks <- gintervals.neighbors(mmc_peaks, as.data.frame(reik_peaks), maxdist = 0, mindist = 0, maxneighbors = 1)
peaks_reik_only <- dplyr::filter(reik_peaks, !(peak_name %in% nei_reik_mcc_peaks[,5]))
peaks_mmc_only <- dplyr::filter(mmc_peaks, !(peak_name %in% nei_mmc_reik_peaks[,4]))
peaks_both <- gintervals.canonic(bind_rows(nei_reik_mcc_peaks[,1:3], nei_mmc_reik_peaks[,1:3]))
peaks_both$peak_name <- peak_names(peaks_both, tad_based = F)
peaks_all <- bind_rows(peaks_reik_only, peaks_both, peaks_mmc_only) %>% gintervals.canonic

print('here2')

scc <- scc_read('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/data/frag_reads/')

print('here3')

scmat <- scc_extract(scc = scc, intervals = peaks_all)
