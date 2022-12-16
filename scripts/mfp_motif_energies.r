### Get mfp motif energies because prego isn't working on amosbase

devtools::load_all("~/src/mcATAC/")
library(prego)
gset_genome('mm10')
load('./output/sequence_modeling/motifs_to_model_proximal_atac.rda')
amd <- all_motif_datasets()
mca <- readRDS('./output/mcatac/mmcortex_mcatac_feat_peaks_test.rds')
vm <- gextract_pwm(intervals = dplyr::select(mca@peaks, chrom, start, end, peak_name), dataset = amd, motifs = mfp)
vmm <- subset(vm, select = -c(chrom, start, end, peak_name))
rownames(vmm) <- vm$peak_name
vmmn <- apply(vmm, 2, function(x) {y <- x; y[x < -1e+6] <- min(x[x > -1e+6], na.rm = T); return(y)})
save(vmmn, file = './output/sequence_modeling/mmcortex_mfp_energies.rda')