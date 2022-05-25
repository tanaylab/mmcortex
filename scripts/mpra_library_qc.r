library(misha)
library(metacell)
# scdb_init('scdb')

setwd('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex')
gsetroot('/home/aviezerl/mm10')

ENH_LIB_PATH = './data/mpra_library_take_3.tsv'
MCMD_PATH = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv'
PEAK_MC_PATH = './data/pl_cort_peak_mc_smoothed_mg.rds'
TRACKS_PATH = './data/motif_tracks_selected_k=64_n=64.txt'
DELS_MOTIFS_PATH = './data/motif_log_sum_exp_ENCODE_dELS_266bp_mm10.rds'
PLOTS_OUT_PATH = './figs/MPRA_QC_plots/'
TADS_PATH = './data/npc_tads_egc.rds'

enh_lib = readr::read_tsv(ENH_LIB_PATH)

mcmd = readr::read_tsv(MCMD_PATH)

cust_st_ord = c('Oligodendrocytes','Astrocytes','NSC','IPC_cyc', 'IPC','IPC_late','iCPN/CfuPN',
                'iCPN_early','iCPN_late','CPN_L2-3','CPN_L5_6','iCfuPN','SCPN','CthPN')
cust_mc_ord_st = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))))
cust_mc_ord_st_ord_md = unlist(lapply(cust_st_ord, function(st) setNames(sort(mcmd$mc[mcmd$st == st]), 
                                                                  rep(st, length(which(mcmd$st == st))))[order(mcmd$mean_day[mcmd$st == st])]))

native_seqs_to_sample = do.call('rbind', 
                                lapply(unique(enh_lib$type), function(typ) {
    return(do.call('rbind', lapply(unique(enh_lib$pattern), function(patt) {
        df_samp = dplyr::filter(enh_lib, type == typ & pattern == patt & is.na(source))
        return(df_samp[sort(sample(nrow(df_samp), min(nrow(df_samp), 5))),])
    })))
})
)
ns_rn = gsub(' ', '', apply(dplyr::select(dplyr::mutate(native_seqs_to_sample, start = start - 117, end = end + 117), chrom, start, end), 1, paste0, collapse = '-'))
pmc = readRDS(PEAK_MC_PATH)
pmc_norm = t(apply(pmc, 1, function(x) x/sum(x)))
ord = cust_mc_ord_st_ord_md

options(repr.plot.width = 14)

purrr::walk(sample(seq_along(ns_rn), 6), function(x,n,i) {
    barplot(pmc_norm[x[[i]],ord], col = mcmd$color[ord], main = paste(c(n[i,], x[[i]])))
}, x = ns_rn, n = native_seqs_to_sample[,c('type', 'pattern')])

trk_s = unlist(read.delim(TRACKS_PATH))
dels_motifs = readRDS(DELS_MOTIFS_PATH)
dels_motifs = dels_motifs[,c('chrom', 'start', 'end', trk_s[trk_s %in% colnames(dels_motifs)])]
motif_ecdfs = apply(subset(dels_motifs, select = -c(chrom, start, end)), 2, function(x) ecdf(x))
dels_to_plot = enh_lib[enh_lib$source %in% native_seqs_to_sample$seq_name,]
dels_motifs_raw = gextract(trk_s[trk_s %in% colnames(dels_motifs)], intervals = dels_to_plot, iterator = 10)
dels_motifs_raw = dels_motifs_raw[with(dels_motifs_raw, order(chrom, start, end)),]
dels_motifs = tgs_matrix_tapply(t(subset(dels_motifs_raw, select = -c(chrom, start, end, intervalID))), dels_motifs_raw$intervalID, function(x) log(sum(exp(x))))
rownames(dels_motifs) = enh_lib$seq_name[enh_lib$source %in% native_seqs_to_sample$seq_name]
seqs_motifs_raw = gextract(trk_s[trk_s %in% colnames(dels_motifs)], intervals = native_seqs_to_sample, iterator = 10)
seqs_motifs_raw = seqs_motifs_raw[with(seqs_motifs_raw, order(chrom, start, end)),]
seqs_motifs = tgs_matrix_tapply(t(subset(seqs_motifs_raw, select = -c(chrom, start, end, intervalID))), seqs_motifs_raw$intervalID, function(x) log(sum(exp(x))))
rownames(seqs_motifs) = native_seqs_to_sample$seq_name
seqs_q = sapply(1:ncol(seqs_motifs), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = subset(seqs_motifs))
rownames(seqs_q) = rownames(seqs_motifs)
dels_q = sapply(1:ncol(dels_motifs), function(x, i) {motif_ecdfs[[i]](x[,i])}, x = dels_motifs)
rownames(dels_q) = rownames(dels_motifs)
seqs_m_z = apply(seqs_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})
dels_m_z = apply(dels_q, 2, function(x) {v = -log2(1-x); v[is.infinite(v)] = max(v[!is.infinite(v)]); return(v)})
seqs_dels_map = as.data.frame(do.call('rbind', lapply(native_seqs_to_sample$seq_name, function(x) which(dels_to_plot$source == x))))
rownames(seqs_dels_map) = native_seqs_to_sample$seq_name
colnames(seqs_dels_map) = 1:ncol(seqs_dels_map)

plot_motif_pair_energy = function(ind1, ind2) {
    vec1 = seqs_motifs[ind1,]
    vec2 = dels_motifs[ind2,]
    plot(vec1, vec2, pch = 16, cex = 2, main = rownames(dels_motifs)[[ind2]],
    xlab = 'Mean energies per motif of native enhancer',
        ylab = 'Mean energies per motif of shadow enhancer')
    abline(0,1,lty='dashed',col='red')
#     if (ind2%%3 == 2) title(rownames(seqs_motifs)[[ind1]])
}

plot_motif_pair_lq = function(ind1, ind2) {
    vec1 = seqs_m_z[ind1,]
    vec2 = dels_m_z[ind2,]
    plot(vec1, vec2, pch = 16, cex = 2, 
        ylab = paste0('-log2(1-q) per motif of ', purrr::map(stringr::str_split(rownames(dels_motifs)[[ind2]], '\\.'), 4), ' shadow enh.'),
        xlab = paste0('-log2(1-q) per motif of native enh.'))
    abline(0,1, lty = 'dashed', col = 'darkgreen', lwd = 2)
    abline(-4,1, lty = 'dashed', col = 'red', lwd = 2)
    abline(4,1, lty = 'dashed', col = 'red', lwd = 2)
    if (ind2%%3 == 2) title(rownames(seqs_motifs)[[ind1]])
}

plot_motif_row = function(i, path) {
    png(paste0(path,native_seqs_to_sample$seq_name[[i]], '_motif_energy.png'), 
                            w = 3000, h = 1000, res = 150)
    par(mfrow = c(1,ncol(seqs_dels_map)), mar = c(5,5,6,4), cex.lab = 2, cex.main = 2)
    sapply(seqs_dels_map[i,], function(j) plot_motif_pair_energy(i, j))
    dev.off()
    png(paste0(path, native_seqs_to_sample$seq_name[[i]], '_motif_log_q.png'), 
                            w = 3000, h = 1000, res = 150)
    par(mfrow = c(1,ncol(seqs_dels_map)), mar = c(5,5,6,4), cex.lab = 2, cex.main = 2)
    sapply(seqs_dels_map[i,], function(j) plot_motif_pair_lq(i, j))
    dev.off()
}


tads = readRDS(TADS_PATH)

if (!dir.exists(PLOTS_OUT_PATH)) dir.create(PLOTS_OUT_PATH)

n = native_seqs_to_sample[,c('type', 'pattern')]
grp_vec = sort(sapply(1:3, rep, nrow(mcmd)))
plot_everything = function(i) {
    st = native_seqs_to_sample$type[[i]]
    st_path = paste0(PLOTS_OUT_PATH, st, '/')
    if (!dir.exists(st_path)) dir.create(st_path)
    pat = native_seqs_to_sample$pattern[[i]]
    pat_path = paste0(st_path, pat, '/')
    if (!dir.exists(pat_path)) dir.create(pat_path)
    enh_path = paste0(pat_path, ns_rn[[i]], '/')
    if (!dir.exists(enh_path)) dir.create(enh_path)
    
    ## Plot native enhancer accessibility across metacells
    png(paste0(enh_path, '/', paste0(c(n[i,], ns_rn[[i]]), collapse = '_'), '.png'), w = 2400, h = 1000, res = 150)
    barplot(pmc_norm[ns_rn[[i]],ord], col = mcmd$color[ord], main = paste(c(n[i,], ns_rn[[i]])))
    dev.off()
    
    ## Plot motif energies
    plot_motif_row(i, enh_path)
    
    ## Boxplot TADs
    tad_cat_vec = as.numeric(sapply(seqs_dels_map[i,], function(k) {
        tadID = gintervals.neighbors(dels_to_plot[k,], tads[,1:5])$tadID; 
        return(tads[tads$tadID == tadID,tail(1:ncol(tads), nrow(mcmd))])
    }))
    png(paste0(enh_path, '/', paste0(c(native_seqs_to_sample$seq_name[i], 'TAD_boxplot'), collapse = '_'), '.png'), 
        w = 1000, h = 1000, res = 150)
    boxplot(tad_cat_vec ~ grp_vec, 
            xlab = 'TAD type', ylab = 'Sum RNA per metacell',
            names = purrr::map(stringr::str_split(dels_to_plot$seq_name[unlist(seqs_dels_map[i,])], '\\.'), 4))
    title(paste(c('TADs of shadow enhancers of ', native_seqs_to_sample$seq_name[i])))
    dev.off()
}

sapply(1:nrow(native_seqs_to_sample), plot_everything)

