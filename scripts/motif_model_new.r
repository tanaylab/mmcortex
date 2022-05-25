library(dplyr)
library(parallel)
library(misha)
library(tgutil)
library(glmnet)
gsetroot('/home/aviezerl/mm10')
SEED = 1337
set.seed(SEED)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
setwd(wd)


prox_enh = unlist(read.delim('./data/mmcortex_proximal_enhancers_mm10.txt', header = F))
dist_enh = unlist(read.delim('./data/mmcortex_distal_enhancers_mm10.txt', header = F))
mcmd = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mcmd_pl_cort.tsv')
pmc = readRDS('./data/pl_cort_peak_mc_smoothed_mg.rds')
pmc_mot_filt = readRDS('./data/peak_mc_motifs_200bp_filt.rds')
ks_results = readRDS('./data/pmc_vs_rg_ks_results_200bp.rds')

print('Imported everything')

pmc_filt = pmc[rownames(pmc) %in% c(prox_enh, dist_enh),]
pmc_norm = t(apply(pmc_filt, 1, function(x) x/sum(x)))
colnames(pmc_norm) = 1:ncol(pmc_norm)

# names(ks_results) = head(tail(colnames(pmc_mot_filt), -3), -1)
# ks_d_mat = sapply(ks_results, function(x) sapply(x, function(y) y$statistic))
# # colnames(ks_d_mat) = colnames(subset(pmc_mot_filt, select = -c(chrom, start, end, intervalID)))
# colnames(ks_d_mat) = grep('chrom|start|end', colnames(pmc_mot_filt), inv=T, v=T)
# rownames(ks_d_mat) = 1:nrow(ks_d_mat)

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  return(data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  ))
}

# qcut = cut(log2(apply(ks_d_mat, 2, function(x) mean(x))), breaks = seq(-7.5, -1, l=40))
# covvec = log2(apply(ks_d_mat, 2, function(x) sd(x)/mean(x)))
# tq = table(qcut)

# mcs = c(sample(mcmd$mc[mcmd$st == 'NSC'], 1),
#         sample(mcmd$mc[mcmd$st == 'IPC'], 1),
#         sample(mcmd$mc[mcmd$st == 'CPN_L5_6'], 1))

# mcs = c(23, 122, 434)
# param_grid = expand.grid(list('mean_thresh' = c(-15, -10, -5), 'q_thresh' = c(0, 0.25, 0.5,0.75), 'ALPH' = c(0,1), 'mci' = mcs, 'cis_f' = FALSE))
# print(head(param_grid))

grid_search_motif_alphas = function(i) {
    print(round(i/nrow(param_grid), digits = 2))
    mt = param_grid[i,'mean_thresh']
    qt = param_grid[i,'q_thresh']
    ALPH = param_grid[i,'ALPH']
    cf = param_grid[i,'cis_f']
    mci = param_grid[i,'mci']
    motif_select = unlist(sapply(tail(unique(names(tq)), mt), 
                             function(nm) colnames(ks_d_mat)[which(covvec > quantile(covvec[qcut == nm], qt) & qcut == nm)]))
    if (cf) {motif_select = grep('cis_bp', motif_select, inv=T, v=T)}
    pmc_mot_ff = subset(pmc_mot_filt, select = motif_select)

    rownames(pmc_mot_ff) = gsub(' ', '', apply(pmc_mot_filt[,1:3], 1, paste0, collapse = '-'))
    inds = 1:nrow(pmc_mot_ff)
    test_inds = sample(inds, round(0.1*length(inds)))
    train_inds = inds[!(inds %in% test_inds)]
    x_train = as.matrix(pmc_mot_ff[train_inds,])
    x_test = as.matrix(pmc_mot_ff[test_inds,])
    y_train = pmc_norm[train_inds,mci]
    y_test = pmc_norm[test_inds,mci]
    lambdas <- 10**seq(-10, -5, by = .5)
    cv_ridge <- cv.glmnet(x_train, y_train, alpha = ALPH, lambda = lambdas)
    optimal_lambda <- cv_ridge$lambda.min
    ridge_reg_opt = glmnet(x_train, y_train, lambda = optimal_lambda, alpha = ALPH, family = 'gaussian')
    predictions_train <- predict(ridge_reg_opt, s = optimal_lambda, newx = x_train)
    predictions_test <- predict(ridge_reg_opt, s = optimal_lambda, newx = x_test)

    ev = eval_results(y_test, predictions_test, x_test)
    return(list('q_thresh' = qt, 'mean_thresh' = mt, 'alpha' = ALPH, 'mci' = mci, 
                'eval_results' = ev, 'cis_f' = cf, 'optimal_lambda' = optimal_lambda))
}

# prm_check = mclapply(1:nrow(param_grid), grid_search_motif_alphas,  mc.cores = 24)

# saveRDS(prm_check, './data/motif_model_param_grid_search_take_2.rds')

params = list('mean_thresh' = -12, 'q_thresh' = 0.35, 'ALPH' = 0.5, 'cis_f' = FALSE)
mt = params[['mean_thresh']]
qt = params[['q_thresh']]
ALPH = params[['ALPH']]
cf = params[['cis_f']]


# motif_select = unlist(sapply(tail(unique(names(tq)), mt), 
#                             function(nm) colnames(ks_d_mat)[which(covvec > quantile(covvec[qcut == nm], qt) & qcut == nm)]))
# if (cf) {motif_select = grep('cis_bp', motif_select, inv=T, v=T)}

# motif_select = grep('kmers\\.{')

# motif_select = c('kmers.C_tmp', 'kmers.CG_tmp')
print('Getting GC content...')
pmc_gc_content = gextract('(seqG + seqC)/2', intervals = pmc_mot_filt[,1:3], 
        iterator = pmc_mot_filt[,1:3], colnames = 'gc_content')
pmc_gc_content = pmc_gc_content[order(as.numeric(pmc_gc_content$intervalID)),]

print('Check if gc content vec is ordered:')
print(which(pmc_gc_content$intervalID != sort(pmc_gc_content$intervalID)))

# nei_pmc_cg2 = readRDS('./data/CGs_in_peak_mc_coords.rds')
nei_pmc_cg2 = readRDS('./data/CGs_in_peak_mc_coords_200bp.rds')
cpg_content = table(nei_pmc_cg2$intervalID1)
cpg_content = setNames(as.numeric(cpg_content), names(cpg_content))
cpg_content[as.character(1:nrow(pmc_mot_filt))[!(1:nrow(pmc_mot_filt) %in% as.numeric(names(cpg_content)))]] = 0
cpg_content = cpg_content[order(as.numeric(names(cpg_content)))]

motif_select = unlist(read.delim('./data/motifs_to_model_250.txt', header = F))
motif_select_match = unique(unlist(lapply(motif_select, function(x) grep(glue::glue("{substr(x, 1, 35)}|{x}_tmp_200bp"), colnames(pmc_mot_filt), v=T))))

# pmc_gc_rn = gsub(' ', '', apply(pmc_gc_content[,1:3], 1, paste0, collapse = '-'))
# motif_select = grep('kmers', colnames(pmc_mot_filt), v=T)
pmc_mot_ff = subset(pmc_mot_filt, select = motif_select_match)
# pmc_mot_ff = pmc_mot_filt
rownames(pmc_mot_ff) = gsub(' ', '', apply(pmc_mot_filt[,1:3], 1, paste0, collapse = '-'))
pmc_mot_ff$gc_content = pmc_gc_content$gc_content
pmc_mot_ff$cpg_content = cpg_content

# pmf_m = apply(subset(pmc_mot_ff, select = -c(chrom, start, end)), 2, mean, na.rm=T)
# pmf_sd = apply(subset(pmc_mot_ff, select = -c(chrom, start, end)), 2, sd, na.rm=T)
# pmc_mot_stzd = sapply(grep('chrom|start|end$', colnames(pmc_mot_ff), inv=T), function(x,i) {(x[,i] - pmf_m[[i-3]])/pmf_sd[[i-3]]}, x = pmc_mot_ff)
pmc_mot_stzd = apply(pmc_mot_ff, 2, function(x) (x - mean(x))/sd(x))

fit_mc = function(mci) {
    print(mci)
    inds = 1:nrow(pmc_mot_stzd)
    test_inds = sample(inds, round(0.1*length(inds)))
    train_inds = inds[!(inds %in% test_inds)]
    x_train = as.matrix(pmc_mot_stzd[train_inds,])
    x_test = as.matrix(pmc_mot_stzd[test_inds,])
    y_train = pmc_norm[train_inds,mci]
    y_test = pmc_norm[test_inds,mci]
    # lambdas <- 10**seq(-10, -7, by = 1)
    lambda <- 10**-10
    # cv_ridge <- cv.glmnet(x_train, y_train, alpha = ALPH, lambda = lambdas)
    # optimal_lambda <- cv_ridge$lambda.min
    ridge_reg_opt = glmnet(x_train, y_train, lambda = lambda, alpha = ALPH, family = 'gaussian')
    predictions_train <- predict(ridge_reg_opt, s = lambda, newx = x_train)
    predictions_test <- predict(ridge_reg_opt, s = lambda, newx = x_test)
    ev = eval_results(y_test, predictions_test, x_test)
    return(list('mci' = mci, 'Rsquare' = ev$Rsquare, #'optimal_lambda' = optimal_lambda, 
          'model' = ridge_reg_opt, 'train_inds' = train_inds, 'test_inds' = test_inds))
}

mc_models = mclapply(1:ncol(pmc_norm), fit_mc, mc.cores = 40)

saveRDS(mc_models, glue::glue('./data/mc_atac_models/stdzd_vars_200bp_motif_select_250_qt={qt}_mt={mt}_ALPH={ALPH}.rds'))