library(glmnet)
library(pheatmap)
devtools::load_all("~/src/mcATAC/")
library(prego)
doMC::registerDoMC(60)
library(matrixStats)
setwd("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex")

## Get the data
load('./output/sequence_modeling/data_list.rda')
xdata <- data_list[['xdata']]
ydata <- data_list[['ydata']]
xtrain <- data_list[['xtrain']]
ytrain <- data_list[['ytrain']]
xtest <- data_list[['xtest']]
ytest <- data_list[['ytest']]
test_inds <- data_list[['test_inds']]
# motifs_mod <- rownames(pltmt_rep)[p_plt$tree_row$order[1:nrow(pltmt_rep)]]
motifs_mod <- data_list[['motifs_mod']]

train_inds <- (1:ncol(xdata))[!(1:ncol(xdata) %in% test_inds)]

misha.ext::gset_genome('mm10')
load("~/raid/proj/mmcortex/output/mcatac/peak_indices_var_atac_clust.rda")
aaa <- readRDS("/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/output/mcatac/mmc_atac_mcls.rds")
mcl_all_norm <- t(t(aaa$mcl_all)/Matrix::colSums(aaa$mcl_all))

seq_coords <- misha.ext::convert_10x_peak_names_to_misha_intervals(rownames(mcl_all_norm))
seq_coords$peak_name <- peak_names(seq_coords, tad_based = F)
coords_exp <- dplyr::mutate(seq_coords, start = start - 100, end = end + 100) %>% 
                select(chrom, start, end) 
seqs_all <- unlist(sapply(1:ceiling(nrow(coords_exp)/1000), function(i) gseq.extract(coords_exp[(1+(i-1)*1000):(min(i*1000, nrow(coords_exp))),])))

# Identify promoter-proximal peaks from each group
tss <- gintervals.load('intervs.global.tss')
nei_sp_atac_prom <- gintervals.neighbors(seq_coords[match(sp_varp_atac, seq_coords$peak_name),], tss, maxneighbors = 1, mindist = -1e3, maxdist = 1e3)
sp_varp_atac_prom <- sp_varp_atac[sp_varp_atac %in% nei_sp_atac_prom$peak_name]
sp_varp_atac_dist <- sp_varp_atac[!(sp_varp_atac %in% nei_sp_atac_prom$peak_name)]

int_select <- coords_exp[seq_coords$peak_name %in% sp_varp_atac_dist,]
rownames(int_select) <- seq_coords$peak_name[seq_coords$peak_name %in% sp_varp_atac_dist]
int_select2 <- int_select[order(match(rownames(int_select), sp_varp_atac_dist)),]
# mex_ag <- rownames(res_all)[rowMaxs(res_all) >= 0.25]
load('./output/sequence_modeling/motifs_for_lasso.rda')
vals <- gextract_pwm(intervals = int_select2,motifs = names(motifs_to_take), prior = 0.01)
vm <- t(as.matrix(subset(vals, select = -c(chrom, start, end))))
vmn <- vm - rowMaxs(vm)
vmnz <- t(apply(vmn, 1, function(x) (x - mean(x))/sd(x)))
run_model <- function(motifs = motifs_mod, amcs = nsc_amc, alpha = 0, prego_mat = NULL,
                        xdata = xdata, ydata = ydata, train_inds = train_inds, test_inds = test_inds) {
    ## get data
    # xdata <- vm[motifs,]
    # xdata <- xdata
    # ydata <- mcl_all_norm[sp_varp_atac_dist,amcs]
    print('here0')
    xdata <- as.matrix(xdata)
    xtrain <- t(xdata[,train_inds])
    print('here01')
    xtest <- t(xdata[,test_inds])
    print('here02')
    ydata <- as.matrix(ydata)
    ytrain <- ydata[train_inds,]
    print('here03')
    ytest <- ydata[test_inds,]
    print(dim(ytest))
    print(dim(xtest))
    ## Run model
    if (!is.null(prego_mat) && nrow(prego_mat) == ncol(xdata)) {
        print('here1')
        xtrain <- bind_cols(xtrain, prego_mat[which(!(1:nrow(prego_mat) %in% test_inds)),])
        xtest <- bind_cols(xtest, prego_mat[test_inds,])
        print('here2')
        # res_glm <- plyr::llply(1:ncol(ytrain), function(i) {
        res_glm <- plyr::llply(1:ncol(prego_mat), function(i) {
                                        nnz_inds <- which(ytrain[,i] != 0)
                                        cv.glmnet(as.matrix(xtrain[nnz_inds,c(1:nrow(xdata), nrow(xdata)+i)]), 
                                                            ytrain[nnz_inds,i], 
                                                            intercept = F, 
                                                            alpha = alpha)
        }, .parallel = T)
        print('here3')
        optimal_lambdas <- lapply(res_glm, function(x) x$lambda.min)
        predictions_train <- do.call('cbind', plyr::llply(1:ncol(prego_mat), function(i) 
                                        predict(res_glm[[i]], 
                                                s = optimal_lambdas[[i]], 
                                                newx = as.matrix(xtrain[,c(1:nrow(xdata), nrow(xdata)+i)])), 
                                                    .parallel=T))
        print('hereA')
        predictions_test <- do.call('cbind', plyr::llply(1:ncol(ytest), function(i) predict(res_glm[[i]], 
                                                                s = optimal_lambdas[[i]], 
                                                                newx = as.matrix(xtest[,c(1:nrow(xdata), nrow(xdata)+i)])), 
                                                                .parallel=T))
        print('hereB')
    } else {
        res_glm <- plyr::llply(1:ncol(ytrain), function(i) {
            nnz_inds <- which(ytrain[,i] != 0)
            cv.glmnet(xtrain[nnz_inds,], ytrain[nnz_inds,i], intercept = F, alpha = alpha)
        }
                    , .parallel = T)
        optimal_lambdas <- lapply(res_glm, function(x) x$lambda.min)
        predictions_train <- do.call('cbind', plyr::llply(1:ncol(ytrain), function(i) 
                                        predict(res_glm[[i]], s = optimal_lambdas[[i]], newx = xtrain), 
                                                    .parallel=T))
        print('here6')
        print(dim(xtest))
        print(dim(ytest))
        predictions_test <- do.call('cbind', plyr::llply(1:ncol(ytest), function(i) predict(res_glm[[i]], 
                                                                    s = optimal_lambdas[[i]], 
                                                                    newx = xtest), .parallel=T))
    }
    ## Evaluate model
    eval_results <- function(true, predicted, df) {
        SSE <- sum((predicted - true)^2)
        SST <- sum((true - mean(true))^2)
        R_square <- 1 - SSE / SST
        RMSE = sqrt(SSE/nrow(df))
        # Model performance metrics
        data.frame(RMSE = RMSE, Rsquare = R_square)
    }
    print('here5')
    res_df <- do.call('rbind', plyr::llply(1:ncol(ytrain), 
                    function(i) eval_results(ytrain[,i], predictions_train[[i]], xtrain)))
    
    print('here4')
    
    print('here7')
    coef_mat <- as.matrix(do.call('cbind', lapply(res_glm, coefficients)))
    colnames(coef_mat) <- colnames(ytrain)
    resl <- list(predictions_train = predictions_train, predictions_test = predictions_test,
                     res_df = res_df, coef_mat = coef_mat, ytest = ytest, ytrain = ytrain)
    return(resl)
}
load('./output/mcatac/atac_metacell_order_annotation.rda')
load('./output/mcatac/ann_colors.rda')
resl_all <- run_model(motifs = rownames(xdata), 
                        amcs =colnames(mcl_all_norm),
                        xdata = vmnz,
                        ydata = mcl_all_norm[sp_varp_atac_dist,colnames(mcl_all_norm)], 
                        alpha = 1,
                        train_inds = train_inds,
                        test_inds = test_inds)

## Add dinuc data
dinuc <- readr::read_tsv("/net/mraid14/export/tgdata/users/aviezerl/proj/motif_reg/output/mmcortex_clust_dinuc.tsv")
dn <- as.matrix(subset(dinuc, select = -c(chrom, start, end, clust)))
rownames(dn) <- misha.ext::convert_misha_intervals_to_10x_peak_names(dinuc)
xdata_dn_cent <- rbind(xdata, t(dn[,grep('center', colnames(dn))]))
xdata_dn_cent_reg <- rbind(xdata, t(dn[,grep('center|region', colnames(dn))]))
xdata_dn_all <- rbind(xdata, t(dn))
resl_cent <- run_model(motifs = rownames(xdata_dn_cent), 
                        amcs = colnames(mcl_all_norm), 
                        xdata = xdata_dn_cent, 
                        ydata = mcl_all_norm[sp_varp_atac_dist,colnames(mcl_all_norm)], alpha = 1,
                        train_inds = train_inds,
                        test_inds = test_inds)
resl_cent_reg <- run_model(motifs = rownames(xdata_dn_cent_reg), 
                        amcs = colnames(mcl_all_norm), 
                        xdata = xdata_dn_cent_reg, 
                        ydata = mcl_all_norm[sp_varp_atac_dist,colnames(mcl_all_norm)], alpha = 1,
                        train_inds = train_inds,
                        test_inds = test_inds)
resl_dn_all <- run_model(motifs = rownames(xdata_dn_all), 
                        amcs = colnames(mcl_all_norm), 
                        xdata = xdata_dn_all, 
                        ydata = mcl_all_norm[sp_varp_atac_dist,colnames(mcl_all_norm)], alpha = 1,
                        train_inds = train_inds,
                        test_inds = test_inds)
resl_dn_only <- run_model(motifs = colnames(dn), 
                        amcs = colnames(mcl_all_norm), 
                        xdata = t(dn), 
                        ydata = mcl_all_norm[sp_varp_atac_dist,colnames(mcl_all_norm)], alpha = 1,
                        train_inds = train_inds,
                        test_inds = test_inds)

p_motif_and_cent <- pheatmap(resl_cent_reg$coef_mat, breaks = c(seq(min(resl_cent_reg$coef_mat), 0, l = 50), 
            seq(1e-10, max(resl_cent_reg$coef_mat), l=51)), color = colorRampPalette(c('blue3', 'white', 'red3'))(100), 
            annotation_col = row_annot, annotation_colors = ann_colors,
            main = "Coefficient matrix for motifs and dinucs in center and region domains")
save_pheatmap(p_motif_and_cent, './output/sequence_modeling/coef_mat_motifs_and_dinuc_center_region.png', h = 1800, w = 2800, res= 200)

p_dn_only <- pheatmap(resl_dn_only$coef_mat, breaks = c(seq(min(resl_dn_only$coef_mat), 0, l = 50), 
            seq(1e-10, max(resl_dn_only$coef_mat), l=51)), color = colorRampPalette(c('blue3', 'white', 'red3'))(100), 
            annotation_col = row_annot, annotation_colors = ann_colors,
            main = "Coefficient matrix for dinucs (all domains)")
save_pheatmap(p_dn_only, './output/sequence_modeling/coef_mat_dinuc_only.png', h = 1800, w = 2800, res= 200)

p_dn_all <- pheatmap(resl_dn_all$coef_mat, breaks = c(seq(min(resl_dn_all$coef_mat), 0, l = 50), 
                seq(1e-10, max(resl_dn_all$coef_mat), l=51)), color = colorRampPalette(c('blue3', 'white', 'red3'))(100), 
                    annotation_col = row_annot, annotation_colors = ann_colors,
            main = "Coefficient matrix for motifs and dinucs in all domains")
save_pheatmap(p_dn_all, './output/sequence_modeling/coef_mat_motifs_and_dinuc_all.png', h = 1800, w = 2800, res= 200)

png('./output/sequence_modeling/figs/spearman_predic_dinucs_only.png', width =2800, height = 1000, res = 200)
par(cex.lab = 1, mar = c(5,5,4,2))
barplot(diag(tgs_cor(resl_dn_only$predictions_test, as.matrix(ytest), spearman=T))[amc_ord], 
            col = ann_colors[['cell_type']][row_annot$cell_type[amc_ord]], ylab=  'Spearman', 
            main = 'Spearman of predicitions vs accessibility\nonlydinucs (all domains)')
dev.off()
png('./output/sequence_modeling/figs/spearman_predic_all_motifs_and_dinucs_center_and_region.png', width =2800, height = 1000, res = 200)
par(cex.lab = 1, mar = c(5,5,4,2))
barplot(diag(tgs_cor(resl_cent_reg$predictions_test, as.matrix(ytest), spearman=T))[amc_ord], 
            col = ann_colors[['cell_type']][row_annot$cell_type[amc_ord]], ylab=  'Spearman', 
            main = 'Spearman of predicitions vs accessibility\n20 motifs and all dinucs, "center" and "region" domains')
dev.off()
png('./output/sequence_modeling/figs/spearman_predic_all_motifs_and_dinucs.png', width =2800, height = 1000, res = 200)
par(cex.lab = 1, mar = c(5,5,4,2))
barplot(diag(tgs_cor(resl_dn_all$predictions_test, as.matrix(ytest), spearman=T))[amc_ord], 
            col = ann_colors[['cell_type']][row_annot$cell_type[amc_ord]], ylab=  'Spearman', 
            main = 'Spearman of predicitions vs accessibility\n20 motifs and all dinucs and domains')
dev.off()



motifs_to_plot <- setdiff(rownames(resl_all[['coef_mat']]), '(Intercept)')
p_lasso_all <- pheatmap(resl_all[['coef_mat']][motifs_to_plot,amc_ord], cluster_cols = F, 
                    annotation_col = row_annot,
                annotation_colors = ann_colors,
                colorRampPalette(c('blue2', 'white', 'red2'))(100), 
                breaks = c(seq(pmin(0,min(resl_all[['coef_mat']])), 0, l=50), 
                            seq(1e-10,pmax(0, max(resl_all[['coef_mat']])),l=51)))
save_pheatmap(p_lasso_all, 'output/sequence_modeling/figs/lasso_all_amc_coef_mat.png', width = 1800, height = 1300, res = 100)
ytest <- mcl_all_norm[sp_varp_atac_dist(test_inds), ]
cor_pred_list <- lapply(list(resl_cent, resl_cent_reg, resl_dn_all), function(x) tgs_cor(as.matrix(x$predictions_test), as.matrix(ytest), spearman=T))
cor_pred <- tgs_cor(as.matrix(resl_all$predictions_test), as.matrix(ytest), spearman=F)
pheatmap(cor_pred, cluster_rows = F, cluster_cols = F, 
                colorRampPalette(c('blue2', 'white', 'red2'))(100), 
                # breaks = c(seq(pmin(0,min(cor_pred)), 0, l=50), 
                #             seq(1e-10,pmax(0, max(cor_pred)),l=51))
                )
png('./output/sequence_modeling/figs/lasso_spearman_barplot.png', width = 1500, height = 500, res = 100)
barplot(diag(tgs_cor(as.matrix(resl_all$predictions_test), 
                as.matrix(resl_all$ytest), spearman=T))[amc_ord], 
                col = ann_colors$cell_type[match(row_annot[amc_ord,], names(ann_colors$cell_type))], 
                ylab = 'Spearman between predicted and observed accessibility')
purrr::walk(seq(0,0.6,0.1), function(x) lines(c(1,220), c(x,x), lty = 2, lwd = 1.5, col = 'grey'))
dev.off()

cor(diag(tgs_cor(as.matrix(resl_all$predictions_test), 
                as.matrix(resl_all$ytest), spearman=T)), colSds(resl_all$ytest), method = 'spearman')


library(gridExtra)
amd <- all_motif_datasets()
motifs_to_plot <- setdiff(rownames(resl_all[['coef_mat']]), '(Intercept)')
motifs_to_plot_ord <- motifs_to_plot[order(match(motifs_to_plot, 
                                rownames(resl_all[['coef_mat']])[p_lasso_all$tree_row$order]))]
motif_plots <- lapply(motifs_to_plot_ord, function(nm) prego::plot_pssm_logo(filter(amd, motif == nm), title = nm))
names(motif_plots) <- motifs_to_plot_ord
gall <- marrangeGrob(grobs = motif_plots, ncol =1, nrow = length(motif_plots))
ggsave('./output/sequence_modeling/figs/lasso_motif_logos.png', gall, height = 40, width = 8, units = 'in')



seqs_train <- seqs_all[match(sp_varp_atac_dist[!(1:length(sp_varp_atac_dist) %in% test_inds)], rownames(aaa$mcl_all_norm))]
seqs_test <- seqs_all[match(sp_varp_atac_dist[test_inds], rownamesaaa$(mcl_all_norm))]

res_mat <- resl_all$predictions_train - resl_all$ytrain
colnames(res_mat) <- colnames(aaa$mcl_all)


load('./output/mcatac/pk_col_annot_atac.rda')
p_resmat <- pheatmap(t(res_mat[,amc_ord]), cluster_rows = F,
            cluster_cols = F, show_colnames = F, show_rownames = F, 
            annotation_row = row_annot, annotation_colors = ann_colors, 
            annotation_col = pk_col_annot_atac,
            # breaks = c(seq(min(res_mat), 0, l=50), seq(1e-10, max(res_mat), l=51)), 
                    color = rev(colorRampPalette(c('blue2', 'lightgoldenrodyellow', 'red3'))(100))
                    )
save_pheatmap(p_resmat, 'output/sequence_modeling/figs/lasso_residuals_mat.png', width = 1800, height = 1300, res = 100)

# prego_on_res_train <- plyr::llply(1:ncol(res_mat), function(i) prego::regress_pwm(sequences = seqs_train, 
#                                                                 response = res_mat[,i], 
#                                                                 min_kmer_cor=0.01, 
#                                                                 unif_prior=0.1, 
#                                                                 score_metric='r2', 
#                                                                 final_metric='ks', parallel = T), .parallel=T)
# saveRDS(prego_on_res_train, './output/sequence_modeling/prego_train_on_lasso_residual_results.rds')
prego_on_res_train <- readRDS('./output/sequence_modeling/prego_train_on_lasso_residual_results.rds')
prego_train_pssms <- lapply(prego_on_res_train, function(x) x$pssm)
train_pssm_df <- do.call('rbind', lapply(seq_along(prego_train_pssms), 
                                    function(i) relocate(mutate(prego_train_pssms[[i]], 
                                    motif = paste0('train_amc_', i)), motif, .before = pos)))

train_reg_energies <- prego::gextract_pwm(intervals = coords_exp[match(sp_varp_atac_dist, rownames(mcl_all_norm)),], 
            dataset = train_pssm_df, prior = 0.01, parallel = T)
tre_mat <- as.matrix(subset(train_reg_energies, select = -c(chrom, start, end)))
tre_mat <- tre_mat[,order(as.numeric(gsub('train_amc_','',colnames(tre_mat))))]
rownames(tre_mat) <- misha.ext::convert_misha_intervals_to_10x_peak_names(train_reg_energies[,1:3])
saveRDS(tre_mat, './output/sequence_modeling/energies_of_motifs_from_lasso_residual_regression.rds')



resl_prego <- run_model(motifs = rownames(xdata), 
                                                amcs =colnames(mcl_all_norm), 
                                                alpha = 1, 
                                                prego_mat = tre_mat)
p_lasso_prego <- pheatmap(resl_prego[['coef_mat']][,amc_ord], cluster_cols = F, 
                    annotation_col = row_annot,
                annotation_colors = ann_colors,
                colorRampPalette(c('blue2', 'white', 'red2'))(100), 
                breaks = c(seq(pmin(0,min(resl_all[['coef_mat']])), 0, l=50), 
                            seq(1e-10,pmax(0, max(resl_all[['coef_mat']])),l=51)))

motifs_to_plot <- setdiff(rownames(resl_prego[['coef_mat']]), '(Intercept)')
motifs_to_plot_ord <- motifs_to_plot[order(match(motifs_to_plot, 
                                rownames(resl_prego[['coef_mat']])[p_lasso_prego$tree_row$order]))]
motif_plots <- lapply(seq_along(prego_train_pssms), function(i) prego::plot_pssm_logo(prego_train_pssms[[i]],
                     title = colnames(mcl_all)[[i]]))
names(motif_plots) <- colnames(mcl_all)
nsc_amc_inds <- which(row_annot == 'NSC')
gall <- marrangeGrob(grobs = motif_plots[nsc_amc_inds], ncol =7, nrow = 6)
ggsave('./output/sequence_modeling/figs/prego_res_lasso_nsc_motif_logos.png', gall, height = 40, width = 40, units = 'in')





cor_pred_prego <- tgs_cor(as.matrix(resl_prego$predictions_test), as.matrix(ytest), spearman=F)