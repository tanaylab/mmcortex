# first plot clustered correlation matrix



generate_flow_clusters = function() {
  
  excluded_mcs = NULL
  
  mct_id = "sing_emb_wt10"
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  heatmap_gset_id = "sing_emb_wt10_bs500f"
  n_cluster = 65
  fig_dir = "figs/test.flow_cluster"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  mc_cl_fn = "data/test.flow_clust_mc.txt"
  
  mc = scdb_mc(mc_id)
  mc_included = c(1:ncol(mc@e_gc))
  mc_included = setdiff(mc_included,excluded_mcs)
  
  
  marks = c("Eomes","Mesp1","T","Snai1","Utf1")
  
  mc_cl = my_clust_mc_by_flow(mct_id = mct_id,
                              mat_id = mat_id,
                              heatmap_gset_id = heatmap_gset_id,
                              n_cluster = n_cluster,
                              fig_dir = fig_dir,
                              mc_cl_fn = mc_cl_fn,
                              mc_included = mc_included,
                              marks = marks)
  
  
  mcell_mc_plot_marks(mc_id = mc_id,gset_id = heatmap_gset_id,mat_id = mat_id,fig_fn = sprintf("%s/heatmap_marker_genes.png",fig_dir),mc_ord = mc_cl$mc)
  
  return(mc_cl)
}

my_clust_mc_by_flow = function(mct_id, mat_id,heatmap_gset_id,n_cluster,fig_dir,mc_cl_fn,marks = NULL,cls_ord = NULL,mc_clust = NULL,mc_included = NULL) {
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  
  
  fclst = my_mctnetwork_plot_cormat(mct_id = mct_id,mat_id = mat_id,gset_id = heatmap_gset_id, K = n_cluster,
                                  fn =  sprintf("%s/flow_clust_mc.png",fig_dir),h=3000,w=3000, text_cex=0.6,
                                  marks = marks,mc_clust = mc_clust,mc_included = mc_included)
  
  
  mc_to_cl = data.frame(mc = as.numeric(names(fclst)), cls = fclst)
  if(0) {
    if(is.null(mc_clust)) {
      write.table(x = mc_to_cl,file = mc_cl_fn,sep = "\t")
    }
  }

  return(mc_to_cl)
}


my_mctnetwork_plot_cormat = function(mct_id, mat_id, gset_id, K, fn, w=1600, h=1600, text_cex=1, cls_pre_order=NULL, marks=NULL,mc_clust = NULL,mc_included = NULL)
{
  mct = scdb_mctnetwork(mct_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
  }
  clst_flows = my_mctnetwork_clust_flows(mct_id, K,mc_included = mc_included)
  
  
  cmat = clst_flows$cmat
  fclst = clst_flows$clust
  if(!is.null(mc_clust)) {
    mc_clust = mc_clust[order(as.numeric(names(mc_clust)))]
    fclst = mc_clust
  }
  
  mc_ord = clst_flows$hc$order
  if(!is.null(cls_pre_order)) {
    cls_rank = 1:K
    names(cls_rank) = as.character(1:K)
    cls_rank[as.character(cls_pre_order)] = 1:K
    mc_mean_age = apply(mct@mc_t, 1, function(x) {return(mean(x*(1:length(x))/sum(x))) })
    mc_ord = order(-cls_rank[fclst]*1000+mc_mean_age)
  }
  
  mc = scdb_mc(mct@mc_id)
  
  shades = colorRampPalette(c("darkblue", "blue","white", "red", "yellow"))(1000)
  png(fn, w=w, h=h)
  n_mc = nrow(cmat)
  
  layout(matrix(c(1,2),nrow=2),heights=c(nrow(cmat)*9+50, 300))
  par(mar=c(0,5,4,5))
  image(cmat[mc_ord, mc_ord], zlim=c(-1,1), col=shades, yaxt='n', xaxt='n')
  N = length(mc_ord)
  
  mc_x = 1:length(mc_ord)
  names(mc_x) = names(fclst)
  mc_x[mc_ord] = 1:N
  for(i in 1:K) {
    abline(h=max(-0.5+mc_x[fclst == i])/(N-1))
    abline(v=max(-0.5+mc_x[fclst == i])/(N-1))
  }
  cl_x = tapply(mc_x, fclst, mean)
  cl_max = tapply(mc_x, fclst, max)
  
  mtext(1:K, side = 3, at=cl_x/N, las=1, cex=1.5)
  
  mtext((1:length(mc_ord))[mc_ord], side = 2, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
  mtext((1:length(mc_ord))[mc_ord], side = 4, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
  par(mar=c(3,5,0,5))
  image(as.matrix(as.numeric(names(fclst))[mc_ord],nrow=1), col=mc@colors[min(as.numeric(names(fclst))):max(as.numeric(names(fclst)))], yaxt='n', xaxt='n')
  mtext(names(fclst)[mc_ord], side = 1, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
  dev.off()
  
  if(!is.null(marks)) {
    egc = log2(mc@e_gc+1e-5)
    dir = sub(".png","_marks",fn)
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    for(g in marks) {
      print(g)
      png(sprintf("%s/%s.png", dir, g),w=1000,h=300)
      plot(1:length(mc_ord), egc[g, names(fclst)[mc_ord]], pch=19, col=mc@colors[as.numeric(names(fclst)[mc_ord])], ylab=g,xaxt='n')
      mtext(1:K, at=cl_x,side=1, las=2)
      
      abline(v=cl_max+0.5)
      grid()
      dev.off()
    }
  }
  mark_fn = sub("png", "mat.png", fn)
  #mcell_mc_plot_marks(mct@mc_id, gset_id=gset_id, mat_id = mat_id, fig_fn=mark_fn, mc_ord = mc_ord, plot_cells=T)
  return(fclst[mc_ord])
}	


my_mctnetwork_clust_flows = function(mct_id, K,mc_included = NULL) {
  mct = scdb_mctnetwork(mct_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
  }
  mc = scdb_mc(mct@mc_id)
  

  
  fls = mctnetwork_get_flow_mat(mct,-1)
  if(!is.null(mc_included)) {
   fls = fls[mc_included,mc_included]
   rownames(fls) = mc_included
   colnames(fls) = mc_included
  }
  #fls3 = fls %*% fls %*% fls
  fls3 = fls %*% fls
  cr3 = tgs_cor(rbind(fls3, t(fls3)))
  diag(cr3) = 0
  cr3_2 = tgs_cor(cr3)
  diag(cr3_2) = 0
  
  hc = hclust(tgs_dist(t(cr3_2)), "ward.D")
  
  clf = cutree(hc, K)
  return(list(cmat = cr3_2, clust = clf, hc=hc))
}


