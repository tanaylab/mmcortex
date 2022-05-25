# first plot clustered correlation matrix
mct_id = "sing_emb_wt10_cap040"
mat_id = "sing_emb_wt10"
heatmap_gset_id = "sing_emb_wt10_bs500f"
n_cluster = 65
fig_dir = "figs/wt10.network"
mc_cl_fn = "data/flow_clust_mc.txt"



plot_reannotated_mc = function() {
  mct_id = "sing_emb_wt10_cap040"
  mat_id = "sing_emb_wt10"
  heatmap_gset_id = "sing_emb_wt10_bs500f"
  n_cluster = 65
  fig_dir = "figs/wt10.network"
  mc_cl_fn = "data/flow_clust_mc.txt"
  marks = read.table("config/flow_clust_marks.txt",sep = "\t",stringsAsFactors = F)
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  cls_ord = read.table("data/wt10.cluster_annotation/cluster_order.txt",sep = "\t",stringsAsFactors = F)
  mc_cluster_annot = read.table("data/wt10.cluster_annotation/mc_cluster_order.txt",sep = "\t",stringsAsFactors = F)
  cls_ord = cls_ord$x
  mc_clust = mc_cluster_annot$new_cluster_id
  names(mc_clust) = mc_cluster_annot$mc_id
  
  mc_ord = my_mctnetwork_plot_cormat(mct_id = mct_id,mat_id = mat_id,gset_id = heatmap_gset_id, K = n_cluster,
                                     fn =  sprintf("%s/flow_clust_mc.png",fig_dir),h=3000,w=3000, text_cex=0.6, 
                                     cls_pre_order=rev(cls_ord), marks = marks,mc_clust = mc_clust)
  
  
  mctnetwork_plot_net(mct_id = mct_id,fn = sprintf("%s/net.png",fig_dir),mc_ord = mc_ord,w = 2000,h=4000,dx_back = 0,dy_ext = 0,plot_mc_ids = F)
  
  
  mctnetwork_plot_net(mct_id = mct_id,fn = sprintf("%s/net_large.png",fig_dir),mc_ord = mc_ord,w = 3000,h=6000,dx_back = 0,dy_ext = 0,plot_mc_ids = T)
}


clust_mc_by_flow = function(mct_id, mat_id,heatmap_gset_id,n_cluster,fig_dir,mc_cl_fn,marks = NULL,cls_ord = NULL,mc_clust = NULL) {
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  #mc_clust = mctnetwork_clust_flows(mct_id = mct_id,K = n_cluster)
  
  mc_ord = my_mctnetwork_plot_cormat(mct_id = mct_id,mat_id = mat_id,gset_id = heatmap_gset_id, K = n_cluster,
                                  fn =  sprintf("%s/flow_clust_mc.png",fig_dir),h=3000,w=3000, text_cex=0.6, 
                                  cls_pre_order=rev(cls_ord), marks = marks,mc_clust = mc_clust)
  
  
  mc_to_cl = data.frame(mc = c(1:length(mc_clust$clust)), cls = mc_clust$clust)
  write.table(x = mc_to_cl,file = mc_cl_fn,sep = "\t")
  return(mc_to_cl)
}


my_mctnetwork_plot_cormat = function(mct_id, mat_id, gset_id, K, fn, w=1600, h=1600, text_cex=1, cls_pre_order=NULL, marks=NULL,mc_clust = NULL)
{
  mct = scdb_mctnetwork(mct_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
  }
  clst_flows = mctnetwork_clust_flows(mct_id, K)
  
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
  names(mc_x) = 1:N
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
  image(as.matrix(mc_ord,nrow=1), col=mc@colors, yaxt='n', xaxt='n')
  mtext((1:length(mc_ord))[mc_ord], side = 1, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
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
      plot(1:length(mc_ord), egc[g, mc_ord], pch=19, col=mc@colors[mc_ord], ylab=g,xaxt='n')
      mtext(1:K, at=cl_x,side=1, las=2)
      
      abline(v=cl_max+0.5)
      grid()
      dev.off()
    }
  }
  mark_fn = sub("png", "mat.png", fn)
  mcell_mc_plot_marks(mct@mc_id, gset_id=gset_id, mat_id = mat_id, fig_fn=mark_fn, mc_ord = mc_ord, plot_cells=T)
  return(mc_ord)
}	