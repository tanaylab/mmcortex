amos_peak_calling_function <- function(trk, percentiles = 0.99, canonical = FALSE) {
    q <- gquantiles(trk, percentiles = percentiles, iterator = 20)
    atac_hits = gscreen(glue::glue("{trk}>{q}"))
	hits_v = gextract(trk, intervals=atac_hits, iterator=20)
	hits_v = hits_v[order(-hits_v[,trk]),]
	hits_peak = hits_v[!duplicated(hits_v$intervalID),]
	hits_peak$start = hits_peak$start - 140
	hits_peak$end = hits_peak$end + 140
	hits_peak = gintervals.force_range(hits_peak)
	hits_peak = gintervals.canonic(hits_peak)
	hits_v = gextract(trk, intervals=hits_peak, iterator=20)
	hits_v = hits_v[order(-hits_v[,trk]),]
	hits_peak = hits_v[!duplicated(hits_v$intervalID),]
	hits_peak$start = hits_peak$start - 140
	hits_peak$end = hits_peak$end + 140
	hits_peak = hits_peak[order(hits_peak$intervalID),]
    if (canonical) {
        hits_peak_can <- gintervals.canonical(hits_peak)
        return(hits_peak_can)
    }
    return(hits_peak)
}

proximal_chromatin_activity <- function(peaks_of_interest,
                                        background_peaks,
                                        mc_sel,
                                        # day,
                                        mca,
                                        mat_rna,
                                        mc_rna,
                                        tss,
                                        tads,
                                        d_puncture = 1e+3,
                                        peak_clusters = NULL,
                                        pred = NULL,
                                        d_proximity_atac = 5e+4, 
                                        d_proximity_rna = 5e+5, 
                                        eps = 1e-5,
                                        restrict_to_tads = TRUE) {
    if (!(tibble::has_name(peaks_of_interest, 'peak_name') & tibble::has_name(background_peaks, 'peak_name'))) {
        stop('Peaks of interest or background peaks do not have peak name')
    }
    # if (length(which(!(peaks_of_interest$peak_name %in% background_peaks$peak_name))) > 0) {
    #     stop('Not all peaks of interest are in background peak set')
    # }
    peaks_of_interest <- peaks_of_interest[!duplicated(peaks_of_interest$peak_name),]
    background_peaks <- background_peaks[!duplicated(background_peaks$peak_name),]
    nei_peaks_peaks <- gintervals.neighbors(peaks_of_interest, background_peaks, maxdist = d_proximity_atac, mindist = -d_proximity_atac, maxneighbors = 1e+6)
    colnames(nei_peaks_peaks)[grep('peak_name', colnames(nei_peaks_peaks))] <- c('peak_name_1', 'peak_name_2')
    nei_peaks_peaks <- nei_peaks_peaks[nei_peaks_peaks$peak_name_1 != nei_peaks_peaks$peak_name_2 & abs(nei_peaks_peaks$dist) >= d_puncture,]

    nei_bg_peaks_tads <- gintervals.neighbors(tads, background_peaks, maxdist = 0, mindist = 0, maxneighbors = 1e+4)
    nei_fg_peaks_tads <- gintervals.neighbors(tads, peaks_of_interest, maxdist = 0, mindist = 0, maxneighbors = 1e+4)
    marginal_peaks_umis <- Matrix::rowSums(mca@mat[background_peaks$peak_name,mc_sel])
    nei_peaks_peaks$tads1 <- nei_fg_peaks_tads$tad_name[match(nei_peaks_peaks[,4], nei_fg_peaks_tads$peak_name)]
    nei_peaks_peaks$tads2 <- nei_bg_peaks_tads$tad_name[match(nei_peaks_peaks[,8], nei_bg_peaks_tads$peak_name)]
    if (restrict_to_tads) {
        nei_peaks_peaks <- nei_peaks_peaks[which(nei_peaks_peaks$tads1 == nei_peaks_peaks$tads2),]
    }
    nei_peaks_peaks$atac_umis_2 <- marginal_peaks_umis[match(nei_peaks_peaks$peak_name_2, names(marginal_peaks_umis))]
    tad_atac_umis <- tapply(nei_peaks_peaks$atac_umis_2, nei_peaks_peaks$peak_name_1, sum)
    
    npgf_sum_atac <- tad_atac_umis
    
    npgf_sum_atac[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% names(npgf_sum_atac))]] <- 0
    npgf_sum_atac[is.na(npgf_sum_atac)] <- 0
    
    npgf_sum_atac <- npgf_sum_atac[peaks_of_interest$peak_name]
    

    
    nei_peaks_genes <- gintervals.neighbors(peaks_of_interest, tss, maxdist = d_proximity_rna, mindist = -d_proximity_rna, maxneighbors = 1e+4)
    
    # print(head(nei_peaks_peaks))
    peaks_tads <- gintervals.neighbors(nei_peaks_genes[,1:4], tads, maxdist = 0, mindist = 0, maxneighbors = 1)
    genes_tads <- gintervals.neighbors(dplyr::rename(nei_peaks_genes[,c(5:7, grep('geneSymbol', colnames(nei_peaks_genes)))], 
                                            chrom = chrom1, start = start1, end = end1),
                                       tads, maxdist = 0, mindist = 0, maxneighbors = 1)
    nei_peaks_genes$peak_tad <- peaks_tads$tad_name
    nei_peaks_genes$gene_tad <- genes_tads$tad_name
    if (restrict_to_tads) {
        npgf <- dplyr::select(nei_peaks_genes[which(nei_peaks_genes$peak_tad == nei_peaks_genes$gene_tad),], 
                      c(1:4, grep('geneSymbol|peak_tad|gene_tad', colnames(nei_peaks_genes))))
    } else {
        npgf <- dplyr::select(nei_peaks_genes, 
                      c(1:4, grep('geneSymbol|peak_tad|gene_tad', colnames(nei_peaks_genes))))
    }
    
    genes_nei_umis <- Matrix::rowSums(mat_rna@mat[unique(npgf$geneSymbol[npgf$geneSymbol %in% rownames(mat_rna@mat)]),
                                                      names(mc_rna@mc[mc_rna@mc %in% mc_sel])])
    npgf$rna_umis <- genes_nei_umis[npgf$geneSymbol]
    npgf$rna_umis[is.na(npgf$rna_umis)] <- 0
    npgf_sum_rna <- tapply(npgf$rna_umis, npgf$peak_name, function(x) sum(x))
    npgf_sum_rna[peaks_of_interest$peak_name[!(peaks_of_interest$peak_name %in% npgf$peak_name)]] <- 0
    npgf_sum_rna[is.na(npgf_sum_rna)] <- 0
    npgf_sum_rna <- npgf_sum_rna[peaks_of_interest$peak_name]
    
    
    sum_umi_df <- dplyr::left_join(tibble::enframe(npgf_sum_atac, name = 'peak_name', value = 'atac_umi'), 
                                   tibble::enframe(npgf_sum_rna, name = 'peak_name', value = 'rna_umi'), by = 'peak_name')
    # r_vec <- egc_by_day_n[,day]
    r_vec <- Matrix::rowSums(mca@mat[,mc_sel])
    r_vec <- log2(1e-5 + r_vec/sum(r_vec))
    sum_umi_df$r <- r_vec[sum_umi_df$peak_name]
    if(!is.null(peak_clusters)) {
        sum_umi_df$cluster <- as.numeric(names(peak_clusters)[match(sum_umi_df$peak_name, peak_clusters)])
    }
    if (!is.null(pred)) {sum_umi_df$pred <- pred[sum_umi_df$peak_name]}
    return(sum_umi_df)
}

transform_x_to_scale <- function(x, new_range = c(0,1)) {
    mnmx <- quantile(x, c(0,1), na.rm = T)
    xlin <- (x - min(x, na.rm = T))/(diff(mnmx))
    xnew <- min(new_range)+diff(new_range)*xlin
    return(xnew)
}

boxplot_vec <- function(xvec, yvec, nm, num_bins = 7, bins = NULL, ylab = '', xlab = '', xlim = NULL, ylim = NULL, xaxt = 's', yaxt = 's', col = 'gray', show_text = TRUE, text_cex = 1.2, text_y_factor = 1, text_y_quantile = .99) {
    xvec_rng <- c(min(xvec, na.rm = T),max(xvec, na.rm = T))
    if (is.null(bins)) {
        cvc <- seq(xvec_rng[[1]]-diff(xvec_rng)/1e+3, xvec_rng[[2]], l = num_bins)
    } else {
        cvc <- bins
    }
    mtfc <- cut(xvec, breaks = cvc)
    boxplot(yvec ~ mtfc, main = nm, ylab = ylab, xlab  = xlab, ylim = ylim, xlim = xlim, xaxt = xaxt, yaxt = yaxt, col = col)
    if (show_text) {
        text(1:length(levels(mtfc)), text_y_factor*rep(quantile(yvec, text_y_quantile, na.rm = T), length(levels(mtfc))), labels = table(mtfc), col = 'red', cex = text_cex)
    }
    
}

vioplot_vec <- function(xvec, yvec, nm, ylab = '', xlab = '') {
    xvec_rng <- c(min(xvec, na.rm = T),max(xvec, na.rm = T))
    cvc <- seq(xvec_rng[[1]], xvec_rng[[2]], l = 7)
    mtfc <- cut(xvec, breaks = cvc)
    vioplot::vioplot(yvec ~ mtfc, main = nm, ylab = ylab, xlab  = xlab)
    text(1:length(levels(mtfc)), rep(quantile(yvec, 0.75), length(levels(mtfc))), labels = table(mtfc), col = 'red', cex = 1.2)
}

multintersect <- function(...) {
    args <- list(...)
    if (length(args) < 2) {stop('Must supply at least two arguments')}
    temp <- args[[1]]
    for (j in 2:length(args)) {
        temp <- intersect(temp, args[[j]])
    }
    print(glue::glue('Length of arguments was {paste0(unlist(lapply(args, length)), collapse = ", ")}'))
    print(glue::glue('Length of intersection is {length(temp)}'))
    return(temp)
}

multunion <- function(...) {
    args <- list(...)
    if (length(args) < 2) {stop('Must supply at least two arguments')}
    temp <- args[[1]]
    for (j in 2:length(args)) {
        temp <- union(temp, args[[j]])
    }
    print(glue::glue('Length of arguments was {paste0(unlist(lapply(args, length)), collapse = ", ")}'))
    print(glue::glue('Length of union is {length(temp)}'))
    return(temp)
}

jacc <- function(x,y) {
    return(length(intersect(x,y))/length(union(x,y)))
}

save_pheatmap_pdf <- function(x, filename, width=2500, height=2500) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree,fontsize = NULL) {
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree, y_factor){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    print(head(full.range))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    print(head(selected.range))
    target_d <- seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = round(y_factor*sum(d.select)))
    # new_d <- selected.range + y_factor*(target_d - selected.range)
    new_d <- target_d
    return(unit(new_d, "npc"))
  }
  if (!is.null(fontsize)) {
      # new.label$gp = grid::gpar(fontsize = fontsize)
      y_factor <- fontsize/16
  } else {
      y_factor <- 1
  }
  print(y_factor)
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "",
                               y_factor = y_factor)
  new.flag <- grid::segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions


  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
  )
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid::grid.newpage()
  grid::grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

get_mc_cc = function(mat_id, mc_id, mc2d_id, nm, plot_mc2d = FALSE, mc2d_png_path = NULL) {
  ## function from Markus
  m_0 = 0.0025
  s_0 = 0.001
  m_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr",
                "Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgo2a",
                "Sgo2b","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)

  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  tot  = Matrix::colSums(m@mat)
  s_tot = Matrix::colSums(m@mat[s_genes,])
  m_tot = Matrix::colSums(m@mat[m_genes,])
  s_score = s_tot/tot
  m_score = m_tot/tot

  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  mc_cc = tibble::rownames_to_column(data.frame(mc_cc))

  colnames(mc_cc) = c('mc', 'cc_score')

  mc_cc = dplyr::mutate(mc_cc, mc = as.numeric(mc))
  mc2d = scdb_mc2d(mc2d_id)

  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  if (plot_mc2d == TRUE) {
    if (is.null(mc2d_png_path)) {
        mc2d_png_path <- paste0('./', mc2d_id, '.png')
        print(glue::glue('mc2d_png_path not given, saving to {mc2d_png_path}'))
    } 
    png(mc2d_png_path, w=800, h=800)
    plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.4, col=ifelse(f[names(mc2d@sc_x)], "lightgray","black"))
    points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=2.5, bg=shades[101 - mc_cc$cc_score])
    dev.off()
  }
  return(mc_cc)
}

plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=seq(1,101,l=6))
{
  if (!is.null(fig_fn)) {
      png(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  text(19, show_vals_ind,cex = 2, labels=signif(vals[show_vals_ind], 1), pos=4)
  if (!is.null(fig_fn)) {
    dev.off()
  }
}

rowMaxmins <- function(mat) {setNames(matrixStats::rowMaxs(mat) - matrixStats::rowMins(mat), rownames(mat))}

get_genes_specific_to_mcs <- function(legc, mc_pos = NULL, mc_neg = NULL, cl_vec = NULL) {
    if (!is.null(mc_pos) && is.null(mc_neg)) {
        cl_vec <- ifelse(colnames(legc) %in% mc_pos, 1, 0)
    } else if (!is.null(mc_pos) && !is.null(mc_neg)) {
        if (!(length(intersect(mc_pos, mc_neg)) == 0)) {
            stop('mc_pos and mc_neg intersect')
        }
        legc <- legc[,c(mc_pos, mc_neg)]
        cl_vec <- c(rep(1, length(mc_pos)), rep(0, length(mc_neg)))
    }
    legc_avg <- t(tgs_matrix_tapply(legc, cl_vec, mean))
    if (ncol(legc_avg) == 2) {
        diffs <- matrixStats::rowDiffs(legc_avg)
        rownames(diffs) <- rownames(legc_avg)
        return(diffs[order( diffs[,1], decreasing = T),])
    } else {
        diffs <- t(plyr::laply(1:ncol(legc_avg), function(i) legc_avg[,i] - matrixStats::rowMaxs(legc_avg[,-i]), .parallel = T))
        colnames(diffs) <- colnames(legc_avg)
        return(lapply(1:ncol(diffs), function(i) {
            df <- diffs[which(diffs[,i] > 0.1),]
            return(df[order(df[,i], decreasing = T),])
        }))
    }
}


mctnetwork_plot_net_YSh = function(mct_id,flow_id, fn, 
						mc_ord = NULL, colors_ordered = NULL,
                        plot_pdf = FALSE,
						propogate=NULL,
						mc_t_score = NULL,
						edge_w_scale=5e-4, 
						flow_thresh = 1e-4,
						w = 2000,h = 2000,
						mc_cex = 0.5,
						dx_back = 0.15, dy_ext = 0.4,
						sigmoid_edge = F, grad_col_edge = F,
						plot_mc_ids = F,miss_color_thresh = 0.5,
            show_over_under_flow = F,
						func_deform=NULL,
                        bg = "white",
						score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000))
{
	if(!is.null(propogate) | !is.null(mc_t_score)) {
		dx_back = 0
		dy_ext = 0
	}
	
  mct = scdb_mctnetwork(mct_id)
  mcf = scdb_mctnetflow(flow_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  net = mct@network
  net$flow = mcf@edge_flows
	mc = scdb_mc(mct@mc_id)
	if(is.null(mc)) {
		stop("cannot find mc object ", mct@mc_id, " matching the mc id in the mctnetwork object! db mismatch? recreate objects?")
	}
	if(is.null(mcf@edge_flows) | sum(mcf@edge_flows)==0) {
		stop("flows seems not to be initialized in mct id ", mct_id, " maybe rerun the mincost algorithm?")
	}
	names(mc@colors) = as.character(1:length(mc@colors))
	
#color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
#order MCs by type, mean age
  if(is.null(mc_ord)) {
	 if(is.null(colors_ordered)) {
			stop("specify either mc_ord or color ord when plotting mctnet network")
	 }
	 mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id,flow_id, colors_ordered)
  } else {
    mc_rank = rep(-1,length(mc_ord))
    mc_rank[mc_ord] = c(1:length(mc_ord))
    names(mc_rank) = as.character(1:length(mc_rank))
  }
  
  mc_rank["-2"] = 0
  mc_rank["-1"] = length(mc_rank)/2
  
  #add growth mc

	f= net$flow > flow_thresh
	nn = net[f,]
	x1 = nn$time1
	x2 = nn$time2
	y1 = as.numeric(mc_rank[as.character(nn$mc1)])
	y2 = as.numeric(mc_rank[as.character(nn$mc2)])

	x1 = ifelse(nn$type1 == "growth", x1 + 0.3, x1)
	x2 = ifelse(nn$type2 == "growth", x2 + 0.3, x2)
	x1 = ifelse(nn$type1 == "norm_b" | nn$type1 == "extend_b",x1-dx_back,x1)
	x2 = ifelse(nn$type2 == "norm_b" | nn$type2 == "extend_b",x2-dx_back,x2)
	y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
	y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
	y2 = ifelse(nn$type2 == "sink", NA, y2)
	y1 = ifelse(nn$type1 == "growth", y2+2.5, y1)
	y2 = ifelse(nn$type2 == "growth", y1+2.5, y2)
	y1 = ifelse(nn$type1 == "extend_b" | nn$type1 == "extend_f",y1+dy_ext, y1)
	y2 = ifelse(nn$type2 == "extend_f" | nn$type2 == "extend_b",y2+dy_ext, y2)

	if(!is.null(func_deform)) {
		x1 = ifelse(nn$type1 == "growth", NA, x1)
		x2 = ifelse(nn$type2 == "growth", NA, x2)
		y1 = ifelse(nn$type1 == "growth", NA, y1)
		y2 = ifelse(nn$type2 == "growth", NA, y2)
		min_x = min(c(x1,x2),na.rm=T)
		min_y = min(c(y1,y2),na.rm=T)
		range_x = (max(c(x1,x2),na.rm=T)-min_x)
		range_y = (max(c(y1,y2),na.rm=T)-min_y)
		ax = c((x1-min_x)/range_x, (x2-min_x)/range_x)
		ay = c((y1-min_y)/range_y, (y2-min_y)/range_y)
		atag = func_deform(ax, ay)
		n1 = length(x1)
		n = length(ax)
		x1 = atag[[1]][1:n1]
		x2 = atag[[1]][(n1+1):n]
		y1 = atag[[2]][1:n1]
		y2 = atag[[2]][(n1+1):n]
		
if(0) {
		xy1 = func_deform((x1-min_x)/range_x, (y1-min_y)/range_y)
		xy2 = func_deform((x2-min_x)/range_x, (y2-min_y)/range_y)
		x1 = xy1[[1]]
		x2 = xy2[[1]]
		y1 = xy1[[2]]
		y2 = xy2[[2]]
}
	}

	nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)	
	f_overflow = nn$type2=="extend_f" & nn$cost > 100
	f_underflow = nn$type2 == "norm_f" & nn$cost < -100
#  nn$flow/(1e-8 + nn$capacity) < miss_color_thresh

    if(plot_pdf) {
        pdf(fn,width = w/1e+2,height =h/1e+2,useDingbats = F)
    } else {
        png(fn, width = w,height = h,bg = bg)
    }
  
	if(is.null(propogate) & is.null(mc_t_score)) {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		mc_rgb = col2rgb(mc@colors)/256
		f = nn$mc1>0 & nn$mc2 > 0
		m1 = as.numeric(nn$mc1[f])
		m2 = as.numeric(nn$mc2[f])
		seg_df = data.frame(x1 = x1[f], y1=y1[f], dx=x2[f]-x1[f], dy=y2[f]-y1[f], 
									r1 = mc_rgb["red",m1],
									r2 = mc_rgb["red",m2],
									g1 = mc_rgb["green",m1],
									g2 = mc_rgb["green",m2],
									b1 = mc_rgb["blue",m1],
									b2 = mc_rgb["blue",m2])

		for(alpha in seq(0,0.98,0.02)) {
			beta = alpha
			beta5 = alpha+0.02
			if(sigmoid_edge) {
			beta = plogis(alpha,loc=0.5,scale=0.1)
			beta5 = plogis(alpha+0.02,loc=0.5,scale=0.1)
			}
			sx1 = seg_df$x1+alpha*seg_df$dx
			sx2 = seg_df$x1+(alpha+0.02)*seg_df$dx
			sy1 = seg_df$y1+beta*seg_df$dy
			sy2 = seg_df$y1+beta5*seg_df$dy
			alpha_col = ifelse(grad_col_edge, alpha,0)
			rgb_r = seg_df$r2*alpha_col+seg_df$r1*(1-alpha_col)
			rgb_g = seg_df$g2*alpha_col+seg_df$g1*(1-alpha_col)
			rgb_b = seg_df$b2*alpha_col+seg_df$b1*(1-alpha_col)
			cols = rgb(rgb_r, rgb_g, rgb_b)
			segments(sx1, sy1, sx2, sy2, 
				col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
				lwd=pmin(nn$flow/edge_w_scale, 10))
		}
#		segments(x1,y1,x2,y2, 
#				col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
#				lwd=pmin(nn$flow/edge_w_scale, 10))
		# f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
		# 							lwd=pmin(nn$flow[f]/edge_w_scale, 10))
		# f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
	  #        					lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,10))
#		points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=1)
	} else if(!is.null(mc_t_score)) {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		mc_t_score = pmax(mc_t_score, quantile(mc_t_score,0.03))
		mc_t_score = pmin(mc_t_score, quantile(mc_t_score,0.97))
		mc_t_score = mc_t_score-min(mc_t_score)
		mc_t_score = mc_t_score/max(mc_t_score)
		mc_t_score = floor(1+999*mc_t_score)
		f = nn$mc1>0 & nn$mc2>0 & nn$time1>0
		max_mc = nrow(mc_t_score)
		m1 = as.numeric(nn$mc1[f]) 
		m2 = as.numeric(nn$mc2[f]) 
		score1 = rep(1, nrow(nn))
		score2 = rep(1, nrow(nn))
		score1[f] = mc_t_score[m1+(nn[f,"time1"]-1)*max_mc]
		score2[f] = mc_t_score[m2+(nn[f,"time2"]-1)*max_mc]

		seg_df = data.frame(x1 = x1, y1=y1, dx=x2-x1, dy=y2-y1, 
														score1= score1, dscore=score2-score1)

		for(alpha in seq(0,0.95,0.05)) {
			x1 = seg_df$x1+alpha*seg_df$dx
			x2 = seg_df$x1+(alpha+0.05)*seg_df$dx
			y1 = seg_df$y1+alpha*seg_df$dy
			y2 = seg_df$y1+(alpha+0.05)*seg_df$dy
			cols = score_shades[floor(seg_df$score1 + (alpha+0.025)*seg_df$dscore)]
			segments(x1, y1, x2, y2, 
				col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
				lwd=pmin(nn$flow/edge_w_scale, 10))
		}
		rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
					col = mc@colors[nn$mc1], border=NA)
		rect(seg_df$x2-0.12,seg_df$y2-0.5, seg_df$x2+0.12, seg_df$y2+0.5,
					col = mc@colors[nn$mc2], border=NA)
	} else {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		max_time = length(propogate)
		m1 = as.numeric(nn$mc1) 
		m2 = as.numeric(nn$mc2) 
		max_m = ncol(propogate[[1]])
		prop_flow = rep(0, nrow(nn))
		for(t in 1:max_time) {
			f = (nn$time1 == t) & nn$mc1>0 & nn$mc2>0
			prop_flow[f] = propogate[[t]][m1[f]+max_m*(m2[f]-1)]
		}
		segments(x1,y1,x2,y2, 
				col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
				lwd=pmin(prop_flow/edge_w_scale, 10))
		points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=m_cex)
	}

	if(plot_mc_ids) {
		f1 = nn$type1!="growth" 
		text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
#	  text(c(x1[f1],x2[f2]),c(y1[f1],y2[f2]),labels = c(nn$mc1[f1],nn$mc2[f2]), cex=1)
	}
	
	dev.off()
}

mctnetwork_plot_propagation = function(mct_id, time, p_anchor, fn, 
						mc_ord = NULL, colors_ordered = NULL,
						edge_w_scale=5e-4, 
						w = 2000,h = 2000, mc_cex=1)
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	probs = mctnetwork_propogate_from_t(mct, time, p_anchor)
	mctnetwork_plot_net(mct_id, fn=fn, 
					propogate = probs$step_m,
					mc_ord = mc_ord, colors_ordered = colors_ordered,
					edge_w_scale = edge_w_scale, 
					w = w, h = h, mc_cex = mc_cex)
}

mm_mctnetwork_plot_net = function(mct_id, flow_id, fn, 
                                  mc_ord = NULL, colors_ordered = NULL,
                                  propogate=NULL,
                                  mc_t_score = NULL,
                                  edge_w_scale=5e-4, 
                                  w = 2000,h = 2000,
                                  mc_cex = 0.5,
                                  dx_back = 0.15, dy_ext = 0.4,
                                  sigmoid_edge = F, grad_col_edge = F,
                                  plot_mc_ids = F,miss_color_thresh = 0.5,
                                  func_deform=NULL,
                                  plot_background_as_grey = F,
                                  bg_col = "gray90",
                                  score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000),
                                  bg_scale = 1,
                                  fr_scale = 2,
                                  max_lwd = 10,
                                  plot_pdf = FALSE,
                                  show_over_under_flow = F,
                                  show_axes = T,
                                  bg = "white")
{
  if(!is.null(propogate) | !is.null(mc_t_score)) {
    dx_back = 0
    dy_ext = 0
  }
  mct = scdb_mctnetwork(mct_id)
  mcf = scdb_mctnetflow(flow_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  net = mct@network
  net$flow = mcf@edge_flows
  mc = scdb_mc(mct@mc_id)
  if(is.null(mc)) {
    stop("cannot find mc object ", mct@mc_id, " matching the mc id in the mctnetwork object! db mismatch? recreate objects?")
  }
  if(is.null(mcf@edge_flows) | sum(mcf@edge_flows)==0) {
    stop("flows seems not to be initialized in mct id ", mct_id, " maybe rerun the mincost algorithm?")
  }
  names(mc@colors) = as.character(1:length(mc@colors))
  
  #color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
  #order MCs by type, mean age
  if(is.null(mc_ord)) {
    if(is.null(colors_ordered)) {
      stop("specify either mc_ord or color ord when plotting mctnet network")
    }
    mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id, flow_id, colors_ordered)
  } else {
    mc_rank = rep(-1,length(mc_ord))
    mc_rank[mc_ord] = c(1:length(mc_ord))
    names(mc_rank) = as.character(1:length(mc_rank))
  }
  
  mc_rank["-2"] = 0
  mc_rank["-1"] = length(mc_rank)/2
  
  #add growth mc
  
  f= net$flow > 1e-4
  nn = net[f,]
  x1 = nn$time1
  x2 = nn$time2
  y1 = as.numeric(mc_rank[as.character(nn$mc1)])
  y2 = as.numeric(mc_rank[as.character(nn$mc2)])
  
  x1 = ifelse(nn$type1 == "growth", x1 + 0.3, x1)
  x2 = ifelse(nn$type2 == "growth", x2 + 0.3, x2)
  x1 = ifelse(nn$type1 == "norm_b" | nn$type1 == "extend_b",x1-dx_back,x1)
  x2 = ifelse(nn$type2 == "norm_b" | nn$type2 == "extend_b",x2-dx_back,x2)
  y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
  y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
  y2 = ifelse(nn$type2 == "sink", NA, y2)
  y1 = ifelse(nn$type1 == "growth", y2+2.5, y1)
  y2 = ifelse(nn$type2 == "growth", y1+2.5, y2)
  y1 = ifelse(nn$type1 == "extend_b" | nn$type1 == "extend_f",y1+dy_ext, y1)
  y2 = ifelse(nn$type2 == "extend_f" | nn$type2 == "extend_b",y2+dy_ext, y2)
  
  if(!is.null(func_deform)) {
    min_x = min(c(x1,x2),na.rm=T)
    min_y = min(c(y1,y2),na.rm=T)
    range_x = (max(c(x1,x2),na.rm=T)-min_x)
    range_y = (max(c(y1,y2),na.rm=T)-min_y)
    xy1 = func_deform((x1-min_x)/range_x, (y1-min_y)/range_y)
    xy2 = func_deform((x2-min_x)/range_x, (y2-min_y)/range_y)
    x1 = xy1[[1]]
    x2 = xy2[[1]]
    y1 = xy1[[2]]
    y2 = xy2[[2]]
  }
  
  nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)
  
  if(plot_pdf) {
    pdf(fn,width = w,height =h,useDingbats = F)
  } else {
    png(fn, width = w,height = h,bg = bg)
  }
  
  f_overflow = nn$type2=="norm_f" & nn$cost > 100
  f_underflow = nn$type2 == "norm_f" & nn$cost < -100
  #  nn$flow/(1e-8 + nn$capacity) < miss_color_thresh
  if(is.null(propogate) & is.null(mc_t_score)) {
    
    if(show_axes) {
      plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    } else {
      plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex,axes = F,xlab = "",ylab = "")
    }
    
    mc_rgb = col2rgb(mc@colors)/256
    f = nn$mc1>0 & nn$mc2 > 0
    m1 = as.numeric(nn$mc1[f])
    m2 = as.numeric(nn$mc2[f])
    seg_df = data.frame(x1 = x1[f], y1=y1[f], dx=x2[f]-x1[f], dy=y2[f]-y1[f], 
                        r1 = mc_rgb["red",m1],
                        r2 = mc_rgb["red",m2],
                        g1 = mc_rgb["green",m1],
                        g2 = mc_rgb["green",m2],
                        b1 = mc_rgb["blue",m1],
                        b2 = mc_rgb["blue",m2])
    
    for(alpha in seq(0,0.98,0.02)) {
      beta = alpha
      beta5 = alpha+0.02
      if(sigmoid_edge) {
        beta = plogis(alpha,loc=0.5,scale=0.1)
        beta5 = plogis(alpha+0.02,loc=0.5,scale=0.1)
      }
      sx1 = seg_df$x1+alpha*seg_df$dx
      sx2 = seg_df$x1+(alpha+0.02)*seg_df$dx
      sy1 = seg_df$y1+beta*seg_df$dy
      sy2 = seg_df$y1+beta5*seg_df$dy
      alpha_col = ifelse(grad_col_edge, alpha,0)
      rgb_r = seg_df$r2*alpha_col+seg_df$r1*(1)
      rgb_g = seg_df$g2*alpha_col+seg_df$g1*(1)
      rgb_b = seg_df$b2*alpha_col+seg_df$b1*(1)
      cols = rgb(rgb_r, rgb_g, rgb_b)
      segments(sx1, sy1, sx2, sy2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, max_lwd))
    }
    if(show_over_under_flow) {
      f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
                               lwd=pmin(nn$flow[f]/edge_w_scale, max_lwd))
      f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
                                lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,max_lwd))
    }
    
    points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=1)
    if(show_over_under_flow) {
      f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
                               lwd=pmin(nn$flow[f]/edge_w_scale, max_lwd))
      f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
                                lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,max_lwd))
    }
    
  } else if(!is.null(mc_t_score)) {
    plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    mc_t_score = pmax(mc_t_score, quantile(mc_t_score,0.03))
    mc_t_score = pmin(mc_t_score, quantile(mc_t_score,0.97))
    mc_t_score = mc_t_score-min(mc_t_score)
    mc_t_score = mc_t_score/max(mc_t_score)
    mc_t_score = floor(1+999*mc_t_score)
    f = nn$mc1>0 & nn$mc2>0 & nn$time1>0
    max_mc = nrow(mc_t_score)
    m1 = as.numeric(nn$mc1[f]) 
    m2 = as.numeric(nn$mc2[f]) 
    score1 = rep(1, nrow(nn))
    score2 = rep(1, nrow(nn))
    score1[f] = mc_t_score[m1+(nn[f,"time1"]-1)*max_mc]
    score2[f] = mc_t_score[m2+(nn[f,"time2"]-1)*max_mc]
    
    seg_df = data.frame(x1 = x1, y1=y1, dx=x2-x1, dy=y2-y1, 
                        score1= score1, dscore=score2-score1)
    
    for(alpha in seq(0,0.95,0.05)) {
      x1 = seg_df$x1+alpha*seg_df$dx
      x2 = seg_df$x1+(alpha+0.05)*seg_df$dx
      y1 = seg_df$y1+alpha*seg_df$dy
      y2 = seg_df$y1+(alpha+0.05)*seg_df$dy
      cols = score_shades[floor(seg_df$score1 + (alpha+0.025)*seg_df$dscore)]
      segments(x1, y1, x2, y2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, max_lwd))
    }
    rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
         col = mc@colors[nn$mc1], border=NA)
    rect(seg_df$x2-0.12,seg_df$y2-0.5, seg_df$x2+0.12, seg_df$y2+0.5,
         col = mc@colors[nn$mc2], border=NA)
  } else {
    
    #mc_col_bg = alpha(c(c("gray","gray"),mc@colors),0.01)
    #names(mc_col_bg) = c(c("-2","-1"),as.character(c(1:length(mc@colors))))
    
    
    plot(c(x1,x2), c(y1,y2), pch=19, col=bg_col,cex=mc_cex*bg_scale*0.2,axes = F,xlab = "",ylab = "")
    
    segments(x1,y1,x2,y2,
             col=bg_col,
             lwd=pmin(nn$flow/edge_w_scale, max_lwd)*bg_scale)
    
    #plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    max_time = length(propogate)
    m1 = as.numeric(nn$mc1) 
    m2 = as.numeric(nn$mc2) 
    max_m = ncol(propogate[[1]])
    prop_flow = rep(0, nrow(nn))
    for(t in 1:max_time) {
      f = (nn$time1 == t) & nn$mc1>0 & nn$mc2>0
      prop_flow[f] = propogate[[t]][m1[f]+max_m*(m2[f]-1)]
    }
    
    prop_flow[prop_flow/edge_w_scale < 1] = 0
    
    segments(x1,y1,x2,y2, 
             col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
             lwd=pmin(prop_flow/edge_w_scale, max_lwd)*fr_scale)
    points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex*rep(pmin(prop_flow/edge_w_scale, max_lwd),2)*0.1*fr_scale)
  }
  
  if(plot_mc_ids) {
    f1 = nn$type1!="growth" 
    text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
    #	  text(c(x1[f1],x2[f2]),c(y1[f1],y2[f2]),labels = c(nn$mc1[f1],nn$mc2[f2]), cex=1)
  }
  
  dev.off()
}