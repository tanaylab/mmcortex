### Utility functions for plots

#' Plot mc+cell graph using pre-defined mc colorization
#'
#' @param mc2d_id mc2d object to plot
#' @param legend_pos position of legend
#' @param plot_edges plot edges between metacells (true by default)
#' @param min_edge_l (defulat 0) length of edges that are consider long
#' @param edge_w width of long edges
#' @param short_edge_w with of short edges
#' @param show_mcid should metacell id be plotted
#' @param cell_outline should single cell be drawn with outline
#' @param sc_cex size of single cell points
#' @param fig_fn override the default file name in the scfig directory
#' @param filt_mc optionally a factor defining which MCs to plot
#' @param colors can provide color vector to overide the default mc@colors
#'
#' @export
my_mcell_mc2d_plot = function(mc2d_id, device = NULL, legend_pos="topleft", plot_edges=T, min_edge_l=1, edge_w = 1, short_edge_w=0, show_mcid = T, cell_outline=F, colors=NULL, fig_fn = NULL, fn_suf="", sc_cex=1, filt_mc=NULL)
{
    if (!is.null(device)) {
        if (device %in% c('png', 'svg', 'pdf')) {
            tgconfig::set_param(param = 'mc_plot_device', value = device, package = 'metacell')
        } else {
            stop('device must be one of png, svg, pdf')
        }
    }
    mcp_2d_height = tgconfig::get_param("mcell_mc2d_height", package = 'metacell')
	mcp_2d_width = tgconfig::get_param("mcell_mc2d_width", package = 'metacell')
	mcp_2d_plot_key = tgconfig::get_param("mcell_mc2d_plot_key", package = 'metacell')
	mcp_2d_cex = tgconfig::get_param("mcell_mc2d_cex", package = 'metacell')
	mcp_2d_legend_cex = tgconfig::get_param("mcell_mc2d_legend_cex", package = 'metacell')
	mc2d = scdb_mc2d(mc2d_id)
	if(is.null(mc2d)) {
		stop("missing mc2d when trying to plot, id ", mc2d_id)
	}
	mc = scdb_mc(mc2d@mc_id)
	if(is.null(mc)) {
		stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
	}
	if(!is.null(filt_mc)) {
		f_sc = filt_mc[mc@mc[names(mc2d@sc_x)]]
		mc2d@sc_x[!f_sc] = NA
		mc2d@sc_y[!f_sc] = NA
		mc2d@mc_x[!filt_mc] = NA
		mc2d@mc_y[!filt_mc] = NA
	}
	# if(is.null(fig_fn)) {
	# 	fig_fn = scfigs_fn(paste(mc2d_id,fn_suf,sep=""), ifelse(plot_edges, "2d_graph_proj", "2d_proj")) 
	# }
	# .plot_start(fig_fn, w=mcp_2d_width, h = mcp_2d_height)
	# png(fig_fn, width = mcp_2d_width, height = mcp_2d_height)
	if(is.null(colors)) {
		cols = mc@colors
	} else {
		cols = colors
	}
	cols[is.na(cols)] = "gray"
	if(cell_outline) {
		plot(mc2d@sc_x, mc2d@sc_y, pch=21, bg=cols[mc@mc[names(mc2d@sc_x)]], cex=sc_cex, lwd=0.5,
        xlab = 'UMAP 1', 
        ylab = 'UMAP 2',
        xaxt = 'n',
        yaxt = 'n')
	} else {
		plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=cols[mc@mc[names(mc2d@sc_x)]], cex=sc_cex,
        xlab = 'UMAP 1', 
        ylab = 'UMAP 2',
        xaxt = 'n',
        yaxt = 'n')
	}
	fr = mc2d@graph$mc1
	to = mc2d@graph$mc2
	if (plot_edges) {
		dx = mc2d@mc_x[fr]-mc2d@mc_x[to]
		dy = mc2d@mc_y[fr]-mc2d@mc_y[to]
		f = sqrt(dx*dx+dy*dy) > min_edge_l
		segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
					lwd=ifelse(f, edge_w, short_edge_w))
	}
	points(mc2d@mc_x, mc2d@mc_y, cex= mcp_2d_cex, col="black", pch=21, bg=cols)
	if(show_mcid) {
		text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)
	}

# 	if(nrow(mc@color_key)!=0 & mcp_2d_plot_key) {
# 		key = mc@color_key[ mc@color_key$color %in% mc@colors, ]
# #		if(nrow(key!=0)) {
# 		if(!is.null(key) & is.vector(key) & nrow(key) != 0) {
# #group	gene	color	priority	T_fold
# 		gmark = tapply(key$gene, key$group, paste, collapse=", ")
# 		gcol = unique(data.frame(col=key$color, group=key$group))
# 		rownames(gcol) = gcol$group
# 		if(is.vector(gmark)) {
# 			gmark = gmark[order(names(gmark))]
# 		}
# 		if(legend_pos == "panel") {
# 			dev.off()
# 			fig_fn = sub(".png", ".2d_proj_legend.png", fig_fn)
# 			png(fig_fn, width = 600, height= length(gmark)*40+400)
# 			plot.new()
# 			legend_pos = "topleft"
# 		}
# 		legend(legend_pos,
# 				legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)),
# 				pch=19, cex=mcp_2d_legend_cex,
# 				col=as.character(gcol[names(gmark), 'col']), bty='n')
# 		}
# 	}

# 	dev.off()
}

