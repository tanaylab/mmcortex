library("Matrix")

embdyn_cc_plots = function(mat_id, mc_id,  mc2d_id=mc_id, tag=mc_id,m_0 = 0.01,s_0 = 0.005)
{
  # output: function returns names of cells which are growth-arrested, i.e. lye below a threshold line
  
	m_genes = c("Mki67","Cenpf","Top2a","Smc4;SMC4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr;IHABP","Cenpa;Cenp-a","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgol2","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")

	s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1

	m = scdb_mat(mat_id)
	mc = scdb_mc(mc_id)
	mc2d = scdb_mc2d(mc2d_id)
	
	s_genes = intersect(rownames(mc@mc_fp), s_genes)
	m_genes = intersect(rownames(mc@mc_fp), m_genes)

	tot  = colSums(m@mat)
	s_tot = colSums(m@mat[s_genes,])
	m_tot = colSums(m@mat[m_genes,])

	s_score = s_tot/tot
	m_score = m_tot/tot

	if(!dir.exists(sprintf("figs/cc_%s",tag))) {
		dir.create(sprintf("figs/cc_%s", tag))
	}

	# m_score = y_0 (1 - s_score/x_0) 
	
  f = (m_score < m_0 * (1- s_score/s_0))

  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  
	
	# Next try to fit a smooth spline
	
  p_coldens =densCols(x = s_score,y = m_score,colramp = colorRampPalette(c("lightgray","blue3", "red", "yellow")))
  #p_coldens =densCols(x = s_score,y = m_score,colramp = colorRampPalette(c("white",rev(inferno(5)))))
  #p_coldens =densCols(x = s_score,y = m_score)
  	
  png(sprintf("figs/cc_%s/cc_coldens_scores.png", tag), w=600, h=600)
  plot(s_score, m_score, pch=19, main = "S phase vs M phase UMIs",cex=0.4,
       xlab = "S phase score",ylab = "M phase score",
       col = p_coldens)
  dev.off()
  
  png(sprintf("figs/cc_%s/cc_scores.png", tag), w=600, h=600)
  plot(s_score, m_score, pch=19, main = "S phase vs M phase UMIs",cex=0.1,
       xlab = "S phase score",ylab = "M phase score")
  points(s_score[f], m_score[f], pch=19, cex=0.1, col="darkred")
  dev.off()
  shades = colorRampPalette(c("white","lightgray","lightblue","blue", "red", "yellow"))
  png(sprintf("figs/cc_%s/cc_sm_scores.png", tag), w=600, h=600)
  
  smoothScatter(s_score, m_score, colramp=shades, pch=19, cex=0.1,
                main = "S phase vs M phase UMIs",xlab = "S phase score",ylab = "M phase score")
  abline(a = m_0,b = - m_0/s_0)
  dev.off()
  
  
  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  png(sprintf("figs/cc_%s/2d_cc.png", tag), w=800, h=800)
  plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.4, col=ifelse(f[names(mc2d@sc_x)], "black", "lightgray"))
  points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=2.5, bg=shades[mc_cc])
  dev.off()
  
  png(sprintf("figs/cc_%s/bars_cc.png", tag), w=3600, h=1500)
  barplot(mc_cc, col=mc@colors, las=2, cex.names=0.7)
  dev.off()
  
	#return(colnames(m@mat)[f])
}
