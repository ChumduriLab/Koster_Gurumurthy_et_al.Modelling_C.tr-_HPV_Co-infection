#########################################################
# assumes that dge_results is a data.table or data.frame object
dedup_DGE_results <- function(dge_results, by="pvalue", logfc_col = "logFC", pval_col = "P.Value", adjp_col = "adj.P.Val", id_col = "GeneSymbol", probeid_col="ProbeName") {
  probes = dge_results[[probeid_col]]
  symbol = dge_results[[id_col]]
  logfc = dge_results[[logfc_col]]
  pval = dge_results[[pval_col]]
  adjp = dge_results[[adjp_col]]
  
  probes_by_symbol = tapply(probes, symbol, function(x) x)
  if (by=="logfc") {
    max_fc = function(x) {which(abs(x) == max(abs(x),na.rm=T))}
    max_fc_ind_per_symbol = lapply(split(logfc,symbol), max_fc)
    sel_probes = mapply( function(x,y) x[y], probes_by_symbol, max_fc_ind_per_symbol)
  } else if (by == "pvalue") {
    minp = function(x) {which(x == min(x,na.rm=T))}
    minp_ind_per_symbol = lapply(split(pval,symbol), minp)
    sel_probes = mapply( function(x,y) x[y], probes_by_symbol, minp_ind_per_symbol)    
  }
  
  return( dge_results[dge_results[[probeid_col]] %in% unlist(sel_probes), ] )
}

suppressMessages(require(RColorBrewer))
RedWhtBlue_orig = brewer.pal(n = 7, name = "RdBu")
RedGyBlue = RedWhtBlue_orig
RedGyBlue[4] = "grey90"
RedGreyBluePal = colorRampPalette(rev(RedGyBlue))(100)

RedGyBlue2 = brewer.pal(n = 7, name = "RdYlBu")
RedGyBlue2[4] = "grey90"
RedGreyBluePal2 = colorRampPalette(rev(RedGyBlue2))(100)

RedGreyBluePal3 = colorRampPalette(c(brewer.pal(n = 7, name = "RdYlBu")[7],  "grey90" ,brewer.pal(n = 7, name = "RdYlBu")[1]))(100)

RedWhtBlue_orig = brewer.pal(n = 7, name = "RdBu")
RedGyBlue_light = RedWhtBlue_orig
RedGyBlue_light[4] = "grey90"
RedGyBlue_light[1] = "red"
RedGyBlue_light[7] = "blue"
RedGreyBlueLightPal = colorRampPalette(rev(RedGyBlue_light))(100)

darker_rgrb = apply(round(t(col2rgb(RedGyBlue_light)) / 1.5), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
RedGreyBlueDarkPal = colorRampPalette(rev(darker_rgrb))(100)


suppressMessages(library(org.Hs.eg.db))
entrez_to_symbol = function(entrezids) {
  tmp = as.character(entrezids)
  tmp = tmp[!is.na(tmp)]
  tmp_df= dcast(select(org.Hs.eg.db, keys=tmp, keytype="ENTREZID", columns="SYMBOL"), ENTREZID ~ ., value.var="SYMBOL", fun.aggregate = function(x) paste(x[!is.na(x)], collapse=",") )
  colnames(tmp_df) = c("EntrezID", "GeneSymbol")  
  return(tmp_df)
}

symbol_to_entrez = function(symbols) {
  tmp = as.character(symbols)
  tmp = tmp[!is.na(tmp)]
  tmp_df= dcast(select(org.Hs.eg.db, keys=tmp, keytype="SYMBOL", columns="ENTREZID"), SYMBOL ~ ., value.var="ENTREZID", fun.aggregate = function(x) paste(x[!is.na(x)], collapse=",") )
  colnames(tmp_df) = c("GeneSymbol", "EntrezID")  
  return(tmp_df)
}

probename_to_symbol_agilent = function(probenames) {
  tmp = normalized$genes[normalized$genes$ProbeName %in% probenames, c("ProbeName", "GeneSymbol")]
  rownames(tmp) = tmp$ProbeName
  return(tmp[probenames,])
}

probename_to_entrez_agilent = function(probenames) {
  tmp = normalized$genes[normalized$genes$ProbeName %in% probenames, c("ProbeName", "EntrezID")]
  rownames(tmp) = tmp$ProbeName
  return(tmp[probenames,])
}



















# Draws heatmap, expression barplot and fold change barplots (both all and FDR< 10%).
# optionally stores TIFF image of final FC barplot
# sel_genes must be a list of Entrez gene IDs
Undiff_vs_Diff_plots = function(sel_genes, sel_samples, title, tiff=F, comparisons_for_fc = c("ecto_2D_vs_whole_all", "ecto_early_vs_whole_unpaired")) {
  
  sel_genes = as.character(sel_genes)
  
  ### 1. Heatmap
  norm_mat_sel_genes = normalized$E[normalized$genes$EntrezID %in% sel_genes, sel_samples]
  rownames(norm_mat_sel_genes) = normalized$genes[normalized$genes$EntrezID %in% sel_genes,"GeneSymbol"]
  colnames(norm_mat_sel_genes) = paste(ed[colnames(norm_mat_sel_genes), "ShortName"], ed[colnames(norm_mat_sel_genes), "patient_or_mouse_id"], ed[colnames(norm_mat_sel_genes), "Passage"], sep="_")
  norm_mat_sel_genes_sorted = norm_mat_sel_genes[,order(colnames(norm_mat_sel_genes))]
  nmat_scaled = t(scale(t(norm_mat_sel_genes_sorted)))
  nmat_scaled[is.na(nmat_scaled)] = 0
  nmat_scaled[abs(nmat_scaled) > 2] <- sign(nmat_scaled[abs(nmat_scaled) > 2]) * 2
  breaks_new = c(-2, seq(-1,1,2/98), 2)
  pheatmap(nmat_scaled, cluster_rows = T, cluster_cols=F,  scale="none", main = title, breaks=breaks_new, show_rownames=T, fontsize_col=9)
  
  #### Barplot
  if(length(unique(rownames(norm_mat_sel_genes_sorted))) < 40) {
    nmsg_ts = melt(norm_mat_sel_genes)
    nmsg_ts$Organoid.Type = ed_MB194_shortname[as.character(nmsg_ts$Var2),"Organoid.Type"]
    nmsg_ts$Tissue.Type = ed_MB194_shortname[as.character(nmsg_ts$Var2),"Tissue.Type"]
    #nmsg_ts$Celltype = ifelse(nmsg_ts$Tissue.Type=="Endo", "Endocervix Whole", paste(nmsg_ts$Tissue.Type, nmsg_ts$Organoid.Type))
    nmsg_ts$Celltype = paste(nmsg_ts$Tissue.Type, nmsg_ts$Organoid.Type)
    nmsg_agg <- aggregate(nmsg_ts$value,
                          by = list(Celltype = nmsg_ts$Celltype, Gene = nmsg_ts$Var1),
                          FUN = function(x) c(mean = mean(x, na.rm=T), sd = sd(x, na.rm=T),
                                              n = sum(!is.na(x))))
    nmsg_agg_fixed = cbind(nmsg_agg[,1:2], nmsg_agg$x)
    nmsg_agg_fixed$se = nmsg_agg_fixed$sd/sqrt(nmsg_agg_fixed$n)
    nmsg_agg_fixed = nmsg_agg_fixed[naturalorder(as.character(nmsg_agg_fixed$Gene)),]
    nmsg_agg_fixed$Gene = factor(as.character(nmsg_agg_fixed$Gene), levels=naturalsort(unique(as.character(nmsg_agg_fixed$Gene))))
    
    limits <- aes(ymax = mean + se, ymin=mean - se)
    dodge <- position_dodge(width=0.9)
    p = ggplot(nmsg_agg_fixed, aes(fill=Celltype, y=mean, x=Gene)) + geom_bar(position=dodge, stat="identity") + geom_errorbar(limits, position=dodge, width=0.25) + theme(axis.text.x=element_text(angle=-90)) + ylab("Mean expression level") + labs(title=title)
    plot(p)
  }
  ####################################################################
  #### Fold change barplot
  for (c in comparisons_for_fc) {
    tmp = all_results[[c]]
    tmp_sel = subset(tmp, as.character(EntrezID) %in% sel_genes)
    tmp_sel_sorted = tmp_sel[order(tmp_sel$logFC, decreasing=F),]
    par(mar=c(5,6,4,2))
    bcol= ifelse(tmp_sel_sorted$adj.P.Val>0.1, "grey", ifelse(tmp_sel_sorted$logFC > 0, "red","blue"))
    barplot(tmp_sel_sorted$logFC, names.arg = tmp_sel_sorted$GeneSymbol, las=2, cex.names=0.8, ylab="Log2 Fold change", col=bcol, main=paste(title,"[",c,"]"), border=NA)
    legend("topleft", legend=c("Significantly downregulated", "Not significant (FDR>=10%)", "Significant upregulated"), fill=c("blue","grey","red"))
    
    ####################################################################
    # FC barplot, significant genes only, dedupped
    
    tmp_sel = subset(tmp, as.character(EntrezID) %in% sel_genes & adj.P.Val < 0.1)
    if (nrow(tmp_sel) > 0) {
      tmp_sel_dedup = dedup_DGE_results(tmp_sel, by="pvalue", logfc_col = "logFC", pval_col = "P.Value", adjp_col = "adj.P.Val", id_col = "GeneSymbol", probeid_col="ProbeName") 
      tmp_sel_sorted = tmp_sel_dedup[order(tmp_sel_dedup$logFC, decreasing=F),]
      
      par(mar=c(12,6,4,2))
      bcol= ifelse(tmp_sel_sorted$adj.P.Val>0.1, "grey", ifelse(tmp_sel_sorted$logFC > 0, "red","blue"))
      barplot(tmp_sel_sorted$logFC, names.arg = tmp_sel_sorted$GeneSymbol, las=2, cex.names=1, cex.axis = 2, ylab="", col=bcol, main=paste(title,"[",c,"]"), ylim=c(-7, 7))
      mtext("Log2 Fold change", side=2, cex=2, line=3.4)
      #legend("topleft", legend=c("Significantly downregulated", "Significantly upregulated"), fill=c("blue","red"), cex=1.3)
    }
  }
  ####################################################################
  if(tiff) {
    par_orig = par()
    image_file = paste(result_folder, paste(title,"_FC_Barplot_MB194.tiff",sep="") ,sep="/")
    tiff(image_file, width=12*100, height=8 * 100)
    par(mar=c(14,10,4,2))
    barplot(tmp_sel_sorted$logFC, names.arg = tmp_sel_sorted$GeneSymbol, las=2, cex.names=3, cex.axis = 3, ylab="", col=bcol, main=title, ylim=c(-5, 10))
    mtext("log2 Fold change", side=2, cex=3, line=4.5)
    dev.off()
    par(par_orig)
  }
}
























