#Last edit: Jan-11, 2022
#Author: Pon Ganish Prakash
#FileName: External_CIN_Analysis.R

#### TO analyze the Microarray dataset of laser-captured epithelium from 128 cervical tissue specimens      ######
### representing all clinically recognized                                                                  ######
### stages of cervical intra epithelial neoplasias (CIN) (Johan A. den Boon et.al 2015, PMID: 26056290)     ######
### The data can be downloaded from the GEO (accession no. GSE63514)                                        ######

# For more information, please refer to the methods section of this manuscript.






############################### #Section-01  ###############################
#This section includes chunks for QC, Pre-processing and Differential gene expression analysis of CIN Microarray data set

# Preamble
rm(list=ls())
library(affy)
library(oligo)
library(dplyr)
library(gcrma)
library(limma)
library(viridis)
#library(pheatmap)
library(ComplexHeatmap)
library(xtable)
library(knitr)
library(AnnotationDbi)
library(hgu133plus2.db)
library(readxl)
library(reshape2)
library(xtable)
library(openxlsx)
suppressPackageStartupMessages({library("maEndToEnd")})
celpath = "GSE63514_Filtered_CIN/"







# Input raw files and collecting metadata
list = list.files(celpath,full.names=TRUE)
raw_data <- oligo::read.celfiles(list) #From OLIGO
raw_data@phenoData
head(Biobase::pData(raw_data))
ph = raw_data@phenoData


#View the PhenoData
head(Biobase::pData(raw_data),100)
Biobase::exprs(raw_data)[1:5, 1:5]
filenames <- sampleNames(raw_data)



#Add additional info  to the metadata
sampleNames <- sub(".gz", "", filenames)
sampleNames <- sub(".CEL", "", sampleNames)
sampleNames <- sub("_U133_Plus2", "", sampleNames)

pData(raw_data)$filenames <- filenames
sampleNames(raw_data) = sampleNames
rownames(Biobase::pData(raw_data))
sample_size = rep(c('Normal','CIN_1','CIN_2','CIN_3'),times=c(24,14,22,40)) # count of samples
pData(raw_data)$Condition <- sample_size
pData(raw_data) %>%head(50)

raw_data@phenoData
head(Biobase::pData(raw_data),50)
tail(Biobase::pData(raw_data),50)
Biobase::exprs(raw_data)[1:5, 1:5]
raw_data@annotation



#Quality check and Filtering

#Visualize samples in PCA space Before Normalization
exp_raw <- Biobase::exprs(raw_data)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],Condition = pData(raw_data)$Condition)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Condition, colour = Condition), size = 2) +
  ggtitle("PCA plot: Raw Expression for 100 samples") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = c('red','darkgreen','orange','blue'))+theme_minimal()


#visualize Probe Intensities
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for 100 samples:raw data")






#Background correction and Normalization using RMA (robust multi-chip average) algorithm

#checking Relative Log Expression data
no_norm_eset <- oligo::rma(raw_data, normalize = FALSE)
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(no_norm_eset)))
RLE_data <- sweep(Biobase::exprs(no_norm_eset), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = "red", 
                                   angle = 60, size = 6.5, hjust = 1))
dev.off()



#RMA Norm
eset_norm <- oligo::rma(raw_data)
#my.affinity.info = compute.affinities.local(data)
#data.gcrma = gcrma(data,affinity.info=my.affinity.info)
eset_calib <- Biobase::exprs(eset_norm)
PCA <- prcomp(t(eset_calib), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Condition = Biobase::pData(eset_norm)$Condition)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Condition, colour = Condition), size = 2) +
  ggtitle("PCA plot: calibrated-normalized data for 100 samples") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = c('red','darkgreen','orange','blue'))+theme_minimal()




annotation_for_heatmap <- data.frame(Condition = sample_size)
row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))
dists <- as.matrix(dist(t(eset_calib), method = "manhattan"))
rownames(dists) <- row.names(pData(eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlGn"))(255))
colnames(dists) <- NULL
diag(dists) <- NA
ann_colors <- list(Condition = c('Normal' = "burlywood3", 'CIN_1' = "chartreuse4",'CIN_2' = "magenta",'CIN_3' = "gold"))
pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,cluster_rows = F,cluster_cols = F,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Correlation heatmap for the calibrated samples")




#Visualize after RMA Norm.
oligo::boxplot(eset_norm,target = "core", main = "Boxplot of log2-intensitites after RMA norm")
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(eset_norm)))
RLE_data <- sweep(Biobase::exprs(eset_norm), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = "red", 
                                   angle = 60, size = 6.5, hjust = 1))
dev.off()





#Filtering the outliers based on intensity values, applying a cutoff
eset_medians <- rowMedians(Biobase::exprs(eset_norm))
hist_res <- hist(eset_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensity values", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
man_threshold <- 6 # Try and play around with cutoffs
abline(v = man_threshold, col = "coral4", lwd = 2)
no_of_samples <- table(paste0(pData(eset_norm)$Condition)) #for referring count post-selection
no_of_samples 
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)


#Now filter out all transcripts that do not have intensities greater than our selected threshold
eset_manfiltered <- subset(eset_norm, idx_man_threshold)





#Mapping Gene names, symbols and Entrez ID
featureNames(eset_manfiltered) #To see all the manufacturer probe names
eset_anno <- AnnotationDbi::select(hgu133plus2.db, #This is latest annotation, fits well with affy website excel files!!!!
                                       keys = (featureNames(eset_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME","ENTREZID"),
                                       keytype = "PROBEID")

eset_anno_res <- eset_anno #stash it in a dummy variable
eset_anno <- subset(eset_anno, !is.na(SYMBOL))



#To get duplicate probe counts
anno_grouped <- group_by(eset_anno, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
probe_stats <- anno_filtered 
nrow(probe_stats)



#Excluding all the duplicates from our data
ids_to_exlude <- (featureNames(eset_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)
eset_final <- subset(eset_manfiltered, !ids_to_exlude) #Excluded duplicates, final object
validObject(eset_final)




#Adding Feature data with newly mapped gene ID information to the object
head(eset_anno)
fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- left_join(fData(eset_final), eset_anno)
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID 
validObject(eset_final)
eset_final_exp <- Biobase::exprs(eset_final) #Exp values of the final object



#Stashing data
#save(raw_data, eset_final, eset_final_exp, file = 'MainPPQC_data.Rdata')





#DGE Analysis using Limma package
individual <- as.character(Biobase::pData(eset_final)$index)
disease <-factor(as.character(Biobase::pData(eset_final)$Condition), levels = c('Normal','CIN_1','CIN_2','CIN_3'))
design_eset_normal <- model.matrix(~ 0  + disease)
colnames(design_eset_normal) <- c('Normal','CIN_1','CIN_2','CIN_3') 
design_eset_normal
head(design_eset_normal[, 1:2])


cont.matrix <- makeContrasts(CIN_1vsNormal=CIN_1-Normal, #Formula for contrasting all as 2 condition
                             CIN_2vsNormal=CIN_2-Normal,CIN_3vsNormal=CIN_3-Normal,
                             CIN_1vsCIN_2=CIN_1-CIN_2,CIN_1vsCIN_3=CIN_1-CIN_3,CIN_2vsCIN_3=CIN_2-CIN_3,
                             levels=design_eset_normal)
fit <- lmFit(eset_final, design = design_eset_normal)
#cont.matrix <- makeContrasts(CIN_1vsNormal=CIN_1-Normal,CIN_1vs_BothCIN2_3 = CIN_1-(CIN_2+CIN_3)/2, #Formula for contrasting 3 condition
# CIN_2vsNormal=CIN_2-Normal,CIN_3vsNormal=CIN_3-Normal,
# CIN_1vsCIN_2=CIN_1-CIN_2,CIN_1vsCIN_3=CIN_1-CIN_3,CIN_2vsCIN_3=CIN_2-CIN_3,
#levels=design_eset_normal)


fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
all_designs = list()
all_fits = list()
all_results = list()
for (c in colnames(cont.matrix)) {
  res = topTable(fit2, coef=c, adjust="BH", number=nrow(eset_final))
  all_results[[c]] = res
}
str(all_results)
all_results$CIN_1vsNormal %>% head()



all_fits[["global_model"]] = fit2
all_designs[["global_model"]] = list(samples=rownames(eset_final), details=list(Disease = disease), 
                                     design = design_eset_normal, global=T, cont.matrix = cont.matrix)

str(all_fits)
str(all_designs)

all_target_conditions = names(all_results)
par(mfrow=c(2,5))
for (tc in all_target_conditions) {
  r = all_results[[tc]]
  plot(r$logFC, -log10(r$adj.P.Val),xlab="log2 Fold Change",ylab="-log10(adj. p-val)", ylim=c(0,max(2,max(-log10(r$adj.P.Val),na.rm=T))))
  title(main=tc, sub=paste("(",nrow(subset(r, adj.P.Val < 0.05))," signif. DE genes)",sep="") )
  abline(h=-log10(0.05),col="red")
  
}

signif_genes_list = lapply(all_results, function(x) {subset(x,P.Value<0.05)$SYMBOL})
combinations = combn(names(signif_genes_list),2)
pairwise_intersections = apply(combinations,2, function(x) { intersect(signif_genes_list[[x[1]]], signif_genes_list[[x[2]]])} )
d = as.data.frame(t(combinations))
d$overlap = lapply(pairwise_intersections, length)
comp_names = names(signif_genes_list)
d_mat = matrix(NA, nrow=length(comp_names), ncol=length(comp_names))
colnames(d_mat) = comp_names
rownames(d_mat) = comp_names


tmp = apply(d, 1, function(x) {x = as.vector(unlist(x)); d_mat[ x[1],x[2] ] <<- as.numeric(x[3]); d_mat[ x[2],x[1] ] <<- as.numeric(x[3]) })
rm(tmp)
ll = unlist(lapply(signif_genes_list, length))
for (n in names(ll)) d_mat[n, n] <- ll[n]
totals_m = matrix(rep(ll, times=nrow(d_mat)), nrow(d_mat), ncol(d_mat))
d_mat_prop = d_mat/totals_m
d$pairwise_union_lengths = apply(combinations,2, function(x) { length(union(signif_genes_list[[x[1]]], signif_genes_list[[x[2]]]))} )
d$jaccard_coef = as.numeric(d$overlap)/d$pairwise_union_lengths


d_t = d[,c(1,2,3,5)]
colnames(d_t) = c("Condition_1", "Condition_2", "Overlap of DE genes", "Jaccard coefficient")
pheatmap(log10(d_mat+1))


all_DE_results_tmp = list()
for (tc in all_target_conditions) {
  tmp = all_results[[tc]]
  tmp$condition = tc
  all_DE_results_tmp[[tc]] = tmp
}
all_DE_results_ts = do.call(rbind, all_DE_results_tmp)
all_DE_results_ts$DE_class = ifelse(all_DE_results_ts$P.Value >0.05, "n.s.", ifelse(all_DE_results_ts$logFC > 0,"Up","Down"))
agg_fun = function(x) ifelse("Down" %in% x, "Down",ifelse("Up" %in% x, "Up","n.s."))
all_DE_results_sw = dcast(all_DE_results_ts, SYMBOL ~ condition, value.var="DE_class", fun.aggregate=agg_fun)
write.xlsx(all_DE_results_ts, 'All_DE_CIN_Results.xlsx', showNA = TRUE) #This has the proper details for each comparisons
final_de = all_DE_results_ts
final_de[final_de == "" | final_de == "  "] <- NA  # Replace blank & space by NA
length(which(is.na(all_DE_results_ts$SYMBOL))) #out of 2,88,561--45178 gene names are NA




#Export DGE results as excel files
result_folder = "./DGE_CIN/"
all_output_txt_files = list()
all_output_excel_files = list()
output_file_prefix = paste(result_folder,"Differential_expression_results_", sep="/")
selected_cols = c("PROBEID", "SYMBOL", "GENENAME", "ENTREZID","logFC","AveExpr","t","P.Value","adj.P.Val")
for (tc in all_target_conditions) {
  filename = paste(output_file_prefix, tc, ".txt", sep="" )
  write.table(all_results[[tc]][,selected_cols], file= filename, row.names=F , sep="\t", dec=".", quote = F)
  all_output_txt_files[[paste("DGE",tc)]] = filename
}



#Stash
#save(all_DE_results_ts, all_DE_results_sw, file = 'DEG_Results.Rdata')




#Visualize results: Generating HeatMap 
str(all_results)
all_results$CIN_1vsNormal
cin1_up = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_1vsNormal' & all_DE_results_ts$DE_class == 'Up',]
cin1_up$Signature = 'Up-CIN1_vs_Normal'
dim(cin1_up) 
cin2_up = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_2vsNormal' & all_DE_results_ts$DE_class == 'Up',]
cin2_up$Signature = 'Up-CIN2_vs_Normal'
dim(cin2_up) 
cin3_up = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_3vsNormal' & all_DE_results_ts$DE_class == 'Up',]
cin3_up$Signature = 'Up-CIN3_vs_Normal'
dim(cin3_up) 



cin1_down = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_1vsNormal' & all_DE_results_ts$DE_class == 'Down',]
cin1_down$Signature = 'Down-CIN1_vs_Normal'
dim(cin1_down) 
cin2_down = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_2vsNormal' & all_DE_results_ts$DE_class == 'Down',]
cin2_down$Signature = 'Down-CIN2_vs_Normal'
dim(cin2_down) 
cin3_down = all_DE_results_ts[all_DE_results_ts$condition == 'CIN_3vsNormal' & all_DE_results_ts$DE_class == 'Down',]
cin3_down$Signature = 'Down-CIN3_vs_Normal'
dim(cin3_down) 




cin_curated = rbind(cin1_up, cin2_up, cin3_up,cin1_down,cin2_down,cin3_down)
dim(cin_curated) 
cin_curated = subset(cin_curated, select = c('PROBEID','Signature'))
dim(cin_curated) 
cin_curated_final = cin_curated %>% distinct(PROBEID, .keep_all= TRUE)
dim(cin_curated_final)  
table(cin_curated_final$Signature)

cin_curated_final$Signature = factor(cin_curated_final$Signature, levels = c("Up-CIN1_vs_Normal", 
                                                                             "Up-CIN2_vs_Normal",
                                                                             "Up-CIN3_vs_Normal",
                                                                             "Down-CIN1_vs_Normal",
                                                                             "Down-CIN2_vs_Normal",
                                                                             "Down-CIN3_vs_Normal"), ordered = TRUE)


str(cin_curated_final)
cin_curated_final = cin_curated_final[complete.cases(cin_curated_final), ]
rownames(cin_curated_final) = cin_curated_final$PROBEID


#Change scaling considering normal samples (Healthy) as baseline
sel_mat = eset_final_exp[cin_curated_final$PROBEID,]
dim(sel_mat) 
na_rows = apply(is.na(sel_mat), 1, sum)>0
sel_mat = sel_mat[!na_rows,]
dim(sel_mat)
baseline_samples = colnames(eset_final_exp)[1:24]
baseline = rowMeans(sel_mat[, baseline_samples], na.rm = T)
sd_row = apply(sel_mat, 1, sd, na.rm=T)
sel_mat_std = sweep(sweep(sel_mat, 1, baseline, "-"), 1, sd_row, "/") 



annotation_for_heatmap <- data.frame(Condition = pData(raw_data)$Condition)
row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))

cin_curated_final2 = subset(cin_curated_final, select = c(Signature))
hmcol <- rev(plasma(100))
ann_colors <- list("Condition" = c('Normal' = "burlywood3", 'CIN_1' = "orange",'CIN_2' = "red",'CIN_3' = "darkred"),
                   "Signature" = c("Up-CIN1_vs_Normal" = 'cadetblue1', 
                                   "Up-CIN2_vs_Normal" = '#6BAED6',
                                   "Up-CIN3_vs_Normal" = '#08519C',
                                   "Down-CIN1_vs_Normal" = 'seagreen1',
                                   "Down-CIN2_vs_Normal" = 'limegreen',"Down-CIN3_vs_Normal" = 'darkgreen'))




breaks_new =c(0, seq(4,10,4/96),11,12,15)
pheatmap(sel_mat, annotation_row = cin_curated_final2, #Absolute Expression
         annotation_col = annotation_for_heatmap,color = plasma(150),breaks = breaks_new,
         annotation_colors = ann_colors,clustering_distance_rows = "correlation",
         legend = TRUE, show_rownames = F, show_colnames = F, use_raster = T,
         cluster_rows = F,cluster_cols =F,gaps_col = c(24,38),
         main = "Abs Exp: DE Genes-CIN Dataset")

dev.off()

breaks_new =c(-2.1, seq(-2,2,4/96) ,2.1)
ComplexHeatmap::pheatmap(sel_mat_std,   #Relative Expression
         annotation_col = annotation_for_heatmap,color = colorRampPalette(c("dodgerblue", "white","deeppink"))(length(breaks_new)),breaks = breaks_new,
         annotation_colors = ann_colors,clustering_distance_rows = "correlation",
         annotation_row = cin_curated_final2,
         legend = TRUE, show_rownames = F, show_colnames = F,
         cluster_rows = F,cluster_cols = F,gaps_col = c(24,38,60), use_raster = T,gaps_row = c(5604),
         main = "Z-scores (control as Baseline):  DE Genes-CIN Dataset")

dev.off()













############################### #Section-02  ###############################
#This section includes chunks for generating figures 3-i,3-j and 3-k in the Manuscript 


#For Fig.3-i:
#Stemness and differentitation related signatures were collected from excel files (Supplementary Table-7) of DEG (p-value <= 0.05 and log2 FC <−1 or >1) 
#between ectocervical stem cell and differentiated organoids as described in:
#Chumduri, C., Gurumurthy, R.K., Berger, H. et al. Opposing Wnt signals regulate cervical squamocolumnar homeostasis and 
#emergence of metaplasia. Nat Cell Biol (2021). https://doi.org/10.1038/s41556-020-00619-0



#I have Split the excel files into two based on logFC > 1 and logFC < -1
ncb_filtered_early = read_excel('Ncb_early_markers_Log2FC_1.xlsx', skip =1) #Filtered for logFC >1 < -1
ncb_filtered_diff = read_excel('Ncb_Diff_markers_Log2Fc_1.xlsx', skip =1) #Filtered for logFC >1 < -1

#ncb_filtered_early = ncb_filtered_early[ncb_filtered_early[, "logFC"] > 1.5, ]
#ncb_filtered_diff = ncb_filtered_diff[ncb_filtered_diff[, "logFC"] < -1.5, ]


ncb_early_genes = ncb_filtered_early$GeneSymbol
ncb_diff_genes = ncb_filtered_diff$GeneSymbol
ncb_early_genes %>% length() #2060
ncb_diff_genes %>% length() #1779

final_diff = na.omit(ncb_diff_genes)
final_diff%>%length() #1779

final_early = na.omit(ncb_early_genes)
final_early%>% length() #2060

ncb_combined = as.data.frame(final_early)
ncb_combined$Signature = 'Stem cell Profile'
colnames(ncb_combined)[1] <- "GeneSymbols"
ncb_early = ncb_combined

ncb_combined = as.data.frame(final_diff)
ncb_combined$Signature = 'Differentiated Profile'
colnames(ncb_combined)[1] <- "GeneSymbols"
ncb_diff = ncb_combined

dim(ncb_early) #2060 x 2
dim(ncb_diff) #1779 x 2
ncb_curated = rbind(ncb_early, ncb_diff)
dim(ncb_curated) #3839 x 2


ncb_heatmap = ncb_curated %>% distinct(GeneSymbols, .keep_all= TRUE)
dim(ncb_heatmap) #3329 x2 #Perfect!!!
ncb_heatmap$Signature = factor(ncb_heatmap$Signature, levels = c("Stem cell Profile","Differentiated Profile"),
                               ordered = TRUE)
str(ncb_heatmap)
rownames(ncb_heatmap) = ncb_heatmap$GeneSymbols





#Checking for the Expression patterns of the selected geneset in Organoid dataset across different infection conditions
load('MB194_micro_array_preprocessed_data.Rdata')
sel_genes = unique(subset(normalized$genes,normalized$genes$GeneSymbol %in% ncb_heatmap$GeneSymbols)$GeneSymbol) #1687 (1435)
analysis_name = "3D"

sel_comparisons = c("E-3D_Ecto_noE6E7_5d_p.i._vs_UI", "F-3D_Ecto_E6E7_5d_p.i._vs_UI","M_3D_Ecto_5dpi_E6E7_vs_UInoE6E7",
                    "K_3D_Ecto_UI_E6E7_vs_noE6E7","L_3D_Ecto_5dpi_E6E7_vs_noE6E7")

baseline_condition = "3D-UI-5d-noE6E7-Ecto"
sel_cells = "3D"
sel_timepoint = "5d"
ctr_inf_effect_comp = "E-3D_Ecto_noE6E7_5d_p.i._vs_UI"
HPV_effect_comp = "K_3D_Ecto_UI_E6E7_vs_noE6E7"

description_text = "3D Ecto cells (NI, 5d p.i. C.tr, 5d p.i. C.tr + E6/E7, E6/E7 alone)"

#To perform z-scoring
get_ge_mat_standardized_gene <- function(genes, sel_samples = NULL, baseline_cond=baseline_condition) {
  sel_probe_flag = normalized$genes$GeneSymbol %in% genes
  tmp_mat = normalized$E[sel_probe_flag,]
  tmp_mat_avg = avereps(tmp_mat, normalized$genes$GeneSymbol[sel_probe_flag])
  
  if(is.null(sel_samples)) sel_samples = rownames(subset(ed, Site== "Ecto" & cells==sel_cells & Time.p.i. %in% sel_timepoint  & 
                                                           Patient == "hc63"))
  
  sel_mat = tmp_mat_avg[, sel_samples]
  na_rows = apply(is.na(sel_mat), 1, sum)>0
  sel_mat = sel_mat[!na_rows,]
  
  baseline_samples = rownames(subset(ed[sel_samples,], Condition == baseline_cond ))
  baseline = rowMeans(sel_mat[, baseline_samples], na.rm = T)
  sd_row = apply(sel_mat, 1, sd, na.rm=T)
  sel_mat_std = sweep(sweep(sel_mat, 1, baseline, "-"), 1, sd_row, "/") 
  ed_sel = ed[sel_samples,]
  return(list(sel_mat_std, ed_sel))
}
tmp = get_ge_mat_standardized_gene(sel_genes)
sel_mat_std = tmp[[1]]
ed_sel = tmp[[2]]



ed_sel$HPV = factor(ed_sel$HPV, levels=c("noE6E7", "E6E7"))
ed_sel$Infection = factor(ed_sel$Infection, levels=c("UI","Ctr"))
colorder = rownames(ed_sel[order(ed_sel$HPV,ed_sel$Infection) ,])
col_anno= ed_sel[, c( "Infection","HPV")]
col_anno = col_anno[match(colorder, rownames(col_anno)),]
anno_colors = list("HPV"=c("E6E7"="#FF0000","noE6E7"="#00A08A"), "Infection" = 
                     c("UI" = "plum","Ctr" = "magenta"),
                   "cells" = c("3D" = "gray"),
                   "Signature" = c("Differentiated Profile" = "darkgreen",
                                   "Stem cell Profile" = "lawngreen"))


target = rownames(sel_mat_std)
dim(sel_mat_std)
dim(ncb_heatmap)
ncb_heatmap = ncb_heatmap[ncb_heatmap$GeneSymbols %in% target,] #remove excess genes as we will be using avereps the genes count will not match
dim(ncb_heatmap)
row_anno = subset(ncb_heatmap, select = (Signature))
row_anno$Signature = factor(row_anno$Signature, 
                            levels = c("Stem cell Profile", 
                                       "Differentiated Profile"), ordered = TRUE)


sel_mat_std_new <- sel_mat_std[rownames(ncb_heatmap),,drop=FALSE]
anno_colors = list("HPV"=c("E6E7"="brown","noE6E7"="black"), 
                   "Infection" = c("UI" = 'yellow2',"Ctr" = "#FF0000"),
                   "Signature" = c("Differentiated Profile" = "darkgreen",
                                   "Stem cell Profile" = "limegreen"))

#breaks_new =c(min(sel_mat_std_new),-1.1, seq(-1,2,4/96) ,2.1, max(sel_mat_std_new))
breaks_new =c(-2,-1.4, seq(-1,1,4/96), 2)# FIXED!!
res_corr = ComplexHeatmap::pheatmap(sel_mat_std_new[, colorder], annotation_col = col_anno, scale="none", 
                                    show_rownames = F,gaps_col = c(6), annotation_row = row_anno,cluster_rows =F,
                                    show_colnames = F, cluster_cols = F, breaks=breaks_new, 
                                    main=paste(analysis_name, "only: NCB set- stemness & Differentiation"), 
                                    clustering_distance_rows = "correlation", 
                                    annotation_colors = anno_colors, fontsize_row = 12, gaps_row = 1686,
                                    border_color = NA,color = colorRampPalette(c("dodgerblue", "white",
                                                                                 "deeppink"))(length(breaks_new)),use_raster = TRUE)
dev.off()











#For Fig.3-j:
#Genes that were enriched for the Gene ontolgy (GO) Terms: 'cornification', 'epidermal cell differentiation',
#'keratinization','keratinocyte differentiation','epidermis development' and 'skin development', were collected in order to
#check for their expression levels under Ctr, HPV and Co-infection conditions 
#(ref. Suppl.Fig 3h or Suppl. Table-4 in this manuscript for GO Terms)
#Also to understand the gene expression patterns in terms of synergistic effect established during infections.
my_data <- read_excel('Supplementary Table 4.xlsx', skip =1)
newdata = subset(my_data, my_data$Description %in% c('cornification', 'epidermal cell differentiation',
                                                     'keratinization','keratinocyte differentiation',
                                                     'epidermis development','skin development'))

cornified = subset(newdata, newdata$Description == 'cornification')
cornified_1 = cornified[1,]
cornified_genes_1 = strsplit(cornified_1$GeneSymbols, ",")[[1]]
length(cornified_genes_1) #15
cornified_2 = cornified[2,]
cornified_genes_2 = strsplit(cornified_2$GeneSymbols, ",")[[1]]
length(cornified_genes_2) #12

epid_cell_diff = subset(newdata, newdata$Description == 'epidermal cell differentiation')
epid_cell_diff_1 = epid_cell_diff[1,]
epid_cell_diff_genes_1 = strsplit(epid_cell_diff_1$GeneSymbols, ",")[[1]]
length(epid_cell_diff_genes_1) #22
epid_cell_diff_2 = epid_cell_diff[2,]
epid_cell_diff_genes_2 = strsplit(epid_cell_diff_2$GeneSymbols, ",")[[1]]
length(epid_cell_diff_genes_2) #15

keratinization = subset(newdata, newdata$Description == 'keratinization')
keratinization_1 = keratinization[1,]
keratinization_genes_1 = strsplit(keratinization_1$GeneSymbols, ",")[[1]]
length(keratinization_genes_1) #17
keratinization_2 = keratinization[2,]
keratinization_genes_2 = strsplit(keratinization_2$GeneSymbols, ",")[[1]]
length(keratinization_genes_2) #14


krt_diff = subset(newdata, newdata$Description == 'keratinocyte differentiation')
krt_diff_1 = krt_diff[1,]
krt_diff_genes_1 = strsplit(krt_diff_1$GeneSymbols, ",")[[1]]
length(krt_diff_genes_1) #20
krt_diff_2 = krt_diff[2,]
krt_diff_genes_2 = strsplit(krt_diff_2$GeneSymbols, ",")[[1]]
length(krt_diff_genes_2) #15

skin = subset(newdata, newdata$Description == 'skin development')
skin_1 = skin[1,]
skin_genes_1 = strsplit(skin_1$GeneSymbols, ",")[[1]]
length(skin_genes_1)  #21 (ns.ns.down)
skin_3 = skin[3,]
skin_genes_3 = strsplit(skin_3$GeneSymbols, ",")[[1]]
length(skin_genes_3)  #19 (up.up.down)


epi = subset(newdata, newdata$Description == 'epidermis development')
epi_1 = epi[1,]
epi_genes_1 = strsplit(epi_1$GeneSymbols, ",")[[1]]
length(epi_genes_1) #24
epi_2 = epi[2,]
epi_genes_2 = strsplit(epi_2$GeneSymbols, ",")[[1]]
length(epi_genes_2) #21


go_combined = c(epi_genes_1, epi_genes_2, cornified_genes_1, cornified_genes_2, epid_cell_diff_genes_1, epid_cell_diff_genes_2,
                keratinization_genes_1, keratinization_genes_2,krt_diff_genes_1, krt_diff_genes_2,skin_genes_1, skin_genes_3)

length(go_combined) #215


go_combined #215
isUnique(go_combined)
go_combined_unique = unique(go_combined)
go_combined_unique %>% length() #48



#additional step to remove NA from go combined markers
sel_genes = unique(subset(normalized$genes,normalized$genes$GeneSymbol %in% go_combined_unique)$GeneSymbol) #48
tmp = get_ge_mat_standardized_gene(sel_genes)
sel_mat_std = tmp[[1]]
ed_sel = tmp[[2]]
ed_sel

ed_sel$HPV = factor(ed_sel$HPV, levels=c("noE6E7", "E6E7"))
ed_sel$Infection = factor(ed_sel$Infection, levels=c("UI","Ctr"))
colorder = rownames(ed_sel[order(ed_sel$HPV,ed_sel$Infection) ,])
col_anno= ed_sel[, c( "Infection","HPV")]
col_anno = col_anno[match(colorder, rownames(col_anno)),]
anno_colors = list("HPV"=c("E6E7"="red","noE6E7"="black"), "Infection" = c("UI" = 'gray',"Ctr" = "orange"))

roworder = order(rownames(sel_mat_std))
go_combined_unique
handcurated_order = c("PAFAH1B1","INHBA","LCE3D","CDSN","KRT16","LCE3E","ITGA6","LAMC2",
                      "CASP14","KRT77","KRTAP22-1",
                      "LOR","PCSK6","COL1A1","CTSV","HES5","PI3",
                      "CST6","KLK12","KLK5","IVL","SPRR1B","KRT75","PTHLH","LCE1C","PKP1","DSC2","ZDHHC21",
                      "KRT78",
                      "COL17A1","FLG","ACER1","KRT26","PAX6","TGM5","DCT","KRT3","RPTN","LIPK","COL5A2","KRTAP19-8",
                      "KLK13","WNT5A","CSTA","KRT73","LIPN","EDA2R","KRT85")


#Re-ordering manually to retain patterns for visualization
sel_mat_std_ordered_handcurated = sel_mat_std[match(handcurated_order, rownames(sel_mat_std)),]
anno_colors = list("HPV"=c("E6E7"="brown","noE6E7"="black"), "Infection" = c("UI" = 'yellow2',"Ctr" = "#FF0000"))

breaks_new =c(-2.1, seq(-2,2,4/96) ,2.1)
breaks_new
res_corr = ComplexHeatmap::pheatmap(sel_mat_std_ordered_handcurated[, colorder], 
                    annotation_col = col_anno, scale="none", show_rownames = T,
                    show_colnames = F, cluster_cols = F, breaks=breaks_new, cluster_rows = F,
                    main=paste(analysis_name, "only: GO Term Markers for Corni, diff..Arranged"), 
                    annotation_colors = anno_colors, fontsize_row = 12, border_color = NA,gaps_col = 6,
                    use_raster = T,color = colorRampPalette(c("dodgerblue", "white",
                                                              "deeppink"))(length(breaks_new)))

dev.off()













#For Fig.3-k:
#Collecting differentially regulated genes between hCEcto E6E7 organoids with or without Ctr infection and 
#to check for their relative expression levels in different stages of CIN samples
load('MB194_micro_array_preprocessed_data.Rdata')
load("DGE_analysis_image.Rdata")



#1. Collecting our DE data of interest for UI_E6E7 vs UI_noE6E7(Effect of E6E7)
de_set_oi  = all_results$K_3D_Ecto_UI_E6E7_vs_noE6E7
de_set_oi$Class = ifelse(de_set_oi$adj.P.Val>0.05, "n.s.", 
                         ifelse(de_set_oi$logFC > 0,"Up","Down"))

table(de_set_oi$Class)
topgenes = de_set_oi[de_set_oi[, "P.Value"] < 0.05, ] #all sig.DE genes 14327
topups = topgenes[topgenes[, "logFC"] > 1.5, ]
dim(topups) #up reg. 903
topdowns = topgenes[topgenes[, "logFC"] < -1.5, ] 
dim(topdowns) #down reg.545
toPrint <- subset(de_set_oi,P.Value < 0.05) #All significant DE Genes
up_genes_de_3d_final = topups$GeneSymbol #903
down_genes_de_3d_final = topdowns$GeneSymbol #545


#This is very important for probeid-gene annotations!!
tmp = unique(eset_final@featureData@data[, c("PROBEID","SYMBOL")])
pn2gs = as.data.frame(tapply(tmp$SYMBOL, tmp$PROBEID, function(x) paste(x, collapse=",")))
colnames(pn2gs) = "GeneSymbol"
eset_final #exp matrix
dim(eset_final)
eset_final
a = eset_final@featureData@data
dim(a)


#extract_genes = a$PROBEID[a$SYMBOL %in% up_genes_de_3d_final] #up
new_set = a[match(up_genes_de_3d_final, a$SYMBOL),]
new_set_up = new_set[complete.cases(new_set), ]
extract_genes_up = new_set_up$PROBEID #699

#extract_genes = a$PROBEID[a$SYMBOL %in% down_genes_de_3d_final] #down
new_set = a[match(down_genes_de_3d_final, a$SYMBOL),]
new_set_down = new_set[complete.cases(new_set), ]
extract_genes_down = new_set_down$PROBEID #326



#2. Collecting our DE data of interest for Ctr.5dpi_E6E7_vs_noE6E7(Effect of Co-infection)
de_set_oi2  = all_results$L_3D_Ecto_5dpi_E6E7_vs_noE6E7
de_set_oi2$Class = ifelse(de_set_oi2$adj.P.Val>0.05, "n.s.", 
                          ifelse(de_set_oi2$logFC > 0,"Up","Down"))

table(de_set_oi2$Class)
topgenes2 = de_set_oi2[de_set_oi2[, "P.Value"] < 0.05, ] #all sig.DE genes 15683
topups2 = topgenes2[topgenes2[, "logFC"] > 1.5, ]
dim(topups2) #up reg. 1141
topdowns2 = topgenes2[topgenes2[, "logFC"] < -1.5, ] 
dim(topdowns2) #down reg.788
toPrint2 <- subset(de_set_oi2,P.Value < 0.05) #All significant DE Genes (15,679)
up_genes_de_3d_final2 = topups2$GeneSymbol #1141
down_genes_de_3d_final2 = topdowns2$GeneSymbol #788



#extract_genes = a$PROBEID[a$SYMBOL %in% up_genes_de_3d_final] #up
new_set2 = a[match(up_genes_de_3d_final2, a$SYMBOL),] #1141
new_set_up2 = new_set2[complete.cases(new_set2), ] #901
extract_genes_up2 = new_set_up2$PROBEID #901

#extract_genes = a$PROBEID[a$SYMBOL %in% down_genes_de_3d_final] #down
new_set2 = a[match(down_genes_de_3d_final2, a$SYMBOL),] #788
new_set_down2 = new_set2[complete.cases(new_set2), ] #448
extract_genes_down2 = new_set_down2$PROBEID #448




#3. Collecting our DE data of interest for Ctr.5dpi_E6E7_vs_noE6E7(Effect of Chlamydia infection)
de_set_oi3  = all_results$`E-3D_Ecto_noE6E7_5d_p.i._vs_UI`
de_set_oi3$Class = ifelse(de_set_oi3$adj.P.Val>0.05, "n.s.", 
                          ifelse(de_set_oi3$logFC > 0,"Up","Down"))

table(de_set_oi3$Class)

topgenes3 = de_set_oi3[de_set_oi3[, "P.Value"] < 0.05, ] #all sig.DE genes 15683
topups3 = topgenes3[topgenes3[, "logFC"] > 1.5, ]
dim(topups3) #up reg. 716
topdowns3 = topgenes3[topgenes3[, "logFC"] < -1.5, ] 
dim(topdowns3) #down reg.292
toPrint3 <- subset(de_set_oi3,P.Value < 0.05) #All significant DE Genes (15,679)
up_genes_de_3d_final3 = topups3$GeneSymbol #716
down_genes_de_3d_final3 = topdowns3$GeneSymbol #292


#extract_genes = a$PROBEID[a$SYMBOL %in% up_genes_de_3d_final] #up
new_set3 = a[match(up_genes_de_3d_final3, a$SYMBOL),] 
new_set_up3 = new_set3[complete.cases(new_set3), ] #473
extract_genes_up3 = new_set_up3$PROBEID #473

#extract_genes = a$PROBEID[a$SYMBOL %in% down_genes_de_3d_final] #down
new_set3 = a[match(down_genes_de_3d_final3, a$SYMBOL),] #292
new_set_down3 = new_set3[complete.cases(new_set3), ] #178
extract_genes_down3 = new_set_down3$PROBEID #178




#Combined:UP-DOWN
comb_genes_up = as.data.frame(extract_genes_up)
comb_genes_up$Signature = 'UP_UI_E6E7_vs_noE6E7'
colnames(comb_genes_up)[1] <- "PROBEID"
comb_genes_down = as.data.frame(extract_genes_down)
comb_genes_down$Signature = 'DOWN_UI_E6E7_vs_noE6E7'
colnames(comb_genes_down)[1] <- "PROBEID"
comb_genes_up2 = as.data.frame(extract_genes_up2)
comb_genes_up2$Signature = 'UP_5dpi_E6E7_vs_noE6E7'
colnames(comb_genes_up2)[1] <- "PROBEID"
comb_genes_down2 = as.data.frame(extract_genes_down2)
comb_genes_down2$Signature = 'DOWN_5dpi_E6E7_vs_noE6E7'
colnames(comb_genes_down2)[1] <- "PROBEID"
comb_genes_up3 = as.data.frame(extract_genes_up3)
comb_genes_up3$Signature = 'UP_5dpi_noE6E7_vs_UI_noE6E7'
colnames(comb_genes_up3)[1] <- "PROBEID"
comb_genes_down3 = as.data.frame(extract_genes_down3)
comb_genes_down3$Signature = 'DOWN_5dpi_noE6E7_vs_UI_noE6E7'
colnames(comb_genes_down3)[1] <- "PROBEID"

dim(comb_genes_down) 
dim(comb_genes_up) 
dim(comb_genes_down2) 
dim(comb_genes_up2) 
dim(comb_genes_down3) 
dim(comb_genes_up3) 


comb_curated2 = rbind(comb_genes_up, comb_genes_down, comb_genes_up2, comb_genes_down2,comb_genes_up3, comb_genes_down3)
dim(comb_curated2) #3046x2
table(comb_curated2$Signature)
comb_heatmap2 = comb_curated2 %>% distinct(PROBEID, .keep_all= TRUE)
dim(comb_heatmap2) #1766 x2 #Perfect!!!
table(comb_heatmap2$Signature)


comb_heatmap2$Signature = factor(comb_heatmap2$Signature, levels = c("UP_UI_E6E7_vs_noE6E7",
                                                                     "DOWN_UI_E6E7_vs_noE6E7",
                                                                     "UP_5dpi_E6E7_vs_noE6E7",
                                                                     "DOWN_5dpi_E6E7_vs_noE6E7",
                                                                     "UP_5dpi_noE6E7_vs_UI_noE6E7",
                                                                     "DOWN_5dpi_noE6E7_vs_UI_noE6E7"),
                                 ordered = TRUE)
str(comb_heatmap2)
comb_heatmap2 = comb_heatmap2[complete.cases(comb_heatmap2), ]
rownames(comb_heatmap2) = comb_heatmap2$PROBEID
table(comb_heatmap2$Signature)


comb_heatmap2$Layer = 'Layer'
comb_heatmap2$Layer[1:898] = 'E6E7_vs_noE6E7'
comb_heatmap2$Layer[899:1399] = 'Ctr_5dpi_E6E7_vs_noE6E7'
comb_heatmap2$Layer[1400:1766] = 'Ctr_5dpi_noE6E7_vs_UI_noE6E7'


comb_heatmap2$Class = 'Class'
comb_heatmap2$Class[1:600] = 'Up'
comb_heatmap2$Class[601:898] = 'Down'
comb_heatmap2$Class[899:1182] = 'Up'
comb_heatmap2$Class[1183:1399] = 'Down'
comb_heatmap2$Class[1400:1686] = 'Up'
comb_heatmap2$Class[1687:1766] = 'Down'


comb_heatmap2$Layer = factor(comb_heatmap2$Layer, levels = c("E6E7_vs_noE6E7","Ctr_5dpi_E6E7_vs_noE6E7",
                                                             "Ctr_5dpi_noE6E7_vs_UI_noE6E7"),
                             ordered = TRUE)

comb_heatmap2$Class = factor(comb_heatmap2$Class, levels = c("Up","Down"),
                             ordered = TRUE)


sel_mat = eset_final_exp[comb_heatmap2$PROBEID,]
dim(sel_mat)

na_rows = apply(is.na(sel_mat), 1, sum)>0
sel_mat = sel_mat[!na_rows,]
dim(sel_mat)

baseline_samples = colnames(eset_final_exp)[1:24]
baseline = rowMeans(sel_mat[, baseline_samples], na.rm = T)
sd_row = apply(sel_mat, 1, sd, na.rm=T)
sel_mat_std = sweep(sweep(sel_mat, 1, baseline, "-"), 1, sd_row, "/") 



rownames(Biobase::pData(raw_data))
sample_size = rep(c('Normal','CIN_1','CIN_2','CIN_3'),times=c(24,14,22,40))
pData(raw_data)$Condition <- sample_size



annotation_for_heatmap <- data.frame(Condition = sample_size)
row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))

annotation_for_heatmap$Condition = factor(annotation_for_heatmap$Condition, levels = c("Normal","CIN_1",
                                                                                       "CIN_2","CIN_3"),
                                          ordered = TRUE)



comb_heatmap3 = subset(comb_heatmap2, select = c(Class, Layer))

hmcol <- rev(plasma(100))
brewer.pal(9, 'Blues')
ann_colors <- list("Condition" = c('Normal' = "#C6DBEF", 'CIN_1' = "#9ECAE1",'CIN_2' = "#3182BD",'CIN_3' = "#08519C"),
                   "Layer" = c('E6E7_vs_noE6E7' = "turquoise3",'Ctr_5dpi_E6E7_vs_noE6E7'='#FF0000',
                               'Ctr_5dpi_noE6E7_vs_UI_noE6E7' = 'black'),
                   "Class" = c('Up' = "limegreen",
                               'Down' = "orange"))


breaks_new =c(-1.1, seq(-1,1,4/96) ,1.1)
ComplexHeatmap::pheatmap(sel_mat_std,
                         annotation_col = annotation_for_heatmap,color = colorRampPalette(c("dodgerblue", "white",
                                                                                            "deeppink1"))(length(breaks_new)),
                         breaks = breaks_new,annotation_colors = ann_colors,clustering_distance_rows = "correlation",
                         annotation_row = comb_heatmap3,#right_annotation = ha,
                         legend = TRUE, show_rownames = F, show_colnames = F,row_names_side = "right",
                         cluster_rows = F,cluster_cols = F,gaps_col = c(24,38,60), #gaps_row = c(596,889,1171),
                         use_raster = TRUE,
                         main = "DE GENES-Ectocervix organoids (Effect of Ctr, E6E7 and Coinf.) Mapped to CIN Dataset")

dev.off()



############# END #################################









