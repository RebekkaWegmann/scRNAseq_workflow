## ---- eval = F-----------------------------------------------------------------------------------------------------------------------------------
## install.packages("BiocManager")


## ---- eval = F-----------------------------------------------------------------------------------------------------------------------------------
## install.packages("devtools")
## 
## # on Windows, you might need to run this in case devtools does not find Rtools 3.5 (only works for older versions of devtools though):
## library(devtools)
## assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
## find_rtools() # is TRUE now


## ---- eval=F-------------------------------------------------------------------------------------------------------------------------------------
## # from CRAN
## install.packages("Rtsne") #Rtsne v. 0.15
## install.packages("ggplot2") # ggplot2 v. 3.2.2
## install.packages("data.table") # data.table v. 1.12.8
## install.packages("RColorBrewer") # RColorBrewer v. 1.1-2
## install.packages("mvoutlier") # mvoutlier 2.0.9, required by some functions from scater
## 
## devtools::install_github("bwlewis/irlba") # irlba 2.3.3, optional. Make sure you use irlba > 2.3.1, older versions contain a bug that results in unreliable output!
## 
## # from Bioconductor
## BiocManager::install("scater", version="3.11") # scater v. 1.16.2
## BiocManager::install("scran", version="3.11") # scran v. 1.16.0


## ---- eval=F-------------------------------------------------------------------------------------------------------------------------------------
## # ensembldb 2.12.1 and EnsDb.Hsapiens.v86 2.99.0
## BiocManager::install(c("ensembldb","EnsDb.Hsapiens.v86"), version="3.11")
## # org.Hs.eg.db v. 3.11.4
## BiocManager::install("org.Hs.eg.db", version="3.11")


## ---- eval=F-------------------------------------------------------------------------------------------------------------------------------------
## # ensembldb 2.4.1 and EnsDb.Mmusculus 2.99.0
## biocLite(c("ensembldb","EnsDb.Mmusculus.v79"))
## # org.Mm.eg.db v. 3.4.1
## biocLite("org.Mm.eg.db")


## ----eval = F------------------------------------------------------------------------------------------------------------------------------------
## BiocManager::install("M3Drop", version="3.11")


## ---- eval = F-----------------------------------------------------------------------------------------------------------------------------------
## install.packages("cluster") # cluster v. 2.0.6
## install.packages("dendextend") # dendextend v. 1.5.2
## install.packages("Ckmeans.1d.dp") # Ckmeans.1d.dp v. 4.2.1
## install.packages("dbscan") # DBSCAN v. 1.1-1
## install.packages(c("mclust","mvtnorm")) # mclust v. 5.4 and mvtnorm 1.0.6
## install.packages("dynamicTreeCut") # dynamicTreeCut v. 1.63-1
## 
## BiocManager::install("SC3", version="3.11") # SC3 v. 1.4.2
## install.packages("MCL") #MCL 1.0


## tar xzf mcl-14-137.tar.gz

## cd mcl-14-137

## ./configure --prefix=$HOME/local

## make install


## ----eval = F------------------------------------------------------------------------------------------------------------------------------------
## BiocManager::install("limma", version="3.11") # limma v. 3.32.5


## ---- message = FALSE----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
graphics.off()

wd = ".."

#Directory where input files are stored
input_dir = file.path(wd,"example_data")

#Directory where code is stored
code_dir = file.path(wd,"code")

#where to save the output data?
out_data_dir = file.path(wd,"example_output")
if(!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

#where to save the produced plots?
plotdir = file.path(wd,"example_plots")
if(!dir.exists(plotdir)) {dir.create(plotdir)}

#use extrenal MCL? ONLY POSSIBLE ON UNIX
use_external_mcl=FALSE

#path to your MCL binary (only needed if use_external_mcl = TRUE)
#mcl_path = "~/local/bin/mcl"

set.seed(17) #to make tSNE plots reproducible

source(file.path(code_dir,"scRNASeq_pipeline_functions.R"))

# loading the libraries that are required throughout the analysis
library_calls()


## ----message=T,eval=T----------------------------------------------------------------------------------------------------------------------------
# Read the dataset
counts = read.delim(file.path(input_dir,"pbmc_example_counts.txt"),sep="\t")

# if we have genes that are not expressed in any cell, discard them
keep_feature = !gene_filter_by_feature_count(counts,0)
counts = counts[keep_feature,]

# make a table of metadata (e.g. batch, cell type annotation,treatment,...)
# Here, we do not have any such information, so we just give each cell a name
# Note that the rownames of this table have to correspond to the column names
# of the count matrix.
cell_annotation = data.frame(cell_idx=paste0("pbmc_C",seq(dim(counts)[2])))
rownames(cell_annotation) = colnames(counts)

# get gene annotations from ensembldb
# optional: get gene descriptions from org.Hs.eg.db (slows down the process a lot!)
# NOTE: the output of get_gene_annotations is a data.table sorted by gene identifier.
#       This means the genes are no longer in the same order as in the count matrix!
gene_annotations = get_gene_annotations(rownames(counts),organism = "human",get_descriptions = F)

# convert this to feature metadata
feature_annotation = as.data.frame(gene_annotations)
rownames(feature_annotation) = gene_annotations$gene_id
feature_annotation = feature_annotation[rownames(counts),]

#construct SingleCellExperiment object
sce =SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts=log2(as.matrix(counts)+1)),
                          colData = cell_annotation,
                          rowData = feature_annotation)



## ----eval=T--------------------------------------------------------------------------------------------------------------------------------------
#calculate QC metrics
sce = calculate_QC_metrics(sce)

# assign cell cycle phase (based on the method from scran)
# because PBMCs are post-mitotic, most cells should be assigned to G0/G1 phase
cc = annotate_cell_cycle(sce)
sce$cell_cycle_phase = cc$phases
sce$cc_scores = cc$scores

#save the SCESet
save(sce,file=file.path(out_data_dir,"sce_raw.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
load(file.path(out_data_dir,"sce_raw.RData"))

p1 = plotHighestExprs(sce,feature_names_to_plot = "symbol", exprs_values = "counts") 
print(p1)



## ------------------------------------------------------------------------------------------------------------------------------------------------
p1.2 = custom_plotHighestExprs(sce,feature_names_to_plot = "symbol",as_percentage=F,exprs_values = "counts", colour_cells_by=c("detected"))
p1.2 = p1.2 + xlab("Expression [raw counts]")
print(p1.2)

# to save a plot, use the ggsave function:
# ggsave(p1.2, file = "saved_example_plot.pdf",height=7,width=7)


## ----warning=F-----------------------------------------------------------------------------------------------------------------------------------
p2 = plotRowData(sce, x="mean", y="detected") 
p2 = p2+xlab("Mean expression [raw counts]")+ggtitle("Mean expression versus detection rate") + scale_x_log10()
print(p2)
# Check total number of zeroes
t = table(counts(sce)==0)
print(t/sum(t)*100)


## ------------------------------------------------------------------------------------------------------------------------------------------------
vars = c("sum","detected","cell_cycle_phase")
v = getVarianceExplained(sce, variables = vars)
p4 = plotExplanatoryVariables(sce, variables=vars)
print(p4 + ggtitle("Percentage of explained variance"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
min_genes = 9.0 #minimum number of features (genes) per cell [log2]
min_UMI = 10.5  #minimum total UMIs / cell [log2]
mt_threshold = 9 #Maximum percentage of mitochondrial genes


## ------------------------------------------------------------------------------------------------------------------------------------------------
plot_RNA_QC(sce, min_genes = min_genes, min_UMI = min_UMI)
plot_MT_QC(sce,mt_threshold)



## ------------------------------------------------------------------------------------------------------------------------------------------------
sce$keep_manual = (   !cell_filter_by_feature_count(counts(sce),2^min_genes) &
                      !cell_filter_by_total_UMI(counts(sce),2^min_UMI) &                                !cell_filter_by_mt_content(sce$subsets_MT_percent,mt_threshold))

table(sce$keep_manual)

sce_clean = sce[,sce$keep_manual]


## ------------------------------------------------------------------------------------------------------------------------------------------------
qc_stats = quickPerCellQC(sce, percent_subsets="subsets_MT_percent")

table(qc_stats$discard, !sce$keep_manual)


## ------------------------------------------------------------------------------------------------------------------------------------------------
n_th = 1
min_counts = 2

keep_feature = !(gene_filter_by_feature_count(counts(sce_clean),n_th, min_counts))
sce_clean = sce_clean[keep_feature,]

sce_clean = calculate_QC_metrics(sce_clean,subset = list(MT=which(rowData(sce_clean)$chr=="MT")))
save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
p = my_plot_PCA(counts = log2(counts(sce_clean)+1),
                scale_pca = T, center_pca = T, return_pca = F, use_irlba=F,
                color = colData(sce_clean)[,"sum",drop=F])
p = p+ggtitle("PCA on raw log2(counts)")
print(p)


## ------------------------------------------------------------------------------------------------------------------------------------------------
p = my_plot_tSNE(counts = log2(counts(sce_clean)+1),
                 is_distance = F, scale_pca = F, n_comp = 50, return_tsne=F,
                 color = colData(sce_clean)[,"sum",drop=F])
p = p+ggtitle("t-SNE on raw log2(counts)")
print(p)


## ----message=F-----------------------------------------------------------------------------------------------------------------------------------
# normalize data
sce_clean = normalize_counts(sce_clean,method = "scran")

# normalized values are automatically log2-transformed and
# stored in the norm_exprs slot of the SCESet
norm_exprs(sce_clean)[1:5,1:5]

save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
# PCA of normalized values
sum_norm = data.frame(sum_expression = colSums(norm_exprs(sce_clean)))
p1 = my_plot_PCA(counts = norm_exprs(sce_clean),
                 color=sum_norm,
                 return_pca = F, scale_pca = T, center_pca = T)
p1 = p1+ggtitle("PCA on normalized counts")
print(p1)

#tSNE of normalized values
p2 = my_plot_tSNE(counts = norm_exprs(sce_clean),
                  color=sum_norm,
                  return_tsne = F, is_distance = F)
p2 = p2 + ggtitle("t-SNE (50 PCs) on normalized counts")
print(p2)


## ----message=F-----------------------------------------------------------------------------------------------------------------------------------
cd20 = t(norm_exprs(sce_clean["ENSG00000156738",]))
colnames(cd20) = "CD20"

go_id = "GO:0002376" 
ens_go = GO_to_gene(go_id)
info_GO = rownames(sce_clean)%in%ens_go
table(info_GO)

p = my_plot_PCA(counts = norm_exprs(sce_clean[info_GO,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - GO:0002376 features",
                color = cd20)
print(p)


## ------------------------------------------------------------------------------------------------------------------------------------------------
info_HVG = info.genes(2^norm_exprs(sce_clean)-1,PLOT=T,qcv=0.25,pv=.1,q=.5,minBiolDisp = 0) 
table(info_HVG)
p = my_plot_PCA(counts = norm_exprs(sce_clean[info_HVG,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - HVG features",
                color = cd20)
print(p)


## ------------------------------------------------------------------------------------------------------------------------------------------------
info_NBdrop = run_DANB(counts(sce_clean),method = "NBDrop",save_plot=F, perc_genes=10) 
info_NBdisp = run_DANB(counts(sce_clean),method = "NBDisp",save_plot=F, perc_genes = 10) 
table(info_NBdrop,info_NBdisp)

p = my_plot_PCA(counts = norm_exprs(sce_clean[info_NBdrop,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - NBDrop features",
                color = cd20)
print(p)
p = my_plot_PCA(counts = norm_exprs(sce_clean[info_NBdisp,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - NBDisp features",
                color = cd20)
print(p)


## ---- eval=T-------------------------------------------------------------------------------------------------------------------------------------
sce_info = sce_clean[info_NBdrop,]
dim(sce_info)

# tSNE map of the cleaned data
# note that by setting return_tsne = T, we can obtain the t-SNE object for later use
set.seed(17) #set seed again, some packages internally set a different one
tsne_info = my_plot_tSNE(counts = norm_exprs(sce_info),
                         scale_pca = F, n_comp = 50, return_tsne=T)$tsne


## ----eval=T--------------------------------------------------------------------------------------------------------------------------------------
save(sce_info, file = file.path(out_data_dir,"sce_info.RData"))
save(tsne_info,file = file.path(out_data_dir,"tsne_info.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
#load the data we need
load(file.path(out_data_dir,"sce_info.RData"))
load(file.path(out_data_dir,"tsne_info.RData"))


## ---- fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------------------------------------------
b_cell = t(norm_exprs(sce_info["ENSG00000156738",]))
colnames(b_cell) = "B-cell"

monocyte = data.frame(Monocyte = colSums(norm_exprs(sce_info)[which(rowData(sce_info)$symbol %in% c('CD14','LYZ','FCGR3A','MS4A7')),]))

t_cell = data.frame(`T-cell` = colSums(norm_exprs(sce_info)[which(rowData(sce_info)$symbol %in% c('CD3E','CD3D','CD3G')),]))

nk_cell = data.frame(`NK cell` = colSums(norm_exprs(sce_info)[which(rowData(sce_info)$symbol %in% c('GNLY','NKG7')),]))

# Make plots
# Note that by providing the tsne input variable instead of counts,
# we can use an existing t-SNE calculation for plotting

p1 = my_plot_tSNE(tsne = tsne_info, color = b_cell, title = "B-cell marker expression")  
p2 = my_plot_tSNE(tsne = tsne_info, color = monocyte, title = "Monocyte marker expression")
p3 = my_plot_tSNE(tsne = tsne_info, color = t_cell, title = "T-cell marker expression")
p4 = my_plot_tSNE(tsne = tsne_info, color = nk_cell, title = " NK cell marker expression")
ggmultiplot(p1,p2,p3,p4,cols=2)


## ------------------------------------------------------------------------------------------------------------------------------------------------
assignment = data.table(tsne1 = tsne_info$Y[,1], tsne2 = tsne_info$Y[,2],cell_type = 'T-cell')
assignment[tsne1 < -10 ,cell_type:='Monocyte']
assignment[tsne2 < -10 ,cell_type:='B-cell']
assignment[tsne2 > 15 & tsne1 > 5,cell_type:='NK Cell']

sce_info$cell_type = assignment$cell_type

p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"cell_type",drop=F])
print(p+labs(title="t-SNE on informative genes",subtitle = "Colored by manual cell annotation"))



## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------
library(SC3)


## ----eval=T--------------------------------------------------------------------------------------------------------------------------------------
library(SC3)
#newer versions of sc3 require that there is a feature_symbol attribute for each gene. Let's quick fix and put the symbol slot there:
rowData(sce_info)$feature_symbol = rowData(sce_info)$symbol 
sce_info = sc3_prepare(sce_info, n_cores = 4)


## ----eval=T--------------------------------------------------------------------------------------------------------------------------------------
sce_info = sc3_estimate_k(sce_info)
str(metadata(sce_info)$sc3)


## ----eval=T--------------------------------------------------------------------------------------------------------------------------------------
sce_info = sc3(sce_info, ks = 5, biology = TRUE, n_cores = 4)


## ---- eval = F-----------------------------------------------------------------------------------------------------------------------------------
## sce_info = sc3(sce_info, ks = 6, biology = F, n_cores = 8)
## sce_info = sc3_run_svm(sce_info)
## 
## sce_info@sc3$svm_train_inds = NULL
## sce_info = sc3_calc_biology(sce_info, k=c(8,13), regime = "marker")
## 
## # to visualize the markers, use my modified fuction:
## 
## # change the plotted gene names to symbol for better readability
## plot_sce = sce_info
## rownames(plot_sce) = rowData(plot_sce)$symbol
## custom_sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90)
## 
## rm(plot_sce)


## ---- out.width="110%",fig.height=4--------------------------------------------------------------------------------------------------------------
sc3_plot_consensus(sce_info, k=5)


## ---- out.width="110%",fig.height=4--------------------------------------------------------------------------------------------------------------
sc3_plot_expression(sce_info, k = 5, show_pdata = c("cell_type"))


## ---- out.width="110%",fig.height=6--------------------------------------------------------------------------------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info[!is.na(rowData(sce_info)$feature_symbol),]
rownames(plot_sce) = rowData(plot_sce)$symbol

sc3_plot_markers(plot_sce, k=5, p.val = 0.1, auroc = 0.80, show_pdata = c("cell_type")) 

# for hybrid SVM approach:
# custom_sc3_plot_markers(plot_sce, k=5, p.val = 0.01, auroc = 0.90) #need to fix for new version


## ------------------------------------------------------------------------------------------------------------------------------------------------
assignment = data.table(clust = sce_info$sc3_5_clusters, cell_type = sce_info$cell_type)
assignment[clust == 1, cell_type:= 'CD4 T-cell']
assignment[clust == 2, cell_type:= 'Monocyte']
assignment[clust==3,cell_type:='B-cell']
assignment[clust==4,cell_type:='NK cell']
assignment[clust == 5, cell_type:= 'CD8 T-cell']

sce_info$SC3_assignment = assignment$cell_type
p = my_plot_tSNE(tsne = tsne_info,
                 color = colData(sce_info)[,"SC3_assignment",drop=F],
                 shape = colData(sce_info)[,"cell_type",drop=F],
                 title = "SC3 assignment",
                 show_proportions = T)
print(p)
save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
pca = my_plot_PCA(counts = norm_exprs(sce_info),return_pca=T)$pca
screeplot(pca,type = 'lines')


## ------------------------------------------------------------------------------------------------------------------------------------------------
dist_eucl = dist.gen(pca$x[,1:6],method='euclidean')
hfit = hclust(dist_eucl,method="average")
plot(hfit, labels = F, main = "hclust on euclidean distance in PCA space")


## ---- message =F---------------------------------------------------------------------------------------------------------------------------------
library(dynamicTreeCut)
groups_hclust_eucl = cutreeDynamic(hfit, distM = as.matrix(dist_eucl), deepSplit = 0, minClusterSize = 5, maxCoreScatter = 0.70, minGap = 0.25)


## ------------------------------------------------------------------------------------------------------------------------------------------------
library(cluster)
si = silhouette(groups_hclust_eucl,dist_eucl)
plot(si, col = "darkgray", border=NA, main = "Silhouette plot for hclust (euclidean in PCA space)")


## ------------------------------------------------------------------------------------------------------------------------------------------------
sce_info$hclust_sil = si[,3]
p = my_plot_tSNE(tsne = tsne_info,
                 color = colData(sce_info)[,"hclust_sil",drop=F],
                 shape = colData(sce_info)[,"cell_type",drop=F],
                 title = "tSNE colored by silhouette width")
print(p+scale_color_distiller(type="div", palette = "RdBu"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
sce_info$hclust_eucl = as.factor(groups_hclust_eucl)
table(sce_info$SC3_assignment,sce_info$hclust_eucl)


## ------------------------------------------------------------------------------------------------------------------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info
rownames(plot_sce) = rowData(plot_sce)$symbol
monocyte_markers = c('CD14','LYZ','FCGR3A','MS4A7','FCER1A','CST3')
p = plotExpression(plot_sce,features = monocyte_markers, x="hclust_eucl",
               colour_by = "hclust_eucl")
print(p)


## ----fig.width = 8-------------------------------------------------------------------------------------------------------------------------------
sce_info$hclust_eucl = factor(sce_info$hclust_eucl, levels = levels(sce_info$hclust_eucl), labels = c("CD4 T-cell","CD14++/CD16- Monocyte","B-cell","CD8 T-cell","NK cell","CD14+/CD16++ Monocyte", "Dendritic cell"))
p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"hclust_eucl",drop=F],
                 shape = colData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (euclidean) assignment",
                 show_proportions = T)
print(p)


## ------------------------------------------------------------------------------------------------------------------------------------------------
library(dendextend)
dist_pearson = dist.gen(t(norm_exprs(sce_info)),method = "pearson")
hfit = hclust(dist_pearson,method="average")
groups_hclust_pearson = cutree(hfit, k=6)
sce_info$hclust_pearson = as.factor(groups_hclust_pearson)
hfit2 = color_branches(hfit, k=6)
hfit2 = hfit2 %>% set("labels", rep("",dim(sce_info)[2]))
plot(hfit2, main = "hclust on Pearson correlation")


## ----fig.width=8---------------------------------------------------------------------------------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"hclust_pearson",drop=F],
                 shape = colData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (Pearson) assignment",
                 show_proportions = T)
print(p)


## ---- message = F--------------------------------------------------------------------------------------------------------------------------------
sim = 1-as.matrix(dist_eucl)/max(dist_eucl)
adj_corr = build_adjacency_matrix(mat = sim,cutoff="auto",is_similarity=T)


## ---- message = F--------------------------------------------------------------------------------------------------------------------------------
# use R MCL:
groups_MCL = MCLcell.clust(adj_corr, use_external_mcl = FALSE)
#use external MCL (UNIX only)
#groups_MCL = MCLcell.clust(adj_corr, use_external_mcl = TRUE, mcl_path = mcl_path)
sce_info$MCL = as.factor(groups_MCL)
table(sce_info$MCL,sce_info$hclust_eucl)


## ---- fig.width = 8------------------------------------------------------------------------------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"MCL",drop=F],
                 shape = colData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=auto) assignment",
                 show_proportions = T)
print(p)


## ----fig.height=6,fig.width = 9------------------------------------------------------------------------------------------------------------------
adj_corr = build_adjacency_matrix(mat = sim,cutoff=0.8,is_similarity=T)
groups_MCL = MCLcell.clust(adj_corr,mcl_path = mcl_path)
sce_info$MCL2 = as.factor(groups_MCL)
table(sce_info$MCL2,sce_info$hclust_eucl)

p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"MCL2",drop=F],
                 shape = colData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=0.8) assignment",
                 show_proportions = T)
p = p + scale_color_manual(values = brewer.pal(10,"Paired"))
print(p)


## ------------------------------------------------------------------------------------------------------------------------------------------------
DBSCAN_groups = run_dbscan(dist_eucl,eps="auto",min_pts=7,tol=0.02)


## ----fig.height=6,fig.width=8--------------------------------------------------------------------------------------------------------------------
sce_info$DBSCAN = as.factor(DBSCAN_groups)
p = my_plot_tSNE(tsne = tsne_info, color = colData(sce_info)[,"DBSCAN",drop=F],
                 shape = colData(sce_info)[,"hclust_eucl",drop=F],
                 title = "DBSCAN assignment",show_proportions = T)
print(p)



## ---- eval = T-----------------------------------------------------------------------------------------------------------------------------------
save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
sce_clean$SC3_assignment = sce_info$SC3_assignment
sce_clean$hclust_eucl = sce_info$hclust_eucl


## ------------------------------------------------------------------------------------------------------------------------------------------------
de_wilcox = run_wilcoxon_test(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                           fc_cutoff = 1, alpha = 0.05)
dim(de_wilcox)
head(de_wilcox)


## ------------------------------------------------------------------------------------------------------------------------------------------------
table(de_wilcox$DE_flag)


## ------------------------------------------------------------------------------------------------------------------------------------------------
mean_exprs_dt = data.table(gene_id = rownames(sce_clean),mean_exprs = rowMeans(norm_exprs(sce_clean)))
de_wilcox = merge(de_wilcox,mean_exprs_dt, by = "gene_id")
names(de_wilcox)


## ------------------------------------------------------------------------------------------------------------------------------------------------
p = generic_scatterplot(de_wilcox, x_col = "mean_exprs", y_col = "log2fc",
                        color = "DE_flag")
print(p+ggtitle('MA plot for Wilcoxon test'))


## ------------------------------------------------------------------------------------------------------------------------------------------------
de_wilcox[,log10_pval:=-log10(adj_pval)]
p = generic_scatterplot(de_wilcox, x_col = "log2fc", y_col = "log10_pval",
                        color = "DE_flag")
print(p+ggtitle('Volcano plot for Wilcoxon test'))


## ----message=F-----------------------------------------------------------------------------------------------------------------------------------
library(DT)
top_DE = de_wilcox[DE_flag==TRUE][order(log2fc,decreasing = T)]$gene_id[1:20]
gene_table = get_gene_annotations(top_DE,get_descriptions = T)
datatable(gene_table, caption = "Top 20 upregulated genes in B-cells")


## ------------------------------------------------------------------------------------------------------------------------------------------------
de_limma_trend = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "trend", fc_cutoff = 1, alpha = 0.05, count_thr = 1, pct=50)

de_limma_trend = merge(de_limma_trend, mean_exprs_dt)

p = generic_scatterplot(de_limma_trend, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag")

print(p + ggtitle("MA plot for limma-trend"))



## ---- eval = F, message = F----------------------------------------------------------------------------------------------------------------------
## library(SC3)
## 
## # add the slots the calc biology function checks
## dummy = list(`0`=0)
## sce_clean@metadata$sc3$consensus = dummy
## sce_clean@metadata$sc3$n_cores = 4
## rowData(sce_clean)$sc3_gene_filter = rep(TRUE,dim(sce_clean)[1]) #do not filter out any genes
## rowData(sce_clean)$feature_symbol = rowData(sce_clean)$symbol
## 
## # add whatever clustering assignment you want to the sc3_0_clusters slot
## sce_clean$sc3_0_clusters = as.factor(sce_info$cell_type)
## 
## sce_clean = sc3_calc_biology(sce_clean, k = 0)


## ---- out.width="110%",fig.height=6, eval=F------------------------------------------------------------------------------------------------------
## # change the plotted gene names to symbol for better readability
## plot_sce = sce_clean[!is.na(rowData(sce_clean)$symbol),]
## rownames(plot_sce) = rowData(plot_sce)$symbol
## custom_sc3_plot_markers(plot_sce, k=0, p.val = 0.01, auroc = 0.90, show_pdata='sc3_0_clusters')


## ------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()


## ------------------------------------------------------------------------------------------------------------------------------------------------
annotate_cell_cycle = function(sce, organism = "human", gene.names = rownames(sce)){
  if(organism == "human"){
    hs.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    assigned = cyclone(sce, pairs=hs.pairs, gene.names = gene.names)} else if (organism == "mouse"){
      mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
      assigned = cyclone(sce, pairs=mm.pairs, gene.names = gene.names)
    } else {stop("Organism has to be human or mouse.")}
  
  return(assigned)
}

