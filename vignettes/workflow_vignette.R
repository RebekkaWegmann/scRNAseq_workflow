## ---- eval = F-----------------------------------------------------------
## install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.5/bioc")
## library(BiocInstaller)

## ---- eval = F-----------------------------------------------------------
## install.packages("devtools")

## ---- eval=F-------------------------------------------------------------
## # from CRAN
## install.packages("Rtsne") #Rtsne v. 0.13
## install.packages("ggplot2") # ggplot2 v. 2.2.1
## install.packages("data.table") # data.table v. 1.10.4
## install.packages("RColorBrewer") # RColorBrewer v. 1.1-2
## install.packages("mvoutlier") # mvoutlier 2.0.8, required by some functions from scater
## 
## devtools::install_github("bwlewis/irlba") # irlba 2.3.2, optional. Make sure you use irlba > 2.3.1, older versions contain a bug that results in unreliable output!
## 
## # from Bioconductor
## biocLite("scater") # scater v. 1.4.0
## biocLite("scran") # scran v. 1.4.5

## ---- eval=F-------------------------------------------------------------
## # ensembldb 2.0.4 and EnsDb.Hsapiens.v79 2.1.0
## biocLite(c("ensembldb","EnsDb.Hsapiens.v79"))
## # org.Hs.eg.db v. 3.4.1
## biocLite("org.Hs.eg.db")

## ---- eval=F-------------------------------------------------------------
## # ensembldb 2.0.4 and EnsDb.Mmusculus 2.1.0
## biocLite(c("ensembldb","EnsDb.Mmusculus.v75"))
## # org.Mm.eg.db v. 3.4.1
## biocLite("org.Mm.eg.db")

## ----eval = F------------------------------------------------------------
## # M3Drop version 3.05.00 (Note: This is still under active development, please let me know if a new version breaks my functions...)
## devtools::install_github("tallulandrews/M3D")

## ---- eval = F-----------------------------------------------------------
## install.packages("cluster") # cluster v. 2.0.6
## install.packages("dendextend") # dendextend v. 1.5.2
## install.packages("Ckmeans.1d.dp") # Ckmeans.1d.dp v. 4.2.1
## install.packages("dbscan") # DBSCAN v. 1.1-1
## install.packages(c("mclust","mvtnorm")) # mclust v. 5.4 and mvtnorm 1.0.6
## install.packages("dynamicTreeCut") # dynamicTreeCut v. 1.63-1
## 
## biocLite("SC3") # SC3 v. 1.4.2
## devtools::install_github("satijalab/seurat") #Seurat 2.0.1
## bioclite("pcaMethods") # pcaMethods 1.68.0, dependency of pcaReduce
## devtools::install_github('JustinaZ/pcaReduce') #pcaReduce 1.0

## tar xzf mcl-14-137.tar.gz

## ----eval = F------------------------------------------------------------
## biocLite("limma") # limma v. 3.32.5
## biocLite("MAST") # MAST v. 1.2.1

## ---- message = FALSE----------------------------------------------------
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

#path to your MCL binary
mcl_path = "~/local/bin/mcl"

set.seed(17) #to make tSNE plots reproducible

source(file.path(code_dir,"scRNASeq_pipeline_functions.R"))

# loading the libraries that are required throughout the analysis
library_calls()

## ----message=F,eval=T----------------------------------------------------
# Read the dataset
counts = read.delim(file.path(input_dir,"pbmc_example_counts.txt"),sep="\t")

# if we have genes that are not expressed in any cell, discard them
keep_feature = !gene_filter_by_feature_count(counts,0)
counts = counts[keep_feature,]

# make a table of metadata (e.g. batch, cell type annotation,treatment,...)
# Here, we do not have any such information, so we just give each cell a name
# Note that the rownames of this table have to correspond to the column names
# of the count matrix.
annot = data.frame(cell_idx=paste0("pbmc_C",seq(dim(counts)[2])))
rownames(annot) = colnames(counts)
pd = new("AnnotatedDataFrame", data=annot)

# get gene annotations from ensembldb
# optional: get gene descriptions from org.Hs.eg.db (slows down the process a lot!)
# NOTE: the output of get_gene_annotations is a data.table sorted by gene identifier.
#       This means the genes are no longer in the same order as in the count matrix!
geneName = get_gene_annotations(rownames(counts),organism = "human",get_descriptions = F)

# convert this to feature metadata
fd_table = as.data.frame(geneName)
rownames(fd_table) = geneName$gene_id
fd_table = fd_table[rownames(counts),]
fd = new("AnnotatedDataFrame",data=fd_table)

#construct SCESet
sce = newSCESet(countData = counts, phenoData=pd,featureData = fd)


## ----eval=T--------------------------------------------------------------
#calculate QC metrics
sce = calculateQCMetrics(sce,feature_controls = list(MT=which(fData(sce)$chr=="MT")))

# assign cell cycle phase (based on the method from scran)
# because PBMCs are post-mitotic, most cells should be assigned to G0/G1 phase
cc = annotate_cell_cycle(sce)
sce$cell_cycle_phase = cc$phases
sce$cc_scores = cc$scores

#save the SCESet
save(sce,file=file.path(out_data_dir,"sce_raw.RData"))

## ------------------------------------------------------------------------
load(file.path(out_data_dir,"sce_raw.RData"))

p1 = plotQC(sce, type="high",feature_names_to_plot = "symbol") 
print(p1)


## ------------------------------------------------------------------------
p1.2 = custom_plotHighestExprs(sce,feature_names_to_plot = "symbol")
p1.2 = p1.2 + xlab("Expression [raw counts, log2]")
print(p1.2)

# to save a plot, use the ggsave function:
ggsave(p1.2, file = "saved_example_plot.pdf",height=7,width=7)

## ----warning=F-----------------------------------------------------------
p2 = plotQC(sce, type = "exprs-freq-vs-mean") 
p2 = p2+xlab("Mean expression [raw counts, log2]")+ggtitle("Mean expression versus detection rate")
print(p2)
# Check total number of zeroes
t = table(counts(sce)==0)
print(t/sum(t)*100)

## ------------------------------------------------------------------------
p3 = plotQC(sce, type="find", variable="total_features") 
print(p3 + ggtitle("Correlation of principal components with total detected features."))

## ------------------------------------------------------------------------
p3.2 = plotQC(sce, type="find", variable="cell_cycle_phase")
print(p3.2+ggtitle("Correlation of principal components with cell cycle phase."))

## ------------------------------------------------------------------------
vars = c("total_counts","total_features","cell_cycle_phase")
p4 = plotQC(sce, type="expl", variables=vars)
print(p4 + ggtitle("Percentage of explained variance"))

## ------------------------------------------------------------------------
min_genes = 9.0 #minimum number of features (genes) per cell [log2]
min_UMI = 10.5  #minimum total UMIs / cell [log2]
mt_threshold = 9 #Maximum percentage of mitochondrial genes

## ------------------------------------------------------------------------
plot_RNA_QC(sce, min_genes = min_genes, min_UMI = min_UMI)
plot_MT_QC(sce,mt_threshold)


## ------------------------------------------------------------------------
sce$keep_manual = (   !cell_filter_by_feature_count(counts(sce),2^min_genes) &
                      !cell_filter_by_total_UMI(counts(sce),2^min_UMI) &                                !cell_filter_by_mt_content(sce$pct_counts_feature_controls_MT,mt_threshold))

table(sce$keep_manual)

sce_clean = sce[,sce$keep_manual]

## ------------------------------------------------------------------------
n_th = 1
min_counts = 2

keep_feature = !(gene_filter_by_feature_count(counts(sce_clean),n_th, min_counts))
sce_clean = sce_clean[keep_feature,]

## ------------------------------------------------------------------------

sce_clean = calculateQCMetrics(sce_clean,feature_controls = list(MT=which(fData(sce_clean)$chr=="MT")))

# The variables used to detect outliers
vars = c( "pct_counts_top_100_features", 
          "total_features", "pct_counts_feature_controls", 
           "log10_counts_endogenous_features", 
           "log10_counts_feature_controls")

sce_clean = plotPCA(sce_clean,
                    size_by = "total_features", 
                    pca_data_input = "pdata",
                    selected_variables = vars,
                    detect_outliers = TRUE,
                    return_SCESet = TRUE)

table(sce_clean$outlier)

#here, we remove the outliers
sce_clean = sce_clean[,!sce_clean$outlier]

# Finally, we again remove genes that are not expressed
keep_feature = !(gene_filter_by_feature_count(counts(sce_clean),n_th,min_counts))
sce_clean = sce_clean[keep_feature,]

save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))

## ------------------------------------------------------------------------
p = my_plot_PCA(counts = log2(counts(sce_clean)+1),
                scale_pca = T, center_pca = T, return_pca = F, use_irlba=F,
                color = pData(sce_clean)[,"total_features",drop=F])
p = p+ggtitle("PCA on raw log2(counts)")
print(p)

## ------------------------------------------------------------------------
p = my_plot_tSNE(counts = log2(counts(sce_clean)+1),
                 is_distance = F, scale_pca = F, n_comp = 50, return_tsne=F,
                 color = pData(sce_clean)[,"total_features",drop=F])
p = p+ggtitle("t-SNE on raw log2(counts)")
print(p)

## ----message=F-----------------------------------------------------------
# normalize data
sce_clean = normalize_counts(sce_clean,method = "scran")

# normalized values are automatically log2-transformed and
# stored in the norm_exprs slot of the SCESet
norm_exprs(sce_clean)[1:5,1:5]

save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))

## ------------------------------------------------------------------------
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

## ----message=F-----------------------------------------------------------
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

## ------------------------------------------------------------------------
info_HVG = info.genes(2^norm_exprs(sce_clean)-1,PLOT=T,qcv=0.25,pv=.1,q=.5,minBiolDisp = 0) 
table(info_HVG)
p = my_plot_PCA(counts = norm_exprs(sce_clean[info_HVG,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - HVG features",
                color = cd20)
print(p)

## ------------------------------------------------------------------------
info_NBdrop = run_DANB(counts(sce_clean),method = "NBDrop",save_plot=F, cutoff = 0.1) 
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

## ---- eval=T-------------------------------------------------------------
sce_info = sce_clean[info_NBdrop,]
dim(sce_info)

# tSNE map of the cleaned data
# note that by setting return_tsne = T, we can obtain the t-SNE object for later use
tsne_info = my_plot_tSNE(counts = norm_exprs(sce_info),
                         scale_pca = F, n_comp = 50, return_tsne=T)$tsne

## ----eval=T--------------------------------------------------------------
save(sce_info, file = file.path(out_data_dir,"sce_info.RData"))
save(tsne_info,file = file.path(out_data_dir,"tsne_info.RData"))

## ------------------------------------------------------------------------
#load the data we need
load(file.path(out_data_dir,"sce_info.RData"))
load(file.path(out_data_dir,"tsne_info.RData"))

## ---- fig.align='center', fig.width=12, fig.height=8---------------------
b_cell = t(norm_exprs(sce_info["ENSG00000156738",]))
colnames(b_cell) = "B-cell"

monocyte = data.frame(Monocyte = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('CD14','LYZ','FCGR3A','MS4A7')),]))

t_cell = data.frame(`T-cell` = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('CD3E','CD3D','CD3G')),]))

nk_cell = data.frame(`NK cell` = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('GNLY','NKG7')),]))

# Make plots
# Note that by providing the tsne input variable instead of counts,
# we can use an existing t-SNE calculation for plotting

p1 = my_plot_tSNE(tsne = tsne_info, color = b_cell, title = "B-cell marker expression")  
p2 = my_plot_tSNE(tsne = tsne_info, color = monocyte, title = "Monocyte marker expression")
p3 = my_plot_tSNE(tsne = tsne_info, color = t_cell, title = "T-cell marker expression")
p4 = my_plot_tSNE(tsne = tsne_info, color = nk_cell, title = " NK cell marker expression")
ggmultiplot(p1,p2,p3,p4,cols=2)

## ------------------------------------------------------------------------
assignment = data.table(tsne1 = tsne_info$Y[,1], tsne2 = tsne_info$Y[,2],cell_type = 'T-cell')
assignment[tsne1 < -10 ,cell_type:='B-cell']
assignment[tsne1 > 5 ,cell_type:='Monocyte']
assignment[tsne2 < -17 & tsne1 > -1,cell_type:='NK Cell']

sce_info$cell_type = assignment$cell_type

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"cell_type",drop=F])
print(p+labs(title="t-SNE on informative genes",subtitle = "Colored by manual cell annotation"))


## ---- echo=F-------------------------------------------------------------
library(SC3)

## ----eval=T--------------------------------------------------------------
library(SC3)
sce_info = sc3_prepare(sce_info, ks = 2:10, n_cores = 4)

## ----eval=T--------------------------------------------------------------
sce_info = sc3_estimate_k(sce_info)
sce_info@sc3$k_estimation

## ----eval=T--------------------------------------------------------------
sce_info = sc3(sce_info, ks = 6, biology = TRUE, n_cores = 4)

## ---- eval = F-----------------------------------------------------------
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
## rownames(plot_sce) = fData(plot_sce)$symbol
## custom_sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90)
## 
## rm(plot_sce)

## ---- out.width="110%",fig.height=4--------------------------------------
sc3_plot_consensus(sce_info, k=6)

## ---- out.width="110%",fig.height=4--------------------------------------
sc3_plot_expression(sce_info, k = 6, show_pdata = c("cell_type"))

## ---- out.width="110%",fig.height=6--------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info
rownames(plot_sce) = fData(plot_sce)$symbol
sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90, show_pdata = c("cell_type"))

# for hybrid SVM approach:
# custom_sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90)

## ------------------------------------------------------------------------
assignment = data.table(clust = sce_info$sc3_6_clusters, cell_type = sce_info$cell_type)
assignment[clust == 1, cell_type:= 'T helper cell']
assignment[clust == 2, cell_type:= 'Monocyte']
assignment[clust==3,cell_type:='?']
assignment[clust==4,cell_type:='CD8+ T-cell']
assignment[clust==5,cell_type:='NK cell']
assignment[clust == 6, cell_type:= 'B-cell']

sce_info$SC3_assignment = assignment$cell_type
p = my_plot_tSNE(tsne = tsne_info,
                 color = pData(sce_info)[,"SC3_assignment",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "SC3 assignment",
                 show_proportions = T)
print(p)
save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))

## ------------------------------------------------------------------------
pca = my_plot_PCA(counts = norm_exprs(sce_info),return_pca=T)$pca
screeplot(pca,type = 'lines')

## ------------------------------------------------------------------------
dist_eucl = dist.gen(pca$x[,1:6],method='euclidean')
hfit = hclust(dist_eucl,method="average")
plot(hfit, labels = F, main = "hclust on euclidean distance in PCA space")

## ---- message =F---------------------------------------------------------
library(dynamicTreeCut)
groups_hclust_eucl = cutreeDynamic(hfit, distM = as.matrix(dist_eucl), deepSplit = 0, minClusterSize = 5, maxCoreScatter = 0.70, minGap = 0.25)

## ------------------------------------------------------------------------
library(cluster)
si = silhouette(groups_hclust_eucl,dist_eucl)
plot(si, col = "darkgray", border=NA, main = "Silhouette plot for hclust (euclidean in PCA space)")

## ------------------------------------------------------------------------
sce_info$hclust_sil = si[,3]
p = my_plot_tSNE(tsne = tsne_info,
                 color = pData(sce_info)[,"hclust_sil",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "tSNE colored by silhouette width")
print(p+scale_color_distiller(type="div", palette = "RdBu"))

## ------------------------------------------------------------------------
sce_info$hclust_eucl = as.factor(groups_hclust_eucl)
table(sce_info$SC3_assignment,sce_info$hclust_eucl)

## ------------------------------------------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info
rownames(plot_sce) = fData(plot_sce)$symbol
monocyte_markers = which(fData(sce_info)$symbol %in% c('CD14','LYZ','FCGR3A','MS4A7'))
p = plotExpression(plot_sce,features = monocyte_markers, x="hclust_eucl",
               colour_by = "hclust_eucl")
print(p)

## ----fig.width = 8-------------------------------------------------------
sce_info$hclust_eucl = factor(sce_info$hclust_eucl, levels = levels(sce_info$hclust_eucl), labels = c("T-helper cell","CD14++/CD16- Monocyte","B-cell","CD8+ T-cell","NK cell","CD14+/CD16++ Monocyte"))
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"hclust_eucl",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (euclidean) assignment",
                 show_proportions = T)
print(p)

## ------------------------------------------------------------------------
library(dendextend)
dist_pearson = dist.gen(t(norm_exprs(sce_info)),method = "pearson")
hfit = hclust(dist_pearson,method="average")
groups_hclust_pearson = cutree(hfit, k=6)
sce_info$hclust_pearson = as.factor(groups_hclust_pearson)
hfit2 = color_branches(hfit, k=6)
hfit2 = hfit2 %>% set("labels", rep("",dim(sce_info)[2]))
plot(hfit2, main = "hclust on Pearson correlation")

## ----fig.width=8---------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"hclust_pearson",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (Pearson) assignment",
                 show_proportions = T)
print(p)

## ------------------------------------------------------------------------
library(pcaReduce)
expr_mat_info = t(norm_exprs(sce_info))
assignment_info = PCAreduce(expr_mat_info, nbt = 5, q = 6, method = "M")

## ----fig.width=12,fig.height=6-------------------------------------------
table(assignment_info[[1]][,4])

plots = list()
for(i in 1:5){
  df = data.frame(apply(assignment_info[[i]],c(1,2),as.factor))
  p = my_plot_tSNE(tsne=tsne_info, color = df[,'cl_id.2',drop=F],
                   shape = pData(sce_info)[,"cell_type",drop=F],
                   title = paste0("pcaReduce, run ", i))
  plots[[i]] = p
}

ggmultiplot(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],cols=3)


## ---- message = F--------------------------------------------------------
sim = 1-as.matrix(dist_eucl)/max(dist_eucl)
adj_corr = build_adjacency_matrix(mat = sim,cutoff="auto",is_similarity=T)

## ---- message = F--------------------------------------------------------
groups_MCL = MCLcell.clust(adj_corr,mcl_path = mcl_path)
sce_info$MCL = as.factor(groups_MCL)
table(sce_info$MCL,sce_info$hclust_eucl)

## ---- fig.width = 8------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"MCL",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=auto) assignment",
                 show_proportions = T)
print(p)

## ----fig.height=6,fig.width = 9------------------------------------------
adj_corr = build_adjacency_matrix(mat = sim,cutoff=0.8,is_similarity=T)
groups_MCL = MCLcell.clust(adj_corr,mcl_path = mcl_path)
sce_info$MCL2 = as.factor(groups_MCL)
table(sce_info$MCL2,sce_info$hclust_eucl)

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"MCL2",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=0.7) assignment",
                 show_proportions = T)
p = p + scale_color_manual(values = brewer.pal(10,"Paired"))
print(p)

## ----eval = T------------------------------------------------------------
library(Seurat)
seurat_assignments = seurat_clustering(sce_info,vars.to.regress=NULL,res=1.2)
sce_info$seurat = as.factor(seurat_assignments)

## ----fig.width=8, eval = T-----------------------------------------------
table(sce_info$seurat,sce_info$hclust_eucl)
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"seurat",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "Seurat assignment",show_proportions = T)
print(p)

## ------------------------------------------------------------------------
DBSCAN_groups = run_dbscan(dist_eucl,eps="auto",min_pts=7,tol=0.005)

## ----fig.height=6,fig.width=8--------------------------------------------
sce_info$DBSCAN = as.factor(DBSCAN_groups)
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"DBSCAN",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "DBSCAN assignment",show_proportions = T)
print(p)


## ----fig.width=7,fig.height=6--------------------------------------------
si = silhouette(DBSCAN_groups, dist_eucl)
plot(si, col = "darkgray", border=NA, main = "Silhouette plot for dbscan (euclidean in PCA space)")

## ------------------------------------------------------------------------
boot = fpc::clusterboot(data = dist_eucl, clustermethod = fpc::dbscanCBI, eps = 5.53,
                         MinPts = 7,B = 100, distances = T, count=F, method = "dist",
                        bootmethod = "boot")

dt = melt(data.table(cluster = c(0:4), boot$bootresult),id.vars="cluster",val="jaccard index",var="boot number")
dt = dt[cluster!=0]
p = ggplot(dt, aes(x=as.factor(cluster),y=`jaccard index`)) + geom_boxplot() +theme_bw()
print(p + ggtitle("Jaccard similarity between cluster assignment on the full data and on 100 bootstrap samples"))

## ------------------------------------------------------------------------
library(mclust)
mclust_assignment = gmm_main(pc_scores = pca$x,n_comp=6, k=1:8,n_cores=1)

## ------------------------------------------------------------------------
sce_info$mclust = as.factor(mclust_assignment)
table(sce_info$mclust,sce_info$hclust_eucl)

## ---- fig.width = 8------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"mclust",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "Mclust assignment",show_proportions = T)
print(p)

## ---- fig.width=8--------------------------------------------------------
mclust_model = gmm_main(pc_scores = pca$x,n_comp=6,do_crossval = F, best_k = 6, return_model = TRUE)
sce_info$mclust_forced = as.factor(mclust_model$classification)

table(sce_info$mclust_forced,sce_info$hclust_eucl)

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"mclust_forced",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "Mclust assignment")
print(p)

## ----message=F,results="hide"--------------------------------------------
mclust_boot = MclustBootstrap(mclust_model, nboot=500)

## ------------------------------------------------------------------------
summary(mclust_boot, what = "ci")

## ----fig.width=8, fig.height = 4-----------------------------------------
par(mfrow=c(1,6))
plot(mclust_boot, what = "pro")

## ---- fig.width=8,fig.height=12------------------------------------------
par(mfrow=c(6,6))
plot(mclust_boot, what = "mean")
par(mfrow=c(6,6))
plot(mclust_boot, what = "var")
par(mfrow=c(1,1))

## ---- eval = F-----------------------------------------------------------
## save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))

## ---- warning=F----------------------------------------------------------
sce_clean$cell_type = sce_info$cell_type
cellsius_out = cellsius_main(sce_clean, group_id = "cell_type", min_n_cells = 5,
                                 verbose =T, min_fc = 2, organism = "human", iter=0, 
                                 max_perc_cells = 50, fc_between_cutoff = 1, mcl_path = mcl_path)

## ------------------------------------------------------------------------
cellsius_out

## ------------------------------------------------------------------------
unique(cellsius_out[cell_idx=="AAAGCAAGTCAAGCGA"]$sub_cluster)

## ------------------------------------------------------------------------
unique(cellsius_out[sub_cluster == "T-cell_1_1"][,c('gene_id','symbol','description')])

## ------------------------------------------------------------------------
cellsius_out[,length(unique(cell_idx)),by="sub_cluster"]

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out, tsne=tsne_info)

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out[main_cluster=="T-cell"], tsne=tsne_info)

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out[sub_cluster=="T-cell_1_1"], tsne=tsne_info)

## ---- fig.width=8,fig.height=6-------------------------------------------
marker_idx = which(rownames(sce_clean)%in%cellsius_out[sub_cluster=='T-cell_1_1']$gene_id)

plotlist = list()
colors = t(norm_exprs(sce_clean)[marker_idx,,drop=F])
colnames(colors) =  fData(sce_clean)$symbol[marker_idx]

for(i in seq(dim(colors)[2])){
  plotlist[[i]] = my_plot_tSNE(tsne = tsne_info,
                      color = colors[,i,drop=F],
                      alpha = 0.8, title = colnames(colors)[i]) + scale_color_distiller(palette = 'RdYlBu')
}

ggmultiplot(plotlist[[1]],plotlist[[2]],plotlist[[3]], plotlist[[4]],cols = 2)

## ---- eval = F-----------------------------------------------------------
## library(shiny)
## launch_marker_vis_app(tsne = tsne_info, sce=sce_clean, marker_idx = marker_idx)

## ------------------------------------------------------------------------
final_clusters = cellsius_final_cluster_assignment(cellsius_out, sce_clean, group_id = "cell_type", min_n_genes=3)

table(final_clusters$cluster)

p = my_plot_tSNE(tsne = tsne_info, color = final_clusters, title = "CellSIUS cluster assignment")
print(p)

## ------------------------------------------------------------------------
sce_clean$SC3_assignment = sce_info$SC3_assignment
sce_clean$hclust_eucl = sce_info$hclust_eucl

## ------------------------------------------------------------------------
de_wilcox = run_wilcoxon_test(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                           fc_cutoff = 1, alpha = 0.05)
dim(de_wilcox)
head(de_wilcox)

## ------------------------------------------------------------------------
table(de_wilcox$DE_flag)

## ------------------------------------------------------------------------
mean_exprs_dt = data.table(gene_id = rownames(sce_clean),mean_exprs = fData(sce_clean)$mean_exprs)
de_wilcox = merge(de_wilcox,mean_exprs_dt, by = "gene_id")
names(de_wilcox)

## ------------------------------------------------------------------------
p = generic_scatterplot(de_wilcox, x_col = "mean_exprs", y_col = "log2fc",
                        color = "DE_flag")
print(p+ggtitle('MA plot for Wilcoxon test'))

## ------------------------------------------------------------------------
de_wilcox[,log10_pval:=-log10(adj_pval)]
p = generic_scatterplot(de_wilcox, x_col = "log2fc", y_col = "log10_pval",
                        color = "DE_flag")
print(p+ggtitle('Volcano plot for Wilcoxon test'))

## ----message=F-----------------------------------------------------------
library(DT)
top_DE = de_wilcox[DE_flag==TRUE][order(log2fc,decreasing = T)]$gene_id[1:20]
gene_table = get_gene_annotations(top_DE,get_descriptions = T)
datatable(gene_table, caption = "Top 20 upregulated genes in B-cells")

## ------------------------------------------------------------------------
de_limma_voom = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "voom", fc_cutoff = 1, alpha = 0.05, count_thr = 1,pct=50)

## ------------------------------------------------------------------------
dim(de_limma_voom)
names(de_limma_voom)

## ------------------------------------------------------------------------
table(de_limma_voom$DE_flag)
table(de_wilcox[de_limma_voom$gene_id,on="gene_id"]$DE_flag,de_limma_voom$DE_flag)

## ------------------------------------------------------------------------
de_limma_voom = merge(de_limma_voom,mean_exprs_dt, by = "gene_id")
names(de_limma_voom)

p = generic_scatterplot(de_limma_voom, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag")
print(p+ggtitle('MA plot for limma-voom'))

## ------------------------------------------------------------------------

voom_no_filt = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "voom", fc_cutoff = 1, alpha = 0.05, count_thr = 0)
voom_no_filt = merge(voom_no_filt, mean_exprs_dt)
p = generic_scatterplot(voom_no_filt, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag")

print(p+ggtitle("MA plot for limma-voom, no filtering"))

table(voom_no_filt[DE_flag==TRUE]$gene_id%in%de_limma_voom[DE_flag==TRUE]$gene_id)

## ------------------------------------------------------------------------
example_dt = copy(voom_no_filt)
example_dt[,logFC:=logFC-1*1/(mean_exprs+1)]
example_dt[,DE_flag := adj.P.Val < 0.05 & abs(logFC)>1]
p = generic_scatterplot(example_dt, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag") + geom_hline(yintercept = 0)
print(p+ggtitle("An example of an MA plot for limma gone wrong"))


## ------------------------------------------------------------------------
de_limma_trend = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "trend", fc_cutoff = 1, alpha = 0.05, count_thr = 1, pct=50)

de_limma_trend = merge(de_limma_trend, mean_exprs_dt)
de_limma_trend[,voom_overlap:=de_limma_trend$DE_flag == de_limma_voom$DE_flag]

p = generic_scatterplot(de_limma_trend, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag", shape = "voom_overlap")

print(p + ggtitle("MA plot for limma-trend"))

table(de_limma_trend$DE_flag,de_limma_voom$DE_flag)
table(de_wilcox[de_limma_trend$gene_id,on="gene_id"]$DE_flag,de_limma_trend$DE_flag)

# Compare p-values between Wilcoxon test and limma
ggplot(data.table(wilcox = de_wilcox[de_limma_trend$gene_id,on="gene_id"]$pval,
                  limmatrend = de_limma_trend$P.Val),
       aes(x=-log10(wilcox),y=-log10(limmatrend)))+geom_point()+theme_bw()+
  ggtitle('Comparison of p-values between limmatrend and Wilcoxon test')

## ---- message = F--------------------------------------------------------
de_MAST = run_MAST(sce_clean,cl_id = "hclust_eucl",cl_ref = "B-cell",norm=T,
                   set_thresh=T,fc_cutoff = 1, alpha=0.05)

## ---- message = F--------------------------------------------------------
de_MAST = run_MAST(sce_clean,cl_id = "hclust_eucl",cl_ref = "B-cell",norm=T,
                   set_thresh=F,fc_cutoff = 1, alpha=0.5)

## ----fig.width=10--------------------------------------------------------
de_MAST = merge(de_MAST, mean_exprs_dt)

p = generic_scatterplot(de_MAST, x_col = "mean_exprs", y_col = "log2FC",
                        color = "DE_flag")
print(p+ggtitle("MA plot for MAST"))

table(de_MAST$DE_flag)
table(de_MAST[de_limma_trend$gene_id,]$DE_flag,de_limma_trend$DE_flag)

# Compare p-values between MAST and limma
ggplot(data.table(MAST = de_MAST[de_limma_trend$gene_id,on="gene_id"]$pval,
                  limmatrend = de_limma_trend$P.Val),
       aes(x=-log10(MAST),y=-log10(limmatrend)))+geom_point()+theme_bw()+
  ggtitle('Comparison of p-values between limmatrend and MAST')

## ----fig.width=10, eval = T----------------------------------------------
merged = merge(de_MAST,
               de_limma_trend[, c("gene_id","DE_flag"),with=F],
               by = "gene_id",all=T)
merged = merge(merged, de_wilcox[,c("gene_id","DE_flag"),with=F],by="gene_id",all=T)

setnames(merged,which(grepl("DE_flag",names(merged))),paste0("DE_flag_",c(1:3)))
merged[,overlap:=factor(DE_flag_1==T & DE_flag_2==T & DE_flag_3==T, labels =c("Not DE or not found by all methods","Overlap of all methods"))]

p = generic_scatterplot(merged, x_col = "mean_exprs", y_col = "log2FC", color = "overlap" )
print(p+ggtitle("Overlap between all methods tested") + ylab("log2FC (MAST)"))

## ---- eval = T, message = F----------------------------------------------
library(SC3)

# add the slots the calc biology function checks
dummy = list(`0`=0)
sce_clean@sc3$consensus = dummy
sce_clean@sc3$n_cores = 4
fData(sce_clean)$sc3_gene_filter = rep(TRUE,dim(sce_clean)[1]) #do not filter out any genes

# add whatever clustering assignment you want to the sc3_0_clusters slot
sce_clean$sc3_0_clusters = as.factor(sce_info$cell_type) 

sce_clean = sc3_calc_biology(sce_clean, k = 0)

## ---- out.width="110%",fig.height=6--------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_clean[!is.na(fData(sce_clean)$symbol),]
rownames(plot_sce) = fData(plot_sce)$symbol
custom_sc3_plot_markers(plot_sce, k=0, p.val = 0.01, auroc = 0.90, show_pdata='sc3_0_clusters')

## ------------------------------------------------------------------------
sessionInfo()

## ------------------------------------------------------------------------
library_calls = function(){
  library(Rtsne) 
  library(ggplot2)
  library(data.table)
  library(scater)
  library(scran)
  library(RColorBrewer)
}

## ------------------------------------------------------------------------
get_gene_annotations = function(gene_list,v=F,get_descriptions=T,organism = "human"){
  library(ensembldb)
  if(organism == "human"){
    library(EnsDb.Hsapiens.v79)
    edb = EnsDb.Hsapiens.v79
  } else if(organism == "mouse"){
    library(EnsDb.Mmusculus.v75)
    edb = EnsDb.Mmusculus.v75
  } else {stop("Currently, annotation is only available for organisms human or mouse.")}
  
  my_get = function(x,db,v){
    out = tryCatch({get(x,db)[1]},
                   error = function(cond){
                     if(v) {message(cond)}
                     return("")},
                   finally = {})
    return(out)
  }
  
  gene_info = ensembldb::genes(edb,filter = list(GeneIdFilter(gene_list)))
  gene_info_dt = data.table(gene_id = names(gene_info),
                            chr = as.character(seqnames(gene_info)),
                            symbol = make.unique(gene_info$symbol),
                            gene_biotype = gene_info$gene_biotype)
  
  geneName = data.table(gene_id = gene_list)
  geneName = merge(geneName,gene_info_dt,by="gene_id",all=T)
  
  if(get_descriptions & organism == "human"){
    library(org.Hs.eg.db)
    geneName[,eg:=my_get(gene_id,db=org.Hs.egENSEMBL2EG,v),by="gene_id"]
    geneName[,description:=my_get(eg,db=org.Hs.egGENENAME,v),by="eg"]
  } else if(get_descriptions){
    library(org.Mm.eg.db)
    geneName[,eg:=my_get(gene_id,db=org.Mm.egENSEMBL2EG,v),by="gene_id"]
    geneName[,description:=my_get(eg,db=org.Mm.egGENENAME,v),by="eg"]
  }
  return(geneName)
}

## ------------------------------------------------------------------------

get_gene_annotations_biomart = function(gene_list){
  library(biomaRt)
  # Load the organism-specific biomart
  ensembl  =  biomaRt::useEnsembl(
    biomart = 'ensembl', 
    dataset = paste0('hsapiens', '_gene_ensembl'),
    version = 83
  )
  #
  geneName  =  biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 
                                            "entrezgene", "description","chromosome_name"), 
                             filters = 'ensembl_gene_id',
                             values = gene_list, 
                             mart = ensembl)
  
  description = lapply(seq(length(geneName$description)),function(i){
    strsplit(geneName$description[i],"[[]")[[1]][1]
  })
  description = (unlist(description))
  geneName$description = description
  colnames(geneName) = c('gene_id','symbol','eg','description','chr')
  geneName = data.table(geneName)
  setkey(geneName,'gene_id')
  geneName = unique(geneName) #remove duplicate entrez gene identifiers
  geneName[,symbol:=make.unique(symbol)]
  #save(geneName,file="data/output/geneName.RData")
  return(geneName)
}

## ------------------------------------------------------------------------
#_____________
# Gene filters
#___________

# Filtering out genes that are not expressed at a minimum of min_counts in at least n_th cells
gene_filter_by_feature_count = function(counts, n_th , min_counts = 1){
  discard = rowSums(counts>=min_counts) <= n_th
  message(paste('Flagged', length(which(discard)), 'genes.'))
  return(discard)
}

#______________
# Cell filters
#______________

# Filtering out cells with fewer than n_th detected features
cell_filter_by_feature_count = function(counts, n_th){
  discard = colSums(counts>0) <= n_th
  message(paste('Flagged', length(which(discard)), 'cells.'))
  return(discard)
}

# Filtering out cells with fewer than n_th total UMI counts
cell_filter_by_total_UMI = function(counts, n_th){
  discard = colSums(counts) <= n_th
  message(paste('Flagged', length(which(discard)), 'cells.'))
  return(discard)
}

# Filtering out cells with high mitochondrial gene content
calc_mt_content = function(counts, geneName){
  mt = rownames(counts[rownames(counts) %in% rownames(geneName[geneName$chromosome_name=="MT",]),])
  mt_not = setdiff(rownames(counts),rownames(counts[mt,]))
  counts_mt = counts[mt,]
  mt_genes.amount = 100/colSums(counts)*colSums(counts[mt,])
  return(mt_genes.amount)
}

cell_filter_by_mt_content = function(mt_genes.amount, t){
  discard = mt_genes.amount > t
  message(paste('Flagged', length(which(discard)),'cells.'))
  return(discard)
}

## ------------------------------------------------------------------------
# Plotting QC based on RNA amount detected per cell
plot_RNA_QC = function(input_sce, min_genes, min_UMI){
  par(mfrow=c(1,3))
  hist(log2(input_sce$total_features),xlab="log2[ # detected genes per cell]", main='', cex.axis=1.5,n=100)
  abline(v=min_genes,col=2)
  hist(log2(input_sce$total_counts),xlab="log2 [# of UMIs per cell]", main='', cex.axis=1.5,n=100)
  abline(v=min_UMI,col=2)
  plot(log2(input_sce$total_features),log2(input_sce$total_counts),xlab="log2[ # detected features per cell]",ylab="log2 [total counts per cell]", cex.axis=1.5)
  abline(v=min_genes,col=2)
  abline(h=min_UMI,col=2)
}

# Plotting mitochondrial gene QC
plot_MT_QC = function(sce, t){
  par(mfrow=c(1,1))
  mt_genes.amount = sce$pct_counts_feature_controls_MT 
  #mt gene content per cell
  plot(seq(length(mt_genes.amount)),(mt_genes.amount),pch=10,cex=1.5,col="gray",main=""
       , cex.axis=2,cex.lab=1.5,xlab="cell index",ylab="Ratio of MT-genes[%]")
  points(seq(length(mt_genes.amount))[mt_genes.amount<t],mt_genes.amount[mt_genes.amount<t],pch=10,cex=1.5)
  abline(h=t,col="red",lty=2,cex=5)
  
  #UMI vs no. genes colored by mt gene content
  plotPhenoData(sce, aes_string(x = "log2(total_features)",
                                y = "log2(total_counts)",
                                colour = "pct_counts_feature_controls_MT"))+
    xlab("Total detected features [log2]") + ylab("Total counts [log2]")+
    ggtitle("Total features vs. total counts, colored by MT content")
}


## ------------------------------------------------------------------------
annotate_cell_cycle = function(sce, organism = "human", gene.names = rownames(sce)){
  if(organism == "human"){
    hs.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    assigned = cyclone(sce, pairs=hs.pairs, gene.names = gene.names)} else if (organism == "mouse"){
      mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
      assigned = cyclone(sce, pairs=mm.pairs, gene.names = gene.names)
    } else {stop("Organism has to be human or mouse.")}
  
  return(assigned)
}

## ------------------------------------------------------------------------
normalize_counts = function(sce,method = "scran"){
  #calculate size factors according to method
  switch(method,
         "TC" ={sce = normalizeExprs(sce, method="none",return_norm_as_exprs=T)},
         "RLE" = {sce = normalizeExprs(sce, method="RLE",return_norm_as_exprs=T)},
         "TMM" = {sce = normalizeExprs(sce, method="TMM",return_norm_as_exprs=T)},
         "UQ" = {sce = normalizeExprs(sce, method="upperquartile",return_norm_as_exprs=T)}, 
         "scran" = {clusters = quickCluster(sce)
         sizes = seq(20,100,5)
         if(min(table(clusters))>max(sizes)){
           sce = computeSumFactors(sce,clusters = clusters,sizes=sizes)
         } else{
           message("Clustering of cells failed, using global scaling factors")
           sce = computeSumFactors(sce)
           if(any(sizeFactors(sce) < 0)) {
             warning("Negative size factors generated. Most likely, this is due to some cells having very low total feature counts. Consider using more stringent QC cutoffs.")
           }
         }
         sce = scater::normalize(sce, return_norm_as_exprs=T)}
  )
  return(sce)  
}

## ------------------------------------------------------------------------
info.genes = function(x,PLOT=F,qcv=.3,pv=.05,q=.95,minBiolDisp=0.1, perc_genes = NULL){
  
  if(!(is.null(perc_genes)|is.null(pv))){
    stop("Please provide either pv or perc_genes, not both!")
  }
  library(statmod)
  
  # calculate mean, variance and CV
  means  = rowMeans(x)
  vars  =  (apply(x,1,var))
  cv2 = vars/means^2
  
  # exclude genes with very low mean
  minMeanForFit  =  unname( quantile( means[ which( cv2 > qcv) ], q ) )#is the 95% of the means with a dispersion greter than 0.3
  useForFit  =  means >= minMeanForFit
  
  #fit model
  fit  =  glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] ) # linear fit
  fit$coefficients
  
  a0 = unname(fit$coefficients["a0"])
  a1 = unname(fit$coefficients["a1tilde"]) #we assume xi = mean(technical size factors) = 0
  minBiolDisp = minBiolDisp^2 #squared minimum biological variance
  psia1theta = a1 # this is the term psi+a1*theta that appears in the formula for omega
  # again, we assume that the technical sf = 0 and mean ratio of all size factors = 1
  m  =  ncol(x)
  cv2th  =  a0 + minBiolDisp + a0 * minBiolDisp #a0 adjusted for min biol variation
  testDenom  =  ( means * psia1theta + means^2 * cv2th ) / ( 1 + cv2th/m ) #omega
  
  pval  =  pchisq(vars*(m-1)/testDenom,df=m-1,lower.tail=F) #Chi^2 distribution
  adj.pval  =  sort(p.adjust(pval,"fdr"))
  
  if(!is.null(pv)){
    info = adj.pval < pv
  } else {
    info = adj.pval < adj.pval[as.integer(perc_genes/100*length(adj.pval))]
  }
  
  if(PLOT){
    if(min(means)<=0) ps = .1-min(means)
    if(min(means)>0)  ps = 0
    xg  =  1000^(seq( min(log(means+ps)), max(log(means+ps)), length.out=5000 ))
    
    vfit  =  a1/xg + a0 #estimated technical variation
    vfit_biol = psia1theta/xg + a0 + minBiolDisp # expected total variation
    
    xlabel = "log[ Average normalized read count]"
    smoothScatter(log(means+ps),log(cv2),xlab=xlabel,ylab="log[ Squared coefficient of variation (CV^2)]")
    points(log(means+ps),log(cv2),col="gray")
    # lines(log(xg[which(vfit>0)]+ps), log(vfit[which(vfit>0)]), col="black", lwd=3,lty=2,ylim=range(cv2) )
    lines(log(xg[which(vfit_biol>0)]+ps), log(vfit_biol[which(vfit_biol>0)]), col="black", lwd=3,ylim=range(cv2) )
    lines(log(xg[which(vfit_biol>0)]+ps),log(vfit_biol[which(vfit_biol>0)] * qchisq(0.975,m-1)/(m-1)),lty=2,col="black")
    lines(log(xg[which(vfit_biol>0)]+ps),log(vfit_biol[which(vfit_biol>0)] * qchisq(0.025,m-1)/(m-1)),lty=2,col="black")
    points(log(means[names(which(info)[TRUE])]+ps),log(cv2[names(which(info)[TRUE])]),col=2,cex=.5)
    
  }
  
  return(info)
}


## ------------------------------------------------------------------------
run_DANB = function(counts,save_plot=T,method="NBDrop",cutoff=NULL, perc_genes = NULL){
  
  if((is.null(perc_genes)&is.null(cutoff))){
    stop("Please provide exactly one of either cutoff or perc_genes")
  }
  
  if(!(is.null(perc_genes)|is.null(cutoff))){
    stop("Please provide either cutoff or perc_genes, not both!")
  }
  
  library(M3Drop) #note that you require version > 2.0 (not on bioconductor yet)
  
  if(!method%in%c("NBDrop","NBDisp")){
    stop("Invalid method selected. Please choose one of \"NBDrop\", \"NBDisp\"")
  }
  
  # fitting the dropout model (DANB)
  fit = NBumiFitModel(counts) #this fits a DANB model
  fit$mus = t(sapply(fit$vals$tjs, function (x) x * fit$vals$tis/fit$vals$total))
  
  size_coeffs = NBumiFitDispVsMean(fit, suppress.plot = T)#get coefcients of mean-dispersion fit
  smoothed_size = exp(size_coeffs[1] + size_coeffs[2] * log(fit$vals$tjs/fit$vals$nc)) #predicted dispersions per gene
  size_mat = matrix(rep(smoothed_size, times = fit$vals$nc), ncol = fit$vals$nc, byrow = FALSE)
  exp_ps  =  (1 + fit$mus/size_mat)^(-size_mat) #these are the fitted values per cell and gene
  exp_tot = rowSums(exp_ps) #this is the total predicted molecules per gene
  
  plot_dt = data.table(Dropout_rate = fit$vals$djs/fit$vals$nc,
                       expression = fit$vals$tjs/fit$vals$nc,
                       sizes = fit$sizes,
                       pred_sizes = smoothed_size,
                       predicted_dropouts=exp_tot/fit$vals$nc)
  
  if(method=="NBDrop"){
    NBumiCheckFitFS(counts,fit,suppress.plot = T) #check whether the fitted droout rates are well correlated with observed ones (i.e. number of zeroes)
    pvals = NBumiFeatureSelectionCombinedDrop(fit) #ranks genes by difference from expected dropout rates
    
    if(is.null(cutoff)){ cutoff = sort(pvals)[as.integer(perc_genes/100*length(pvals))]}
    
    info = rownames(counts)%in%names(pvals[which(pvals<cutoff)]) #select the top features based on dropout rates
    plot_dt[,info:=info]
    
    p = ggplot(plot_dt,aes(x=log10(expression),y=Dropout_rate)) +
      geom_point(aes(color=info),alpha=0.7,size=2) + 
      geom_line(aes(x=log10(expression),y=predicted_dropouts),color=colors()[30],size=1.2)+
      theme_bw() + xlab("Expression [log10]") + ylab("Dropout Rate")+
      ggtitle("Dropout rate vs. Expression")+ theme(text = element_text(size=17))+
      scale_color_manual(values = colors()[c(226,32)],name="is_outlier")
    print(p)
    if(save_plot){
      ggsave(p,file=file.path(plotdir,"NBdrop_plot.pdf"),height=6,width=8)
    }
  } else {
    resids = NBumiFeatureSelectionHighVar(fit) #ranks genes by difference from fitted mean-dispersion power law
    if(is.null(cutoff)){cutoff = perc_genes/100}
    info = rownames(counts)%in%names(resids[which(resids<quantile(resids,cutoff))])
    plot_dt[,info:=info]
    plot_dt[,est_cv2:=(1+expression/sizes)/expression] #predicted variance according to DANB model
    plot_dt[,ave_cv2:=(1+expression/pred_sizes)/expression] #predicted variance according to linear fit of dispersions
    
    p = ggplot(plot_dt,aes(x=log10(expression),y=log10(est_cv2))) +
      geom_point(aes(color=info),alpha=0.7,size=2) + 
      geom_line(aes(x=log10(expression),y=log10(ave_cv2)),color=colors()[30],size=1.2)+
      theme_bw() + xlab("Expression [log10]") + ylab("Estimated CV^2 [log10]")+
      ggtitle("Mean - Dispersion Relationship")+ theme(text = element_text(size=17))+
      scale_color_manual(values = colors()[c(226,32)],name="is_outlier")
    print(p)
    if(save_plot){
      ggsave(p,file=file.path(plotdir,"NBdisp_plot.pdf"),height=6,width=8)
    }
  }
  return(info)
}

## ------------------------------------------------------------------------
GO_to_gene = function(go_id, organism = "human"){
  if(organism == "human"){
    library(org.Hs.eg.db)
    gene_eg = get(go_id,org.Hs.egGO2ALLEGS) # retrieve all entrez gene identifiers mapping to go_id
    gene_ens = unlist(lapply(gene_eg, function(x) get(x,org.Hs.egENSEMBL))) #convert to ensembl
  } else if(organism == "mouse"){
    library(org.Mm.eg.db)
    gene_eg = get(go_id,org.Mm.egGO2ALLEGS) # retrieve all entrez gene identifiers mapping to go_id
    gene_ens = unlist(lapply(gene_eg, function(x) get(x,org.Mm.egENSEMBL))) #convert to ensembl
  } else {stop("Organism has to be human or mouse.")}
  return(gene_ens)
}

## ------------------------------------------------------------------------
build_adjacency_matrix = function(mat,cutoff="auto", is_similarity = F){
  library(Ckmeans.1d.dp)
  if(!is_similarity){
    message("Computing cell Pearson correlation coefficient")
    corr.cells = cor(mat,method="pearson")
  } else {corr.cells = mat}
  
  adj.corr = corr.cells
  
  if(cutoff=="auto"){
    # we find the best correlation cutoff by looking for a "valley"
    # in the histogram of correlations. This function attempts to set the
    # cutoff automatically, but might not always succeed...
    
    # if there are more than 500 cells, randomly sample 500 correlations
    if(dim(corr.cells)[1]>500){
      idx = sample(seq(dim(corr.cells)[1]),size=500)
    } else {idx = seq(dim(corr.cells)[1])}
    
    freqs = hist(corr.cells[idx,idx],breaks=dim(corr.cells[idx,idx])[1]/10)
    k1d = Ckmeans.1d.dp(corr.cells,k=2)
    cutoff = max(as.vector(corr.cells)[which(k1d$cluster==1)])
    abline(v=cutoff,col="red")
  } else if (is.numeric(cutoff)){cutoff=cutoff} else {
    stop("Please provide a numeric value for corr.cutoff or set to \"auto\"")
  }
  
  message("Building the adjacency matrix")
  adj.corr[adj.corr<cutoff]=0
  adj.corr[adj.corr>0] = 1
  return(list(adj=adj.corr,cor=corr.cells,cutoff=cutoff))
}

MCLcell.clust=function(adj_list,selfloop=T,mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(igraph)
  
  adj = adj_list$adj
  corr.cells = adj_list$cor
  corr.cutoff = adj_list$cutoff
  
  if(!selfloop) diag(adj)=0 # no self-loop for MCL
  message("Building Graph")
  graphs = get.data.frame( graph.adjacency(adj), what = "edges") # gene connection for graphs
  graphs = data.frame(graphs,CORR=sapply(seq(dim(graphs)[1]), function(i) corr.cells[graphs$from[i],graphs$to[i]] -corr.cutoff))
  
  write.table(graphs, file = "tmp.mcl.inp",row.names=F,col.names=F,sep = " ")
  message("Running MCL")
  system(paste0(mcl_path, " tmp.mcl.inp --abc -o tmp.mcl.out"))
  x = scan("tmp.mcl.out", what="", sep="\n")
  MCL.cells = strsplit(x, "[[:space:]]+")
  MCL.cells = lapply(seq(length(MCL.cells)), function(i){
    tmp = sapply(seq(length(MCL.cells[[i]])),function(j){
      gsub('\"','',MCL.cells[[i]][j])
    })
  })
  system("rm tmp.mcl.inp tmp.mcl.out")
  
  groups.MCL = matrix(rep(-1,dim(corr.cells)[2]),ncol=1)
  rownames(groups.MCL) = colnames(corr.cells)
  for(i in seq(length(MCL.cells))) groups.MCL[MCL.cells[[i]],]=i
  
  #if necessary, collapse all clusters containing only 1 cell to a big "unassigned"
  groups.MCL[groups.MCL %in% names(table(groups.MCL)[which(table(groups.MCL)==1)])] = 0
  
  return(groups.MCL)
}


## ------------------------------------------------------------------------
run_dbscan = function(dist,eps="auto",min_pts,tol=0.01){
  library(dbscan)
  #automatic determination of eps (the "elbow" in the kNNdistplot)
  if(eps=="auto"){
    kNNdist = sort(kNNdist(dist,min_pts))
    i = seq(1,length(kNNdist),as.integer(0.001*length(kNNdist)))
    slope_prev = 100
    for(indx in i){
      slope = kNNdist[indx]/indx
      if(slope_prev>=slope-tol*slope){
        slope_prev = slope
      } else {
        elbow = indx
        break
      }
    }
    eps = kNNdist[elbow]
    print(paste("Epsilon: ",eps))
  } else if(!is.numeric(eps)){
    stop("Please provide a value for eps or set it to \"auto\"")} else {eps=eps}
  
  kNNdistplot(dist,k=min_pts)
  abline(h=eps,col="red")
  res = dbscan(dist,eps = eps,minPts = min_pts)
  return(res$cluster)
}

## ------------------------------------------------------------------------
gmm_main = function(norm_exprs=NULL,pc_scores=NULL,n_comp=10,do_crossval=T, model_type = "VVI",
                     best_k=NULL,tolerance = 1.96,k=1:10,n_cores = 4, return_model = F){
  library(mclust)
  library(parallel)
  
  if(is.null(pc_scores)){
    if(is.null(norm_exprs)){
      stop("Missing expression values. Please provide either a matrix of normalized counts or pre-computed PCA scores.")
    }
    print("Running PCA...")
    pca = pca_stuff(norm_exprs)
    top_pc_scores = pca$x[,1:n_comp]
  } else {
    top_pc_scores = pc_scores[,1:n_comp]
  }
  
  if(do_crossval){
    #fit models with k-fold crossvalidation
    folds = 1:10
    n_folds = length(folds)
    # this randomly determines which samples should be excluded from
    # model fitting during each fold of the cross-validation. 
    idx = 	sample(rep(1:n_folds, length.out = nrow(top_pc_scores))) 
    
    #Set up parallel processing
    library(parallel)
    cl = makeCluster(n_cores) 
    funs = as.character(lsf.str(envir=parent.frame())) #let clusters know about functions in workspace
    clusterExport(cl,funs)
    
    #benchmark
    #time = system.time(parSapply(cl,folds,cross_val,data=testset,idx=idx,structure='VVI',components=k))
    print("Determining number of clusters...")
    likes = parSapply(cl,folds,cross_val,data=top_pc_scores,idx=idx,structure=model_type,components=k)
    stopCluster(cl)
    
    mean_likes = apply(likes,1,function(x) sum(x[which(is.finite(x))])/length(which(x!=0 & is.finite(x))))
    sd_likes = apply(likes,1,function(x) sd(x[which(x!=0 & is.finite(x))]))
    sd_likes = sd_likes[which(!is.na(sd_likes))]
    mean_likes = mean_likes[which(!is.na(sd_likes))]
    best_idx = which(mean_likes==max(mean_likes))
    ci_likes = mean_likes[best_idx]-tolerance*sd_likes[best_idx]
    
    best_k_idx = min(which(mean_likes>=ci_likes)) #smallest numebr of components that
    #fit reasonably well
    best_k = k[best_k_idx]
    mean_likes[best_k_idx]
    col=rep(1,n_folds)
    col[best_k_idx] = 33
    
    # Plot likelihood vs number of clusters
    # This should look like a ROC curve ideally. The best model should have
    # a high likelihood with the smallest no. of clusters, i.e. be the one
    # where the slope decreases. If there is no clear decrease in the
    # slope, this means that pretty much any model fits equally well/bad and that
    # most likely the clustering produces a nonsensical result
    like_dt = data.table(cluster=k,color=as.factor(col),likes)
    like_dt_melt = melt(like_dt,id.vars=c("cluster","color"),val="log-Likelihood",var="fold")
    p = ggplot(like_dt_melt[`log-Likelihood`!=0 & is.finite(`log-Likelihood`)],aes(x=as.factor(cluster),y=`log-Likelihood`,fill=color)) +
      geom_boxplot() + labs(x = "Number of clusters",
                            title = paste0("No. clusters v.s. log-likelihood, ",n_folds,"-fold crossvalidation"),
                            subtitle = paste0(n_comp," principal components used to calcualte model")) + 
      theme_bw() +scale_fill_manual(values=c("white","red"),guide=F)
    
    #ggsave(p,file=file.path(plotdir,paste0("mclust_crossval_",n_comp,".pdf")),height=5,width=7)
    print(p)
    
    print(paste0("Found ",best_k," clusters."))
  } else {best_k = best_k}
  
  if(is.null(best_k)){
    stop("Please provide a value for best_k or set do_crossval = TRUE")
  }
  print("Assigning cells to clusters...")
  #assignments of cells to clusters
  model = calc_model(top_pc_scores,best_k,model_type)
  cluster_assignments = model$classification
  if(return_model){
    return(model)
  } else {
    return(cluster_assignments)}
}


## ------------------------------------------------------------------------
#do pca
pca_stuff = function(log_data_hv,scale_pca=T,center_pca=T){
  pca = prcomp(t(log_data_hv[,-1,with=F]),scale=scale_pca,center=center_pca)
  return(pca)
}

#Fitting GMM with mclust:
##############################################################################
# function to calculate different models
# k = number of compnents
# structure = model structure / constraints. See mclustModelNames for details.

calc_model = function(data,k,structure){
  return(Mclust(data,G=k,modelNames = structure,initialization=
                  list(subset=sample(1:nrow(data),size=as.integer(nrow(data)*4/5)))))
}

#############################################################################
#Functions to calculate log-likelihood out of what mclust returns

# Probability density function for a Gaussian mixture
# Presumes the mixture object has the structure used by mclust

dnormalmix = function(x,mixture,log=FALSE) {
  lambda 	= mixture$parameters$pro
  k 		= length(lambda)
  # Calculate share of likelihood for all data for one component
  like_component = function(x, component) {
    lambda[component] * dmvnorm(
      x,
      mean = mixture$parameters$mean[,component],
      sigma 	= mixture$parameters$variance$sigma[,,component]
    )
  }
  # Create array with likelihood shares from all components over all data
  likes 	= sapply(1:k, like_component ,x = x)
  # Add up contributions from components
  d 	= rowSums(likes)
  if (log) {
    d 	= log(d)
  }
  return(d)
}

# Log likelihood function for a Gaussian mixture, potentially on new data
loglike_normalmix = function(x,mixture) {
  loglike 	= dnormalmix(x, mixture, log = TRUE)
  return(sum(loglike))
}

###############################################################################
#Cross validation things
#Cross validation
#data = input data
#idx = a random sample of folds (e.g. 11432...)
#fold = the current fold
#structure = model structure for mclust (e.g. 'VVI)
#components = a vector containing the numebr of components for which to test models

cross_val = function(fold,data,idx,structure,components){
  #library(mclust,lib.loc = .libPaths()[[2]]) #for the VM
  library(mclust)
  library(mvtnorm)
  like_test = c()
  for(k in components){
    out = tryCatch(
      {
        calc_model(data[which(idx!=fold),],k,structure)
      },
      error=function(cond){
        #try to find another model
        out2 = 0
        counter = 0
        while(out2 == 0 && counter<=5){
          out2 = tryCatch(
            {
              calc_model(data[which(idx!=fold),],k,structure)
            },
            error = return(0),
            warning=return(0),
            finally= {counter = counter +1}
          )
        }
        message('There was an error: \n')
        message(cond)
        write.csv(cond,'gmm.log',append=TRUE)
        return(out2)
      },
      warning=function(cond){
        #try to find another model
        out2 = 0
        counter = 0
        while(out2 == 0 && counter<=10){
          out2 = tryCatch(
            {
              calc_model(data[which(idx!=fold),],k,structure)
            },
            error = return(0),
            warning=return(0),
            finally={counter = counter +1}
          )
        }
        message('There was a warning: \n')
        message(cond)
        write.csv(cond,'gmm.log',append=TRUE)
        return(out2)
      },
      finally={
        message('\n done.')
      }
    )
    if(class(out)=='Mclust'){  
      like_test = append(like_test,loglike_normalmix(data[which(idx==fold),],out))
    }
    else{
      like_test = append(like_test,0)
    }
  }
  return(like_test)
}

## ------------------------------------------------------------------------
seurat_clustering = function(sce,vars.to.regress=NULL,res=0.6,n_comp=10){
  
  library(Seurat)
  #make SEURAT object, scale and optionally regress out confounders
  tmp_seurat = CreateSeuratObject(raw.data = counts(sce))
  tmp_seurat@data = norm_exprs(sce) #add the normalized values
  
  # This next step is a bit of cheating. Seurat expects us to run the complete
  # workflow on the same object and checks whether data have been normalized
  # by checking if object@calc.params$NormalizeData$normalization.method exists.
  # Since we provided normalized values, and do not want to re-run normalization,
  # we just put a dummy value in that slot.
  
  tmp_seurat@calc.params$NormalizeData = list(normalization.method ="dummy")
  
  if(!is.null(vars.to.regress)){
    if(any(!vars.to.regress%in%names(pData(sce)))){
      stop("Variables to regress out have to be column names in pData(sce)")
    }
    tmp_seurat = AddMetadata(object = tmp_seurat, metadata = pData(sce)[,vars.to.regress])
    tmp_seurat = ScaleData(object = tmp_seurat,vars.to.regress=vars.to.regress)
  } else {
    tmp_seurat = ScaleData(object = tmp_seurat)
  }
  
  tmp_seurat = RunPCA(object = tmp_seurat, pc.genes = rownames(sce), do.print = FALSE) 
  tmp_seurat = FindClusters(object = tmp_seurat, reduction.type = "pca", dims.use = 1:n_comp, 
                       resolution = res, print.output = 0, save.SNN = TRUE)
  seurat_assignment = tmp_seurat@ident
  return(seurat_assignment)
}


## ------------------------------------------------------------------------
cellsius_main = function(sce,group_id,min_n_cells=10,verbose = T, min_fc = 2,fc_between_cutoff = 1,
                                     organism = "human", corr_cutoff = NULL, iter=0, max_perc_cells = 50,
                                    mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(Ckmeans.1d.dp)
  
  expr_dt = data.table(gene_id = rownames(sce),norm_exprs(sce))
  expr_dt_melt = melt(expr_dt,id.vars="gene_id",val="expr",var="cell_idx")
  expr_dt_melt = merge(expr_dt_melt,
                       data.table(cell_idx=colnames(sce),main_cluster=as.character(pData(sce)[,group_id])),
                       by="cell_idx")
  
  #Identify genes with significant bimodal distribution
  
  expr_dt_melt[,c("N_cells","within_p","pos0","pos1","Dpos"):=cellsius_find_bimodal_genes(expr,min_n_cells = min_n_cells, max_perc_cells = max_perc_cells),by=c('gene_id','main_cluster')]
  expr_dt_melt[,sig := within_p<100 & Dpos > min_fc] 
  expr_dt_melt[sig==T, within_adj_p:=p.adjust(within_p),by=c('cell_idx')] #correct for multiple testing, only consider genes where test has actually been run
  expr_dt_melt[,sig:=within_adj_p<0.1] 
  expr_dt_melt = expr_dt_melt[gene_id %in% expr_dt_melt[!is.na(sig) & sig==T]$gene_id] 
  
  # If no bimodal gene were found, exit and return NA
  if(dim(expr_dt_melt)[1] == 0){
    print("No genes with bimodal distribution found, returning NA.")
    return(NA)
  }
  # Check whether these genes are specific to the subcluster
  
  for(clust in unique(expr_dt_melt$main_cluster)){
    expr_dt_melt = expr_dt_melt[,paste0(clust,"_",c("p_between","fc")):=cellsius_test_cluster_specificity(
      expr,main_cluster,clust, fc_between_cutoff = fc_between_cutoff),by="gene_id"]
    
    expr_dt_melt[main_cluster==clust,keep:=(expr_dt_melt[main_cluster==clust][[paste0(clust,"_p_between")]] < 0.1)]
  }
  
  expr_dt_melt = expr_dt_melt[keep==TRUE & !is.na(sig)]
  
  # If there are still non-specific genes, discard them (this can happen for
  # very high expressed genes like mitochondrial genes)
  expr_dt_melt[,n_clust_per_gene:=length(unique(main_cluster)),by='gene_id']
  expr_dt_melt = expr_dt_melt[n_clust_per_gene==1]
  expr_dt_melt[,n_clust_per_gene:=NULL]
  
  # Identify correlated gene sets with MCL
  expr_dt_melt = expr_dt_melt[,gene_cluster:=0]
  expr_dt_melt = cellsius_find_gene_sets(expr_dt_melt, corr_cutoff = corr_cutoff, mcl_path=mcl_path)
  
  # discard gene sets that only contain one gene (those are assigned to cluster 0)
  expr_dt_melt = expr_dt_melt[gene_cluster !=0 ]
  
  if(dim(expr_dt_melt)[1] == 0){
    print("No subclusters found, returning NA.")
    return(NA)
  }
  
  # Extract cell subclusters
  expr_dt_melt[,sub_cluster:=main_cluster]
  expr_dt_melt[,mean_expr := mean(expr), by = c('main_cluster','gene_cluster','cell_idx')]
  expr_dt_melt[,sub_cluster:=cellsius_sub_cluster(mean_expr,sub_cluster,gene_cluster, iter=iter),by=c('main_cluster','gene_cluster')]
  
  # Check how many cells belong to the subgroup relative to the total cluster size.
  # If a sub cluster contains more than max_perc_cells cells, discard it.
  clust_list = expr_dt_melt[,list(sub = length(unique(cell_idx))) ,by=c('sub_cluster','main_cluster')]
  clust_list[,tot := sum(sub)/(length(sub_cluster)/2), by= 'main_cluster']
  clust_list = clust_list[grep('_1$',sub_cluster)]
  clust_list[,perc:=sub/tot*100]
  discard_sub_clust = clust_list[perc > max_perc_cells]$sub_cluster
  discard_sub_clust = append(discard_sub_clust,gsub('_1$','_0',discard_sub_clust))
  
  expr_dt_melt = expr_dt_melt[!sub_cluster%in%discard_sub_clust]
  
  # If verbose is TRUE, print a summary of the results
  if(verbose){
    # annotate genes (only if verbose)
    gene_info = get_gene_annotations(unique(expr_dt_melt$gene_id),get_descriptions = T,
                                     organism = organism)
    expr_dt_melt = merge(expr_dt_melt,gene_info, by = 'gene_id')
    
    cellsius_print_summary(expr_dt_melt)
  }
  return(expr_dt_melt)
}

## ------------------------------------------------------------------------
##################################################
# STEP 1: Identify genes with bimodal distribution
##################################################

cellsius_find_bimodal_genes = function(expr, min_n_cells, max_perc_cells){
  
  #skip genes with 0 expression
  if(sum(expr)==0){
    return(list(-1,100,-1,-1,-1))
  }
  # run k-means
  k1d = Ckmeans.1d.dp(expr,k=2)
  # check if a cluster with more than n cells exists
  indx = which(k1d$size>min_n_cells) 
  if(length(indx)>1 ){
    
    # do statistic only if in pos2 cells are less than max_perc_cells% of the total cells in the cluster
    if(k1d$size[2] < round(length(expr)*max_perc_cells/100)){ 
      
      t1=tryCatch({t.test(expr[which(k1d$cluster==2)],y=expr[which(k1d$cluster==1)])},
                   error = function(cond){return(0)},
                   finally={}
      )
      
      if(!is.numeric(t1)){
        
        p1=t1$p.value
        N0=k1d$size[1] # number of cells where the gene is downregulated
        N1=k1d$size[2] # number of cells  where the gene is upregulated
        pos0=k1d$centers[1] 
        pos1=k1d$centers[2]
        Dpos=pos1-pos0
        return(list(N1,p1,pos0,pos1,Dpos))
      } #else {print(paste("ttest failed, dpos = ",pos1-pos0))} # for testing
    }
  }
  # if no cluster was found, return a list of dummy values
  return(list(-1,100,-1,-1,-1))
}

##################################################
# STEP 2: Check whether these genes are specific to one cell subgroup
###################################################

cellsius_test_cluster_specificity = function(exprs, cluster, current_cluster, fc_between_cutoff){
  
  in_clust = which(cluster == current_cluster)
  k1d = Ckmeans.1d.dp(exprs[in_clust],k=2)
  in_subclust = in_clust[which(k1d$cluster==2)]
  
  mean_in = mean(exprs[in_subclust])
  mean_out = mean(exprs[-in_subclust])
  mean_out_nozero = mean(exprs[-in_subclust][exprs[-in_subclust]>0])
  
  # If there are subclusters, but all cells outside the subcluster express 0,
  # set mean_out_nozero to 0
  if(length(in_subclust>0) && !any(exprs[-in_subclust]>0)){mean_out_nozero=0}
  
  fc = mean_in - mean_out
  
  ts = tryCatch({t.test(exprs[in_subclust],exprs[-in_clust])},
                error = function(cond){ return(0)})
  
  if(!is.numeric(ts)){pv = ts$p.value} else {
    #print(paste("ttest failed, fc = ",mean_in-mean_out_nozero)) #for testing only
    pv=999}
  
  if(!is.nan(mean_out_nozero) && (mean_in-mean_out_nozero < fc_between_cutoff)) pv = 999
  return(list(pv,fc))
}

#####################################################
# STEP 3: MCL clustering to find correlated gene sets
#####################################################

cellsius_find_gene_sets = function(expr_dt_melt, corr_cutoff = NULL, min_corr = 0.35, max_corr = 0.5,
                          mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(igraph)
  
  for(clust in unique(expr_dt_melt$main_cluster)){
    
    if(length(unique(expr_dt_melt[main_cluster == clust]$gene_id))==1) { next }
    
    mat = dcast.data.table(expr_dt_melt[main_cluster==clust], gene_id ~ cell_idx,
                           value.var = 'expr')
    mat = mat[rowSums(mat[,-1,with=F])!=0,]
    corr.mat = cor(t(mat[,-1,with=F]))
    dimnames(corr.mat) = list(mat$gene_id,mat$gene_id)
    
    if(is.null(corr_cutoff)){
      corr_cutoff = max(quantile(corr.mat[corr.mat!=1],0.95),min_corr)
      corr_cutoff = min(corr_cutoff, max_corr)}
    adj.corr = corr.mat
    adj.corr[adj.corr<corr_cutoff] = 0
    adj.corr[adj.corr>=corr_cutoff] = 1
    diag(adj.corr) = 0 # no self-loop for MCL
    
    graphs = get.data.frame( graph_from_adjacency_matrix(adj.corr), what = "edges") # gene connection for graphs
    
    # if a graph has no edges (i.e. all genes are uncorrelated), 
    # assign all genes to cluster "0" and go to next iteration
    if(dim(graphs)[1]==0){
      expr_dt_melt = expr_dt_melt[main_cluster == clust, gene_cluster := 0]
      next
    }
    
    graphs = data.frame(graphs,CORR=sapply(seq(dim(graphs)[1]), function(i) corr.mat[graphs$from[i],graphs$to[i]] -corr_cutoff))
    write.table(graphs, file = "tmp.mcl.inp",row.names=F,col.names=F,sep = " ")
    system(paste0(mcl_path, " tmp.mcl.inp --abc -o tmp.mcl.out"))
    x = scan("tmp.mcl.out", what="", sep="\n")
    y = strsplit(x, "[[:space:]]+")
    y = lapply(seq(length(y)), function(i){
      tmp = sapply(seq(length(y[[i]])),function(j){
        gsub('\"','',y[[i]][j])
      })
    })
    
    for(i in seq(length(y))){
      if(length(y[[i]] > 1)){
        expr_dt_melt = expr_dt_melt[main_cluster==clust & gene_id %in% y[[i]],gene_cluster:=i]
      }
    }
  }
  
  return(expr_dt_melt)
}  

############################################
# Step 4: Assign cells to subgroups
############################################

cellsius_sub_cluster = function(mean_expr,sub_cluster,gene_cluster, iter = 0){
  
  k1d = Ckmeans.1d.dp(mean_expr,k=2)$cluster
  cells_sub = (k1d==2)
  
  if(iter == 0){return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))}
  
  # if iter is set higher than 0, a second step of kmeans clustering. 
  # This will remove the lowest peak and can sometimes help to get a more
  # accurate classification.
  
  k1d = Ckmeans.1d.dp(mean_expr[cells_sub],k=2)$cluster
  
  if (max(k1d)>1) {
    cells_sub[cells_sub] = (k1d==2)
    return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))
  }
  return(paste0(sub_cluster,"_",gene_cluster,"_",0))
}

#######################################
# Step 5: Print summary
#######################################

cellsius_print_summary = function(expr_dt_melt){
  cat('--------------------------------------------------------\n',
      'Summary of rare cell types\n',
      '--------------------------------------------------------\n\n')
  for(clust in unique(expr_dt_melt$main_cluster)){
    
    if(!any(expr_dt_melt[main_cluster==clust]$gene_cluster!=0)){
      next
    }
    cat('Main cluster: ', clust,  '\n', '---------------\n')
    subclusts = unique(expr_dt_melt[main_cluster==clust & gene_cluster!=0][order(gene_cluster)]$gene_cluster)
    for(subclust in subclusts){
      
      cat('Subcluster: ', subclust, '\n',
          'Number of cells: ', 
          length(unique(expr_dt_melt[main_cluster==clust & 
                                       sub_cluster == paste(clust,subclust,1,sep="_")]$cell_idx)),
          '\n Marker genes: \n')
      
      print(unique(expr_dt_melt[main_cluster==clust & gene_cluster == subclust][,c("gene_id","symbol","description")]))
      cat('\n\n')
      
    }
  }
}

############################################
# OPTIONAL: Final assignment to unique clusters
# Note: This is different from the previous subcluster asignemnt, where a cell can potentially be
# a member of multiple subgroups.
############################################

cellsius_final_cluster_assignment = function(rare, sce, group_id, min_n_genes = 3){
  
  rare[,n_genes:=length(unique(gene_id)),by='sub_cluster']
  
  assignments = data.table(cell_idx = colnames(sce), pData(sce)[,group_id])
  names(assignments) = c('cell_idx', 'group')
  assignments$group = as.character(assignments$group)
  assignments = merge(assignments, rare[n_genes>=min_n_genes,c('cell_idx','main_cluster','sub_cluster')],by='cell_idx',all=T)
  assignments = unique(assignments)
  
  final_assignment = function(main_cluster,sub_cluster){
    
    if(length(sub_cluster)==1){
      if(is.na(sub_cluster) || grepl("0$",sub_cluster)){
        out = main_cluster
      } else {
        out = gsub('_\\d$','',sub_cluster)
      }
    } else {
      subclusts = gsub('_\\d$', '',sub_cluster[grepl("1$",sub_cluster)])
      out = paste(subclusts,collapse='-')
      if(out == ''){out = main_cluster}
    }
    return(out)
  }
  
  assignments[,final:=final_assignment(group,sub_cluster),by="cell_idx"]
  assignments = unique(assignments[,c('cell_idx','final')])
  
  out = data.frame(cluster = as.character(assignments$final), row.names = assignments$cell_idx)
  out = out[colnames(sce),,drop=F]

  return(out)
}


## ------------------------------------------------------------------------
# Visualize output of rare cell type algorithm on tSNE map
# tsne = RTsne object 
# rare = output of the rare cell types algorithm (a data.table)
plot_rare_cells = function(tsne,rare){
  
  tsne_dt = data.table(tSNE1 = tsne$Y[,1], tSNE2 = tsne$Y[,2], cell_idx = rownames(tsne$Y))
  tsne_dt = merge(tsne_dt, rare[,c('cell_idx','main_cluster','sub_cluster')],
                  by = c('cell_idx'), all = T)
  tsne_dt[is.na(main_cluster),main_cluster:='Other']
  tsne_dt[main_cluster == 'Other',sub_cluster:='none']
  tsne_dt[grepl('_0$',sub_cluster),sub_cluster:= 'none']
  
  setkey(tsne_dt, 'cell_idx')
  tsne_dt = unique(tsne_dt)
  
  rc_cols = brewer.pal(10,"Spectral")[rep(c(1,9,7,2,6,10,3,8),3)]
  
  p = ggplot(tsne_dt, aes(x = tSNE1, y= tSNE2)) + 
    geom_point(color = "darkgray", alpha = 0.5, size = 1.5)+
    theme_bw() + theme(text = element_text(size = 15))
  p = p + geom_point(data = tsne_dt[sub_cluster!='none'], aes(x=tSNE1, y=tSNE2, color = sub_cluster))+
    scale_color_manual(values = rc_cols) + guides(color = guide_legend(title = 'Subcluster'))
  return(p)                                                                     
}

## ------------------------------------------------------------------------
run_wilcoxon_test = function(sce,cl_id,cl_ref,
                             alpha=0.05,fc_cutoff=0.5,pseudocount=NULL){
  
  if(!is.null(pseudocount)){warning("The pseudocount argument is deprecated and will be removed!")}
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  in_clust = pData(sce)[,cl_id] %in% cl_ref
  pvals = apply(norm_exprs(sce),1,
                function(x) wilcox.test(x=x[in_clust],y=x[!in_clust])$p.value)
  fc = apply(norm_exprs(sce),1,
             function(x) mean(x[in_clust])-mean(x[!in_clust]))
  #log2((mean(2^x[in_clust]-1)+pseudocount)/(mean(2^x[!in_clust]-1)+pseudocount))
  
  out = data.table(gene_id = rownames(sce), pval = pvals, log2fc = fc)
  out[,adj_pval:=p.adjust(pval,method="fdr")]
  out[,DE_flag:=(adj_pval<alpha & abs(log2fc)>fc_cutoff)]
  
  return(out)
}

## ------------------------------------------------------------------------
run_limma = function(sce,cl_id,cl_ref,
                    alpha=0.05,fc_cutoff=0.5,count_thr=1,pct=50, method = "trend"){
  
  if(!method %in% c("voom","trend")){
    stop("Method has to be either \"voom\" or \"trend\".")
  }
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  res = as.factor(pData(sce)[,cl_id] %in% cl_ref)
  design = model.matrix(~0+res)
  colnames(design) = c("Rest","Cluster")
  rownames(design) = colnames(sce)
  
  # filter out genes not detected at a count of count-thr in at least
  # 50% of cells in at least one cluster. Outlier cells (clusters with only one cell) are ignored.
  clust_sizes = table(pData(sce)[,cl_id])
  clusts = names(clust_sizes[which(clust_sizes>1)])
  
  keep_mat = matrix(rep(NA,dim(sce)[1]*length(clusts)),ncol = length(clusts))
  for(i in seq(length(clusts))){
    keep_mat[,i] = rowSums(counts(sce)[,pData(sce)[,cl_id]==clusts[i]]>=count_thr)>=pct/100*length(which(pData(sce)[,cl_id]==clusts[i]))  
  }
  keep = apply(keep_mat, 1, function(x) any(x))
  
  #convert to DGEList and filter
  dge = convertTo(sce[keep,],"edgeR")
  contrast_matrix = limma::makeContrasts(Cluster-Rest, levels = design)
  
  if(method == "voom"){
    # transforming counts
    voomt = limma::voom(dge,plot=T,design = design)
    
    #do differential expression analysis on voom transformed data
    fit = limma::lmFit(voomt, design)
    fit2 = limma::contrasts.fit(fit,contrast_matrix)
    fit2 = limma::eBayes(fit2)
  } else {
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=1)
    fit = limma::lmFit(logCPM, design)
    fit2 = limma::contrasts.fit(fit,contrast_matrix)
    fit2 = limma::eBayes(fit2, trend = T)
  }
  
  diff_exp = limma::topTable(fit2,adjust="BH",number = dim(dge)[1])
  out = data.table(gene_id = rownames(diff_exp),diff_exp)
  out[,DE_flag:=as.factor(adj.P.Val < alpha & abs(logFC) > fc_cutoff)]
  
  return(out)
}

## ------------------------------------------------------------------------
run_MAST = function(sce,cl_id,cl_ref,n_cores = 8,nbins=10,min_per_bin=30,
                    alpha=0.05,fc_cutoff=0.5,norm=F,set_thresh=T){
  
  library(MAST)
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  options(mc.cores = n_cores)
  if(norm){
    sca = FromMatrix(norm_exprs(sce), pData(sce), fData(sce))} else {
      sca = FromMatrix(log2(counts(sce)+1), pData(sce), fData(sce))
    }
  # adaptive thresholding
  # note how the threshold moves with median expression
  if(set_thresh){
    message("Calculating expression thresholds...\n
            Check the MAST_theresholds plot. If there are no clear bimodal\n
            distributions, the thresholds are likely to be nonsense.\n
            If that is the case, re-run this function setting set_thresh = F")
    thres = thresholdSCRNACountMatrix(assay(sca), nbins = nbins, min_per_bin = min_per_bin)
    if(!any(thres$cutpoint!=0)){message("All cut points are zero. Try using a different
                                        value of nbins and min_per_bin or set set_thresh=FALSE")}
    par(mfrow=c(nbins%%4+1,4))
    plot(thres)
    dev.copy2pdf(file = file.path(plotdir,"MAST_thresholds.pdf"),width=8,height=3*(nbins%%4+1))
    par(mfrow=c(1,1))
    assays(sca) = list(thresh=thres$counts_threshold, counts=assay(sca))
    }
  
  cond=factor(colData(sca)[,cl_id]==cl_ref)
  cond=relevel(cond,"FALSE")
  colData(sca)$condition=cond
  # calculate the cellular detection rate as no. detected features / no. total features
  # and center it 
  colData(sca)$cngeneson = scale(pData(sce)$total_features/dim(sce)[1],scale=F)
  # fit model (note that this will take time for large datasets!!)
  message("Fitting models...")
  zlmCond = zlm(~condition + cngeneson, sca)
  summaryCond = summary(zlmCond, doLRT='conditionTRUE') 
  #extract the results as a data.table
  summaryDt = summaryCond$datatable
  fcHurdle = merge(summaryDt[contrast=='conditionTRUE' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionTRUE' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  names(fcHurdle)[1:3] = c("gene_id","pval","log2FC")
  fcHurdle[,DE_flag:=as.factor(fdr<alpha & abs(log2FC)>fc_cutoff)]
  return(fcHurdle)
}


## ------------------------------------------------------------------------
generic_scatterplot = function(dt, x_col, y_col, 
                        color = NULL, shape = NULL, size = NULL, alpha = 0.8, abs_size = 2){
  
  plot_dt = dt[,c(x_col,y_col),with=F]
  
  #by default, do not show any legends
  show_col = F
  show_shape = F
  show_size = F
  continuous = F
  discrete = F
  
  if(!is.null(color)){
    plot_dt[,c(color):=dt[[color]]]
    if(!is.numeric(dt[[color]])){
      show_col = guide_legend(color)
      discrete = T
    } else{
      show_col = guide_colorbar(color)
      continuous = T
    }
  } else {
    color = "dummy_color"
    plot_dt[,dummy_color:=factor(1)]
  }

  if(!is.null(shape)){
    show_shape = guide_legend(shape)
    plot_dt[,c(shape):=dt[[shape]]]
  } else {
    shape = "dummy_shape"
    plot_dt[,dummy_shape:=factor(1)]
  }
  
  if(!is.null(size)){
    show_size = guide_legend(size)
    plot_dt[,c(size):=dt[[size]]]
  } else {
    size ="dummy_size"
    plot_dt[,dummy_size:= 1]
  }
  
  p = ggplot(na.omit(plot_dt), aes_string(x=paste0("`",x_col,"`"),y=paste0("`",y_col,"`")))
  if(size == "dummy_size"){
    p = p+geom_point(aes_string(color=paste0("`",color,"`"),
                                shape=paste0("`",shape,"`"),
                                size=paste0("`",size,"`")),
                     alpha=alpha,size=abs_size)
  } else{
    p = p+ geom_point(aes_string(color=paste0("`",color,"`"),
                                 shape=paste0("`",shape,"`"),
                                 size=paste0("`",size,"`")),
                      alpha=alpha)}
  
  p = p + guides(color = show_col, shape = show_shape, size = show_size)+
    theme_bw()+theme(text=element_text(size=15))
  
  if(continuous){
    p = p+scale_color_continuous(low = "lightgray", high = "darkblue")
  }
    
  if(discrete){
    if(length(unique(na.omit(plot_dt[[color]])))==2){
      p = p+scale_color_manual(values = c("lightgray","darkred"))
    } else{p = p + scale_color_brewer(type = "qual",palette = 2)}
  }
    
  if(color == "dummy_color"){
    p = p + scale_color_manual(values = "darkgrey")
  }
  
  if(shape!="dummy_shape"){
    if(length(unique(na.omit(plot_dt[[shape]])))>9) {stop("Too many different shapes provided. Please provide a maximum of 9 shapes, otherwise, your plot will be a mess.")}
     p = p + scale_shape_manual(values = c(15,16,17,18,8,0,1,2,3))
  }
  
  return(p)
} 

## ------------------------------------------------------------------------
my_plot_PCA = function(counts=NULL,pca=NULL, alpha = 0.7, scale_pca = T, center_pca=T,comp=c(1,2),
                       color=NULL,size=NULL,shape=NULL,return_pca=F,title="PCA",abs_size=2, use_irlba=F){
  
  if(is.null(counts)&is.null(pca)){
    message('Please provide either a count matrix or pre-computed 
            principal component analysis as a prcomp object.')
  } else  if(is.null(pca)){
    if(use_irlba){
      library(irlba)
      if(packageDescription("irlba")$Version <= 2.3){
        stop("Please update irlba to version 2.3.2 (github). There is a bug in versions < 2.3 which results in unreliable output.")
      }
      pca = prcomp_irlba(t(counts),n=max(comp),center=center_pca,scale. = scale_pca)
      rownames(pca$x) = colnames(counts)} else{
        pca<-prcomp(t(counts), scale = scale_pca, center = center_pca)}
  } 
  
  pca_1.2<-cbind(pca$x[,comp[1]],pca$x[,comp[2]])
  
  if(!use_irlba){
    sdev1<-round((pca$sdev[comp[1]]**2)/sum(pca$sdev **2)*100,2)
    sdev2<-round((pca$sdev[comp[2]]**2)/sum(pca$sdev **2)*100,2)
  } 
  
  pca_1.2 = data.table(pca_1.2)
  names(pca_1.2) = paste0("PC",comp)
  
  if(!is.null(color)){
    pca_1.2[,color:=color[rownames(pca$x),]]
    setnames(pca_1.2,"color",colnames(color))
    color = colnames(color)
  }
  
  if(!is.null(size)){
    pca_1.2[,size:=size[rownames(pca$x),]]
    setnames(pca_1.2, "size", colnames(size))
    size = colnames(size)
  }
  
  if(!is.null(shape)){
    pca_1.2[,shape:=shape[rownames(pca$x),]]
    setnames(pca_1.2, "shape", colnames(shape))
    shape = colnames(shape)
  }
  
  p = generic_scatterplot(pca_1.2, x_col = paste0("PC",comp[1]), y_col = paste0("PC",comp[2]), color = color,
                          size = size, shape = shape, abs_size = abs_size, alpha = alpha)
  if(use_irlba){p = p+ggtitle(title)} else {
    p = p + xlab(paste("PC ",comp[1], " [",sdev1, "%]",sep="")) +
      ylab(paste("PC ", comp[2], " [",sdev2, "%]",sep="")) + ggtitle(title)
  }
  
  if(return_pca){
    return(list(plot = p, pca = pca))
  } else{
    return(p)
  }
  
}


## ------------------------------------------------------------------------
plot_pca_loadings = function(pca, comp = 1){
  
  loading_dt = data.table(pca$rotation[,comp])
  names(loading_dt) = "loading"
  loading_dt[,gene_id:=rownames(pca$rotation)]
  loading_dt = loading_dt[order(loading,decreasing=T)]
  
  n = length(loading_dt$loading)
  loading_dt = loading_dt[append(c(1:15),seq(n-15,n,1)),]
  
  loading_dt$gene_id = factor(loading_dt$gene_id, levels = loading_dt$gene_id)
  
  p = ggplot(loading_dt, aes(x=loading, y=gene_id))+
    geom_point()+theme_bw()
  
  return(p)
}


## ------------------------------------------------------------------------
my_plot_tSNE = function(counts=NULL,tsne=NULL,alpha = 0.7, color=NULL,abs_size = 2,
                        size=NULL,shape=NULL,return_tsne=F,is_distance=F,show_proportions=F,
                        n_comp = 50, scale_pca=F, use_irlba = F, title="tSNE"){
  
  if(is.null(counts)&is.null(tsne)){
    message('Please provide either a count matrix or pre-computed 
            tSNE map as an Rtsne object.')
  } else  if(is.null(tsne)){
    if(use_irlba){
      library(irlba)
      if(packageDescription("irlba")$Version <= 2.3){
        stop("Please update irlba to version 2.3.2 (github). There is a bug in versions < 2.3 which results in unreliable output.")
      }
      pca = prcomp_irlba(t(counts),n=n_comp,center=T,scale. = scale_pca)
      rownames(pca$x) = colnames(counts)
      tsne = Rtsne(pca$x,pca = F, initial_dims = n_comp, is_distance = F)
      rownames(tsne$Y) = colnames(counts)
    } else{
      tsne<-Rtsne(t(counts),initial_dims=n_comp,pca=T,
                  is_distance=is_distance,pca_scale=scale_pca)
      rownames(tsne$Y) = colnames(counts)
    }
  } 
  
  tsne_1.2 = data.table(tsne$Y)
  names(tsne_1.2) = c("tSNE1","tSNE2")
  
  if(!is.null(color)){
    if(show_proportions){
      if(is.numeric(color[[colnames(color)]])){stop("Proportions can only be calculated for discrete variables. Please provide your color scale as character or factor.")}
      color[[colnames(color)]] = paste0(color[[colnames(color)]]," [",
                                        round(calc_percentage(color[[colnames(color)]]),2),"%]")
    }
    tsne_1.2[,color:= color[rownames(tsne$Y),]]
    setnames(tsne_1.2,"color",colnames(color))
    color = colnames(color)
  }
  
  if(!is.null(size)){
    tsne_1.2[,size:=size[rownames(tsne$Y),]]
    setnames(tsne_1.2, "size", colnames(size))
    size = colnames(size)
  }
  
  if(!is.null(shape)){
    tsne_1.2[,shape:=shape[rownames(tsne$Y),]]
    setnames(tsne_1.2, "shape", colnames(shape))
    shape = colnames(shape)
  }
  
  p = generic_scatterplot(tsne_1.2, x_col = "tSNE1", y_col = "tSNE2", color = color,
                          size = size, shape = shape, abs_size = abs_size, alpha = alpha)
  p = p+ggtitle(title)
  
  if(return_tsne){
    return(list(plot = p, tsne = tsne))
  } else{
    return(p)
  }
  
  }


## ------------------------------------------------------------------------
make_pairs_plot = function(input_mat,main=""){
  ## puts histograms on the diagonal
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, 
         col = "lightgray", ...)
  }
  
  ## put (absolute) correlations on the upper panels,
  ## with size proportional to the correlations.
  
  panel.cor <- function(x, y, digits = 2, 
                        prefix = "", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="na.or.complete"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  pairs(input_mat,upper.panel=panel.cor,diag.panel=panel.hist,main=main)
}

## ------------------------------------------------------------------------
ggmultiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

