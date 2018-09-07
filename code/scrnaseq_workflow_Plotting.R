#####################################################
#
# scRNASeq pipeline functions
#
# Plotting and dimensionality reduction
# _______________________
#
# This script contains helper functions to make PCA and tSNE plots, some other miscellaneous ggplot-wrappers and modified versions of some SC3 and scater plotting functions.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

##############################
# General dimensionality reduction and plotting
#############################

#plot PCA. Need to submit either the counts or a pre-calculated PCA
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

# Visualize the top PC loadings

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

#plot tSNE. Need to submit either the counts or a pre-calculated Rtsne object
# note thet for discrete color labels, the colors argument must be a factor
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

# Convenience function to calculate percentages of a vector of discrete variables
calc_percentage = function(x){
  perc = table(x)/sum(table(x))*100
  idx = sapply(x, function(x) which(names(perc)==x)) 
  return(as.numeric(perc[idx]))
}


# Putting multiple ggplots on one image
# The function is from here: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
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
#_________________________________________________
# Make a pairs plot (e.g. to compare size factors)
# Input:
# - input_mat = a matrix of columns that you want to plot against each other
# 
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

# Generic scatterplot using data.table and ggplot2
#__________________________________

# Convenience function to plot 2 columns of a data.table against each other
# dt =  data table
# x_col = name of the column to lot on x axis
# y_col = name of column to plot on y axis
# color = name of column to use as color values. If numeric, a continuous color scheme will be applied, if character or factor, colors will be discrete.
# shape = name of column to use for shape
# size = name of column to use for size
# alpha = transparency

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

# Function to launch the marker_vis shiny app (pre-alpha version and very buggy...)
launch_marker_vis_app = function(tsne,sce,marker_idx){
  plot_dt = data.table(tSNE1=tsne$Y[,1],tSNE2=tsne$Y[,2],
                       t(norm_exprs(sce)[marker_idx,,drop=F]))
  names(plot_dt)[-c(1,2)] = fData(sce)$symbol[marker_idx]
  assign("plot_dt",plot_dt,.GlobalEnv)
  runApp(file.path(code_dir,'marker_vis_app'))
}

# The following functions are modified from scater / SC3

## Custom version of the sc3_plot_markers function, taken from the SC3 package
custom_sc3_plot_markers = function (object, k, auroc = 0.85, p.val = 0.01, show_pdata = NULL, order_dend = F) 
{
  if (is.null(object@sc3$consensus)) {
    warning(paste0("Please run sc3_consensus() first!"))
    return(object)
  }
  if(order_dend){
    hc <- object@sc3$consensus[[as.character(k)]]$hc
  }
  dataset <- get_processed_dataset(object)
  if (!is.null(object@sc3$svm_train_inds)) {
    dataset <- dataset[, object@sc3$svm_train_inds]
  }
  add_ann_col <- FALSE
  ann <- NULL
  if (!is.null(show_pdata)) {
    ann <- SC3:::make_col_ann_for_heatmaps(object, show_pdata)
    if (!is.null(ann)) {
      add_ann_col <- TRUE
      rownames(ann) <- colnames(dataset)
    }
  }
  markers <- SC3:::organise_marker_genes(object, k, p.val, auroc)
  markers <- SC3:::markers_for_heatmap(markers)
  row.ann <- data.frame(Cluster = factor(markers[, 1], levels = unique(markers[, 
                                                                               1])))
  if(order_dend){
    order = hc$order
  } else { order = order(object[[paste('sc3',k,'clusters',sep='_')]])}
  
  groups = object[[paste('sc3',k,'clusters',sep='_')]][order]
  gaps = sapply(unique(groups), function(y) max(which(groups==y)))
  
  rownames(row.ann) <- rownames(markers)
  do.call(pheatmap::pheatmap, c(list(dataset[rownames(markers), order,
                                             drop = FALSE], show_colnames = FALSE, cluster_rows = FALSE, 
                                     cluster_cols = FALSE, gaps_col = gaps,
                                     annotation_row = row.ann, 
                                     annotation_names_row = FALSE, gaps_row = which(diff(markers[, 
                                                                                                 1]) != 0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
}


