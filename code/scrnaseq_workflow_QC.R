#####################################################
#
# scRNASeq pipeline functions
#
# PART II: Quality control
# _______________________
#
# This script contains all functions called from the Quality Control section of the workflow.
# 
# Authors:
#   Rebekka Wegmann (wegmann@imsb.biol.ethz.ch)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

######################
# Quality Control
######################

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

#__________________
# Cell cycle phase annotation
#__________________

#using the scran package, which internally calls cyclone (default)

annotate_cell_cycle = function(sce, organism = "human", gene.names = rownames(sce)){
  if(organism == "human"){
    hs.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    assigned = cyclone(sce, pairs=hs.pairs, gene.names = gene.names)} else if (organism == "mouse"){
      mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
      assigned = cyclone(sce, pairs=mm.pairs, gene.names = gene.names)
    } else {stop("Organism has to be human or mouse.")}
  
  return(assigned)
}

#___________
# QC visualization
#___________

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
  mt_genes.amount = sce$pct_counts_MT 
  #mt gene content per cell
  plot(seq(length(mt_genes.amount)),(mt_genes.amount),pch=10,cex=1.5,col="gray",main=""
       , cex.axis=2,cex.lab=1.5,xlab="cell index",ylab="Ratio of MT-genes[%]")
  points(seq(length(mt_genes.amount))[mt_genes.amount<t],mt_genes.amount[mt_genes.amount<t],pch=10,cex=1.5)
  abline(h=t,col="red",lty=2,cex=5)
  
  #UMI vs no. genes colored by mt gene content
  plotColData(sce, x = "log10_total_features_by_counts",
                                y = "log10_total_counts",
                                colour_by = "pct_counts_MT")+
    xlab("Total detected features [log10]") + ylab("Total counts [log10]")+
    ggtitle("Total features vs. total counts, colored by MT content")
}

