#####################################################
#
# scRNASeq pipeline functions
#
# PART I: Setup
# _______________________
#
# This script contains heleper fucntions for reading data, gene annotation and loading required libraries.
# 
# Authors:
#   Rebekka Wegmann (wegmann@imsb.biol.ethz.ch)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

################
# Library Calls
################

# note that these are only what is needed throughout the analysis
# as some of the other packages load a LOT of depndencies, you might have
# to unload them / restart R after using them and before being able to
# load any new packages

library_calls = function(){
  library(Rtsne) 
  library(ggplot2)
  library(data.table)
  library(scater)
  library(scran)
  library(RColorBrewer)
}

################################
# Reading data & Gene annotation
################################

# Helper function to read count matrix from csv. 

read_data = function(infile,experiment_id){
  counts  =  read.delim(infile, header=T, stringsAsFactors=FALSE)
  # matrix count format
  rownames(counts) = counts$Gene.Id
  counts = counts[,-1]
  colnames(counts) =  paste0(experiment_id,"_C",seq(dim(counts)[2]),sep="")
  return(counts)
}


#----------------------------------
# Annotating genes using ensembldb (the advantage is that it is faster than biomaRt, will
# return exactly one entry per gene including non-coding ones, and uses always
# the same version (regardless of bioconductor version), so this is used
# by default)

get_gene_annotations = function(gene_list,v=F,get_descriptions=T,organism = "human"){
  library(ensembldb)
  if(organism == "human"){
    library(EnsDb.Hsapiens.v86)
    edb = EnsDb.Hsapiens.v86
  } else if(organism == "mouse"){
    library(EnsDb.Mmusculus.v79)
    edb = EnsDb.Mmusculus.v79
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

get_gene_symbols_in_order = function(gene_list){
  # set get_descriptions to T if wanted
  gene_annot <- get_gene_annotations(gene_list, get_descriptions = F)
  setkey(gene_annot, 'gene_id')
  gene_annot <- gene_annot[gene_list]
  gene_names <- gene_annot$symbol
  return(gene_names)
}

#-----------------------------------
# Annotating genes using biomaRt

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


####################
# utility
####################

# merge two sce objects by adding columns to each other and fill in missing genes
# the combined object will contain only assays and metadata present in both sces
# if filter_genes is true, onlygenes detected in both experiments will be kept

merge_sce = function(sce1,sce2, filter_genes = T){
  if(any(colnames(sce1) %in% colnames(sce2))){
    stop('Column names must be unique!')
  }
  all_rownames = unique(c(rownames(sce1), rownames(sce2)))
  all_assays = intersect(names(assays(sce1)), names(assays(sce2)))
  all_coldata = intersect(names(colData(sce1)), names(colData(sce2)))
  all_rowdata = intersect(names(rowData(sce1)), names(rowData(sce2)))
  
  c = rbind(colData(sce1)[,all_coldata], colData(sce2)[,all_coldata])
  r = unique(rbind(rowData(sce1)[,all_rowdata], rowData(sce2)[,all_rowdata]))[all_rownames,]
  
  sce_combined = SingleCellExperiment(colData = c, rowData = r)
  
  
  for(i in 1:length(all_assays)){
    tmp = matrix(nrow = length(all_rownames), ncol = dim(sce1)[2] + dim(sce2)[2], data=0)
    rownames(tmp) = all_rownames
    colnames(tmp) = c(colnames(sce1), colnames(sce2))
    tmp[rownames(sce1), colnames(sce1)] = assays(sce1)[[which(assayNames(sce1) == all_assays[i])]]
    tmp[rownames(sce2), colnames(sce2)] = assays(sce2)[[which(assayNames(sce2) == all_assays[i])]]
    assays(sce_combined)[[i]] = tmp
  }
  names(assays(sce_combined)) = all_assays
  
  if(filter_genes){
    genes_in_both = intersect(rownames(sce1), rownames(sce2))
    sce_combined = sce_combined[genes_in_both,]
  }
  return(sce_combined)
  
}
