#####################################################
#
# scRNASeq pipeline functions
#
# PART III: Normalization
# _______________________
#
# This script contains a wrapper function to different normalizations, based on the scran R package.
# 
# Authors:
#   Rebekka Wegmann (wegmann@imsb.biol.ethz.ch)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

#########################################
# Normalizations
# Input:
# - sce = an SCESet
# - method = the method you want to use to normalize. Choices are:
#     TC = total count normalization (multiply this by 10^6 to get CPMs)
#     UQ = upperquartile
#     RLE = relative log-expression, as in DESeq2
#     TMM = trimmed mean of M-values, as in edgeR
#     scran (default) = Lun sum factors, implemented in scran package
# Output:
# - the normalized expression values in the norm_exprs slot of the SCESet
#########################################

normalize_counts = function(sce,method = "scran"){
  #calculate size factors according to method
  switch(method,
         "TC" ={sizeFactors(sce) = sce$sum/mean(sce$sum)},
         "RLE" = {sf =  edgeR::calcNormFactors(counts(sce), method = "RLE")
                  sizeFactors(sce) =sf*sce$total_counts/mean(sf*sce$total_counts)},
         "TMM" = {sf =  edgeR::calcNormFactors(counts(sce), method = "TMM")
                  sizeFactors(sce) =sf*sce$total_counts/mean(sf*sce$total_counts)},
         "UQ" = {sf =  edgeR::calcNormFactors(counts(sce), method = "upperquartile")
                  sizeFactors(sce) =sf*sce$total_counts/mean(sf*sce$total_counts)}, 
         "scran" = {
           clusters = quickCluster(sce)
           sizes = seq(20, 100, 5)
           if (min(table(clusters)) > max(sizes)) {
             sce = computeSumFactors(sce, clusters = clusters, sizes = sizes)
           } else{
             message("Clustering of cells failed, using global scaling factors")
             sce = computeSumFactors(sce)
             if (any(sizeFactors(sce) < 0)) {
               warning(
                 "Negative size factors generated. Most likely, this is due to some cells having very low total feature counts. Consider using more stringent QC cutoffs."
               )
             }
           }
         })
  
         sce = scater::logNormCounts(sce)
         norm_exprs(sce) = logcounts(sce)
  
  return(sce)  
}
