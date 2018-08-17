#' OS data
#'
#' Pre-processed gene expression data from high-grade osteosarcoma biopsies and sample characteristics were downloaded from GEO (GSE42352). We converted nuIDs to gene symbols using the annotation platform GPL10295. For genes with duplicate gene symbols, we selected the gene with the highest variance. Finally, we subsetted the data to the 53 patients for which 5 year metastasis information was available.
#'
#' @docType data
#' @keywords datasets
#' @name OSdata
#' @usage data(OSdata)
#' @format Variable "exp": Data frame containing expression data for 24998 genes and 53 samples. Variable "targets": Data frame containing information on whether patients developed metastases within 5 years or not, 53 samples and 2 columns.
NULL
