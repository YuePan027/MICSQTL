#' Example protein expression data
#'
#' A subset of proteomic data PsychENCODE
#'
#' @format ## `protein_data`
#' A data frame with 2242 rows (protein) and 100 columns (sample).
"protein_data"


#' Example SNP data
#'
#' A data set stores the information of genetic variants at each location.
#'
#' @format ## `SNP_data`
#' A sparse matrix with 344388 rows (SNP) from one chromosome and 100 columns (sample, should match the sample in `protein_data`).
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2) for a sample at a genetic location.
"SNP_data"


#' Protein annotation data
#'
#' Annotations of each protein from `protein_data`
#'
#' @format ## `anno_protein`
#' A data frame with 2242 rows and 3 columns:
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start}{Start position}
#'   \item{Symbol}{Unique symbol for each protein}
#'   ...
#' }
"anno_protein"


#' SNP annotation data
#'
#' Annotations of each SNP from `SNP_data`
#'
#' @format ## `anno_SNP`
#' A data frame with 344388 rows and 3 columns:
#' \describe{
#'   \item{CHROM}{Chromosome}
#'   \item{POS}{Position}
#'   \item{ID}{Unique identifier for each SNP}
#'   ...
#' }
"anno_SNP"


#' Signature matrix
#'
#' A signature matrix which serves as a reference of known cellular signatures.
#'
#' @format ## `ref_data`
#' A data frame with 2242 rows (proteins) and 4 columns (cell types).
"ref_data"


#' Example gene expression data
#'
#' A subset of RNA-seq data PsychENCODE
#'
#' @format ## `gene_data`
#' A data frame with 2000 rows (genes) and 100 columns (sample).
"gene_data"


#' Metadata
#'
#' Metadata with interested phenotype.
#'
#' @format ## `meta`
#' A data frame with 100 rows (sample) and 2 columns (disease status and gender).
"meta"


#' Deconvoluted proportion based on gene expression
#'
#' Cell-type proportion estimates from gene expression
#'
#' @format ## `prop_gene`
#' A data frame of cell-type proportion estimates for each sample based on RNA-seq expression.
"prop_gene"
