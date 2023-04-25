#' Example data
#'
#' The example input files are organized as a `SummarizedExperiment `object.
#' 
#' #' A `SummarizedExperiment `object with following example data:
#'   - protein_data: A subset of proteomic data PsychENCODE with 2242 rows
#' (protein) and 127 columns (sample).
#' 
#'   - anno_protein: A data frame with 2242 rows and 4 columns
#' (Chr, Start, End, Symbol) as annotations of each protein from `protein_data`.
#' 
#'   - ref_protein: A signature matrix with 2242 rows (protein) and
#' 4 columns (cell types), which serves as a reference of known
#' cellular signatures.
#' 
#'   - gene_data: A data frame with 2867 rows (genes) and
#' 127 columns (sample).
#' 
#'   - ref_gene: A signature matrix with 4872 rows (genes) and 5 columns
#' (cell types), which serves as a reference of known cellular signatures.
#' 
#'   - SNP_data: A sparse matrix with 2000 rows (SNP), which stores the
#' information of genetic variants at each location from one chromosome
#' and 127 columns (sample, should match the sample in `protein_data`).
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2)
#' for a sample at a genetic location.
#' 
#'   - anno_SNP: A data frame with 2000 rows and 3 columns (CHROM, POS, ID),
#' which stores Annotations of each SNP from `SNP_data`.
#' 
#'   - meta: A data frame with 127 rows (sample) and 2 columns
#' (disease status and gender) as metadata.
#'
#' @name se
#' @docType data
#' @keywords datasets
"se"