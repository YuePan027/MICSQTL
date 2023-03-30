#' Estimation of cellular composition using cross-source information
#'
#' This function returns a `SummarizedExperiment` object including cross-source cell-type proportion estimates for each sample.
#'
#' This is a function developed to implement cell-type proportion deconvolution using initial proportion generated from a different source.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression data frame contained in `counts` slot.
#' @param ini_prop Initial cell proportion generated from a different source (e.g. gene expression).
#' It should be a data frame with rows as matched samples as in protein expression data, and columns as cell types.
#' The sum of cell proportions for each samples should be 1.
#' @param mrk_prot Marker proteins used for cell proportion estimates.
#'
#' @return A `SummarizedExperiment`. The cross-source cell-type proportion estimates for each sample will be stored as an element (`cross_prop`) in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @import TCA
#'
#' @export
#'
#'
#'
#'
cross_prop <- function(se, ini_prop, mrk_prot){
    if (!all(colnames(assay(se)) == rownames(ini_prop))){
        stop("Samples in protein_data do not match that in initial proportion")
    }
    mrk_prot <- intersect(rownames(assay(se)), mrk_prot)
    tca_res <- TCA::tca(X = assay(se)[mrk_prot,],
                   W = ini_prop,
                   refit_W = TRUE,
                   refit_W.sparsity = length(mrk_prot))
    cross_prop <- tca_res$W
    se@metadata$cross_prop <- cross_prop
    return(se)
}
