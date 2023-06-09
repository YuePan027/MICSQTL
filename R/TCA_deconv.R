#' Estimation of cellular expression in high-throughput data
#' from heterogeneous
#'  tissues
#'
#' This function returns a `SummarizedExperiment` object including
#' deconvoluted cellular sample-wise expression.
#'
#' This is a function developed to implement deconvolution for
#' cell-type-specific expression per bulk sample using `TCA`.
#'
#' @param se A `SummarizedExperiment` object with bulk expression data frame
#' contained in `counts` slot, and
#' cell-type proportion estimates for each sample contained as an element
#' (`prop`) in `metadata` slot.
#' @param test A character string indicate which data set to estimate the
#' cellular expression.
#' If not "array", it must contained in `metadata` slot.
#' @param prop A matrix of cellular composition. It has to be provided
#' unless `test` is set to be "array".
#'
#' @return A `SummarizedExperiment`. The results after `TCA` deconvolution
#' will be stored
#' as an element (`TCA_deconv`) in `metadata` slot.
#' It is a list with the length of the number of cell types (as in `prop`
#' in `metadata` slot).
#' Each element stores the deconvoluted protein expression matrix per bulk
#' sample within that cell type.
#'
#' @import SummarizedExperiment
#' @importFrom methods slot
#'
#' @export
#'
#' @examples
#' data(se)
#' se <- deconv(se, source = "protein", method = "nnls")
#' se <- feature_filter(se,
#'     target_protein = c("ABCA1", "ABCA2"),
#'     filter_method = c("allele", "distance"), filter_allele = 0.15,
#'     filter_geno = 0.05, ref_position = "TSS"
#' )
#' se <- TCA_deconv(se, prop = methods::slot(se, "metadata")$prop)
#'
TCA_deconv <- function(se, test = "array", prop) {
    if (test == "array") {
        res_tca <- TCA::tca(
            X = as.matrix(assay(se)),
            W = prop,
            refit_W = FALSE,
            verbose = FALSE
        )
        res_full <- TCA::tensor(tca.mdl = res_tca, X = as.matrix(assay(se)))
        names(res_full) <- colnames(prop)
        methods::slot(se, "metadata")$TCA_deconv <- res_full
    } else {
        res_tca <- TCA::tca(
            X = methods::slot(se, "metadata")[[test]],
            W = prop,
            refit_W = FALSE,
            verbose = FALSE
        )
        res_full <- TCA::tensor(
            tca.mdl = res_tca,
            X = methods::slot(se, "metadata")[[test]]
        )
        names(res_full) <- colnames(prop)
        methods::slot(se, "metadata")$TCA_deconv2 <- res_full
    }

    return(se)
}
