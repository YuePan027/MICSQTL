#' Cell-type-specific differential expression (csQTL)
#'
#' This function returns a `SummarizedExperiment` object including csQTL
#' proteins based on samples' genotype.
#'
#' This is a function developed to implement cell-type-specific differential
#' expression using `TOAST`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression data
#' frame contained in `counts` slot.
#' The information from genetic variants should be stored in a P
#' (the number of SNP) by N (the number of samples,
#' should match the sample in `counts` slot) matrix contained as an element
#' (`SNP_data`) in `metadata` slot.
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2)
#' for a sample at a genetic location.
#' The annotations of these SNP should be stored as an element (`anno_SNP`)
#' in `metadata` slot.
#' It should include at least the following columns: "CHROM" (which chromosome
#' the SNP is on),
#' "POS" (position of that SNP) and "ID" (a unique identifier for each SNP,
#' usually a combination of chromosome and
#' its position).
#' The information on cellular composition is required and stored as `prop` in
#' `metadata` slot.
#' It is an N (the number of samples, should match the sample in `counts` slot)
#' by K (the number of cell types) matrix.
#' This can be obtained by running `deconv()` before any filtering steps, or
#' use any source data directly.
#' @param BPPARAM For applying `bplapply`.
#'
#' @return A `SummarizedExperiment`. The csQTL results will be stored as an
#' element (`TOAST_output`) in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @import TOAST
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @importFrom S4Vectors metadata
#'
#' @export
#'
#' @examples
#' data(se)
#' se <- deconv(se, source = "protein", method = "nnls")
#' target_protein <- c("ABCA1")
#' se <- feature_filter(se,
#'     target_protein = target_protein,
#'     filter_method = c("allele", "distance"), filter_allele = 0.15,
#'     filter_geno = 0.05, ref_position = "TSS"
#' )
#' se <- csQTL(se)
#'
csQTL <- function(se, BPPARAM = bpparam()) {
    prop <- metadata(se)$prop
    assay(se) <- as.data.frame(assay(se))
    protein_tab <- as.data.frame(metadata(se)$target_dat)
    SNP_dat <- metadata(se)$SNP_data
    SNP_ID <- metadata(se)$anno_SNP$ID

    Res_TOAST <- lapply(
        seq_len(length(metadata(se)$choose_SNP_list)),
        function(x) {
            test_protein <- metadata(se)$choose_SNP_list[[x]]
            protein_name <-
                names(metadata(se)$choose_SNP_list)[x]
            message("csQTL test for protein ", protein_name, "\n")

            res <- bplapply(seq_len(length(test_protein)), function(i) {
                design <- data.frame(factor(SNP_dat[test_protein[i], ]))
                Design_out <- makeDesign(design, prop)
                fitted_model <- fitModel(
                    Design_out,
                    as.matrix(protein_tab[protein_name, , drop = FALSE])
                )
                res_table <- csTest(fitted_model,
                    coef = colnames(design),
                    cell_type = NULL,
                    verbose = FALSE
                )
                res_table[["SNP"]] <- SNP_ID[test_protein[i]]
                return(res_table)
            }, BPPARAM = BPPARAM) 

            res_df <- data.frame(matrix(unlist(lapply(
                seq_len(ncol(prop)),
                function(i) {
                    unlist(map(map(res, i), "p_value"))
                }
            )), ncol = ncol(prop)))
            colnames(res_df) <- colnames(prop)
            res_df$protein <- protein_name
            res_df$SNP <- unlist(map(res, "SNP"))
            return(res_df[, c("protein", "SNP", colnames(prop)), drop = FALSE])
        }
    )
    metadata(se)$TOAST_output <- Res_TOAST
    return(se)
}
