#' Cell-type-specific differential expression (csQTL)
#'
#' This function returns a `SummarizedExperiment` object including csQTL proteins based on samples' genotype.
#'
#' This is a function developed to implement cell-type-specific differential expression using `TOAST`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression data frame contained in `counts` slot.
#' The information from genetic variants should be stored in a P (the number of SNP) by N (the number of samples, should match the sample in `counts` slot) matrix contained as an element (`SNP_data`) in `metadata` slot.
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2) for a sample at a genetic location.
#' The annotations of these SNP should be stored as an element (`anno_SNP`) in `metadata` slot.
#' It should include at least the following columns: "CHROM" (which chromosome the SNP is on),
#' "POS" (position of that SNP) and "ID" (a unique identifier for each SNP, usually a combination of chromosome and its position).
#' The information on cellular composition is required and stored as `prop` in `metadata` slot.
#' It is an N (the number of samples, should match the sample in `counts` slot) by K (the number of cell types) matrix.
#' This can be obtained by running `deconv()` before any filtering steps, or use any source data directly.
#' @param BPPARAM For applying `bplapply`.
#'
#' @return A `SummarizedExperiment`. The csQTL results will be stored as an element (`TOAST_output`) in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @import TOAST
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#'
#' @export
#'
#'
csQTL <- function(se, BPPARAM = bpparam()){

  assay(se) <- as.data.frame(assay(se))

  prop <- se@metadata$prop
  protein_tab <- as.data.frame(se@metadata$target_dat)
  SNP_dat <- se@metadata$SNP_data
  SNP_ID <- se@metadata$anno_SNP$ID

  Res_TOAST <- lapply(1:length(se@metadata$choose_SNP_list), function(x){
    test_protein <- se@metadata$choose_SNP_list[[x]]
    protein_name <- names(se@metadata$choose_SNP_list)[x]
    cat("csQTL test for protein", protein_name, "\n")

    res <- bplapply(1:length(test_protein), function(i){
      design <- data.frame(factor(SNP_dat[test_protein[i],]))#data.frame(factor(t(SNP_dat)[, test_protein[i]]))
      Design_out <- makeDesign(design, prop)
      fitted_model <- fitModel(Design_out,
                               as.matrix(protein_tab[protein_name,]))
      res_table <- csTest(fitted_model,
                          coef = colnames(design),
                          cell_type = NULL)
      res_table[["SNP"]] <- SNP_ID[test_protein[i]]
      return(res_table)
    },BPPARAM=BPPARAM) %>%
      suppressWarnings() %>%
      suppressMessages()

    res_df <- data.frame(matrix(unlist(lapply(1:ncol(prop), function(i){
      unlist(map(map(res,i), "p_value"))
    })), ncol = ncol(prop)))
    colnames(res_df) <- colnames(prop)
    res_df$protein <- protein_name
    res_df$SNP <- unlist(purrr::map(res,"SNP"))
    # res_df$min_p <- apply(res_df[, 1:ncol(se@metadata$prop)], 1, FUN = min)
    # idx4 <- match(res_df$SNP[which(res_df$min_p < filter_TOAST)], anno_SNP$ID)
    return(res_df[, c("protein", "SNP", colnames(prop))])
  })
  #se@metadata$choose_SNP_list <- map(Res_TOAST, 2)
  #res_TOAST <- map(Res_TOAST, 1)
  se@metadata$TOAST_output <- Res_TOAST
  return(se)
}
