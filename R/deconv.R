#' Estimation of cellular composition in high-throughput data from
#' heterogeneous tissues
#'
#' This function returns a `SummarizedExperiment` object including cell-type
#' proportion estimates for each sample.
#'
#' This is a function developed to implement cell-type proportion deconvolution
#' using either single or cross sources.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression
#' data frame contained in `counts` slot,
#' a bulk transcript expression data frame (`gene_data`) contained in
#' `metadata` slot, and
#' a "signature matrix" which serves as a reference of known cellular signatures
#'  contained as an element start with
#' `sig` (`ref_protein` or `ref_gene` depends on `source`) in `metadata` slot.
#' Note that the 'signature matrix' should only include markers that have been
#' demonstrated to be useful in previous
#' literature to ensure reliable results.
#' @param source A character string denotes which molecular profiles to be
#' deconvoluted. The setting of `proteins`
#' or `transcript` means single-source
#' deconvolution with source-specific signature matrix, while 'cross' means
#' proteome deconvolution based on
#' pre-estimated transcriptome proportion.
#' @param method A character string denotes which deconvolution method to use.
#' In the current version, only `nnls` is supported for single source.
#' @param iter If TRUE, iterate estimated proportion until converged.
#' @param Max_iter The maximum number of iterations if not converged.
#' @param Diff_max The convergence metric for each iteration, 
#' defined as the maximum difference from the previous iteration.
#'
#' @return A `SummarizedExperiment`. The cell-type proportion estimates for
#' each sample will be stored as an
#' element start with `prop` in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @importFrom magrittr "%>%"
#' @importFrom nnls nnls
#' @importFrom TCA tca
#' @importFrom S4Vectors metadata
#' 
#'
#' @export
#'
#' @examples
#' data(se)
#' se <- deconv(se, source = "protein", method = "nnls")
#'
deconv <- function(se,
                   source = "protein",
                   method = "nnls",
                   iter = FALSE,
                   Max_iter = 20, Diff_max = 0.015) {
    if (source == "protein") {
        assay(se) <- as.data.frame(assay(se))
        in_use <- intersect(
            rownames(assay(se)),
            rownames(metadata(se)$ref_protein)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        protein_sub <- as.data.frame(assay(se)[in_use, , drop = FALSE])
        ref_protein <-
            metadata(se)$ref_protein[in_use, , drop = FALSE]
        
        if (method == "nnls") {
            decon_nnls <- apply(
                protein_sub, 2,
                function(y) nnls(as.matrix(ref_protein), y)$x
            ) %>% t()
            prop <- decon_nnls / rowSums(decon_nnls)
            rownames(prop) <- colnames(protein_sub)
            colnames(prop) <- colnames(ref_protein)
        }
    }
    if (source == "transcript") {
        in_use <- intersect(
            rownames(metadata(se)$gene_data),
            rownames(metadata(se)$ref_gene)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        gene_sub <-
            as.data.frame(metadata(se)$gene_data[in_use, ,
                drop = FALSE
            ])
        ref_gene <-
            metadata(se)$ref_gene[in_use, , drop = FALSE]
    
        if (method == "nnls") {
            decon_nnls <- apply(gene_sub, 2, function(y) {
                nnls::nnls(as.matrix(ref_gene), y)$x
            }) %>% t()
            prop <- decon_nnls / rowSums(decon_nnls)
            rownames(prop) <- colnames(gene_sub)
            colnames(prop) <- colnames(ref_gene)
        }
    }
    if (source == "cross") {
        in_use <- intersect(
            rownames(metadata(se)$gene_data),
            rownames(metadata(se)$ref_gene)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        gene_sub <-
            as.data.frame(metadata(se)$gene_data[in_use, ,
                drop = FALSE
            ])
        ref_gene <-
            metadata(se)$ref_gene[in_use, , drop = FALSE]
        result <- metadata(se)$prop_gene
        ini_prop <- result[, seq_len(ncol(ref_gene))]
        mrk_prot <- intersect(rownames(assay(se)), in_use)
        tca_res <- tca(
            X = assay(se)[mrk_prot, , drop = FALSE],
            W = ini_prop,
            refit_W = TRUE,
            refit_W.sparsity = length(mrk_prot)
        )
        prop <- tca_res$W
    }
    metadata(se)$prop <- prop
    if(iter){
        TCA_iter <- function(ini_prop, max_iter, diff_max){
          iter_idx <- 1
          while (iter_idx < max_iter){
            tca_res <- TCA::tca(
              X = assay(se)[mrk_prot, , drop = FALSE],
              W = ini_prop,
              refit_W = TRUE,
              refit_W.sparsity = length(mrk_prot)
            )
            diff <- abs(ini_prop - tca_res$W)
            if(max(diff) < diff_max){
              return(tca_res$W)
            } else{
              iter_idx <- iter_idx + 1
              ini_prop <- tca_res$W
            }
          }
        }
        Res <- TCA_iter(prop, Max_iter, Diff_max)
        metadata(se)$prop <- Res
    }
    return(se)
}
