#' Estimation of cellular composition in high-throughput data from
#' heterogeneous tissues
#'
#' This function returns a `SummarizedExperiment` object including cell-type
#' proportion estimates for each sample.
#'
#' This is a function developed to implement cell-type proportion deconvolution
#' using either `CIBERSORT` or `nnls`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression
#' data frame contained in `counts` slot,
#' a bulk transcript expression data frame (`gene_data`) contained in
#' `metadata` slot, and
#' a "signature matrix" which serves as a reference of known cellular signatures
#'  contained as an element start with
#' `sig` (`sig_protein` or `sig_gene` depends on `source`) in `metadata` slot.
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
#' In the current version, only
#' `CIBERSORT` or `nnls` is supported.
#'
#' @return A `SummarizedExperiment`. The cell-type proportion estimates for
#' each sample will be stored as an
#' element start with `prop` (depends on `source`) in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @importFrom magrittr "%>%"
#' @importFrom methods slot
#'
#' @export
#'
#' @examples
#' data(se)
#' se <- deconv(se, source = "protein", method = "cibersort")
#'
deconv <- function(se,
                   source = "protein",
                   method = "cibersort") {
    if (source == "protein") {
        assay(se) <- as.data.frame(assay(se))
        in_use <- intersect(
            rownames(assay(se)),
            rownames(methods::slot(se, "metadata")$sig_protein)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        protein_sub <- as.data.frame(assay(se)[in_use, , drop = FALSE])
        sig_protein <- 
            methods::slot(se, "metadata")$sig_protein[in_use, , drop = FALSE]
        if (method == "cibersort") {
            result <- CIBERSORT(
                sig_matrix = sig_protein,
                mixture_file = as.data.frame(protein_sub),
                perm = 0, QN = TRUE,
                absolute = FALSE, abs_method = "sig.score"
            )
            prop <- result[, seq_len(ncol(sig_protein))]
        }
        if (method == "nnls") {
            decon_nnls <- apply(
                protein_sub, 2,
                function(y) nnls::nnls(as.matrix(sig_protein), y)$x
            ) %>% t()
            prop <- decon_nnls / rowSums(decon_nnls)
            rownames(prop) <- colnames(protein_sub)
            colnames(prop) <- colnames(sig_protein)
        }
        methods::slot(se, "metadata")$prop <- prop
    }
    if (source == "transcript") {
        in_use <- intersect(
            rownames(methods::slot(se, "metadata")$gene_data),
            rownames(methods::slot(se, "metadata")$sig_gene)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        gene_sub <- 
            as.data.frame(methods::slot(se, "metadata")$gene_data[in_use, ,
            drop = FALSE
        ])
        sig_gene <- 
            methods::slot(se, "metadata")$sig_gene[in_use, , drop = FALSE]
        if (method == "cibersort") {
            result <- CIBERSORT(
                sig_matrix = sig_gene,
                mixture_file = as.data.frame(gene_sub),
                perm = 0, QN = TRUE,
                absolute = FALSE, abs_method = "sig.score"
            )
            prop <- result[, seq_len(ncol(sig_gene))]
        }
        if (method == "nnls") {
            decon_nnls <- apply(gene_sub, 2, function(y) {
                nnls::nnls(as.matrix(sig_gene), y)$x
            }) %>% t()
            prop <- decon_nnls / rowSums(decon_nnls)
            rownames(prop) <- colnames(gene_sub)
            colnames(prop) <- colnames(sig_gene)
        }
        methods::slot(se, "metadata")$prop <- prop
    }
    if (source == "cross") {
        in_use <- intersect(
            rownames(methods::slot(se, "metadata")$gene_data),
            rownames(methods::slot(se, "metadata")$sig_gene)
        )
        if (any(length(in_use) == 0)) {
            stop("None of the feaures in 'signature matrix' exist in bulk
                 expression data.")
        }
        gene_sub <- 
            as.data.frame(methods::slot(se, "metadata")$gene_data[in_use, ,
            drop = FALSE
        ])
        sig_gene <- 
            methods::slot(se, "metadata")$sig_gene[in_use, , drop = FALSE]
        result <- CIBERSORT(
            sig_matrix = sig_gene,
            mixture_file = as.data.frame(gene_sub),
            perm = 0, QN = TRUE,
            absolute = FALSE, abs_method = "sig.score"
        )
        ini_prop <- result[, seq_len(ncol(sig_gene))]
        mrk_prot <- intersect(rownames(assay(se)), in_use)
        tca_res <- TCA::tca(
            X = assay(se)[mrk_prot, ],
            W = ini_prop,
            refit_W = TRUE,
            refit_W.sparsity = length(mrk_prot)
        )
        prop <- tca_res$W
        methods::slot(se, "metadata")$prop <- prop
    }
    return(se)
}
