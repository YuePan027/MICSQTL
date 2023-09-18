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
#' `metadata` slot. A "signature matrix" functions as a reference containing 
#' known cellular signatures (either `ref_protein` or `ref_gene` as an element 
#' in the `metadata` slot) may be necessary for certain `source`and `method` 
#' options. To ensure the reliability of the results obtained, we strongly 
#' recommend that the "signature matrix" should exclusively comprise markers 
#' that have been previously validated in the literature. 
#' @param source A character string denotes which molecular profiles to be
#' deconvoluted. The setting of `proteins`
#' or `transcript` means single-source
#' deconvolution with source-specific signature matrix, while `cross` means
#' proteome deconvolution based on transcriptome-proteome with matched samples. 
#' @param method A character string specifies the deconvolution method to be 
#' employed. In the current version, only `nnls` is supported for 
#' single-source deconvolution. 
#' For cross-source deconvolution, 'Joint' or 'TCA' are valid options. 
#' If 'Joint', an external reference containing cell counts 
#' in a similar tissue type (typically obtainable from small-scale single-cell 
#' or flow cytometry experiments) is necessary. 
#' If 'TCA', an input of pre-estimated transcriptome 
#' proportions, denoted as `prop_gene` as an element in the `metadata` slot, 
#' is required. This input can be derived from single-source deconvolution 
#' using `nnls` included in this package, or from an external source.
#' @param use_refactor A numeric value indicate the number of proteins included
#' for proportion estimates based on refactor values. Note that `ajive_decomp`
#' with `refactor_loading = TRUE` required if this method applied. 
#' If NULL, then all proteins included in assay will be used.
#' @param Step A numeric vector indicates the step size in projected gradient 
#' descent for cell count fraction parameter and cell size parameter, 
#' respectively. Only valid if `method = Joint`.
#' @param Eps A numeric value indicates the convergence criteria for projected 
#' gradient descent. Only valid if `method = Joint`.
#' @param Iter A numeric value indicates the maximum iteration time for 
#' projected gradient descent. Only valid if `method = Joint`.
#' @param cell_counts A matrix containing cell counts across multiple subjects, 
#' where subjects are represented as rows and cell types as columns. Each entry 
#' (i, j) in the matrix indicates the count of cells belonging to the ith 
#' subject and jth cell type. Only required if `method = Joint`.
#'
#' @return A `SummarizedExperiment`. The cell-type proportion estimates for
#' each sample will be stored as an
#' element start with `prop` in `metadata` slot. If `method = Joint`, then the 
#' cellular fractions obtained from proteomics and transcriptomics are stored 
#' in the `prop` and `prop2`, respectively, within the `metadata` slot.
#' 
#'
#' @import SummarizedExperiment
#' @importFrom magrittr "%>%"
#' @importFrom nnls nnls
#' @importFrom TCA tca
#' @importFrom S4Vectors metadata
#' @importFrom dirmult rdirichlet dirmult
#' @importFrom stats rnorm sd
#' 
#'
#' @export
#'
#' @examples
#' data(se)
#' se <- deconv(se, source = "protein", method = "nnls", use_refactor = NULL)
#'
deconv <- function(se,
                   source = "protein",
                   method = c("nnls", "Joint", "TCA"),
                   use_refactor = c(1000, NULL),
                   Step = c(10^(-6), 10^(-6)),
                   Eps = 10^(-5),
                   Iter = 500, cell_counts) {
    if (source == "protein") {
        assay(se) <- as.data.frame(assay(se))
        in_use <- intersect(
            rownames(assay(se)),
            rownames(metadata(se)$ref_protein)
        )
        if(!is.null(use_refactor)){
            in_use <- intersect(in_use, 
                                names(se@metadata$refactor[[1]])
                                [order(se@metadata$refactor[[1]],
                                       decreasing = FALSE)[1:use_refactor]])
        }
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
        metadata(se)$prop <- prop
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
        metadata(se)$prop <- prop
    }
    if (source == "cross") {
        if(method == "Joint"){
            in_prot <- rownames(assay(se))
            if(!is.null(use_refactor)){
                if (use_refactor > 1500 | use_refactor < 500) {
                    warning("Selecting too few or too many features is not 
                            recommended as it can result in invalid outcomes.")
                }
                in_prot <- names(se@metadata$refactor[[1]])
                in_prot <- in_prot[order(se@metadata$refactor[[1]], 
                                         decreasing = FALSE)[1:use_refactor]]
            }
            n_sample <- ncol(assay(se)[in_prot, ])
            K <- ncol(cell_counts)
            n_feature1 <- nrow(assay(se)[in_prot, ])
            n_feature2 <- nrow(metadata(se)$gene_data)
            my_res <- lapply(1:n_sample, function(i){
                set.seed(1234)
                X1 <- matrix(rnorm(n = K*n_feature1, 
                                   mean = mean(as.vector(
                                       as.matrix(assay(se)[in_prot, ]))), 
                                   sd = sd(as.vector(
                                       as.matrix(assay(se)[in_prot, ])))), 
                             nrow = n_feature1)
                set.seed(1234)
                X2 <- matrix(rnorm(n = K*n_feature2, 
                                   mean = mean(as.vector(as.matrix(
                                       metadata(se)$gene_data))), 
                                   sd = sd(as.vector(as.matrix(
                                       metadata(se)$gene_data)))),
                             nrow = n_feature2)
                Y1 <- as.matrix(assay(se)[in_prot, ])[,i]
                Y2 <- as.matrix(metadata(se)$gene_data)[,i]
                health.cell.prop.alpha <- ini_prep(cell_counts)
                subj.var <- 0.04 # larger -> larger subject-level proportion var
                set.seed(i)
                ini_p <- as.vector(rdirichlet(
                    1,(health.cell.prop.alpha*((1-subj.var)/subj.var))/
                        sum(health.cell.prop.alpha)))
                ini_s <- list(rep(1,K), rep(1,K))
                res <- MICSQTL_optim(Y1, Y2,
                                     ini_p, ini_s,
                                     X1 = X1, X2 = X2,
                                     step_p = Step[1],
                                     step_s = Step[2],
                                     eps = Eps,
                                     iter = Iter)
                return(res[c("prop1", "prop2")])
                
            })
            prop <- do.call(rbind, lapply(my_res, 
                                          function(sub_res) sub_res$prop1))
            rownames(prop) <- colnames(assay(se)[in_prot, ])
            colnames(prop) <- colnames(cell_counts)
            metadata(se)$prop <- prop
            
            prop2 <- do.call(rbind, lapply(my_res, 
                                          function(sub_res) sub_res$prop2))
            rownames(prop2) <- colnames(metadata(se)$gene_data)
            colnames(prop2) <- colnames(cell_counts)
            metadata(se)$prop2 <- prop2
        } else if(method == "TCA"){
            result <- metadata(se)$prop_gene
            ini_prop <- result[, seq_len(ncol(result))]
            in_prot <- rownames(assay(se))
            if(!is.null(use_refactor)){
                in_prot <- names(se@metadata$refactor[[1]])
                in_prot <- in_prot[order(se@metadata$refactor[[1]], 
                                         decreasing = FALSE)[1:use_refactor]]
            }
            prot_data <- 
                tca_res <- TCA::tca(
                    X = assay(se)[in_prot, ],
                    W = ini_prop,
                    refit_W = TRUE,
                    refit_W.sparsity = nrow(assay(se)[in_prot, ])
                )
            prop <- tca_res$W
            metadata(se)$prop <- prop
        }
    }

    return(se)
}
