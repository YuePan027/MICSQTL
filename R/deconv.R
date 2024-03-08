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
#' data frame contained in assay,
#' a bulk transcript expression data frame (`gene_data`) contained in
#' `metadata` slot (can be rescaled using either log or MinMax transformations,
#' but need to be consistent between two bulk data). A "signature matrix" 
#' functions as a reference containing known cellular signatures 
#' (either `ref_protein` or `ref_gene` as an element in the `metadata` slot) 
#' may be necessary for certain `source`and `method` options. 
#' To ensure the reliability of the results obtained, we strongly recommend 
#' that the "signature matrix" should exclusively comprise markers that have 
#' been previously validated in the literature. 
#' @param source A character string denotes which molecular profiles to be
#' deconvoluted. The setting of `proteins`or `transcript` means single-source
#' deconvolution with source-specific signature matrix, while `cross` means
#' proteome deconvolution based on transcriptome-proteome with matched samples. 
#' @param method A character string specifies the deconvolution method to be 
#' employed. In the current version, only 'nnls' is supported for single-source 
#' deconvolution. Besides the bulk data, `ref_protein` or `ref_gene` is 
#' required for protein or gene deconvolution, respectively. 
#' For cross-source deconvolution, `JNMF` or `TCA` are valid options. 
#' If `JNMF`, an external reference containing cell counts 
#' in a similar tissue type (typically obtainable from small-scale single-cell 
#' or flow cytometry experiments) is necessary if `pinit = "rdirichlet"`; 
#' A "signature matrix" is required for other methods.
#' If `TCA`, an input of pre-estimated transcriptome 
#' proportions, denoted as `prop_gene` as an element in the `metadata` slot, 
#' is required. This input can be derived from single-source deconvolution 
#' using `nnls` included in this package, or from an external source.
#' @param use_refactor A numeric value indicate the number of proteins included
#' for proportion estimates based on refactor values. Note that `ajive_decomp`
#' with `refactor_loading = TRUE` required if this method applied. 
#' If NULL, then all proteins included in assay will be used.
#' @param Step A numeric vector indicates the step size in projected gradient 
#' descent for cell count fraction parameter and cell size parameter, 
#' respectively. Only valid if `method = JNMF`.
#' @param Eps A numeric value indicates the convergence criteria for projected 
#' gradient descent. Only valid if `method = JNMF`.
#' @param Iter A numeric value indicates the maximum iteration time for 
#' projected gradient descent. Only valid if `method = JNMF`.
#' @param cell_counts A matrix containing cell counts across multiple subjects, 
#' where subjects are represented as rows and cell types as columns. Each entry 
#' (i, j) in the matrix indicates the count of cells belonging to the ith 
#' subject and jth cell type. Only required if `method = JNMF` and 
#' `pinit = "rdirichlet"`.
#' @param pinit Accepts either a numeric matrix or a character indicating the 
#' method for initializing initial values for cellular fraction. If `pinit` is a
#' numeric matrix (pre-estimated transcriptome proportions using other methods), 
#' each row represents the cellular fraction for each sample 
#' across various cell types. The resulting cellular fraction will match the 
#' cell types defined in `pinit`. Alternatively, `pinit` can be generated using 
#' either the `rdirichlet` or `nnls` method.
#' @param ref_pnl A "signature matrix" functions as a reference containing 
#' known cellular signatures. It is optional. If provided, the initial values 
#' for purified data will be generated based on `ref_pnl`. Otherwise, 
#' the initial values for purified data will be generated using a normal 
#' distribution based on bulk data. Please note that the input signature matrix 
#' should have the same rescaling transformation as the bulk 
#' transcriptomes/proteomic.
#'
#' @return A `SummarizedExperiment`. The cell-type proportion estimates for each
#' sample are stored as elements starting with prop in the metadata slot. 
#' If `method = JNMF`, then the cellular fractions obtained from proteomics and 
#' transcriptomics are stored in the `prop` and `prop2` elements, respectively, 
#' within the metadata slot. The purified data is stored in a list with the 
#' same length as the number of subjects (the number of columns in the assay). 
#' For subject i, the purified protein expression data can be obtained by 
#' accessing `se_sim@metadata$purified[[i]][["X1"]]`, and similarly, 
#'  the purified transcript expression data can be obtained by accessing 
#' `se_sim@metadata$purified[[i]][["X2"]`].
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
                   source = "cross",
                   method = c("nnls", "JNMF", "TCA"),
                   use_refactor = c(1000, NULL),
                   Step = c(10^(-8), 10^(-6)),
                   Eps = 10^(-4),
                   Iter = 500, 
                   cell_counts = NULL,
                   pinit = "nnls", 
                   ref_pnl = NULL) {
    if (source == "protein") {
        assay(se) <- as.data.frame(assay(se))
        in_use <- intersect(
            rownames(assay(se)),
            rownames(metadata(se)$ref_protein)
        )
        if(!is.null(use_refactor)){
            in_use <- intersect(in_use, 
                                names(se@metadata$refactor[[1]])
                                [order(se@metadata$refactor[[1]], decreasing = 
                                           FALSE)[seq_len(use_refactor)]])
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
        if (!method %in% c("JNMF", "TCA", "nnls")) {
            stop("At the current version, 
                 only nnls, JNMF, and TCA are supported methods.")
        }
        if(method == "JNMF"){
            in_prot <- rownames(assay(se))
            in_rna <- rownames(metadata(se)$gene_data)
            if(!is.null(use_refactor)){
                if (use_refactor > 1500 | use_refactor < 500) {
                    warning("Selecting too few or too many features is not 
                            recommended as it can result in invalid outcomes.")
                }
                in_prot <- names(se@metadata$refactor[[1]])
                in_prot <- in_prot[
                    order(se@metadata$refactor[[1]], 
                          decreasing = FALSE)[seq_len(use_refactor)]]
            }
            n_sample <- ncol(assay(se)[in_prot, , drop = FALSE])
            if(is.null(cell_counts)){
                if(is.null(ref_pnl)){
                    K <- ncol(pinit)
                }
                if(!is.null(ref_pnl)){
                    K <- ncol(ref_pnl)
                    in_prot <- intersect(in_prot, rownames(ref_pnl))
                    in_rna <- intersect(in_rna, rownames(ref_pnl))
                }
            } else{
                K <- ncol(cell_counts)
                if(!is.null(ref_pnl)){
                    in_prot <- intersect(in_prot, rownames(ref_pnl))
                    in_rna <- intersect(in_rna, rownames(ref_pnl))
                }
            }
            n_feature1 <- length(in_prot)
            n_feature2 <- length(in_rna)
            my_res <- lapply(seq_len(n_sample), function(i){
                message("sample", i)
                if(!is.null(ref_pnl)){
                    X1 <- ref_pnl[in_prot,, drop = FALSE]
                    X2 <- ref_pnl[in_rna,, drop = FALSE]
                } else{
                    X1 <- matrix(rnorm(n = K*n_feature1, 
                                       mean = mean(as.vector(
                                           as.matrix(assay(se)[in_prot,, 
                                                               drop = FALSE]))), 
                                       sd = sd(as.vector(
                                           as.matrix(assay(
                                               se)[in_prot,, drop = FALSE ])))), 
                                 nrow = n_feature1)
                    X2 <- matrix(rnorm(n = K*n_feature2, 
                                       mean = mean(as.vector(as.matrix(
                                           metadata(se)$gene_data))), 
                                       sd = sd(as.vector(as.matrix(
                                           metadata(se)$gene_data)))),
                                 nrow = n_feature2) 
                }
                Y1 <- as.matrix(assay(se)
                                [in_prot, , drop = FALSE])[,i, drop = TRUE]
                Y2 <- as.matrix(metadata(se)$gene_data
                                [in_rna, ,drop = FALSE])[,i, drop = TRUE]
                if(is.matrix(pinit)){
                    ini_p <- pinit[i,, drop = TRUE]
                } else if(pinit == "rdirichlet"){
                    if (is.null(cell_counts)) {
                        stop("cell_counts is required for rdirichlet method.")
                    }
                    health.cell.prop.alpha <- ini_prep(cell_counts)
                    subj.var <- 0.04 #larger:larger subject-level proportion var
                    ini_p <- as.vector(rdirichlet(
                        1,(health.cell.prop.alpha*((1-subj.var)/subj.var))/
                            sum(health.cell.prop.alpha)))
                } else if (pinit == "nnls"){
                    if (is.null(ref_pnl)) {
                        stop("ref_pnl is required for nnls method.")
                    }
                    result <- nnls::nnls(ref_pnl[in_rna, , drop = FALSE],
                                         metadata(se)$gene_data[in_rna,i])   
                    ini_p <- result$x/sum(result$x)
                }
                ini_s <- list(rep(1,K), rep(1,K))
                if(!is.null(cell_counts)){
                    colnames(X1) <- colnames(cell_counts)
                    colnames(X2) <- colnames(cell_counts)
                } else{
                    colnames(X1) <- colnames(ref_pnl)
                    colnames(X2) <- colnames(ref_pnl)
                }
                rownames(X1) <- in_prot
                rownames(X2) <- in_rna
                res <- MICSQTL_optim(Y1, Y2,
                                     ini_p, ini_s,
                                     X1 = X1, X2 = X2,
                                     step_p = Step[1],
                                     step_s = Step[2],
                                     eps = Eps,
                                     iter = Iter)
                return(res[c("prop1", "prop2", "X1", "X2", 
                             "ini_p", "s1", "s2")])
            })
            prop <- do.call(rbind, lapply(my_res, 
                                          function(sub_res) sub_res$prop1))
            rownames(prop) <- colnames(assay(se)[in_prot, , drop = FALSE])
            prop2 <- do.call(rbind, lapply(my_res, 
                                           function(sub_res) sub_res$prop2))
            rownames(prop2) <- colnames(metadata(se)$gene_data)
            
            if(!is.null(cell_counts)){
                colnames(prop) <- colnames(cell_counts)
                colnames(prop2) <- colnames(cell_counts)
            } else{
                colnames(prop) <- colnames(ref_pnl)
                colnames(prop2) <- colnames(ref_pnl)
            }
            metadata(se)$prop <- prop
            metadata(se)$prop2 <- prop2
            
            prop_init <- do.call(rbind, lapply(my_res, 
                                               function(sub_res) sub_res$ini_p))
            rownames(prop_init) <- colnames(metadata(se)$gene_data)
            if(!is.null(cell_counts)){
                colnames(prop_init) <- colnames(cell_counts)
            } else{
                colnames(prop_init) <- colnames(ref_pnl)
            }
            metadata(se)$prop_init <- prop_init
            metadata(se)$purified <- my_res
        } else if(method == "TCA"){
            result <- metadata(se)$prop_gene
            ini_prop <- result[, seq_len(ncol(result))]
            in_prot <- rownames(assay(se))
            if(!is.null(use_refactor)){
                in_prot <- names(se@metadata$refactor[[1]])
                in_prot <- in_prot[
                    order(se@metadata$refactor[[1]], 
                          decreasing = FALSE)[seq_len(use_refactor)]]
            }
            prot_data <- 
                tca_res <- TCA::tca(
                    X = assay(se)[in_prot, , drop = FALSE],
                    W = ini_prop,
                    refit_W = TRUE,
                    refit_W.sparsity = nrow(assay(se)[in_prot, , drop = FALSE])
                )
            prop <- tca_res$W
            metadata(se)$prop <- prop
        }
    }
    return(se)
}