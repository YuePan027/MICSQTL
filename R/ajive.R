#' Integrative analysis for modes of joint variation
#'
#' This function returns a `SummarizedExperiment` object including results from
#' AJIVE (Angle based Joint and Individual Variation Explained),
#' an integrative analysis tool, and common normalized scored based on it.
#' This implemented the AJIVE algorithm 
#' (see R version of AJIVE package (https://github.com/idc9/r_jive)
#' or detailed manuscript at https://arxiv.org/pdf/1704.02060.pdf)
#' an integrative analysis method.
#' It is useful when there are multiple data matrices measured on the same set
#' of samples.
#' It decomposes each data matrix as three parts: (1) Joint variation across
#' data types
#' (2) Individual structured variation for each data type and (3) Residual
#' noise.
#' Common normalized scores are one of the desirable output to explore the
#' joint behavior that is shared by different data sources.
#'
#'
#' @param se A `SummarizedExperiment` object with bulk expression data frame
#' contained in `counts` slot.
#' In addition, data measured on the same set of samples in `metadata` slot is
#' required.
#' @param ini_rank A vector with each element corresponds to a initial signal
#' rank for AJIVE decomposition.
#' Please refer to the original paper (Feng, Qing, et al. "Angle-based joint
#' and individual variation explained." Journal of multivariate analysis 166
#' (2018): 241-265.)
#' on the choice of initial ranks.
#' @param test A character string indicate which data are used as a secondary
#' data block measured on the same set of samples.
#' @param use_marker If TRUE, only markers contained in `ref_gene` are used.
#' @param level A character string indicate if the integrative analysis should
#' be done at cell type specific level and
#' which cell type should be used. By default, the integrative analysis is done
#' at bulk level.
#' @param plot If TRUE, visualization on common normalized scores across 
#' different data sources will be stored as an element (`cns_plot`) in 
#' `metadata` slot.
#' @param score A character of variable name indicating which common
#' normalized score is used for boxplot and ridge plot (valid if plot = TRUE).
#' @param group_var A character of variable name indicating which variable is
#' used as the group variable to compare common normalized scores 
#' (valid if plot = TRUE).
#' @param scatter A logical value indicating whether to make scatter plot or
#' not. Only valid when the joint rank is at least two (valid if plot = TRUE).
#' @param scatter_x A character of variable name indicating which common
#' normalized scores on horizontal axis (valid if plot = TRUE).
#' @param scatter_y A character of variable name indicating which common
#' normalized scores on vertical axis (valid if plot = TRUE).
#' @param refactor_loading A logical value indicating whether to output the 
#' refactor joint loadings for each data source. A smaller value indicates 
#' features with the highest variance explained within the joint space.
#'
#' @return A `SummarizedExperiment`. The results from AJIVE will be stored as
#' an element (`ajive_res`) in `metadata` slot.
#' The generated common normalized scores will be stored as an element (`cns`)
#' in `metadata` slot. The visualization of the common normalized scores will 
#' be stored as an element (`cns_plot`) in `metadata` slot.
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @import ggridges
#' @importFrom ggpubr ggarrange
#'
#' @export
#'
#' @examples
#' data(se)
#' metadata(se)$gene_data <- metadata(se)$gene_data[seq_len(100),] #reduce time
#' se <- ajive_decomp(se, use_marker = TRUE)
#'
ajive_decomp <- function(se, ini_rank = c(20, 20), test = "gene_data",
                         use_marker = FALSE,
                         level = "bulk",
                         plot = FALSE,
                         score = "cns_1", group_var = "disease",
                         scatter = FALSE, scatter_x, scatter_y,
                         refactor_loading = FALSE) {
    if(is.null(metadata(se)$ajive_res)){
        if (level == "bulk") {
            dat1 <- assay(se)
            if (use_marker) {
                in_use <- intersect(
                    rownames(metadata(se)$gene_data),
                    rownames(metadata(se)$ref_gene)
                )
                dat2 <- metadata(se)[[test]][in_use, , drop = FALSE]
            } else {
                dat2 <- metadata(se)[[test]]
            }
        } else {
            dat1 <- metadata(se)$TCA_deconv[[level]]
            dat2 <- metadata(se)$TCA_deconv2[[level]]
        }
        
        if (!all(colnames(dat1) == colnames(dat2))) {
            stop("Samples in the first data do not match that 
                 in the second data")
        }
        
        blocks_test <- list(
            scale(t(dat1), center = TRUE, scale = FALSE),
            scale(t(dat2), center = TRUE, scale = FALSE)
        )
        ajive_res <- ajive(blocks_test, ini_rank)
        metadata(se)$ajive_res <- ajive_res
        cns <- get_common_normalized_scores(ajive_res)
        metadata(se)$cns <- cns
        joint_rank <- ajive_res$joint_rank
        metadata(se)$joint_rank <- joint_rank
    }
    if(plot){
        cns <- metadata(se)$cns
        colnames(cns) <- paste("cns", seq_len(ncol(cns)), sep = "_")
        df <- data.frame(cns, metadata(se)$meta)
        
        p1 <- ggplot(
            df,
            aes(
                x = .data[[group_var]], y = .data[[score]],
                fill = .data[[group_var]]
            )
        ) +
            geom_point(
                position = position_jitterdodge(
                    jitter.width = 0.1,
                    dodge.width = 0.7
                ),
                aes(fill = .data[[group_var]], color = .data[[group_var]]),
                pch = 21, alpha = 0.5
            ) +
            geom_boxplot(lwd = 0.7, outlier.shape = NA) +
            theme_classic() +
            ylab(score)
        
        p2 <- ggplot(
            df,
            aes(
                x = .data[[score]], y = .data[[group_var]],
                fill = .data[[group_var]]
            )
        ) +
            geom_density_ridges(
                aes(
                    point_color = .data[[group_var]],
                    point_fill = .data[[group_var]],
                    point_shape = .data[[group_var]]
                ),
                alpha = .2, point_alpha = 1, jittered_points = TRUE
            ) +
            scale_point_color_hue(l = 40) +
            theme_classic() +
            xlab(score)
        
        if (ncol(cns) < 2 & scatter) {
            message("The scatter plot is not supported since joint rank
                is smaller than 2")
        }
        
        if (ncol(cns) < 2 | !scatter) {
            plot <- ggarrange(p1, p2,
                             ncol = 2, nrow = 1,
                             common.legend = TRUE
            )
        } else if (scatter) {
            p3 <- ggplot(
                df,
                aes(
                    x = .data[[scatter_x]], y = .data[[scatter_y]],
                    color = .data[[group_var]]
                )
            ) +
                geom_point() +
                scale_point_color_hue(l = 40) +
                theme_classic()
            plot <- ggarrange(p1, p2, p3,
                             ncol = 2, nrow = 2,
                             common.legend = TRUE
            )
        }
        metadata(se)$cns_plot <- plot
    }
    if(refactor_loading){
        refactor <- lapply(seq_len(2), function(i){
            loadings <- get_block_loadings(metadata(se)$ajive_res, i, "joint")
            if(i == 1){
                data <- as.matrix(assay(se))
            } else{
                data <- as.matrix(metadata(se)$gene_dat)
            }
            protein_app <- loadings %*%
                t(loadings) %*% data  
            row_dist1 <- rowSums((data- protein_app)^2)
            names(row_dist1) <- rownames(data)
            return(row_dist1)
        })
        metadata(se)$refactor <- refactor
    }
    return(se)
}
