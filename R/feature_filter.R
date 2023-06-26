#' Feature filtering 
#'
#' This function returns a `SummarizedExperiment` object including SNPs used to
#' test for each protein in downstream analysis.
#'
#' This is a function developed to filter unwanted proteins or SNPs with less
#' variation among samples for downstream analysis.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression data
#' frame contained in `counts` slot.
#' Annotations on each row (protein) should be stored in rowData() with protein
#' symbol as row names
#' The first column should be a character vector indicating which chromosome
#' each protein is on.
#' A "Start" column with numeric values indicating the start position on that
#' chromosome and
#' a "Symbol" column as a unique name for each protein is also required.
#' The information from genetic variants should be stored in a P
#' (the number of SNP) by N
#' (the number of samples, should match the sample in `counts` slot) matrix
#' contained as an element
#' (`SNP_data`) in `metadata` slot.
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2)
#' for a sample at a genetic location.
#' The annotations of these SNP should be stored as an element (`anno_SNP`)
#' in `metadata` slot.
#' It should include at least the following columns: "CHROM"
#' (which chromosome the SNP is on),
#' "POS" (position of that SNP) and "ID" (a unique identifier for each SNP,
#' usually a combination of
#' chromosome and its position).
#' @param target_protein A character vector contains proteins names that will
#' be used for downstream analysis.
#' By default, all proteins in `counts` slot will be used.
#' @param target_SNP A character vector contains SNP IDs that will be used
#' for downstream analysis.
#' If not provided, all SNPs will be used for further filtering.
#' @param filter_method A character string denotes which filtering method
#' will be used to filter out unrelated SNPs.
#' If "allele", then the minor allele frequency below argument `filter_allele`
#' will be filtered out.
#' If "distance", then only cis-acting SNPs for each protein
#' (defined as SNPs on the same chromosome and
#' within 1M base pair (bp) range of that protein) will be included for
#' downstream analysis.
#' if "null", then the same SNPs will be used for each protein.
#' @param filter_allele A numeric value denotes the threshold for minor
#' allele frequency. Only works when `filter_method`
#' contains "allele".
#' @param filter_geno A numeric value denotes the threshold for minimum
#' genotype group proportion.
#' Only works when `filter_method` contains "allele".
#' @param ref_position A character string denotes the reference position
#' on protein when `filter_method` contains "distance",
#' where "TSS" refers to transcription start site, and "genebody" refers
#' to the middle point of "Start" and "End" position.
#' @param BPPARAM For applying `bplapply`.
#'
#' @return A `SummarizedExperiment`. The results after filtering will be
#' stored as an element
#' (`choose_SNP_list`) in `metadata` slot.
#' `choose_SNP_list` is a list with the length of the number of proteins
#' for downstream analysis.
#' Each element stores the index of SNPs to be tested for corresponding
#' protein.
#' The proteins with no SNPs correspond to it will be removed from the
#' returned list.
#'
#' @import SummarizedExperiment
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @importFrom S4Vectors metadata
#'
#' @export
#'
#' @examples
#' data(se)
#' target_protein <- rowData(se)[rowData(se)$Chr == 9, ][seq_len(3), "Symbol"]
#' se <- feature_filter(se,
#'     target_protein = target_protein,
#'     filter_method = c("allele", "distance"),
#'     filter_allele = 0.15, filter_geno = 0.05, ref_position = "TSS"
#' )
#'
feature_filter <- function(se,
                           target_protein = NULL,
                           target_SNP = NULL,
                           filter_method = c("allele", "distance", "null"),
                           filter_allele = 0.25,
                           filter_geno = 0.05,
                           ref_position = c("TSS", "genebody"),
                           BPPARAM = bpparam()) {
    assay(se) <- as.data.frame(assay(se))

    if (!all(colnames(assay(se)) ==
        colnames(metadata(se)$SNP_data))) {
        stop("Samples in protein_data do not match that in SNP_data")
    }

    if (!is.null(target_protein)) {
        metadata(se)$target_dat <- 
          assay(se)[target_protein, , drop = FALSE]
    } else {
        metadata(se)$target_dat <- assay(se)
        target_protein <- rownames(assay(se))
    }

    if (!is.null(target_SNP)) {
        metadata(se)$SNP_data <-
            metadata(se)$SNP_data[
                match(target_SNP, metadata(se)$anno_SNP$ID), ,
                drop = FALSE
            ]
        metadata(se)$anno_SNP <-
            metadata(se)$anno_SNP[
                match(target_SNP, metadata(se)$anno_SNP$ID), ,
                drop = FALSE
            ]
    }

    choose_SNP_list <-
        rep(list(seq_len(nrow(metadata(se)$SNP_data))),
            each = length(target_protein)
        )

    if ("allele" %in% filter_method) {
        prop <- apply(metadata(se)$SNP_data, 1, prop_allele)
        prop2 <- apply(metadata(se)$SNP_data, 1, prop_geno)
        idx <- which(pmin(prop, 1 - prop) > filter_allele & prop2 > filter_geno)
        metadata(se)$SNP_data <-
            metadata(se)$SNP_data[idx, , drop = FALSE]
        metadata(se)$anno_SNP <-
            metadata(se)$anno_SNP[idx, , drop = FALSE]
        choose_SNP_list <-
            rep(list(seq_len(nrow(metadata(se)$SNP_data))),
                each = length(target_protein)
            )
    }
    if ("distance" %in% filter_method) {
        anno_SNP <- metadata(se)$anno_SNP
        choose_SNP_list <- bplapply(seq_len(length(target_protein)),
            function(i) {
                message(
                    "Filter SNP based on distance for protein ",
                    target_protein[i], "\n"
                )
                df <-
                    rowData(se)[which(rowData(se)$Symbol ==
                        target_protein[i]), , drop = FALSE]
                if (ref_position == "TSS") {
                    df$ref_pos <- df$Start
                } else if (ref_position == "genebody") {
                    df$ref_pos <- (df$Start + df$end) / 2
                }
                idx2 <- which((abs(as.numeric(anno_SNP$POS) - df$ref_pos) <
                    1000000) & anno_SNP$CHROM == 
                      as.character(as.data.frame(df[1,1,drop=FALSE])))
                return(idx2)
            },
            BPPARAM = BPPARAM
        )
        names(choose_SNP_list) <- target_protein
        idx3 <- which(unlist(lapply(choose_SNP_list, length)) != 0)
        metadata(se)$target_dat <-
            metadata(se)$target_dat[idx3, , drop = FALSE]
        choose_SNP_list <- choose_SNP_list[idx3]
    }

    metadata(se)$choose_SNP_list <- choose_SNP_list
    return(se)
}
