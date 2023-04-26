data(se)
set.seed(1234)
se@metadata$SNP_data <- 
    se@metadata$SNP_data[, sample(ncol(se@metadata$SNP_data))]

test_that("feature filter function gives error for unmatched samples", {
    expect_error(
        feature_filter(se,
            target_protein = c("ABCA1"),
            filter_method = c("allele")
        ),
        "Samples in protein_data do not match that in SNP_data"
    )
})


data(se)
se@metadata$anno_SNP <- se@metadata$anno_SNP[1:100, ]
test_that("feature filter function gives error for unmatched annotation", {
    expect_error(
        feature_filter(se,
            target_protein = c("ABCA1"),
            filter_method = c("allele")
        ),
        "SNPs contained in annotation data frame `anno_SNP` must match the SNPs
        in `SNP_data`"
    )
})
