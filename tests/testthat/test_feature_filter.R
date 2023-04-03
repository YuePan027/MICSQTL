se <- SummarizedExperiment(assays = list(counts = mcQTL::protein_data),
                           rowData = mcQTL::anno_protein)
se@metadata <- list(SNP_data = mcQTL::SNP_data[,rev(seq_len(ncol(mcQTL::SNP_data)))], anno_SNP = mcQTL::anno_SNP)
test_that("feature filter function gives error for unmatched samples", {
    expect_error(feature_filter(se, target_protein = c("Protein_5"),
                                filter_method = c("allele")), 
                 "Samples in protein_data do not match that in SNP_data.")
})


se@metadata <- list(SNP_data = mcQTL::SNP_data, anno_SNP = mcQTL::anno_SNP[1:100,])
test_that("feature filter function gives error for unmatched annotation", {
    expect_error(feature_filter(se, target_protein = c("Protein_5"),
                                filter_method = c("allele")), 
                 "SNPs contained in annotation data frame `anno_SNP` must match the SNPs in `SNP_data`.")
})