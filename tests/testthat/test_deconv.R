N <- 50 # 50 samples
P <- 100 # 500 features
data <- matrix(rnorm(N*P, N, P), ncol = N)
colnames(data) <- paste("subject", 1:ncol(data), sep = "_")
rownames(data) <- paste("feature", 1:nrow(data), sep = "_")
se <- SummarizedExperiment(assays = list(counts = data))
metadata(se) <- list(sig_matrix = MICSQTL::ref_data)


test_that("deconv function gives error for invalid inputs", {
    expect_error(deconv(se, "newmethod"), "Only 'cibersort' or 'nnls' is valid in current version.")
    expect_error(deconv(se, "cibersort"), "None of the feaures in 'signature matrix' exist in bulk expression data.")
})
