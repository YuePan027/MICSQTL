% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCA_deconv.R
\name{TCA_deconv}
\alias{TCA_deconv}
\title{Estimation of cellular expression in high-throughput data
from heterogeneous tissues}
\usage{
TCA_deconv(se, test = "array", prop)
}
\arguments{
\item{se}{A `SummarizedExperiment` object with bulk expression data frame
contained in `counts` slot, and
cell-type proportion estimates for each sample contained as an element
(`prop`) in `metadata` slot.}

\item{test}{A character string indicate which data set to estimate the
cellular expression.
If not "array", it must contained in `metadata` slot.}

\item{prop}{A matrix of cellular composition. It has to be provided
unless `test` is set to be "array".}
}
\value{
A `SummarizedExperiment`. The results after `TCA` deconvolution
will be stored
as an element (`TCA_deconv`) in `metadata` slot.
It is a list with the length of the number of cell types (as in `prop`
in `metadata` slot).
Each element stores the deconvoluted protein expression matrix per bulk
sample within that cell type.
}
\description{
This function returns a `SummarizedExperiment` object including
deconvoluted cellular sample-wise expression.
}
\details{
This is a function developed to implement deconvolution for
cell-type-specific expression per bulk sample using `TCA`.
}
\examples{
data(se)
se <- deconv(se, source = "protein", method = "nnls")
se <- feature_filter(se,
    target_protein = c("ABCA1", "ABCA2"),
    filter_method = c("allele", "distance"), filter_allele = 0.15,
    filter_geno = 0.05, ref_position = "TSS"
)
se <- TCA_deconv(se, prop = metadata(se)$prop)

}
