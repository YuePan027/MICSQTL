## Allele proportion
##
## This function calculates the allele frequency.
##
prop_allele <- function(x) {
    prop <- (sum(x == 0) * 2 + sum(x == 1)) / (2 * length(x))
    return(prop)
}


## Genotype proportion

prop_geno <- function(x) {
    prop <- min(table(x) / length(x))
    return(prop)
}


cross_prop <- function(se, ini_prop, mrk_prot) {
    if (!all(colnames(assay(se)) == rownames(ini_prop))) {
        stop("Samples in protein_data do not match that in initial proportion")
    }
    mrk_prot <- intersect(rownames(assay(se)), mrk_prot)
    tca_res <- TCA::tca(
        X = assay(se)[mrk_prot, ],
        W = ini_prop,
        refit_W = TRUE,
        refit_W.sparsity = length(mrk_prot)
    )
    cross_prop <- tca_res$W
    methods::slot(se, "metadata")$cross_prop <- cross_prop
    return(se)
}


CoreAlg <- function(X, y, absolute, abs_method) {
    # try different values of nu
    svn_itor <- 3

    res <- function(i) {
        if (i == 1) {
            nus <- 0.25
        }
        if (i == 2) {
            nus <- 0.5
        }
        if (i == 3) {
            nus <- 0.75
        }
        model <- svm(X, y,
            type = "nu-regression", kernel = "linear",
            nu = nus, scale = FALSE
        )
        model
    }

    # if (Sys.info()["sysname"] == "Windows") {
    out <- mclapply(seq_len(svn_itor), res, mc.cores = 1)
    # } else {
    # out <- mclapply(seq_len(svn_itor), res, mc.cores = svn_itor)
    # }

    nusvm <- rep(0, svn_itor)
    corrv <- rep(0, svn_itor)

    # do cibersort
    t <- 1
    while (t <= svn_itor) {
        weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights < 0)] <- 0
        w <- weights / sum(weights)
        u <- sweep(X, MARGIN = 2, w, "*")
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- stats::cor(k, y)
        t <- t + 1
    }

    # pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]

    # get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q < 0)] <- 0
    if (!absolute || abs_method == "sig.score") w <- (q / sum(q))
    # relative space (returns fractions)
    if (absolute && abs_method == "no.sumto1") w <- q
    # absolute space (returns scores)

    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]

    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}

# do permutations

doPerm <- function(perm, X, Y, absolute, abs_method) {
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()

    while (itor <= perm) {
        # print(itor)

        # random mixture
        yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])

        # standardize mixture
        yr <- (yr - mean(yr)) / stats::sd(yr)

        # run CIBERSORT core algorithm
        result <- CoreAlg(X, yr, absolute, abs_method)

        mix_r <- result$mix_r

        # store correlation
        if (itor == 1) {
            dist <- mix_r
        } else {
            dist <- rbind(dist, mix_r)
        }

        itor <- itor + 1
    }
    newList <- list("dist" = dist)
}




#' CIBERSORT is an analytical tool developed by Newman et al. to provide an
#' estimation of the abundances of member cell types in a mixed cell population,
#'  using gene expression data.

#' @param sig_matrix  Cell type GEP barcode matrix: row 1 = sample labels;
#' column 1 = gene symbols; no missing values; default =LM22.txt download
#' from CIBERSORT (https://cibersort.stanford.edu/runcibersort.php)
#' @param mixture_file  GEP matrix: row 1 = sample labels; column
#' 1 = gene symbols; no missing values
#' @param perm Set permutations for statistical analysis
#' (at least 100 permutations recommended).
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param absolute  Run CIBERSORT in absolute mode (default = FALSE)
#' - note that cell subsets will be scaled by their absolute levels and will
#' not be represented as fractions
#' (to derive the default output, normalize absolute levels such that they
#' sum to 1 for each mixture sample)
#' - the sum of all cell subsets in each mixture sample will be added to the
#' ouput ('Absolute score').
#' If LM22 is used, this score will capture total immune content.
#' @param abs_method  if absolute is set to TRUE, choose method: 'no.sumto1'
#' or 'sig.score'
#' - sig.score = for each mixture sample, define S as the median expression
#' level of all genes in the signature matrix divided by the median expression
#' level of all genes in the mixture. Multiple cell subset fractions by S.
#' - no.sumto1 = remove sum to 1 constraint
#' @author Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' @return cibersrot with immune cell fractions
#'
#' @import e1071
#' @import parallel
#' @import preprocessCore
#' @import tidyverse
#' @importFrom glue glue


CIBERSORT <- function(sig_matrix, mixture_file, perm, QN = TRUE, absolute,
                      abs_method = "sig.score") {
    if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0) {
        stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
    }

    if (absolute && abs_method != "no.sumto1" && abs_method != "sig.score") {
        stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
    }

    # read in data
    X <- sig_matrix

    Y <- tibble::rownames_to_column(mixture_file, var = "symbol")
    # to prevent crashing on duplicated gene symbols,
    # add unique numbers to identical names
    dups <- dim(Y)[1] - length(unique(Y[, 1]))
    if (dups > 0) {
        warning(glue::glue(dups, " duplicated gene symbol(s)
                           found in mixture file!", sep = ""))
        rownames(Y) <- make.names(Y[, 1], unique = TRUE)
    } else {
        rownames(Y) <- Y[, 1]
    }
    # Y <- Y[,-1]
    ###################################
    X <- data.matrix(X)
    Y <- data.matrix(Y)
    colname <- colnames(Y)[-1]
    Y <- as.matrix(Y[, -1], ncol = ncol(Y) - 1)
    colnames(Y) <- colname

    # order
    X <- X[order(rownames(X)), ]
    Y <- as.matrix(Y[order(rownames(Y)), ])
    colnames(Y) <- colname

    P <- perm # number of permutations

    # anti-log if max < 50 in mixture file
    if (max(Y) < 50) {
        Y <- 2^Y
    }

    # quantile normalization of mixture file
    # library(preprocessCore)
    if (QN == TRUE) {
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }

    # store original mixtures
    Yorig <- Y
    Ymedian <- max(stats::median(Yorig), 1)

    # intersect genes
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    YintX <- Ygns %in% Xgns
    Y <- as.matrix(Y[YintX, ])
    colnames(Y) <- colname
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY, ]

    # standardize sig matrix
    X <- (X - mean(X)) / stats::sd(as.vector(X))

    # empirical null distribution of correlation coefficients
    if (P > 0) {
        nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)
    }

    header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
    if (absolute) {
        header <- c(
            header,
            paste("Absolute score (", abs_method, ")",
                sep = ""
            )
        )
    }

    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999

    if (mixtures == 1) {
        y <- Y[, itor]

        # standardize mixture
        y <- (y - mean(y)) / stats::sd(y)

        # run SVR core algorithm
        result <- CoreAlg(X, y, absolute, abs_method)

        # get results
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse

        if (absolute && abs_method == "sig.score") {
            w <- w * stats::median(Y[, itor]) / Ymedian
        }

        # calculate p-value
        if (P > 0) {
            pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
        }

        # print output
        out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
        if (absolute) out <- c(out, sum(w))
        output <- matrix(out, nrow = 1)
    } else {
        # iterate through mixtures
        while (itor <= mixtures) {
            y <- Y[, itor]

            # standardize mixture
            y <- (y - mean(y)) / stats::sd(y)

            # run SVR core algorithm
            result <- CoreAlg(X, y, absolute, abs_method)

            # get results
            w <- result$w
            mix_r <- result$mix_r
            mix_rmse <- result$mix_rmse

            if (absolute && abs_method == "sig.score") {
                w <- w * stats::median(Y[, itor]) / Ymedian
            }

            # calculate p-value
            if (P > 0) {
                pval <- 1 - (which.min(abs(nulldist - mix_r)) /
                    length(nulldist))
            }

            # print output
            out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
            if (absolute) out <- c(out, sum(w))
            if (itor == 1) {
                output <- out
            } else {
                output <- rbind(output, out)
            }

            itor <- itor + 1
        }
    }



    # save results
    # write.table(rbind(header,output), file="CIBERSORT-Results.txt",
    # sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    # return matrix object containing all results
    obj <- rbind(header, output)
    obj <- obj[, -1]
    obj <- obj[-1, ]
    if (mixtures == 1) {
        obj <- matrix(as.numeric(obj), nrow = 1)
    } else {
        obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
    }
    rownames(obj) <- colnames(Y)
    if (!absolute) {
        colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
    } else {
        colnames(obj) <- c(
            colnames(X), "P-value", "Correlation", "RMSE",
            paste("Absolute score (", abs_method, ")", sep = "")
        )
    }
    obj
}








# Below are from ajive (https://github.com/idc9/r_jive)
## The singular value threshold.
##
## Computes the singluar value theshold for the data matrix (half way
# between the rank and rank + 1 singluar value).

get_sv_threshold <- function(singular_values, rank) {
    .5 * (singular_values[rank] + singular_values[rank + 1])
}

## Computes the joint scores.
##
## Estimate the joint rank with the wedin bound, compute the signal scores SVD,
# double check each joint component.
##

get_joint_scores <- function(blocks, block_svd, initial_signal_ranks,
                             sv_thresholds,
                             n_wedin_samples = 1000, n_rand_dir_samples = 1000,
                             joint_rank = NA) {
    if (is.na(n_wedin_samples) & is.na(n_rand_dir_samples) &
        is.na(joint_rank)) {
        stop("at least one of n_wedin_samples, n_rand_dir_samples,
             or joint_rank must not be NA",
            call. = FALSE
        )
    }

    K <- length(blocks)
    n_obs <- dim(blocks[[1]])[1]

    # SVD of the signal scores matrix -----------------------------------------
    signal_scores <- list()
    for (k in seq_len(K)) {
        signal_scores[[k]] <-
            block_svd[[k]][["u"]][, seq_len(initial_signal_ranks[k])]
    }

    M <- do.call(cbind, signal_scores)
    M_svd <- get_svd(M, rank = min(initial_signal_ranks))


    # estimate joint rank with wedin bound and random direction bound ----------

    rank_sel_results <- list()
    rank_sel_results[["obs_svals"]] <- M_svd[["d"]]

    if (is.na(joint_rank)) {
        # maybe comptue wedin bound
        if (!is.na(n_wedin_samples)) {
            block_wedin_samples <- matrix(NA, K, n_wedin_samples)

            for (k in seq_len(K)) {
                block_wedin_samples[k, ] <- get_wedin_bound_samples(
                    X = blocks[[k]],
                    SVD = block_svd[[k]],
                    signal_rank = initial_signal_ranks[k],
                    num_samples = n_wedin_samples
                )
            }

            wedin_samples <- K - colSums(block_wedin_samples)
            wedin_svsq_threshold <- stats::quantile(wedin_samples, .05)

            rank_sel_results[["wedin"]] <- list(
                block_wedin_samples = block_wedin_samples,
                wedin_samples = wedin_samples,
                wedin_svsq_threshold = wedin_svsq_threshold
            )
        } else {
            wedin_svsq_threshold <- NA
        }

        # maybe compute random direction bound
        if (!is.na(n_rand_dir_samples)) {
            rand_dir_samples <-
                get_random_direction_bound(
                    n_obs = n_obs,
                    dims = initial_signal_ranks,
                    num_samples = n_rand_dir_samples
                )
            rand_dir_svsq_threshold <- stats::quantile(rand_dir_samples, .95)

            rank_sel_results[["rand_dir"]] <- list(
                rand_dir_samples = rand_dir_samples,
                rand_dir_svsq_threshold = rand_dir_svsq_threshold
            )
        } else {
            rand_dir_svsq_threshold <- NA
        }

        overall_sv_sq_threshold <- max(wedin_svsq_threshold,
            rand_dir_svsq_threshold,
            na.rm = TRUE
        )
        joint_rank_estimate <- sum(M_svd[["d"]]^2 > overall_sv_sq_threshold)

        rank_sel_results[["overall_sv_sq_threshold"]] <- overall_sv_sq_threshold
        rank_sel_results[["joint_rank_estimate"]] <- joint_rank_estimate
    } else { # user provided joint rank
        joint_rank_estimate <- joint_rank
        rank_sel_results[["joint_rank_estimate"]] <- joint_rank
    }


    # estimate joint score space ------------------------------------

    if (joint_rank_estimate >= 1) {
        joint_scores <- M_svd[["u"]][, seq_len(joint_rank_estimate),
            drop = FALSE
        ]

        # reconsider joint score space ------------------------------------
        # remove columns of joint_scores that have a
        # trivial projection from one of the data matrices

        to_remove <- c()
        for (k in seq_len(K)) {
            for (j in seq_len(joint_rank_estimate)) {
                score <- t(blocks[[k]]) %*% joint_scores[, j]
                sv <- norm(score)

                if (sv < sv_thresholds[[k]]) {
                    glue::glue(paste("removing column", j))
                    to_remove <- c(to_remove, j)
                    break
                }
            }
        }
        to_keep <- setdiff(seq_len(joint_rank_estimate), to_remove)
        joint_rank <- length(to_keep)
        joint_scores <- joint_scores[, to_keep, drop = FALSE]
    } else {
        joint_scores <- NA
    }


    list(joint_scores = joint_scores, rank_sel_results = rank_sel_results)
}

## Computes the final JIVE decomposition.
##
## Computes X = J + I + E for a single data block and the respective SVDs.
##

get_final_decomposition <- function(X, joint_scores,
                                    sv_threshold, full = TRUE) {
    jive_decomposition <- list()
    jive_decomposition[["individual"]] <-
        get_individual_decomposition(X, joint_scores, sv_threshold, full)
    jive_decomposition[["joint"]] <-
        get_joint_decomposition(X, joint_scores, full)


    if (full) {
        jive_decomposition[["noise"]] <-
            X - (jive_decomposition[["joint"]][["full"]] +

                jive_decomposition[["individual"]][["full"]])
    } else {
        jive_decomposition[["noise"]] <- NA
    }

    jive_decomposition
}

## Computes the individual matix for a data block.
##

get_individual_decomposition <- function(X, joint_scores, sv_threshold,
                                         full = TRUE) {
    if (any(is.na(joint_scores))) {
        indiv_decomposition <- get_svd(X)
    } else {
        X_orthog <- (diag(dim(X)[1]) - joint_scores %*% t(joint_scores)) %*% X
        indiv_decomposition <- get_svd(X_orthog)
    }


    indiv_rank <- sum(indiv_decomposition[["d"]] > sv_threshold)

    indiv_decomposition <- truncate_svd(
        decomposition = indiv_decomposition,
        rank = indiv_rank
    )

    if (full) {
        indiv_decomposition[["full"]] <- svd_reconstruction(indiv_decomposition)
    } else {
        indiv_decomposition[["full"]] <- NA
    }

    indiv_decomposition[["rank"]] <- indiv_rank
    indiv_decomposition
}

## Computes the joint matix for a data block.
##
##

get_joint_decomposition <- function(X, joint_scores, full = TRUE) {
    if (any(is.na(joint_scores))) {
        joint_decomposition <- list(full = NA, rank = 0, u = NA, d = NA, v = NA)
        return(joint_decomposition)
    }
    joint_rank <- dim(joint_scores)[2]
    J <- joint_scores %*% t(joint_scores) %*% X

    joint_decomposition <- get_svd(J, joint_rank)

    if (full) {
        joint_decomposition[["full"]] <- J
    } else {
        joint_decomposition[["full"]] <- NA
    }

    joint_decomposition[["rank"]] <- joint_rank
    joint_decomposition
}


## Singluar Value Decomposition.
##
## Returns a possibly truncated SVD of a data matrix.
##
## Wraps the svd function. Removes the extra singluar values if a truncated
# svd is computed.

get_svd <- function(X, rank = NULL) {
    # SVD <- get_svd(X, rank=2)

    if (is.null(rank)) {
        svd(X)
    } else if (rank == 0) {
        # TODO: what to do
        decomposition <- list()
        decomposition[["u"]] <- matrix(0, ncol = 1, nrow = dim(X)[1])
        decomposition[["d"]] <- 0
        decomposition[["v"]] <- matrix(0, ncol = 1, nrow = dim(X)[2])
        decomposition
    } else {
        decomposition <- svd(X, nu = rank, nv = rank)
        decomposition[["d"]] <- decomposition[["d"]][seq_len(rank)]
        decomposition
    }
}



## Truncates an SVD.
##
## Removes columns from the U, D, V matrix computed form an SVD.
##

truncate_svd <- function(decomposition, rank) {
    if (rank == 0) {
        n <- dim(decomposition[["u"]])[1]
        d <- dim(decomposition[["v"]])[1]
        decomposition[["u"]] <- matrix(0, ncol = 1, nrow = n)
        decomposition[["d"]] <- 0
        decomposition[["v"]] <- matrix(0, ncol = 1, nrow = d)
    } else {
        decomposition[["u"]] <-
            decomposition[["u"]][, seq_len(rank), drop = FALSE]
        decomposition[["d"]] <-
            decomposition[["d"]][seq_len(rank)]
        decomposition[["v"]] <-
            decomposition[["v"]][, seq_len(rank), drop = FALSE]
    }

    decomposition
}


## Reconstruces the original matrix from its SVD.
##
## Computes UDV^T to get the approximate (or full) X matrix.
##

svd_reconstruction <- function(decomposition) {
    # decomposition rank -- need to truncated singluar values
    r <- dim(decomposition[["u"]])[2]

    decomposition[["u"]] %*%
        diag(decomposition[["d"]][seq_len(r)], nrow = r, ncol = r) %*%
        t(decomposition[["v"]])
}



get_random_direction_bound <- function(n_obs, dims, num_samples = 1000) {
    n_blocks <- length(dims)
    rand_dir_samples <- rep(0, num_samples)
    for (s in seq_len(num_samples)) {
        rand_subspaces <- list()
        for (b in seq(n_blocks)) {
            X <- matrix(
                stats::rnorm(n_obs * dims[b], mean = 0, sd = 1),
                n_obs, dims[b]
            )
            U <- get_svd(X)[["u"]]

            rand_subspaces[[b]] <- U
        }
        M <- do.call(cbind, rand_subspaces)
        M_svd <- get_svd(M, rank = min(dims))

        rand_dir_samples[s] <- M_svd[["d"]][1]^2
    }

    rand_dir_samples
}



get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples = 1000) {
    # resample for U and V
    U_perp <- SVD[["u"]][, -(seq_len(signal_rank))]
    U_sampled_norms <- wedin_bound_resampling(
        X = X,
        perp_basis = U_perp,
        right_vectors = FALSE,
        num_samples = num_samples
    )

    V_perp <- SVD[["v"]][, -(seq_len(signal_rank))]
    V_sampled_norms <- wedin_bound_resampling(
        X = X,
        perp_basis = V_perp,
        right_vectors = TRUE,
        num_samples = num_samples
    )

    sigma_min <- SVD[["d"]][signal_rank]
    wedin_bound_samples <-
        mapply(
            function(u, v) min(max(u, v) / sigma_min, 1)^2,
            U_sampled_norms, V_sampled_norms
        )

    wedin_bound_samples
}


## Resampling procedure for the wedin bound
##

wedin_bound_resampling <- function(X, perp_basis, right_vectors,
                                   num_samples = 1000) {
    rank <- dim(perp_basis)[2]
    resampled_norms <- rep(0, num_samples)

    for (s in seq_len(num_samples)) {
        sampled_col_index <- sample.int(
            n = dim(perp_basis)[2],
            size = rank,
            replace = TRUE
        )


        perp_resampled <- perp_basis[, sampled_col_index]

        if (right_vectors) {
            resampled_projection <- X %*% perp_resampled
        } else {
            resampled_projection <- t(perp_resampled) %*% X
        }

        # operator L2 norm
        resampled_norms[s] <- norm(resampled_projection,
            type = "2"
        )
    }

    resampled_norms
}

## Returns the common normalized scores.
##
## Represents the joint signal among the data blocks.
##

get_common_normalized_scores <- function(ajive_output) {
    ajive_output[["joint_scores"]]
}

## Angle based Joint and Individual Variation Explained
##
## Computes the JIVE decomposition.
##
##
##
ajive <- function(blocks, initial_signal_ranks, full = TRUE,
                  n_wedin_samples = 1000, n_rand_dir_samples = 1000,
                  joint_rank = NA) {
    K <- length(blocks)

    if (K < 2) {
        stop("ajive expects at least two data matrices.")
    }

    if (sum(vapply(blocks, function(X) any(is.na(X)), logical(1))) > 0) {
        stop("Some of the blocks has missing data --
             ajive expects full data matrices.")
    }

    # TODO: should we give the option to center the data?
    # if(center){
    #     blocks <- lapply(blocks,function(X) scale(X, center=T, scale=FALSE))
    # }

    # step 1: initial signal space extraction --------------------------------
    # initial estimate of signal space with SVD

    block_svd <- list()
    sv_thresholds <- rep(0, K)
    for (k in seq_len(K)) {
        block_svd[[k]] <- get_svd(blocks[[k]])

        sv_thresholds[k] <- get_sv_threshold(
            singular_values = block_svd[[k]][["d"]],
            rank = initial_signal_ranks[k]
        )
    }


    # step 2: joint sapce estimation -------------------------------------------

    out <- get_joint_scores(blocks, block_svd, initial_signal_ranks,
        sv_thresholds,
        n_wedin_samples = n_wedin_samples,
        n_rand_dir_samples = n_rand_dir_samples,
        joint_rank = joint_rank
    )
    joint_rank_sel_results <- out$rank_sel_results
    joint_scores <- out$joint_scores

    joint_rank <- out[["rank_sel_results"]][["joint_rank_estimate"]]
    # dim(joint_scores)[2]


    # step 3: final decomposition ----------------------------------------------

    block_decomps <- list()
    for (k in seq_len(K)) {
        block_decomps[[k]] <- get_final_decomposition(
            X = blocks[[k]],
            joint_scores = joint_scores,
            sv_threshold = sv_thresholds[k]
        )
    }

    jive_decomposition <- list(block_decomps = block_decomps)
    jive_decomposition[["joint_scores"]] <- joint_scores
    jive_decomposition[["joint_rank"]] <- joint_rank

    jive_decomposition[["joint_rank_sel"]] <- joint_rank_sel_results
    jive_decomposition
}
