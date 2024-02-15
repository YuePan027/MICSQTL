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
        X = assay(se)[mrk_prot, , drop = FALSE],
        W = ini_prop,
        refit_W = TRUE,
        refit_W.sparsity = length(mrk_prot)
    )
    cross_prop <- tca_res$W
    metadata(se)$cross_prop <- cross_prop
    return(se)
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

get_block_loadings <- function(ajive_output, k, type){
    if(! type  %in% c('joint', 'individual')){
        stop('type must be: joint or individual')
    }
    
    ajive_output$block_decomps[[k]][[type]][['v']]
}

grad_p <- function(X1, Y1, X2, Y2, p, s1, s2){
    grad <- rowSums(t(as.vector(X1 %*% (p*s1) - Y1) * X1) * s1) +
        rowSums(t(as.vector(X2 %*% (p*s2) - Y2) * X2) * s2)
    return(grad)
}

grad_si <- function(X,Y,p,s){
    grad <- rowSums(t(as.vector(X %*% (p*s) - Y) * X) * p)
    return(grad)
}

grad_Xi <- function(X,Y,p,s){
    grad <- (X %*% (p*s) - Y) %*% matrix(p*s, nrow = 1)
    return(grad)
}

ini_prep <- function(cell_counts){
    est <- dirmult(cell_counts)
    cell.prop.alpha <- est$gamma
    return(cell.prop.alpha)
}


MICSQTL_optim <- function(Y1, Y2,
                           ini_p, ini_s,
                           X1, 
                           X2,
                           step_p,
                           step_s,
                           eps,
                           iter){
    p <- ini_p
    s1 <- ini_s[[1]]
    s2 <- ini_s[[2]]
    eps_t <- 0.2
    iter_t <- 1
    X1_init <- X1
    X2_init <- X2
    res <- list()
    while(eps_t > eps & iter_t < iter){
        p_u <- pmax(p - step_p * grad_p(X1, Y1, X2, Y2, p, s1, s2), 0)
        s1_u <- pmax(s1 - step_s * grad_si(X1, Y1, p, s1), 0)
        s2_u <- pmax(s2 - step_s * grad_si(X2, Y2, p, s2), 0)  
        X1_u <- X1 - 0.1 * grad_Xi(X1, Y1, p, s1)
        X2_u <- X2 - 0.1 * grad_Xi(X2, Y2, p, s2)
        
        if(min(Y1) >= 0){
            X1_u[X1_u < 0] <- 0
        } 
        if(min(Y2) >=0){
            X2_u[X2_u < 0] <- 0
        }
    
        eps_t <- max(abs(p_u - p),
                     abs(p*s1 - p_u*s1_u),
                     abs(p*s2 - p_u*s2_u),
                     max(abs(X1_u - X1)),
                     max(abs(X2_u - X2)))
        p <- p_u
        s1 <- s1_u
        s2 <- s2_u
        max_diff_X1 <- max(abs(X1_u - X1))
        max_diff_X2 <- max(abs(X2_u - X2))
        X1 <- X1_u
        X2 <- X2_u
        iter_t <- iter_t + 1
    }
    if(iter_t >= iter & eps_t > eps){
        warning("Max iteration reached without convergence")
    }
    cat(paste0("iter = ", iter_t, " eps = ", eps_t, "\n"))
    res <- list(iter = iter_t,
                L = sum(as.vector(X1 %*% (p*s1) - Y1)^2) +
                    sum(as.vector(X2 %*% (p*s2) - Y2)^2),
                p = p,
                s1 = s1,
                s2 = s2,
                eps = eps_t,
                X1 = X1,
                X2 = X2,
                max_diff_X1 = max_diff_X1,
                max_diff_X2 = max_diff_X2,
                prop1 = (p*s1) / sum(p*s1),
                prop2 = (p*s2) / sum(p*s2),
                ini_p = ini_p,
                X1_init = X1_init,
                X2_init = X2_init
    )
    return(res)
}