#' Compute gamma
#'
#' @param X A numeric vector.
#' @param R A numeric vector.
#' @param k A numeric parameter.
#' @param delta_max Maximum for the range.
#' @param n_deltas A number of steps.
#' @param compute_gamma  A logical value
#' @return The index of the first sign change, or NA if no sign change occurs.
#' @export
gmm_test_delta_sweep <- function(X, R, delta_max, k, n_deltas = 100,
                                 compute_gamma = TRUE) {
  # X, R: n x 1 orthogonal vectors
  # delta_max: maximum delta value to test
  # n_deltas: number of delta values to test
  # compute_gamma: if TRUE, estimate gamma via OLS for each delta

  n <- length(X)
#
#   # Check orthogonality
#   if (abs(sum(X * R)) > 1e-10) {
#     warning("X and R are not perfectly orthogonal. Inner product = ", sum(X * R))
#   }

  # Create delta grid
  deltas <- seq(0.01, delta_max, length.out = n_deltas)  # Start from 0.01 to avoid S=X

  # Storage for results
  results <- data.frame(
    delta = deltas,
    gamma_hat = NA,
    J_stat = NA,
    df = NA,
    p_value = NA,
    sigma2_hat = NA,
    condition_number = NA,
    S_variance = NA
  )
  k=-1
  for (i in 1:length(deltas)) {
    delta <- deltas[i]

    # Construct instrument: S = X - delta * R
    S <- X -k*delta * R

    # Estimate gamma via OLS: E = X - gamma * S + error
    if (compute_gamma) {
      # OLS: gamma_hat = (S'S)^(-1) S'X
      gamma_hat <- sum(S * X) / sum(S^2)
    } else {
      gamma_hat <- 1  # Fixed gamma
    }

    # Compute OLS residual: E = X - gamma_hat * S
    E <- X - gamma_hat * S

    # For univariate case, convert to matrix form for GMM test
    E_matrix <- matrix(E, ncol = 1)
    S_matrix <- matrix(S, ncol = 1)

    # Run GMM test
    tryCatch({
      gmm_result <- gmm_test_homoskedasticity_robust(E_matrix, S_matrix)

      results$gamma_hat[i] <- gamma_hat
      results$J_stat[i] <- gmm_result$J
      results$df[i] <- gmm_result$df
      results$p_value[i] <- gmm_result$p_value
      results$sigma2_hat[i] <- gmm_result$sigma2_hat
      results$condition_number[i] <- gmm_result$condition_number
      results$S_variance[i] <- var(S)

    }, error = function(e) {
      # Handle numerical issues
      results$gamma_hat[i] <<- gamma_hat
      results$J_stat[i] <<- NA
      results$p_value[i] <<- NA
      results$S_variance[i] <<- var(S)
    })
  }

  return(results)
}
