#' Robust version of GMM test for the delta sweep
#'
#' @param E A numeric vector.
#' @param S A numeric vector.
#' @param regularization a tolerance
#' @return The index of the first sign change, or NA if no sign change occurs.
#' @export
gmm_test_homoskedasticity_robust <- function(E, S, regularization = 1e-8) {
  n <- nrow(E)
  p <- ncol(E)
  q <- ncol(S)

  # Estimate common variance under null
  sigma2_hat <- mean(rowSums(E^2)) / p

  # Compute moment conditions: vec(e_i e_i' - sigma^2 I) ⊗ s_i
  moments <- matrix(0, n, p^2 * q)
  I_p <- diag(p)

  for (i in 1:n) {
    ee_matrix <- outer(E[i, ], E[i, ]) - sigma2_hat * I_p
    ee_vec <- as.vector(ee_matrix)
    moments[i, ] <- kronecker(ee_vec, S[i, ])
  }

  # Sample moments
  m_bar <- colMeans(moments)
  S_hat <- var(moments)

  # Regularization for numerical stability
  eigenvals <- eigen(S_hat, only.values = TRUE)$values
  min_eigenval <- min(eigenvals)
  max_eigenval <- max(eigenvals)

  if (min_eigenval < regularization) {
    S_hat <- S_hat + diag(regularization, nrow(S_hat))
    min_eigenval <- min_eigenval + regularization
  }

  # Compute J-statistic
  J_stat <- n * t(m_bar) %*% solve(S_hat, m_bar)
  df <- p^2 * q
  p_value <- 1 - pchisq(J_stat, df)

  return(list(
    J = as.numeric(J_stat),
    df = df,
    p_value = p_value,
    sigma2_hat = sigma2_hat,
    condition_number = max_eigenval / min_eigenval
  ))
}
