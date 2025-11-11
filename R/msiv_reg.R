#' This function performs Instrumental Variable Regression using multiple synthetic IVs
#'
#' The synthetic IVs are identified using the dual tendency DT condition for coplanar regression vectors.
#' The DT is formulated for both the homoscedastic and heteroscedastic cases.
#' @param data A data frame.
#' @param Y Name of the dependent variable.
#' @param X Character vector of endogenous variable names.
#' @param H Character vector of exogenous variable names.
#' @param reps Number of bootstrap loops.
#'
#' @return A list containing IV regression results and generated SIVs.
#' @export
#' @importFrom stats lm resid cov cor pf rnorm as.formula predict fitted complete.cases ecdf
#' @importFrom sandwich vcovHC
#' @importFrom AER ivreg
#' @examples
#' df <- wooldridge::mroz
#' data <- df[complete.cases(df), ]
#' result <- msiv_reg(data, "hours", c("lwage", "educ"),
#'                    c("age", "kidslt6", "kidsge6", "nwifeinc"), reps=5)
msiv_reg <- function(data, Y, X, H, reps = 5) {

  rad2deg <- function(rad) (rad * 180) / pi

  # Safe formula builder
  make_formula <- function(lhs, rhs_vars) {
    rhs_vars <- rhs_vars[rhs_vars != ""]
    rhs_str <- paste(rhs_vars, collapse = " + ")
    if (lhs != "") {
      as.formula(paste(lhs, "~", rhs_str))
    } else {
      as.formula(paste("~", rhs_str))
    }
  }

  # Two-sample Anderson-Darling statistic
  ad2_stat <- function(x, y) {
    if (!is.numeric(x) || !is.numeric(y)) stop("Inputs must be numeric vectors")
    if (length(x) < 2 || length(y) < 2) stop("Samples must have at least two observations")
    n <- length(x)
    m <- length(y)
    z <- sort(c(x, y))
    z <- z[-length(z)]
    ecdf_x <- ecdf(x)
    ecdf_y <- ecdf(y)
    W <- rank(z) / (n + m)
    eps <- .Machine$double.eps
    stat <- (n * m / (n + m)^2) * sum((ecdf_x(z) - ecdf_y(z))^2 / pmax((1 - W) * W, eps))
    return(stat)
  }

  check_sign_change <- function(x) {
    signs <- sign(x)
    sign_changes <- any(diff(signs) != 0, na.rm = TRUE)
    return(as.integer(sign_changes))
  }

  check_initial_abs_increase <- function(x) {
    abs_x <- abs(x)
    initial_increasing <- all(diff(abs_x[1:min(20, length(abs_x))]) >= 0)
    return(as.integer(initial_increasing))
  }

  find_first_sign_change <- function(x) {
    sign_changes <- which(diff(sign(x)) != 0)
    if (length(sign_changes) > 0) {
      return(sign_changes[1] + 1)
    } else {
      return(NA)
    }
  }

  # Ensure variables exist in data
  if (!all(c(Y, X, H) %in% names(data))) {
    stop("Some variables are missing in the dataset.")
  }

  if (reps > 10) reps <- 10
  if (reps < 2) reps <- 2

  N <- nrow(data)
  siv_list <- list()
  signk <- numeric(length(X))

  # Loop for each endogenous variable
  for (g in seq_along(X)) {
    xname <- X[g]
    rhs_y <- c(H, setdiff(X, xname))

    # Factoring out effects of other exogenous variables
    formula_y <- make_formula(Y, rhs_y)
    fity <- lm(formula_y, data = data)
    y <- resid(fity)

    formula_x <- make_formula(xname, rhs_y)
    fitx <- lm(formula_x, data = data)
    x <- resid(fitx)

    # Saving transformed x and y
    data$x <- (x - mean(x))
    data$y <- (y - mean(y))
    x0 <- x
    y0 <- y

    # Generating vector orthogonal to x
    fity <- lm(y0 ~ x0, data = data)
    V <- resid(fity)
    V <- (V - mean(V)) / sd(V) * sd(x0)
    data$R <- V

    # Determining the sign of cor(x,u)
    signc <- matrix(0, ncol = 5, nrow = 2)
    signc[, 1] <- c(1, -1)

    for (j in 1:2) {
      k <- signc[j, 1]

      dd <- 3
      d <- 0.01
      delt <- 0.01
      i <- 1
      dl <- round(dd / delt)
      m1 <- numeric(dl)
      m2 <- numeric(dl)

      while (d < dd) {
        data$siv <- (data$x - k * d * data$R)
        rls <- lm(data$x ~ data$siv, data = data)
        s.rls <- summary(rls, vcov. = function(x) vcovHC(x, type = "HC1"), diagnostics = TRUE)
        data$ev21 <- resid(s.rls)
        m1[i] <- cov(data$ev21^2, data$siv)
        m2[i] <- cor(data$ev21^2, data$siv)
        d <- d + delt
        i <- i + 1
      }

      signc[j, 3] <- check_initial_abs_increase(m1)
      signc[j, 5] <- check_initial_abs_increase(m2)

      if (signc[j, 3] != 1) {
        index <- find_first_sign_change(m1)
        if (!is.na(index) && index <= length(m1)) {
          m <- m1[index:length(m1)]
        } else {
          m <- m1
        }
      } else {
        m <- m1
      }

      signc[j, 2] <- check_sign_change(m1)
      signc[j, 4] <- check_sign_change(m2)
    }

    ch <- rowSums(signc[, -1])
    k <- signc[which.max(ch), 1]
    signk[g] <- k

    if (k != 0) {
      # Bootstrap sampling loop
      d0i <- numeric(reps)
      d0ri <- numeric(reps)
      d0rni <- numeric(reps)
      S <- round(N * 0.999)
      l <- 1

      while (l < reps) {
        set.seed(3 * l)
        data_sample <- data[sample(1:N, S, replace = TRUE), ]

        dd <- 3
        d <- 0.01
        delt <- 0.01
        i <- 1
        dl <- round(dd / delt)
        m1 <- numeric(dl)
        dv2 <- numeric(dl)
        x4 <- numeric(dl)
        vvar <- numeric(dl)
        st <- numeric(dl)

        while (d < dd) {
          data_sample$siv <- (data_sample$x - k * d * data_sample$R)
          rls <- lm(x ~ siv, data = data_sample)
          s.rls <- summary(rls, vcov. = function(x) vcovHC(x, type = "HC1"), diagnostics = TRUE)
          data_sample$ev21 <- resid(s.rls)

          # FGLS
          ehatsq <- resid(rls)^2
          sighatsq.rls <- lm(log(ehatsq) ~ siv, data = data_sample)
          data_sample$vari <- sqrt(exp(fitted(sighatsq.rls)))
          vvar[i] <- var(data_sample$vari)
          fgls <- lm(x ~ siv, weights = 1 / data_sample$vari, data = data_sample)
          data_sample$ev22 <- resid(fgls)
          m1[i] <- cor(data_sample$ev21^2, data_sample$siv)

          n <- length(data_sample$ev21)
          ssr1 <- sum((predict(lm((ev21^2) ~ siv, data = data_sample)) - mean(data_sample$ev21^2))^2)
          sse1 <- sum(data_sample$ev21^2)
          x1 <- (ssr1 / 2) / (sse1 / n^2)^2

          ssr2 <- sum((predict(lm((ev22^2) ~ siv, data = data_sample)) - mean(data_sample$ev22^2))^2)
          sse2 <- sum(data_sample$ev22^2)
          x2 <- (ssr2 / 2) / (sse2 / n^2)^2

          x3 <- x1 / x2
          dv2[i] <- pf(x3, df1 = 1, df2 = 1, lower.tail = TRUE)
          st[i] <- d

          # Non-parametric
          samp1 <- (predict(lm((data_sample$ev21^2) ~ data_sample$siv, data = data_sample)))^2
          samp2 <- (predict(lm((data_sample$ev22^2) ~ data_sample$siv, data = data_sample)))^2
          ad0 <- ad2_stat(x = samp1, y = samp2)
          x4[i] <- 1 - ad0
          d <- d + delt
          i <- i + 1
        }

        # Simple delta estimate
        #d0i[l] <- which.min(abs(m1)) * delt
        rvec <- data_sample$R
        xvec <- data_sample$x
        # Run the enhanced delta sweep
        results <- gmm_test_delta_sweep(xvec, rvec, k=k, delta_max = dd, n_deltas = 200)
        # Find interesting points
        min_j_idx <- which.min(results$J_stat)
        d0i[l] <- results$delta[min_j_idx]
        # Parametric heteroscedastic
        d0ri[l] <- which.min(dv2) * delt

        # Non-parametric heteroscedastic
        d0rni[l] <- which.min(x4) * delt

        l <- l + 1
      }

      # Random shocks for SIVs
      v1 <- rnorm(N, 0, sd(x))
      v2 <- rnorm(N, 0, sd(x))
      v3 <- rnorm(N, 0, sd(x))

      # Simple homogeneous case
      d0i <- d0i[complete.cases(d0i)]
      d0m <- mean(d0i)
      siv_list[[paste0("siv1", g)]] <- (data$x - k * d0m * data$R)
      siv_list[[paste0("siv1b", g)]] <- (data$x - k * d0m * data$R + v1)

      # Parametric heterogeneous case
      d0ri <- d0ri[complete.cases(d0ri)]
      d0rm <- mean(d0ri)
      siv_list[[paste0("siv2", g)]] <- (data$x - k * d0rm * data$R)
      siv_list[[paste0("siv2b", g)]] <- (data$x - k * d0rm * data$R + v2)

      # Non-parametric heterogeneous case
      d0rni <- d0rni[complete.cases(d0rni)]
      d0rnm <- mean(d0rni)
      siv_list[[paste0("siv3", g)]] <- (data$x - k * d0rnm * data$R)
      siv_list[[paste0("siv3b", g)]] <- (data$x - k * d0rnm * data$R + v3)

    } else {
      print("NO endogeneity problem. All SIV estimates are the same as OLS")
      d0m <- 0.001
      d0rm <- 0.001
      d0rnm <- 0.001

      siv_list[[paste0("siv1", g)]] <- (data$x - k * d0m * data$R)
      siv_list[[paste0("siv2", g)]] <- (data$x - k * d0rm * data$R)
      siv_list[[paste0("siv3", g)]] <- (data$x - k * d0rnm * data$R)
    }
  }

  # Add SIVs to data
  for (name in names(siv_list)) {
    data[[name]] <- siv_list[[name]]
  }

  # Extract SIV variable names
  vars1 <- names(siv_list)[grepl("^siv1[0-9]*$", names(siv_list))]
  vars2 <- names(siv_list)[grepl("^siv2[0-9]*$", names(siv_list))]
  vars3 <- names(siv_list)[grepl("^siv3[0-9]*$", names(siv_list))]

  # Build formulas
  formula <- make_formula(Y, c(X, H))
  instruments1 <- make_formula("", c(vars1, H))
  instruments2 <- make_formula("", c(vars2, H))
  instruments3 <- make_formula("", c(vars3, H))

  # Run IV regressions
  iv1 <- AER::ivreg(formula, instruments1, data = data)
  iv2 <- AER::ivreg(formula, instruments2, data = data)
  iv3 <- AER::ivreg(formula, instruments3, data = data)

  return(list(
    IV1 = iv1,
    IV2 = iv2,
    IV3 = iv3,
    siv_list = siv_list,
    signk = signk
  ))
}
