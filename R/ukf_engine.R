# =============================================================================
# ukf_engine.R
# Core Unscented Kalman Filter (UKF) implementation.
#
# Functions:
#   propagate_model()  -- RK4 propagation of sigma points through the ODE
#   UKF_dT()           -- Streamlined UKF prediction + update for one time step
#   UKF_blend()        -- Full pass over a time series, calling UKF_dT at each step
#
# Bugs fixed vs. original:
#   - pracma dependency made explicit (loaded in main script)
#   - Noise injection matrix dimensions generalized from hardcoded N_y=2
#   - Chi-square computed against clean smoothed data, not against noisy x
#   - Cholesky jitter loop extracted to helper for DRY code
#   - forcePositive clamp applied to ALL parameter rows uniformly
#   - Vectorized sigma-point covariance accumulation (no inner for-loops)
# =============================================================================

source("R/constants.R")


# -----------------------------------------------------------------------------
# Internal helper: make a square matrix positive-definite via adaptive jitter,
# falling back to Matrix::nearPD if jitter alone is insufficient.
# Returns the stabilized matrix.
.stabilize_pd <- function(M, label = "M") {
  n       <- nrow(M)
  jitter  <- 0
  success <- FALSE
  first   <- TRUE

  while (!success && jitter <= UKF_CONSTANTS$JITTER_MAX) {
    Mj      <- M + diag(jitter, n)
    eigvals <- eigen(Mj, symmetric = TRUE, only.values = TRUE)$values
    cond    <- max(eigvals) / (min(eigvals) + .Machine$double.eps)

    if (all(eigvals > UKF_CONSTANTS$EIGVAL_MIN) &&
        cond < UKF_CONSTANTS$COND_NUM_MAX) {
      M       <- Mj
      success <- TRUE
    } else {
      jitter <- if (first) { first <- FALSE; UKF_CONSTANTS$JITTER_INIT } else jitter * 10
    }
  }

  if (!success) {
    if (!requireNamespace("Matrix", quietly = TRUE))
      stop(label, ": stabilization failed and package 'Matrix' is unavailable.")
    M <- as.matrix(Matrix::nearPD(M)$mat)
  }
  M
}


# -----------------------------------------------------------------------------
#' propagate_model
#'
#' Propagates the augmented state matrix through the ODE for one dT step using
#' a fixed-step RK4 integrator with nn = round(dT/dt) sub-steps.
#'
#' @param t_dummy  Scalar dummy time (ODE models do not use explicit time).
#' @param ode_model  Function with signature ode_model(t, x, p).
#' @param dt   Sub-step size (should be 0.1 * dT or smaller).
#' @param dT   Data time step (interval between fMRI volumes).
#' @param N_p  Number of model parameters (top rows of augmented state).
#' @param x    Augmented state matrix (N_p + N_y) x N_sigma.
#' @return     Propagated augmented state matrix, same dimensions as x.
#' @export
propagate_model <- function(t_dummy, ode_model, dt, dT, N_p, x) {
  stopifnot(is.numeric(t_dummy), length(t_dummy) == 1,
            is.matrix(x), N_p >= 1, dt > 0, dT > 0, dt <= dT)

  nn <- max(1L, round(dT / dt))

  p <- x[seq_len(N_p), , drop = FALSE]
  y <- x[(N_p + 1):nrow(x), , drop = FALSE]

  # Warn if parameters are in a potentially stiff regime
  if (any(abs(p) > UKF_CONSTANTS$STIFFNESS_WARN)) {
    warning("propagate_model: |parameter| > ", UKF_CONSTANTS$STIFFNESS_WARN,
            ". Consider reducing dt for numerical stability.")
  }

  for (n in seq_len(nn)) {
    k1 <- dt * ode_model(t_dummy, y,           p)
    k2 <- dt * ode_model(t_dummy, y + k1 / 2,  p)
    k3 <- dt * ode_model(t_dummy, y + k2 / 2,  p)
    k4 <- dt * ode_model(t_dummy, y + k3,       p)
    y  <- y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }

  rbind(p, y)
}


# -----------------------------------------------------------------------------
#' UKF_dT
#'
#' Streamlined Unscented Kalman Filter for a single time step dT.
#' Generates 2*N_x sigma points (no central point), propagates them through the
#' ODE, then performs the Kalman measurement update.
#'
#' @param t_dummy   Scalar dummy time variable.
#' @param ode_model ODE model function(t, x, p).
#' @param xhat      Augmented state vector (N_x x 1) from previous step.
#' @param Pxx       Augmented covariance matrix (N_x x N_x).
#' @param y_obs     Observed state vector (N_y x 1) at current time step.
#' @param N_p       Number of model parameters.
#' @param N_y       Number of observed state variables.
#' @param R         Observation noise covariance (N_y x N_y).
#' @param dt        RK4 sub-step size.
#' @param dT        Data time step.
#' @param R_scale   Scalar SD of observation noise.
#' @param Q_scale   Scalar SD of process noise on parameters.
#' @param forcePositive  Logical; clamp parameters to > PARAM_MIN if TRUE.
#' @return List: xhat (N_x x 1), Pxx (N_x x N_x), K (N_x x N_y).
#' @export
UKF_dT <- function(t_dummy, ode_model, xhat, Pxx, y_obs,
                   N_p, N_y, R, dt, dT, R_scale, Q_scale,
                   forcePositive = FALSE) {

  N_x     <- N_p + N_y
  N_sigma <- 2L * N_x

  # ---- Sigma points --------------------------------------------------------
  Pxx   <- .stabilize_pd(Pxx, "Pxx")
  S     <- tryCatch(t(chol(N_x * Pxx, pivot = TRUE)),
                    error = function(e) {
                      Pxx_pd <- as.matrix(Matrix::nearPD(N_x * Pxx)$mat)
                      t(chol(Pxx_pd, pivot = TRUE))
                    })

  Xa <- matrix(xhat, nrow = N_x, ncol = N_sigma) + cbind(S, -S)

  # ---- Propagate -----------------------------------------------------------
  X      <- propagate_model(t_dummy, ode_model, dt, dT, N_p, Xa)
  xtilde <- rowMeans(X)                          # (N_x x 1)

  # ---- Predicted state covariance (vectorized) -----------------------------
  dX  <- X - xtilde                              # (N_x x N_sigma) deviations
  Pxx <- (dX %*% t(dX)) / N_sigma

  # Process noise on parameter block only
  Q               <- matrix(0, N_x, N_x)
  Q[1:N_p, 1:N_p] <- Q_scale^2 * diag(N_p) #old version: Q[1:N_p, 1:N_p] <- Q_scale * diag(N_p)
  Pxx             <- Pxx + Q

  # ---- Observation prediction (vectorized) ---------------------------------
  Y      <- X[(N_p + 1):(N_p + N_y), , drop = FALSE]
  ytilde <- rowMeans(Y)                          # (N_y x 1)

  dY  <- Y - ytilde                              # (N_y x N_sigma) deviations
  Pyy <- R + (dY %*% t(dY)) / N_sigma
  Pyy <- .stabilize_pd(Pyy, "Pyy")

  # Cross-covariance (vectorized)
  Pxy <- (dX %*% t(dY)) / N_sigma               # (N_x x N_y)

  # ---- Kalman gain ---------------------------------------------------------
  K <- tryCatch(Pxy %*% solve(Pyy),
                error = function(e) Pxy %*% MASS::ginv(Pyy))

  # ---- State update --------------------------------------------------------
  xhat <- xtilde + K %*% (y_obs - ytilde)

  if (forcePositive)
    xhat[1:N_p] <- pmax(UKF_CONSTANTS$PARAM_MIN, xhat[1:N_p])

  # ---- Covariance update ---------------------------------------------------
  Pxx <- Pxx - K %*% t(Pxy)
  Pxx <- .stabilize_pd(Pxx, "Pxx post-update")

  list(xhat = xhat, Pxx = Pxx, K = K)
}


# -----------------------------------------------------------------------------
#' UKF_blend
#'
#' Applies UKF_dT to every time point in ts_data (one forward pass).
#' The chi-square goodness-of-fit is computed against the clean smoothed data,
#' not against the noise-injected copy used internally.
#'
#' @param t_dummy     Scalar dummy time (NOT a time vector).
#' @param ts_data     Matrix (T x (1 + N_y)): col 1 = time, cols 2:(N_y+1) = signals.
#' @param ode_model   ODE model function(t, x, p).
#' @param N_p         Number of model parameters.
#' @param N_y         Number of observed state variables.
#' @param param_guess Initial parameter vector of length N_p.
#' @param dt          RK4 sub-step size.
#' @param dT          Data time step.
#' @param R_scale     Observation noise scale.
#' @param Q_scale     Process noise scale.
#' @param forcePositive  Logical; clamp parameters positive if TRUE.
#' @param seeded      Logical; set seed for reproducible noise injection.
#' @return List: param_est, xhat, error, chisq.
#' @export
UKF_blend <- function(t_dummy, ts_data, ode_model, N_p, N_y,
                      param_guess, dt, dT,
                      R_scale = 0.3, Q_scale = 0.015,
                      forcePositive = FALSE, seeded = FALSE) {

  # --- Input validation -----------------------------------------------------
  stopifnot(
    "t_dummy must be a scalar" = is.numeric(t_dummy) && length(t_dummy) == 1,
    "ts_data must be a matrix" = is.matrix(ts_data),
    "ts_data must have N_y + 1 columns" = ncol(ts_data) == N_y + 1,
    "param_guess length must equal N_p" = length(param_guess) == N_p,
    "dt must be less than dT" = dt < dT,
    "R_scale must be positive" = R_scale > 0,
    "Q_scale must be positive" = Q_scale > 0
  )

  num_time  <- nrow(ts_data)
  N_x       <- N_p + N_y
  clean_obs <- t(ts_data[, -1, drop = FALSE])    # (N_y x T) — used for chi-sq

  # --- Initialize storage ---------------------------------------------------
  xhat   <- matrix(0, nrow = N_x, ncol = num_time)
  Pxx    <- vector("list", num_time)
  Ks     <- vector("list", num_time)
  errors <- matrix(0, nrow = N_x, ncol = num_time)

  for (i in seq_len(num_time))
    Pxx[[i]] <- matrix(0, N_x, N_x)
  for (i in seq_len(num_time))
    Ks[[i]]  <- matrix(0, N_x, N_y)

  # --- Augmented initial state: [params (repeated); observed data] ----------
  z  <- matrix(param_guess, nrow = N_p, ncol = num_time)
  y0 <- clean_obs                                # (N_y x T)
  x  <- rbind(z, y0)

  xhat[, 1] <- x[, 1]

  R        <- (R_scale)^2 * cov(t(x[(N_p + 1):(N_p + N_y), ]))
  # old version: Pxx[[1]] <- pracma::blkdiag(Q_scale * diag(N_p), R)
  Pxx[[1]] <- pracma::blkdiag(Q_scale^2 * diag(N_p), R) 

  # --- Noise-injected observations for UKF measurement updates --------------
  # FIX: dimensions now driven by N_y, not hardcoded to 2
  if (seeded) set.seed(1L)
  noise <- pracma::sqrtm(R)$B %*% matrix(rnorm(N_y * num_time),
                                          nrow = N_y, ncol = num_time)
  y_noisy <- clean_obs + noise

  # --- Forward pass ---------------------------------------------------------
  UKF_kstep <- list(xhat = xhat[, 1], Pxx = Pxx[[1]], K = Ks[[1]])

  for (k in 2:num_time) {
    UKF_kstep <- tryCatch(
      UKF_dT(t_dummy, ode_model,
             xhat[, k - 1], Pxx[[k - 1]], y_noisy[, k],
             N_p, N_y, R, dt, dT, R_scale, Q_scale, forcePositive),
      error = function(e) {
        message(sprintf("UKF_blend: step %d skipped — %s", k, conditionMessage(e)))
        UKF_kstep   # carry forward previous estimate
      }
    )
    xhat[, k]   <- UKF_kstep$xhat
    Pxx[[k]]    <- UKF_kstep$Pxx
    Ks[[k]]     <- UKF_kstep$K
    errors[, k] <- sqrt(pmax(0, diag(Pxx[[k]])))
  }

  param_estimated <- xhat[1:N_p, num_time, drop = FALSE]

  # FIX: chi-square against CLEAN data, not noise-injected copy
  resid <- clean_obs - xhat[(N_p + 1):(N_p + N_y), , drop = FALSE]
  chisq <- mean(resid^2)

  list(
    param_est = param_estimated,
    xhat      = xhat,
    error     = t(errors[1:N_p, num_time, drop = FALSE]),
    chisq     = chisq
  )
}
