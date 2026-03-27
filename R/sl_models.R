# =============================================================================
# sl_models.R
# Stuart-Landau oscillator models for UKF parameter estimation on BOLD signals.
#
# The Stuart-Landau (SL) model is the normal form of a supercritical Hopf
# bifurcation, making it the minimal model for capturing the transition between
# damped (subcritical, a < 0) and sustained (supercritical, a > 0) oscillation.
#
# THESIS HYPOTHESIS (MDD criticality):
#   Healthy resting-state BOLD operates near criticality (a ≈ 0), while MDD
#   drives the system subcritical (a < 0).  rtfMRI-nf should shift a toward
#   zero (toward criticality) in the active group vs. sham.
#
# IDENTIFIABILITY (why SL where coupled pendulum failed):
#   - Coupled pendulum: coupling k modulates frequency splitting → undetectable
#     in 1/f BOLD spectrum (proven flat chi² surface, 2026-03 experiments).
#   - Stuart-Landau: coupling K modulates AMPLITUDE variance → detectable from
#     amplitude envelope, not spectral peak → identifiable from resting BOLD.
#
# MODEL (complex form, z = x + iy):
#   dz/dt = (a + iω)z  −  |z|²·z  +  K(z_j − z_i)
#
# MODEL (Cartesian real form for UKF):
#   dx/dt =  a·x − ω·y − (x²+y²)·x  [+ K(x_j − x_i) for coupled]
#   dy/dt =  ω·x + a·y − (x²+y²)·y  [+ K(y_j − y_i) for coupled]
#
# ANALYTIC SIGNAL INPUT:
#   Both x and y are observable via the Hilbert transform:
#   x(t) = BOLD signal (real part)
#   y(t) = Hilbert{BOLD}(t) (imaginary part, 90° phase-shifted)
#   r(t) = sqrt(x²+y²) = instantaneous amplitude envelope
#
# PARAMETERS AND UNITS (all in data-normalised units, TR = 2 s):
#   a  : bifurcation parameter (dimensionless)
#          a < 0 → subcritical damped oscillation (MDD hypothesis)
#          a = 0 → Hopf bifurcation (criticality)
#          a > 0 → supercritical limit cycle, amplitude r* = sqrt(a)
#   ω  : natural angular frequency (rad / TR)
#          BOLD band 0.01–0.10 Hz → ω ∈ [0.126, 1.257] rad/TR
#          ω_init from spectral peak of BOLD signal
#   K  : diffusive coupling strength (dimensionless, ≥ 0)
#          K = 0 → uncoupled oscillators
#          K > 0 → amplitude synchronisation (detectable from BOLD)
#
# STAGE STRUCTURE (mirrors two-stage LHO approach):
#   Stage 1  make_sl_single()           N_p=2 [a, ω],  N_y=2 [x, y]
#            → estimates (a_i, ω_i) per ROI independently
#   Stage 2  make_sl_paired_fixed_aw()  N_p=1 [K],     N_y=4 [x1,y1,x2,y2]
#            → estimates K_{ij} for ROI pairs, fixing a/ω from Stage 1
#
# BOUNDS SUMMARY:
#   a  : [-A_MAX,  A_MAX]   A_MAX = 2.0 (no physiological prior; wide)
#   ω  : [ω_MIN,   ω_MAX]   from BOLD band limits (see SL_BOUNDS below)
#   K  : [0,       K_MAX]   K_MAX = 5.0 (positive coupling only)
#
# =============================================================================


# --- Physiological bounds (TR = 2 s, BOLD band 0.01–0.10 Hz) ----------------
SL_BOUNDS <- list(
  A_MIN   = -2.0,                              # bifurcation lower bound
  A_MAX   =  2.0,                              # bifurcation upper bound
  OM_MIN  =  2 * pi * 0.01 * 2.0,             # ≈ 0.1257 rad/TR
  OM_MAX  =  2 * pi * 0.10 * 2.0,             # ≈ 1.2566 rad/TR
  K_MIN   =  0.0,
  K_MAX   =  5.0
)

# Default initial guesses ---------------------------------------------------
SL_INIT <- list(
  A_INIT  = 0.0,        # start at criticality — let UKF determine direction
  OM_INIT = 2 * pi * 0.032 * 2.0,   # ≈ dominant BOLD freq ~0.032 Hz
  K_INIT  = 0.1
)


# =============================================================================
# STAGE 1 — make_sl_single
# =============================================================================

# -----------------------------------------------------------------------------
#' make_sl_single
#'
#' Single-region Stuart-Landau oscillator for UKF Stage 1.
#' Estimates bifurcation parameter a and natural frequency ω from one ROI's
#' BOLD analytic signal.
#'
#' EQUATIONS:
#'   dx/dt = a·x − ω·y − (x²+y²)·x
#'   dy/dt = ω·x + a·y − (x²+y²)·y
#'
#' PARAMETERS: p[1,] = a, p[2,] = ω   (N_p = 2)
#' STATE:       x[1,] = x (BOLD),  x[2,] = y (Hilbert)   (N_y = 2)
#' OBSERVED:    both rows (ts_data cols: [time, x_BOLD, y_Hilbert])
#'
#' IDENTIFIABILITY:
#'   a → controls amplitude envelope decay/growth (τ = 1/|a| TRs)
#'   ω → controls oscillation phase rate (detectable from y(t) lead/lag)
#'   For a < 0 (MDD): signal decays, amplitude < 1/e at t ≈ 1/|a| TRs
#'   For a > 0:       limit cycle, amplitude → sqrt(a) asymptotically
#'
#' @param  none (pure factory — all parameters estimated by UKF)
#' @return ODE function with signature (t, x, p)
#' @export
make_sl_single <- function() {
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    # Clamp inside ODE to catch sigma-point excursions beyond iter-level bounds
    a   <- pmax(SL_BOUNDS$A_MIN,  pmin(SL_BOUNDS$A_MAX,  p[1, ]))
    om  <- pmax(SL_BOUNDS$OM_MIN, pmin(SL_BOUNDS$OM_MAX, p[2, ]))

    xr  <- x[1, ]
    xi  <- x[2, ]
    r2  <- pmin(xr^2 + xi^2, 25.0)   # clamp |z|² — prevents cubic blow-up

    dxr <- a * xr - om * xi - r2 * xr
    dxi <- om * xr + a * xi - r2 * xi

    rbind(dxr, dxi)
  }
}

sl_single_model <- make_sl_single()


# =============================================================================
# STAGE 2 — make_sl_paired_fixed_aw
# =============================================================================

# -----------------------------------------------------------------------------
#' make_sl_paired_fixed_aw
#'
#' Paired Stuart-Landau oscillator for UKF Stage 2.
#' Both a and ω are fixed per region from Stage 1; only coupling K is free.
#'
#' EQUATIONS (diffusive coupling on both real and imaginary parts):
#'   dx₁/dt = a₁·x₁ − ω₁·y₁ − (x₁²+y₁²)·x₁  +  K·(x₂−x₁)
#'   dy₁/dt = ω₁·x₁ + a₁·y₁ − (x₁²+y₁²)·y₁   +  K·(y₂−y₁)
#'   dx₂/dt = a₂·x₂ − ω₂·y₂ − (x₂²+y₂²)·x₂   +  K·(x₁−x₂)
#'   dy₂/dt = ω₂·x₂ + a₂·y₂ − (x₂²+y₂²)·y₂    +  K·(y₁−y₂)
#'
#' PARAMETERS: p[1,] = K   (N_p = 1)
#' STATE:       [x₁, y₁, x₂, y₂]  (N_y = 4)
#' OBSERVED:    all four rows (ts_data cols: [time, x1, y1, x2, y2])
#' FIXED:       a₁, ω₁, a₂, ω₂ closed over from Stage 1 estimates
#'
#' IDENTIFIABILITY OF K:
#'   K controls amplitude covariance: Cov(r₁, r₂) increases with K.
#'   Unlike frequency-domain coupling (pendulum k), this is detectable from
#'   amplitude envelope statistics regardless of spectral resolution.
#'   Bowl depth in chi²(K) scales as K², not K (signal-to-noise favourable).
#'
#' @param a1_fixed   Bifurcation parameter for region 1 (from Stage 1)
#' @param a2_fixed   Bifurcation parameter for region 2 (from Stage 1)
#' @param om1_fixed  Natural frequency for region 1, rad/TR (from Stage 1)
#' @param om2_fixed  Natural frequency for region 2, rad/TR (from Stage 1)
#' @return ODE function with signature (t, x, p)
#' @export
make_sl_paired_fixed_aw <- function(a1_fixed, om1_fixed, a2_fixed, om2_fixed) {
  force(a1_fixed); force(om1_fixed)
  force(a2_fixed); force(om2_fixed)

  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    K   <- pmax(SL_BOUNDS$K_MIN, pmin(SL_BOUNDS$K_MAX, p[1, ]))

    x1  <- x[1, ]
    y1  <- x[2, ]
    x2  <- x[3, ]
    y2  <- x[4, ]

    r1_sq <- pmin(x1^2 + y1^2, 25.0)
    r2_sq <- pmin(x2^2 + y2^2, 25.0)

    dx1 <- a1_fixed * x1 - om1_fixed * y1 - r1_sq * x1 + K * (x2 - x1)
    dy1 <- om1_fixed * x1 + a1_fixed * y1 - r1_sq * y1 + K * (y2 - y1)
    dx2 <- a2_fixed * x2 - om2_fixed * y2 - r2_sq * x2 + K * (x1 - x2)
    dy2 <- om2_fixed * x2 + a2_fixed * y2 - r2_sq * y2 + K * (y1 - y2)

    rbind(dx1, dy1, dx2, dy2)
  }
}


# =============================================================================
# PREPROCESSING HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
#' hilbert_analytic
#'
#' Compute the analytic signal via FFT-based Hilbert transform.
#' Returns a list with real (original) and imaginary (Hilbert) components.
#'
#' @param x  Numeric vector (bandpass-filtered BOLD signal, zero-mean)
#' @return   List: $real (= x), $imag (90°-shifted), $amplitude, $phase
#' @export
hilbert_analytic <- function(x) {
  n <- length(x)

  # If the input contains no finite samples, return NA vectors so callers
  # can detect and skip degenerate ROIs without producing NaNs from fft.
  if (!any(is.finite(x))) {
    na_vec <- rep(NA_real_, n)
    return(list(
      real = na_vec,
      imag = na_vec,
      amplitude = na_vec,
      phase = na_vec,
      amp_scale = NA_real_
    ))
  }

  # For robustness, require all samples to be finite. If there are some
  # non-finite entries, treat as degenerate and return NA vectors instead
  # of attempting an FFT on invalid data.
  if (any(!is.finite(x))) {
    na_vec <- rep(NA_real_, n)
    return(list(
      real = na_vec,
      imag = na_vec,
      amplitude = na_vec,
      phase = na_vec,
      amp_scale = NA_real_
    ))
  }

  X <- fft(x)

  # One-sided analytic signal construction (Marple 1999):
  h <- numeric(n)
  if (n %% 2 == 0) {
    h[1]       <- 1
    h[2:(n/2)] <- 2
    h[n/2 + 1] <- 1
  } else {
    h[1]            <- 1
    h[2:((n + 1)/2)] <- 2
  }

  Z  <- Re(fft(X * h, inverse = TRUE)) / n
  Zh <- Im(fft(X * h, inverse = TRUE)) / n

  # Normalise to unit mean amplitude so the SL cubic term |z|²z operates
  # in its intended regime (r ≈ 1).  amp_scale preserved for back-transform.
  amp_mean <- mean(sqrt(x^2 + Zh^2), na.rm = TRUE)
  if (!is.finite(amp_mean) || amp_mean <= 0) amp_mean <- .Machine$double.eps
  x_n  <- x / amp_mean
  Zh_n <- Zh / amp_mean

  list(
    real      = x_n,
    imag      = Zh_n,
    amplitude = sqrt(x_n^2 + Zh_n^2),
    phase     = atan2(Zh_n, x_n),
    amp_scale = amp_mean
  )
}


# -----------------------------------------------------------------------------
#' prepare_sl_single_input
#'
#' Build ts_data matrix for Stage 1 SL UKF from one ROI column.
#' Computes Hilbert transform to get imaginary component.
#'
#' @param smoothed_df  Data frame with Time column + ROI columns
#' @param roi_name     Column name of the ROI to use
#' @return             Matrix (T x 3): [time, x_BOLD, y_Hilbert]
#' @export
prepare_sl_single_input <- function(smoothed_df, roi_name) {
  stopifnot(roi_name %in% colnames(smoothed_df))
  time_col <- colnames(smoothed_df)[1]
  x_bold   <- smoothed_df[[roi_name]]
  ha       <- hilbert_analytic(x_bold)
  as.matrix(data.frame(
    time    = smoothed_df[[time_col]],
    x_bold  = ha$real,
    y_hilb  = ha$imag
  ))
}


# -----------------------------------------------------------------------------
#' prepare_sl_paired_input
#'
#' Build ts_data matrix for Stage 2 paired SL UKF from two ROI columns.
#' Both ROIs get the Hilbert treatment → 4 observed channels.
#'
#' @param smoothed_df  Data frame with Time column + ROI columns
#' @param roi1         Column name of ROI 1
#' @param roi2         Column name of ROI 2
#' @return             Matrix (T x 5): [time, x1, y1, x2, y2]
#' @export
prepare_sl_paired_input <- function(smoothed_df, roi1, roi2) {
  stopifnot(roi1 %in% colnames(smoothed_df),
            roi2 %in% colnames(smoothed_df))
  time_col <- colnames(smoothed_df)[1]
  ha1      <- hilbert_analytic(smoothed_df[[roi1]])
  ha2      <- hilbert_analytic(smoothed_df[[roi2]])
  as.matrix(data.frame(
    time = smoothed_df[[time_col]],
    x1   = ha1$real,  y1 = ha1$imag,
    x2   = ha2$real,  y2 = ha2$imag
  ))
}


# -----------------------------------------------------------------------------
#' omega_init_from_spectrum
#'
#' Estimate ω initialisation value from the dominant spectral peak of a BOLD
#' signal, clamped to the physiological range [OM_MIN, OM_MAX].
#'
#' @param x   Numeric vector (BOLD signal)
#' @param TR  Repetition time in seconds (default 2.0)
#' @return    ω_init in rad/TR
#' @export
# -----------------------------------------------------------------------------
#' omega_from_phase
#'
#' Estimate natural frequency ω from the mean instantaneous frequency of the
#' analytic signal.  More robust than spectral peak for 1/f BOLD signals, where
#' the dominant peak is always at the lowest BOLD band frequency (0.01 Hz).
#'
#' Method: unwrap the instantaneous phase from the Hilbert transform, compute
#' the mean phase increment per TR, then clamp to [OM_MIN, OM_MAX].
#'
#' @param ha   Output of hilbert_analytic() — list with $phase field
#' @param TR   Repetition time in seconds (default 2.0)
#' @return     ω estimate in rad/TR, clamped to physiological BOLD range
#' @export
omega_from_phase <- function(ha, TR = 2.0) {
  phase_unwrap <- function(ph) {
    # Simple unwrap: accumulate phase increments
    d <- diff(ph)
    d <- d - 2 * pi * round(d / (2 * pi))
    c(ph[1], ph[1] + cumsum(d))
  }

  ph_unwrapped <- phase_unwrap(ha$phase)
  # Mean phase increment per sample = mean instantaneous frequency (rad/TR)
  dph          <- diff(ph_unwrapped)
  om_est       <- median(dph[dph > 0])   # median of positive increments only

  if (is.na(om_est) || !is.finite(om_est))
    om_est <- SL_BOUNDS$OM_MIN   # fallback

  max(SL_BOUNDS$OM_MIN, min(SL_BOUNDS$OM_MAX, om_est))
}

#' omega_from_phase_v2
#' Alternative ω estimation from phase increments with outlier rejection.
#' Outliers can arise from noise or phase slips in the Hilbert transform, which
#' can bias the median increment.  This version excludes increments outside the
#' 10th–90th percentile range before taking the median.
#' In practice, this may yield more stable ω estimates for noisy BOLD signals.
#' Note: the original omega_from_phase() is simpler and may be sufficient in many cases; this is an optional refinement.
# -----------------------------------------------------------------------------
#' @param ha   Output of hilbert_analytic() — list with $phase field
#' @param TR   Repetition time in seconds (default 2.0)
#' @return     ω estimate in rad/TR, clamped to physiological BOLD range
#' @export
omega_from_phase_v2 <- function(ha, TR=2.0) {
  unwrap <- function(ph) { d<-diff(ph); d<-d-2*pi*round(d/(2*pi)); c(ph[1],ph[1]+cumsum(d)) }
  dph  <- diff(unwrap(ha$phase))
  dph_c <- dph[dph >= quantile(dph,0.10,na.rm=TRUE) & dph <= quantile(dph,0.90,na.rm=TRUE)]
  om   <- median(dph_c, na.rm=TRUE)
  if (is.na(om)||!is.finite(om)||om<=0) om <- SL_BOUNDS$OM_MIN
  max(SL_BOUNDS$OM_MIN, min(SL_BOUNDS$OM_MAX, om))
}


# -----------------------------------------------------------------------------
#' omega_init_from_spectrum
#'
#' Spectral-peak based ω initialisation (kept for reference/comparison).
#' NOTE: for 1/f BOLD signals this returns 0.01 Hz (lower bound) in ~90% of
#' ROIs because the dominant power is always at the lowest frequency bin.
#' Prefer omega_from_phase() for Stage 1 initialisation.
#'
#' @param x   Numeric vector (BOLD signal)
#' @param TR  Repetition time in seconds (default 2.0)
#' @return    ω_init in rad/TR
#' @export
omega_init_from_spectrum <- function(x, TR = 2.0) {
  n    <- length(x)
  sp   <- Mod(fft(x - mean(x)))[2:(n %/% 2 + 1)]^2
  freqs <- (1:(n %/% 2)) / (n * TR)

  idx_band <- which(freqs >= 0.01 & freqs <= 0.10)
  if (length(idx_band) == 0) idx_band <- seq_along(freqs)

  f_peak   <- freqs[idx_band[which.max(sp[idx_band])]]
  om_est   <- 2 * pi * f_peak * TR
  max(SL_BOUNDS$OM_MIN, min(SL_BOUNDS$OM_MAX, om_est))
}


# -----------------------------------------------------------------------------
#' make_sl_single_fixed_om
#'
#' Single-region Stuart-Landau with ω FIXED from instantaneous phase estimate.
#' Only bifurcation parameter a is estimated by the UKF (N_p = 1).
#'
#' WHY: In Stage 1 with N_p=2 [a, ω], ω is driven to the lower bound (0.01 Hz)
#' in ~90% of BOLD ROIs because BOLD has a 1/f spectrum — the UKF has no
#' gradient to move ω away from the boundary.  Fixing ω from the analytic
#' signal phase (omega_from_phase) removes that degree of freedom and gives
#' the UKF a well-conditioned 1D problem for a.
#'
#' EQUATIONS (same as make_sl_single, ω closed over):
#'   dx/dt = a·x − ω_fixed·y − (x²+y²)·x
#'   dy/dt = ω_fixed·x + a·y − (x²+y²)·y
#'
#' PARAMETERS: p[1,] = a   (N_p = 1)
#' STATE:       [x, y]     (N_y = 2)
#' FIXED:       om_fixed from omega_from_phase()
#'
#' @param om_fixed  Natural frequency in rad/TR (from omega_from_phase)
#' @return ODE function with signature (t, x, p)
#' @export
make_sl_single_fixed_om <- function(om_fixed) {
  force(om_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a   <- pmax(SL_BOUNDS$A_MIN, pmin(SL_BOUNDS$A_MAX, p[1, ]))
    om  <- om_fixed

    xr  <- x[1, ]
    xi  <- x[2, ]
    r2  <- pmin(xr^2 + xi^2, 25.0)

    dxr <- a * xr - om * xi - r2 * xr
    dxi <- om * xr + a * xi - r2 * xi

    rbind(dxr, dxi)
  }
}
