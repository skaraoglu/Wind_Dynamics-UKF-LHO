# =============================================================================
# models.R
# Consolidated ODE models for UKF parameter estimation.
#
# Contains ALL models used across chapters 3-5:
#   - LHO two-stage pipeline (PRIMARY for wind — Chapter 5)
#   - LHO comparison models (symmetric/asymmetric coupled)
#   - Nonlinear pendulum models (Chapter 4 legacy, for comparison)
#   - Stuart-Landau models (Chapter 4 primary, secondary for wind)
#
# All models follow the interface:
#   ode_model(t, x, p)  ->  matrix of first derivatives (same dims as x)
# where:
#   t  : scalar dummy time (models have no explicit time dependence)
#   x  : (N_y x N_sigma) matrix of state variables (or (N_y x 1) vector)
#   p  : (N_p x N_sigma) matrix of parameters     (or (N_p x 1) vector)
#
# STATE LAYOUT for paired second-order models:
#   x = [theta1, theta1_dot, theta2, theta2_dot]  (N_y = 4)
#   Returns d/dt [theta1, theta1_dd, theta2, theta2_dd]
#            = [theta1_dot, accel1, theta2_dot, accel2]
#
# STATE LAYOUT for single second-order models:
#   x = [theta, theta_dot]  (N_y = 2)
#   Returns d/dt [theta, theta_dd] = [theta_dot, accel]
#
# STATE LAYOUT for SL models (complex analytic signal):
#   x = [x_real, y_imag]  (N_y = 2, single) or [x1,y1,x2,y2] (N_y = 4, paired)
# =============================================================================


# #############################################################################
#
#   SECTION 1: LHO TWO-STAGE PIPELINE (PRIMARY FOR WIND)
#
#   Stage 1: make_lho_single()    → per-station gamma (N_p=1)
#   Stage 2: make_lho_fixed_ab()  → inter-station k   (N_p=1)
#
#   Advantages over SL for wind:
#     - No Hilbert transform required (works on [position, velocity])
#     - Linear dynamics: UKF sigma-point propagation is EXACT
#     - N_p=1 per stage: strictly 1D optimisation, no drift
#     - Physics match: damped harmonic oscillator is correct model
#       for wind speed perturbations around the detrended mean
#
# #############################################################################


# -----------------------------------------------------------------------------
# make_lho_single  (Stage 1 — single-station)
#
# Fits only gamma to a single station's wind speed time series.
# No coupling term. Simplest possible identifiable model.
#
# EQUATIONS:
#   theta_dd = -gamma * theta - 2*zeta*sqrt(gamma) * theta_dot
#
# PARAMETERS: p[1,] = gamma  (squared natural frequency, rad^2/day^2)
# STATE:      [theta, theta_dot]  (N_y = 2)
# FIXED:      zeta (damping ratio, closed over)
#
# ts_data: [time, theta, theta_dot]  (3 columns)
#   theta     = detrended, z-scored wind speed
#   theta_dot = diff(theta) / dT  (numerical velocity)
#
# @param zeta_fixed  Damping ratio (default 0.1, lightly underdamped)
# @return ODE function (t, x, p)
# @export
make_lho_single <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    zeta       <- zeta_fixed
    theta      <- x[1, ]
    theta_dot  <- x[2, ]

    c_damp     <- 2 * zeta * sqrt(pmax(gamma, 0))
    theta_ddot <- -gamma * theta - c_damp * theta_dot

    rbind(theta_dot, theta_ddot)
  }
}

lho_single_model <- make_lho_single(0.1)


# -----------------------------------------------------------------------------
# make_lho_fixed_ab  (Stage 2 — coupled, frequencies fixed from Stage 1)
#
# Both natural frequencies a (station 1) and b (station 2) are fixed from
# Stage 1 estimates.  Only coupling k is free.  This gives a strictly 1D
# chi-sq surface in k, which is necessary and sufficient for identifiability
# via the normal-mode splitting mechanism.
#
# EQUATIONS:
#   theta1_dd = -a*theta1 - 2*zeta*sqrt(a)*theta1_dot - k*(theta1 - theta2)
#   theta2_dd = -b*theta2 - 2*zeta*sqrt(b)*theta2_dot + k*(theta1 - theta2)
#
# PARAMETERS: p[1,] = k  (coupling strength, rad^2/day^2)
# STATE:      [theta1, theta1_dot, theta2, theta2_dot]  (N_y = 4)
# FIXED:      a, b (from Stage 1), zeta (closed over)
#
# ts_data: [time, theta1, theta1_dot, theta2, theta2_dot]  (5 columns)
#
# Normal mode frequencies (for a=b=gamma):
#   omega_minus = sqrt(gamma)          (in-phase mode)
#   omega_plus  = sqrt(gamma + 2*k)    (out-of-phase mode)
#   k is recovered from mode splitting: k = (omega_plus^2 - omega_minus^2) / 2
#
# @param a_fixed     Squared frequency for station 1 (from Stage 1)
# @param b_fixed     Squared frequency for station 2 (from Stage 1)
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
# @export
make_lho_fixed_ab <- function(a_fixed, b_fixed, zeta_fixed = 0.1) {
  force(a_fixed); force(b_fixed); force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    k          <- p[1, ]
    a          <- a_fixed
    b          <- b_fixed
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(max(a, 0))
    cb <- 2 * zeta * sqrt(max(b, 0))

    theta1_ddot <- -a * theta1 - ca * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -b * theta2 - cb * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}


# #############################################################################
#
#   SECTION 2: LHO COMPARISON MODELS
#
#   Used for model-selection diagnostics (Dchi2 tests) and sensitivity
#   analyses. Not the primary analysis pipeline.
#
# #############################################################################


# -----------------------------------------------------------------------------
# make_lho_fixed_zeta  (symmetric coupled LHO, gamma shared)
#
# Both stations share gamma (symmetric assumption a = b = gamma).
# Coupling k estimated jointly with gamma.  N_p = 2.
#
# EQUATIONS:
#   theta1_dd = -gamma*theta1 - 2*zeta*sqrt(gamma)*theta1_dot - k*(theta1 - theta2)
#   theta2_dd = -gamma*theta2 - 2*zeta*sqrt(gamma)*theta2_dot + k*(theta1 - theta2)
#
# PARAMETERS: p = [gamma, k]  (N_p = 2)
# STATE:      [theta1, theta1_dot, theta2, theta2_dot]  (N_y = 4)
#
# IDENTIFIABILITY:
#   Imposing a = b = gamma breaks the structural ridge (a-d, b-d, k+d).
#   gamma and k are recovered from the two normal mode frequencies:
#     gamma = omega_minus^2,  k = (omega_plus^2 - omega_minus^2) / 2
#
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
# @export
make_lho_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    k          <- p[2, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    c_damp <- 2 * zeta * sqrt(pmax(gamma, 0))

    theta1_ddot <- -gamma * theta1 - c_damp * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -gamma * theta2 - c_damp * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

coupled_osc_model_lho_fz <- make_lho_fixed_zeta(0.1)


# -----------------------------------------------------------------------------
# make_lho_abk_fixed_zeta  (asymmetric coupled LHO)
#
# Allows different natural frequencies per station: a != b.
# Used for Dchi2 comparison: if chi2_lho - chi2_abk > threshold, the
# symmetric assumption is violated.
#
# PARAMETERS: p = [a, b, k]  (N_p = 3)
# WARNING: Has identifiability ridge (a-d, b-d, k+d).
#          Use only for model-selection, not production coupling estimation.
#
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
# @export
make_lho_abk_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a          <- p[1, ]
    b          <- p[2, ]
    k          <- p[3, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(pmax(a, 0))
    cb <- 2 * zeta * sqrt(pmax(b, 0))

    theta1_ddot <- -a * theta1 - ca * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -b * theta2 - cb * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

coupled_osc_model_lho_abk_fz <- make_lho_abk_fixed_zeta(0.1)


# #############################################################################
#
#   SECTION 3: NONLINEAR PENDULUM MODELS (Chapter 4 legacy)
#
#   Included for completeness and cross-model comparison.
#   Identical physics to LHO but with sin(theta) nonlinearity retained.
#   For z-scored data, sin(theta) ~ theta to 5 significant figures,
#   so these should give nearly identical results to the LHO variants.
#
# #############################################################################


# -----------------------------------------------------------------------------
# make_glk2_fixed_zeta  (nonlinear symmetric pendulum, zeta fixed)
#
# EQUATIONS:
#   theta1_dd = -gamma*sin(theta1) - 2*zeta*sqrt(gamma)*theta1_dot
#               - k*(sin(theta1) - sin(theta2))
#   theta2_dd = -gamma*sin(theta2) - 2*zeta*sqrt(gamma)*theta2_dot
#               + k*(sin(theta1) - sin(theta2))
#
# PARAMETERS: p = [gamma, k]  (N_p = 2)
#
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
# @export
make_glk2_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    k          <- p[2, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    c_damp <- 2 * zeta * sqrt(pmax(gamma, 0))

    theta1_ddot <- -gamma * sin(theta1) - c_damp * theta1_dot - k * (sin(theta1) - sin(theta2))
    theta2_ddot <- -gamma * sin(theta2) - c_damp * theta2_dot + k * (sin(theta1) - sin(theta2))

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

coupled_osc_model_glk2_fz <- make_glk2_fixed_zeta(0.1)


# -----------------------------------------------------------------------------
# make_abk2_fixed_zeta  (nonlinear asymmetric pendulum, zeta fixed)
#
# PARAMETERS: p = [a, b, k]  (N_p = 3)
#
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
# @export
make_abk2_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a          <- p[1, ]
    b          <- p[2, ]
    k          <- p[3, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(pmax(a, 0))
    cb <- 2 * zeta * sqrt(pmax(b, 0))

    theta1_ddot <- -a * sin(theta1) - ca * theta1_dot - k * (sin(theta1) - sin(theta2))
    theta2_ddot <- -b * sin(theta2) - cb * theta2_dot + k * (sin(theta1) - sin(theta2))

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

coupled_osc_model_abk2_fz <- make_abk2_fixed_zeta(0.1)


# #############################################################################
#
#   SECTION 4: STUART-LANDAU MODELS (Chapter 4 primary)
#
#   The SL model is the normal form of a supercritical Hopf bifurcation.
#   Requires Hilbert-transformed analytic signal as input.
#   Included for cross-domain comparison with Chapter 4 BOLD results.
#
#   For wind data, the SL is a SECONDARY analysis — the LHO pipeline is
#   the primary model because wind physics match damped oscillation, not
#   Hopf bifurcation. The SL is retained to enable direct comparison of
#   bifurcation parameter 'a' between wind and BOLD MDD (a ~ -0.29).
#
# #############################################################################

# NOTE: SL_BOUNDS is defined in wind_config.R (or sl_models.R for BOLD).
# It must be loaded before these functions are called.

# -----------------------------------------------------------------------------
# make_sl_single  (single-station SL, free a and omega)
#
# EQUATIONS:
#   dx/dt = a*x - omega*y - (x^2+y^2)*x
#   dy/dt = omega*x + a*y - (x^2+y^2)*y
#
# PARAMETERS: p = [a, omega]  (N_p = 2)
# STATE:      [x_real, y_hilbert]  (N_y = 2)
# INPUT:      Hilbert analytic signal (bandpass-filtered)
#
# @return ODE function (t, x, p)
# @export
make_sl_single <- function() {
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a   <- pmax(SL_BOUNDS$A_MIN,  pmin(SL_BOUNDS$A_MAX,  p[1, ]))
    om  <- pmax(SL_BOUNDS$OM_MIN, pmin(SL_BOUNDS$OM_MAX, p[2, ]))

    xr  <- x[1, ]
    xi  <- x[2, ]
    r2  <- pmin(xr^2 + xi^2, 25.0)

    dxr <- a * xr - om * xi - r2 * xr
    dxi <- om * xr + a * xi - r2 * xi

    rbind(dxr, dxi)
  }
}

sl_single_model <- make_sl_single()


# -----------------------------------------------------------------------------
# make_sl_single_fixed_om  (single-station SL, omega fixed from phase estimate)
#
# PARAMETERS: p = [a]  (N_p = 1)
# FIXED:      omega (closed over from omega_from_phase)
#
# @param om_fixed  Natural frequency in rad/day (from omega_from_phase)
# @return ODE function (t, x, p)
# @export
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


# -----------------------------------------------------------------------------
# make_sl_paired_fixed_aw  (paired SL, a and omega fixed per station)
#
# PARAMETERS: p = [K]  (N_p = 1, diffusive coupling)
# STATE:      [x1, y1, x2, y2]  (N_y = 4)
# FIXED:      a1, om1, a2, om2 from Stage 1
#
# @param a1_fixed, om1_fixed  Station 1 bifurcation and frequency
# @param a2_fixed, om2_fixed  Station 2 bifurcation and frequency
# @return ODE function (t, x, p)
# @export
make_sl_paired_fixed_aw <- function(a1_fixed, om1_fixed, a2_fixed, om2_fixed) {
  force(a1_fixed); force(om1_fixed)
  force(a2_fixed); force(om2_fixed)

  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    K   <- pmax(SL_BOUNDS$K_MIN, pmin(SL_BOUNDS$K_MAX, p[1, ]))

    x1  <- x[1, ]; y1  <- x[2, ]
    x2  <- x[3, ]; y2  <- x[4, ]

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
# MODEL INVENTORY SUMMARY
# =============================================================================
# Section 1 — LHO Two-Stage Pipeline (PRIMARY):
#   make_lho_single(zeta)              N_p=1 [gamma]     N_y=2  Stage 1
#   make_lho_fixed_ab(a, b, zeta)      N_p=1 [k]         N_y=4  Stage 2
#
# Section 2 — LHO Comparison:
#   make_lho_fixed_zeta(zeta)          N_p=2 [gamma, k]  N_y=4  symmetric
#   make_lho_abk_fixed_zeta(zeta)      N_p=3 [a, b, k]   N_y=4  asymmetric
#
# Section 3 — Nonlinear Pendulum (legacy):
#   make_glk2_fixed_zeta(zeta)         N_p=2 [gamma, k]  N_y=4  nonlinear sym
#   make_abk2_fixed_zeta(zeta)         N_p=3 [a, b, k]   N_y=4  nonlinear asym
#
# Section 4 — Stuart-Landau (Ch.4 cross-domain):
#   make_sl_single()                   N_p=2 [a, omega]   N_y=2  Hilbert input
#   make_sl_single_fixed_om(om)        N_p=1 [a]          N_y=2  Hilbert input
#   make_sl_paired_fixed_aw(a1,o1,a2,o2) N_p=1 [K]        N_y=4  Hilbert input
# =============================================================================

cat("[models.R] Loaded: LHO pipeline (make_lho_single, make_lho_fixed_ab),\n")
cat("           LHO comparison (lho_fz, lho_abk_fz),\n")
cat("           Nonlinear pendulum (glk2_fz, abk2_fz),\n")
cat("           Stuart-Landau (sl_single, sl_single_fixed_om, sl_paired_fixed_aw)\n")
