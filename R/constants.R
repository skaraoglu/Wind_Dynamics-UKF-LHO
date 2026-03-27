# =============================================================================
# constants.R
# Core UKF numerical constants. Domain-independent.
# =============================================================================

UKF_CONSTANTS <- list(
  # --- Cholesky / positive-definiteness stabilisation ---
  JITTER_INIT    = 1e-8,       # initial jitter added to covariance matrix
  JITTER_MAX     = 1e-2,       # maximum jitter before falling back to nearPD
  EIGVAL_MIN     = 1e-12,      # minimum eigenvalue for positive-definiteness
  COND_NUM_MAX   = 1e10,       # maximum condition number

  # --- Parameter estimation ---
  PARAM_MIN      = 1e-8,       # minimum parameter value when forcePositive = TRUE
  PARAM_TOL_DEFAULT = 1e-3,    # default L2 convergence tolerance
  MAXSTEPS_DEFAULT  = 1000,    # default maximum iterations
  CHISQ_PLATEAU_TOL = 1e-6,   # chi-sq plateau detection tolerance


  # --- Numerical safety ---
  STIFFNESS_WARN = 100         # warn when |parameter| exceeds this value
)
