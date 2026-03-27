# =============================================================================
# wind_config.R
# Wind-specific configuration overriding BOLD-specific bounds from sl_models.R
#
# DOMAIN ADAPTATION NOTES:
#   This file is sourced AFTER sl_models.R to override SL_BOUNDS and SL_INIT
#   for daily wind speed time series. The SL model equations are unchanged —
#   only the parameter bounds and time constants are adapted.
#
# TIME UNITS:
#   In the BOLD pipeline, time is measured in TRs (1 TR = 2 s).
#   In the wind pipeline, time is measured in DAYS (1 step = 1 day).
#   Therefore ω is in rad/day, dT = 1.0 day, dt = 0.1 day.
#
# FREQUENCY BANDS:
#   Wind speed exhibits oscillatory structure at multiple timescales:
#     Annual cycle:     T = 365 d  →  ω = 2π/365 ≈ 0.0172 rad/day
#     Semi-annual:      T = 182 d  →  ω = 2π/182 ≈ 0.0345 rad/day
#     Monthly:          T ~30 d    →  ω ≈ 0.209 rad/day
#     Synoptic weather: T = 3-10 d →  ω = 0.63 – 2.09 rad/day
#   We set OM_MIN to capture the annual cycle and OM_MAX to capture
#   synoptic-scale weather fluctuations (3-day period).
#
# BIFURCATION PARAMETER:
#   a bounds remain at [-2, 2]. Wind oscillations are expected to be
#   near-critical or supercritical (a ≥ 0) for dominant seasonal modes,
#   unlike BOLD which was deeply subcritical (a ≈ -0.29 in MDD).
#
# COUPLING:
#   K bounds remain at [0, 5]. Two nearby stations (~50 km apart) should
#   exhibit detectable atmospheric coupling through shared large-scale
#   pressure systems modulated by local topography.
# =============================================================================


# --- Override SL bounds for wind (units: rad/day, dimensionless) -------------
# v2 FIX: bounds now match the SYNOPTIC band (3-15 day period), not the
# annual-to-synoptic range.  The annual cycle is externally forced and must
# be removed by bandpass filtering before the SL model sees the data.
# The OM bounds will be overridden per-band in the notebook when running
# intraseasonal or broadband analyses.
SL_BOUNDS <- list(
  A_MIN   = -2.0,                             # bifurcation lower bound
  A_MAX   =  2.0,                             # bifurcation upper bound
  OM_MIN  =  2 * pi / 15,                     # 15-day period ≈ 0.419 rad/day
  OM_MAX  =  2 * pi / 3,                      # 3-day period  ≈ 2.094 rad/day
  K_MIN   =  0.0,
  K_MAX   =  5.0
)

# --- Override SL initial guesses for wind ------------------------------------
SL_INIT <- list(
  A_INIT  = 0.0,                              # start at criticality
  OM_INIT = 2 * pi / 7,                       # ~7-day period (mid-synoptic)
  K_INIT  = 0.1
)

# --- Wind-specific UKF tuning -----------------------------------------------
WIND_UKF <- list(
  dT       = 1.0,           # 1 day between observations
  dt       = 0.1,           # RK4 sub-step = 0.1 day (10 sub-steps per day)
  R_scale  = 0.3,           # observation noise scale (start same as BOLD)
  Q_scale  = 0.015,         # process noise scale (start same as BOLD)
  t_dummy  = 0.0            # scalar dummy time (models are autonomous)
)

# --- Coupled pendulum bounds for wind (LHO / gLk models) --------------------
# gamma = squared natural frequency (rad²/day²)
# For annual cycle: gamma = (2π/365)² ≈ 2.96e-4
# For synoptic:     gamma = (2π/3)²   ≈ 4.39
PENDULUM_BOUNDS <- list(
  GAMMA_MIN = 1e-4,          # below annual cycle
  GAMMA_MAX = 5.0,           # above synoptic scale
  K_MIN     = 0.0,
  K_MAX     = 5.0
)

# --- Seasonal segmentation --------------------------------------------------
SEASONS <- list(
  winter = c(12, 1, 2),      # Dec-Jan-Feb
  spring = c(3, 4, 5),       # Mar-Apr-May
  summer = c(6, 7, 8),       # Jun-Jul-Aug
  autumn = c(9, 10, 11)      # Sep-Oct-Nov
)

# --- Station metadata --------------------------------------------------------
STATIONS <- list(
  aralik = list(
    code    = 18195,
    name    = "Aralik",
    prefix  = "a"
  ),
  igdir  = list(
    code    = 17763,
    name    = "Igdir Havalimani",
    prefix  = "h"
  )
)

# --- Signal types to analyse -------------------------------------------------
SIGNAL_TYPES <- c("avg_spd", "max_spd")

cat("[wind_config] SL_BOUNDS overridden for wind timescales.\n")
cat(sprintf("  OM range: [%.4f, %.4f] rad/day (period: [%.1f, %.1f] days)\n",
            SL_BOUNDS$OM_MIN, SL_BOUNDS$OM_MAX,
            2*pi/SL_BOUNDS$OM_MAX, 2*pi/SL_BOUNDS$OM_MIN))
cat(sprintf("  A  range: [%.1f, %.1f]\n", SL_BOUNDS$A_MIN, SL_BOUNDS$A_MAX))
cat(sprintf("  K  range: [%.1f, %.1f]\n", SL_BOUNDS$K_MIN, SL_BOUNDS$K_MAX))
cat(sprintf("  dT = %.1f day, dt = %.2f day\n", WIND_UKF$dT, WIND_UKF$dt))
