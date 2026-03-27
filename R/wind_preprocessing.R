# =============================================================================
# wind_preprocessing.R  (v2 — BANDPASS FIX)
#
# CRITICAL FIX (v2):
#   v1 omitted band-pass filtering.  hilbert_analytic() requires NARROWBAND
#   input (docstring: "@param x Numeric vector (bandpass-filtered BOLD signal)").
#   The BOLD pipeline worked because AFNI pre-applied a 0.01-0.10 Hz bandpass.
#   For wind, we must replicate this step explicitly.
# =============================================================================

# --- BAND-PASS FILTER (FFT-based, zero-phase) --------------------------------
bandpass_fft <- function(x, dT = 1.0, f_low = 1/15, f_high = 1/3) {
  n <- length(x)
  if (n < 10) return(x)
  mu <- mean(x, na.rm = TRUE)
  x_c <- x - mu
  na_mask <- is.na(x_c)
  x_c[na_mask] <- 0
  X <- fft(x_c)
  freqs <- (0:(n-1)) / (n * dT)
  freqs[freqs > 1/(2*dT)] <- freqs[freqs > 1/(2*dT)] - 1/dT
  freqs <- abs(freqs)
  mask <- (freqs >= f_low) & (freqs <= f_high)
  mask[1] <- FALSE
  X[!mask] <- 0 + 0i
  x_f <- Re(fft(X, inverse = TRUE)) / n
  x_f[na_mask] <- NA
  x_f
}

bandpass_quality_check <- function(x_raw, x_filtered) {
  x_raw <- x_raw[!is.na(x_raw)]
  x_filtered <- x_filtered[!is.na(x_filtered)]
  p_raw <- sum(x_raw^2)
  if (p_raw == 0) return(0)
  sum(x_filtered^2) / p_raw
}

# --- PREDEFINED FREQUENCY BANDS -----------------------------------------------
WIND_BANDS <- list(
  synoptic = list(
    name = "synoptic", f_low = 1/15, f_high = 1/3,
    label = "3-15 day (synoptic weather)",
    om_min = 2 * pi / 15, om_max = 2 * pi / 3
  ),
  intraseasonal = list(
    name = "intraseasonal", f_low = 1/60, f_high = 1/15,
    label = "15-60 day (intraseasonal)",
    om_min = 2 * pi / 60, om_max = 2 * pi / 15
  ),
  broadband = list(
    name = "broadband", f_low = 1/60, f_high = 1/3,
    label = "3-60 day (combined)",
    om_min = 2 * pi / 60, om_max = 2 * pi / 3
  )
)

# --- DATA LOADING --------------------------------------------------------------
load_wind_data <- function(path = "data/wind.csv") {
  stopifnot(file.exists(path))
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (grepl("/", df$date[1])) {
    df$date <- as.Date(df$date, format = "%d/%m/%Y")
  } else {
    df$date <- as.Date(df$date)
  }
  df <- df[order(df$date), ]
  df$day_index   <- seq_len(nrow(df))
  df$year        <- as.integer(format(df$date, "%Y"))
  df$month       <- as.integer(format(df$date, "%m"))
  df$day_of_year <- as.integer(format(df$date, "%j"))
  cat(sprintf("[preprocessing] Loaded %d rows from %s to %s\n",
              nrow(df), min(df$date), max(df$date)))
  df
}

# --- MISSING DATA ---------------------------------------------------------------
summarise_missing <- function(df) {
  cols <- c("a_max_dir","a_max_spd","a_avg_dir","a_avg_spd",
            "h_max_dir","h_max_spd","h_avg_dir","h_avg_spd")
  cols <- intersect(cols, colnames(df))
  cat("\n--- Missing Value Summary ---\n")
  miss_df <- data.frame(
    column = cols, n_total = nrow(df),
    n_missing = sapply(cols, function(c) sum(is.na(df[[c]]))),
    pct_missing = sapply(cols, function(c) round(100*mean(is.na(df[[c]])),1)),
    stringsAsFactors = FALSE)
  print(miss_df, row.names = FALSE)
  invisible(miss_df)
}

interpolate_gaps <- function(x, max_gap_days = 7) {
  n <- length(x)
  if (!any(is.na(x))) return(x)
  nas <- is.na(x); runs <- rle(nas); idx <- 1
  for (i in seq_along(runs$lengths)) {
    gap_end <- idx + runs$lengths[i] - 1
    if (runs$values[i] && runs$lengths[i] <= max_gap_days) {
      left  <- if (idx > 1) x[idx - 1] else NA
      right <- if (gap_end < n) x[gap_end + 1] else NA
      if (!is.na(left) && !is.na(right)) {
        x[idx:gap_end] <- seq(left, right, length.out = runs$lengths[i]+2)[2:(runs$lengths[i]+1)]
      } else if (!is.na(left)) { x[idx:gap_end] <- left
      } else if (!is.na(right)) { x[idx:gap_end] <- right }
    }
    idx <- gap_end + 1
  }
  x
}

# --- FULL PIPELINE: detrend -> z-score -> BANDPASS -> ready for Hilbert --------
preprocess_wind_signal <- function(x, max_gap_days = 7, band = "synoptic") {
  raw <- x; n <- length(x)
  n_na_before <- sum(is.na(x))
  x <- interpolate_gaps(x, max_gap_days)
  n_interpolated <- n_na_before - sum(is.na(x))
  valid <- !is.na(x)
  if (sum(valid) < 10) {
    warning("fewer than 10 valid points")
    return(list(signal=x, signal_unfiltered=x, raw=raw, trend_coefs=c(0,0),
                mean=0, sd=1, band_power_frac=NA, n_interpolated=0,
                n_remaining_na=sum(is.na(x)), band_info=NULL))
  }
  t_vec <- seq_len(n)
  fit_df <- data.frame(y = x[valid], t = t_vec[valid])
  fit <- lm(y ~ t, data = fit_df)
  trend_at_valid <- predict(fit, newdata = data.frame(t = t_vec[valid]))
  x_dt <- x; x_dt[valid] <- x[valid] - trend_at_valid
  mu <- mean(x_dt[valid], na.rm = TRUE); sigma <- sd(x_dt[valid], na.rm = TRUE)
  if (sigma < 1e-10) sigma <- 1
  x_z <- x_dt; x_z[valid] <- (x_dt[valid] - mu) / sigma

  # *** BANDPASS FILTER (if band is specified) ***
  # band = NULL: skip filtering (for LHO path, which works on raw signal)
  # band = "synoptic" etc: apply bandpass (for SL path, which needs Hilbert)
  if (!is.null(band)) {
    band_info <- WIND_BANDS[[band]]
    if (is.null(band_info)) { band_info <- WIND_BANDS$synoptic }
    x_f <- bandpass_fft(x_z, dT = 1.0, f_low = band_info$f_low, f_high = band_info$f_high)
    bpf <- bandpass_quality_check(x_z[valid], x_f[valid])
  } else {
    x_f <- x_z        # no filtering — return z-scored signal directly
    band_info <- NULL
    bpf <- 1.0        # all power retained
  }

  list(signal = x_f, signal_unfiltered = x_z, raw = raw,
       trend_coefs = coef(fit), mean = mu, sd = sigma,
       band_power_frac = bpf, n_interpolated = n_interpolated,
       n_remaining_na = sum(is.na(x_f)), band_info = band_info)
}

# --- ts_data BUILDERS (now receive bandpass-filtered input) --------------------
prepare_wind_single_input <- function(signal, day_index) {
  valid <- !is.na(signal)
  if (sum(valid) < length(signal)) signal[!valid] <- 0
  ha <- hilbert_analytic(signal)
  as.matrix(data.frame(time = day_index, x_real = ha$real, y_hilb = ha$imag))
}

prepare_wind_paired_input <- function(signal1, signal2, day_index) {
  signal1[is.na(signal1)] <- 0; signal2[is.na(signal2)] <- 0
  ha1 <- hilbert_analytic(signal1); ha2 <- hilbert_analytic(signal2)
  as.matrix(data.frame(time = day_index, x1 = ha1$real, y1 = ha1$imag,
                        x2 = ha2$real, y2 = ha2$imag))
}

# --- SEGMENTATION ---------------------------------------------------------------
segment_by_season <- function(df) {
  segments <- list()
  for (yr in unique(df$year)) {
    for (sname in names(SEASONS)) {
      months <- SEASONS[[sname]]
      sub <- df[df$year == yr & df$month %in% months, ]
      if (nrow(sub) >= 30) {
        label <- sprintf("%s_%d", sname, yr)
        sub$day_index <- seq_len(nrow(sub))
        segments[[label]] <- sub
      }
    }
  }
  cat(sprintf("[preprocessing] Created %d seasonal segments\n", length(segments)))
  segments
}

segment_by_year <- function(df) {
  segments <- list()
  for (yr in unique(df$year)) {
    sub <- df[df$year == yr, ]
    sub$day_index <- seq_len(nrow(sub))
    segments[[as.character(yr)]] <- sub
  }
  cat(sprintf("[preprocessing] Created %d yearly segments\n", length(segments)))
  segments
}

# --- OMEGA ESTIMATION -----------------------------------------------------------
omega_from_spectrum_wind <- function(x, dT = 1.0) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 20) return(SL_BOUNDS$OM_MIN)
  sp <- Mod(fft(x - mean(x)))[2:(n %/% 2 + 1)]^2
  freqs <- (1:(n %/% 2)) / (n * dT)
  f_peak <- freqs[which.max(sp)]
  om_est <- 2 * pi * f_peak * dT
  max(SL_BOUNDS$OM_MIN, min(SL_BOUNDS$OM_MAX, om_est))
}
