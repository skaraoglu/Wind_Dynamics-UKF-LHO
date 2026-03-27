# =============================================================================
# wind_logging.R
# Detailed logging system for the wind dynamics experiment.
#
# Provides a structured logging framework that writes per-fit results,
# summary statistics, and diagnostic information to CSV and text logs.
# Mirrors the verbose logging from the BOLD experiment pipeline.
# =============================================================================


# --- Initialise global log storage -------------------------------------------
WIND_LOG <- new.env(parent = emptyenv())
WIND_LOG$log_dir     <- "logs"
WIND_LOG$start_time  <- Sys.time()
WIND_LOG$entries     <- list()
WIND_LOG$entry_count <- 0L
WIND_LOG$master_log  <- NULL


# -----------------------------------------------------------------------------
#' init_logging
#'
#' Initialise the logging system. Creates log directory, opens master log CSV.
#'
#' @param log_dir  Directory for log files (default: "logs")
#' @return Invisible NULL
#' @export
init_logging <- function(log_dir = "logs") {
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  WIND_LOG$log_dir    <- log_dir
  WIND_LOG$start_time <- Sys.time()
  WIND_LOG$entries    <- list()
  WIND_LOG$entry_count <- 0L

  # Master log CSV: one row per fit
  master_path <- file.path(log_dir, "master_log.csv")
  WIND_LOG$master_log <- master_path

  header <- data.frame(
    timestamp     = character(),
    experiment    = character(),
    hypothesis    = character(),
    station       = character(),
    signal_type   = character(),
    segment       = character(),
    model         = character(),
    N_p           = integer(),
    N_y           = integer(),
    T_obs         = integer(),
    param_names   = character(),
    param_values  = character(),
    chisq         = numeric(),
    chisq_best    = numeric(),
    iterations    = integer(),
    param_norm    = numeric(),
    converged     = logical(),
    elapsed_sec   = numeric(),
    notes         = character(),
    stringsAsFactors = FALSE
  )
  write.csv(header, master_path, row.names = FALSE)

  # Human-readable session log
  session_path <- file.path(log_dir, "session_log.txt")
  cat(sprintf("═══════════════════════════════════════════════════════════════\n"),
      file = session_path)
  cat(sprintf("  WIND DYNAMICS EXPERIMENT LOG\n"), file = session_path, append = TRUE)
  cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = session_path, append = TRUE)
  cat(sprintf("═══════════════════════════════════════════════════════════════\n\n"),
      file = session_path, append = TRUE)

  cat(sprintf("[logging] Initialised. Master log: %s\n", master_path))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' log_fit
#'
#' Log a single UKF fit result to the master CSV and session log.
#'
#' @param result       Output of iterative_param_optim() or similar
#' @param experiment   String label (e.g., "H1_single_station")
#' @param hypothesis   String (e.g., "H1")
#' @param station      Station name (e.g., "Aralik")
#' @param signal_type  Signal type (e.g., "avg_spd")
#' @param segment      Segment label (e.g., "full", "winter_2016")
#' @param model        Model name (e.g., "SL_single", "SL_paired", "LHO")
#' @param N_p          Number of parameters
#' @param N_y          Number of observed states
#' @param T_obs        Number of time points
#' @param param_names  Character vector of parameter names
#' @param elapsed_sec  Elapsed time in seconds
#' @param notes        Optional notes string
#' @return Invisible NULL
#' @export
log_fit <- function(result, experiment, hypothesis, station, signal_type,
                    segment, model, N_p, N_y, T_obs, param_names,
                    elapsed_sec = NA, notes = "") {

  WIND_LOG$entry_count <- WIND_LOG$entry_count + 1L

  # Extract param values
  pvals <- if (!is.null(result$par)) result$par else result$param_est
  pvals_str  <- paste(round(as.numeric(pvals), 6), collapse = ";")
  pnames_str <- paste(param_names, collapse = ";")

  # Convergence check
  converged <- !is.null(result$param_norm) && result$param_norm < 1e-3

  row <- data.frame(
    timestamp     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    experiment    = experiment,
    hypothesis    = hypothesis,
    station       = station,
    signal_type   = signal_type,
    segment       = segment,
    model         = model,
    N_p           = N_p,
    N_y           = N_y,
    T_obs         = T_obs,
    param_names   = pnames_str,
    param_values  = pvals_str,
    chisq         = if (!is.null(result$chisq)) tail(result$chisq, 1) else NA,
    chisq_best    = if (!is.null(result$value)) result$value else NA,
    iterations    = if (!is.null(result$steps)) result$steps else NA,
    param_norm    = if (!is.null(result$param_norm)) result$param_norm else NA,
    converged     = converged,
    elapsed_sec   = elapsed_sec,
    notes         = notes,
    stringsAsFactors = FALSE
  )

  # Append to master CSV
  write.table(row, WIND_LOG$master_log,
              sep = ",", row.names = FALSE, col.names = FALSE,
              quote = TRUE, append = TRUE)

  # Append to in-memory list
  WIND_LOG$entries[[WIND_LOG$entry_count]] <- row

  # Human-readable session log
  session_path <- file.path(WIND_LOG$log_dir, "session_log.txt")
  cat(sprintf("\n── Fit #%d [%s] ──────────────────────────────────────────\n",
              WIND_LOG$entry_count, hypothesis),
      file = session_path, append = TRUE)
  cat(sprintf("  Experiment : %s\n", experiment), file = session_path, append = TRUE)
  cat(sprintf("  Station    : %s | Signal: %s | Segment: %s\n",
              station, signal_type, segment), file = session_path, append = TRUE)
  cat(sprintf("  Model      : %s  (N_p=%d, N_y=%d, T=%d)\n",
              model, N_p, N_y, T_obs), file = session_path, append = TRUE)
  for (i in seq_along(param_names)) {
    cat(sprintf("  %-12s: %g\n", param_names[i], as.numeric(pvals)[i]),
        file = session_path, append = TRUE)
  }
  cat(sprintf("  chi²       : %g  (best: %g)\n",
              row$chisq, row$chisq_best), file = session_path, append = TRUE)
  cat(sprintf("  Iterations : %s | Converged: %s | Time: %.1f s\n",
              row$iterations, row$converged, ifelse(is.na(elapsed_sec), 0, elapsed_sec)),
      file = session_path, append = TRUE)
  if (nchar(notes) > 0) {
    cat(sprintf("  Notes      : %s\n", notes), file = session_path, append = TRUE)
  }

  # Console output
  cat(sprintf("[%s] %s | %s/%s/%s | params=[%s] | chi²=%.4f | %d iters | %.1fs\n",
              hypothesis, model, station, signal_type, segment,
              pvals_str, row$chisq_best, row$iterations,
              ifelse(is.na(elapsed_sec), 0, elapsed_sec)))

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' log_section
#'
#' Write a section header to the session log and console.
#'
#' @param title  Section title string
#' @export
log_section <- function(title) {
  session_path <- file.path(WIND_LOG$log_dir, "session_log.txt")
  header <- sprintf("\n\n%s\n  %s\n%s\n",
                    strrep("═", 65), title, strrep("═", 65))
  cat(header, file = session_path, append = TRUE)
  cat(header)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' log_message
#'
#' Write an informational message to session log and console.
#'
#' @param msg  Message string
#' @export
log_message <- function(msg) {
  session_path <- file.path(WIND_LOG$log_dir, "session_log.txt")
  line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, file = session_path, append = TRUE)
  cat(line)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' log_chi2_surface
#'
#' Log a chi² surface scan for coupling identifiability analysis.
#'
#' @param K_grid      Numeric vector of K values scanned
#' @param chisq_vals  Corresponding chi² values
#' @param label       Descriptive label
#' @param hypothesis  Hypothesis tag (e.g., "H3")
#' @export
log_chi2_surface <- function(K_grid, chisq_vals, label, hypothesis = "H3") {
  df <- data.frame(K = K_grid, chisq = chisq_vals)
  fname <- file.path(WIND_LOG$log_dir,
                     sprintf("chi2_surface_%s_%s.csv", hypothesis,
                             gsub("[^a-zA-Z0-9]", "_", label)))
  write.csv(df, fname, row.names = FALSE)

  # Summary
  k_best <- K_grid[which.min(chisq_vals)]
  bowl_depth <- max(chisq_vals) - min(chisq_vals)
  log_message(sprintf("Chi² surface [%s]: K_best=%.4f, bowl_depth=%.4f, range=[%.4f, %.4f]",
                      label, k_best, bowl_depth, min(chisq_vals), max(chisq_vals)))
}


# -----------------------------------------------------------------------------
#' get_master_log
#'
#' Read and return the master log as a data frame.
#'
#' @return Data frame of all logged fits
#' @export
get_master_log <- function() {
  if (file.exists(WIND_LOG$master_log)) {
    read.csv(WIND_LOG$master_log, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
}


# -----------------------------------------------------------------------------
#' print_summary
#'
#' Print a summary of the experiment to console and session log.
#'
#' @export
print_summary <- function() {
  elapsed <- as.numeric(difftime(Sys.time(), WIND_LOG$start_time, units = "mins"))
  ml <- get_master_log()

  log_section(sprintf("EXPERIMENT SUMMARY  (%.1f minutes elapsed)", elapsed))
  log_message(sprintf("Total fits logged: %d", nrow(ml)))

  if (nrow(ml) > 0) {
    log_message(sprintf("Convergence rate: %.1f%% (%d/%d)",
                        100 * mean(ml$converged, na.rm = TRUE),
                        sum(ml$converged, na.rm = TRUE), nrow(ml)))
    log_message(sprintf("Mean chi²: %.4f  (SD: %.4f)",
                        mean(ml$chisq_best, na.rm = TRUE),
                        sd(ml$chisq_best, na.rm = TRUE)))

    # Per-hypothesis summary
    for (h in unique(ml$hypothesis)) {
      sub <- ml[ml$hypothesis == h, ]
      log_message(sprintf("  %s: %d fits, mean chi²=%.4f, converged=%d/%d",
                          h, nrow(sub),
                          mean(sub$chisq_best, na.rm = TRUE),
                          sum(sub$converged, na.rm = TRUE), nrow(sub)))
    }
  }
  invisible(NULL)
}
