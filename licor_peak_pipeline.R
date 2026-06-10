# LICOR CO2 injection peak extraction pipeline
# Purpose:
#   1) Read a LICOR .txt output file
#   2) Detect sequential CO2 injection peaks
#   3) Extract peak maximum CO2 and instantaneous H2O at the peak maximum
#   4) Calculate baseline-corrected CO2 area under each peak
#   5) Calculate average H2O across each detected peak width
#   6) Export a polished CSV summary

# ---------------------------
# User settings
# ---------------------------
#Where your LiCOR output text file is located; make sure it corresponds to your file path
input_file  <- "Day1_T0_LICOR_GW_SP.txt"

#Where you want this data to be saved; make sure it corresponds to your file system.
output_file <- "LICOR_CO2_peak_summary.csv"

# These are cleaned versions of the LICOR headers.
co2_col <- "CO2_(umol_mol-1)"
h2o_col <- "H2O_(mmol_mol-1)"

# Peak detection parameters. Tune these if needed.
smooth_window <- 5        # odd integer; light smoothing for noisy CO2 trace
threshold_mad <- 4        # higher = fewer peaks; lower = more sensitive
min_prominence <- 8       # minimum CO2 rise above local baseline, in umol mol-1
min_peak_gap_sec <- 8     # minimum time between separate injections, in seconds
edge_fraction <- 0.08     # peak edges where signal falls to baseline + 8% of peak height

# TRUE = integrate CO2 above the local baseline, usually best for injection peaks.
# FALSE = integrate raw CO2 values.
baseline_correct_area <- TRUE

# ---------------------------
# Helper functions
# ---------------------------
clean_names_licor <- function(x) {
  x <- gsub("CO2", "CO2", x, fixed = TRUE)
  x <- gsub("CO2", "CO2", x, fixed = TRUE)
  x <- gsub("CO₂", "CO2", x, fixed = TRUE)
  x <- gsub("H₂O", "H2O", x, fixed = TRUE)
  x <- gsub("µ", "u", x, fixed = TRUE)
  x <- gsub("⁻¹", "-1", x, fixed = TRUE)
  x <- gsub("°", "deg", x, fixed = TRUE)
  x <- gsub(":", "", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9_()./-]", "", x)
  x
}

roll_mean_centered <- function(x, k = 5) {
  if (k <= 1) return(x)
  if (k %% 2 == 0) stop("smooth_window must be an odd integer.")
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

trapz <- function(x, y) {
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

# ---------------------------
# Read and prepare data
# ---------------------------
raw <- read.delim(
  input_file,
  skip = 1,                 # first line is LICOR metadata, not the header
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Drop empty trailing columns if present.
raw <- raw[, colSums(!is.na(raw) & raw != "") > 0, drop = FALSE]
names(raw) <- clean_names_licor(names(raw))

if (!all(c(co2_col, h2o_col) %in% names(raw))) {
  stop(
    "Could not find required columns. Available columns are:\n",
    paste(names(raw), collapse = "\n")
  )
}

raw$datetime <- as.POSIXct(
  paste(raw$`System_Date_(Y-M-D)`, raw$`System_Time_(hms)`),
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

# LICOR output may round times to whole seconds, causing duplicate timestamps.
# For integration, use evenly spaced elapsed time if timestamps are duplicated.
if (any(duplicated(raw$datetime))) {
  total_elapsed <- as.numeric(max(raw$datetime, na.rm = TRUE) - min(raw$datetime, na.rm = TRUE), units = "secs")
  raw$elapsed_sec <- seq(0, total_elapsed, length.out = nrow(raw))
} else {
  raw$elapsed_sec <- as.numeric(raw$datetime - raw$datetime[1], units = "secs")
}

raw$co2 <- as.numeric(raw[[co2_col]])
raw$h2o <- as.numeric(raw[[h2o_col]])

# Smooth CO2 only for detection. Use raw CO2 for reported maxima and area.
raw$co2_smooth <- roll_mean_centered(raw$co2, smooth_window)
raw$co2_smooth[is.na(raw$co2_smooth)] <- raw$co2[is.na(raw$co2_smooth)]

# Dynamic baseline / threshold.
baseline_global <- stats::median(raw$co2_smooth, na.rm = TRUE)
noise_mad <- stats::mad(raw$co2_smooth, constant = 1, na.rm = TRUE)
threshold <- baseline_global + threshold_mad * noise_mad

# ---------------------------
# Find local maxima above threshold
# ---------------------------
co2s <- raw$co2_smooth
previous_value <- c(-Inf, head(co2s, -1))
next_value <- c(tail(co2s, -1), -Inf)

candidate_idx <- which(co2s > threshold & co2s >= previous_value & co2s > next_value)

# Keep the tallest candidate within each minimum gap window.
if (length(candidate_idx) > 0) {
  kept <- integer(0)
  for (idx in candidate_idx) {
    if (length(kept) == 0) {
      kept <- c(kept, idx)
    } else {
      time_gap <- raw$elapsed_sec[idx] - raw$elapsed_sec[tail(kept, 1)]
      if (time_gap < min_peak_gap_sec) {
        if (raw$co2[idx] > raw$co2[tail(kept, 1)]) kept[length(kept)] <- idx
      } else {
        kept <- c(kept, idx)
      }
    }
  }
  peak_idx <- kept
} else {
  peak_idx <- integer(0)
}

# ---------------------------
# Determine peak widths and summarize each peak
# ---------------------------
results <- data.frame()

for (i in seq_along(peak_idx)) {
  pk <- peak_idx[i]

  left_search_start <- if (i == 1) 1 else peak_idx[i - 1]
  right_search_end <- if (i == length(peak_idx)) nrow(raw) else peak_idx[i + 1]
  local_region <- raw$co2_smooth[left_search_start:right_search_end]
  local_baseline <- stats::quantile(local_region, probs = 0.15, na.rm = TRUE)

  peak_height <- raw$co2_smooth[pk] - local_baseline
  if (is.na(peak_height) || peak_height < min_prominence) next

  edge_level <- local_baseline + edge_fraction * peak_height

  left <- pk
  while (left > 1 && raw$co2_smooth[left] > edge_level) left <- left - 1

  right <- pk
  while (right < nrow(raw) && raw$co2_smooth[right] > edge_level) right <- right + 1

  idx <- left:right

  y_for_area <- raw$co2[idx]
  if (baseline_correct_area) y_for_area <- pmax(raw$co2[idx] - local_baseline, 0)

  max_row <- idx[which.max(raw$co2[idx])]

  results <- rbind(
    results,
    data.frame(
      `sample number` = nrow(results) + 1,
      `maximum co2 (umol_mol-1)` = raw$co2[max_row],
      `Instantaneous H2O` = raw$h2o[max_row],
      `CO2 curve area` = trapz(raw$elapsed_sec[idx], y_for_area),
      `Average H2O` = mean(raw$h2o[idx], na.rm = TRUE),
      check.names = FALSE
    )
  )
}

write.csv(results, output_file, row.names = FALSE)

message("Detected ", nrow(results), " peaks.")
message("Saved summary to: ", normalizePath(output_file, mustWork = FALSE))
