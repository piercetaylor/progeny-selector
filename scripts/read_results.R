#!/usr/bin/env Rscript
# Reads results.csv (docs/adr/0016, contract 1.4.0) with explicit column types.
#
# FIXED_COLUMNS mirrors progeny_selector.io.export.FIXED_COLUMNS exactly: the 34
# documented columns plus the trailing `token_profile` (35 names), which is the
# header written when there are no rows. Dynamic per-target, per-avoid-locus and
# per-chromosome columns, when present, fall between `results_schema` and
# `token_profile` and are read by readr's guess (`.default = col_guess()`).
suppressPackageStartupMessages(library(readr))

FIXED_COLUMNS <- c(
  "rank_overall",
  "rank_in_family",
  "sample_id",
  "line_name",
  "family_id",
  "generation",
  "passes_filters",
  "exclusion_reason",
  "composite_score",
  "foreground_all_pass",
  "avoid_all_pass",
  "rpp_total",
  "rpp_carrier",
  "rpp_noncarrier",
  "expected_rpp",
  "drag_total_est",
  "drag_total_max",
  "drag_unit",
  "ibs_rp",
  "ibs_donor",
  "missing_rate",
  "het_rate",
  "expected_het",
  "qc_flags",
  "role",
  "n_informative_called",
  "frac_a",
  "frac_h",
  "frac_b",
  "background_model",
  "background_unit",
  "rank_mode",
  "assembly",
  "results_schema",
  "token_profile"
)

# NA is missing (readr's default na = c("", "NA") would also treat an empty cell as
# missing, which is wrong here: empty text in exclusion_reason, qc_flags and notes is
# a genuine empty string, docs/adr/0016). na = "NA" alone keeps that distinction.
read_results <- function(path) {
  readr::read_csv(
    path,
    col_types = cols(
      rank_overall = col_integer(),
      rank_in_family = col_integer(),
      sample_id = col_character(),
      line_name = col_character(),
      family_id = col_character(),
      generation = col_character(),
      passes_filters = col_logical(),
      exclusion_reason = col_character(),
      composite_score = col_double(),
      foreground_all_pass = col_logical(),
      avoid_all_pass = col_logical(),
      rpp_total = col_double(),
      rpp_carrier = col_double(),
      rpp_noncarrier = col_double(),
      expected_rpp = col_double(),
      drag_total_est = col_double(),
      drag_total_max = col_double(),
      drag_unit = col_character(),
      ibs_rp = col_double(),
      ibs_donor = col_double(),
      missing_rate = col_double(),
      het_rate = col_double(),
      expected_het = col_double(),
      qc_flags = col_character(),
      role = col_character(),
      n_informative_called = col_integer(),
      frac_a = col_double(),
      frac_h = col_double(),
      frac_b = col_double(),
      background_model = col_character(),
      background_unit = col_character(),
      rank_mode = col_character(),
      assembly = col_character(),
      results_schema = col_character(),
      token_profile = col_character(),
      .default = col_guess()
    ),
    na = "NA",
    show_col_types = FALSE
  )
}

# Dynamic column families read by guess:
#   target_<id>_status, avoid_<id>_status   -> character (pass/fail/unknown)
#   drag_<id>_left_max, drag_<id>_right_max -> double
#   rpp_<chrom>                             -> double
#   recomb_<id>_left, recomb_<id>_right     -> logical

main <- function(args) {
  if (length(args) < 1) {
    stop("usage: read_results.R <results.csv>")
  }
  path <- args[[1]]
  df <- read_results(path)

  prefix_len <- length(FIXED_COLUMNS)
  actual_prefix <- names(df)[seq_len(min(prefix_len, ncol(df)))]
  if (!identical(actual_prefix, FIXED_COLUMNS[seq_along(actual_prefix)]) ||
      length(actual_prefix) < prefix_len) {
    expected <- FIXED_COLUMNS
    mismatch <- NA_integer_
    for (i in seq_along(expected)) {
      got <- if (i <= length(actual_prefix)) actual_prefix[[i]] else NA_character_
      if (!identical(got, expected[[i]])) {
        mismatch <- i
        break
      }
    }
    stop(sprintf(
      "results.csv fixed prefix mismatch at column %d: expected %s, got %s",
      mismatch,
      expected[[mismatch]],
      if (mismatch <= length(actual_prefix)) actual_prefix[[mismatch]] else "<missing>"
    ))
  }

  if (nrow(df) > 0) {
    schema <- df[["results_schema"]]
    if (!all(startsWith(schema, "1."))) {
      stop("results.csv results_schema is not a 1.x version in every row")
    }
  }

  locus_status_cols <- grep("^(target|avoid)_.*_status$", names(df), value = TRUE)

  cat(sprintf("rows: %d\n", nrow(df)))
  cat(sprintf("passing: %d\n", sum(df[["passes_filters"]], na.rm = TRUE)))
  cat("per-locus status columns:\n")
  for (col in locus_status_cols) {
    cat(sprintf("  %s\n", col))
  }
}

if (!interactive()) {
  main(commandArgs(trailingOnly = TRUE))
}
