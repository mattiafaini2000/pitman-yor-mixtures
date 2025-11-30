# Report panels and standalone exploratory settings. seed = NULL keeps the
# standalone files' unspecified RNG state; 123 is the report setup seed.
prior_study_config <- function(variant = "report_dp") {
  settings <- list(
    variant = variant, sample_size = 1000L, replications = 300L,
    concentration = 1, discounts = 0, base_mean = 0, base_sd = 1,
    seed = NULL, sampling = "sample_int", study = "rank_size",
    mean_mode = "zero_filled", max_rank = 50L, min_occurrences = 0L,
    min_tail_rank = 5L, include_dp_exponential = FALSE,
    max_size = 50L, min_tail_size = 3L
  )
  if (variant %in% c("report_dp", "report_py")) {
    settings$discounts <- if (variant == "report_dp") 0 else 0.3
    settings$seed <- 123L
    settings$sampling <- "sample"
    settings$study <- "count_mass"
    settings$mean_mode <- "conditional"
    settings$max_rank <- NULL
    settings$min_occurrences <- 150L
  } else if (variant %in% c("rank_50", "rank_120", "cluster_sizes")) {
    settings$discounts <- seq(0, 0.9, length.out = 10)
    if (variant == "rank_120") {
      settings$max_rank <- 120L
      settings$include_dp_exponential <- TRUE
    }
    if (variant == "cluster_sizes") settings$study <- "cluster_size"
  } else if (variant == "synthetic_observations") {
    settings$sample_size <- 500L
    settings$discounts <- 0.3
    settings$replications <- 1L
    settings$sampling <- "sample"
    settings$study <- "observations"
    settings$base_variance <- 0.5
    settings$base_sd <- sqrt(0.5)
    settings$observation_sd <- 0.1
    settings$unused_source_sigma2 <- 1
  } else {
    stop("Choose report_dp, report_py, rank_50, rank_120, cluster_sizes, or synthetic_observations.")
  }
  settings
}
