# The report's finite-sample expected K_i, for i = 1, ..., n.
expected_cluster_count <- function(n, discount, concentration) {
  validate_prior_parameters(n, discount, concentration)
  observations <- seq_len(n)
  if (discount == 0) {
    return(cumsum(concentration / (concentration + observations - 1)))
  }
  concentration / discount * (
    cumprod((concentration + discount + observations - 1) /
              (concentration + observations - 1)) - 1
  )
}

# Source-panel asymptotic references. The DP panel explicitly puts K_1 = 1.
asymptotic_cluster_count <- function(n, discount, concentration,
                                     dp_first_point = TRUE) {
  validate_prior_parameters(n, discount, concentration)
  observations <- seq_len(n)
  if (discount == 0) {
    values <- concentration * log(observations)
    if (dp_first_point) values[1] <- 1
    return(values)
  }
  gamma(concentration + 1) / (discount * gamma(concentration + discount)) *
    observations^discount
}

# Expected masses in stick-breaking (appearance) order, not size-ranked masses.
expected_stick_breaking_mass <- function(n_blocks, discount, concentration) {
  validate_prior_parameters(n_blocks, discount, concentration)
  indices <- seq_len(n_blocks)
  if (discount == 0) {
    return(1 / (concentration + 1) *
             (concentration / (concentration + 1))^(indices - 1))
  }
  residual_products <- 1
  if (n_blocks > 1L) {
    previous <- seq_len(n_blocks - 1L)
    residual_products <- c(
      1, cumprod((concentration + previous * discount) /
                   (1 + concentration + (previous - 1) * discount))
    )
  }
  (1 - discount) / (1 + concentration + (indices - 1) * discount) *
    residual_products
}

# Formula for E[K_{n,r}] from the report. Rising factorials use the ordinary
# unit increment; n and each block size must be positive integers, size <= n.
expected_blocks_by_size <- function(n, sizes = seq_len(n), discount,
                                    concentration) {
  validate_prior_parameters(n, discount, concentration)
  if (any(!is.finite(sizes)) || any(sizes < 1 | sizes > n | sizes != as.integer(sizes))) {
    stop("sizes must be integer block sizes between 1 and n.")
  }
  log_rising_factorial <- function(start, terms) {
    if (terms == 0L) return(0)
    lgamma(start + terms) - lgamma(start)
  }
  vapply(sizes, function(size) {
    exp(lchoose(n, size) + log_rising_factorial(1 - discount, size - 1) +
          log_rising_factorial(concentration + discount, n - size) -
          log_rising_factorial(concentration + 1, n - 1))
  }, numeric(1))
}

# A reference with fixed slope -1/d and an amplitude fitted to the same tail.
# The optional d = 0 curve separately fits amplitude * exp(-rank).
fit_rank_reference <- function(rank_summary, discount, min_tail_rank = 5L,
                               include_dp_exponential = FALSE) {
  tail <- rank_summary[rank_summary$rank >= min_tail_rank &
                         rank_summary$mean_proportion > 0, , drop = FALSE]
  fitted <- rep(NA_real_, nrow(rank_summary))
  log_amplitude <- NA_real_
  reference <- "none"
  if (discount > 0 && nrow(tail) > 0L) {
    log_amplitude <- mean(log(tail$mean_proportion) + (1 / discount) * log(tail$rank))
    fitted <- exp(log_amplitude) * rank_summary$rank^(-1 / discount)
    reference <- "fixed_slope_power_law"
  } else if (discount == 0 && include_dp_exponential && nrow(tail) > 0L) {
    log_amplitude <- mean(log(tail$mean_proportion) + tail$rank)
    fitted <- exp(log_amplitude) * exp(-rank_summary$rank)
    reference <- "fixed_rate_exponential"
  }
  data.frame(rank_summary, fitted_reference = fitted,
             log_amplitude = log_amplitude, reference = reference,
             min_tail_rank = min_tail_rank)
}

# The cluster-size illustration fixes slope -(1+d) after truncation, and fits
# its log amplitude using observed positive frequencies at size >= min_tail_size.
fit_cluster_size_reference <- function(size_summary, discount, min_tail_size = 3L) {
  tail <- size_summary[size_summary$size >= min_tail_size &
                         size_summary$probability > 0, , drop = FALSE]
  log_amplitude <- if (nrow(tail) > 0L) {
    mean(log(tail$probability) + (1 + discount) * log(tail$size))
  } else NA_real_
  data.frame(size_summary,
             fitted_reference = exp(log_amplitude) * size_summary$size^(-(1 + discount)),
             log_amplitude = log_amplitude, min_tail_size = min_tail_size)
}
