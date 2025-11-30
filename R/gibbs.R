# Standardize observations once; retain the transformation for interpretation.
standardize_observations <- function(observations) {
  observations <- as.numeric(observations)
  if (length(observations) < 2L || any(!is.finite(observations))) {
    stop("At least two finite observations are required.")
  }
  center <- mean(observations)
  scale <- stats::sd(observations)
  if (scale <= 0) stop("Observations must have positive sample variance.")
  list(values = (observations - center) / scale, center = center, scale = scale)
}

validate_mixture_arguments <- function(observations, niter, concentration,
                                       discount, base_variance,
                                       variance_shape, variance_rate) {
  if (length(observations) < 2L || any(!is.finite(observations))) {
    stop("At least two finite observations are required.")
  }
  if (length(niter) != 1L || niter < 1L || niter != as.integer(niter)) {
    stop("niter must be a positive integer.")
  }
  if (!is.finite(concentration) || concentration <= 0 ||
      !is.finite(discount) || discount < 0 || discount >= 1) {
    stop("This implementation supports concentration > 0 and 0 <= discount < 1.")
  }
  if (any(!is.finite(c(base_variance, variance_shape, variance_rate))) ||
      any(c(base_variance, variance_shape, variance_rate) <= 0)) {
    stop("Base variance and inverse-gamma parameters must be positive.")
  }
}

# R reference sampler on already standardized observations. All iterations are
# returned; posterior helpers discard burn-in exactly once.
py_sampler_r <- function(standardized_observations, niter = 1000L,
                         concentration = 1, discount = 0.3, base_mean = 0,
                         base_variance = 1, variance_shape = 2,
                         variance_rate = 1) {
  observations <- as.numeric(standardized_observations)
  validate_mixture_arguments(observations, niter, concentration, discount,
                             base_variance, variance_shape, variance_rate)
  if (!requireNamespace("plyr", quietly = TRUE)) stop("Install plyr to use the R sampler.")
  observation_count <- length(observations)
  theta_store <- matrix(0, nrow = niter, ncol = observation_count)
  sigma2_store <- numeric(niter)
  k_store <- numeric(niter)
  last_allocation_probabilities <- vector("list", niter)
  theta <- stats::rnorm(observation_count, 0, 1)
  sigma2 <- 1

  for (iteration in seq_len(niter)) {
    for (observation_index in seq_len(observation_count)) {
      other_theta <- theta[-observation_index]
      observation <- observations[observation_index]
      cluster_counts <- plyr::count(other_theta)
      other_cluster_means <- cluster_counts$x
      other_cluster_sizes <- cluster_counts$freq
      other_cluster_count <- length(other_cluster_means)
      allocation_weights <- numeric(other_cluster_count + 1L)
      reference_point <- 0
      conditional_mean <- base_variance / (base_variance + sigma2) * observation +
        sigma2 / (base_variance + sigma2) * base_mean
      conditional_variance <- (1 / base_variance + 1 / sigma2)^(-1)
      # New-cluster-first ordering and the original exponentiation are retained.
      allocation_weights[1L] <- exp(
        log(concentration + other_cluster_count * discount) -
          log(concentration + observation_count - 1) +
          stats::dnorm(observation, reference_point, sqrt(sigma2), log = TRUE) +
          stats::dnorm(reference_point, base_mean, sqrt(base_variance), log = TRUE) -
          stats::dnorm(reference_point, conditional_mean, sqrt(conditional_variance), log = TRUE)
      )
      allocation_weights[2:(other_cluster_count + 1L)] <- exp(
        log(other_cluster_sizes - discount) - log(concentration + observation_count - 1) +
          stats::dnorm(observation, other_cluster_means, sqrt(sigma2), log = TRUE)
      )
      allocation_probabilities <- allocation_weights / sum(allocation_weights)
      selected_cluster <- sample.int(other_cluster_count + 1L, size = 1L,
                                     prob = allocation_probabilities)
      if (selected_cluster == 1L) {
        theta[observation_index] <- stats::rnorm(1, conditional_mean,
                                               sqrt(conditional_variance))
      } else {
        theta[observation_index] <- other_cluster_means[selected_cluster - 1L]
      }
    }
    # pi contains only the final observation's vector in this iteration.
    last_allocation_probabilities[[iteration]] <- allocation_probabilities
    theta_store[iteration, ] <- theta
    posterior_shape <- variance_shape + observation_count / 2
    posterior_rate <- variance_rate + 1 / 2 * sum((observations - theta)^2)
    sigma2 <- 1 / stats::rgamma(1, shape = posterior_shape, rate = posterior_rate)
    sigma2_store[iteration] <- sigma2

    # Occupied-mean acceleration follows storage and the shared-variance update.
    cluster_counts <- plyr::count(theta)
    cluster_means <- cluster_counts$x
    cluster_count <- length(cluster_means)
    for (cluster_index in seq_len(cluster_count)) {
      member_indices <- which(theta == cluster_means[cluster_index])
      cluster_sample_mean <- mean(observations[member_indices])
      cluster_size <- length(member_indices)
      mean_variance <- (cluster_size / sigma2 + 1 / base_variance)^(-1)
      mean_location <- mean_variance *
        (base_mean / base_variance + cluster_size * cluster_sample_mean / sigma2)
      theta[member_indices] <- stats::rnorm(1, mean_location, sqrt(mean_variance))
    }
    k_store[iteration] <- cluster_count
  }
  list(theta_store = theta_store, sigma2_store = sigma2_store, k_store = k_store,
       pi = last_allocation_probabilities, data_scaled = observations)
}

# Compile only on an explicit call after R starts with a project-local tempdir.
load_cpp_sampler <- function(root = project_root()) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install Rcpp to use the C++ backend.")
  checked_project_path(tempdir(), root)
  cache_directory <- ensure_project_directory(".cache/rcpp", root)
  temporary_directory <- ensure_project_directory(".cache/tmp", root)
  previous_environment <- Sys.getenv(c("TMP", "TEMP", "TMPDIR"), unset = NA_character_)
  on.exit({
    for (environment_name in names(previous_environment)) {
      previous_value <- previous_environment[[environment_name]]
      if (is.na(previous_value)) Sys.unsetenv(environment_name) else {
        do.call(Sys.setenv, setNames(list(previous_value), environment_name))
      }
    }
  }, add = TRUE)
  Sys.setenv(TMP = temporary_directory, TEMP = temporary_directory, TMPDIR = temporary_directory)
  sampler_environment <- new.env(parent = globalenv())
  Rcpp::sourceCpp(file = file.path(root, "src", "gibbs.cpp"),
                  env = sampler_environment, cacheDir = cache_directory,
                  rebuild = FALSE, showOutput = FALSE)
  sampler_environment$py_sampler_cpp
}

# Fit raw observations with one centering/scaling step. R and Rcpp retain their
# own RNG routines; no agreement of seeded trajectories is assumed.
fit_py_mixture <- function(observations, niter = 1000L, concentration = 1,
                           discount = 0.3, base_mean = 0, base_variance = 1,
                           variance_shape = 2, variance_rate = 1,
                           backend = c("r", "rcpp"), cpp_sampler = NULL,
                           root = project_root()) {
  backend <- match.arg(backend)
  standardized <- standardize_observations(observations)
  validate_mixture_arguments(standardized$values, niter, concentration, discount,
                             base_variance, variance_shape, variance_rate)
  sampler <- py_sampler_r
  if (backend == "rcpp") {
    sampler <- if (is.null(cpp_sampler)) load_cpp_sampler(root) else cpp_sampler
  }
  fit <- sampler(standardized_observations = standardized$values, niter = niter,
                 concentration = concentration, discount = discount,
                 base_mean = base_mean, base_variance = base_variance,
                 variance_shape = variance_shape, variance_rate = variance_rate)
  fit$center <- standardized$center
  fit$scale <- standardized$scale
  fit$backend <- backend
  fit$parameters <- list(concentration = concentration, discount = discount,
                          base_mean = base_mean, base_variance = base_variance,
                          variance_shape = variance_shape, variance_rate = variance_rate)
  fit
}
