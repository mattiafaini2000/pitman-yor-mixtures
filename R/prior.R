# Sequential prior partitions. Concentration and discount correspond to c and d.
validate_prior_parameters <- function(n, discount, concentration) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (length(discount) != 1L || !is.finite(discount) ||
      discount < 0 || discount >= 1) {
    stop("discount must be in [0, 1).")
  }
  if (length(concentration) != 1L || !is.finite(concentration) ||
      concentration <= 0) {
    stop("This implementation requires positive concentration.")
  }
  invisible(NULL)
}

# Return latent draws, appearance-order allocations, sizes, and the K_i path.
# base_draw is called once for each newly occupied block, including the first.
simulate_partition <- function(n, discount, concentration, base_draw,
                               sampling = c("sample", "sample_int")) {
  validate_prior_parameters(n, discount, concentration)
  sampling <- match.arg(sampling)
  if (!is.function(base_draw)) stop("base_draw must be a function.")

  latent_values <- numeric(n)
  allocations <- numeric(n)
  cluster_count <- numeric(n)
  latent_values[1] <- base_draw()
  allocations[1] <- 1
  cluster_count[1] <- 1
  occupied_values <- latent_values[1]
  cluster_sizes <- 1

  if (n > 1L) {
    for (observation in seq.int(2L, n)) {
      # The first category creates a new block; subsequent categories retain
      # appearance order. Do not move the base draw ahead of this choice.
      probabilities <- c(
        concentration + discount * cluster_count[observation - 1L],
        cluster_sizes - discount
      ) / (concentration + observation - 1L)
      if (sampling == "sample") {
        selected <- sample(seq_along(probabilities), 1L, prob = probabilities)
      } else {
        selected <- sample.int(length(probabilities), 1L, prob = probabilities)
      }

      if (selected == 1L) {
        latent_values[observation] <- base_draw()
        occupied_values <- c(occupied_values, latent_values[observation])
        cluster_sizes <- c(cluster_sizes, 1)
        cluster_count[observation] <- cluster_count[observation - 1L] + 1
        allocations[observation] <- cluster_count[observation]
      } else {
        block <- selected - 1L
        latent_values[observation] <- occupied_values[block]
        allocations[observation] <- block
        cluster_sizes[block] <- cluster_sizes[block] + 1
        cluster_count[observation] <- cluster_count[observation - 1L]
      }
    }
  }

  list(latent_values = latent_values, allocations = allocations,
       cluster_sizes = cluster_sizes, cluster_count = cluster_count)
}

# Replications remain sequential; RNG state is controlled by the entry point.
simulate_prior_replicates <- function(n, discount, concentration, replications,
                                      base_draw, sampling = "sample") {
  if (length(replications) != 1L || !is.finite(replications) ||
      replications < 1 || replications != as.integer(replications)) {
    stop("replications must be a positive integer.")
  }
  lapply(seq_len(replications), function(replication) {
    simulate_partition(n, discount, concentration, base_draw, sampling)
  })
}

# Finite-sample proportions, not draws of the infinite random probability mass.
# Conditional means use only replications where the rank exists; zero_filled
# means give an absent rank proportion zero before averaging all replications.
summarize_rank_proportions <- function(partitions, n,
                                       ordering = c("ranked", "appearance"),
                                       mode = c("conditional", "zero_filled"),
                                       max_rank = NULL, min_occurrences = 0L) {
  ordering <- match.arg(ordering)
  mode <- match.arg(mode)
  proportions <- lapply(partitions, function(partition) {
    sizes <- partition$cluster_sizes
    if (ordering == "ranked") sizes <- sort(sizes, decreasing = TRUE)
    sizes / n
  })
  available_rank <- max(lengths(proportions))
  if (is.null(max_rank)) max_rank <- available_rank
  if (max_rank < 1 || max_rank != as.integer(max_rank)) {
    stop("max_rank must be a positive integer.")
  }
  rows <- lapply(seq_len(max_rank), function(rank) {
    exists <- lengths(proportions) >= rank
    present <- vapply(proportions[exists], function(values) values[rank], numeric(1))
    mean_proportion <- if (mode == "conditional") mean(present) else {
      complete <- numeric(length(partitions))
      complete[exists] <- present
      mean(complete)
    }
    data.frame(rank = rank, mean_proportion = mean_proportion,
               occurrences = sum(exists), replications = length(partitions),
               ordering = ordering, mean_mode = mode)
  })
  summary <- do.call(rbind, rows)
  summary[summary$occurrences >= min_occurrences, , drop = FALSE]
}

# Count every occupied block across every replication before normalizing.
# This is a pooled distribution over blocks, not a mean of per-partition ratios.
pooled_cluster_size_distribution <- function(partitions, n) {
  counts <- integer(n)
  for (partition in partitions) {
    size_counts <- table(partition$cluster_sizes)
    observed_sizes <- as.integer(names(size_counts))
    counts[observed_sizes] <- counts[observed_sizes] + as.integer(size_counts)
  }
  observed_sizes <- which(counts > 0)
  data.frame(size = observed_sizes, pooled_count = counts[observed_sizes],
             probability = counts[observed_sizes] / sum(counts),
             total_blocks = sum(counts))
}

# The supplied synthetic example generated observations with sd = 0.1;
# its separate sigma2 = 1 assignment did not control these observations.
simulate_normal_observations <- function(n, discount, concentration,
                                        base_mean = 0, base_variance = 0.5,
                                        observation_sd = 0.1) {
  partition <- simulate_partition(
    n, discount, concentration,
    function() stats::rnorm(1L, mean = base_mean, sd = sqrt(base_variance))
  )
  observations <- stats::rnorm(n, mean = partition$latent_values,
                              sd = observation_sd)
  list(observations = observations, partition = partition,
       base_mean = base_mean, base_variance = base_variance,
       observation_sd = observation_sd)
}
