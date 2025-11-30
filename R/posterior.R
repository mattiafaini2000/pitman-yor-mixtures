# Burn-in is a count of initial stored iterations, not a sampler setting.
posterior_iterations <- function(theta_store, burn_in = 100L) {
  if (!is.matrix(theta_store) || length(burn_in) != 1L || burn_in < 0L ||
      burn_in != as.integer(burn_in) || burn_in >= nrow(theta_store)) {
    stop("burn_in must leave at least one stored matrix row.")
  }
  seq.int(burn_in + 1L, nrow(theta_store))
}

# Historical Gaussian-mixture summaries use undiscounted counts and omit the
# base-measure term. The denominator is explicit; coordinates stay standardized.
historical_mixture_density <- function(theta_store, sigma2_store, grid,
                                       concentration = 1, discount = 0.3,
                                       burn_in = 100L,
                                       denominator = c("n_plus_concentration", "n_plus_discount")) {
  denominator <- match.arg(denominator)
  kept_iterations <- posterior_iterations(theta_store, burn_in)
  if (length(sigma2_store) != nrow(theta_store)) stop("Variance and theta rows must align.")
  observation_count <- ncol(theta_store)
  normalization <- observation_count + if (denominator == "n_plus_concentration") {
    concentration
  } else discount
  density_store <- matrix(0, nrow = length(kept_iterations), ncol = length(grid))
  for (saved_index in seq_along(kept_iterations)) {
    iteration <- kept_iterations[saved_index]
    cluster_counts <- plyr::count(theta_store[iteration, ])
    for (grid_index in seq_along(grid)) {
      density_store[saved_index, grid_index] <- sum(
        cluster_counts$freq / normalization *
          stats::dnorm(grid[grid_index], cluster_counts$x, sqrt(sigma2_store[iteration]))
      )
    }
  }
  list(grid = grid, density = colMeans(density_store), density_store = density_store,
       denominator = denominator, burn_in = burn_in)
}

# Exact theta equality encodes common membership; both source calculations are
# retained. Dahl's representative is chosen only among retained sampled rows.
representative_partition <- function(theta_store, burn_in = 100L,
                                      calculation = c("outer", "scalar")) {
  calculation <- match.arg(calculation)
  kept_iterations <- posterior_iterations(theta_store, burn_in)
  saved_labels <- theta_store[kept_iterations, , drop = FALSE]
  saved_count <- nrow(saved_labels)
  observation_count <- ncol(saved_labels)
  similarity <- matrix(0, observation_count, observation_count)
  if (calculation == "outer") {
    for (saved_index in seq_len(saved_count)) {
      labels <- saved_labels[saved_index, ]
      similarity <- similarity + outer(labels, labels, FUN = "==")
    }
  } else {
    for (saved_index in seq_len(saved_count)) {
      for (first_index in seq_len(observation_count)) {
        for (second_index in seq_len(observation_count)) {
          similarity[first_index, second_index] <- similarity[first_index, second_index] +
            ifelse(saved_labels[saved_index, first_index] ==
                     saved_labels[saved_index, second_index], 1, 0)
        }
      }
    }
  }
  similarity <- similarity / saved_count
  squared_distance <- numeric(saved_count)
  for (saved_index in seq_len(saved_count)) {
    if (calculation == "outer") {
      adjacency <- outer(saved_labels[saved_index, ], saved_labels[saved_index, ], FUN = "==")
      difference <- similarity - adjacency
      squared_distance[saved_index] <- sum(difference * difference)
    } else {
      for (first_index in seq_len(observation_count)) {
        for (second_index in seq_len(observation_count)) {
          squared_distance[saved_index] <- squared_distance[saved_index] +
            (similarity[first_index, second_index] -
               ifelse(saved_labels[saved_index, first_index] ==
                        saved_labels[saved_index, second_index], 1, 0))^2
        }
      }
    }
  }
  selected_iteration <- kept_iterations[which.min(squared_distance)]
  theta <- theta_store[selected_iteration, ]
  list(iteration = selected_iteration, theta = theta,
       cluster = match(theta, unique(theta)), similarity = similarity,
       squared_distance = squared_distance, retained_iterations = kept_iterations)
}

rank_sizes_partition <- function(theta) {
  sizes <- tabulate(match(theta, unique(theta)))
  sort(sizes, decreasing = TRUE)
}

# Earlier report exploration: use the last stored iteration, round its centers,
# then match each original center to the nearest rounded value (first tie).
rounded_last_partition <- function(theta_store, digits = 3L) {
  if (!is.matrix(theta_store) || nrow(theta_store) < 1L) stop("Supply stored theta rows.")
  last_theta <- theta_store[nrow(theta_store), ]
  rounded_means <- unique(round(last_theta, digits))
  assignment <- vapply(last_theta, function(value) {
    which.min(abs(rounded_means - value))
  }, integer(1))
  list(iteration = nrow(theta_store), centers = rounded_means, cluster = assignment)
}

# Missing ranks receive zero before averaging, including iterations with fewer
# occupied clusters. Return proportions of the observed sample, not random masses.
posterior_rank_size <- function(theta_store, burn_in = 100L) {
  kept_iterations <- posterior_iterations(theta_store, burn_in)
  rank_samples <- lapply(kept_iterations, function(iteration) {
    rank_sizes_partition(theta_store[iteration, ]) / ncol(theta_store)
  })
  maximum_rank <- max(lengths(rank_samples))
  rank_matrix <- matrix(0, length(rank_samples), maximum_rank)
  for (saved_index in seq_along(rank_samples)) {
    rank_matrix[saved_index, seq_along(rank_samples[[saved_index]])] <- rank_samples[[saved_index]]
  }
  data.frame(rank = seq_len(maximum_rank), proportion = colMeans(rank_matrix))
}

# Exploratory impurity (1 - sum(p^2)) and entropy on supplied partitions.
# The global quantities use joint contingency proportions, not an accuracy score.
partition_impurity <- function(reference_partition, estimated_partition) {
  if (length(reference_partition) != length(estimated_partition) ||
      length(reference_partition) == 0L || anyNA(reference_partition) ||
      anyNA(estimated_partition)) stop("Supply two complete aligned partition vectors.")
  by_reference <- lapply(unique(reference_partition), function(reference_label) {
    subset_estimated <- estimated_partition[reference_partition == reference_label]
    proportions <- table(subset_estimated) / length(subset_estimated)
    data.frame(reference = as.character(reference_label),
               impurity = 1 - sum(proportions^2),
               entropy = -sum(proportions * log(proportions + 1e-12)))
  })
  contingency <- table(reference_partition, estimated_partition)
  joint_proportions <- as.vector(contingency) / sum(contingency)
  list(by_reference = do.call(rbind, by_reference), contingency = contingency,
       joint_entropy = -sum(joint_proportions * log(joint_proportions + 1e-12)),
       joint_impurity = 1 - sum(joint_proportions^2))
}
