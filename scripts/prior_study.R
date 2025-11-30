# Run from the repository root: Rscript --vanilla scripts/prior_study.R rank_120
source("R/io.R")
source("R/prior.R")
source("R/theory.R")
source("R/prior_plotting.R")
source("config/prior_study.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 2L) stop("Usage: prior_study.R [variant] [output_directory]")
variant <- if (length(arguments) >= 1L) arguments[1] else "report_dp"
settings <- prior_study_config(variant)
output_directory <- if (length(arguments) >= 2L) arguments[2] else {
  file.path("outputs", "prior", variant)
}
output_directory <- ensure_project_directory(output_directory)
if (!is.null(settings$seed)) set.seed(settings$seed)

if (settings$study == "observations") {
  study <- simulate_normal_observations(
    settings$sample_size, settings$discounts, settings$concentration,
    settings$base_mean, settings$base_variance, settings$observation_sd
  )
  saveRDS(list(configuration = settings, study = study),
          project_output_path(file.path(output_directory, "study.rds")))
} else {
  all_partitions <- lapply(settings$discounts, function(discount) {
    simulate_prior_replicates(
      settings$sample_size, discount, settings$concentration,
      settings$replications,
      function() stats::rnorm(1L, settings$base_mean, settings$base_sd),
      settings$sampling
    )
  })
  # The standalone rank studies use one rank truncation across the full grid,
  # then zero-pad every discount/replication to that common rank bound.
  rank_bound <- settings$max_rank
  if (settings$study == "rank_size") {
    available_rank <- max(vapply(all_partitions, function(partitions) {
      max(lengths(lapply(partitions, `[[`, "cluster_sizes")))
    }, numeric(1)))
    rank_bound <- min(settings$max_rank, available_rank)
  }
  studies <- lapply(seq_along(settings$discounts), function(index) {
    discount <- settings$discounts[index]
    partitions <- all_partitions[[index]]
    if (settings$study == "count_mass") {
      count_summary <- data.frame(
        observations = seq_len(settings$sample_size),
        empirical = rowMeans(do.call(cbind, lapply(partitions, `[[`, "cluster_count"))),
        exact = expected_cluster_count(settings$sample_size, discount, settings$concentration),
        asymptotic = asymptotic_cluster_count(settings$sample_size, discount, settings$concentration)
      )
      mass_summary <- do.call(rbind, lapply(c("appearance", "ranked"), function(ordering) {
        summary <- summarize_rank_proportions(
          partitions, settings$sample_size, ordering, settings$mean_mode,
          settings$max_rank, settings$min_occurrences
        )
        summary$stick_breaking_reference <- expected_stick_breaking_mass(
          max(lengths(lapply(partitions, `[[`, "cluster_sizes"))),
          discount, settings$concentration
        )[summary$rank]
        summary
      }))
      list(discount = discount, counts = count_summary, mass = mass_summary)
    } else if (settings$study == "rank_size") {
      summary <- summarize_rank_proportions(
        partitions, settings$sample_size, "ranked", settings$mean_mode,
        rank_bound, settings$min_occurrences
      )
      summary <- fit_rank_reference(summary, discount, settings$min_tail_rank,
                                    settings$include_dp_exponential)
      summary$discount <- discount
      summary
    } else {
      summary <- pooled_cluster_size_distribution(partitions, settings$sample_size)
      summary <- summary[summary$size <= settings$max_size, , drop = FALSE]
      summary <- fit_cluster_size_reference(summary, discount, settings$min_tail_size)
      summary$discount <- discount
      summary
    }
  })
  saveRDS(list(configuration = settings, summaries = studies),
          project_output_path(file.path(output_directory, "study.rds")))
  if (settings$study == "count_mass") {
    write_project_csv(studies[[1]]$counts, file.path(output_directory, "cluster_counts.csv"))
    write_project_csv(studies[[1]]$mass, file.path(output_directory, "block_proportions.csv"))
    figures <- list(cluster_counts = plot_prior_cluster_counts(studies[[1]]$counts),
                    block_proportions = plot_prior_mass(studies[[1]]$mass))
  } else {
    summaries <- do.call(rbind, studies)
    write_project_csv(summaries, file.path(output_directory, "summary.csv"))
    figures <- if (settings$study == "rank_size") {
      list(rank_size = plot_prior_rank_reference(summaries, settings$include_dp_exponential))
    } else list(cluster_sizes = plot_prior_cluster_size_reference(summaries))
  }
  for (figure_name in names(figures)) {
    ggplot2::ggsave(project_output_path(file.path(output_directory, paste0(figure_name, ".pdf"))),
                    plot = figures[[figure_name]], width = 7, height = 4.5)
  }
}
