# Plot summary tables only; these helpers do not simulate partitions.
plot_prior_cluster_counts <- function(count_summary) {
  ggplot2::ggplot(count_summary, ggplot2::aes(x = observations)) +
    ggplot2::geom_line(ggplot2::aes(y = empirical, colour = "Empirical")) +
    ggplot2::geom_line(ggplot2::aes(y = exact, colour = "Exact expectation")) +
    ggplot2::geom_line(ggplot2::aes(y = asymptotic, colour = "Asymptotic reference"),
                       linetype = "dashed") +
    ggplot2::scale_colour_manual(values = c("Empirical" = "blue",
                                            "Exact expectation" = "red",
                                            "Asymptotic reference" = "green")) +
    ggplot2::labs(x = "Number of observations", y = "Expected occupied blocks", colour = NULL) +
    ggplot2::theme_minimal()
}

plot_prior_mass <- function(mass_summary) {
  ggplot2::ggplot(mass_summary, ggplot2::aes(x = rank)) +
    ggplot2::geom_line(ggplot2::aes(y = mean_proportion), colour = "blue") +
    ggplot2::geom_point(ggplot2::aes(y = mean_proportion), colour = "blue") +
    ggplot2::geom_line(ggplot2::aes(y = stick_breaking_reference),
                       colour = "red", linetype = "dashed") +
    ggplot2::facet_wrap(~ordering) +
    ggplot2::labs(x = "Block index", y = "Mean finite-sample block proportion",
                  caption = "Red: stick-breaking mass reference in appearance order") +
    ggplot2::theme_minimal()
}

plot_prior_rank_reference <- function(rank_summary, include_dp_exponential = FALSE) {
  rank_summary$discount_label <- factor(round(rank_summary$discount, 2))
  plot <- ggplot2::ggplot(rank_summary, ggplot2::aes(x = rank, colour = discount_label)) +
    ggplot2::geom_line(ggplot2::aes(y = mean_proportion)) +
    ggplot2::geom_point(ggplot2::aes(y = mean_proportion), alpha = 0.5, size = 1) +
    ggplot2::geom_line(
      data = rank_summary[rank_summary$discount > 0, , drop = FALSE],
      ggplot2::aes(y = fitted_reference), linetype = "dashed", na.rm = TRUE
    ) +
    ggplot2::scale_x_log10(limits = c(4, NA)) +
    ggplot2::labs(x = "Block rank (1 = largest)", y = "Mean finite-sample block proportion",
                  colour = "Discount") + ggplot2::theme_minimal()
  if (include_dp_exponential) {
    plot <- plot + ggplot2::geom_line(
      data = rank_summary[rank_summary$discount == 0, , drop = FALSE],
      ggplot2::aes(y = fitted_reference), linetype = "dotted", linewidth = 0.7,
      na.rm = TRUE
    ) + ggplot2::scale_y_log10(limits = c(1e-8, 1))
  } else {
    plot <- plot + ggplot2::scale_y_log10(limits = c(NA, 1))
  }
  plot
}

plot_prior_cluster_size_reference <- function(size_summary) {
  size_summary$discount_label <- factor(round(size_summary$discount, 2))
  ggplot2::ggplot(size_summary, ggplot2::aes(x = size, colour = discount_label)) +
    ggplot2::geom_line(ggplot2::aes(y = probability)) +
    ggplot2::geom_point(ggplot2::aes(y = probability), alpha = 0.5, size = 1) +
    ggplot2::geom_line(ggplot2::aes(y = fitted_reference), linetype = "dashed") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Occupied block size", y = "Pooled block-size frequency",
                  colour = "Discount") + ggplot2::theme_minimal()
}
