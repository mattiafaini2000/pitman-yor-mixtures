galaxies_config <- list(
  niter = 1000L, burn_in = 100L, seed = NULL,
  grid = seq(-4, 4, length.out = 1000L),
  output_directory = "outputs/galaxies"
)

sampler_profiles <- list(
  report = list(concentration = 1, discount = 0.3, base_mean = 0,
                 base_variance = 1, variance_shape = 2, variance_rate = 1),
  cpp_defaults = list(concentration = 1, discount = 0.3, base_mean = 0,
                       base_variance = 1, variance_shape = 1, variance_rate = 1),
  cpp_lab_variant = list(concentration = 1, discount = 0.9, base_mean = 0,
                          base_variance = 1, variance_shape = 1, variance_rate = 1)
)

sensitivity_profiles <- list(
  discount_density = list(parameter = "discount", values = seq(0, 0.9, length.out = 10L),
                           summary = "density", burn_in = 100L,
                           denominator = "n_plus_discount"),
  discount_clusters = list(parameter = "discount", values = seq(0, 0.9, length.out = 10L),
                            summary = "partition", burn_in = 100L),
  concentration_density = list(parameter = "concentration", values = seq(0.5, 1.5, length.out = 10L),
                                summary = "density", burn_in = 100L,
                                denominator = "n_plus_concentration"),
  concentration_clusters = list(parameter = "concentration", values = seq(0.5, 1.5, length.out = 10L),
                                 summary = "partition", burn_in = 200L)
)

posterior_rank_config <- list(
  dataset = "normal_synthetic", sample_size = 1000L, seed = 123L,
  niter = 1000L, burn_in = 100L, concentration = 1,
  discounts = seq(0, 0.9, length.out = 10L),
  output_directory = "outputs/posterior_rank_size"
)
