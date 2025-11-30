# The historical exploration uses a synthetic normal sample, not galaxies.
source("R/io.R")
source("R/gibbs.R")
source("R/posterior.R")
source("config/galaxies.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 1L) stop("Usage: posterior_rank_size.R [r|rcpp]")
backend <- if (length(arguments)) arguments[1L] else "rcpp"
if (!backend %in% c("r", "rcpp")) stop("Unknown backend.")
settings <- posterior_rank_config
set.seed(settings$seed)
observations <- stats::rnorm(settings$sample_size, 0, 1)
parameters <- sampler_profiles$cpp_defaults
parameters$concentration <- settings$concentration
cpp_sampler <- if (backend == "rcpp") load_cpp_sampler() else NULL
curves <- vector("list", length(settings$discounts))
for (discount_index in seq_along(settings$discounts)) {
  parameters$discount <- settings$discounts[discount_index]
  # Synthetic observations were passed directly to C++, without standardization.
  sampler <- if (backend == "rcpp") cpp_sampler else py_sampler_r
  fit <- do.call(sampler, c(list(standardized_observations = observations,
                                 niter = settings$niter), parameters))
  curves[[discount_index]] <- posterior_rank_size(fit$theta_store, settings$burn_in)
}
maximum_rank <- max(vapply(curves, nrow, integer(1)))
curve_table <- do.call(rbind, lapply(seq_along(curves), function(discount_index) {
  curve <- curves[[discount_index]]
  data.frame(rank = seq_len(maximum_rank),
             proportion = c(curve$proportion, rep(0, maximum_rank - nrow(curve))),
             discount = settings$discounts[discount_index])
}))
output_directory <- ensure_project_directory(settings$output_directory)
write_project_csv(curve_table, file.path(output_directory, "rank_size.csv"))
saveRDS(list(curves = curve_table, observations = observations, settings = settings,
             parameters = sampler_profiles$cpp_defaults, backend = backend),
         project_output_path(file.path(output_directory, "rank_size.rds")))
