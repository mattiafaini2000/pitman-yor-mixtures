# Each mode retains its own summary definition and burn-in.
source("R/io.R")
source("R/gibbs.R")
source("R/posterior.R")
source("config/galaxies.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 2L) stop("Usage: parameter_sensitivity.R [study] [r|rcpp]")
study <- if (length(arguments) >= 1L) arguments[1L] else "discount_density"
backend <- if (length(arguments) >= 2L) arguments[2L] else "rcpp"
if (!study %in% names(sensitivity_profiles) || !backend %in% c("r", "rcpp")) stop("Unknown study or backend.")
settings <- sensitivity_profiles[[study]]
parameters <- sampler_profiles$cpp_defaults
if (!is.null(galaxies_config$seed)) set.seed(galaxies_config$seed)
cpp_sampler <- if (backend == "rcpp") load_cpp_sampler() else NULL
summaries <- vector("list", length(settings$values))
for (parameter_index in seq_along(settings$values)) {
  parameters[[settings$parameter]] <- settings$values[parameter_index]
  fit <- do.call(fit_py_mixture, c(list(observations = MASS::galaxies,
                                       niter = galaxies_config$niter, backend = backend,
                                       cpp_sampler = cpp_sampler), parameters))
  if (settings$summary == "density") {
    density <- historical_mixture_density(fit$theta_store, fit$sigma2_store,
                                           galaxies_config$grid, parameters$concentration,
                                           parameters$discount, settings$burn_in,
                                           settings$denominator)
    summaries[[parameter_index]] <- data.frame(x = density$grid, density = density$density,
                                               parameter = settings$parameter,
                                               value = settings$values[parameter_index])
  } else {
    partition <- representative_partition(fit$theta_store, settings$burn_in)
    summaries[[parameter_index]] <- data.frame(observation = seq_along(fit$data_scaled),
                                               x = fit$data_scaled, cluster = partition$cluster,
                                               parameter = settings$parameter,
                                               value = settings$values[parameter_index],
                                               selected_iteration = partition$iteration)
  }
}
summary_table <- do.call(rbind, summaries)
output_directory <- ensure_project_directory(file.path("outputs/sensitivity", study))
write_project_csv(summary_table, file.path(output_directory, "summary.csv"))
saveRDS(list(summary = summary_table, settings = settings, sampler_parameters = sampler_profiles$cpp_defaults,
             niter = galaxies_config$niter, backend = backend,
             center = fit$center, scale = fit$scale),
         project_output_path(file.path(output_directory, "summary.rds")))
