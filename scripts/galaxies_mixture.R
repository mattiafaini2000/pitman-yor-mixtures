# Run from the repository root: Rscript --vanilla scripts/galaxies_mixture.R r report
source("R/io.R")
source("R/gibbs.R")
source("R/posterior.R")
source("config/galaxies.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 2L) stop("Usage: galaxies_mixture.R [r|rcpp] [report|cpp_defaults|cpp_lab_variant]")
backend <- if (length(arguments) >= 1L) arguments[1L] else "r"
profile_name <- if (length(arguments) >= 2L) arguments[2L] else "report"
if (!backend %in% c("r", "rcpp") || !profile_name %in% names(sampler_profiles)) stop("Unknown backend or profile.")
settings <- galaxies_config
parameters <- sampler_profiles[[profile_name]]
if (!is.null(settings$seed)) set.seed(settings$seed)
fit <- do.call(fit_py_mixture, c(list(observations = MASS::galaxies,
                                     niter = settings$niter, backend = backend), parameters))
density <- historical_mixture_density(fit$theta_store, fit$sigma2_store, settings$grid,
                                       concentration = parameters$concentration,
                                       discount = parameters$discount, burn_in = settings$burn_in)
partition <- representative_partition(fit$theta_store, settings$burn_in)
output_directory <- ensure_project_directory(settings$output_directory)
saveRDS(list(fit = fit, burn_in = settings$burn_in, profile = profile_name,
             density = density, partition = partition),
         project_output_path(file.path(output_directory, "mixture.rds")))
write_project_csv(data.frame(x = density$grid, density = density$density),
                   file.path(output_directory, "density.csv"))
write_project_csv(data.frame(observation = seq_along(fit$data_scaled),
                             x = fit$data_scaled, cluster = partition$cluster),
                   file.path(output_directory, "partition.csv"))
