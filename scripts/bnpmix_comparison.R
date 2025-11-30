arguments <- commandArgs(trailingOnly = TRUE)
profile <- if (length(arguments)) arguments[[1]] else "location"
source("R/io.R")
source("R/bnpmix.R")
root <- project_root()
observations <- MASS::galaxies
center <- mean(observations)
scale_value <- stats::sd(observations)
standardized <- (observations - center) / scale_value
fit <- fit_bnpmix(standardized, profile)
result <- list(fit = fit, center = center, scale = scale_value,
               profile = profile, settings = bnpmix_profiles(standardized)[[profile]],
               vi_partition = bnpmix_partition(fit, "VI"))
saveRDS(result, project_output_path(file.path("outputs/bnpmix", paste0(profile, ".rds")), root))
