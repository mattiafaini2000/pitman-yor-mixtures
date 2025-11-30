# Optional package examples, separate from the custom shared-variance sampler.
bnpmix_profiles <- function(standardized_observations) {
  observations <- standardized_observations
  margin <- 0.1 * stats::sd(observations)
  local_grid <- seq(min(observations) - margin, max(observations) + margin,
                    length.out = 100)
  list(
    location = list(
      mcmc = list(niter = 1000, nburn = 100, method = "MAR", model = "L",
                  nupd = 1000, print_message = TRUE, hyper = FALSE),
      prior = list(strength = 1, discount = 0.9, m0 = 0, s20 = 1, a0 = 1, b0 = 1),
      output = list(grid = seq(-4, 4, length.out = 1000),
                    out_type = "FULL", out_param = TRUE)),
    location_scale_short = list(
      mcmc = list(niter = 2000, nburn = 1000, method = "MAR", model = "LS", hyper = TRUE),
      prior = list(m0 = 0, k0 = 0.2, a0 = 3, b0 = 2),
      output = list(grid = local_grid, out_type = "FULL", out_param = TRUE)),
    location_scale_long = list(
      mcmc = list(niter = 5000, nburn = 2500, method = "MAR", model = "LS", hyper = TRUE),
      prior = list(m0 = 0, k0 = 0.2, a0 = 2, b0 = 2,
                   m1 = 0, s21 = stats::var(observations), tau1 = 1, zeta1 = 1,
                   a1 = 1, b1 = 1),
      output = list(grid = local_grid, out_type = "FULL", out_param = TRUE))
  )
}

# Input is already standardized; package burn-in is applied within PYdensity.
fit_bnpmix <- function(standardized_observations, profile = "location") {
  profiles <- bnpmix_profiles(standardized_observations)
  if (!profile %in% names(profiles)) stop("Unknown BNPmix profile.")
  settings <- profiles[[profile]]
  do.call(BNPmix::PYdensity, c(list(y = standardized_observations), settings))
}

# VI/Binder package partitions are separate from the sampled-candidate Dahl summary.
bnpmix_partition <- function(fit, loss = c("VI", "Binder")) {
  BNPmix::partition(fit, dist = match.arg(loss))
}
