# Pitman-Yor mixtures

A statistics project by **Leonardo Arrigoni and Mattia Faini** on Pitman-Yor processes, prior random partitions, and Gaussian mixture density estimation and clustering.

Read the full [English report](pitman_yor_report_en.pdf), *The Pitman-Yor Process and Pitman-Yor Process Mixture Models*.

## Report summary

The report introduces the Pitman-Yor process through stick-breaking, its predictive distribution and its exchangeable partition probability function. The discount parameter extends the Dirichlet-process special case: the expected number of occupied clusters grows logarithmically for the Dirichlet process and at a power-law rate for positive discount. The project examines both the order in which clusters appear and their ordering by size, which describe different aspects of the random partition.

The reported prior studies compare empirical cluster counts and block proportions with analytical expressions and asymptotic references, using 300 replications in the illustrated settings. Increasing discount produces more occupied clusters and slower decay of size-ranked proportions. Separate rank-size and cluster-size illustrations use fitted reference amplitudes; their curves should be read as empirical comparisons rather than independent theoretical validation.

For the standardized `MASS::galaxies` observations, the report studies a Gaussian location mixture with shared variance using marginal Gibbs sampling. Its historical figures show broadly similar density curves across the illustrated discount and concentration grids, while representative partitions become more fragmented as these parameters increase. Partitions are selected from sampled states using squared distance to the posterior similarity matrix. The figures summarize runs of 1,000 Gibbs iterations; they do not establish convergence or independent Monte Carlo replications.

## Code and dependencies

| Location | Contents |
| --- | --- |
| `R/` | Prior simulation, formulas, R sampler, posterior summaries, plotting and file handling |
| `src/gibbs.cpp` | Optional Rcpp sampler and original-name compatibility wrapper |
| `scripts/`, `config/` | Explicit studies and their numerical settings |
| `results/reference/` | Five unchanged historical summary tables, not complete chains |
| `pitman_yor_report_en.pdf` | Complete English report |

Prior simulation and formulas use base R and `stats`. The R sampler and historical density summaries require `plyr`; the prior plotting entry point requires `ggplot2`. Galaxies examples use `MASS`, supplied with standard R installations. The optional C++ backend requires `Rcpp` and an existing compiler toolchain; the package comparison requires `BNPmix`. [r_dependencies.csv](r_dependencies.csv) records these roles without asserting historical package versions.

The dataset is accessed as `MASS::galaxies`, an 82-element velocity vector attributed in the [MASS documentation](https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/galaxies.html) to Roeder (1990) and Postman, Huchra and Geller (1986). No external data download is required. The saved tables contain standardized coordinates, density summaries and representative labels; their exact generating configurations and full posterior chains are unavailable.

## Example commands

Run from the checkout root with dependencies already installed. Set temporary paths before starting R, particularly for the optional Rcpp backend:

```powershell
New-Item -ItemType Directory -Force .cache/tmp | Out-Null
$env:TMP = $env:TEMP = $env:TMPDIR = (Resolve-Path .cache/tmp).Path
Rscript --vanilla scripts/prior_study.R report_dp
Rscript --vanilla scripts/prior_study.R rank_120
Rscript --vanilla scripts/galaxies_mixture.R r report
Rscript --vanilla scripts/parameter_sensitivity.R discount_density rcpp
Rscript --vanilla scripts/posterior_rank_size.R rcpp
Rscript --vanilla scripts/bnpmix_comparison.R location
```

Generated tables, study objects and prior plots go under `outputs/`. Configuration distinguishes the report's variance prior, C++ defaults, burn-in counts and density normalizations. [Method notes](docs/methods.md) explain these settings and the other named variants. Sourcing reusable modules does not start a study or compile C++.
