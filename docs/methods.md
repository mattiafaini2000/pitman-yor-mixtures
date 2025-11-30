# Method notes

**Parameters and prior partitions**

Public functions use `concentration` for the source's `c`, and `discount` for `d`. The reported process has `0 <= d < 1` and `c > -d`; executable simulation and mixture helpers restrict concentration to positive values. At step `i`, the new-block category comes first with weight `c + d*K`, followed by existing blocks with weights `n_j - d`, divided by `c + i - 1`. The base draw occurs only after the new-block category is selected, including the initialization at `i = 1`. Sample size one returns that initialized partition directly. The named `sample` and `sample_int` variants retain the respective source sampling calls.

`simulate_partition()` returns latent values, appearance-order allocations, final block sizes and the occupied-count path. Appearance labels are not size ranks. `summarize_rank_proportions()` distinguishes conditional averages over replicates where a rank exists from averages that fill missing ranks with zero. The report's count/mass panels retain ranks occurring at least 150 times; later rank-size studies use zero filling and a common truncation across the discount grid. `pooled_cluster_size_distribution()` counts blocks across all replications and normalizes once, rather than averaging individually normalized distributions.

`R/theory.R` separates exact expected occupied counts, expected stick-breaking masses and expected counts of blocks of a specified size from fitted rank/size reference curves. For positive discount, rank references fix slope `-1/d`; cluster-size references fix slope `-(1+d)`. Their amplitudes are estimated from the same empirical tails, beginning at rank 5 or size 3. The `rank_120` variant handles `d = 0` with a separately fitted exponential; `rank_50` has no Dirichlet reference curve. These are fitted illustrations. Direct product/gamma arithmetic and finite truncation remain subject to their usual range limitations.

| Prior study | Settings and summary |
| --- | --- |
| `report_dp`, `report_py` | `n=1000`, 300 replications, `c=1`, `d=0` or `0.3`; conditional count/mass panels; seed 123 |
| `rank_50`, `rank_120` | Ten discounts from 0 to 0.9; zero-filled means through at most 50 or 120 ranks |
| `cluster_sizes` | Pooled block frequencies through size 50; fitted tail amplitude |
| `synthetic_observations` | `n=500`, `c=1`, `d=0.3`, base variance 0.5 and observation standard deviation 0.1 |

Standalone settings without a supplied seed retain `seed = NULL`. Starting a report profile with seed 123 does not recover the original document's shared RNG stream through all chunks. The synthetic configuration records the unused source variable `sigma2 = 1` separately from the actual observation standard deviation.

**Mixture inference and stored states**

The principal model is `X_i | theta_i, sigma2 ~ N(theta_i, sigma2)`, with a Pitman-Yor prior on the locations, Gaussian base measure and inverse-gamma prior on the shared variance. R performs centering and sample-standard-deviation scaling once, and returns both transformation values. The C++ core takes its supplied coordinates directly. The `report` profile uses variance shape 2 and rate 1; `cpp_defaults` uses shape 1 and rate 1. Rates passed to R's gamma draw become reciprocal scales in the Rcpp gamma call. Variance parameters are not standard deviations.

Both implementations retain allocation updates, the shared-variance update and occupied-mean acceleration. `theta_store` is written before acceleration, while `sigma2_store` follows its variance update. All iterations are retained; posterior helpers remove the first `burn_in` rows once. The compatibility argument `nburn` is unused within the C++ sampler. `pi` contains only the last observation's allocation-probability vector in each sweep. Exact shared theta values encode membership. The backends use their respective ordering and sampling routines, without a claim of identical seeded trajectories. Exponentiation before probability normalization can underflow for extreme inputs.

`load_cpp_sampler()` compiles only when explicitly called, with a local `cacheDir`; R must already have started with a temporary directory inside the checkout. No compiled objects are distributed. The third original lab script's external comparison objects are not supplied; the sampler profile preserves its model parameters without inventing those objects.

**Posterior summaries**

The historical density helper averages sums of occupied Gaussian components weighted by undiscounted counts. It omits the base-measure/new-cluster contribution and is not a normalized full Pitman-Yor predictive density. `discount_density` retains normalization by `n + d`; `concentration_density` and the main mixture illustration use `n + c`. Coordinates and densities remain on the standardized scale. Converting a density to original units would require the scale Jacobian.

| Sensitivity study | Grid | Discarded initial iterations |
| --- | --- | ---: |
| `discount_density`, `discount_clusters` | Ten discounts from 0 to 0.9, concentration 1 | 100 |
| `concentration_density` | Ten concentrations from 0.5 to 1.5, discount 0.3 | 100 |
| `concentration_clusters` | Same concentration grid | 200 |

The posterior similarity matrix averages pairwise equality over retained theta rows. Dahl's representative minimizes the sum of squared distances to that matrix among those same sampled rows; `which.min` retains the first minimum. Scalar-loop and outer-product calculations are available. A separate helper retains the earlier report's rounded final-row partition. Cluster labels and colours across separate fits have no matched identity.

Posterior rank-size exploration uses the original synthetic `rnorm(1000, 0, 1)` sample, seed 123, without adding standardization. It is a separate dataset from galaxies. Exploratory comparison utilities accept explicit reference/estimated labels: `1 - sum(p^2)` is impurity/diversity, and joint-contingency entropy is not a normalized clustering-accuracy score.

**Optional package examples and reference artifacts**

`R/bnpmix.R` retains separate location and location-scale profiles. The location example uses `L`, `MAR`, `hyper=FALSE`, strength 1, discount 0.9 and 1,000 iterations with 100 burn-in. The two `LS` examples have `hyper=TRUE`, 2,000/5,000 iterations and 1,000/2,500 burn-in. Their variance prior is specified with the documented `b0` field. The [BNPmix manual](https://cran.r-project.org/web/packages/BNPmix/BNPmix.pdf) identifies `MAR` as the marginal sampler and documents these fields; this API inspection does not establish a historical package version. VI/Binder package partitions are separate estimators from the sampled-candidate Dahl partition. No package agreement or execution is claimed.

The five historical RDS files were inspected only for table schemas: `dens.RDS` and `dens_c.RDS` have 10,000 rows of `x,density,d` or `x,density,c`; `df_data.RDS` has 82 `x` values; `plot_df.RDS` and `plot_df_c.RDS` have 820 rows of `x,d,cluster` or `x,c,cluster`. They are summary snapshots, not full chains. Exact generation settings remain unestablished. Original asset hashes, mappings and source locations are recorded in local, uncommitted inventories.
