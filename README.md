# A Bayesian flexible hypothesis test for multiple comparisons

**Authors**  
Ingrid Guevara, Vanda Inacio, Antonio Lijoi, and Luis Gutiérrez

This repository contains R examples for Bayesian multiple comparisons of treatment distributions using Dirichlet process mixtures. Each example simulates data for three treatments and estimates posterior probabilities for partitions of those treatments. Treatments in the same block share a distribution under that partition.

## Repository files

| File | Purpose |
| --- | --- |
| [`dpm_univariate_normal.R`](dpm_univariate_normal.R) | Simulates univariate normal mixtures and plots the true and estimated densities with pointwise 50% credible bands. |
| [`dpm_multivariate_normal.R`](dpm_multivariate_normal.R) | Simulates bivariate normal mixtures and displays one faceted plot with the simulated observations and the true and estimated density contours. |
| [`partition_proposals.R`](partition_proposals.R) | Provides the shared partition prior, log-Stirling numbers, move/split/merge proposals, and likelihood and plotting helpers. |

## Run the examples

Install the required packages in R.

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "gsl", "mvtnorm", "MASS", "CholWishart"))
```

Open R or RStudio with the repository root as your working directory. Run either example with the corresponding command below. The examples simulate their own data, and the seeds make the simulation and MCMC draws reproducible in the same R environment.

```r
set.seed(123)
source("dpm_univariate_normal.R")

set.seed(123)
source("dpm_multivariate_normal.R")
```

Each example defaults to 10,000 MCMC iterations and uses the second half of the chain for the posterior model table and predictive densities. To try a shorter run, change `numMCMC` in the **MCMC sampler** section.

## Read the output

- The console prints the posterior probabilities of the sampled partitions. A label such as `1,2,1` groups treatments 1 and 3 together and places treatment 2 in a separate block.
- The plot opens in the interactive R graphics window. Treatments share a color when they belong to the same block of the highest-posterior partition. Line type distinguishes the true densities from the estimates.
- The estimated densities average over the retained MCMC draws, so treatments shown in the same color can still have different estimated curves.
