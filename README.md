# BePhyNE

**Bayesian environmental Phylogenetic Niche Estimation in R**

BePhyNE jointly estimates species' environmental niches and their evolution across
a phylogeny. It combines presence–absence or presence–pseudoabsence records,
environmental predictors, and a phylogenetic tree in a Bayesian hierarchical
model. Related species inform one another's niche estimates through a shared
evolutionary model.

For each species and environmental predictor, BePhyNE fits a symmetric, unimodal
response curve described by an **optimum**, **breadth**, and **tolerance**
(curve height). The package provides functions for MCMC fitting, posterior
summaries, response-curve plots, predictive evaluation, and ancestral niche
reconstruction through continuous stochastic mapping.

[Full guide](vignettes/BePhyNE.pdf) ·
[Vignette source](vignettes/BePhyNE.Rmd) ·
[Example data](data/Pleth_data_vignette.csv) ·
[Report an issue](https://github.com/sean-mchugh/BePhyNE/issues)

## The model

BePhyNE connects two parts of an analysis:

- **Environmental niche estimation:** occurrence records inform species-specific
  response curves along environmental gradients.
- **Phylogenetic inference:** niche optimum and log breadth evolve together under
  a multivariate Brownian-motion model, with estimated ancestral root values,
  evolutionary rates, and covariance.

In the current fitting workflow, each predictor has its own evolutionary model.
For example, precipitation optimum and log breadth can covary, but precipitation
optimum and temperature optimum have no estimated evolutionary covariance.
Predictors contribute additively on the logit scale, without interaction terms.
Tolerance is estimated with a species-level prior and is not part of the
two-trait Brownian-motion model used here.

These are model-based estimates of environmental response. Their interpretation
depends on the occurrence sample, how absences or pseudoabsences were chosen, the
predictors, the phylogeny, and the priors. Curve height should not automatically
be interpreted as prevalence or physiological tolerance.

## Installation

Use a recent version of R, at least **R 4.1**: the current source contains the
native `|>` pipe, although `DESCRIPTION` still lists an older minimum version.
Install the development package from GitHub:

```r
install.packages("devtools")  # Once, if needed
devtools::install_github("sean-mchugh/BePhyNE", build_vignettes = FALSE)

library(BePhyNE)
```

Package dependencies are listed in [DESCRIPTION](DESCRIPTION). The installation
above skips rebuilding the vignette; the [PDF guide](vignettes/BePhyNE.pdf) is
available directly in this repository. Optional packages for stochastic mapping
are listed below.

## Input data

You need a phylogeny of class `phylo` with branch lengths and a data frame with
one row per occurrence or absence record. The example data have these columns:

| Column | Contents |
| --- | --- |
| `species` | Species name matching a tip label in `tree$tip.label` exactly |
| `PA` | `1` for presence; `0` for absence or pseudoabsence |
| `bio12` | Annual precipitation at the sampled location |
| `bio1` | Annual mean temperature at the sampled location |

Each environmental predictor occupies its own numeric column. Coordinates and
other columns can remain in the data frame; select the predictors to fit with
`env_preds`. Check predictor units, missing values, and species names before
fitting. The example below assumes that the modeled species have usable presence
and absence records, with no missing values in those records.

`format_BePhyNE_data()` orders the species to match the tree and standardizes the
predictors by default. Keep its returned `scale` object: it records the centering
and scaling needed to label response curves in the original predictor units.
You can supply an existing scaling object through `scale_atr`, for example when
using a common environmental reference dataset.

## Example: plethodontid salamanders

This example follows the vignette using eastern North American plethodontid
salamanders and two climate predictors. It uses the packaged tree and downloads
the example CSV from this repository.

**The short chain below demonstrates the workflow. It is not a sufficient basis
for biological inference.** Choose run length, burn-in, and tuning after
examining convergence and effective sample sizes.

### 1. Prepare the data

```r
library(BePhyNE)
set.seed(123)

data("ENA_Pleth_Tree", package = "BePhyNE")
tree <- ENA_Pleth_Tree

pa_data <- read.csv(
  "https://raw.githubusercontent.com/sean-mchugh/BePhyNE/main/data/Pleth_data_vignette.csv",
  stringsAsFactors = FALSE
)

env_preds <- c("bio12", "bio1")
Npred <- length(env_preds)

data_obj <- format_BePhyNE_data(
  pa_data = pa_data,
  tree = tree,
  sp_col = "species",
  occ_col = "PA",
  env_preds = env_preds
)

scale_atr <- data_obj$scale
sets <- separate.data(data_obj$data, ratio = 0.5)
```

From a local copy of the repository, you can read
`"data/Pleth_data_vignette.csv"` instead of the URL. To use your own observations,
replace the CSV and load your tree with, for example,
`ape::read.tree("your_tree.tre")`.

`separate.data()` randomly splits each species' presence and absence records
between training and prediction sets. Here, half are used for fitting and half
for evaluation. Check that both sets contain enough records of both classes.
For analyses requiring spatially independent evaluation, prepare an appropriate
spatial split yourself.

### 2. Set priors, starting values, and proposal tuning

```r
Prior_scale <- make_all_priors(
  N = Npred,
  tips = length(tree$tip.label),
  heights_mean_by_sp = 0.95,
  heights_sd_by_sp = 0.15,
  r = 2,
  p = 1,
  plot = FALSE
)

startPars <- get_starting_values(
  Prior_scale = Prior_scale,
  tree = tree,
  data = sets$training,
  reps_before_POE = 1000
)

move_details <- make_tuning(tree = tree, pred = Npred)
```

This uses the constructor's default root and evolutionary priors and explicitly
sets the tolerance-prior arguments. These are starting choices for the example,
not universal recommendations. Review `?make_all_priors` and choose priors suited
to your scaled predictors and tree. Here, `N` is the number of environmental
predictors; `r = 2` is the number of evolving niche traits, optimum and log breadth.

Starting-value generation searches for usable likelihoods and can require many
draws. `reps_before_POE` controls when it switches from full draws to replacing
problematic species' starting values; it is not an overall runtime limit.

### 3. Run the MCMC and save the log

```r
out_dir <- "bephyne-output"
dir.create(out_dir, showWarnings = FALSE)
filename <- file.path(out_dir, "plethodontid-demo")

# Use a fresh prefix so the log read below belongs to this run.
if (file.exists(paste0(filename, ".pars.log"))) {
  stop("Choose a new filename prefix before running another chain.")
}

iterations <- 1000

fit <- BePhyNE_MCMC(
  tree = tree,
  pa_data = sets$training,
  Prior_scale = Prior_scale,
  startPars = startPars,
  move_details = move_details,
  iterations = iterations,
  trim_freq = 10,
  write2file = TRUE,
  append2existingfile = FALSE,
  filename = filename
)
```

The sampler writes `bephyne-output/plethodontid-demo.pars.log`, retaining states
every `trim_freq` iterations. File output is required for the log-reading steps
below. Save the tree, predictor order, scaling object, priors, seed, and package
version with your analysis so the results can be interpreted later.

### 4. Summarize and plot the posterior

```r
logdf <- read_BePhyNE_log(paste0(filename, ".pars.log"))

# Illustrative only: choose burn-in from chain diagnostics.
burnin <- floor(0.2 * iterations)
posterior <- logdf[logdf$Iteration > burnin, , drop = FALSE]

log_summary <- summarize_logdf(posterior, HPD_prob = 0.95)

head(log_summary$ESS)
head(log_summary$HPD)
log_summary$median_parlist$traits[[1]][[1]]

plot_summary_ridgeplot(
  tree = tree,
  log_summary = log_summary,
  predictor_names = env_preds,
  scale_atr = scale_atr
)
```

The summary includes effective sample sizes (`ESS`), highest posterior density
intervals (`HPD`), and posterior parameter summaries (`median_parlist`). Raw log
values and `HPD` intervals use the sampler's parameterization, including log
breadth and transformed tolerance. The trait matrices in `median_parlist` contain
back-transformed response parameters.

This example keeps the summary in standardized predictor units for prediction
and passes `scale_atr` to the plotting function to label the environmental axes.
Keep that distinction when comparing curves or supplying new prediction data.
The ridge plot displays curves based on posterior median parameters.

Before interpreting results, inspect traces, compare independently initialized
chains, assess effective sample sizes, and examine sensitivity to priors and
proposal tuning. A completed run alone does not establish convergence.

### 5. Evaluate predictions

```r
prediction_stats <- AUC_posterior_median(
  log_summary = log_summary,
  pa_data = sets$predicting
)

plot_AUC_treebarplot(tree, prediction_stats)
```

These functions calculate and plot species-level AUC statistics using the
posterior median response parameters and the held-out observations. This
evaluation does not integrate over posterior uncertainty; its interpretation
also depends on how the test records and pseudoabsences were sampled.

## Optional: ancestral niche reconstruction

Continuous stochastic mapping uses Bruce Stagg Martin's
[`contsimmap`](https://github.com/bstaggmartin/contsimmap) and
[`evorates`](https://github.com/bstaggmartin/evorates) packages, as described in the
vignette. Install them separately for this step:

```r
devtools::install_github("bstaggmartin/evorates")
devtools::install_github("bstaggmartin/contsimmap")

library(evorates)
library(contsimmap)

set.seed(456)
maps <- make_simmaps_BePhyNE(
  tree = tree,
  logdf = posterior,
  char_names = c("bio12_optimum", "bio12_breadth",
                 "bio1_optimum", "bio1_breadth"),
  nsims = 10
)

plot_BePhyNE_simmap(maps)
```

Supply the retained posterior log, with burn-in already removed. Character names
follow predictor order, with optimum then breadth for each predictor. Ten maps
are a demonstration; use more posterior samples for an analysis, keeping `nsims`
no larger than the number of available rows. This plotting call leaves values
on the fitted predictor scale.

## Documentation and help

Read the [guide](vignettes/BePhyNE.pdf) for the longer walkthrough and use R's help
pages for argument details:

```r
?format_BePhyNE_data
?make_all_priors
?BePhyNE_MCMC
?summarize_logdf
```

The [vignette source](vignettes/BePhyNE.Rmd) records the original workflow and
contains local paths and some older examples. The examples above use the current
exported function interfaces. For questions or bug reports, open a
[GitHub issue](https://github.com/sean-mchugh/BePhyNE/issues) with a small example,
the exact error, and `sessionInfo()`.

## Reference

McHugh, S. W., Espíndola, A., White, E., and Uyeda, J. C. (2022).
[Jointly Modeling Species Niche and Phylogenetic Model in a Bayesian Hierarchical
Framework](https://doi.org/10.1101/2022.07.06.499056). *bioRxiv*, preprint.

When reporting an analysis, also record the BePhyNE version and the GitHub commit
used to install it.

## License

The `License` field in [DESCRIPTION](DESCRIPTION) is currently a placeholder;
a software license has not yet been specified in this package.
