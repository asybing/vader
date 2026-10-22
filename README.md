# VADER

**Vedotin ADC Dose Escalation in R**

VADER is an R-based clinical trial simulation framework for evaluating dose-escalation strategies for vedotin antibody-drug conjugates (ADCs). It links patient-level bodysize simulation, population pharmacokinetic (PopPK) simulation, an exposure–dose limiting toxicity (DLT) model, and adaptive dose-escalation rules within a single workflow.

The package is intended to support simulation of individual dose-escalation trials and Monte Carlo evaluation of operating characteristics across repeated trials. For most users, the two primary functions are:

- `simde()` — simulate one dose-escalation trial.
- `simde_mc()` — run repeated dose-escalation trials and summarize performance across replicates.

All other package functions should generally be treated as internal helper functions that support these two user-facing workflows.

## Scientific context

VADER was developed to evaluate how accumulated clinical and pharmacokinetic knowledge from the vedotin ADC platform can inform first-in-human dose escalation. The underlying workflow was designed to:

1. define the candidate dose levels and trial design;
2. simulate individual patient characteristics;
3. simulate individual ADC and/or MMAE exposure using a PopPK model;
4. simulate DLT occurrence from a prespecified exposure–DLT relationship;
5. update dose-level summaries after each newly enrolled participant;
6. apply the selected dose-escalation rules;
7. continue enrollment until a stopping criterion is met; and
8. summarize the selected maximum tolerated dose (MTD) or maximum administered dose (MAD), enrollment, DLTs, and other trial-level outcomes.

The framework can be used to compare 3+3, Bayesian optimal interval (BOIN), and Bayesian logistic regression model (BLRM) designs.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("asybing/vader")
```

Load the package:

```r
library(vader)
```

## Recommended workflow

### 1. Define the simulation scenario

Before running a simulation, specify the major components of the virtual trial:

- candidate dose levels
- dose units;
- dosing schedule (ii and addl);
- DLT assessment window;
- drug–antibody ratio (DAR);
- popPK model;
- PK exposure metric;
- target toxicity rate;
- maximum cohort size;
- maximum total trial enrollment;
- dose-escalation design (3+3, BOIN, or BLRM);
- use of an informative prior (only needed for BLRM).

### 2. Simulate a single trial with `simde()`

Use `simde()` when you want to inspect the conduct and outcome of one virtual dose-escalation trial.

```r
trial <- simde(
  # candidate doses
  # dose units
  # dosing schedule
  # DAR
  # PopPK model
  # exposure–DLT model
  # escalation design
  # cohort and trial limits
)

trial
```

`simde()` performs the following cycle:

1. enroll the next participant or cohort;
2. simulate individual body-size characteristics or other required covariates;
3. simulate individual exposure using the PopPK model;
4. simulate DLT occurrence using the exposure–DLT relationship;
5. update the cumulative data at each dose level;
6. classify dose levels according to the selected design;
7. determine whether to move **up**, move **down**, **stay**, or **stop** (udss); and
8. repeat until the trial meets a stopping rule.

A single-trial simulation is useful for checking whether the model inputs and escalation behavior are clinically reasonable before running a larger Monte Carlo experiment.

### 3. Run repeated trials with `simde_mc()`

Use `simde_mc()` to repeat the trial simulation under the same scenario and evaluate the design's operating characteristics.

```r
set.seed(123)

mc_results <- simde_mc(
  # the same core scenario inputs used by simde()
  # number of simulation replicates
)

mc_results
```

Typical operating characteristics of interest include:

- probability of selecting the model-defined target or “true” MTD;
- total number of participants enrolled;
- number of participants experiencing a DLT;
- number of participants enrolled above the 'true' MTD.

## `simde()` versus `simde_mc()`

| Function | Use it when you want to... | Typical output level |
|---|---|---|
| `simde()` | simulate and inspect one virtual dose-escalation trial | participant, cohort, dose-level, and trial results |
| `simde_mc()` | evaluate a design across many repeated virtual trials | replicate-level results and aggregate operating characteristics |

A practical development sequence is:

```text
Define scenario -> run simde() -> inspect behavior -> refine inputs -> run simde_mc() -> compare designs
```

## Helper functions

The remaining functions in VADER support the internal steps required by `simde()` and `simde_mc()`, such as patient-covariate generation, PopPK simulation, exposure calculation, DLT generation, dose-level summarization, tolerability assessment, escalation decisions, and result formatting.

Users generally should not need to call these helper functions directly. Keeping analyses centered on `simde()` and `simde_mc()` provides a more stable, reproducible interface and reduces the risk of bypassing required workflow steps.

Advanced users and package developers may inspect helper-function documentation or source code when extending or validating the package, but helper functions should not be considered the primary public API unless explicitly documented otherwise.

## Dose-escalation logic

VADER supports a common simulation architecture across multiple designs while allowing the decision rules to differ by approach.

- **3+3:** escalation and stopping decisions are based on the number of observed DLTs among participants treated at the current dose. Optional backfill may allow additional enrollment at tolerable dose levels.
- **BOIN:** the observed toxicity rate is compared with prespecified lower and upper decision boundaries to determine escalation, continued enrollment, or de-escalation.
- **BLRM:** dose-level toxicity probabilities are updated as clinical data accumulate. Escalation decisions may incorporate an informative or weakly informative prior and overdose control criteria.
