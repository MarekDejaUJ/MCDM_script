# MCDM_script
A Multi-criteria decision model for AI implementation in academic libraries

# AI readiness in academic libraries
**Fuzzy Delphi → Bayesian BWM → Sensitivity → Credal rankings**

This repository contains data and R code for a ten-step workflow that turns first-round expert judgments into uncertainty-aware importance weights and a compromise ranking. Stage 1 performs a single-round fuzzy Delphi to screen indicators while preserving dissent. Stage 2 estimates panel weights with Bayesian Best–Worst Method (BBWM) using an informative Dirichlet prior centered on the fuzzy-Delphi consensus, then generates sensitivity scenarios and credal dominance graphs. Stage 3 in the paper uses fuzzy VIKOR; inputs for that step are produced here. Methods draw on Hashemi Petrudi et al. (2022); Kuo & Chen (2008); Murray et al. (1985); Li, Wang, & Rezaei (2020); Munim et al. (2022); Saner et al. (2022); Opricovic (2011). Tooling: tidyverse (Wickham et al., 2019); rstan (Guo et al., 2015); bayesplot (Gabry et al., 2019); loo (Vehtari et al., 2017); tidybayes (Kay, 2024).

---

## Repository contents (flat layout)

- `ai.csv` — raw first-round expert ratings (wide; one row per expert; columns `c1_v1 … c5_v31` and `id`).  
- `group2.csv` — expert → institutional archetype mapping used for archetype aggregation.  
- `fuzzy_delphi.R` — Stage 1 implementation; produces `delphi_stage.rds` and `fuzzy.csv`.  
- `bwm_hier.R` — Stage 2 (pooled Bayesian BWM by default; optional hierarchical variant toggled inside); produces `bbwm_stage.rds`.  
- `sensitivity.R` — Stage 2a stress tests main-criteria tilts; produces `sensitivity_stage.rds`.  
- `credal.R` — Stage 2b credal dominance graphs from posterior draws; writes `main_credal_ranking.png`.  
- `rhat.R` — quick R-hat/ESS extraction from stored Stan fits.  
- `bwm_dirichlet_lognormal.stan` — pooled BBWM with Delphi-centered Dirichlet prior and lognormal ratio likelihood.  
- `delphi_stage.rds`, `bbwm_stage.rds`, `sensitivity_stage.rds` — staged artifacts for downstream steps.  
- `code_AI_MCDM.Rproj`, `LICENSE` (CC0), `README.md`.

> Reverse coding: selected C4 items are transformed as `5 − x` before fuzzy aggregation so higher is uniformly “better”; see `fuzzy_delphi.R` for the exact columns.

---

## How to run (quick start)

```r
# setwd to repo root or open the .Rproj
source("fuzzy_delphi.R")     # → delphi_stage.rds
source("bwm_hier.R")         # → bbwm_stage.rds  (pooled BBWM by default)
source("sensitivity.R")      # → sensitivity_stage.rds + printed tables
source("credal.R")           # → main_credal_ranking.png
source("rhat.R")             # → R-hat / ESS summaries from stored fits
