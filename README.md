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
```

To switch to the hierarchical expert-level BBWM, edit `bwm_hier.R` to point to the hierarchical Stan file variant (commented in the script). Expect ~10× longer runtime; pooled BBWM remains the default for speed and reproducibility.

Downstream scripts expect pooled parameter name `alpha` (hierarchical exposes the group mean as `w_group`; the scripts handle the name internally).

---

## Research protocol ↔ repo steps (traceability)
### Stage 1: Fuzzy Delphi (screening)

1) Treat each 1–4 Likert response as a degenerate TFN; aggregate per indicator with min (pessimistic), geometric mean (center), max (optimistic) (Hashemi Petrudi et al., 2022; Kuo & Chen, 2008; Murray et al., 1985).
2) Defuzzify each TFN by geometric aggregation $a_i=(l_im_iu_i)^{1/3}$; retain by rule $a_i \ge 0.8\bar{a}$ (modal-point alternative also possible) (Hashemi Petrudi et al., 2020; Wang et al., 2014).
Implemented in: `fuzzy_delphi.R` → `delphi_stage.rds` with W, df_group (l,m,u), df_defuzz (a), sel_vec, and meta.

### Stage 2: Bayesian BWM (weighting)
3) From retained indicators, compute best→others (BO) and others→worst (OW) ratios per expert (Hashemi Petrudi et al., 2022).
4) Center the Dirichlet prior on the normalized Delphi centers $\tilde{m}$ (replacing generic $w^*$) (Munim et al., 2022).
5) Fit pooled BBWM at the domain level (five main criteria) with lognormal likelihood on BO/OW ratios and Dirichlet prior $\alpha \sim Dirichlet(c\tilde{m})$.
6) Fit pooled BBWM within each domain on retained indicators; obtain local posteriors.
7) Compute global indicator weights as $w_{global}=w_{domain}\times w_{local}$; keep posterior draws for uncertainty.
Implemented in: `bwm_hier.R` (pooled by default) → `bbwm_stage.rds` with fit_global, fits_local, summaries, and weights_global_sum (Li et al., 2020; Munim et al., 2022; Saner et al., 2022).

### Stage 2a: Sensitivity (policy tilt)
8) Vary one main criterion (e.g., C1) over a grid 0.10…0.90; proportionally renormalize the others; recompute global indicator weights and re-rank (Talib et al., 2019; Hashemi Petrudi et al., 2022).
Implemented in: `sensitivity.R` → `sensitivity_stage.rds` and wide tables with original vs scenario weights and ranks.

### Stage 2b: Credal rankings (partial orders)
9) Convert posterior draws to pairwise exceedance probabilities $P(w_i > w_j)$; draw directed edges above a confidence threshold to form dominance digraphs (Li et al., 2020; Munim et al., 2022).
Implemented in: `credal.R` → `main_credal_ranking.png` (+ domain-level graphs if desired).

### Stage 3: Fuzzy VIKOR (compromise ranking)
10) Build an archetype × indicator fuzzy decision matrix (min–gmean–max per cell), draw weight vectors from BBWM posteriors, compute fuzzy $S, R,$ and $Q$ and defuzzify; report distributions and the compromise solution with acceptable-advantage and acceptable-stability checks (Opricovic, 2011).
Inputs produced here; implement with your preferred FuzzyMCDM/VIKOR routine if needed.

---

## Data and staged artifacts
- `ai.csv` — raw Likert ratings (1–4). `id` identifies experts; indicators `c1_v1 … c5_v31`.
- `group2.csv` — archetype labels per expert (`id`, `type`).
- `delphi_stage.rds` — list with `W` (experts × retained indicators), `df_group` (TFNs), `df_defuzz` (consensus a), `sel_vec`, `meta`.
- `bbwm_stage.rds` — list with `fit_global`, `fits_local`, `alpha_global_sum`, `alpha_locals_sum`, `weights_global_sum`, `domains`, `indicators_local`.
- `sensitivity_stage.rds` — tables of scenario main-weights and indicator weights/ranks.

All scripts use relative paths; staged `.rds` files allow resuming without recomputation. Seeds are set where applicable.

---

## Software requirements
- R ≥ 4.x
- **CRAN**: `tidyverse`, `rstan`, `bayesplot`, `loo`, `tidybayes`, `igraph`, `tidygraph`, `ggraph`, `ggrepel`
- **Stan toolchain** configured for `rstan` (see platform-specific setup).
- **Parallel chains**: `options(mc.cores = parallel::detectCores())`. Tune `adapt_delta` and `max_treedepth` in scripts if you see divergent transitions or depth saturations.

---

## OSF connection and citation
This GitHub repository is connected to the Open Science Framework via the GitHub add-on. The OSF project hosts a public copy of these materials and can be registered to mint a DOI for archival citation. Add the OSF DOI here when issued. If you archive a GitHub release on Zenodo, add that DOI here as the software citation. In manuscripts, cite the paper, the OSF dataset DOI, and the Zenodo/GitHub software DOI as appropriate.

**Citation template**
Author(s). Year. AI readiness in academic libraries: fuzzy Delphi, Bayesian BWM, sensitivity, and credal rankings. OSF, DOI: [OSF DOI]. Software release, DOI: [Zenodo DOI].

---

## License
All files are released under Creative Commons CC0 1.0 Universal (public-domain dedication). You may copy, modify, distribute, and perform the work, including for commercial purposes, without permission. Scholarly attribution is appreciated.

---

## References (methods)
Hashemi Petrudi et al., 2020; Hashemi Petrudi et al., 2022; Kuo & Chen, 2008; Li, Wang, & Rezaei, 2020; Murray et al., 1985; Munim et al., 2022; Opricovic, 2011; Saner et al., 2022; Wickham et al., 2019; Guo et al., 2015; Gabry et al., 2019; Vehtari et al., 2017; Kay, 2024.
