# MCDM for AI readiness
A Multi-criteria decision model for AI implementation in academic libraries

# AI readiness in academic libraries
**Fuzzy Delphi → Bayesian BWM → Sensitivity → Credal rankings**

This repository contains data and R code for a ten-step workflow that turns first-round expert judgments into uncertainty-aware importance weights and a compromise ranking. Stage 1 performs a single-round fuzzy Delphi to screen indicators while preserving dissent. Stage 2 estimates panel weights with Bayesian Best–Worst Method (BBWM) using an informative Dirichlet prior centered on the fuzzy-Delphi consensus, then generates sensitivity scenarios and credal dominance graphs. Stage 3 in the paper uses fuzzy VIKOR; inputs for that step are produced here. Methods draw on Hashemi Petrudi et al. (2022); Kuo & Chen (2008); Murray et al. (1985); Li, Wang, & Rezaei (2020); Munim et al. (2022); Saner et al. (2022); Opricovic (2011). Tooling: tidyverse (Wickham et al., 2019); rstan (Guo et al., 2015); bayesplot (Gabry et al., 2019); loo (Vehtari et al., 2017); tidybayes (Kay, 2024).

---

## Repository contents (flat layout)

- `ai.csv` — raw first-round expert ratings (wide; one row per expert; columns `c1_v1 … c5_v31` and `id`).  
- `group2.csv` — expert → institutional archetype mapping used for archetype aggregation.  
- `fuzzy_delphi.R` — Stage 1 implementation; produces `delphi_stage.rds`.  
- `bwm_hier.R` — Stage 2 (pooled Bayesian BWM by default; optional hierarchical variant toggled inside); produces `bbwm_stage.rds`.  
- `sensitivity.R` — Stage 2a stress tests main-criteria tilts; produces `sensitivity_stage.rds`.  
- `credal.R` — Stage 2b credal dominance graphs from posterior draws; writes `main_credal_ranking.png`.  
- `rhat.R` — quick R-hat/ESS extraction from stored Stan fits.
- `fuzzy_VIKOR.R` — Stage 3; produces `vikor_stage.rds` in `results_vikor` folder.
- `bwm_dirichlet_lognormal.stan` — pooled BBWM with Delphi-centered Dirichlet prior and lognormal ratio likelihood.  
- `delphi_stage.rds`, `bbwm_stage.rds`, `sensitivity_stage.rds` — staged artifacts for downstream steps.  
- `code_AI_MCDM.Rproj`, `LICENSE` (CC0), `README.md`.

> Reverse coding: selected C4 items are transformed as `5 − x` before fuzzy aggregation so higher is uniformly “better” in case of preperdness; see `fuzzy_delphi.R` for the exact columns.

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
-  **Fuzzy VIKOR**: `FuzzyMCDM` v1.1 (archived). Install via:
  ```r
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_version("FuzzyMCDM", version = "1.1", repos = "https://cran.r-project.org")
  # or:
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/FuzzyMCDM/FuzzyMCDM_1.1.tar.gz",
    repos = NULL, type = "source"
  )
  library(FuzzyMCDM)
```
Or via Github: https://github.com/cran/FuzzyMCDM
---

## OSF connection and citation
This GitHub repository is connected to the Open Science Framework via the GitHub add-on. The OSF project hosts a public copy of these materials and can be registered to mint a DOI for archival citation. 

**Citation template**
Deja, M. (2025, September 9). A Multi-criteria decision model for AI implementation in academic libraries. Retrieved from osf.io/2rnkv


---

## License
All files are released under Creative Commons CC0 1.0 Universal (public-domain dedication). You may copy, modify, distribute, and perform the work, including for commercial purposes, without permission. Scholarly attribution is appreciated.

---

## References (methods)
Hashemi Petrudi, S. H., Ghomi, H., & Mazaheriasad, M. (2022). An Integrated Fuzzy Delphi and Best Worst Method (BWM) for performance measurement in higher education. Decision Analytics Journal, 4, 100121. https://doi.org/10.1016/j.dajour.2022.100121
Hashemi Petrudi, S. H., Tavana, M., & Abdi, M. (2020). A comprehensive framework for analyzing challenges in humanitarian supply chain management: A case study of the Iranian Red Crescent Society. International Journal of Disaster Risk Reduction, 42, 101340. https://doi.org/10.1016/j.ijdrr.2019.101340
Kuo, Y.-F., & Chen, P.-C. (2008). Constructing performance appraisal indicators for mobility of the service industries using Fuzzy Delphi Method. Expert Systems with Applications, 35(4), 1930–1939. https://doi.org/10.1016/j.eswa.2007.08.068
Li, L., Wang, X., & Rezaei, J. (2020). A Bayesian Best-Worst Method-Based Multicriteria Competence Analysis of Crowdsourcing Delivery Personnel. Complexity, 2020, 1–17. https://doi.org/10.1155/2020/4250417
Munim, Z. H., Balasubramaniyan, S., Kouhizadeh, M., & Ullah Ibne Hossain, N. (2022). Assessing blockchain technology adoption in the Norwegian oil and gas industry using Bayesian Best Worst Method. Journal of Industrial Information Integration, 28, 100346. https://doi.org/10.1016/j.jii.2022.100346
Murray, T. J., Pipino, L. L., & Van Gigch, J. P. (1985). A pilot study of fuzzy set modification of Delphi*. Human Systems Management, 5(1), 76–80. https://doi.org/10.3233/HSM-1985-5111
Opricovic, S. (2011). Fuzzy VIKOR with an application to water resources planning. Expert Systems with Applications, 38(10), 12983–12990. https://doi.org/10.1016/j.eswa.2011.04.097
Opricovic, S., & Tzeng, G.-H. (2004). Compromise solution by MCDM methods: A comparative analysis of VIKOR and TOPSIS. European Journal of Operational Research, 156(2), 445–455. https://doi.org/10.1016/S0377-2217(03)00020-1
Saner, H. S., Yucesan, M., & Gul, M. (2022). A Bayesian BWM and VIKOR-based model for assessing hospital preparedness in the face of disasters. Natural Hazards, 111(2), 1603–1635. https://doi.org/10.1007/s11069-021-05108-7
Talib, F., Asjad, M., Attri, R., Siddiquee, A. N., & Khan, Z. A. (2019). Ranking model of total quality management enablers in healthcare establishments using the best-worst method. The TQM Journal, 31(5), 790–814. https://doi.org/10.1108/TQM-04-2019-0118
