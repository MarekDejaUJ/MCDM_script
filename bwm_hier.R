# bwm_hier.R — hierarchical Bayesian BWM (global + local) with fuzzy‑Delphi priors

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(rstan)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 0) Config
config <- list(
  delphi_rds     = "delphi_stage.rds",
  stan_file      = "bwm_dirichlet_lognormal.stan",
  out_rds        = "bbwm_stage.rds",
  seed           = 20250908,
  chains         = 4,
  iter           = 6000,
  warmup         = 1000,
  adapt_delta    = 0.95,
  max_treedepth  = 12,
  c_prior        = NA_real_   # if NA, will default to #experts
)

# 1) Load fuzzy‑Delphi artifacts
stopifnot(file.exists(config$delphi_rds))
D <- readRDS(config$delphi_rds)
W               <- D$W                     # experts × indicators
meta            <- D$meta                  # criterion, domain, reverse, direction
sel_vec         <- D$sel_vec               # retained indicators
df_group        <- D$df_group              # per‑criterion panel TFNs (l,m,u)
normalize       <- function(x) x / sum(x)
gmean           <- function(x) exp(mean(log(x)))
set.seed(config$seed)

# 2) Subset to retained criteria (keeps column order from sel_vec)
W_cons <- W[, sel_vec, drop = FALSE]
stopifnot(ncol(W_cons) == length(sel_vec))

# 3) Domains and lookups
domains <- paste0("c", 1:5)
crit_dom <- tibble(criterion = sel_vec) %>% left_join(meta[, c("criterion","domain")], by = "criterion")

# 4) Domain‑level fuzzy priors (geometric mean of m within each retained domain)
df_fuzzy_main <- df_group %>%
  filter(criterion %in% sel_vec) %>%
  mutate(domain = str_extract(criterion, "^(c[1-5])")) %>%
  group_by(domain) %>%
  summarize(m_fuzzy = gmean(m), .groups = "drop") %>%
  arrange(match(domain, domains)) %>%
  mutate(m_fuzzy = normalize(m_fuzzy))

# 5) Build W_main (items × experts = 5 × J): per domain, average each expert's scores across its selected indicators
W_main <- map(domains, function(dom){
  cols <- crit_dom %>% filter(domain == dom) %>% pull(criterion)
  if (length(cols) == 0) return(rep(NA_real_, nrow(W_cons)))
  rowMeans(W_cons[, cols, drop = FALSE], na.rm = TRUE)
}) %>%
  set_names(domains) %>%
  bind_cols() %>%
  as.matrix() %>%
  t()

# 6) Helper to prepare BWM ratios (items × experts → list)
make_bwm <- function(mat_items_by_experts){
  M <- mat_items_by_experts
  M[is.na(M)] <- 0
  M[M <= 0] <- 1e-2
  m <- nrow(M); J <- ncol(M)
  best  <- apply(M, 2, which.max)
  worst <- apply(M, 2, which.min)
  bos <- ows <- matrix(NA_real_, J, m)
  for (j in seq_len(J)) for (i in seq_len(m)){
    bos[j, i] <- M[best[j], j] / M[i, j]
    ows[j, i] <- M[i, j]       / M[worst[j], j]
  }
  list(m = m, J = J, bos = bos, ows = ows, best = best, worst = worst)
}

# 7) Write Stan model if missing (separate file)
stan_model_code <- "
// bwm_dirichlet_lognormal.stan — BBWM with Dirichlet center and lognormal ratios
// Data are built from best-to-others and others-to-worst ratios

data {
  int<lower=1> m;               // number of items
  int<lower=1> J;               // number of experts
  matrix[J, m] bos;             // Best-to-Other ratios
  matrix[J, m] ows;             // Other-to-Worst ratios
  int<lower=1,upper=m> best[J]; // index of best item per expert
  int<lower=1,upper=m> worst[J];// index of worst item per expert
  vector<lower=0>[m] m_fuzzy;   // normalized fuzzy-Delphi centers (sum to 1)
  real<lower=0> c;              // Dirichlet concentration parameter
}
parameters {
  simplex[m] alpha;             // item weights
  real<lower=0> sigma_b;        // noise for best-to-other
  real<lower=0> sigma_w;        // noise for other-to-worst
}
model {
  // Prior
  alpha ~ dirichlet(c * m_fuzzy);
  sigma_b ~ normal(0, 1);
  sigma_w ~ normal(0, 1);
  // Likelihood
  for (j in 1:J) {
    for (i in 1:m) {
      bos[j, i] ~ lognormal(log(alpha[best[j]] / alpha[i]), sigma_b);
      ows[j, i] ~ lognormal(log(alpha[i] / alpha[worst[j]]),        sigma_w);
    }
  }
}
"
if (!file.exists(config$stan_file)) writeLines(stan_model_code, config$stan_file)

# 8) Fit global (5 domains)
prep_global <- make_bwm(W_main)
stopifnot(prep_global$m == length(domains))
J <- prep_global$J
c_prior <- if (is.na(config$c_prior)) J else config$c_prior

data_global <- c(prep_global, list(m_fuzzy = df_fuzzy_main$m_fuzzy, c = c_prior))
fit_global <- rstan::stan(
  file    = config$stan_file,
  data    = data_global,
  chains  = config$chains,
  iter    = config$iter,
  warmup  = config$warmup,
  control = list(adapt_delta = config$adapt_delta, max_treedepth = config$max_treedepth),
  seed    = config$seed
)

Aglob <- as.matrix(fit_global, pars = "alpha")
alpha_global_summary <- tibble(
  domain = domains,
  mean   = colMeans(Aglob),
  median = apply(Aglob, 2, median),
  q05    = apply(Aglob, 2, quantile, 0.05),
  q95    = apply(Aglob, 2, quantile, 0.95)
)

# 9) Fit locals per domain
fits_local <- list()
indicators_local <- list()
alpha_locals_summary <- map_dfr(domains, function(dom){
  cols <- crit_dom %>% filter(domain == dom) %>% pull(criterion)
  if (!length(cols)) return(tibble(domain = character(), indicator = character(), mean_local = numeric(), median_local = numeric(), q05 = numeric(), q95 = numeric()))
  W_dom <- W_cons[, cols, drop = FALSE]        # experts × m_dom
  M_dom <- t(W_dom)                             # items × experts
  prep  <- make_bwm(M_dom)
  fuzzy_raw <- df_group %>% filter(criterion %in% cols) %>% arrange(match(criterion, cols)) %>% pull(m)
  fuzzy_norm <- normalize(fuzzy_raw)
  data_d <- c(prep, list(m_fuzzy = fuzzy_norm, c = c_prior))
  fit_d <- rstan::stan(
    file    = config$stan_file,
    data    = data_d,
    chains  = config$chains,
    iter    = config$iter,
    warmup  = config$warmup,
    control = list(adapt_delta = 0.90, max_treedepth = config$max_treedepth),
    seed    = config$seed + match(dom, domains)
  )
  fits_local[[dom]] <<- fit_d
  indicators_local[[dom]] <<- cols
  Adom <- as.matrix(fit_d, pars = "alpha")
  tibble(
    domain       = dom,
    indicator    = cols,
    mean_local   = colMeans(Adom),
    median_local = apply(Adom, 2, median),
    q05          = apply(Adom, 2, quantile, 0.05),
    q95          = apply(Adom, 2, quantile, 0.95)
  )
})

# 10) Compose global indicator weights
weights_global <- alpha_locals_summary %>%
  left_join(alpha_global_summary %>% select(domain, mean_main = mean, median_main = median), by = "domain") %>%
  mutate(
    global_mean   = mean_local   * mean_main,
    global_median = median_local * median_main
  ) %>%
  arrange(desc(global_mean))

# 11) Save artifacts for later stages (credal, sensitivity, VIKOR)
saveRDS(list(
  config             = config,
  domains            = domains,
  crit_dom           = crit_dom,
  sel_vec            = sel_vec,
  fit_global         = fit_global,
  fits_local         = fits_local,
  indicators_local   = indicators_local,
  alpha_global_sum   = alpha_global_summary,
  alpha_locals_sum   = alpha_locals_summary,
  weights_global_sum = weights_global
), file = config$out_rds)

# 12) Quick console snapshots
print(alpha_global_summary)
print(weights_global, n = Inf, width = Inf)
