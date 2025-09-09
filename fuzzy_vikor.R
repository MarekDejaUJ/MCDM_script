# fuzzy_vikor.R — posterior‑weighted fuzzy VIKOR

suppressPackageStartupMessages({
  library(tidyverse)
  library(FuzzyMCDM)
})

# if (!requireNamespace("FuzzyMCDM", quietly = TRUE)) {
#   install.packages("FuzzyMCDM")
# }

# 0) Config
config <- list(
  delphi_rds = "delphi_stage.rds",      # from fuzzy_delphi.R
  bbwm_rds   = "bbwm_stage.rds",        # from bwm_hier.R
  groups_csv = "group2.csv",            # columns: id; type; level
  groups_delim = ";",
  out_dir    = "results_vikor",
  n_sims     = 1000,
  v          = 0.5,                      # trade-off parameter (0.5 = balance S and R)
  seed       = 20250908
)

set.seed(config$seed)

# 1) Load artifacts
D <- readRDS(config$delphi_rds)
B <- readRDS(config$bbwm_rds)

sel_cons   <- D$sel_vec                 # retained indicators (ordered)
meta       <- D$meta %>% filter(criterion %in% sel_cons) %>% arrange(match(criterion, sel_cons))
df_fdm     <- D$df_fdm                  # expert × criterion TFNs (degenerate)

domains          <- B$domains
fit_global       <- B$fit_global
fits_local       <- B$fits_local
indicators_local <- B$indicators_local

# 2) Grouping info for alternatives
if (file.exists(config$groups_csv)) {
  df_group_raw <- readr::read_delim(config$groups_csv, delim = config$groups_delim, show_col_types = FALSE)
} else {
  stop("groups_csv not found: ", config$groups_csv)
}
stopifnot(all(c("id","type") %in% names(df_group_raw)))

alts <- unique(df_group_raw$type)

# 3) Decision fuzzy matrix (rows = alts; columns = 3*|criteria|, ordered by sel_cons)
#    Each cell is a TFN assembled as (min l, gmean m, max u) across experts in that alt

gmean <- function(x) exp(mean(log(x)))

fdm_sel <- df_fdm %>%
  filter(criterion %in% sel_cons) %>%
  left_join(df_group_raw %>% select(id, type), by = c("expert" = "id"))

n_crit <- length(sel_cons)

decision_fuzzy <- matrix(NA_real_, nrow = length(alts), ncol = n_crit * 3,
                         dimnames = list(alts, rep(sel_cons, each = 3)))

for (alt in alts) {
  sub <- fdm_sel %>% filter(type == alt)
  for (j in seq_along(sel_cons)) {
    cri <- sel_cons[j]
    x <- sub %>% filter(criterion == cri)
    decision_fuzzy[alt, (3*j-2):(3*j)] <- c(min(x$l), gmean(x$m), max(x$u))
  }
}

# 4) Cost/benefit vector from metadata (default 'max'; override by editing meta$direction)
cb <- meta %>% pull(direction)
stopifnot(all(cb %in% c("max","min")))

# 5) Posterior draws of weights (compose global × local)
A_global <- as.matrix(fit_global, pars = "alpha")
colnames(A_global) <- domains

A_local <- purrr::map(domains, function(dom) {
  mat <- as.matrix(fits_local[[dom]], pars = "alpha")
  colnames(mat) <- indicators_local[[dom]]
  mat
}) %>% setNames(domains)

w_sims <- matrix(NA_real_, nrow = config$n_sims, ncol = n_crit, dimnames = list(NULL, sel_cons))

for (s in seq_len(config$n_sims)) {
  ig <- sample.int(nrow(A_global), 1)
  for (dom in domains) {
    if (!dom %in% names(A_local)) next
    il   <- sample.int(nrow(A_local[[dom]]), 1)
    inds <- indicators_local[[dom]]
    keep <- inds[inds %in% sel_cons]
    if (!length(keep)) next
    w_sims[s, match(keep, sel_cons)] <- A_global[ig, dom] * A_local[[dom]][il, keep]
  }
  w_sims[s, ] <- w_sims[s, ] / sum(w_sims[s, ])
}

# 6) Run fuzzy VIKOR across posterior weight draws
vikor_runs <- vector("list", config$n_sims)
for (s in seq_len(config$n_sims)) {
  w_tfn <- rep(w_sims[s, ], each = 3)
  res <- FuzzyVIKOR(decision_fuzzy, w_tfn, cb, config$v)
  res$sim <- s
  res$alt <- res$Alternatives
  vikor_runs[[s]] <- res %>% select(sim, alt, Def_S, Def_R, Def_Q, Ranking)
}

vikor_posterior <- bind_rows(vikor_runs)

# 7) Summaries and decision diagnostics
vikor_summary <- vikor_posterior %>%
  group_by(alt) %>%
  summarize(
    mean_S = mean(Def_S), q05_S = quantile(Def_S, 0.05), q50_S = median(Def_S), q95_S = quantile(Def_S, 0.95),
    mean_R = mean(Def_R), q05_R = quantile(Def_R, 0.05), q50_R = median(Def_R), q95_R = quantile(Def_R, 0.95),
    mean_Q = mean(Def_Q), q05_Q = quantile(Def_Q, 0.05), q50_Q = median(Def_Q), q95_Q = quantile(Def_Q, 0.95),
    mean_rank = mean(Ranking),
    p_best = mean(Ranking == 1),
    p_top2 = mean(Ranking <= 2),
    p_top3 = mean(Ranking <= 3),
    .groups = "drop"
  ) %>%
  arrange(mean_Q)

# Acceptable-advantage (C1) and acceptable-stability (C2) using median Q, S, R
crisp <- vikor_summary %>% arrange(q50_Q)
J <- nrow(crisp)
DQ <- 1 / (J - 1)
Adv <- (crisp$q50_Q[2] - crisp$q50_Q[1]) / (crisp$q50_Q[J] - crisp$q50_Q[1])
C1 <- Adv > DQ
best_alt <- crisp$alt[1]
C2 <- (crisp$q50_S[crisp$alt == best_alt] == min(crisp$q50_S)) || (crisp$q50_R[crisp$alt == best_alt] == min(crisp$q50_R))
vikor_decision <- tibble(
  best_alt = best_alt,
  DQ = DQ,
  Adv = Adv,
  C1 = C1,
  C2 = C2
)

# 8) Save artifacts
if (!dir.exists(config$out_dir)) dir.create(config$out_dir, recursive = TRUE)
readr::write_csv(vikor_summary,   file.path(config$out_dir, "vikor_summary.csv"))
readr::write_csv(vikor_posterior, file.path(config$out_dir, "vikor_posterior.csv"))
readr::write_csv(vikor_decision,  file.path(config$out_dir, "vikor_decision.csv"))

saveRDS(list(
  config          = config,
  sel_cons        = sel_cons,
  alts            = alts,
  decision_fuzzy  = decision_fuzzy,
  cb              = cb,
  w_sims          = w_sims,
  vikor_posterior = vikor_posterior,
  vikor_summary   = vikor_summary,
  vikor_decision  = vikor_decision
), file = file.path(config$out_dir, "vikor_stage.rds"))

# 9) Console snapshot
print(vikor_summary, n = Inf, width = Inf)
print(vikor_decision)
