# fuzzy_delphi.R

# 0) Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# 1) Config
config <- list(
  input_ai_csv      = "ai.csv",         # id, c1_v1 … c5_v31
  output_fuzzy_csv  = "fuzzy.csv",      # long TFNs for inspection / archival
  max_scale         = 4,                 # Likert upper bound (1–4)
  selection_rule    = "conservative"     # default selection to carry forward
)

# 2) Helpers
gmean <- function(x) exp(mean(log(x)))
normalize <- function(x) x / sum(x)

# 3) Load and orient raw data (reverse C4 at intake)
df_bwm_raw <- readr::read_csv(config$input_ai_csv, show_col_types = FALSE)
stopifnot("id" %in% names(df_bwm_raw))

criteria <- setdiff(names(df_bwm_raw), "id")
meta <- tibble(criterion = criteria) %>%
  mutate(
    domain    = str_extract(criterion, "^(c[1-5])"),
    reverse   = domain == "c4",           # reverse AI-as-threat (C4)
    direction = "max"                      # placeholder for VIKOR; update later per indicator if needed
  ) %>%
  arrange(criterion)

# apply reverse transform: 1↔4 on C4
rev_vars <- meta %>% filter(reverse) %>% pull(criterion)
if (length(rev_vars)) {
  df_bwm <- df_bwm_raw %>%
    mutate(across(all_of(rev_vars), ~ (config$max_scale + 1) - .x))
} else {
  df_bwm <- df_bwm_raw
}

# 4) Numeric matrix W (experts × indicators)
W <- df_bwm %>% select(all_of(c("id", meta$criterion))) %>% select(-id) %>% as.matrix()
stopifnot(ncol(W) == nrow(meta))
colnames(W) <- meta$criterion

# 5) Long TFNs (degenerate per expert)
df_fdm <- df_bwm %>%
  pivot_longer(
    cols      = all_of(meta$criterion),
    names_to  = "criterion",
    values_to = "score"
  ) %>%
  rename(expert = id) %>%
  mutate(l = score, m = score, u = score) %>%
  select(expert, criterion, l, m, u)

readr::write_csv(df_fdm, config$output_fuzzy_csv)

# 6) Panel TFN by indicator (min, gmean, max)
df_group <- df_fdm %>%
  group_by(criterion) %>%
  summarize(
    l = min(l),
    m = gmean(m),
    u = max(u),
    .groups = "drop"
  ) %>%
  arrange(criterion)

# 7) Defuzzify and thresholds (Hashemi Petrudi et al.; geometric aggregation)
df_defuzz <- df_group %>% mutate(a = (l * m * u)^(1/3))

thr_strict        <- 0.8 * config$max_scale
thr_mean          <- mean(df_defuzz$a)
thr_moderate      <- (thr_strict + thr_mean) / 2
thr_conservative  <- 0.8 * thr_mean

selected_strict       <- df_defuzz %>% filter(a >= thr_strict) %>% arrange(desc(a))
selected_mean         <- df_defuzz %>% filter(a >= thr_mean) %>% arrange(desc(a))
selected_moderate     <- df_defuzz %>% filter(a >= thr_moderate) %>% arrange(desc(a))
selected_conservative <- df_defuzz %>% filter(a >= thr_conservative) %>% arrange(desc(a))

selection_summary <- tibble(
  approach   = c("strict", "mean", "moderate", "conservative"),
  threshold  = c(thr_strict, thr_mean, thr_moderate, thr_conservative),
  n_selected = c(nrow(selected_strict), nrow(selected_mean), nrow(selected_moderate), nrow(selected_conservative)),
  criteria   = list(selected_strict$criterion, selected_mean$criterion, selected_moderate$criterion, selected_conservative$criterion)
)

# 8) Carry forward the configured selection
sel_tbl <- switch(
  config$selection_rule,
  strict        = selected_strict,
  mean          = selected_mean,
  moderate      = selected_moderate,
  conservative  = selected_conservative,
  selected_conservative
)
sel_vec <- sel_tbl$criterion

# 9) Save stage artifacts for reuse
saveRDS(list(
  config      = config,
  meta        = meta,
  W           = W,
  df_fdm      = df_fdm,
  df_group    = df_group,
  df_defuzz   = df_defuzz,
  thresholds  = list(strict = thr_strict, mean = thr_mean, moderate = thr_moderate, conservative = thr_conservative),
  selections  = selection_summary,
  sel_tbl     = sel_tbl,
  sel_vec     = sel_vec
), file = "delphi_stage.rds")

# 10) Minimal console outputs for quick checks
print(selection_summary)
print(sel_tbl %>% select(criterion, a) %>% head(10))
