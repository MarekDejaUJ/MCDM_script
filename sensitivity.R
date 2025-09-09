# sensitivity.R — ranked output with baseline + scenarios (weights + ranks)

suppressPackageStartupMessages({
  library(tidyverse)
})

# 0) Config
config <- list(
  bbwm_rds       = "bbwm_stage.rds",
  out_rds        = "sensitivity_stage.rds",
  out_dir        = "results_sensitivity",
  stat_main      = "median",   # "mean" or "median"
  stat_local     = "median",   # "mean" or "median"
  domain_to_vary = "c1",
  grid           = seq(0.1, 0.9, by = 0.1)
)

# 1) Load BBWM artifacts
B <- readRDS(config$bbwm_rds)

domains          <- B$domains
alpha_global_sum <- B$alpha_global_sum  # columns: domain, mean, median, q05, q95
alpha_locals_sum <- B$alpha_locals_sum  # columns: domain, indicator, mean_local, median_local, q05, q95

# 2) Column selectors
get_main_col  <- function(stat)  switch(stat, mean = "mean",         median = "median",         stop("stat_main must be 'mean' or 'median'"))
get_local_col <- function(stat)  switch(stat, mean = "mean_local",   median = "median_local",   stop("stat_local must be 'mean' or 'median'"))

main_col  <- get_main_col(config$stat_main)
local_col <- get_local_col(config$stat_local)
stopifnot(main_col %in% names(alpha_global_sum))
stopifnot(local_col %in% names(alpha_locals_sum))

# 3) Baseline vectors (normalized main weights)
alpha_main <- alpha_global_sum %>%
  arrange(match(domain, domains)) %>%
  transmute(domain, w = .data[[main_col]]) %>%
  { setNames(.$w, .$domain) }
alpha_main <- alpha_main / sum(alpha_main)

local_tbl <- alpha_locals_sum %>%
  transmute(domain, indicator, weight_local = .data[[local_col]])

# 4) Baseline global weights (no scenario tilt)
baseline_tbl <- local_tbl %>%
  mutate(original_weight = as.numeric(alpha_main[domain]) * weight_local) %>%
  arrange(desc(original_weight)) %>%
  mutate(original_rank = row_number()) %>%
  select(indicator, original_weight, original_rank, domain)

# 5) Reallocation helper
reallocate_main <- function(alpha_main, vary_domain, w_new){
  stopifnot(vary_domain %in% names(alpha_main))
  others <- setdiff(names(alpha_main), vary_domain)
  others_sum <- sum(alpha_main[others])
  out <- alpha_main
  out[vary_domain] <- w_new
  out[others] <- alpha_main[others] / others_sum * (1 - w_new)
  out / sum(out)
}

# 6) Build sensitivity scenarios (long)
sensitivity_long <- purrr::map_dfr(config$grid, function(w_new){
  new_main <- reallocate_main(alpha_main, config$domain_to_vary, w_new)
  tibble(
    domain          = names(new_main),
    weight_main_new = as.numeric(new_main),
    scenario        = paste0("S", sprintf("%d", round(w_new * 10)))
  ) %>%
    left_join(local_tbl, by = "domain") %>%
    mutate(weight_global_new = weight_main_new * weight_local) %>%
    select(scenario, domain, indicator, weight_local, weight_main_new, weight_global_new)
})

# 7) Rank within each scenario
ranked <- sensitivity_long %>%
  group_by(scenario) %>%
  arrange(desc(weight_global_new), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 8) Wide tables
# 8a) Main weights per scenario (for reporting)
table_main <- sensitivity_long %>%
  distinct(scenario, domain, weight_main_new) %>%
  arrange(scenario, domain) %>%
  pivot_wider(names_from = domain, values_from = weight_main_new)

# 8b) Indicator-level: include baseline, then scenarios (both weight and rank)
wide_scenarios <- ranked %>%
  select(scenario, indicator, weight_global_new, rank) %>%
  arrange(indicator, scenario) %>%
  pivot_wider(names_from = scenario, values_from = c(weight_global_new, rank))

final_table <- baseline_tbl %>%
  select(indicator, original_weight) %>%
  left_join(wide_scenarios, by = "indicator") %>%
  arrange(desc(original_weight))

# 9) Save artifacts
if (!dir.exists(config$out_dir)) dir.create(config$out_dir, recursive = TRUE)
readr::write_csv(table_main,  file.path(config$out_dir, paste0("table_main_",  config$domain_to_vary, ".csv")))
readr::write_csv(final_table, file.path(config$out_dir, paste0("table_ranked_", config$domain_to_vary, ".csv")))

saveRDS(list(
  config            = config,
  alpha_main        = alpha_main,
  local_tbl         = local_tbl,
  baseline_tbl      = baseline_tbl,
  sensitivity_long  = sensitivity_long,
  ranked            = ranked,
  table_main        = table_main,
  final_table       = final_table
), file = config$out_rds)

# 10) Console snapshots
print(table_main)
print(final_table, n = Inf, width = Inf)
