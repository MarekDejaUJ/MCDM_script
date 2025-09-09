# After the script has been run and the models have been fitted

# 1. Get R-hat for the global model
rhat_global <- summary(fit_global)$summary[, "Rhat"]
print("R-hat for Global Criteria:")
print(rhat_global)

# 2. Get R-hat for each local model
rhat_locals <- lapply(fits_local, function(fit) {
  summary(fit)$summary[, "Rhat"]
})

# Print the R-hat values for each local model
print("R-hat for Local Indicators (by domain):")
print(rhat_locals)
