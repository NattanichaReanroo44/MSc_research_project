# ============================================================
# SENSITIVITY: MAINTENANCE EXPENDITURE INSTEAD OF DISASTER EXPENDITURE
# ============================================================

# Standardize maintain_expen, keeping the same variable name
bayes_aft_data_maintain <- bayes_aft_data %>%
  mutate(maintain_expen = as.numeric(scale(maintain_expen)))

# ============================================================
# FORMULAS
# ============================================================

# M8: with magnitude and depth
formula_m8_maintain <- bf(
  recovery_y | cens(censor_type, recovery_y2) ~
    magnitude + depth + ntl_decline_ratio +
    gdp2015 + population + maintain_expen +
    (1 | earthquake_id)
)

# M12: without magnitude and depth
formula_m12_maintain <- bf(
  recovery_y | cens(censor_type, recovery_y2) ~
    ntl_decline_ratio + log_gdp + population + maintain_expen +
    (1 | earthquake_id)
)

# ============================================================
# FIT BOTH WEIBULL RANDOM-EFFECT MODELS
# ============================================================

m8_maintain_result <- fit_bayes_aft(
  formula = formula_m8_maintain,
  model_name = "M8_maintain",
  family_name = "weibull",
  random_effect = TRUE,
  data = bayes_aft_data_maintain
)

m12_maintain_result <- fit_bayes_aft(
  formula = formula_m12_maintain,
  model_name = "M12_maintain",
  family_name = "weibull",
  random_effect = TRUE,
  data = bayes_aft_data_maintain
)

maintain_results <- list(
  m8_maintain_result,
  m12_maintain_result
)

# ============================================================
# FIT STATUS
# ============================================================

maintain_fit_status <- map_dfr(
  maintain_results,
  ~ tibble(
    Model = .x$Model,
    Family = .x$Family,
    Structure = .x$Structure,
    Successful = .x$Successful,
    Error = .x$Error
  )
)

print(maintain_fit_status, n = Inf, width = Inf)

# ============================================================
# CONVERGENCE
# ============================================================

maintain_convergence <- map_dfr(
  maintain_results,
  extract_convergence
)

print(maintain_convergence, n = Inf, width = Inf)

# ============================================================
# ROBUST LOO / ELPD
# ============================================================

maintain_loo_results <- map(
  maintain_results,
  function(x) {
    if (!x$Successful) return(NULL)
    
    compute_loo_robust(
      fit = x$Fit,
      model_name = x$Model,
      family_name = x$Family,
      structure = x$Structure
    )
  }
)

maintain_loo_results <- compact(maintain_loo_results)

maintain_loo_summary <- map_dfr(
  maintain_loo_results,
  ~ tibble(
    Model = .x$Model,
    Family = .x$Family,
    Structure = .x$Structure,
    LOO_method = .x$Method,
    ELPD_LOO = .x$ELPD_LOO,
    SE_ELPD = .x$SE_ELPD,
    P_LOO = .x$P_LOO,
    LOOIC = .x$LOOIC,
    Max_Pareto_k = .x$Max_Pareto_k,
    N_k_over_0_7 = .x$N_k_over_0_7,
    N_k_over_1 = .x$N_k_over_1
  )
) %>%
  arrange(desc(ELPD_LOO))

print(maintain_loo_summary, n = Inf, width = Inf)

# ============================================================
# DIRECT PAIRED LOO COMPARISON
# ============================================================

maintain_loo_objects <- setNames(
  map(maintain_loo_results, "LOO"),
  map_chr(maintain_loo_results, "Model")
)

maintain_loo_compare <- loo::loo_compare(maintain_loo_objects)

maintain_loo_compare

maintain_loo_compare_table <- as.data.frame(maintain_loo_compare) %>%
  rownames_to_column("Model") %>%
  as_tibble() %>%
  rename(
    Delta_ELPD = elpd_diff,
    SE_Delta_ELPD = se_diff
  )

print(maintain_loo_compare_table, n = Inf, width = Inf)

# ============================================================
# COEFFICIENTS
# ============================================================

maintain_coefficients <- map_dfr(
  maintain_results,
  function(x) {
    if (!x$Successful) return(NULL)
    
    posterior_summary <- as.data.frame(
      posterior_summary(x$Fit, pars = "^b_")
    )
    
    posterior_summary %>%
      rownames_to_column("Variable") %>%
      as_tibble() %>%
      mutate(
        Model = x$Model,
        Variable = sub("^b_", "", Variable),
        Time_Ratio = exp(Estimate)
      ) %>%
      select(
        Model,
        Variable,
        Estimate,
        Est.Error,
        Q2.5,
        Q97.5,
        Time_Ratio
      )
  }
)

print(maintain_coefficients, n = Inf, width = Inf)

# ============================================================
# MAINTAIN_EXPEN EFFECT ONLY
# ============================================================

maintain_effect <- maintain_coefficients %>%
  filter(Variable == "maintain_expen") %>%
  mutate(
    Lower_Time_Ratio = exp(Q2.5),
    Upper_Time_Ratio = exp(Q97.5)
  ) %>%
  select(
    Model,
    Estimate,
    Est.Error,
    Q2.5,
    Q97.5,
    Time_Ratio,
    Lower_Time_Ratio,
    Upper_Time_Ratio
  )

print(maintain_effect, n = Inf, width = Inf)


# ============================================================
# POSTERIOR-PREDICTIVE MAE/RMSE
# ============================================================

maintain_mae_rmse <- map_dfr(
  maintain_results,
  ~ calculate_bayes_mae_rmse(
    x = .x,
    data = bayes_aft_data_maintain
  )
)

print(maintain_mae_rmse, n = Inf, width = Inf)

# ============================================================
# JOIN MAE/RMSE WITH LOO/ELPD
# ============================================================

maintain_performance <- maintain_loo_summary %>%
  left_join(
    maintain_mae_rmse %>%
      select(
        Model,
        Family,
        Structure,
        MAE,
        RMSE,
        Median_Error
      ),
    by = c(
      "Model",
      "Family",
      "Structure"
    )
  ) %>%
  arrange(desc(ELPD_LOO))

print(maintain_performance, n = Inf, width = Inf)

# ============================================================
# CLEAN REPORTING TABLE
# ============================================================

maintain_performance_table <- maintain_performance %>%
  transmute(
    Model,
    Family,
    Structure,
    `ELPD-LOO` = round(ELPD_LOO, 2),
    `SE ELPD` = round(SE_ELPD, 2),
    MAE = round(MAE, 2),
    RMSE = round(RMSE, 2),
    `Median error` = round(Median_Error, 2),
    LOOIC = round(LOOIC, 2),
    `LOO method` = LOO_method,
    `Max Pareto k` = round(Max_Pareto_k, 2)
  )

print(maintain_performance_table, n = Inf, width = Inf)

# ============================================================
# BEST BY ELPD
# ============================================================

best_maintain_elpd <- maintain_performance %>%
  filter(is.finite(ELPD_LOO)) %>%
  slice_max(
    ELPD_LOO,
    n = 1,
    with_ties = FALSE
  )

best_maintain_elpd

# ============================================================
# BEST BY MAE
# ============================================================

best_maintain_mae <- maintain_performance %>%
  filter(is.finite(MAE)) %>%
  slice_min(
    MAE,
    n = 1,
    with_ties = FALSE
  )

best_maintain_mae

# ============================================================
# BEST BY RMSE
# ============================================================

best_maintain_rmse <- maintain_performance %>%
  filter(is.finite(RMSE)) %>%
  slice_min(
    RMSE,
    n = 1,
    with_ties = FALSE
  )

best_maintain_rmse


# ============================================================
# SAVE SENSITIVITY RESULTS
# ============================================================

bayesian_maintain_sensitivity_results <- list(
  data = bayes_aft_data_maintain,
  m9_result = m9_maintain_result,
  m12_result = m12_maintain_result,
  fit_status = maintain_fit_status,
  convergence = maintain_convergence,
  loo_results = maintain_loo_results,
  loo_summary = maintain_loo_summary,
  loo_comparison = maintain_loo_compare_table,
  coefficients = maintain_coefficients,
  maintain_effect = maintain_effect
)

saveRDS(
  bayesian_maintain_sensitivity_results,
  file = "bayesian_maintain_expenditure_sensitivity.rds"
)
