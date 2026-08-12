# ============================================================
# PACKAGES
# ============================================================

install.packages('cmprsk')

library(survival)
library(cmprsk)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(knitr)
library(kableExtra)


# ============================================================
# TRUE SIMULATION PARAMETERS
# ============================================================

# Data: all_simulation_scenarios
true_parameters <- c(X1 = 0.50, X2 = 0.35)


# ============================================================
# EXTRACT AFT COEFFICIENT METRICS
# ============================================================

extract_aft_metrics <- function(
    fit,
    model_name,
    simulation_id,
    scenario,
    true_beta = true_parameters) {
  
  model_summary <- summary(fit)$table
  coefficient_names <- intersect(names(true_beta), rownames(model_summary))
  
  if (length(coefficient_names) == 0) {
    stop("No matching coefficients were found.")
  }
  
  estimate <- model_summary[coefficient_names, "Value"]
  standard_error <- model_summary[coefficient_names, "Std. Error"]
  lower_95 <- estimate - 1.96 * standard_error
  upper_95 <- estimate + 1.96 * standard_error
  
  tibble(
    simulation_id = simulation_id,
    scenario = scenario,
    model = model_name,
    parameter = coefficient_names,
    true_value = true_beta[coefficient_names],
    estimate = as.numeric(estimate),
    standard_error = as.numeric(standard_error),
    lower_95 = as.numeric(lower_95),
    upper_95 = as.numeric(upper_95),
    bias = estimate - true_beta[coefficient_names],
    squared_error = (estimate - true_beta[coefficient_names])^2,
    covered =
      lower_95 <= true_beta[coefficient_names] &
      upper_95 >= true_beta[coefficient_names]
  )
}


# ============================================================
# COX MEDIAN PREDICTION
# ============================================================

predict_cox_median <- function(fit, newdata) {
  vapply(
    seq_len(nrow(newdata)),
    function(i) {
      survival_curve <- survfit(
        fit,
        newdata = newdata[i, , drop = FALSE]
      )
      
      survival_table <- summary(survival_curve)$table
      
      predicted_median <- if (is.matrix(survival_table)) {
        survival_table[1, "median"]
      } else {
        survival_table["median"]
      }
      
      predicted_median <- as.numeric(predicted_median)
      
      if (
        length(predicted_median) == 0 ||
        !is.finite(predicted_median)
      ) {
        return(NA_real_)
      }
      
      predicted_median
    },
    numeric(1)
  )
}


# ============================================================
# PREDICTION METRICS
# ============================================================

calculate_prediction_metrics <- function(
    observed_data,
    predicted_time,
    risk_score,
    model_name,
    simulation_id,
    scenario) {
  
  true_target_time <- exp(observed_data$mu)
  
  valid_prediction <-
    is.finite(true_target_time) &
    is.finite(predicted_time)
  
  observed_c_index <- tryCatch(
    {
      survival::concordance(
        Surv(
          observed_data$observed_time,
          observed_data$event
        ) ~ risk_score,
        reverse = TRUE
      )$concordance
    },
    error = function(e) NA_real_
  )
  
  true_c_index <- tryCatch(
    {
      survival::concordance(
        Surv(
          observed_data$true_recovery_time,
          rep(1L, nrow(observed_data))
        ) ~ risk_score,
        reverse = TRUE
      )$concordance
    },
    error = function(e) NA_real_
  )
  
  if (!any(valid_prediction)) {
    return(
      tibble(
        simulation_id = simulation_id,
        scenario = scenario,
        model = model_name,
        n_prediction = 0L,
        median_mae = NA_real_,
        median_rmse = NA_real_,
        observed_c_index = observed_c_index,
        true_c_index = true_c_index
      )
    )
  }
  
  prediction_error <-
    predicted_time[valid_prediction] -
    true_target_time[valid_prediction]
  
  tibble(
    simulation_id = simulation_id,
    scenario = scenario,
    model = model_name,
    n_prediction = sum(valid_prediction),
    median_mae = mean(abs(prediction_error)),
    median_rmse = sqrt(mean(prediction_error^2)),
    observed_c_index = observed_c_index,
    true_c_index = true_c_index
  )
}


# ============================================================
# FIT STANDARD SURVIVAL MODELS
# ============================================================

fit_standard_survival_models <- function(
    simulation_data,
    simulation_id,
    scenario,
    true_beta = true_parameters) {
  
  analysis_data <- simulation_data %>%
    mutate(event = as.integer(event))
  
  # Cox proportional hazards model
  cox_fit <- coxph(
    Surv(observed_time, event) ~ X1 + X2,
    data = analysis_data,
    x = TRUE
  )
  
  cox_predicted_median <- predict_cox_median(
    fit = cox_fit,
    newdata = analysis_data
  )
  
  cox_risk_score <- predict(
    cox_fit,
    newdata = analysis_data,
    type = "lp"
  )
  
  cox_prediction_metrics <- calculate_prediction_metrics(
    observed_data = analysis_data,
    predicted_time = cox_predicted_median,
    risk_score = cox_risk_score,
    model_name = "Cox PH",
    simulation_id = simulation_id,
    scenario = scenario
  )
  
  # Lognormal AFT model
  lognormal_fit <- survreg(
    Surv(observed_time, event) ~ X1 + X2,
    data = analysis_data,
    dist = "lognormal"
  )
  
  lognormal_predicted_time <- predict(
    lognormal_fit,
    newdata = analysis_data,
    type = "quantile",
    p = 0.5
  )
  
  # Larger risk score = shorter survival
  lognormal_risk_score <- -predict(
    lognormal_fit,
    newdata = analysis_data,
    type = "lp"
  )
  
  lognormal_prediction_metrics <- calculate_prediction_metrics(
    observed_data = analysis_data,
    predicted_time = lognormal_predicted_time,
    risk_score = lognormal_risk_score,
    model_name = "Lognormal AFT",
    simulation_id = simulation_id,
    scenario = scenario
  )
  
  lognormal_coefficient_metrics <- extract_aft_metrics(
    fit = lognormal_fit,
    model_name = "Lognormal AFT",
    simulation_id = simulation_id,
    scenario = scenario,
    true_beta = true_beta
  )
  
  # Weibull AFT model
  weibull_fit <- survreg(
    Surv(observed_time, event) ~ X1 + X2,
    data = analysis_data,
    dist = "weibull"
  )
  
  weibull_predicted_time <- predict(
    weibull_fit,
    newdata = analysis_data,
    type = "quantile",
    p = 0.5
  )
  
  weibull_risk_score <- -predict(
    weibull_fit,
    newdata = analysis_data,
    type = "lp"
  )
  
  weibull_prediction_metrics <- calculate_prediction_metrics(
    observed_data = analysis_data,
    predicted_time = weibull_predicted_time,
    risk_score = weibull_risk_score,
    model_name = "Weibull AFT",
    simulation_id = simulation_id,
    scenario = scenario
  )
  
  weibull_coefficient_metrics <- extract_aft_metrics(
    fit = weibull_fit,
    model_name = "Weibull AFT",
    simulation_id = simulation_id,
    scenario = scenario,
    true_beta = true_beta
  )
  
  list(
    coefficient_metrics = bind_rows(
      lognormal_coefficient_metrics,
      weibull_coefficient_metrics
    ),
    prediction_metrics = bind_rows(
      cox_prediction_metrics,
      lognormal_prediction_metrics,
      weibull_prediction_metrics
    ),
    fits = list(
      cox = cox_fit,
      lognormal = lognormal_fit,
      weibull = weibull_fit
    )
  )
}


# ============================================================
# SAFE MODEL FITTING
# ============================================================

fit_standard_models_safely <- function(
    simulation_data,
    simulation_id,
    scenario) {
  
  tryCatch(
    {
      result <- fit_standard_survival_models(
        simulation_data = simulation_data,
        simulation_id = simulation_id,
        scenario = scenario
      )
      
      list(
        successful = TRUE,
        error_message = NA_character_,
        coefficient_metrics = result$coefficient_metrics,
        prediction_metrics = result$prediction_metrics
      )
    },
    error = function(e) {
      list(
        successful = FALSE,
        error_message = conditionMessage(e),
        coefficient_metrics = NULL,
        prediction_metrics = NULL
      )
    }
  )
}


# ============================================================
# RUN STANDARD MODELS ACROSS SIMULATIONS
# ============================================================

run_standard_simulation_list <- function(
    dataset_list,
    scenario_name) {
  
  results <- map2(
    dataset_list,
    seq_along(dataset_list),
    ~ fit_standard_models_safely(
      simulation_data = .x,
      simulation_id = .y,
      scenario = scenario_name
    )
  )
  
  coefficient_results <- map_dfr(
    results,
    "coefficient_metrics"
  )
  
  prediction_results <- map_dfr(
    results,
    "prediction_metrics"
  )
  
  failure_results <- imap_dfr(
    results,
    function(result, simulation_id) {
      tibble(
        simulation_id = simulation_id,
        scenario = scenario_name,
        successful = result$successful,
        error_message = result$error_message
      )
    }
  )
  
  list(
    coefficient_results = coefficient_results,
    prediction_results = prediction_results,
    failures = failure_results
  )
}


# ============================================================
# SCENARIO 1: NON-INFORMATIVE CENSORING
# ============================================================

scenario1_results_10 <- run_standard_simulation_list(
  dataset_list = scenario1_10,
  scenario_name = "Non-informative 10%"
)

scenario1_results_25 <- run_standard_simulation_list(
  dataset_list = scenario1_25,
  scenario_name = "Non-informative 25%"
)

scenario1_results_40 <- run_standard_simulation_list(
  dataset_list = scenario1_40,
  scenario_name = "Non-informative 40%"
)


# ============================================================
# SCENARIO 2: INFORMATIVE CENSORING
# ============================================================

scenario2_results_weak <- run_standard_simulation_list(
  dataset_list = scenario2_informative$Weak,
  scenario_name = "Informative weak"
)

scenario2_results_moderate <- run_standard_simulation_list(
  dataset_list = scenario2_informative$Moderate,
  scenario_name = "Informative moderate"
)

scenario2_results_strong <- run_standard_simulation_list(
  dataset_list = scenario2_informative$Strong,
  scenario_name = "Informative strong"
)


# ============================================================
# FINE-GRAY MODEL
# ============================================================

fit_fine_gray_model <- function(
    simulation_data,
    simulation_id,
    scenario) {
  
  covariate_matrix <- model.matrix(
    ~ X1 + X2,
    data = simulation_data
  )[, -1, drop = FALSE]
  
  fine_gray_fit <- cmprsk::crr(
    ftime = simulation_data$observed_time,
    fstatus = simulation_data$event_type,
    cov1 = covariate_matrix,
    failcode = 1,
    cencode = 0
  )
  
  estimates <- fine_gray_fit$coef
  standard_errors <- sqrt(diag(fine_gray_fit$var))
  
  tibble(
    simulation_id = simulation_id,
    scenario = scenario,
    model = "Fine-Gray",
    parameter = names(estimates),
    estimate = as.numeric(estimates),
    standard_error = as.numeric(standard_errors),
    lower_95 = estimates - 1.96 * standard_errors,
    upper_95 = estimates + 1.96 * standard_errors,
    subdistribution_hazard_ratio = exp(estimates)
  )
}


# ============================================================
# RUN FINE-GRAY MODELS
# ============================================================

run_fine_gray_list <- function(dataset_list, scenario_name) {
  map2_dfr(
    dataset_list,
    seq_along(dataset_list),
    function(simulation_data, simulation_id) {
      
      tryCatch(
        fit_fine_gray_model(
          simulation_data = simulation_data,
          simulation_id = simulation_id,
          scenario = scenario_name
        ),
        error = function(e) {
          tibble(
            simulation_id = simulation_id,
            scenario = scenario_name,
            model = "Fine-Gray",
            parameter = NA_character_,
            estimate = NA_real_,
            standard_error = NA_real_,
            lower_95 = NA_real_,
            upper_95 = NA_real_,
            subdistribution_hazard_ratio = NA_real_
          )
        }
      )
    }
  )
}

fine_gray_10 <- run_fine_gray_list(
  scenario3_competing$competing_10,
  "Competing risk 10%"
)

fine_gray_25 <- run_fine_gray_list(
  scenario3_competing$competing_25,
  "Competing risk 25%"
)

fine_gray_40 <- run_fine_gray_list(
  scenario3_competing$competing_40,
  "Competing risk 40%"
)


# ============================================================
# SCENARIO 3: COMPETING RISKS AS CENSORING
# ============================================================

prepare_competing_as_censoring <- function(data) {
  data %>%
    mutate(
      event = as.integer(event_type == 1)
    )
}

scenario3_standard_10 <- run_standard_simulation_list(
  dataset_list = map(
    scenario3_competing$competing_10,
    prepare_competing_as_censoring
  ),
  scenario_name = "Competing risk 10%"
)

scenario3_standard_25 <- run_standard_simulation_list(
  dataset_list = map(
    scenario3_competing$competing_25,
    prepare_competing_as_censoring
  ),
  scenario_name = "Competing risk 25%"
)

scenario3_standard_40 <- run_standard_simulation_list(
  dataset_list = map(
    scenario3_competing$competing_40,
    prepare_competing_as_censoring
  ),
  scenario_name = "Competing risk 40%"
)


# ============================================================
# COMBINE COEFFICIENT RESULTS
# ============================================================

all_coefficient_results <- bind_rows(
  scenario1_results_10$coefficient_results,
  scenario1_results_25$coefficient_results,
  scenario1_results_40$coefficient_results,
  
  scenario2_results_weak$coefficient_results,
  scenario2_results_moderate$coefficient_results,
  scenario2_results_strong$coefficient_results,
  
  scenario3_standard_10$coefficient_results,
  scenario3_standard_25$coefficient_results,
  scenario3_standard_40$coefficient_results
)


# ============================================================
# COEFFICIENT PERFORMANCE SUMMARY
# ============================================================

coefficient_performance_summary <-
  all_coefficient_results %>%
  group_by(
    scenario,
    model,
    parameter
  ) %>%
  summarise(
    n_successful = sum(is.finite(estimate)),
    mean_estimate = mean(estimate, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    estimation_rmse = sqrt(
      mean(squared_error, na.rm = TRUE)
    ),
    coverage_95 = mean(covered, na.rm = TRUE),
    mean_standard_error = mean(
      standard_error,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

coefficient_performance_summary


# ============================================================
# COMBINE PREDICTION RESULTS
# ============================================================

all_prediction_results <- bind_rows(
  scenario1_results_10$prediction_results,
  scenario1_results_25$prediction_results,
  scenario1_results_40$prediction_results,
  
  scenario2_results_weak$prediction_results,
  scenario2_results_moderate$prediction_results,
  scenario2_results_strong$prediction_results,
  
  scenario3_standard_10$prediction_results,
  scenario3_standard_25$prediction_results,
  scenario3_standard_40$prediction_results
)


# ============================================================
# PREDICTION PERFORMANCE SUMMARY
# ============================================================

prediction_performance_summary <-
  all_prediction_results %>%
  group_by(
    scenario,
    model
  ) %>%
  summarise(
    n_successful = sum(is.finite(median_mae)),
    mean_mae = mean(median_mae, na.rm = TRUE),
    sd_mae = sd(median_mae, na.rm = TRUE),
    mean_prediction_rmse = mean(
      median_rmse,
      na.rm = TRUE
    ),
    sd_prediction_rmse = sd(
      median_rmse,
      na.rm = TRUE
    ),
    mean_observed_c_index = mean(
      observed_c_index,
      na.rm = TRUE
    ),
    sd_observed_c_index = sd(
      observed_c_index,
      na.rm = TRUE
    ),
    mean_true_c_index = mean(
      true_c_index,
      na.rm = TRUE
    ),
    sd_true_c_index = sd(
      true_c_index,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

prediction_performance_summary


# ============================================================
# FINE-GRAY PERFORMANCE SUMMARY
# ============================================================

all_fine_gray_results <- bind_rows(
  fine_gray_10,
  fine_gray_25,
  fine_gray_40
)

fine_gray_summary <- all_fine_gray_results %>%
  filter(!is.na(parameter)) %>%
  group_by(
    scenario,
    parameter
  ) %>%
  summarise(
    n_successful = sum(is.finite(estimate)),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    mean_subdistribution_hazard_ratio = mean(
      subdistribution_hazard_ratio,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

fine_gray_summary


# ============================================================
# SAVE SIMULATION RESULTS
# ============================================================

simulation_model_results <- list(
  coefficient_results = all_coefficient_results,
  coefficient_summary = coefficient_performance_summary,
  prediction_results = all_prediction_results,
  prediction_summary = prediction_performance_summary,
  fine_gray_results = all_fine_gray_results,
  fine_gray_summary = fine_gray_summary
)

saveRDS(
  simulation_model_results,
  file = "simulation_model_results.rds"
)


# ============================================================
# COEFFICIENT RESULTS TABLE
# ============================================================

coefficient_main_table <- coefficient_performance_summary %>%
  select(
    scenario,
    model,
    parameter,
    mean_bias,
    estimation_rmse,
    coverage_95
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = c(
      mean_bias,
      estimation_rmse,
      coverage_95
    )
  ) %>%
  mutate(
    across(
      starts_with("mean_bias"),
      ~ round(.x, 3)
    ),
    across(
      starts_with("estimation_rmse"),
      ~ round(.x, 3)
    ),
    across(
      starts_with("coverage_95"),
      ~ round(100 * .x, 1)
    )
  ) %>%
  rename(
    Scenario = scenario,
    Model = model,
    `Bias X1` = mean_bias_X1,
    `Bias X2` = mean_bias_X2,
    `RMSE X1` = estimation_rmse_X1,
    `RMSE X2` = estimation_rmse_X2,
    `Coverage X1 (%)` = coverage_95_X1,
    `Coverage X2 (%)` = coverage_95_X2
  )


# ============================================================
# PREDICTION RESULTS TABLE
# ============================================================

predictive_main_table <- prediction_performance_summary %>%
  mutate(
    `Median MAE, mean (SD)` = sprintf(
      "%.2f (%.2f)",
      mean_mae,
      sd_mae
    ),
    `Median RMSE, mean (SD)` = sprintf(
      "%.2f (%.2f)",
      mean_prediction_rmse,
      sd_prediction_rmse
    ),
    `Observed C-index, mean (SD)` = sprintf(
      "%.3f (%.3f)",
      mean_observed_c_index,
      sd_observed_c_index
    ),
    `True C-index, mean (SD)` = sprintf(
      "%.3f (%.3f)",
      mean_true_c_index,
      sd_true_c_index
    )
  ) %>%
  select(
    Scenario = scenario,
    Model = model,
    `Median MAE, mean (SD)`,
    `Median RMSE, mean (SD)`,
    `Observed C-index, mean (SD)`,
    `True C-index, mean (SD)`
  )

print(
  predictive_main_table,
  n = Inf,
  width = Inf
)


# ============================================================
# FINE-GRAY RESULTS TABLE
# ============================================================

fine_gray_main_table <- fine_gray_summary %>%
  mutate(
    `Estimate, mean (SD)` = sprintf(
      "%.3f (%.3f)",
      mean_estimate,
      sd_estimate
    ),
    `Mean subdistribution HR` = round(
      mean_subdistribution_hazard_ratio,
      3
    )
  ) %>%
  select(
    Scenario = scenario,
    Parameter = parameter,
    `Estimate, mean (SD)`,
    `Mean subdistribution HR`
  )


# ============================================================
# FINAL OUTPUTS
# ============================================================

coefficient_performance_summary
coefficient_main_table
print(predictive_main_table, n = 27)
fine_gray_main_table