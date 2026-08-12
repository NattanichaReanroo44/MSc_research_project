# ============================================================
# PACKAGES
# ============================================================

library(survival)
library(flexsurv)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(brms)

# ============================================================
# DATA PREPARATION
# ============================================================

pga_recovery_econ_clean1 = read_csv('/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv')

data_table <- pga_recovery_econ_clean1 %>%
  mutate(
    population = population / 1e6, # million persons
    gdp2015 = gdp2015 / 1e6, # trillion yen
    dis_expen = dis_expen / 1e6, # billion yen
    maintain_expen = maintain_expen / 1e6,
    restore_expen = restore_expen / 1e6,
    log_gdp = log(gdp2015),
    log_depth = log(depth),
    event = ifelse(first_event == "recovered", 1, 0)
  )

aft_cv_data <- data_table

# ============================================================
# CANDIDATE AFT MODELS
# ============================================================

candidate_formulas <- list(
  M1 = Surv(recovery_days, event) ~ magnitude + depth,
  M2 = Surv(recovery_days, event) ~ magnitude + depth + mean_pga,
  M3 = Surv(recovery_days, event) ~ magnitude + depth + max_pga,
  M4 = Surv(recovery_days, event) ~ magnitude + depth + ntl_decline_ratio,
  M5 = Surv(recovery_days, event) ~ magnitude + depth + gdp2015 + population + dis_expen,
  M6 = Surv(recovery_days, event) ~ magnitude + depth + mean_pga + ntl_decline_ratio,
  M7 = Surv(recovery_days, event) ~ magnitude + depth + max_pga + ntl_decline_ratio,
  M8 = Surv(recovery_days, event) ~ magnitude + depth + ntl_decline_ratio + gdp2015 + population + dis_expen,
  M9 = Surv(recovery_days, event) ~ magnitude + depth + mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen,
  M10 = Surv(recovery_days, event) ~ magnitude + depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen,
  M11 = Surv(recovery_days, event) ~ magnitude + depth + log_gdp + population + dis_expen,
  M12 = Surv(recovery_days, event) ~ magnitude + depth + ntl_decline_ratio + log_gdp + population + dis_expen,
  M13 = Surv(recovery_days, event) ~ magnitude + depth + mean_pga + ntl_decline_ratio + log_gdp + population + dis_expen,
  M14 = Surv(recovery_days, event) ~ magnitude + depth + max_pga + ntl_decline_ratio + log_gdp + population + dis_expen
)

# ============================================================
# HELD-OUT LOG PREDICTIVE DENSITY
# ============================================================

heldout_log_density <- function(fit, test_row) {
  evaluation_time <- test_row$recovery_days[[1]]
  event_status <- test_row$event[[1]]
  
  predicted_hazard <- summary(
    fit,
    newdata = test_row,
    type = "hazard",
    t = evaluation_time,
    ci = FALSE,
    cross = FALSE
  )[[1]]$est[[1]]
  
  predicted_cumhaz <- summary(
    fit,
    newdata = test_row,
    type = "cumhaz",
    t = evaluation_time,
    ci = FALSE,
    cross = FALSE
  )[[1]]$est[[1]]
  
  if (
    !is.finite(predicted_cumhaz) ||
    predicted_cumhaz < 0
  ) {
    stop("Invalid predicted cumulative hazard.")
  }
  
  if (event_status == 1) {
    if (
      !is.finite(predicted_hazard) ||
      predicted_hazard <= 0
    ) {
      stop("Invalid predicted hazard for observed event.")
    }
    
    log(predicted_hazard) - predicted_cumhaz
  } else {
    -predicted_cumhaz
  }
}

# ============================================================
# HELD-OUT MEDIAN PREDICTION
# ============================================================

heldout_median_prediction <- function(fit, test_row) {
  prediction <- predict(
    fit,
    newdata = test_row,
    type = "quantile",
    p = 0.5,
    conf.int = FALSE
  )
  
  if (".pred_quantile" %in% names(prediction)) {
    predicted_median <- prediction$.pred_quantile[[1]]
  } else if ("est" %in% names(prediction)) {
    predicted_median <- prediction$est[[1]]
  } else {
    stop(
      paste(
        "Median prediction column not found.",
        "Available columns:",
        paste(names(prediction), collapse = ", ")
      )
    )
  }
  
  if (
    length(predicted_median) != 1 ||
    is.na(predicted_median) ||
    !is.finite(predicted_median) ||
    predicted_median <= 0
  ) {
    stop("Invalid median recovery-time prediction.")
  }
  
  predicted_median
}

# ============================================================
# LOOCV FOR ONE AFT MODEL
# ============================================================

loocv_one_aft_model <- function(
    formula,
    data,
    distribution,
    model_name) {
  
  n <- nrow(data)
  
  observation_results <- purrr::map_dfr(
    seq_len(n),
    function(i) {
      training_data <- data[-i, , drop = FALSE]
      test_data <- data[i, , drop = FALSE]
      
      fitted_model <- tryCatch(
        flexsurv::flexsurvreg(
          formula = formula,
          data = training_data,
          dist = distribution
        ),
        error = function(e) e
      )
      
      # Model fitting failure
      if (inherits(fitted_model, "error")) {
        return(
          tibble::tibble(
            observation = i,
            event = test_data$event[[1]],
            recovery_days = test_data$recovery_days[[1]],
            predicted_median = NA_real_,
            prediction_error = NA_real_,
            squared_error = NA_real_,
            absolute_error = NA_real_,
            log_predictive_density = NA_real_,
            fit_failed = TRUE,
            log_score_failed = TRUE,
            median_prediction_failed = TRUE,
            error_message = conditionMessage(fitted_model)
          )
        )
      }
      
      score_result <- tryCatch(
        {
          list(
            value = heldout_log_density(
              fit = fitted_model,
              test_row = test_data
            ),
            error = NA_character_
          )
        },
        error = function(e) {
          list(
            value = NA_real_,
            error = paste0(
              "Log-score error: ",
              conditionMessage(e)
            )
          )
        }
      )
      
      prediction_result <- tryCatch(
        {
          list(
            value = heldout_median_prediction(
              fit = fitted_model,
              test_row = test_data
            ),
            error = NA_character_
          )
        },
        error = function(e) {
          list(
            value = NA_real_,
            error = paste0(
              "Median-prediction error: ",
              conditionMessage(e)
            )
          )
        }
      )
      
      observed_time <- test_data$recovery_days[[1]]
      event_status <- test_data$event[[1]]
      predicted_median <- prediction_result$value
      
      # Prediction error is calculated only for observed recoveries
      if (
        event_status == 1 &&
        is.finite(observed_time) &&
        is.finite(predicted_median)
      ) {
        prediction_error <- predicted_median - observed_time
        squared_error <- prediction_error^2
        absolute_error <- abs(prediction_error)
      } else {
        prediction_error <- NA_real_
        squared_error <- NA_real_
        absolute_error <- NA_real_
      }
      
      errors <- c(
        score_result$error,
        prediction_result$error
      )
      
      errors <- errors[
        !is.na(errors) &
          nzchar(errors)
      ]
      
      combined_error <- if (length(errors) == 0) {
        NA_character_
      } else {
        paste(errors, collapse = " | ")
      }
      
      tibble::tibble(
        observation = i,
        event = event_status,
        recovery_days = observed_time,
        predicted_median = predicted_median,
        prediction_error = prediction_error,
        squared_error = squared_error,
        absolute_error = absolute_error,
        log_predictive_density = score_result$value,
        fit_failed = FALSE,
        log_score_failed = !is.finite(score_result$value),
        median_prediction_failed = !is.finite(predicted_median),
        error_message = combined_error
      )
    }
  )
  
  observation_results %>%
    dplyr::mutate(
      model = model_name,
      distribution = distribution,
      .before = 1
    )
}

# ============================================================
# RUN LOOCV FOR ALL AFT MODELS
# ============================================================

aft_distributions <- c(
  Lognormal = "lnorm",
  Loglogistic = "llogis",
  Weibull = "weibull"
)

aft_loocv_predictions <- imap_dfr(
  aft_distributions,
  function(distribution_code, distribution_name) {
    imap_dfr(
      candidate_formulas,
      function(model_formula, model_name) {
        message(
          "LOOCV: ",
          distribution_name,
          " ",
          model_name
        )
        
        loocv_one_aft_model(
          formula = model_formula,
          data = aft_cv_data,
          distribution = distribution_code,
          model_name = model_name
        ) %>%
          mutate(
            distribution = distribution_name
          )
      }
    )
  }
)

# ============================================================
# LOOCV CHECKS
# ============================================================

aft_loocv_predictions %>%
  summarise(
    total_rows = n(),
    valid_medians = sum(is.finite(predicted_median)),
    failed_medians = sum(
      median_prediction_failed,
      na.rm = TRUE
    ),
    valid_event_predictions = sum(
      event == 1 &
        is.finite(predicted_median),
      na.rm = TRUE
    )
  )

aft_loocv_predictions[
  aft_loocv_predictions$distribution == "Weibull",
]

table(aft_loocv_predictions$distribution)

weibull_worst <- aft_loocv_predictions %>%
  filter(
    distribution == "Weibull"
  ) %>%
  arrange(log_predictive_density) %>%
  dplyr::select(
    model,
    observation,
    event,
    recovery_days,
    log_predictive_density,
    fit_failed,
    error_message
  )

weibull_worst %>%
  slice_head(n = 30)

saveRDS(
  aft_loocv_predictions,
  "aft_loocv_predictions.rds"
)

# ============================================================
# LOOCV PERFORMANCE SUMMARY
# ============================================================

aft_loocv_comparison <- aft_loocv_predictions %>%
  group_by(
    distribution,
    model
  ) %>%
  summarise(
    n = n(),
    
    failed_log_predictions = sum(
      !is.finite(log_predictive_density)
    ),
    
    failed_median_predictions = sum(
      median_prediction_failed,
      na.rm = TRUE
    ),
    
    n_events = sum(
      event == 1,
      na.rm = TRUE
    ),
    
    n_events_used_rmse = sum(
      event == 1 &
        is.finite(predicted_median) &
        is.finite(recovery_days),
      na.rm = TRUE
    ),
    
    elpd_loo = {
      valid_scores <- log_predictive_density[
        is.finite(log_predictive_density)
      ]
      
      if (length(valid_scores) == n()) {
        sum(valid_scores)
      } else {
        NA_real_
      }
    },
    
    mean_log_score = {
      valid_scores <- log_predictive_density[
        is.finite(log_predictive_density)
      ]
      
      if (length(valid_scores) == n()) {
        mean(valid_scores)
      } else {
        NA_real_
      }
    },
    
    se_elpd = {
      valid_scores <- log_predictive_density[
        is.finite(log_predictive_density)
      ]
      
      if (
        length(valid_scores) == n() &&
        length(valid_scores) > 1
      ) {
        sqrt(length(valid_scores)) *
          sd(valid_scores)
      } else {
        NA_real_
      }
    },
    
    rmse_event = {
      valid_rmse <-
        event == 1 &
        is.finite(predicted_median) &
        is.finite(recovery_days)
      
      if (sum(valid_rmse) > 0) {
        sqrt(
          mean(
            (
              predicted_median[valid_rmse] -
                recovery_days[valid_rmse]
            )^2
          )
        )
      } else {
        NA_real_
      }
    },
    
    mae_event = {
      valid_mae <-
        event == 1 &
        is.finite(predicted_median) &
        is.finite(recovery_days)
      
      if (sum(valid_mae) > 0) {
        mean(
          abs(
            predicted_median[valid_mae] -
              recovery_days[valid_mae]
          )
        )
      } else {
        NA_real_
      }
    },
    
    mean_error_event = {
      valid_error <-
        event == 1 &
        is.finite(predicted_median) &
        is.finite(recovery_days)
      
      if (sum(valid_error) > 0) {
        mean(
          predicted_median[valid_error] -
            recovery_days[valid_error]
        )
      } else {
        NA_real_
      }
    },
    
    .groups = "drop"
  ) %>%
  filter(
    failed_log_predictions == 0
  ) %>%
  arrange(
    desc(elpd_loo)
  ) %>%
  mutate(
    delta_elpd =
      max(elpd_loo, na.rm = TRUE) -
      elpd_loo,
    negative_mean_log_score =
      -mean_log_score
  )

aft_loocv_comparison

# ============================================================
# BEST LOOCV MODEL
# ============================================================

best_aft_cv_model <- aft_loocv_comparison %>%
  slice_max(
    order_by = elpd_loo,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::select(
    distribution,
    model
  )

best_aft_cv_model

# ============================================================
# POINTWISE ELPD DIFFERENCES
# ============================================================

best_model_scores <- aft_loocv_predictions %>%
  semi_join(
    best_aft_cv_model,
    by = c(
      "distribution",
      "model"
    )
  ) %>%
  dplyr::select(
    observation,
    best_log_predictive_density =
      log_predictive_density
  )

aft_loocv_differences <- aft_loocv_predictions %>%
  filter(
    !is.na(log_predictive_density)
  ) %>%
  left_join(
    best_model_scores,
    by = "observation"
  ) %>%
  mutate(
    pointwise_elpd_difference =
      log_predictive_density -
      best_log_predictive_density
  ) %>%
  group_by(
    distribution,
    model
  ) %>%
  summarise(
    elpd_difference = sum(
      pointwise_elpd_difference
    ),
    se_difference = sqrt(
      n() * var(
        pointwise_elpd_difference
      )
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(elpd_difference)
  )

aft_loocv_differences

# ============================================================
# COMBINE AIC AND LOOCV RESULTS
# ============================================================

aft_model_comparison <- aft_comparison %>%
  dplyr::select(
    distribution,
    model,
    parameters,
    AIC,
    delta_AIC,
    AIC_weight
  ) %>%
  left_join(
    aft_loocv_comparison %>%
      dplyr::select(
        distribution,
        model,
        elpd_loo,
        se_elpd,
        delta_elpd,
        mean_log_score,
        rmse_event,
        mae_event,
        mean_error_event,
        n_events,
        n_events_used_rmse,
        failed_median_predictions
      ),
    by = c(
      "distribution",
      "model"
    )
  ) %>%
  left_join(
    aft_loocv_differences %>%
      dplyr::select(
        distribution,
        model,
        elpd_difference,
        se_difference
      ),
    by = c(
      "distribution",
      "model"
    )
  ) %>%
  arrange(
    desc(elpd_loo)
  )

aft_model_comparison

# ============================================================
# PREDICTIVE PERFORMANCE TABLE
# ============================================================

aft_predictive_table <- aft_model_comparison %>%
  mutate(
    AIC = round(AIC, 1),
    delta_AIC = round(delta_AIC, 2),
    AIC_weight = round(AIC_weight, 3),
    elpd_loo = round(elpd_loo, 1),
    se_elpd = round(se_elpd, 1),
    delta_elpd = round(delta_elpd, 1),
    elpd_difference = round(elpd_difference, 1),
    se_difference = round(se_difference, 1),
    rmse_event = round(rmse_event, 1),
    mae_event = round(mae_event, 1),
    mean_error_event = round(mean_error_event, 1)
  ) %>%
  dplyr::select(
    distribution,
    model,
    AIC,
    delta_AIC,
    AIC_weight,
    elpd_loo,
    se_elpd,
    elpd_difference,
    se_difference,
    rmse_event,
    mae_event,
    mean_error_event,
    n_events_used_rmse,
    failed_median_predictions
  )

print(
  aft_predictive_table,
  n = 42
)

aft_predictive_table %>%
  filter(
    !is.na(rmse_event)
  ) %>%
  arrange(
    rmse_event
  ) %>%
  dplyr::select(
    distribution,
    model,
    rmse_event,
    mae_event,
    mean_error_event,
    n_events_used_rmse,
    elpd_loo
  ) %>%
  print(n = 42)

# ============================================================
# SIMPLIFIED PREDICTIVE TABLE
# ============================================================

aft_predictive_table <- aft_model_comparison %>%
  mutate(
    AIC = round(AIC, 1),
    delta_AIC = round(delta_AIC, 2),
    AIC_weight = round(AIC_weight, 3),
    elpd_loo = round(elpd_loo, 1),
    se_elpd = round(se_elpd, 1),
    rmse_event = round(rmse_event, 1),
    mae_event = round(mae_event, 1)
  ) %>%
  dplyr::select(
    distribution,
    model,
    AIC,
    delta_AIC,
    AIC_weight,
    elpd_loo,
    se_elpd,
    rmse_event,
    mae_event
  ) %>%
  arrange(
    desc(elpd_loo)
  )

print(aft_predictive_table, n = 42)