# ============================================================
# INTERVAL-CENSORED RECOVERY MODELS USING icenReg
#
# Models:
#   1. Log-normal AFT
#   2. Weibull AFT
#   3. Log-logistic AFT
#   4. Semi-parametric interval-censored PH
#
# Candidate specifications:
#   M1-M14
#
# Prediction assessment:
#   Leave-one-out cross-validation (LOOCV)
#
# event = 1 : recovery observed between monthly observations
# event = 0 : subsequent earthquake -> right censoring
# ============================================================


# ============================================================
# PACKAGES
# ============================================================
install.packages(c('icenReg','forcats'))
library(icenReg)
library(survival)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lubridate)
library(ggplot2)
library(forcats)
library(readr)
# ============================================================
# CONSTRUCT INTERVAL DATA
# ============================================================
pga_recovery_econ_clean1 = read_csv('/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv')
data_table <- pga_recovery_econ_clean1 %>%
  mutate(
    population = population / 1e6,        # million persons
    gdp2015 = gdp2015 / 1e6,             # trillion yen (from million yen)
    dis_expen = dis_expen / 1e6,         # billion yen (from thousand yen)
    maintain_expen = maintain_expen / 1e6,
    restore_expen = restore_expen / 1e6
  )
data_table <- data_table %>% mutate(event = ifelse(first_event == "recovered", 1, 0))
colnames(data_table)
interval_data <- data_table %>% mutate(date = as.Date(date),
                                       observed_date = date + recovery_days,
                                       previous_month_date = floor_date(observed_date, unit = "month") %m-% months(1),
                                       # Interval endpoints measured as days since earthquake
                                       recovery_left = case_when(event == 1 ~  as.numeric(previous_month_date - date),
                                                                 # Subsequent earthquake: right censored at its occurrence time
                                                                 event == 0 ~ recovery_days, TRUE ~ NA_real_),
                                       recovery_right = case_when(event == 1 ~ recovery_days, event == 0 ~ Inf, TRUE ~ NA_real_),
                                       recovery_left = pmax(recovery_left, 0),
                                       censoring_type = case_when(
                                         event == 1 ~ "Interval-censored recovery",
                                         event == 0 ~ "Right-censored at subsequent earthquake",
                                         TRUE ~ NA_character_),
                                       log_gdp = log(gdp2015))
interval_analysis_data = interval_data

interval_surv_formula <- function(rhs) {as.formula(
  paste("Surv(recovery_left, recovery_right, type = 'interval2') ~", rhs))}

interval_candidate_formulas <- list(
  M1 = interval_surv_formula("magnitude + depth"),
  M2 = interval_surv_formula("magnitude + depth + mean_pga"),
  M3 = interval_surv_formula("magnitude + depth + max_pga"),
  M4 = interval_surv_formula("magnitude + depth + ntl_decline_ratio"),
  M5 = interval_surv_formula("magnitude + depth + gdp2015 + population + dis_expen"),
  M6 = interval_surv_formula("magnitude + depth + mean_pga + ntl_decline_ratio"),
  M7 = interval_surv_formula("magnitude + depth + max_pga + ntl_decline_ratio"),
  M8 = interval_surv_formula("magnitude + depth + ntl_decline_ratio + gdp2015 + population + dis_expen"),
  M9 = interval_surv_formula("magnitude + depth + mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen"),
  M10 = interval_surv_formula("magnitude + depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen"),
  M11 = interval_surv_formula("magnitude + depth + log_gdp + population + dis_expen"),
  M12 = interval_surv_formula("magnitude + depth + ntl_decline_ratio + log_gdp + population + dis_expen"),
  M13 = interval_surv_formula("magnitude + depth + mean_pga + ntl_decline_ratio + log_gdp + population + dis_expen"),
  M14 = interval_surv_formula("magnitude + depth + max_pga + ntl_decline_ratio + log_gdp + population + dis_expen"))

# ============================================================
# PARAMETRIC AFT FITTING FUNCTION
# ============================================================

fit_ic_par <- function(formula, model_name, distribution, distribution_label, data) {
  tryCatch(
    {
      fit <- icenReg::ic_par(formula = formula, data = data, model = "aft", dist = distribution)
      list(Model = model_name, Distribution = distribution_label,
        Fit = fit, Successful = TRUE, Error = NA_character_)},
    error = function(e) {
      list(Model = model_name, Distribution = distribution_label,
        Fit = NULL, Successful = FALSE, Error = conditionMessage(e))})}

# ============================================================
# LOG-NORMAL
# ============================================================

icen_lognormal_results <- imap(
  interval_candidate_formulas,
  ~ fit_ic_par(formula = .x, model_name = .y, distribution = "lnorm",
    distribution_label = "Log-normal AFT", data = interval_analysis_data))

# ============================================================
# WEIBULL
# ============================================================

icen_weibull_results <- imap(
  interval_candidate_formulas,
  ~ fit_ic_par(formula = .x,  model_name = .y, distribution = "weibull",
    distribution_label = "Weibull AFT", data = interval_analysis_data))

# ============================================================
# LOG-LOGISTIC
# ============================================================

icen_loglogistic_results <- imap(
  interval_candidate_formulas,
  ~ fit_ic_par(formula = .x, model_name = .y, distribution = "loglogistic", 
    distribution_label = "Log-logistic AFT", data = interval_analysis_data))

# ============================================================
# SEMI-PARAMETRIC PH
# ============================================================

fit_ic_ph <- function(formula, model_name, data) {
  tryCatch(
    {
      fit <- icenReg::ic_sp(formula = formula, data = data,  model = "ph", bs_samples = 0)
      list(Model = model_name, Distribution = "Interval PH",  Fit = fit,
           Successful = TRUE, Error = NA_character_)},
    error = function(e) {
      list(
        Model = model_name, Distribution = "Interval PH", Fit = NULL,
        Successful = FALSE, Error = conditionMessage(e))})}


icen_ph_results <- imap(interval_candidate_formulas, ~ fit_ic_ph(
    formula = .x, model_name = .y, data = interval_analysis_data))

# ============================================================
# FIT STATUS
# ============================================================

extract_fit_status <- function(results) {
  imap_dfr(results, function(x, model_name) {
      tibble(
        Model = model_name,
        Distribution = x$Distribution,
        Successful = x$Successful,
        Error = x$Error)})}


interval_fit_status <- bind_rows(
  extract_fit_status(icen_lognormal_results),
  extract_fit_status(icen_weibull_results),
  extract_fit_status(icen_loglogistic_results),
  extract_fit_status(icen_ph_results))


print(interval_fit_status, n = Inf, width = Inf)
interval_fit_status %>% count(Distribution, Successful)

# ============================================================
# SUCCESSFUL FITS
# ============================================================

icen_lognormal_models <-
  icen_lognormal_results %>%
  keep(~ .x$Successful) %>%
  map("Fit")


icen_weibull_models <-
  icen_weibull_results %>%
  keep(~ .x$Successful) %>%
  map("Fit")


icen_loglogistic_models <-
  icen_loglogistic_results %>%
  keep(~ .x$Successful) %>%
  map("Fit")


icen_ph_models <-
  icen_ph_results %>%
  keep(~ .x$Successful) %>%
  map("Fit")

length(icen_lognormal_models)
length(icen_weibull_models)
length(icen_loglogistic_models)
length(icen_ph_models)


# ============================================================
# PARAMETRIC AIC
# ============================================================

extract_ic_par_statistics <- function(fit, model_name, distribution) {
  n_regression <- length(fit$reg_pars)
  n_baseline <- length(fit$baseline)
  k <- n_regression + n_baseline
  ll <- as.numeric(fit$llk)
  tibble(
    Model = model_name,
    Distribution = distribution,
    Parameters = k,
    LogLik = ll,
    AIC = -2 * ll + 2 * k)}

parametric_fit_table <- bind_rows(
  imap_dfr(
    icen_lognormal_models,
    ~ extract_ic_par_statistics(
      fit = .x,
      model_name = .y,
      distribution = "Log-normal AFT")),
  imap_dfr(
    icen_weibull_models,
    ~ extract_ic_par_statistics(
      fit = .x,
      model_name = .y,
      distribution = "Weibull AFT")),
  imap_dfr(
    icen_loglogistic_models,
    ~ extract_ic_par_statistics(
      fit = .x,
      model_name = .y,
      distribution = "Log-logistic AFT"))
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Relative_Likelihood = exp(-0.5 * Delta_AIC),
    AIC_Weight = Relative_Likelihood / sum(Relative_Likelihood),
    AIC_Rank = row_number()) %>% select(-Relative_Likelihood)

print(parametric_fit_table, n = Inf)

# ============================================================
# INTERVAL-AWARE ERROR FUNCTION
# ============================================================

calculate_interval_error <- function(predicted_median, left, right) {
  if (!is.finite(predicted_median)) {
    return(NA_real_)}
  # ----------------------------------------------
  # Interval-censored recovery
  # ----------------------------------------------
  if (is.finite(right)) {
    if (predicted_median < left) {
      return(left - predicted_median)}
    if (predicted_median > right) {
      return(predicted_median - right)}
    # Prediction lies inside observed interval
    return(0)}
  # ----------------------------------------------
  # Right-censored observation
  # ----------------------------------------------
  if (is.infinite(right)) {
    if (predicted_median < left) {
      return(left - predicted_median)}
    # Prediction is compatible with T > left
    return(0)}
  NA_real_}

# ============================================================
# LOOCV PARAMETRIC AFT
# ============================================================

loocv_ic_par <- function(formula, model_name, distribution, distribution_label, data) {
  n <- nrow(data)
  results <-  vector("list", n)
  for (i in seq_len(n)) {
    train_data <- data[-i, , drop = FALSE]
    test_data <- data[i,  , drop = FALSE]
    results[[i]] <- tryCatch(
      {
        fit <- icenReg::ic_par(
          formula = formula,
          data = train_data,
          model = "aft",
          dist = distribution
        )
        predicted_median <- predict(fit, newdata = test_data, type = "response")
        predicted_median <- as.numeric(predicted_median)
        left <- as.numeric(test_data$recovery_left)
        right <- as.numeric(test_data$recovery_right)
        interval_error <- calculate_interval_error(predicted_median = predicted_median,
            left =  left, right = right)
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = distribution_label,
          Recovery_left = left,
          Recovery_right = right,
          Predicted_median = predicted_median,
          Interval_error = interval_error,
          Squared_interval_error = interval_error^2,
          Successful = is.finite(predicted_median) && is.finite(interval_error),
          Error = NA_character_)},
      error = function(e) {
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = distribution_label,
          Recovery_left = as.numeric(test_data$recovery_left),
          Recovery_right = as.numeric(test_data$recovery_right),
          Predicted_median = NA_real_,
          Interval_error = NA_real_,
          Squared_interval_error = NA_real_,
          Successful = FALSE,
          Error = conditionMessage(e))})}
  bind_rows(results)
}

set.seed(2026)

loocv_lognormal <- imap_dfr(interval_candidate_formulas,
  ~ loocv_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "lnorm",
    distribution_label = "Log-normal AFT",
    data = interval_analysis_data))

loocv_weibull <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "weibull",
    distribution_label = "Weibull AFT",
    data = interval_analysis_data))

loocv_loglogistic <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "loglogistic",
    distribution_label = "Log-logistic AFT",
    data = interval_analysis_data))

# ============================================================
# LOOCV INTERVAL PH
# ============================================================

loocv_ic_ph <- function(formula, model_name, data) {
  n <- nrow(data)
  results <- vector("list", n)
  for (i in seq_len(n)) {
    train_data <- data[-i, , drop = FALSE]
    test_data <- data[i, , drop = FALSE]
    results[[i]] <- tryCatch(
      {
        fit <- icenReg::ic_sp(
          formula = formula,
          data = train_data,
          model = "ph",
          bs_samples = 0
        )
        predicted_median <- predict(fit, newdata = test_data, type = "response")
        predicted_median <- as.numeric(predicted_median)
        left <- as.numeric(test_data$recovery_left)
        right <- as.numeric(test_data$recovery_right)
        interval_error <- calculate_interval_error(predicted_median, left, right)
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = "Interval PH",
          Recovery_left = left,
          Recovery_right = right,
          Predicted_median = predicted_median,
          Interval_error = interval_error,
          Squared_interval_error = interval_error^2,
          Successful = is.finite(predicted_median) && is.finite(interval_error),
          Error = NA_character_)},
      error = function(e) {
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = "Interval PH",
          Recovery_left = as.numeric(test_data$recovery_left),
          Recovery_right = as.numeric(test_data$recovery_right),
          Predicted_median = NA_real_,
          Interval_error = NA_real_,
          Squared_interval_error = NA_real_,
          Successful =  FALSE,
          Error = conditionMessage(e)
        )
      }
    )
  }
  
  
  bind_rows(results)
}

loocv_ph <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_ic_ph(
    formula = .x,
    model_name = .y,
    data = interval_analysis_data))

# ============================================================
# ALL LOOCV RESULTS
# ============================================================

loocv_predictions <- bind_rows(
  loocv_lognormal,
  loocv_weibull,
  loocv_loglogistic,
  loocv_ph)


print(loocv_predictions, n = 100)

loocv_predictions %>% group_by(Model, Distribution) %>%
  summarise(N = n(),
            Successful = sum(Successful),
            Failed = sum(!Successful),
            Finite_predictions = sum(is.finite(Predicted_median)),
    .groups = "drop") %>% print(n = Inf)
loocv_predictions %>%
  filter(
    !Successful
  ) %>%
  count(Model, Distribution, Error, sort = TRUE) %>%
  print(n = Inf, width = Inf)

# ============================================================
# LOOCV PERFORMANCE SUMMARY
# ============================================================

loocv_performance <- loocv_predictions %>% group_by(Model, Distribution) %>%
  summarise(N = n(),
    N_successful =sum(Successful & is.finite(Interval_error)),
    Success_rate = N_successful / N,
    LOOCV_MAE = if (any(is.finite(Interval_error))) {
        mean(Interval_error[is.finite(Interval_error)])
      } else {NA_real_},
    LOOCV_RMSE =
      if (any(is.finite(Interval_error))) {sqrt(mean(Interval_error[is.finite(Interval_error)]^2))
      } else {NA_real_},
    Median_Error = if (any(is.finite(Interval_error))) {
        median(Interval_error[is.finite(Interval_error)])
      } else {NA_real_},
    .groups = "drop"
  ) %>%
  arrange(desc(Success_rate),LOOCV_MAE,LOOCV_RMSE) %>%
  mutate(Prediction_Rank =if_else(Success_rate == 1 & is.finite(LOOCV_MAE),
        rank(LOOCV_MAE,ties.method = "first"), NA_integer_))


print(loocv_performance, n = Inf, width = Inf)

# ============================================================
# OVERALL PREDICTIVE COMPARISON
# ============================================================

overall_prediction_table <- loocv_performance %>%
  arrange(desc(Success_rate),LOOCV_MAE,LOOCV_RMSE) %>%
  transmute(Model, Distribution,
    `LOOCV success (%)` = round(100 * Success_rate, 1),
    `LOOCV MAE` = round(LOOCV_MAE, 2),
    `LOOCV RMSE` = round(LOOCV_RMSE, 2),
    `Median error` = round(Median_Error, 2),
    `Prediction rank` = Prediction_Rank
  )


print(overall_prediction_table, n = Inf, width = Inf)

best_predictive_model <- loocv_performance %>%
  filter(Success_rate == 1, is.finite(LOOCV_MAE)) %>%
  slice_min(LOOCV_MAE, n = 1, with_ties = FALSE)

best_predictive_model


best_parametric_aic <- parametric_fit_table %>%
  slice_min(AIC, n = 1, with_ties = FALSE)
best_parametric_aic

loocv_mae_plot <- loocv_performance %>%
  filter(is.finite(LOOCV_MAE)) %>%
  ggplot(
    aes(
      x = reorder(paste(Model, Distribution), LOOCV_MAE),
      y = LOOCV_MAE,
      fill = Distribution)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "LOOCV predictive performance",
    subtitle = "Interval-censored recovery models",
    x = NULL,
    y = "Interval-aware MAE (days)",
    fill = "Model family") +
  theme_bw()

loocv_mae_plot


loocv_rmse_plot <- loocv_performance %>% filter(
    is.finite(LOOCV_RMSE)) %>%
  ggplot(aes(x = reorder(paste(Model, Distribution), LOOCV_RMSE),
      y = LOOCV_RMSE,
      fill = Distribution)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "LOOCV prediction RMSE",
    subtitle = "Interval-censored recovery models",
    x = NULL,
    y = "Interval-aware RMSE (days)",
    fill = "Model family") +
  theme_bw()


loocv_rmse_plot

# ============================================================
# SAVE RESULTS
# ============================================================

interval_icenreg_results <- list(
  data = interval_analysis_data,
  formulas = interval_candidate_formulas,
  fit_status = interval_fit_status,
  lognormal_models = icen_lognormal_models,
  weibull_models = icen_weibull_models,
  loglogistic_models = icen_loglogistic_models,
  interval_ph_models = icen_ph_models,
  parametric_fit_table = parametric_fit_table,
  loocv_predictions = loocv_predictions,
  loocv_performance = loocv_performance,
  parametric_final_table = parametric_final_table,
  overall_prediction_table = overall_prediction_table,
  best_predictive_model = best_predictive_model,
  best_parametric_aic = best_parametric_aic
)

best_predictive_model
best_parametric_aic
saveRDS(interval_icenreg_results, file = "interval_censored_icenreg_results.rds")


# ============================================================
# PARAMETRIC AFT SURVIVAL PROBABILITY
# ============================================================

get_icpar_survival <- function(fit, newdata, time, distribution, eps = 1e-12) {
  if (time <= 0) return(1)
  
  eta <- as.numeric(predict(fit, newdata = newdata, type = "lp"))
  baseline <- fit$baseline
  
  if (distribution == "lnorm") {
    sigma <- exp(as.numeric(baseline[1]))
    z <- (log(time) - eta) / sigma
    surv <- 1 - pnorm(z)
    
  } else if (distribution == "weibull") {
    shape <- exp(as.numeric(baseline[1]))
    scale <- exp(eta)
    surv <- exp(-(time / scale)^shape)
    
  } else if (distribution == "loglogistic") {
    shape <- exp(as.numeric(baseline[1]))
    scale <- exp(eta)
    surv <- 1 / (1 + (time / scale)^shape)
    
  } else {
    stop(paste("Unsupported distribution:", distribution))
  }
  
  pmin(pmax(surv, eps), 1)
}

# ============================================================
# ELPD LOOCV FOR PARAMETRIC INTERVAL-CENSORED AFT
# ============================================================

loocv_elpd_ic_par <- function(formula, model_name, distribution,
                              distribution_label, data, eps = 1e-12) {
  n <- nrow(data)
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    train_data <- data[-i, , drop = FALSE]
    test_data <- data[i, , drop = FALSE]
    
    results[[i]] <- tryCatch(
      {
        fit <- icenReg::ic_par(
          formula = formula,
          data = train_data,
          model = "aft",
          dist = distribution
        )
        
        left <- as.numeric(test_data$recovery_left)
        right <- as.numeric(test_data$recovery_right)
        
        # Interval-censored: P(L < T <= R) = S(L) - S(R)
        if (is.finite(right)) {
          survival_left <- get_icpar_survival(
            fit, test_data, left, distribution, eps
          )
          
          survival_right <- get_icpar_survival(
            fit, test_data, right, distribution, eps
          )
          
          predictive_probability <- survival_left - survival_right
          
          # Right-censored: P(T > L) = S(L)
        } else if (is.infinite(right)) {
          predictive_probability <- get_icpar_survival(
            fit, test_data, left, distribution, eps
          )
          
        } else {
          predictive_probability <- NA_real_
        }
        
        if (is.finite(predictive_probability)) {
          predictive_probability <- max(predictive_probability, eps)
        }
        
        log_predictive_density <- if (is.finite(predictive_probability)) {
          log(predictive_probability)
        } else {
          NA_real_
        }
        
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = distribution_label,
          Recovery_left = left,
          Recovery_right = right,
          Predictive_probability = predictive_probability,
          Log_predictive_density = log_predictive_density,
          Successful = is.finite(log_predictive_density),
          Error = NA_character_
        )
      },
      error = function(e) {
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = distribution_label,
          Recovery_left = as.numeric(test_data$recovery_left),
          Recovery_right = as.numeric(test_data$recovery_right),
          Predictive_probability = NA_real_,
          Log_predictive_density = NA_real_,
          Successful = FALSE,
          Error = conditionMessage(e)
        )
      }
    )
  }
  
  bind_rows(results)
}

# ============================================================
# TEST M1 LOG-NORMAL
# ============================================================

test_elpd <- loocv_elpd_ic_par(
  formula = interval_candidate_formulas$M1,
  model_name = "M1",
  distribution = "lnorm",
  distribution_label = "Log-normal AFT",
  data = interval_analysis_data
)

test_elpd %>% count(Successful, Error)
# ============================================================
# LOG-NORMAL ELPD LOOCV
# ============================================================

elpd_lognormal <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_elpd_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "lnorm",
    distribution_label = "Log-normal AFT",
    data = interval_analysis_data
  )
)
# ============================================================
# WEIBULL ELPD LOOCV
# ============================================================

elpd_weibull <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_elpd_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "weibull",
    distribution_label = "Weibull AFT",
    data = interval_analysis_data
  )
)

# ============================================================
# LOG-LOGISTIC ELPD LOOCV
# ============================================================

elpd_loglogistic <- imap_dfr(
  interval_candidate_formulas,
  ~ loocv_elpd_ic_par(
    formula = .x,
    model_name = .y,
    distribution = "loglogistic",
    distribution_label = "Log-logistic AFT",
    data = interval_analysis_data
  )
)

elpd_loocv_predictions <- bind_rows(
  elpd_lognormal,
  elpd_weibull,
  elpd_loglogistic
)
elpd_loocv_predictions %>% group_by(Model, Distribution) %>%
  summarise(
    N = n(),
    N_successful = sum(Successful),
    N_failed = sum(!Successful),
    .groups = "drop") %>%print(n = Inf)

elpd_loocv_predictions %>% filter(!Successful) %>%
  count(Model, Distribution, Error, sort = TRUE) %>%
  print(n = Inf, width = Inf)
# ============================================================
# ELPD LOOCV SUMMARY
# ============================================================

elpd_loocv_summary <- elpd_loocv_predictions %>%
  group_by(Model, Distribution) %>%
  summarise(
    N = n(),
    N_successful = sum(is.finite(Log_predictive_density)),
    Success_rate = N_successful / N,
    ELPD_LOOCV = if (all(is.finite(Log_predictive_density))) {
      sum(Log_predictive_density)
    } else {
      NA_real_
    },
    Mean_log_predictive_density = mean(
      Log_predictive_density,
      na.rm = TRUE
    ),
    SE_ELPD = if (sum(is.finite(Log_predictive_density)) > 1) {
      sqrt(
        N * var(
          Log_predictive_density,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  arrange(desc(ELPD_LOOCV)) %>%
  mutate(ELPD_Rank = row_number())

print(elpd_loocv_summary, n = Inf, width = Inf)


best_elpd_model <- elpd_loocv_summary %>%
  filter(is.finite(ELPD_LOOCV)) %>%
  slice_max(
    ELPD_LOOCV,
    n = 1,
    with_ties = FALSE
  )

best_elpd_model


best_elpd_model_name <- best_elpd_model$Model
best_elpd_distribution <- best_elpd_model$Distribution


best_elpd_pointwise <- elpd_loocv_predictions %>%
  filter(
    Model == best_elpd_model_name,
    Distribution == best_elpd_distribution
  ) %>%
  select(
    Observation,
    Best_Log_predictive_density = Log_predictive_density
  )


elpd_difference_summary <- elpd_loocv_predictions %>%
  left_join(
    best_elpd_pointwise,
    by = "Observation"
  ) %>%
  mutate(
    Pointwise_ELPD_Difference =
      Log_predictive_density - Best_Log_predictive_density
  ) %>%
  group_by(Model, Distribution) %>%
  summarise(
    Delta_ELPD = sum(
      Pointwise_ELPD_Difference,
      na.rm = TRUE
    ),
    SE_Delta_ELPD = sqrt(
      n() * var(
        Pointwise_ELPD_Difference,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )


elpd_loocv_summary <- elpd_loocv_summary %>%
  left_join(
    elpd_difference_summary,
    by = c("Model", "Distribution")
  ) %>%
  arrange(desc(ELPD_LOOCV)) %>%
  mutate(ELPD_Rank = row_number())

print(elpd_loocv_summary, n = Inf, width = Inf)


# ============================================================
# JOIN MAE/RMSE + ELPD
# ============================================================

combined_loocv_performance <- loocv_performance %>%
  left_join(
    elpd_loocv_summary %>%
      select(
        Model,
        Distribution,
        ELPD_LOOCV,
        SE_ELPD,
        Delta_ELPD,
        SE_Delta_ELPD,
        ELPD_Rank
      ),
    by = c("Model", "Distribution")
  )

print(combined_loocv_performance, n = Inf, width = Inf)


# ============================================================
# FINAL PARAMETRIC COMPARISON
# ============================================================

final_parametric_comparison <- parametric_fit_table %>%
  left_join(
    combined_loocv_performance,
    by = c("Model", "Distribution")
  ) %>%
  arrange(
    desc(ELPD_LOOCV),
    LOOCV_MAE,
    AIC
  )

print(final_parametric_comparison, n = Inf, width = Inf)


# ============================================================
# LOG-DEPTH SENSITIVITY ANALYSIS FOR INTERVAL PH
# Uses existing fit_ic_ph()
# ============================================================

interval_analysis_data_logdepth <- interval_analysis_data %>%
  mutate(
    log_depth = log(depth)
  )

interval_candidate_formulas_logdepth <- list(
  M1 = interval_surv_formula("magnitude + log_depth"),
  M2 = interval_surv_formula("magnitude + log_depth + mean_pga"),
  M3 = interval_surv_formula("magnitude + log_depth + max_pga"),
  M4 = interval_surv_formula("magnitude + log_depth + ntl_decline_ratio"),
  M5 = interval_surv_formula("magnitude + log_depth + gdp2015 + population + dis_expen"),
  M6 = interval_surv_formula("magnitude + log_depth + mean_pga + ntl_decline_ratio"),
  M7 = interval_surv_formula("magnitude + log_depth + max_pga + ntl_decline_ratio"),
  M8 = interval_surv_formula(
    "magnitude + log_depth + ntl_decline_ratio + gdp2015 + population + dis_expen"
  ),
  M9 = interval_surv_formula(
    "magnitude + log_depth + mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen"
  ),
  M10 = interval_surv_formula(
    "magnitude + log_depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen"
  ),
  M11 = interval_surv_formula(
    "magnitude + log_depth + log_gdp + population + dis_expen"
  ),
  M12 = interval_surv_formula(
    "magnitude + log_depth + ntl_decline_ratio + log_gdp + population + dis_expen"
  ),
  M13 = interval_surv_formula(
    "magnitude + log_depth + mean_pga + ntl_decline_ratio + log_gdp + population + dis_expen"
  ),
  M14 = interval_surv_formula(
    "magnitude + log_depth + max_pga + ntl_decline_ratio + log_gdp + population + dis_expen"
  )
)

# ============================================================
# FIT USING YOUR EXISTING FUNCTION
# ============================================================

icen_ph_logdepth_results <- imap(
  interval_candidate_formulas_logdepth,
  ~ fit_ic_ph(
    formula = .x,
    model_name = .y,
    data = interval_analysis_data_logdepth
  )
)

# ============================================================
# FIT STATUS USING EXISTING extract_fit_status()
# ============================================================

icen_ph_logdepth_status <- extract_fit_status(
  icen_ph_logdepth_results
)

print(
  icen_ph_logdepth_status,
  n = Inf,
  width = Inf
)

# ============================================================
# KEEP SUCCESSFUL PH MODELS
# ============================================================

icen_ph_logdepth_models <- icen_ph_logdepth_results %>%
  keep(~ .x$Successful) %>%
  map("Fit")

names(icen_ph_logdepth_models)
length(icen_ph_logdepth_models)

# ============================================================
# LOOCV FOR LOG-DEPTH INTERVAL PH
# ============================================================

loocv_ph_logdepth <- imap_dfr(
  interval_candidate_formulas_logdepth,
  ~ loocv_ic_ph(
    formula = .x,
    model_name = .y,
    data = interval_analysis_data_logdepth
  )
)

loocv_ph_logdepth

# ============================================================
# LOG-DEPTH PH LOOCV PERFORMANCE
# ============================================================

loocv_ph_logdepth_performance <- loocv_ph_logdepth %>%
  group_by(Model, Distribution) %>%
  summarise(
    N = n(),
    N_successful = sum(Successful & is.finite(Interval_error)),
    Success_rate = N_successful / N,
    LOOCV_MAE = if (any(is.finite(Interval_error))) {
      mean(Interval_error[is.finite(Interval_error)])
    } else {
      NA_real_
    },
    LOOCV_RMSE = if (any(is.finite(Interval_error))) {
      sqrt(mean(Interval_error[is.finite(Interval_error)]^2))
    } else {
      NA_real_
    },
    Median_Error = if (any(is.finite(Interval_error))) {
      median(Interval_error[is.finite(Interval_error)])
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  arrange(
    desc(Success_rate),
    LOOCV_MAE,
    LOOCV_RMSE
  ) %>%
  mutate(
    Prediction_Rank = if_else(
      Success_rate == 1 & is.finite(LOOCV_MAE),
      rank(LOOCV_MAE, ties.method = "first"),
      NA_integer_
    )
  )

print(loocv_ph_logdepth_performance, n = Inf, width = Inf)

colnames(final_parametric_comparison)

# ============================================================
# BRIER SCORE FOR BEST MODEL IN EACH PARAMETRIC DISTRIBUTION
# ============================================================

best_each_distribution <- final_parametric_comparison %>%
  filter(is.finite(ELPD_LOOCV)) %>%
  group_by(Distribution) %>%
  slice_max(ELPD_LOOCV, n = 1, with_ties = FALSE) %>%
  ungroup()

best_each_distribution

brier_times <- seq(10, 290, by = 10)

loocv_brier_ic_par <- function(formula, model_name, distribution,
                               distribution_label, data, times) {
  n <- nrow(data)
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    train_data <- data[-i, , drop = FALSE]
    test_data <- data[i, , drop = FALSE]
    
    results[[i]] <- tryCatch(
      {
        fit <- icenReg::ic_par(
          formula = formula,
          data = train_data,
          model = "aft",
          dist = distribution
        )
        
        left <- as.numeric(test_data$recovery_left)
        right <- as.numeric(test_data$recovery_right)
        
        map_dfr(times, function(t) {
          survival_probability <- get_icpar_survival(
            fit = fit,
            newdata = test_data,
            time = t,
            distribution = distribution
          )
          
          predicted_recovery_probability <- 1 - survival_probability
          
          observed_status <- case_when(
            t <= left ~ 0,
            is.finite(right) & t >= right ~ 1,
            TRUE ~ NA_real_
          )
          
          brier <- if (is.finite(observed_status)) {
            (observed_status - predicted_recovery_probability)^2
          } else {
            NA_real_
          }
          
          tibble(
            Observation = i,
            Model = model_name,
            Distribution = distribution_label,
            Time = t,
            Recovery_left = left,
            Recovery_right = right,
            Observed_status = observed_status,
            Predicted_probability = predicted_recovery_probability,
            Brier = brier,
            Evaluable = is.finite(brier)
          )
        })
      },
      error = function(e) {
        tibble(
          Observation = i,
          Model = model_name,
          Distribution = distribution_label,
          Time = times,
          Recovery_left = as.numeric(test_data$recovery_left),
          Recovery_right = as.numeric(test_data$recovery_right),
          Observed_status = NA_real_,
          Predicted_probability = NA_real_,
          Brier = NA_real_,
          Evaluable = FALSE
        )
      }
    )
  }
  
  bind_rows(results)
}

best_lognormal_name <- best_each_distribution %>%
  filter(Distribution == "Log-normal AFT") %>%
  pull(Model)

best_weibull_name <- best_each_distribution %>%
  filter(Distribution == "Weibull AFT") %>%
  pull(Model)

best_loglogistic_name <- best_each_distribution %>%
  filter(Distribution == "Log-logistic AFT") %>%
  pull(Model)

best_lognormal_name
best_weibull_name
best_loglogistic_name

brier_best_lognormal <- loocv_brier_ic_par(
  formula = interval_candidate_formulas[[best_lognormal_name]],
  model_name = best_lognormal_name,
  distribution = "lnorm",
  distribution_label = "Log-normal AFT",
  data = interval_analysis_data,
  times = brier_times
)

brier_best_weibull <- loocv_brier_ic_par(
  formula = interval_candidate_formulas[[best_weibull_name]],
  model_name = best_weibull_name,
  distribution = "weibull",
  distribution_label = "Weibull AFT",
  data = interval_analysis_data,
  times = brier_times
)

brier_best_loglogistic <- loocv_brier_ic_par(
  formula = interval_candidate_formulas[[best_loglogistic_name]],
  model_name = best_loglogistic_name,
  distribution = "loglogistic",
  distribution_label = "Log-logistic AFT",
  data = interval_analysis_data,
  times = brier_times
)

best_model_brier_predictions <- bind_rows(
  brier_best_lognormal,
  brier_best_weibull,
  brier_best_loglogistic
)

best_model_brier_predictions %>%
  group_by(Model, Distribution) %>%
  summarise(
    N = n(),
    N_evaluable = sum(is.finite(Brier)),
    N_missing = sum(!is.finite(Brier)),
    .groups = "drop"
  ) %>%
  print(n = Inf)

best_model_brier_by_time <- best_model_brier_predictions %>%
  group_by(Model, Distribution, Time) %>%
  summarise(
    N_evaluable = sum(is.finite(Brier)),
    Brier = if (any(is.finite(Brier))) {
      mean(Brier[is.finite(Brier)])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

print(best_model_brier_by_time, n = Inf, width = Inf)

best_model_mean_brier <- best_model_brier_by_time %>%
  group_by(Model, Distribution) %>%
  summarise(
    Mean_LOOCV_Brier = mean(Brier, na.rm = TRUE),
    SD_Brier = sd(Brier, na.rm = TRUE),
    Min_Brier = min(Brier, na.rm = TRUE),
    Max_Brier = max(Brier, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Mean_LOOCV_Brier)

best_model_mean_brier

best_each_distribution_with_brier <- best_each_distribution %>%
  left_join(
    best_model_mean_brier,
    by = c("Model", "Distribution")
  ) %>%
  arrange(Mean_LOOCV_Brier)

best_each_distribution_with_brier

best_parametric_summary_table <- best_each_distribution_with_brier %>%
  transmute(
    Model,
    Distribution,
    AIC = round(AIC, 2),
    `LOOCV MAE` = round(LOOCV_MAE, 2),
    `LOOCV RMSE` = round(LOOCV_RMSE, 2),
    `ELPD-LOOCV` = round(ELPD_LOOCV, 2),
    `SE ELPD` = round(SE_ELPD, 2),
    `Mean LOOCV Brier` = round(Mean_LOOCV_Brier, 4),
    `SD Brier` = round(SD_Brier, 4)
  )

print(best_parametric_summary_table, n = Inf, width = Inf)

brier_plot <- ggplot(
  best_model_brier_by_time,
  aes(
    x = Time,
    y = Brier,
    group = Distribution,
    color = Distribution
  )
) +
  geom_line() +
  geom_point() +
  labs(
    title = "LOOCV Brier Scores",
    subtitle = "Best model from each parametric distribution",
    x = "Time since initial earthquake (days)",
    y = "Brier score",
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

brier_plot


# ============================================================
# COEFFICIENTS FOR BEST PARAMETRIC MODELS
# ============================================================

best_lognormal_fit <- icen_lognormal_models[["M5"]]
best_loglogistic_fit <- icen_loglogistic_models[["M12"]]
best_weibull_fit <- icen_weibull_models[["M2"]]

summary(best_lognormal_fit)
summary(best_loglogistic_fit)
summary(best_weibull_fit)

best_model_coefficients <- bind_rows(
  tibble(
    Model = "M5",
    Distribution = "Log-normal AFT",
    Variable = names(best_lognormal_fit$reg_pars),
    Estimate = as.numeric(best_lognormal_fit$reg_pars),
    Time_Ratio = exp(as.numeric(best_lognormal_fit$reg_pars))
  ),
  tibble(
    Model = "M12",
    Distribution = "Log-logistic AFT",
    Variable = names(best_loglogistic_fit$reg_pars),
    Estimate = as.numeric(best_loglogistic_fit$reg_pars),
    Time_Ratio = exp(as.numeric(best_loglogistic_fit$reg_pars))
  ),
  tibble(
    Model = "M2",
    Distribution = "Weibull AFT",
    Variable = names(best_weibull_fit$reg_pars),
    Estimate = as.numeric(best_weibull_fit$reg_pars),
    Time_Ratio = exp(as.numeric(best_weibull_fit$reg_pars))
  )
)

print(best_model_coefficients, n = Inf, width = Inf)
