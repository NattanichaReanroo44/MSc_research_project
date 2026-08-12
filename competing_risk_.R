# ============================================================
# PACKAGES
# ============================================================

install.packages('riskRegression')

library(riskRegression)
library(prodlim)
library(cmprsk)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(patchwork)


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
    restore_expen = restore_expen / 1e6
  )

fine_gray_data <- data_table %>%
  mutate(event = ifelse(first_event == "recovered", 1, 0))

fine_gray_data <- fine_gray_data %>%
  mutate(
    observed_time = as.numeric(recovery_days),
    event_type = case_when(
      event == 1 ~ 1L, # recovery
      event == 0 ~ 2L, # subsequent earthquake
      TRUE ~ NA_integer_
    ),
    log_gdp = log(gdp2015)
  )

fine_gray_data %>% count(event, event_type)


# ============================================================
# DATA SUMMARY
# ============================================================

fine_gray_data %>%
  summarise(
    n_total = n(),
    n_recovery = sum(event_type == 1L, na.rm = TRUE),
    n_subsequent_earthquake = sum(
      event_type == 2L,
      na.rm = TRUE
    ),
    n_missing_status = sum(is.na(event_type)),
    minimum_time = min(observed_time, na.rm = TRUE),
    median_time = median(observed_time, na.rm = TRUE),
    maximum_time = max(observed_time, na.rm = TRUE)
  )

all(fine_gray_data$observed_time > 0, na.rm = TRUE)
table(fine_gray_data$event_type, useNA = "ifany")

model_variables <- c(
  "observed_time",
  "event_type",
  "magnitude",
  "depth",
  "mean_pga",
  "max_pga",
  "ntl_decline_ratio",
  "gdp2015",
  "log_gdp",
  "population",
  "dis_expen"
)

analysis_data <- fine_gray_data


# ============================================================
# CANDIDATE FINE-GRAY MODELS
# ============================================================

fine_gray_candidate_formulas <- list(
  M1 = Hist(observed_time, event_type) ~ magnitude + depth,
  M2 = Hist(observed_time, event_type) ~ magnitude + depth + mean_pga,
  M3 = Hist(observed_time, event_type) ~ magnitude + depth + max_pga,
  M4 = Hist(observed_time, event_type) ~ magnitude + depth + ntl_decline_ratio,
  M5 = Hist(observed_time, event_type) ~ magnitude + depth + gdp2015 + population + dis_expen,
  M6 = Hist(observed_time, event_type) ~ magnitude + depth + mean_pga + ntl_decline_ratio,
  M7 = Hist(observed_time, event_type) ~ magnitude + depth + max_pga + ntl_decline_ratio,
  M8 = Hist(observed_time, event_type) ~ magnitude + depth + ntl_decline_ratio + gdp2015 +
    population + dis_expen,
  M9 = Hist(observed_time, event_type) ~ magnitude + depth + mean_pga + ntl_decline_ratio +
    gdp2015 + population + dis_expen,
  M10 = Hist(observed_time, event_type) ~ magnitude + depth + max_pga + ntl_decline_ratio +
    gdp2015 + population + dis_expen,
  M11 = Hist(observed_time, event_type) ~ magnitude + depth + log_gdp + population + dis_expen,
  M12 = Hist(observed_time, event_type) ~ magnitude + depth + ntl_decline_ratio +
    log_gdp + population + dis_expen,
  M13 = Hist(observed_time, event_type) ~ magnitude + depth + mean_pga + ntl_decline_ratio +
    log_gdp + population + dis_expen,
  M14 = Hist(observed_time, event_type) ~ magnitude + depth + max_pga + ntl_decline_ratio +
    log_gdp + population + dis_expen
)


# ============================================================
# FIT FINE-GRAY MODELS
# ============================================================

fit_one_fine_gray_model <- function(model_formula, model_name, data) {
  tryCatch(
    {
      fitted_model <- riskRegression::FGR(
        formula = model_formula,
        data = data,
        cause = 1
      )
      
      fitted_model$call$formula <- model_formula
      fitted_model$formula <- model_formula
      
      list(
        model_name = model_name,
        formula = model_formula,
        fit = fitted_model,
        successful = TRUE,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        model_name = model_name,
        formula = model_formula,
        fit = NULL,
        successful = FALSE,
        error_message = conditionMessage(e)
      )
    }
  )
}

fine_gray_fit_results <- imap(
  fine_gray_candidate_formulas,
  ~ fit_one_fine_gray_model(
    model_formula = .x,
    model_name = .y,
    data = analysis_data
  )
)


# ============================================================
# MODEL FIT STATUS
# ============================================================

fine_gray_fit_status <- imap_dfr(
  fine_gray_fit_results,
  function(result, model_name) {
    tibble(
      Model = model_name,
      Successful = result$successful,
      Error = result$error_message
    )
  }
)

fine_gray_models <- fine_gray_fit_results %>%
  keep(~ .x$successful) %>%
  map("fit")

names(fine_gray_models)
length(fine_gray_models)

fine_gray_models$M1$call
fine_gray_models$M13$call

summary(fine_gray_models$M1)


# ============================================================
# COEFFICIENTS AND SUBDISTRIBUTION HAZARD RATIOS
# ============================================================

extract_fine_gray_coefficients <- function(fitted_model, model_name) {
  crr_fit <- fitted_model$crrFit
  
  estimates <- as.numeric(crr_fit$coef)
  variable_names <- names(crr_fit$coef)
  standard_errors <- sqrt(diag(crr_fit$var))
  z_values <- estimates / standard_errors
  p_values <- 2 * pnorm(-abs(z_values))
  
  lower_95 <- estimates - 1.96 * standard_errors
  upper_95 <- estimates + 1.96 * standard_errors
  
  tibble(
    Model = model_name,
    Variable = variable_names,
    Estimate = estimates,
    SE = standard_errors,
    Lower_95 = lower_95,
    Upper_95 = upper_95,
    SHR = exp(estimates),
    SHR_Lower_95 = exp(lower_95),
    SHR_Upper_95 = exp(upper_95),
    z_value = z_values,
    p_value = p_values
  )
}

fine_gray_coefficient_results <- imap_dfr(
  fine_gray_models,
  extract_fine_gray_coefficients
)

fine_gray_coefficient_results

fine_gray_coefficient_table <- fine_gray_coefficient_results %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    SHR = round(SHR, 3),
    `95% CI for SHR` = sprintf(
      "%.3f--%.3f",
      SHR_Lower_95,
      SHR_Upper_95
    ),
    `p-value` = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  select(
    Model,
    Variable,
    Estimate,
    SE,
    SHR,
    `95% CI for SHR`,
    `p-value`
  )

fine_gray_coefficient_table


# ============================================================
# PSEUDO-AIC MODEL COMPARISON
# ============================================================

extract_fine_gray_fit <- function(fitted_model, model_name) {
  crr_fit <- fitted_model$crrFit
  
  number_parameters <- length(crr_fit$coef)
  log_pseudo_likelihood <- as.numeric(crr_fit$loglik)
  pseudo_aic <- -2 * log_pseudo_likelihood + 2 * number_parameters
  
  tibble(
    Model = model_name,
    N = nrow(analysis_data),
    Parameters = number_parameters,
    Converged = isTRUE(
      crr_fit$converged
    ),
    Log_Pseudo_Likelihood = log_pseudo_likelihood,
    Pseudo_AIC = pseudo_aic
  )
}

fine_gray_fit_table <- imap_dfr(
  fine_gray_models,
  extract_fine_gray_fit
) %>%
  arrange(Pseudo_AIC) %>%
  mutate(
    Delta_Pseudo_AIC = Pseudo_AIC - min(Pseudo_AIC),
    relative_likelihood = exp(-0.5 * Delta_Pseudo_AIC),
    Pseudo_AIC_Weight = relative_likelihood / sum(relative_likelihood),
    Pseudo_AIC_Rank = row_number()
  ) %>%
  select(-relative_likelihood)

fine_gray_fit_table


# ============================================================
# CROSS-VALIDATION SETTINGS
# ============================================================

prediction_times <- seq(
  from = 10,
  to = 300,
  by = 10
)

prediction_times

set.seed(2026)


# ============================================================
# CROSS-VALIDATED BRIER SCORE AND AUC
# ============================================================

fine_gray_cv_score <- riskRegression::Score(
  object = fine_gray_models,
  formula = Hist(observed_time, event_type) ~ 1,
  data = analysis_data,
  cause = 1,
  times = prediction_times,
  metrics = c("brier", "auc"),
  split.method = "cv5",
  B = 1,
  null.model = TRUE,
  contrasts = FALSE,
  conf.int = FALSE,
  cens.method = "ipcw",
  seed = 2026,
  parallel = "no",
  progress.bar = 3,
  errorhandling = "stop"
)

fine_gray_cv_score
summary(fine_gray_cv_score)


# ============================================================
# BRIER SCORE RESULTS
# ============================================================

cv_brier_table <- fine_gray_cv_score$Brier$score %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_with(tolower)

names(cv_brier_table)

cv_brier_table

if ("score" %in% names(cv_brier_table) &&
    !"brier" %in% names(cv_brier_table)) {
  cv_brier_table <- cv_brier_table %>%
    rename(
      brier = score
    )
}

if ("time" %in% names(cv_brier_table) &&
    !"times" %in% names(cv_brier_table)) {
  cv_brier_table <- cv_brier_table %>%
    rename(
      times = time
    )
}

cv_brier_model_table <- cv_brier_table %>%
  filter(
    !model %in% c(
      "Reference",
      "Null model"
    )
  ) %>%
  arrange(
    times,
    brier
  )

cv_brier_model_table


# ============================================================
# AUC RESULTS
# ============================================================

cv_auc_table <- fine_gray_cv_score$AUC$score %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_with(tolower)

names(cv_auc_table)

cv_auc_table

if ("score" %in% names(cv_auc_table) &&
    !"auc" %in% names(cv_auc_table)) {
  cv_auc_table <- cv_auc_table %>%
    rename(
      auc = score
    )
}

if ("time" %in% names(cv_auc_table) &&
    !"times" %in% names(cv_auc_table)) {
  cv_auc_table <- cv_auc_table %>%
    rename(
      times = time
    )
}

cv_auc_model_table <- cv_auc_table %>%
  filter(
    !model %in% c(
      "Reference",
      "Null model"
    )
  ) %>%
  arrange(
    times,
    desc(auc)
  )

cv_auc_model_table


# ============================================================
# MEAN CROSS-VALIDATION PERFORMANCE
# ============================================================

mean_brier_by_model <- cv_brier_model_table %>%
  group_by(model) %>%
  summarise(
    Mean_CV_Brier = mean(brier, na.rm = TRUE),
    SD_CV_Brier = sd(brier, na.rm = TRUE),
    .groups = "drop"
  )

mean_auc_by_model <- cv_auc_model_table %>%
  group_by(model) %>%
  summarise(
    Mean_CV_AUC = mean(
      auc,
      na.rm = TRUE
    ),
    SD_CV_AUC = sd(
      auc,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ============================================================
# FINAL MODEL COMPARISON
# ============================================================

fine_gray_model_comparison <- fine_gray_fit_table %>%
  rename(model = Model) %>%
  left_join(
    mean_brier_by_model,
    by = "model"
  ) %>%
  left_join(
    mean_auc_by_model,
    by = "model"
  ) %>%
  arrange(
    Mean_CV_Brier,
    desc(Mean_CV_AUC),
    Parameters
  ) %>%
  mutate(
    Prediction_Rank = row_number()
  ) %>%
  rename(
    Model = model
  )

fine_gray_model_comparison

fine_gray_model_comparison_table <-
  fine_gray_model_comparison %>%
  transmute(
    Model,
    Parameters,
    `Pseudo-AIC` = round(Pseudo_AIC, 2),
    `Delta pseudo-AIC` = round(Delta_Pseudo_AIC, 2),
    `Pseudo-AIC weight` = round(Pseudo_AIC_Weight, 3),
    `Mean CV Brier` = round(Mean_CV_Brier, 4),
    `Mean CV AUC` = round(Mean_CV_AUC, 3),
    `Prediction rank` = Prediction_Rank
  )

fine_gray_model_comparison_table


# ============================================================
# BEST MODEL
# ============================================================

best_fine_gray_model_name <-
  fine_gray_model_comparison %>%
  filter(
    is.finite(Mean_CV_Brier)
  ) %>%
  slice_min(
    order_by = Mean_CV_Brier,
    n = 1,
    with_ties = FALSE
  ) %>%
  pull(Model)

best_fine_gray_model_name

best_fine_gray_model <- fine_gray_models[[best_fine_gray_model_name]]

best_fine_gray_coefficients <- fine_gray_coefficient_table %>%
  filter(
    Model == best_fine_gray_model_name
  )

best_fine_gray_coefficients
summary(best_fine_gray_model)


# ============================================================
# BRIER SCORE PLOT
# ============================================================

brier_plot <- ggplot(
  cv_brier_model_table,
  aes(
    x = times,
    y = brier,
    group = model,
    color = model
  )
) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cross-validated Brier Scores",
    subtitle = "Fine-Gray candidate models",
    x = "Time since initial earthquake (days)",
    y = "Brier score"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

brier_plot


# ============================================================
# AUC PLOT
# ============================================================

auc_plot <- ggplot(
  cv_auc_model_table,
  aes(
    x = times,
    y = auc,
    group = model,
    color = model
  )
) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cross-validated Time-dependent AUC",
    subtitle = "Fine-Gray candidate models",
    x = "Time since initial earthquake (days)",
    y = "AUC",
    color = "Model"
  ) +
  theme_bw()

auc_plot


# ============================================================
# CUMULATIVE INCIDENCE
# ============================================================

cumulative_incidence_fit <- cmprsk::cuminc(
  ftime = analysis_data$observed_time,
  fstatus = analysis_data$event_type,
  cencode = 0
)

cumulative_incidence_fit


# ============================================================
# SAVE RESULTS
# ============================================================

fine_gray_all_results <- list(
  data = analysis_data,
  formulas = fine_gray_candidate_formulas,
  fit_status = fine_gray_fit_status,
  fitted_models = fine_gray_models,
  coefficient_results = fine_gray_coefficient_results,
  coefficient_table = fine_gray_coefficient_table,
  fit_statistics = fine_gray_fit_table,
  cross_validation = fine_gray_cv_score,
  brier_results = cv_brier_model_table,
  auc_results = cv_auc_model_table,
  model_comparison = fine_gray_model_comparison,
  model_comparison_table = fine_gray_model_comparison_table,
  best_model_name = best_fine_gray_model_name,
  best_model = best_fine_gray_model,
  best_model_coefficients = best_fine_gray_coefficients,
  cumulative_incidence = cumulative_incidence_fit
)

saveRDS(
  fine_gray_all_results,
  file = "fine_gray_14_candidate_models.rds"
)

write.csv(
  fine_gray_coefficient_table,
  file = "fine_gray_all_coefficient_results.csv",
  row.names = FALSE
)

write.csv(
  fine_gray_model_comparison_table,
  file = "fine_gray_model_comparison.csv",
  row.names = FALSE
)

write.csv(
  best_fine_gray_coefficients,
  file = "fine_gray_best_model_coefficients.csv",
  row.names = FALSE
)


# ============================================================
# FINAL OUTPUTS
# ============================================================

fine_gray_fit_status

fine_gray_model_comparison_table

best_fine_gray_model_name

best_fine_gray_coefficients


# ============================================================
# COMBINED PERFORMANCE PLOT
# ============================================================

combined_plot <- brier_plot + auc_plot + plot_layout(ncol = 2)

combined_plot