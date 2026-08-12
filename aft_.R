# ============================================================
# PACKAGES
# ============================================================

library(survival)
library(dplyr)
library(purrr)
library(broom)

# ============================================================
# DATA PREPARATION
# ============================================================

pga_recovery_econ_clean1 = read_csv(
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv'
)

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
# FIT AFT MODELS
# ============================================================

aft_distributions <- c(
  Weibull = "weibull",
  Lognormal = "lognormal",
  Loglogistic = "loglogistic"
)

aft_models <- map(
  aft_distributions,
  function(dist_name) {
    map(
      candidate_formulas,
      ~ survreg(
        .x,
        data = data_table,
        dist = dist_name,
        x = TRUE,
        y = TRUE
      )
    )
  }
)

# ============================================================
# AIC MODEL COMPARISON
# ============================================================

aft_comparison <- imap_dfr(
  aft_models,
  function(model_set, distribution_name) {
    imap_dfr(
      model_set,
      function(model_object, model_name) {
        tibble(
          distribution = distribution_name,
          model = model_name,
          n = model_object$n,
          parameters = length(coef(model_object)) + 1,
          log_likelihood = as.numeric(logLik(model_object)),
          AIC = AIC(model_object)
        )
      }
    )
  }
) %>%
  group_by(distribution) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC))
  ) %>%
  ungroup() %>%
  arrange(distribution, AIC)

print(aft_comparison, n = 42)

# ============================================================
# AFT COEFFICIENT RESULTS
# ============================================================

aft_results <- imap_dfr(
  aft_models,
  function(model_set, distribution_name) {
    imap_dfr(
      model_set,
      function(model_object, model_name) {
        tidy(
          model_object,
          conf.int = TRUE
        ) %>%
          filter(term != "Log(scale)") %>%
          mutate(
            distribution = distribution_name,
            model = model_name,
            time_ratio = exp(estimate),
            time_ratio_low = exp(conf.low),
            time_ratio_high = exp(conf.high)
          )
      }
    )
  }
) %>%
  dplyr::select(
    distribution,
    model,
    term,
    estimate,
    time_ratio,
    time_ratio_low,
    time_ratio_high,
    p.value
  )

aft_results

# ============================================================
# BEST AFT MODEL
# ============================================================

best_aft_model <- aft_comparison %>%
  slice_min(
    order_by = AIC,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::select(
    distribution,
    model,
    parameters,
    log_likelihood,
    AIC
  )

best_aft_model

best_aft_results <- aft_results %>%
  semi_join(
    best_aft_model,
    by = c("distribution", "model")
  )

best_aft_results

# ============================================================
# BEST MODEL RESULTS TABLE
# ============================================================

best_aft_table <- aft_results %>%
  semi_join(
    best_aft_model,
    by = c("distribution", "model")
  ) %>%
  mutate(
    Predictor = recode(
      term,
      magnitude = "Earthquake magnitude",
      depth = "Earthquake depth",
      mean_pga = "Mean PGA",
      max_pga = "Maximum PGA",
      ntl_decline_ratio = "Nighttime light decline ratio",
      gdp2015 = "GDP",
      log_gdp = "Log GDP",
      population = "Population",
      dis_expen = "Disaster expenditure"
    ),
    `Time ratio` = sprintf("%.3f", time_ratio),
    `95% CI` = sprintf(
      "%.3f--%.3f",
      time_ratio_low,
      time_ratio_high
    ),
    `p-value` = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    Distribution = distribution,
    Model = model,
    Predictor,
    `Time ratio`,
    `95% CI`,
    `p-value`
  )

best_aft_table

# ============================================================
# LOGNORMAL M6 MODEL
# ============================================================

lognormal_m6_aft <- survreg(
  Surv(recovery_days, event) ~
    magnitude +
    depth +
    mean_pga +
    ntl_decline_ratio,
  data = bayes_data,
  dist = "lognormal"
)

summary(lognormal_m6_aft)