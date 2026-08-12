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
    restore_expen = restore_expen / 1e6
  ) %>%
  mutate(
    log_gdp = log(gdp2015),
    log_depth = log(depth),
    event = ifelse(first_event == "recovered", 1, 0)
  )

# ============================================================
# CANDIDATE COX MODELS
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
  M10 = Surv(recovery_days, event) ~ magnitude + depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen
)

# ============================================================
# FIT COX MODELS
# ============================================================

cox_models <- map(
  candidate_formulas,
  ~ coxph(
    .x,
    data = data_table,
    x = TRUE,
    y = TRUE
  )
)

# ============================================================
# MODEL COMPARISON
# ============================================================

cox_comparison <- imap_dfr(
  cox_models,
  ~ tibble(
    model = .y,
    n = .x$n,
    events = .x$nevent,
    parameters = length(coef(.x)),
    log_likelihood = as.numeric(logLik(.x)),
    AIC = AIC(.x),
    concordance = summary(.x)$concordance[1]
  )
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC))
  )

cox_comparison

# ============================================================
# COX COEFFICIENT RESULTS
# ============================================================

cox_results <- imap_dfr(
  cox_models,
  ~ tidy(
    .x,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      model = .y,
      .before = 1
    )
) %>%
  dplyr::select(
    model,
    term,
    hazard_ratio = estimate,
    conf.low,
    conf.high,
    p.value
  )

cox_results

# ============================================================
# PROPORTIONAL HAZARDS TEST
# ============================================================

ph_tests <- imap_dfr(
  cox_models,
  ~ {
    test <- cox.zph(.x)
    
    as.data.frame(test$table) %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        model = .y,
        .before = 1
      )
  }
)

ph_tests

# ============================================================
# CANDIDATE MODELS WITH LOG-TRANSFORMED DEPTH
# ============================================================

candidate_formulas_log <- list(
  M1 = Surv(recovery_days, event) ~ magnitude + log_depth,
  M2 = Surv(recovery_days, event) ~ magnitude + log_depth + mean_pga,
  M3 = Surv(recovery_days, event) ~ magnitude + log_depth + max_pga,
  M4 = Surv(recovery_days, event) ~ magnitude + log_depth + ntl_decline_ratio,
  M5 = Surv(recovery_days, event) ~ magnitude + log_depth + gdp2015 + population + dis_expen,
  M6 = Surv(recovery_days, event) ~ magnitude + log_depth + mean_pga + ntl_decline_ratio,
  M7 = Surv(recovery_days, event) ~ magnitude + log_depth + max_pga + ntl_decline_ratio,
  M8 = Surv(recovery_days, event) ~ magnitude + log_depth + ntl_decline_ratio + gdp2015 + population + dis_expen,
  M9 = Surv(recovery_days, event) ~ magnitude + log_depth + mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen,
  M10 = Surv(recovery_days, event) ~ magnitude + log_depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen
)

# ============================================================
# FIT LOG-DEPTH COX MODELS
# ============================================================

cox_models2 <- map(
  candidate_formulas_log,
  ~ coxph(
    .x,
    data = data_table,
    x = TRUE,
    y = TRUE
  )
)

# ============================================================
# LOG-DEPTH MODEL COMPARISON
# ============================================================

cox_comparison2 <- imap_dfr(
  cox_models2,
  ~ tibble(
    model = .y,
    n = .x$n,
    events = .x$nevent,
    parameters = length(coef(.x)),
    log_likelihood = as.numeric(logLik(.x)),
    AIC = AIC(.x),
    concordance = summary(.x)$concordance[1]
  )
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC))
  )

cox_comparison2

# ============================================================
# LOG-DEPTH COEFFICIENT RESULTS
# ============================================================

cox_results2 <- imap_dfr(
  cox_models2,
  ~ tidy(
    .x,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      model = .y,
      .before = 1
    )
) %>%
  dplyr::select(
    model,
    term,
    hazard_ratio = estimate,
    conf.low,
    conf.high,
    p.value
  )

cox_results2

# ============================================================
# LOG-DEPTH PROPORTIONAL HAZARDS TEST
# ============================================================

ph_tests2 <- imap_dfr(
  cox_models2,
  ~ {
    test <- cox.zph(.x)
    
    as.data.frame(test$table) %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        model = .y,
        .before = 1
      )
  }
)

ph_tests2

# ============================================================
# M5 RESULTS TABLE
# ============================================================

cox_M5 <- cox_results %>%
  filter(model == "M5") %>%
  mutate(
    Predictor = recode(
      term,
      magnitude = "Earthquake magnitude",
      depth = "Earthquake depth",
      gdp2015 = "GDP",
      population = "Population",
      dis_expen = "Disaster expenditure"
    ),
    `Hazard ratio` = sprintf("%.3f", hazard_ratio),
    `95%% CI` = sprintf("%.3f–%.3f", conf.low, conf.high),
    `p-value` = if_else(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    Predictor,
    `Hazard ratio`,
    `95%% CI`,
    `p-value`
  )

cox_M5