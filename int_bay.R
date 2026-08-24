# ============================================================
# BAYESIAN INTERVAL-CENSORED AFT + LOO/ELPD
# ============================================================

library(brms)
library(purrr)
library(tibble)
library(loo)
library(dplyr)


# ============================================================
# STANDARDIZE CONTINUOUS PREDICTORS BUT KEEP ORIGINAL NAMES
# ============================================================

bayes_aft_data <- interval_analysis_data %>%
  mutate(
    earthquake_id = factor(eq_id),
    censor_type = case_when(
      is.finite(recovery_right) ~ "interval",
      is.infinite(recovery_right) ~ "right",
      TRUE ~ NA_character_
    ),
    recovery_y = if_else(recovery_left <= 0, 0.01, recovery_left),
    recovery_y2 = if_else(is.finite(recovery_right), recovery_right, recovery_y),
    
    magnitude = as.numeric(scale(magnitude)),
    depth = as.numeric(scale(depth)),
    mean_pga = as.numeric(scale(mean_pga)),
    max_pga = as.numeric(scale(max_pga)),
    ntl_decline_ratio = as.numeric(scale(ntl_decline_ratio)),
    gdp2015 = as.numeric(scale(gdp2015)),
    log_gdp = as.numeric(scale(log_gdp)),
    population = as.numeric(scale(population)),
    dis_expen = as.numeric(scale(dis_expen))
  ) %>%
  filter(
    !is.na(censor_type),
    is.finite(recovery_y),
    recovery_y > 0
  )

summary(
  bayes_aft_data %>%
    select(
      magnitude,
      depth,
      mean_pga,
      max_pga,
      ntl_decline_ratio,
      gdp2015,
      log_gdp,
      population,
      dis_expen
    )
)

bayes_aft_data$recovery_y2
# ------------------------------------------------------------
# COMMON PRIORS
# ------------------------------------------------------------

common_priors <- c(
  prior(normal(0, 0.3), class = "b"),
  prior(normal(log(90), 1), class = "Intercept")
)

# Earthquake-level random-effect prior
random_effect_prior <- prior(
  student_t(3, 0, 0.3),
  class = "sd",
  group = "earthquake_id"
)

# Lognormal priors
lognormal_fixed_priors <- c(
  common_priors,
  prior(student_t(3, 0, 0.35), class = "sigma")
)

lognormal_random_priors <- c(
  lognormal_fixed_priors,
  random_effect_prior
)

# Weibull priors
weibull_fixed_priors <- c(
  common_priors,
  prior(lognormal(log(1.5), 0.25), class = "shape")
)

weibull_random_priors <- c(
  weibull_fixed_priors,
  random_effect_prior
)

# ------------------------------------------------------------
# FORMULAS
# ------------------------------------------------------------

bayesian_rhs <- list(
  M1 = "magnitude + depth",
  M2 = "magnitude + depth + mean_pga",
  M3 = "magnitude + depth + max_pga",
  M4 = "magnitude + depth + ntl_decline_ratio",
  M5 = "magnitude + depth + gdp2015 + population + dis_expen",
  M6 = "magnitude + depth + mean_pga + ntl_decline_ratio",
  M7 = "magnitude + depth + max_pga + ntl_decline_ratio",
  M8 = "magnitude + depth + ntl_decline_ratio + gdp2015 + population + dis_expen",
  M9 = "magnitude + depth + mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen",
  M10 = "magnitude + depth + max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen",
  M11 = "magnitude + depth + log_gdp + population + dis_expen",
  M12 = "magnitude + depth + ntl_decline_ratio + log_gdp + population + dis_expen",
  M13 = "magnitude + depth + mean_pga + ntl_decline_ratio + log_gdp + population + dis_expen",
  M14 = "magnitude + depth + max_pga + ntl_decline_ratio + log_gdp + population + dis_expen"
)

fixed_formulas <- map(
  bayesian_rhs,
  ~ bf(as.formula(
    paste(
      "recovery_y | cens(censor_type, recovery_y2) ~",
      .x
    )
  ))
)

random_formulas <- map(
  bayesian_rhs,
  ~ bf(as.formula(
    paste(
      "recovery_y | cens(censor_type, recovery_y2) ~",
      .x,
      "+ (1 | earthquake_id)"
    )
  ))
)

# ------------------------------------------------------------
# FIT FUNCTION
# ------------------------------------------------------------

fit_bayes_aft <- function(formula, model_name, family_name, random_effect, data) {
  if (family_name == "lognormal") {
    family_use <- lognormal()
    priors_use <- if (random_effect) lognormal_random_priors else lognormal_fixed_priors
  } else if (family_name == "weibull") {
    family_use <- weibull()
    priors_use <- if (random_effect) weibull_random_priors else weibull_fixed_priors
  } else {
    stop("family_name must be 'lognormal' or 'weibull'")
  }
  
  tryCatch(
    {
      fit <- brm(
        formula = formula,
        data = data,
        family = family_use,
        prior = priors_use,
        chains = 4,
        iter = 8000,
        warmup = 4000,
        cores = min(4, parallel::detectCores()),
        seed = 2026,
        control = list(
          adapt_delta = 0.99,
          max_treedepth = 15
        ),
        save_pars = save_pars(all = TRUE),
        refresh = 0
      )
      
      list(
        Model = model_name,
        Family = family_name,
        Structure = ifelse(random_effect, "Random earthquake", "Fixed"),
        Fit = fit,
        Successful = TRUE,
        Error = NA_character_
      )
    },
    error = function(e) {
      list(
        Model = model_name,
        Family = family_name,
        Structure = ifelse(random_effect, "Random earthquake", "Fixed"),
        Fit = NULL,
        Successful = FALSE,
        Error = conditionMessage(e)
      )
    }
  )
}

# ------------------------------------------------------------
# FIT ALL LOGNORMAL MODELS
# ------------------------------------------------------------

lognormal_fixed_results <- imap(
  fixed_formulas,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "lognormal",
    random_effect = FALSE,
    data = bayes_aft_data
  )
)

lognormal_random_results <- imap(
  random_formulas,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "lognormal",
    random_effect = TRUE,
    data = bayes_aft_data
  )
)

# ------------------------------------------------------------
# FIT ALL WEIBULL MODELS
# ------------------------------------------------------------

weibull_fixed_results <- imap(
  fixed_formulas,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "weibull",
    random_effect = FALSE,
    data = bayes_aft_data
  )
)

weibull_random_results <- imap(
  random_formulas,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "weibull",
    random_effect = TRUE,
    data = bayes_aft_data
  )
)

# ------------------------------------------------------------
# FIT STATUS
# ------------------------------------------------------------

all_fit_results <- c(
  lognormal_fixed_results,
  lognormal_random_results,
  weibull_fixed_results,
  weibull_random_results
)

fit_status <- map_dfr(
  all_fit_results,
  ~ tibble(
    Model = .x$Model,
    Family = .x$Family,
    Structure = .x$Structure,
    Successful = .x$Successful,
    Error = .x$Error
  )
)

print(fit_status, n = Inf, width = Inf)

# ------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
# ------------------------------------------------------------

extract_convergence <- function(x) {
  if (!x$Successful) {
    return(tibble(
      Model = x$Model,
      Family = x$Family,
      Structure = x$Structure,
      Max_Rhat = NA_real_,
      Min_Bulk_ESS = NA_real_,
      Min_Tail_ESS = NA_real_
    ))
  }
  
  s <- posterior::summarise_draws(
    posterior::as_draws_array(x$Fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  
  tibble(
    Model = x$Model,
    Family = x$Family,
    Structure = x$Structure,
    Max_Rhat = max(s$rhat, na.rm = TRUE),
    Min_Bulk_ESS = min(s$ess_bulk, na.rm = TRUE),
    Min_Tail_ESS = min(s$ess_tail, na.rm = TRUE)
  )
}

convergence_table <- map_dfr(
  all_fit_results,
  extract_convergence
)

print(convergence_table, n = Inf, width = Inf)

# ------------------------------------------------------------
# SAFE LOO:
# 1. standard PSIS-LOO
# 2. if problematic k -> moment matching
# 3. if still problematic -> reloo
# ------------------------------------------------------------

compute_loo_robust <- function(fit, model_name, family_name, structure,
                               k_threshold = 0.7) {
  tryCatch(
    {
      loo_initial <- brms::loo(fit)
      k_initial <- loo::pareto_k_values(loo_initial)
      max_k_initial <- max(k_initial, na.rm = TRUE)
      
      loo_final <- loo_initial
      method_used <- "PSIS-LOO"
      
      if (any(k_initial > k_threshold, na.rm = TRUE)) {
        loo_mm <- tryCatch(
          brms::loo(
            fit,
            moment_match = TRUE
          ),
          error = function(e) NULL
        )
        
        if (!is.null(loo_mm)) {
          k_mm <- loo::pareto_k_values(loo_mm)
          loo_final <- loo_mm
          method_used <- "Moment matching"
          
          if (any(k_mm > k_threshold, na.rm = TRUE)) {
            loo_reloo <- tryCatch(
              brms::loo(
                fit,
                reloo = TRUE
              ),
              error = function(e) NULL
            )
            
            if (!is.null(loo_reloo)) {
              loo_final <- loo_reloo
              method_used <- "ReLOO"
            }
          }
        } else {
          loo_reloo <- tryCatch(
            brms::loo(
              fit,
              reloo = TRUE
            ),
            error = function(e) NULL
          )
          
          if (!is.null(loo_reloo)) {
            loo_final <- loo_reloo
            method_used <- "ReLOO"
          }
        }
      }
      
      final_k <- tryCatch(
        loo::pareto_k_values(loo_final),
        error = function(e) rep(NA_real_, nrow(bayes_aft_data))
      )
      
      est <- loo_final$estimates
      
      list(
        Model = model_name,
        Family = family_name,
        Structure = structure,
        LOO = loo_final,
        Method = method_used,
        ELPD_LOO = est["elpd_loo", "Estimate"],
        SE_ELPD = est["elpd_loo", "SE"],
        P_LOO = est["p_loo", "Estimate"],
        LOOIC = est["looic", "Estimate"],
        Max_Pareto_k = suppressWarnings(max(final_k, na.rm = TRUE)),
        N_k_over_0_7 = sum(final_k > 0.7, na.rm = TRUE),
        N_k_over_1 = sum(final_k > 1, na.rm = TRUE),
        Successful = TRUE,
        Error = NA_character_
      )
    },
    error = function(e) {
      list(
        Model = model_name,
        Family = family_name,
        Structure = structure,
        LOO = NULL,
        Method = NA_character_,
        ELPD_LOO = NA_real_,
        SE_ELPD = NA_real_,
        P_LOO = NA_real_,
        LOOIC = NA_real_,
        Max_Pareto_k = NA_real_,
        N_k_over_0_7 = NA_integer_,
        N_k_over_1 = NA_integer_,
        Successful = FALSE,
        Error = conditionMessage(e)
      )
    }
  )
}

# ------------------------------------------------------------
# RUN ROBUST LOO FOR ALL SUCCESSFUL MODELS
# ------------------------------------------------------------

loo_results <- map(
  all_fit_results,
  function(x) {
    if (!x$Successful) {
      return(list(
        Model = x$Model,
        Family = x$Family,
        Structure = x$Structure,
        LOO = NULL,
        Method = NA_character_,
        ELPD_LOO = NA_real_,
        SE_ELPD = NA_real_,
        P_LOO = NA_real_,
        LOOIC = NA_real_,
        Max_Pareto_k = NA_real_,
        N_k_over_0_7 = NA_integer_,
        N_k_over_1 = NA_integer_,
        Successful = FALSE,
        Error = x$Error
      ))
    }
    
    compute_loo_robust(
      fit = x$Fit,
      model_name = x$Model,
      family_name = x$Family,
      structure = x$Structure
    )
  }
)

loo_summary <- map_dfr(
  loo_results,
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
    N_k_over_1 = .x$N_k_over_1,
    Successful = .x$Successful,
    Error = .x$Error
  )
) %>%
  arrange(desc(ELPD_LOO))

print(loo_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# BUILD NAMED LIST OF FINAL LOO OBJECTS
# ------------------------------------------------------------

loo_object_list <- list()

for (x in loo_results) {
  if (x$Successful && !is.null(x$LOO)) {
    nm <- paste(
      x$Model,
      x$Family,
      ifelse(x$Structure == "Fixed", "fixed", "random"),
      sep = "_"
    )
    loo_object_list[[nm]] <- x$LOO
  }
}

# ------------------------------------------------------------
# GLOBAL LOO COMPARISON
# ------------------------------------------------------------

global_loo_compare <- loo::loo_compare(
  loo_object_list
)

global_loo_compare

global_loo_table <- as.data.frame(
  global_loo_compare
) %>%
  rownames_to_column("Model_specification") %>%
  as_tibble() %>%
  rename(
    Delta_ELPD = elpd_diff,
    SE_Delta_ELPD = se_diff
  )

print(global_loo_table, n = Inf, width = Inf)

# ------------------------------------------------------------
# JOIN ABSOLUTE ELPD + DELTA ELPD
# ------------------------------------------------------------

loo_summary <- loo_summary %>%
  mutate(
    Model_specification = paste(
      Model,
      Family,
      ifelse(Structure == "Fixed", "fixed", "random"),
      sep = "_"
    )
  )

bayesian_model_comparison <- global_loo_table %>%
  left_join(
    loo_summary,
    by = "Model_specification"
  ) %>%
  arrange(desc(ELPD_LOO)) %>%
  mutate(ELPD_Rank = row_number())

print(
  bayesian_model_comparison,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# BEST MODEL
# ------------------------------------------------------------

best_bayesian_model <- bayesian_model_comparison %>%
  filter(is.finite(ELPD_LOO)) %>%
  slice_max(
    ELPD_LOO,
    n = 1,
    with_ties = FALSE
  )

best_bayesian_model

# ------------------------------------------------------------
# FIXED VS RANDOM WITHIN EACH MODEL/FAMILY
# ------------------------------------------------------------

fixed_random_comparison <- loo_summary %>%
  select(
    Model,
    Family,
    Structure,
    ELPD_LOO,
    SE_ELPD,
    LOOIC,
    LOO_method,
    Max_Pareto_k
  ) %>%
  tidyr::pivot_wider(
    names_from = Structure,
    values_from = c(
      ELPD_LOO,
      SE_ELPD,
      LOOIC,
      LOO_method,
      Max_Pareto_k
    )
  ) %>%
  mutate(
    Random_minus_Fixed_ELPD =
      `ELPD_LOO_Random earthquake` - ELPD_LOO_Fixed,
    Preferred_structure = case_when(
      Random_minus_Fixed_ELPD > 0 ~ "Random earthquake",
      Random_minus_Fixed_ELPD < 0 ~ "Fixed",
      TRUE ~ "Equal"
    )
  )

print(
  fixed_random_comparison,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# FINAL REPORTING TABLE
# ------------------------------------------------------------

bayesian_loo_table <- bayesian_model_comparison %>%
  transmute(
    Rank = ELPD_Rank,
    Model,
    Family,
    Structure,
    `LOO method` = LOO_method,
    `ELPD-LOO` = round(ELPD_LOO, 2),
    `SE ELPD` = round(SE_ELPD, 2),
    `Delta ELPD` = round(Delta_ELPD, 2),
    `SE Delta ELPD` = round(SE_Delta_ELPD, 2),
    `p-LOO` = round(P_LOO, 2),
    LOOIC = round(LOOIC, 2),
    `Max Pareto k` = round(Max_Pareto_k, 2)
  )

print(
  bayesian_loo_table,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

bayesian_aft_results <- list(
  data = bayes_aft_data,
  fixed_formulas = fixed_formulas,
  random_formulas = random_formulas,
  fit_status = fit_status,
  convergence = convergence_table,
  lognormal_fixed = lognormal_fixed_results,
  lognormal_random = lognormal_random_results,
  weibull_fixed = weibull_fixed_results,
  weibull_random = weibull_random_results,
  loo_results = loo_results,
  loo_summary = loo_summary,
  global_loo_comparison = global_loo_table,
  final_comparison = bayesian_model_comparison,
  fixed_random_comparison = fixed_random_comparison,
  best_model = best_bayesian_model,
  reporting_table = bayesian_loo_table
)

saveRDS(
  bayesian_aft_results,
  file = "bayesian_interval_aft_loo_results.rds"
)


# ============================================================
# BAYESIAN INTERVAL-CENSORED AFT WITHOUT MAGNITUDE/DEPTH
# REUSE EXISTING FUNCTIONS
# ============================================================


bayes_aft_data_noeq <- interval_analysis_data %>%
  mutate(
    earthquake_id = factor(eq_id),
    censor_type = case_when(
      is.finite(recovery_right) ~ "interval",
      is.infinite(recovery_right) ~ "right",
      TRUE ~ NA_character_
    ),
    recovery_y = if_else(recovery_left <= 0, 0.01, recovery_left),
    recovery_y2 = if_else(is.finite(recovery_right), recovery_right, recovery_y),
    mean_pga = as.numeric(scale(mean_pga)),
    max_pga = as.numeric(scale(max_pga)),
    ntl_decline_ratio = as.numeric(scale(ntl_decline_ratio)),
    gdp2015 = as.numeric(scale(gdp2015)),
    log_gdp = as.numeric(scale(log_gdp)),
    population = as.numeric(scale(population)),
    dis_expen = as.numeric(scale(dis_expen))
  ) %>%
  filter(!is.na(censor_type), is.finite(recovery_y), recovery_y > 0)

bayes_aft_data_noeq %>% count(censor_type)

summary(
  bayes_aft_data_noeq %>%
    select(mean_pga, max_pga, ntl_decline_ratio, gdp2015,
           log_gdp, population, dis_expen)
)

bayesian_rhs_noeq <- list(
  M1 = "1",
  M2 = "mean_pga",
  M3 = "max_pga",
  M4 = "ntl_decline_ratio",
  M5 = "gdp2015 + population + dis_expen",
  M6 = "mean_pga + ntl_decline_ratio",
  M7 = "max_pga + ntl_decline_ratio",
  M8 = "ntl_decline_ratio + gdp2015 + population + dis_expen",
  M9 = "mean_pga + ntl_decline_ratio + gdp2015 + population + dis_expen",
  M10 = "max_pga + ntl_decline_ratio + gdp2015 + population + dis_expen",
  M11 = "log_gdp + population + dis_expen",
  M12 = "ntl_decline_ratio + log_gdp + population + dis_expen",
  M13 = "mean_pga + ntl_decline_ratio + log_gdp + population + dis_expen",
  M14 = "max_pga + ntl_decline_ratio + log_gdp + population + dis_expen"
)

fixed_formulas_noeq <- map(
  bayesian_rhs_noeq,
  ~ bf(as.formula(
    paste("recovery_y | cens(censor_type, recovery_y2) ~", .x)
  ))
)

random_formulas_noeq <- map(
  bayesian_rhs_noeq,
  ~ bf(as.formula(
    paste(
      "recovery_y | cens(censor_type, recovery_y2) ~",
      .x,
      "+ (1 | earthquake_id)"
    )
  ))
)

lognormal_fixed_results_noeq <- imap(
  fixed_formulas_noeq,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "lognormal",
    random_effect = FALSE,
    data = bayes_aft_data_noeq
  )
)

lognormal_random_results_noeq <- imap(
  random_formulas_noeq,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "lognormal",
    random_effect = TRUE,
    data = bayes_aft_data_noeq
  )
)

weibull_fixed_results_noeq <- imap(
  fixed_formulas_noeq,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "weibull",
    random_effect = FALSE,
    data = bayes_aft_data_noeq
  )
)

weibull_random_results_noeq <- imap(
  random_formulas_noeq,
  ~ fit_bayes_aft(
    formula = .x,
    model_name = .y,
    family_name = "weibull",
    random_effect = TRUE,
    data = bayes_aft_data_noeq
  )
)

all_fit_results_noeq <- c(
  lognormal_fixed_results_noeq,
  lognormal_random_results_noeq,
  weibull_fixed_results_noeq,
  weibull_random_results_noeq
)

fit_status_noeq <- map_dfr(
  all_fit_results_noeq,
  ~ tibble(
    Model = .x$Model,
    Family = .x$Family,
    Structure = .x$Structure,
    Successful = .x$Successful,
    Error = .x$Error
  )
)

print(fit_status_noeq, n = Inf, width = Inf)

convergence_table_noeq <- map_dfr(
  all_fit_results_noeq,
  extract_convergence
)

print(convergence_table_noeq, n = Inf, width = Inf)

convergence_check_noeq <- convergence_table_noeq %>%
  mutate(
    Rhat_OK = Max_Rhat < 1.01,
    Bulk_ESS_OK = Min_Bulk_ESS >= 400,
    Tail_ESS_OK = Min_Tail_ESS >= 400,
    Convergence_OK = Rhat_OK & Bulk_ESS_OK & Tail_ESS_OK
  )

print(convergence_check_noeq, n = Inf, width = Inf)

loo_results_noeq <- map(
  all_fit_results_noeq,
  function(x) {
    if (!x$Successful) {
      return(list(
        Model = x$Model,
        Family = x$Family,
        Structure = x$Structure,
        LOO = NULL,
        Method = NA_character_,
        ELPD_LOO = NA_real_,
        SE_ELPD = NA_real_,
        P_LOO = NA_real_,
        LOOIC = NA_real_,
        Max_Pareto_k = NA_real_,
        N_k_over_0_7 = NA_integer_,
        N_k_over_1 = NA_integer_,
        Successful = FALSE,
        Error = x$Error
      ))
    }
    
    compute_loo_robust(
      fit = x$Fit,
      model_name = x$Model,
      family_name = x$Family,
      structure = x$Structure
    )
  }
)

loo_summary_noeq <- map_dfr(
  loo_results_noeq,
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
    N_k_over_1 = .x$N_k_over_1,
    Successful = .x$Successful,
    Error = .x$Error
  )
) %>%
  arrange(desc(ELPD_LOO))

print(loo_summary_noeq, n = Inf, width = Inf)

loo_summary_noeq %>%
  select(Model, Family, Structure, LOO_method,
         Max_Pareto_k, N_k_over_0_7, N_k_over_1) %>%
  print(n = Inf, width = Inf)

loo_object_list_noeq <- list()

for (x in loo_results_noeq) {
  if (x$Successful && !is.null(x$LOO)) {
    model_id <- paste(
      x$Model,
      x$Family,
      ifelse(x$Structure == "Fixed", "fixed", "random"),
      "noeq",
      sep = "_"
    )
    loo_object_list_noeq[[model_id]] <- x$LOO
  }
}

global_loo_compare_noeq <- loo::loo_compare(
  loo_object_list_noeq
)

global_loo_compare_noeq

global_loo_table_noeq <- as.data.frame(
  global_loo_compare_noeq
) %>%
  rownames_to_column("Model_specification") %>%
  as_tibble() %>%
  rename(
    Delta_ELPD = elpd_diff,
    SE_Delta_ELPD = se_diff
  )

print(global_loo_table_noeq, n = Inf, width = Inf)

loo_summary_noeq <- loo_summary_noeq %>%
  mutate(
    Model_specification = paste(
      Model,
      Family,
      ifelse(Structure == "Fixed", "fixed", "random"),
      "noeq",
      sep = "_"
    )
  )

bayesian_model_comparison_noeq <- global_loo_table_noeq %>%
  left_join(
    loo_summary_noeq,
    by = "Model_specification"
  ) %>%
  arrange(desc(ELPD_LOO)) %>%
  mutate(ELPD_Rank = row_number())

print(bayesian_model_comparison_noeq, n = Inf, width = Inf)

best_bayesian_model_noeq <- bayesian_model_comparison_noeq %>%
  filter(is.finite(ELPD_LOO)) %>%
  slice_max(ELPD_LOO, n = 1, with_ties = FALSE)

best_bayesian_model_noeq

best_fixed_model_noeq <- loo_summary_noeq %>%
  filter(Structure == "Fixed", is.finite(ELPD_LOO)) %>%
  slice_max(ELPD_LOO, n = 1, with_ties = FALSE)

best_fixed_model_noeq

best_random_model_noeq <- loo_summary_noeq %>%
  filter(Structure == "Random earthquake", is.finite(ELPD_LOO)) %>%
  slice_max(ELPD_LOO, n = 1, with_ties = FALSE)

best_random_model_noeq

best_lognormal_model_noeq <- loo_summary_noeq %>%
  filter(Family == "lognormal", is.finite(ELPD_LOO)) %>%
  slice_max(ELPD_LOO, n = 1, with_ties = FALSE)

best_lognormal_model_noeq

best_weibull_model_noeq <- loo_summary_noeq %>%
  filter(Family == "weibull", is.finite(ELPD_LOO)) %>%
  slice_max(ELPD_LOO, n = 1, with_ties = FALSE)

best_weibull_model_noeq

fixed_random_comparison_noeq <- loo_summary_noeq %>%
  select(
    Model,
    Family,
    Structure,
    ELPD_LOO,
    SE_ELPD,
    LOOIC,
    LOO_method,
    Max_Pareto_k
  ) %>%
  pivot_wider(
    names_from = Structure,
    values_from = c(
      ELPD_LOO,
      SE_ELPD,
      LOOIC,
      LOO_method,
      Max_Pareto_k
    )
  ) %>%
  mutate(
    Random_minus_Fixed_ELPD =
      `ELPD_LOO_Random earthquake` - ELPD_LOO_Fixed,
    Preferred_structure = case_when(
      Random_minus_Fixed_ELPD > 0 ~ "Random earthquake",
      Random_minus_Fixed_ELPD < 0 ~ "Fixed",
      TRUE ~ "Equal"
    )
  ) %>%
  arrange(desc(Random_minus_Fixed_ELPD))

print(fixed_random_comparison_noeq, n = Inf, width = Inf)

bayesian_loo_table_noeq <- bayesian_model_comparison_noeq %>%
  transmute(
    Rank = ELPD_Rank,
    Model,
    Family,
    Structure,
    `LOO method` = LOO_method,
    `ELPD-LOO` = round(ELPD_LOO, 2),
    `SE ELPD` = round(SE_ELPD, 2),
    `Delta ELPD` = round(Delta_ELPD, 2),
    `SE Delta ELPD` = round(SE_Delta_ELPD, 2),
    `p-LOO` = round(P_LOO, 2),
    LOOIC = round(LOOIC, 2),
    `Max Pareto k` = round(Max_Pareto_k, 2)
  )

print(bayesian_loo_table_noeq, n = Inf, width = Inf)

bayesian_aft_results_noeq <- list(
  data = bayes_aft_data_noeq,
  rhs = bayesian_rhs_noeq,
  fixed_formulas = fixed_formulas_noeq,
  random_formulas = random_formulas_noeq,
  fit_status = fit_status_noeq,
  convergence = convergence_table_noeq,
  convergence_check = convergence_check_noeq,
  lognormal_fixed = lognormal_fixed_results_noeq,
  lognormal_random = lognormal_random_results_noeq,
  weibull_fixed = weibull_fixed_results_noeq,
  weibull_random = weibull_random_results_noeq,
  all_fits = all_fit_results_noeq,
  loo_results = loo_results_noeq,
  loo_summary = loo_summary_noeq,
  global_loo_comparison = global_loo_table_noeq,
  final_comparison = bayesian_model_comparison_noeq,
  fixed_random_comparison = fixed_random_comparison_noeq,
  best_overall_model = best_bayesian_model_noeq,
  best_fixed_model = best_fixed_model_noeq,
  best_random_model = best_random_model_noeq,
  best_lognormal_model = best_lognormal_model_noeq,
  best_weibull_model = best_weibull_model_noeq,
  reporting_table = bayesian_loo_table_noeq
)

saveRDS(
  bayesian_aft_results_noeq,
  file = "bayesian_interval_aft_no_magnitude_depth_results.rds"
)


# ============================================================
# BAYESIAN MAE/RMSE FOR BOTH ANALYSES
# ============================================================
# ============================================================
# COMMON MAE/RMSE FUNCTION
# ============================================================

calculate_bayes_mae_rmse <- function(x, data) {
  if (!x$Successful || is.null(x$Fit)) {
    return(tibble(
      Model = x$Model,
      Family = x$Family,
      Structure = x$Structure,
      N = NA_integer_,
      MAE = NA_real_,
      RMSE = NA_real_,
      Median_Error = NA_real_,
      Successful = FALSE,
      Error = x$Error
    ))
  }
  
  tryCatch({
    posterior_pred <- posterior_predict(
      x$Fit,
      newdata = data,
      re_formula = NULL
    )
    
    predicted_median <- apply(posterior_pred, 2, median, na.rm = TRUE)
    
    interval_error <- case_when(
      is.finite(data$recovery_right) & predicted_median < data$recovery_left ~
        data$recovery_left - predicted_median,
      is.finite(data$recovery_right) & predicted_median > data$recovery_right ~
        predicted_median - data$recovery_right,
      is.finite(data$recovery_right) &
        predicted_median >= data$recovery_left &
        predicted_median <= data$recovery_right ~ 0,
      is.infinite(data$recovery_right) & predicted_median < data$recovery_left ~
        data$recovery_left - predicted_median,
      is.infinite(data$recovery_right) & predicted_median >= data$recovery_left ~ 0,
      TRUE ~ NA_real_
    )
    
    tibble(
      Model = x$Model,
      Family = x$Family,
      Structure = x$Structure,
      N = sum(is.finite(interval_error)),
      MAE = mean(interval_error, na.rm = TRUE),
      RMSE = sqrt(mean(interval_error^2, na.rm = TRUE)),
      Median_Error = median(interval_error, na.rm = TRUE),
      Successful = TRUE,
      Error = NA_character_
    )
  },
  error = function(e) {
    tibble(
      Model = x$Model,
      Family = x$Family,
      Structure = x$Structure,
      N = NA_integer_,
      MAE = NA_real_,
      RMSE = NA_real_,
      Median_Error = NA_real_,
      Successful = FALSE,
      Error = conditionMessage(e)
    )
  })
}

# ============================================================
# WITH MAGNITUDE + DEPTH
# ============================================================

set.seed(2026)

bayes_mae_rmse_all <- map_dfr(
  all_fit_results,
  ~ calculate_bayes_mae_rmse(.x, bayes_aft_data)
)

print(bayes_mae_rmse_all, n = Inf, width = Inf)

bayesian_loo_table_full <- bayesian_loo_table %>%
  left_join(
    bayes_mae_rmse_all %>%
      select(Model, Family, Structure, MAE, RMSE, Median_Error),
    by = c("Model", "Family", "Structure")
  ) %>%
  relocate(MAE, RMSE, Median_Error, .after = `SE ELPD`)

print(bayesian_loo_table_full, n = Inf, width = Inf)

# ============================================================
# BEST WITH MAGNITUDE + DEPTH
# ============================================================

best_bayes_elpd <- bayesian_loo_table_full %>%
  filter(is.finite(`ELPD-LOO`)) %>%
  slice_max(`ELPD-LOO`, n = 1, with_ties = FALSE)

best_bayes_mae <- bayesian_loo_table_full %>%
  filter(is.finite(MAE)) %>%
  slice_min(MAE, n = 1, with_ties = FALSE)

best_bayes_rmse <- bayesian_loo_table_full %>%
  filter(is.finite(RMSE)) %>%
  slice_min(RMSE, n = 1, with_ties = FALSE)

best_bayes_elpd
best_bayes_mae
best_bayes_rmse

bayesian_loo_mae_rmse_table <- bayesian_loo_table_full %>%
  mutate(
    MAE = round(MAE, 2),
    RMSE = round(RMSE, 2),
    Median_Error = round(Median_Error, 2)
  ) %>%
  arrange(Rank)

print(bayesian_loo_mae_rmse_table, n = Inf, width = Inf)

# ============================================================
# ADD TO ORIGINAL RESULTS
# ============================================================

bayesian_aft_results$mae_rmse_all <- bayes_mae_rmse_all
bayesian_aft_results$loo_mae_rmse_table <- bayesian_loo_mae_rmse_table
bayesian_aft_results$best_by_elpd <- best_bayes_elpd
bayesian_aft_results$best_by_mae <- best_bayes_mae
bayesian_aft_results$best_by_rmse <- best_bayes_rmse

saveRDS(
  bayesian_aft_results,
  file = "bayesian_interval_aft_loo_results.rds"
)

# ============================================================
# WITHOUT MAGNITUDE + DEPTH
# ============================================================

set.seed(2026)

bayes_mae_rmse_all_noeq <- map_dfr(
  all_fit_results_noeq,
  ~ calculate_bayes_mae_rmse(.x, bayes_aft_data_noeq)
)

print(bayes_mae_rmse_all_noeq, n = Inf, width = Inf)

bayesian_loo_table_full_noeq <- bayesian_loo_table_noeq %>%
  left_join(
    bayes_mae_rmse_all_noeq %>%
      select(Model, Family, Structure, MAE, RMSE, Median_Error),
    by = c("Model", "Family", "Structure")
  ) %>%
  relocate(MAE, RMSE, Median_Error, .after = `SE ELPD`)

print(bayesian_loo_table_full_noeq, n = Inf, width = Inf)

# ============================================================
# BEST WITHOUT MAGNITUDE + DEPTH
# ============================================================

best_bayes_elpd_noeq <- bayesian_loo_table_full_noeq %>%
  filter(is.finite(`ELPD-LOO`)) %>%
  slice_max(`ELPD-LOO`, n = 1, with_ties = FALSE)

best_bayes_mae_noeq <- bayesian_loo_table_full_noeq %>%
  filter(is.finite(MAE)) %>%
  slice_min(MAE, n = 1, with_ties = FALSE)

best_bayes_rmse_noeq <- bayesian_loo_table_full_noeq %>%
  filter(is.finite(RMSE)) %>%
  slice_min(RMSE, n = 1, with_ties = FALSE)

best_bayes_elpd_noeq
best_bayes_mae_noeq
best_bayes_rmse_noeq

bayesian_loo_mae_rmse_table_noeq <- bayesian_loo_table_full_noeq %>%
  mutate(
    MAE = round(MAE, 2),
    RMSE = round(RMSE, 2),
    Median_Error = round(Median_Error, 2)
  ) %>%
  arrange(Rank)

print(bayesian_loo_mae_rmse_table_noeq, n = Inf, width = Inf)

# ============================================================
# ADD TO NOEQ RESULTS
# ============================================================

bayesian_aft_results_noeq$mae_rmse_all <- bayes_mae_rmse_all_noeq
bayesian_aft_results_noeq$loo_mae_rmse_table <- bayesian_loo_mae_rmse_table_noeq
bayesian_aft_results_noeq$best_by_elpd <- best_bayes_elpd_noeq
bayesian_aft_results_noeq$best_by_mae <- best_bayes_mae_noeq
bayesian_aft_results_noeq$best_by_rmse <- best_bayes_rmse_noeq

saveRDS(
  bayesian_aft_results_noeq,
  file = "bayesian_interval_aft_no_magnitude_depth_results.rds"
)

# ============================================================
# COMPARE BEST ELPD MODELS
# ============================================================

best_bayesian_sensitivity_comparison <- bind_rows(
  best_bayes_elpd %>%
    mutate(Analysis = "With magnitude and depth"),
  best_bayes_elpd_noeq %>%
    mutate(Analysis = "Without magnitude and depth")
) %>%
  select(
    Analysis,
    Rank,
    Model,
    Family,
    Structure,
    `ELPD-LOO`,
    `SE ELPD`,
    MAE,
    RMSE,
    Median_Error,
    `Delta ELPD`,
    `SE Delta ELPD`,
    LOOIC,
    `LOO method`,
    `Max Pareto k`
  )

print(best_bayesian_sensitivity_comparison, n = Inf, width = Inf)

saveRDS(
  best_bayesian_sensitivity_comparison,
  file = "bayesian_aft_with_vs_without_earthquake_covariates.rds"
)




