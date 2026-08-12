# ============================================================
# PACKAGES
# ============================================================

library(readr)
library(brms)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
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

data_table <- data_table %>%
  mutate(
    event = ifelse(first_event == "recovered", 1, 0)
  )

bayes_data <- data_table %>%
  mutate(
    earthquake_id = factor(eq_id),
    log_gdp = log(gdp2015),
    censoring = if_else(event == 1, "none", "right"),
    censoring = factor(
      censoring,
      levels = c("none", "right")
    )
  )

table(
  event = bayes_data$event,
  censoring = bayes_data$censoring
)

nlevels(bayes_data$earthquake_id) # 17
table(bayes_data$earthquake_id)
colnames(bayes_data)


# ============================================================
# STANDARDISE CONTINUOUS VARIABLES
# ============================================================

continuous_vars <- c(
  "magnitude",
  "depth",
  "mean_pga",
  "max_pga",
  "gdp2015",
  "ntl_decline_ratio",
  "population",
  "log_gdp",
  "dis_expen",
  "maintain_expen",
  "restore_expen"
)

scaling_parameters <- map_dfr(
  continuous_vars,
  function(variable_name) {
    tibble(
      variable = variable_name,
      mean = mean(
        bayes_data[[variable_name]],
        na.rm = TRUE
      ),
      sd = sd(
        bayes_data[[variable_name]],
        na.rm = TRUE
      )
    )
  }
)

scaling_parameters

zero_variance_variables <- scaling_parameters %>%
  filter(!is.finite(sd) | sd == 0)

zero_variance_variables

bayes_data <- bayes_data %>%
  mutate(
    across(
      all_of(continuous_vars),
      ~ as.numeric(scale(.x))
    )
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
# BUILD BAYESIAN FORMULAS
# ============================================================

extract_formula_rhs <- function(formula_object) {
  if (!inherits(formula_object, "formula")) {
    stop("The supplied object is not a formula.")
  }
  
  paste(
    deparse(formula_object[[3]]),
    collapse = ""
  )
}

make_bayesian_formula <- function(
    formula_object,
    earthquake_effect = c("NoRE", "Random")) {
  
  earthquake_effect <- match.arg(earthquake_effect)
  rhs <- extract_formula_rhs(formula_object)
  
  if (earthquake_effect == "Random") {
    rhs <- paste0(
      rhs,
      " + (1 | earthquake_id)"
    )
  }
  
  formula_text <- paste0(
    "recovery_days | cens(censoring) ~ ",
    rhs
  )
  
  brms::bf(
    as.formula(formula_text)
  )
}

bayesian_formulas <- crossing(
  model = names(candidate_formulas),
  earthquake_effect = c("NoRE", "Random")
) %>%
  mutate(
    formula = map2(
      model,
      earthquake_effect,
      function(model_name, effect_type) {
        make_bayesian_formula(
          formula_object = candidate_formulas[[model_name]],
          earthquake_effect = effect_type
        )
      }
    )
  )

bayesian_formulas %>%
  dplyr::select(
    model,
    earthquake_effect,
    formula
  )

# bayesian_formulas$formula[[28]]


# ============================================================
# MODEL GRID
# ============================================================

bayesian_model_grid <- crossing(
  distribution = c("Lognormal", "Weibull"),
  model = names(candidate_formulas),
  earthquake_effect = c("NoRE", "Random")
) %>%
  left_join(
    bayesian_formulas,
    by = c(
      "model",
      "earthquake_effect"
    )
  ) %>%
  mutate(
    fit_name = paste(
      distribution,
      model,
      earthquake_effect,
      sep = "_"
    )
  ) %>%
  arrange(
    distribution,
    earthquake_effect,
    model
  )

bayesian_model_grid
nrow(bayesian_model_grid)

bayesian_model_grid %>%
  dplyr::select(
    fit_name,
    distribution,
    model,
    earthquake_effect
  ) %>%
  print(n = 56)


# ============================================================
# PRIOR SPECIFICATION
# ============================================================

get_prior(
  bf(
    recovery_days | cens(censoring) ~
      magnitude + depth + (1 | earthquake_id)
  ),
  data = bayes_data,
  family = lognormal()
)

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
  prior(
    student_t(3, 0, 0.35),
    class = "sigma"
  )
)

lognormal_random_priors <- c(
  lognormal_fixed_priors,
  random_effect_prior
)

# Weibull priors
weibull_fixed_priors <- c(
  common_priors,
  prior(
    lognormal(log(1.5), 0.25),
    class = "shape"
  )
)

weibull_random_priors <- c(
  weibull_fixed_priors,
  random_effect_prior
)


# ============================================================
# SELECT MODEL PRIORS
# ============================================================

select_model_priors <- function(
    distribution,
    earthquake_effect) {
  
  distribution <- match.arg(
    distribution,
    choices = c("Lognormal", "Weibull")
  )
  
  earthquake_effect <- match.arg(
    earthquake_effect,
    choices = c("NoRE", "Random")
  )
  
  if (
    distribution == "Lognormal" &&
    earthquake_effect == "NoRE"
  ) {
    return(lognormal_fixed_priors)
  }
  
  if (
    distribution == "Lognormal" &&
    earthquake_effect == "Random"
  ) {
    return(lognormal_random_priors)
  }
  
  if (
    distribution == "Weibull" &&
    earthquake_effect == "NoRE"
  ) {
    return(weibull_fixed_priors)
  }
  
  if (
    distribution == "Weibull" &&
    earthquake_effect == "Random"
  ) {
    return(weibull_random_priors)
  }
  
  stop(
    "Unsupported distribution or earthquake-effect combination."
  )
}


# ============================================================
# VALIDATE PRIORS
# ============================================================

example_random_formula <- make_bayesian_formula(
  candidate_formulas$M6,
  earthquake_effect = "Random"
)

get_prior(
  formula = example_random_formula,
  data = bayes_data,
  family = lognormal()
)

example_nore_formula <- make_bayesian_formula(
  candidate_formulas$M6,
  earthquake_effect = "NoRE"
)

example_random_formula <- make_bayesian_formula(
  candidate_formulas$M6,
  earthquake_effect = "Random"
)

validate_prior(
  prior = select_model_priors(
    distribution = "Lognormal",
    earthquake_effect = "NoRE"
  ),
  formula = example_nore_formula,
  data = bayes_data,
  family = lognormal()
)

validate_prior(
  prior = select_model_priors(
    distribution = "Lognormal",
    earthquake_effect = "Random"
  ),
  formula = example_random_formula,
  data = bayes_data,
  family = lognormal()
)

validate_prior(
  prior = select_model_priors(
    distribution = "Weibull",
    earthquake_effect = "NoRE"
  ),
  formula = example_nore_formula,
  data = bayes_data,
  family = weibull()
)

validate_prior(
  prior = select_model_priors(
    distribution = "Weibull",
    earthquake_effect = "Random"
  ),
  formula = example_random_formula,
  data = bayes_data,
  family = weibull()
)

# All four run without error


# ============================================================
# PRIOR PREDICTIVE CHECK MODELS
# ============================================================

prior_check_grid <- tribble(
  ~check_name, ~model, ~distribution, ~earthquake_effect,
  "Lognormal_M1_NoRE", "M1", "Lognormal", "NoRE",
  "Lognormal_M9_Random", "M9", "Lognormal", "Random",
  "Weibull_M1_NoRE", "M1", "Weibull", "NoRE",
  "Weibull_M9_Random", "M9", "Weibull", "Random"
)

select_family <- function(distribution) {
  switch(
    distribution,
    Lognormal = brms::lognormal(),
    Weibull = brms::weibull(),
    stop("Unknown distribution: ", distribution)
  )
}

fit_prior_check <- function(
    model,
    distribution,
    earthquake_effect,
    check_name) {
  
  model_formula <- make_bayesian_formula(
    formula_object = candidate_formulas[[model]],
    earthquake_effect = earthquake_effect
  )
  
  model_priors <- select_model_priors(
    distribution = distribution,
    earthquake_effect = earthquake_effect
  )
  
  brm(
    formula = model_formula,
    data = bayes_data,
    family = select_family(distribution),
    prior = model_priors,
    sample_prior = "only",
    chains = 4,
    iter = 4000,
    warmup = 1000,
    backend = "rstan",
    seed = 2026,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    ),
    file = file.path(
      "prior_predictive_models",
      check_name
    ),
    file_refit = "on_change",
    refresh = 500
  )
}

dir.create(
  "prior_predictive_models",
  showWarnings = FALSE,
  recursive = TRUE
)

prior_check_fits <- pmap(
  prior_check_grid,
  function(
    check_name,
    model,
    distribution,
    earthquake_effect) {
    
    message(
      "Prior check: ",
      check_name
    )
    
    fit_prior_check(
      model = model,
      distribution = distribution,
      earthquake_effect = earthquake_effect,
      check_name = check_name
    )
  }
)

names(prior_check_fits) <- prior_check_grid$check_name


# ============================================================
# PRIOR PREDICTIVE SUMMARY
# ============================================================

summarise_prior_fit <- function(
    fit,
    model_name,
    ndraws = 1000) {
  
  yrep <- posterior_predict(
    fit,
    ndraws = ndraws
  )
  
  values <- as.vector(yrep)
  
  tibble(
    model = model_name,
    minimum = min(values),
    q01 = quantile(values, 0.01),
    q05 = quantile(values, 0.05),
    median = median(values),
    q95 = quantile(values, 0.95),
    q99 = quantile(values, 0.99),
    maximum = max(values),
    proportion_below_1_day = mean(values < 1),
    proportion_above_5_years = mean(values > 5 * 365),
    proportion_above_10_years = mean(values > 10 * 365)
  )
}

prior_predictive_summary <- imap_dfr(
  prior_check_fits,
  summarise_prior_fit
)

prior_predictive_summary


# ============================================================
# PRIOR PREDICTIVE PLOTS
# ============================================================

prior_predictive_plots <- imap(
  prior_check_fits,
  function(fit, model_name) {
    
    yrep <- posterior_predict(
      fit,
      ndraws = 50
    )
    
    bayesplot::ppc_dens_overlay(
      y = log1p(bayes_data$recovery_days),
      yrep = log1p(yrep)
    ) +
      labs(
        title = model_name,
        x = "log(1 + recovery days)"
      )
  }
)

prior_predictive_plots$Lognormal_M1_NoRE
prior_predictive_plots$Lognormal_M9_Random
prior_predictive_plots$Weibull_M1_NoRE
prior_predictive_plots$Weibull_M9_Random

prior_predictive_plots <- imap(
  prior_check_fits,
  function(fit, model_name) {
    
    yrep <- posterior_predict(
      fit,
      ndraws = 50
    )
    
    bayesplot::ppc_dens_overlay(
      y = bayes_data$recovery_days,
      yrep = yrep
    ) +
      coord_cartesian(
        xlim = c(0, 1000)
      ) +
      labs(
        title = model_name,
        x = "Recovery days"
      )
  }
)

prior_predictive_plots$Lognormal_M1_NoRE
prior_predictive_plots$Lognormal_M9_Random
prior_predictive_plots$Weibull_M1_NoRE
prior_predictive_plots$Weibull_M9_Random


# ============================================================
# CHECK PRIOR PREDICTIVE DRAWS
# ============================================================

yrep <- posterior_predict(
  prior_check_fits$Lognormal_M9_Random,
  ndraws = 5
)

dim(yrep)
summary(as.vector(yrep))
apply(yrep, 1, summary)


# ============================================================
# TEST MODEL
# ============================================================

dir.create(
  "bayesian_test_models",
  showWarnings = FALSE,
  recursive = TRUE
)

test_formula <- make_bayesian_formula(
  candidate_formulas$M6,
  earthquake_effect = "NoRE"
)

test_fit <- brm(
  formula = test_formula,
  data = bayes_data,
  family = lognormal(),
  prior = select_model_priors(
    distribution = "Lognormal",
    earthquake_effect = "NoRE"
  ),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  backend = "rstan",
  seed = 20,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  save_pars = save_pars(all = TRUE),
  file = "bayesian_test_models/Lognormal_M6_NoRE",
  file_refit = "on_change"
)

summary(test_fit)

posterior::summarise_draws(
  posterior::as_draws_df(test_fit),
  posterior::default_convergence_measures()
)

nuts_params(test_fit) %>%
  filter(Parameter == "divergent__") %>%
  summarise(
    divergences = sum(Value)
  )

pp_check(
  test_fit,
  type = "dens_overlay",
  ndraws = 100
)


# ============================================================
# FIT ALL BAYESIAN AFT MODELS
# ============================================================

fit_one_model <- function(
    formula,
    distribution,
    earthquake_effect,
    fit_name) {
  
  dir.create(
    "bayesian_aft_models",
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  brm(
    formula = formula,
    data = bayes_data,
    family = select_family(distribution),
    prior = select_model_priors(
      distribution,
      earthquake_effect
    ),
    chains = 4,
    iter = 8000,
    warmup = 4000,
    backend = "rstan",
    seed = 2026,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    save_pars = save_pars(all = TRUE),
    file = file.path(
      "bayesian_aft_models",
      fit_name
    ),
    file_refit = "on_change",
    refresh = 500
  )
}

bayesian_fits <- pmap(
  bayesian_model_grid,
  function(
    distribution,
    model,
    earthquake_effect,
    formula,
    fit_name) {
    
    message(
      "Fitting: ",
      fit_name
    )
    
    tryCatch(
      fit_one_model(
        formula = formula,
        distribution = distribution,
        earthquake_effect = earthquake_effect,
        fit_name = fit_name
      ),
      error = function(e) {
        structure(
          list(
            fit_name = fit_name,
            error_message = conditionMessage(e)
          ),
          class = "bayesian_fit_error"
        )
      }
    )
  }
)

length(
  list.files("bayesian_aft_models")
)

names(bayesian_fits) <- bayesian_model_grid$fit_name

successful_fits <- bayesian_fits[
  !map_lgl(
    bayesian_fits,
    ~ inherits(.x, "bayesian_fit_error")
  )
]


# ============================================================
# INITIAL PSIS-LOO
# ============================================================

loo_objects <- imap(
  successful_fits,
  function(fit, fit_name) {
    
    message(
      "LOO: ",
      fit_name
    )
    
    loo(
      fit,
      moment_match = FALSE
    )
  }
)


# ============================================================
# PARETO-K DIAGNOSTICS
# ============================================================

pareto_summary <- imap_dfr(
  loo_objects,
  function(loo_object, fit_name) {
    
    k <- loo::pareto_k_values(
      loo_object
    )
    
    tibble(
      fit_name = fit_name,
      max_pareto_k = max(k),
      n_k_above_0_7 = sum(k > 0.7),
      n_k_above_1 = sum(k > 1)
    )
  }
)

pareto_summary

problem_models <- pareto_summary %>%
  filter(
    n_k_above_0_7 > 0
  ) %>%
  pull(fit_name)

problem_models


# ============================================================
# MOMENT MATCHING
# ============================================================

for (i in seq_along(problem_models)) {
  fit_name <- problem_models[i]
  
  message(
    sprintf(
      "[%d/%d] Moment matching: %s",
      i,
      length(problem_models),
      fit_name
    )
  )
  
  loo_objects[[fit_name]] <- loo(
    successful_fits[[fit_name]],
    moment_match = TRUE
  )
}

saveRDS(
  loo_objects,
  file = "loo_objects_moment_matched.rds"
)


# ============================================================
# CHECK AFTER MOMENT MATCHING
# ============================================================

pareto_after_mm <- imap_dfr(
  loo_objects,
  function(loo_object, fit_name) {
    
    k_values <- loo::pareto_k_values(
      loo_object
    )
    
    tibble(
      fit_name = fit_name,
      max_pareto_k = max(
        k_values,
        na.rm = TRUE
      ),
      n_k_above_0_7 = sum(
        k_values > 0.7,
        na.rm = TRUE
      ),
      n_k_above_1 = sum(
        k_values > 1,
        na.rm = TRUE
      )
    )
  }
) %>%
  arrange(
    desc(max_pareto_k)
  )

print(
  pareto_after_mm,
  n = 56
)

remaining_problem_models <- pareto_after_mm %>%
  filter(
    n_k_above_0_7 > 0
  ) %>%
  pull(fit_name)

length(remaining_problem_models)
remaining_problem_models


# ============================================================
# EXACT RELOO
# ============================================================

for (i in seq_along(remaining_problem_models)) {
  fit_name <- remaining_problem_models[i]
  
  message(
    sprintf(
      "[%d/%d] Exact reloo: %s",
      i,
      length(remaining_problem_models),
      fit_name
    )
  )
  
  loo_objects[[fit_name]] <- tryCatch(
    brms::reloo(
      x = successful_fits[[fit_name]],
      loo = loo_objects[[fit_name]],
      k_threshold = 0.7
    ),
    error = function(e) {
      warning(
        "reloo failed for ",
        fit_name,
        ": ",
        conditionMessage(e)
      )
      
      loo_objects[[fit_name]]
    }
  )
}

saveRDS(
  loo_objects,
  file = "loo_objects_final.rds"
)


# ============================================================
# FINAL PARETO-K CHECK
# ============================================================

loo_objects <- readRDS(
  "loo_objects_final.rds"
)

pareto_final <- imap_dfr(
  loo_objects,
  function(x, name) {
    
    k <- loo::pareto_k_values(x)
    
    tibble(
      fit_name = name,
      max_k = max(k),
      n_k_above_0_7 = sum(k > 0.7)
    )
  }
)

pareto_final %>%
  filter(
    n_k_above_0_7 > 0
  )


# ============================================================
# ELPD MODEL COMPARISON
# ============================================================

elpd_table <- imap_dfr(
  loo_objects,
  function(loo_object, fit_name) {
    
    estimates <- loo_object$estimates
    
    tibble(
      fit_name = fit_name,
      elpd_loo = estimates["elpd_loo", "Estimate"],
      se_elpd = estimates["elpd_loo", "SE"],
      p_loo = estimates["p_loo", "Estimate"],
      looic = estimates["looic", "Estimate"],
      se_looic = estimates["looic", "SE"]
    )
  }
) %>%
  arrange(
    desc(elpd_loo)
  )

print(
  elpd_table,
  n = 56
)

loo_comparison <- loo::loo_compare(
  loo_objects
)

loo_comparison_table <- as.data.frame(
  loo_comparison
) %>%
  tibble::rownames_to_column(
    "fit_name"
  ) %>%
  as_tibble()

print(
  loo_comparison_table,
  n = 56,
  width = Inf
)


# ============================================================
# BEST FULL MODEL
# ============================================================

best_fit <- successful_fits[
  ["Lognormal_M13_Random"]
]

epred <- posterior_epred(
  best_fit
)

predicted_days <- apply(
  epred,
  2,
  median
)

prediction_df <- bayes_data %>%
  mutate(
    predicted_days = predicted_days
  ) %>%
  filter(
    event == 1
  )

ggplot(
  prediction_df,
  aes(
    x = recovery_days,
    y = predicted_days
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    x = "Observed recovery time (days)",
    y = "Predicted recovery time (days)"
  ) +
  theme_bw()

pp_check(
  best_fit,
  type = "dens_overlay",
  ndraws = 50
)

pp_check(
  best_fit,
  type = "ecdf_overlay",
  ndraws = 50
)


# ============================================================
# CONDITIONAL EFFECTS FOR FULL MODEL
# ============================================================

newdata <- bayes_data %>%
  summarise(
    magnitude = mean(magnitude),
    depth = mean(depth),
    mean_pga = mean(mean_pga),
    ntl_decline_ratio = mean(ntl_decline_ratio),
    log_gdp = mean(log_gdp),
    population = mean(population),
    dis_expen = mean(dis_expen),
    earthquake_id = earthquake_id[1]
  )

conditional_effects(
  best_fit,
  method = "posterior_epred"
)


# ============================================================
# REDUCED CANDIDATE MODELS
# ============================================================

candidate_formulas_reduced <- map(
  candidate_formulas,
  ~ update(
    .x,
    . ~ . - magnitude - depth
  )
)

candidate_formulas_reduced

imap(
  candidate_formulas_reduced,
  ~ cat(
    .y,
    ": ",
    deparse(.x),
    "\n"
  )
)

reduced_bayesian_formulas <- tibble(
  model = names(candidate_formulas_reduced)
) %>%
  mutate(
    formula = map(
      model,
      ~ make_bayesian_formula(
        formula_object = candidate_formulas_reduced[[.x]],
        earthquake_effect = "Random"
      )
    )
  )

reduced_bayesian_formulas %>%
  select(
    model,
    formula
  ) %>%
  print(n = 14)


# ============================================================
# REDUCED MODEL GRID
# ============================================================

reduced_bayesian_model_grid <- crossing(
  distribution = c("Lognormal", "Weibull"),
  model = names(candidate_formulas_reduced)
) %>%
  left_join(
    reduced_bayesian_formulas,
    by = "model"
  ) %>%
  mutate(
    earthquake_effect = "Random",
    fit_name = paste(
      distribution,
      model,
      "Random_NoMagDepth",
      sep = "_"
    )
  ) %>%
  arrange(
    distribution,
    model
  )

print(
  reduced_bayesian_model_grid %>%
    select(
      fit_name,
      distribution,
      model,
      earthquake_effect
    ),
  n = 28
)


# ============================================================
# FIT REDUCED BAYESIAN MODELS
# ============================================================

dir.create(
  "bayesian_reduced_28_models",
  recursive = TRUE,
  showWarnings = FALSE
)

fit_one_reduced_model <- function(
    formula,
    distribution,
    fit_name) {
  
  model_path <- file.path(
    "bayesian_reduced_28_models",
    fit_name
  )
  
  model_file <- paste0(
    model_path,
    ".rds"
  )
  
  if (file.exists(model_file)) {
    message(
      "Loading existing model: ",
      fit_name
    )
    
    return(
      readRDS(model_file)
    )
  }
  
  message(
    "Fitting reduced model: ",
    fit_name
  )
  
  brm(
    formula = formula,
    data = bayes_data,
    family = select_family(distribution),
    prior = select_model_priors(
      distribution = distribution,
      earthquake_effect = "Random"
    ),
    chains = 4,
    iter = 8000,
    warmup = 4000,
    backend = "rstan",
    seed = 2026,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    save_pars = save_pars(all = TRUE),
    file = model_path,
    file_refit = "on_change",
    refresh = 500
  )
}

reduced_bayesian_fits <- pmap(
  reduced_bayesian_model_grid,
  function(
    distribution,
    model,
    formula,
    earthquake_effect,
    fit_name) {
    
    tryCatch(
      fit_one_reduced_model(
        formula = formula,
        distribution = distribution,
        fit_name = fit_name
      ),
      error = function(e) {
        structure(
          list(
            fit_name = fit_name,
            error_message = conditionMessage(e)
          ),
          class = "bayesian_fit_error"
        )
      }
    )
  }
)

names(reduced_bayesian_fits) <-
  reduced_bayesian_model_grid$fit_name


# ============================================================
# REDUCED MODEL FIT STATUS
# ============================================================

reduced_fit_status <- imap_dfr(
  reduced_bayesian_fits,
  function(fit, fit_name) {
    
    tibble(
      fit_name = fit_name,
      successful = !inherits(
        fit,
        "bayesian_fit_error"
      ),
      error_message = if (
        inherits(
          fit,
          "bayesian_fit_error"
        )
      ) {
        fit$error_message
      } else {
        NA_character_
      }
    )
  }
)

reduced_fit_status %>%
  count(successful)

reduced_fit_status %>%
  filter(
    !successful
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

saveRDS(
  reduced_bayesian_fits,
  "reduced_bayesian_fits_28.rds"
)

saveRDS(
  reduced_bayesian_model_grid,
  "reduced_bayesian_model_grid_28.rds"
)


# ============================================================
# KEEP SUCCESSFUL REDUCED FITS
# ============================================================

successful_reduced_fits <- reduced_bayesian_fits[
  !map_lgl(
    reduced_bayesian_fits,
    ~ inherits(.x, "bayesian_fit_error")
  )
]

successful_index <- !map_lgl(
  reduced_bayesian_fits,
  ~ inherits(.x, "bayesian_fit_error")
)

names(successful_reduced_fits) <-
  reduced_bayesian_model_grid$fit_name[
    successful_index
  ]

names(successful_reduced_fits)
length(successful_reduced_fits)


# ============================================================
# SAMPLING DIAGNOSTICS
# ============================================================

reduced_sampling_diagnostics <- imap_dfr(
  successful_reduced_fits,
  function(fit, fit_name) {
    
    fit_summary <- posterior::summarise_draws(
      posterior::as_draws_df(fit),
      posterior::default_convergence_measures()
    )
    
    nuts <- brms::nuts_params(fit)
    
    tibble(
      fit_name = fit_name,
      max_rhat = max(
        fit_summary$rhat,
        na.rm = TRUE
      ),
      min_bulk_ess = min(
        fit_summary$ess_bulk,
        na.rm = TRUE
      ),
      min_tail_ess = min(
        fit_summary$ess_tail,
        na.rm = TRUE
      ),
      divergences = sum(
        nuts$Value[
          nuts$Parameter == "divergent__"
        ]
      ),
      max_treedepth_hits = sum(
        nuts$Value[
          nuts$Parameter == "treedepth__"
        ] >= 15
      )
    )
  }
)

print(
  reduced_sampling_diagnostics,
  n = 26,
  width = Inf
)

reduced_sampling_diagnostics %>%
  filter(
    max_rhat >= 1.01 |
      divergences > 0 |
      max_treedepth_hits > 0
  )


# ============================================================
# INITIAL PSIS-LOO FOR REDUCED MODELS
# ============================================================

reduced_loo_objects <- imap(
  successful_reduced_fits,
  function(fit, fit_name) {
    
    message(
      "LOO: ",
      fit_name
    )
    
    brms::loo(
      fit,
      moment_match = FALSE
    )
  }
)

saveRDS(
  reduced_loo_objects,
  "reduced_loo_objects_initial.rds"
)


# ============================================================
# REDUCED PARETO-K DIAGNOSTICS
# ============================================================

reduced_pareto_summary <- imap_dfr(
  reduced_loo_objects,
  function(loo_object, fit_name) {
    
    k <- loo::pareto_k_values(
      loo_object
    )
    
    tibble(
      fit_name = fit_name,
      max_pareto_k = max(
        k,
        na.rm = TRUE
      ),
      n_k_above_0_7 = sum(
        k > 0.7,
        na.rm = TRUE
      ),
      n_k_above_1 = sum(
        k > 1,
        na.rm = TRUE
      )
    )
  }
) %>%
  arrange(
    desc(max_pareto_k)
  )

print(
  reduced_pareto_summary,
  n = 26,
  width = Inf
)

reduced_problem_models <- reduced_pareto_summary %>%
  filter(
    n_k_above_0_7 > 0
  ) %>%
  pull(fit_name)

length(reduced_problem_models) # 24
reduced_problem_models


# ============================================================
# MOMENT MATCHING FOR REDUCED MODELS
# ============================================================

for (i in seq_along(reduced_problem_models)) {
  fit_name <- reduced_problem_models[i]
  
  message(
    sprintf(
      "[%d/%d] Moment matching: %s",
      i,
      length(reduced_problem_models),
      fit_name
    )
  )
  
  reduced_loo_objects[[fit_name]] <- tryCatch(
    brms::loo(
      successful_reduced_fits[[fit_name]],
      moment_match = TRUE
    ),
    error = function(e) {
      warning(
        "Moment matching failed for ",
        fit_name,
        ": ",
        conditionMessage(e)
      )
      
      reduced_loo_objects[[fit_name]]
    }
  )
}

saveRDS(
  reduced_loo_objects,
  "reduced_loo_objects_moment_matched.rds"
)


# ============================================================
# CHECK AFTER MOMENT MATCHING
# ============================================================

reduced_pareto_after_mm <- imap_dfr(
  reduced_loo_objects,
  function(loo_object, fit_name) {
    
    k <- loo::pareto_k_values(
      loo_object
    )
    
    tibble(
      fit_name = fit_name,
      max_pareto_k = max(
        k,
        na.rm = TRUE
      ),
      n_k_above_0_7 = sum(
        k > 0.7,
        na.rm = TRUE
      ),
      n_k_above_1 = sum(
        k > 1,
        na.rm = TRUE
      )
    )
  }
) %>%
  arrange(
    desc(max_pareto_k)
  )

reduced_remaining_problem_models <-
  reduced_pareto_after_mm %>%
  filter(
    n_k_above_0_7 > 0
  ) %>%
  pull(fit_name)

reduced_remaining_problem_models


# ============================================================
# EXACT RELOO FOR REDUCED MODELS
# ============================================================

for (i in seq_along(reduced_remaining_problem_models)) {
  fit_name <- reduced_remaining_problem_models[i]
  
  message(
    sprintf(
      "[%d/%d] Exact reloo: %s",
      i,
      length(reduced_remaining_problem_models),
      fit_name
    )
  )
  
  reduced_loo_objects[[fit_name]] <- brms::reloo(
    x = successful_reduced_fits[[fit_name]],
    loo = reduced_loo_objects[[fit_name]],
    k_threshold = 0.7
  )
}

saveRDS(
  reduced_loo_objects,
  "reduced_loo_objects_final.rds"
)


# ============================================================
# REDUCED MODEL ELPD TABLE
# ============================================================

reduced_elpd_table <- imap_dfr(
  reduced_loo_objects,
  function(loo_object, fit_name) {
    
    estimates <- loo_object$estimates
    
    tibble(
      fit_name = fit_name,
      elpd_loo = estimates["elpd_loo", "Estimate"],
      se_elpd = estimates["elpd_loo", "SE"],
      p_loo = estimates["p_loo", "Estimate"],
      looic = estimates["looic", "Estimate"],
      se_looic = estimates["looic", "SE"]
    )
  }
) %>%
  left_join(
    reduced_bayesian_model_grid %>%
      select(
        fit_name,
        distribution,
        model
      ),
    by = "fit_name"
  ) %>%
  arrange(
    desc(elpd_loo)
  )

print(
  reduced_elpd_table,
  n = 26,
  width = Inf
)


# ============================================================
# REDUCED MODEL COMPARISON
# ============================================================

reduced_loo_comparison <- loo::loo_compare(
  reduced_loo_objects
)

reduced_loo_comparison_table <- as.data.frame(
  reduced_loo_comparison
) %>%
  tibble::rownames_to_column(
    "fit_name"
  ) %>%
  as_tibble()

print(
  reduced_loo_comparison_table,
  n = 28,
  width = Inf
)


# ============================================================
# BAYESIAN PREDICTION METRICS
# ============================================================

calculate_bayesian_metrics <- function(fit, data) {
  
  epred <- posterior_epred(
    fit,
    newdata = data,
    re_formula = NULL
  )
  
  if (ncol(epred) != nrow(data)) {
    stop(
      "Predictions do not match the number of rows in data."
    )
  }
  
  predicted_days <- apply(
    epred,
    2,
    median,
    na.rm = TRUE
  )
  
  metric_data <- tibble(
    observation = seq_len(nrow(data)),
    observed_days = data$recovery_days,
    predicted_days = predicted_days,
    event = data$event
  ) %>%
    filter(
      event == 1,
      is.finite(observed_days),
      is.finite(predicted_days)
    ) %>%
    mutate(
      error = predicted_days - observed_days,
      absolute_error = abs(error),
      squared_error = error^2
    )
  
  metrics <- metric_data %>%
    summarise(
      n_events = n(),
      mae_event = mean(absolute_error),
      mse_event = mean(squared_error),
      rmse_event = sqrt(
        mean(squared_error)
      ),
      mean_error_event = mean(error),
      median_absolute_error_event = median(
        absolute_error
      )
    )
  
  list(
    metrics = metrics,
    observation_results = metric_data
  )
}


# ============================================================
# REDUCED MODEL PREDICTION METRICS
# ============================================================

reduced_prediction_metrics <- imap_dfr(
  successful_reduced_fits,
  function(fit, fit_name) {
    
    result <- calculate_bayesian_metrics(
      fit = fit,
      data = bayes_data
    )
    
    result$metrics %>%
      mutate(
        fit_name = fit_name,
        .before = 1
      )
  }
)

reduced_model_results <- reduced_elpd_table %>%
  left_join(
    reduced_prediction_metrics,
    by = "fit_name"
  ) %>%
  arrange(
    desc(elpd_loo)
  )

print(
  reduced_model_results,
  n = 26,
  width = Inf
)


# ============================================================
# FULL MODEL PREDICTION METRICS
# ============================================================

list.files(
  "bayesian_aft_models"
)

model_files <- list.files(
  "bayesian_aft_models",
  pattern = "\\.rds$",
  full.names = TRUE
)

bayesian_prediction_metrics <- map_dfr(
  model_files,
  function(model_file) {
    
    fit_name <- tools::file_path_sans_ext(
      basename(model_file)
    )
    
    message(
      "Calculating metrics: ",
      fit_name
    )
    
    fit <- readRDS(
      model_file
    )
    
    tryCatch(
      {
        result <- calculate_bayesian_metrics(
          fit = fit,
          data = bayes_data
        )
        
        result$metrics %>%
          mutate(
            fit_name = fit_name,
            .before = 1
          )
      },
      error = function(e) {
        tibble(
          fit_name = fit_name,
          n_events = NA_integer_,
          mae_event = NA_real_,
          mse_event = NA_real_,
          rmse_event = NA_real_,
          mean_error_event = NA_real_,
          median_absolute_error_event = NA_real_,
          error_message = conditionMessage(e)
        )
      }
    )
  }
)

nrow(bayesian_prediction_metrics)

print(
  bayesian_prediction_metrics,
  n = 56,
  width = Inf
)


# ============================================================
# BEST REDUCED MODEL
# ============================================================

best_reduced_fit <- successful_reduced_fits[
  ["Lognormal_M9_Random_NoMagDepth"]
]

# Posterior expected recovery-time draws
epred_reduced <- posterior_epred(
  best_reduced_fit,
  newdata = bayes_data,
  re_formula = NULL
)

# Posterior median prediction
predicted_days_reduced <- apply(
  epred_reduced,
  2,
  median,
  na.rm = TRUE
)

prediction_df_reduced <- bayes_data %>%
  mutate(
    predicted_days = predicted_days_reduced
  ) %>%
  filter(
    event == 1,
    is.finite(recovery_days),
    is.finite(predicted_days)
  )


# ============================================================
# REDUCED MODEL PREDICTION PLOT
# ============================================================

observed_predicted_reduced_plot <- ggplot(
  prediction_df_reduced,
  aes(
    x = recovery_days,
    y = predicted_days
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    title = "Observed versus Predicted Recovery Time",
    subtitle = paste(
      "Reduced Lognormal M9 model with",
      "earthquake-level random effects"
    ),
    x = "Observed recovery time (days)",
    y = "Posterior median predicted recovery time (days)"
  ) +
  theme_bw()

observed_predicted_reduced_plot


# ============================================================
# REDUCED MODEL POSTERIOR PREDICTIVE CHECK
# ============================================================

ecdf_check_reduced <- pp_check(
  best_reduced_fit,
  type = "ecdf_overlay",
  ndraws = 50
) +
  coord_cartesian(
    xlim = c(0, 1000)
  ) +
  labs(
    title = "Posterior Predictive ECDF Check",
    subtitle = paste(
      "Reduced Lognormal M9 model excluding",
      "earthquake magnitude and focal depth"
    ),
    x = "Recovery time (days)",
    y = "Empirical cumulative probability"
  ) +
  theme_bw()

ecdf_check_reduced

reduced_side_by_side <-
  observed_predicted_reduced_plot |
  ecdf_check_reduced

reduced_side_by_side


# ============================================================
# REDUCED MODEL CONDITIONAL EFFECTS
# ============================================================

representative_newdata_reduced <- bayes_data %>%
  summarise(
    mean_pga = mean(
      mean_pga,
      na.rm = TRUE
    ),
    ntl_decline_ratio = mean(
      ntl_decline_ratio,
      na.rm = TRUE
    ),
    gdp2015 = mean(
      gdp2015,
      na.rm = TRUE
    ),
    population = mean(
      population,
      na.rm = TRUE
    ),
    dis_expen = mean(
      dis_expen,
      na.rm = TRUE
    ),
    earthquake_id = first(
      earthquake_id
    )
  )

representative_newdata_reduced

reduced_conditional_effects <- conditional_effects(
  best_reduced_fit,
  method = "posterior_epred"
)

reduced_effect_plots <- plot(
  reduced_conditional_effects,
  plot = FALSE,
  ask = FALSE
)

combined_reduced_effects_plot <-
  wrap_plots(
    reduced_effect_plots,
    ncol = 2
  ) +
  plot_annotation(
    title = "Conditional Effects of the Best Reduced Bayesian AFT Model",
    subtitle = paste(
      "Lognormal M9 with earthquake-level random effects;",
      "magnitude and focal depth excluded"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        size = 11
      )
    )
  )

combined_reduced_effects_plot


# ============================================================
# BEST FULL MODEL PREDICTIONS
# ============================================================

best_fit <- readRDS(
  "bayesian_aft_models/Lognormal_M13_Random.rds"
)

# Posterior expected recovery-time draws
epred <- posterior_epred(
  best_fit,
  newdata = bayes_data,
  re_formula = NULL
)

# Posterior median prediction
predicted_days <- apply(
  epred,
  2,
  median,
  na.rm = TRUE
)

prediction_df <- bayes_data %>%
  mutate(
    predicted_days = predicted_days
  ) %>%
  filter(
    event == 1,
    is.finite(recovery_days),
    is.finite(predicted_days)
  )


# ============================================================
# FULL MODEL PREDICTION PLOT
# ============================================================

observed_predicted_plot <- ggplot(
  prediction_df,
  aes(
    x = recovery_days,
    y = predicted_days
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    title = "Observed versus Predicted Recovery Time",
    subtitle = "Lognormal M13 model with earthquake-level random effects",
    x = "Observed recovery time (days)",
    y = "Posterior median predicted recovery time (days)"
  ) +
  theme_bw()


# ============================================================
# FULL MODEL POSTERIOR PREDICTIVE CHECK
# ============================================================

ecdf_check <- pp_check(
  best_fit,
  type = "ecdf_overlay",
  ndraws = 50
) +
  coord_cartesian(
    xlim = c(0, 1000)
  ) +
  labs(
    title = "Posterior Predictive ECDF Check",
    subtitle = "Lognormal M13 model with earthquake-level random effects",
    x = "Recovery time (days)",
    y = "Empirical cumulative probability"
  ) +
  theme_bw()

full_side_by_side <-
  observed_predicted_plot |
  ecdf_check

full_side_by_side