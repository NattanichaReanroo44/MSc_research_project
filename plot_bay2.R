library(brms)
library(patchwork)
# ============================================================
# GET BEST BAYESIAN MODEL FIT
# ============================================================

best_model_name <- best_bayes_elpd_noeq$Model
best_family <- best_bayes_elpd_noeq$Family
best_structure <- best_bayes_elpd_noeq$Structure


if (best_family == "lognormal" && best_structure == "Fixed") {
  best_fit <- lognormal_fixed_results_noeq[[best_model_name]]$Fit
} else if (best_family == "lognormal" && best_structure == "Random earthquake") {
  best_fit <- lognormal_random_results_noeq[[best_model_name]]$Fit
} else if (best_family == "weibull" && best_structure == "Fixed") {
  best_fit <- weibull_fixed_results_noeq[[best_model_name]]$Fit
} else if (best_family == "weibull" && best_structure == "Random earthquake") {
  best_fit <- weibull_random_results_noeq[[best_model_name]]$Fit
}

summary(best_fit)

# ============================================================
# POSTERIOR COEFFICIENTS + TIME RATIOS
# ============================================================

best_coefficients <- posterior_summary(
  best_fit,
  pars = "^b_"
) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>%
  as_tibble() %>%
  mutate(
    Variable = sub("^b_", "", Parameter),
    Time_Ratio = exp(Estimate),
    TR_Lower = exp(Q2.5),
    TR_Upper = exp(Q97.5)
  )

print(best_coefficients, n = Inf, width = Inf)

# ============================================================
# TIME-RATIO FOREST PLOT
# ============================================================

best_time_ratio_plot <- best_coefficients %>%
  filter(Variable != "Intercept") %>%
  ggplot(
    aes(
      x = Time_Ratio,
      y = reorder(Variable, Time_Ratio)
    )
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(
      xmin = TR_Lower,
      xmax = TR_Upper
    ),
    height = 0.15
  ) +
  geom_point(size = 2.5) +
  labs(
    title = "Posterior time ratios",
    subtitle = paste(
      best_model_name,
      best_family,
      best_structure
    ),
    x = "Time ratio (95% credible interval)",
    y = NULL
  ) +
  theme_bw()

best_time_ratio_plot

ggsave(
  "time_ratio.png",
  best_time_ratio_plot,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)


# ============================================================
# TRACE PLOTS
# ============================================================


library(patchwork)

p1 <- plot(
  best_fit,
  variable = c(
    "b_Intercept",
    "b_ntl_decline_ratio",
    "b_log_gdp",
    "b_population"
  ),
  ask = FALSE
)

p2 <- plot(
  best_fit,
  variable = c(
    "b_dis_expen",
    "sd_earthquake_id__Intercept",
    "shape"
  ),
  ask = FALSE
)


top_plot <- wrap_plots(p1, ncol = 1)
bottom_plot <- wrap_plots(p2, ncol = 1)

combined_plot <- top_plot / bottom_plot

print(combined_plot)

ggsave(
  filename = "M12_MCMC_diagnostics.png",
  plot = combined_plot,
  width = 13,
  height = 15,
  dpi = 300
)

# ============================================================
# POSTERIOR PARAMETER DISTRIBUTIONS
# ============================================================

mcmc_plot(
  best_fit,
  type = "areas",
  pars = "^b_",
  prob = 0.95
)

# ============================================================
# CONDITIONAL EFFECTS:
# WITH VS WITHOUT MAGNITUDE/DEPTH
# ============================================================

# Best models
best_full_fit <- bayesian_aft_results$weibull_random[["M8"]]$Fit
best_noeq_fit <- bayesian_aft_results_noeq$weibull_random[["M12"]]$Fit

# ============================================================
# CONDITIONAL EFFECTS
# Population-level effects only
# ============================================================

full_effects <- conditional_effects(
  best_full_fit,
  re_formula = NA,
  robust = TRUE
)

noeq_effects <- conditional_effects(
  best_noeq_fit,
  re_formula = NA,
  robust = TRUE
)

# ============================================================
# HELPER FUNCTION
# ============================================================

extract_conditional_effect <- function(effect_list, variable, model_name) {
  effect_data <- effect_list[[variable]]
  predictor_column <- names(effect_data)[1]
  
  tibble(
    predictor_value = effect_data[[predictor_column]],
    estimate = effect_data$estimate__,
    lower = effect_data$lower__,
    upper = effect_data$upper__,
    variable = variable,
    model = model_name
  )
}

# ============================================================
# VARIABLES
# ============================================================

# Appearing in both M8 and M12
common_variables <- c(
  "ntl_decline_ratio",
  "population",
  "dis_expen"
)

# Appearing only in M8 with magnitude/depth
full_only_variables <- c(
  "magnitude",
  "depth",
  "gdp2015"
)

# Appearing only in M12 without magnitude/depth
noeq_only_variables <- c(
  "log_gdp"
)

# ============================================================
# COMMON EFFECTS
# ============================================================

common_effect_data <- map_dfr(
  common_variables,
  function(variable) {
    bind_rows(
      extract_conditional_effect(
        effect_list = full_effects,
        variable = variable,
        model_name = "With magnitude/depth (M8)"
      ),
      extract_conditional_effect(
        effect_list = noeq_effects,
        variable = variable,
        model_name = "Without magnitude/depth (M12)"
      )
    )
  }
)

# ============================================================
# FULL MODEL ONLY
# ============================================================

full_only_effect_data <- map_dfr(
  full_only_variables,
  function(variable) {
    extract_conditional_effect(
      effect_list = full_effects,
      variable = variable,
      model_name = "With magnitude/depth (M8)"
    )
  }
)

# ============================================================
# NOEQ MODEL ONLY
# ============================================================

noeq_only_effect_data <- map_dfr(
  noeq_only_variables,
  function(variable) {
    extract_conditional_effect(
      effect_list = noeq_effects,
      variable = variable,
      model_name = "Without magnitude/depth (M12)"
    )
  }
)

# ============================================================
# COMBINE
# ============================================================

all_effect_data <- bind_rows(
  common_effect_data,
  full_only_effect_data,
  noeq_only_effect_data
)

all_effect_data %>%
  distinct(variable, model)

# ============================================================
# LABELS
# ============================================================

variable_labels <- c(
  magnitude = "Magnitude",
  depth = "Focal depth",
  mean_pga = "Mean PGA",
  ntl_decline_ratio = "NTL decline ratio",
  gdp2015 = "GDP 2015",
  log_gdp = "Log GDP",
  population = "Population",
  dis_expen = "Disaster expenditure"
)

all_effect_data <- all_effect_data %>%
  mutate(
    variable_label = factor(
      variable_labels[variable],
      levels = variable_labels[
        c(
          "magnitude",
          "depth",
          "mean_pga",
          "ntl_decline_ratio",
          "gdp2015",
          "log_gdp",
          "population",
          "dis_expen"
        )
      ]
    )
  )

# ============================================================
# COMBINED CONDITIONAL-EFFECT PLOT
# ============================================================

combined_conditional_effects_plot <- ggplot(
  all_effect_data,
  aes(
    x = predictor_value,
    y = estimate,
    colour = model,
    fill = model
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ variable_label,
    scales = "free_x",
    ncol = 2
  ) +
  coord_cartesian(
    ylim = c(0, 500)
  ) +
  scale_colour_manual(
    values = c(
      "With magnitude/depth (M8)" = "blue",
      "Without magnitude/depth (M12)" = "red"
    )
  ) +
  scale_fill_manual(
    values = c(
      "With magnitude/depth (M8)" = "blue",
      "Without magnitude/depth (M12)" = "red"
    )
  ) +
  labs(
    title = "Conditional Effects of Covariates With and Without Earthquake Characteristics",
    x = "Standardised covariate value",
    y = "Expected recovery time (days)",
    colour = "Model",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

combined_conditional_effects_plot

ggsave(
  "combined_conditional_effects_plot.png",
  combined_conditional_effects_plot,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300
)

# ============================================================
# MODEL SUMMARIES
# ============================================================

summary(best_full_fit)
summary(best_noeq_fit)

# Compute LOO for each fitted model
loo_M8 <- loo(best_full_fit, moment_match = TRUE)

loo_M12 <- loo(best_noeq_fit, moment_match = TRUE)

# Direct pairwise comparison
loo_compare(list(M12 = loo_M12, M8  = loo_M8))
