# ============================================================
# PACKAGES
# ============================================================

library(brms)
library(dplyr)
library(purrr)
library(ggplot2)


# ============================================================
# LOAD BAYESIAN AFT MODELS
# ============================================================

best_fit <- readRDS(
  "bayesian_aft_models/Lognormal_M13_Random.rds"
)

best_reduced_fit <- readRDS(
  "reduced_bayesian_aft_models/Lognormal_M9_Random_NoMagDepth.rds"
)


# ============================================================
# CONDITIONAL EFFECTS
# ============================================================

# Population-level conditional effects
full_effects <- conditional_effects(
  best_fit,
  method = "posterior_epred",
  re_formula = NA
)

reduced_effects <- conditional_effects(
  best_reduced_fit,
  method = "posterior_epred",
  re_formula = NA
)

names(full_effects)
names(reduced_effects)


# ============================================================
# EXTRACT CONDITIONAL EFFECTS
# ============================================================

extract_conditional_effect <- function(
    effect_list,
    variable,
    model_name) {
  
  if (!variable %in% names(effect_list)) {
    return(NULL)
  }
  
  effect_data <- effect_list[[variable]]
  
  tibble(
    predictor_value = effect_data[[variable]],
    estimate = effect_data$estimate__,
    lower = effect_data$lower__,
    upper = effect_data$upper__,
    variable = variable,
    model = model_name
  )
}


# ============================================================
# MODEL COVARIATES
# ============================================================

# Covariates included in both models
common_variables <- c(
  "mean_pga",
  "ntl_decline_ratio",
  "population",
  "dis_expen"
)

# Covariates included only in the full M13 model
full_only_variables <- c(
  "magnitude",
  "depth",
  "log_gdp"
)

# Covariate included only in the reduced M9 model
reduced_only_variables <- c(
  "gdp2015"
)


# ============================================================
# EXTRACT COMMON MODEL EFFECTS
# ============================================================

common_effect_data <- map_dfr(
  common_variables,
  function(variable) {
    bind_rows(
      extract_conditional_effect(
        effect_list = full_effects,
        variable = variable,
        model_name = "Full M13"
      ),
      extract_conditional_effect(
        effect_list = reduced_effects,
        variable = variable,
        model_name = "Reduced M9"
      )
    )
  }
)


# ============================================================
# EXTRACT MODEL-SPECIFIC EFFECTS
# ============================================================

full_only_effect_data <- map_dfr(
  full_only_variables,
  function(variable) {
    extract_conditional_effect(
      effect_list = full_effects,
      variable = variable,
      model_name = "Full M13"
    )
  }
)

reduced_only_effect_data <- map_dfr(
  reduced_only_variables,
  function(variable) {
    extract_conditional_effect(
      effect_list = reduced_effects,
      variable = variable,
      model_name = "Reduced M9"
    )
  }
)


# ============================================================
# COMBINE CONDITIONAL EFFECTS
# ============================================================

all_effect_data <- bind_rows(
  common_effect_data,
  full_only_effect_data,
  reduced_only_effect_data
)

all_effect_data %>%
  distinct(variable, model)


# ============================================================
# VARIABLE LABELS
# ============================================================

variable_labels <- c(
  magnitude = "Magnitude",
  depth = "Focal depth",
  mean_pga = "Mean PGA",
  ntl_decline_ratio = "NTL decline ratio",
  log_gdp = "Log GDP",
  gdp2015 = "GDP 2015",
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
          "log_gdp",
          "gdp2015",
          "population",
          "dis_expen"
        )
      ]
    )
  )


# ============================================================
# CONDITIONAL EFFECTS PLOT
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
  geom_line(
    linewidth = 1
  ) +
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
      "Full M13" = "blue",
      "Reduced M9" = "red"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Full M13" = "blue",
      "Reduced M9" = "red"
    )
  ) +
  labs(
    title = "Conditional Effects of the Full and Reduced Bayesian AFT Models",
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
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

combined_conditional_effects_plot


# ============================================================
# MODEL SUMMARIES
# ============================================================

summary(best_fit)
summary(best_reduced_fit)