# ============================================================
# OBSERVED INTERVALS VS POSTERIOR PREDICTIONS
# BEST BAYESIAN MODEL
# ============================================================

library(brms)
library(dplyr)
library(ggplot2)
library(grid)

# ============================================================
# IDENTIFY BEST MODEL
# ============================================================

best_model_name <- best_bayes_elpd_noeq$Model
best_family <- best_bayes_elpd_noeq$Family
best_structure <- best_bayes_elpd_noeq$Structure

best_model_name
best_family
best_structure

# ============================================================
# GET BEST FIT OBJECT
# ============================================================

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
# POSTERIOR PREDICTIONS
# ============================================================

set.seed(2026)

best_posterior_pred <- posterior_predict(
  best_fit,
  newdata = bayes_aft_data,
  re_formula = NULL
)

best_predicted_median <- apply(
  best_posterior_pred,
  2,
  median,
  na.rm = TRUE
)

best_predicted_lower <- apply(
  best_posterior_pred,
  2,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

best_predicted_upper <- apply(
  best_posterior_pred,
  2,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

# ============================================================
# CREATE PLOTTING DATA
# ============================================================

best_prediction_data <- bayes_aft_data %>%
  mutate(
    Observation = row_number(),
    Predicted_median = best_predicted_median,
    Predicted_lower = best_predicted_lower,
    Predicted_upper = best_predicted_upper,
    Censoring = case_when(
      is.finite(recovery_right) ~ "Interval-censored",
      is.infinite(recovery_right) ~ "Right-censored",
      TRUE ~ NA_character_
    ),
    Interval_error = case_when(
      is.finite(recovery_right) &
        Predicted_median < recovery_left ~
        recovery_left - Predicted_median,
      is.finite(recovery_right) &
        Predicted_median > recovery_right ~
        Predicted_median - recovery_right,
      is.finite(recovery_right) &
        Predicted_median >= recovery_left &
        Predicted_median <= recovery_right ~ 0,
      is.infinite(recovery_right) &
        Predicted_median < recovery_left ~
        recovery_left - Predicted_median,
      is.infinite(recovery_right) &
        Predicted_median >= recovery_left ~ 0,
      TRUE ~ NA_real_
    )
  )

# ============================================================
# MAE / RMSE FOR BEST MODEL
# ============================================================

best_prediction_performance <- best_prediction_data %>%
  summarise(
    Model = best_model_name,
    Family = best_family,
    Structure = best_structure,
    N = sum(is.finite(Interval_error)),
    MAE = mean(Interval_error, na.rm = TRUE),
    RMSE = sqrt(mean(Interval_error^2, na.rm = TRUE)),
    Median_Error = median(Interval_error, na.rm = TRUE)
  )

print(best_prediction_performance, width = Inf)



# ============================================================
# OBSERVED INTERVALS VS POSTERIOR PREDICTIONS
# EACH OBSERVATION ON ITS OWN LINE
# ============================================================

best_prediction_plot_data <- best_prediction_data %>%
  mutate(
    Sort_time = if_else(
      is.finite(recovery_right),
      (recovery_left + recovery_right) / 2,
      recovery_left
    )
  ) %>%
  arrange(Sort_time) %>%
  mutate(
    Plot_order = row_number(),
    Observation_label = Observation
  )

arrow_length <- 15

best_interval_prediction_plot <- ggplot(best_prediction_plot_data) +
  geom_hline(
    yintercept = seq(
      0.5,
      nrow(best_prediction_plot_data) + 0.5,
      by = 1
    ),
    linewidth = 0.25,
    colour = "grey85"
  ) +
  geom_segment(
    data = best_prediction_plot_data %>%
      filter(is.finite(recovery_right)),
    aes(
      x = recovery_left,
      xend = recovery_right,
      y = Plot_order,
      yend = Plot_order
    ),
    linewidth = 1.3
  ) +
  geom_segment(
    data = best_prediction_plot_data %>%
      filter(is.infinite(recovery_right)),
    aes(
      x = recovery_left,
      xend = recovery_left + arrow_length,
      y = Plot_order,
      yend = Plot_order
    ),
    linewidth = 1.1,
    linetype = "solid",
    arrow = arrow(
      length = unit(0.12, "cm"),
      type = "open"
    )
  ) +
  geom_point(
    aes(
      x = Predicted_median,
      y = Plot_order
    ),
    shape = 21,
    size = 2,
    stroke = 1,
    fill = "white"
  ) +
  scale_y_continuous(
    breaks = best_prediction_plot_data$Plot_order,
    labels = best_prediction_plot_data$Observation_label,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    title = "Observed recovery intervals and posterior predictions",
    subtitle = paste(best_model_name, best_family, best_structure),
    x = "Recovery time (days)",
    y = "Observation",
    caption = "Open circle = posterior median; solid line = observed interval; solid arrow = right-censored"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    axis.title = element_text(
      size = 13
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.ticks.y = element_blank(),
    plot.caption = element_text(
      size = 10,
      hjust = 0
    ),
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

best_interval_prediction_plot

ggsave(
  "best_bayesian_observed_vs_predicted.png",
  best_interval_prediction_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

# ============================================================
# PLOT 2
# POSTERIOR MEDIAN + 95% PREDICTIVE INTERVAL
# ============================================================

best_posterior_interval_plot <- best_prediction_data %>%
  ggplot(
    aes(
      x = reorder(
        as.factor(Observation),
        Predicted_median
      ),
      y = Predicted_median
    )
  ) +
  geom_errorbar(
    aes(
      ymin = Predicted_lower,
      ymax = Predicted_upper
    ),
    width = 0
  ) +
  geom_point(size = 2) +
  coord_flip() +
  labs(
    title = "Posterior predictive recovery times",
    subtitle = paste(
      best_model_name,
      best_family,
      best_structure
    ),
    x = "Observation",
    y = "Predicted recovery time (days)"
  ) +
  theme_bw()

best_posterior_interval_plot

# ============================================================
# PLOT 3
# INTERVAL-AWARE ABSOLUTE ERROR
# ============================================================

best_error_plot <- best_prediction_data %>%
  ggplot(
    aes(
      x = reorder(
        as.factor(Observation),
        Interval_error
      ),
      y = Interval_error
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Interval-aware prediction errors",
    subtitle = paste(
      best_model_name,
      best_family,
      best_structure
    ),
    x = "Observation",
    y = "Absolute error (days)"
  ) +
  theme_bw()

best_error_plot

# ============================================================
# PLOT 4
# PREDICTION VS OBSERVED INTERVAL MIDPOINT
# INTERVAL-CENSORED CASES ONLY
# ============================================================

best_prediction_midpoint_plot <- best_prediction_data %>%
  filter(is.finite(recovery_right)) %>%
  mutate(
    Observed_midpoint =
      (recovery_left + recovery_right) / 2
  ) %>%
  ggplot(
    aes(
      x = Observed_midpoint,
      y = Predicted_median
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(size = 2.3) +
  labs(
    title = "Posterior prediction versus observed recovery interval midpoint",
    subtitle = paste(
      best_model_name,
      best_family,
      best_structure
    ),
    x = "Observed interval midpoint (days)",
    y = "Posterior median predicted recovery time (days)"
  ) +
  theme_bw()

best_prediction_midpoint_plot

# ============================================================
# SAVE INTO RESULT OBJECT
# ============================================================

bayesian_aft_results$best_fit <- best_fit
bayesian_aft_results$best_prediction_data <- best_prediction_data
bayesian_aft_results$best_prediction_performance <- best_prediction_performance
bayesian_aft_results$best_interval_prediction_plot <- best_interval_prediction_plot
bayesian_aft_results$best_posterior_interval_plot <- best_posterior_interval_plot
bayesian_aft_results$best_error_plot <- best_error_plot
bayesian_aft_results$best_prediction_midpoint_plot <- best_prediction_midpoint_plot

saveRDS(
  bayesian_aft_results,
  file = "bayesian_interval_aft_loo_results.rds"
)
