# ============================================================
# PACKAGES
# ============================================================

library(dplyr)
library(ggplot2)
library(survival)
library(splines)
library(e1071)
library(purrr)
library(broom)
library(patchwork)

# ============================================================
# DATA PREPARATION
# ============================================================

pga_recovery_econ_clean1 = read_csv(
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv'
)

pga_recovery_econ_clean2 = pga_recovery_econ_clean1 %>%
  filter(first_event == "recovered")

data_table <- pga_recovery_econ_clean1 %>%
  mutate(
    population = population / 1e6, # million persons
    gdp2015 = gdp2015 / 1e6, # trillion yen
    dis_expen = dis_expen / 1e6, # billion yen
    maintain_expen = maintain_expen / 1e6,
    restore_expen = restore_expen / 1e6
  )

data_table

data <- data_table %>%
  mutate(
    recovery_group = ifelse(
      recovery_days < 100,
      "Short (<150 days)",
      "Long (>=150 days)"
    ),
    log_gdp = log(gdp2015),
    log_population = log(population),
    log_dis_expen = log(dis_expen),
    log_maintain_expen = log(maintain_expen),
    log_restore_expen = log(restore_expen)
  )

data$recovery_group <- factor(
  data$recovery_group,
  levels = c(
    "Short (<150 days)",
    "Long (>=150 days)"
  )
)

# ============================================================
# RECOVERY TIME DISTRIBUTION
# ============================================================

ggplot(
  data,
  aes(recovery_days)
) +
  geom_histogram(bins = 30) +
  geom_vline(
    xintercept = 150,
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    x = "Recovery days",
    y = "Frequency"
  )

table(data$first_event)

# ============================================================
# ECONOMIC VARIABLES
# ============================================================

eco_vars <- c(
  "gdp2015",
  "population",
  "dis_expen",
  "maintain_expen",
  "restore_expen"
)

# ============================================================
# NORMALITY TEST
# ============================================================

normality <- map_df(
  eco_vars,
  ~ {
    test <- shapiro.test(data[[.x]])
    
    data.frame(
      variable = .x,
      W = test$statistic,
      p_value = test$p.value
    )
  }
)

normality

# All p-values < 0.05, indicating evidence against normality

# ============================================================
# WILCOXON RANK-SUM TEST
# ============================================================

wilcox_results <- map_df(
  eco_vars,
  ~ {
    test <- wilcox.test(
      data[[.x]] ~ data$recovery_group
    )
    
    data.frame(
      variable = .x,
      p_value = test$p.value
    )
  }
)

wilcox_results

# ============================================================
# BOXPLOTS
# ============================================================

plots <- lapply(
  eco_vars,
  function(v) {
    ggplot(
      data,
      aes(
        x = recovery_group,
        y = .data[[v]]
      )
    ) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        x = "Recovery group",
        y = v
      )
  }
)

wrap_plots(
  plots,
  ncol = 3
)

# ============================================================
# BOXPLOTS ON LOG SCALE
# ============================================================

plots <- lapply(
  eco_vars,
  function(v) {
    ggplot(
      data,
      aes(
        x = recovery_group,
        y = .data[[v]]
      )
    ) +
      geom_boxplot() +
      scale_y_log10() +
      theme_minimal() +
      labs(
        x = "Recovery group",
        y = paste0(v, " (log scale)")
      )
  }
)

wrap_plots(
  plots,
  ncol = 3
)

# ============================================================
# WILCOXON TEST ON LOG-TRANSFORMED VARIABLES
# ============================================================

wilcox_results_log <- lapply(
  eco_vars,
  function(v) {
    test <- wilcox.test(
      log(data[[v]]) ~ data$recovery_group
    )
    
    data.frame(
      variable = v,
      p_value = test$p.value
    )
  }
) %>%
  bind_rows()

wilcox_results_log