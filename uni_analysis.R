# ============================================================
# PACKAGES
# ============================================================

library(survival)
library(survminer)
library(dplyr)
library(splines)
library(broom)
library(purrr)

# ============================================================
# LOAD DATA
# ============================================================

ntl = read.csv('/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL_clean.csv')
pga_clean = read.csv('/Users/nattanicha/Desktop/MSc_project/pga_clean_selected.csv')
econ_df = read.csv('/Users/nattanicha/Desktop/MSc_project/econ_df.csv')
eq = read.csv('/Users/nattanicha/Desktop/MSc_project/eq_df.csv')
pga_recovery = read.csv('/Users/nattanicha/Desktop/MSc_project/recovery_df.csv')
pga_recovery_econ_clean1 = read_csv('/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv')

# ============================================================
# DATA PREPARATION
# ============================================================

data <- pga_recovery_econ_clean1

data_table <- data %>%
  mutate(
    population = population / 1e6, # million persons
    gdp2015 = gdp2015 / 1e6, # trillion yen
    dis_expen = dis_expen / 1e6, # billion yen
    maintain_expen = maintain_expen / 1e6,
    restore_expen = restore_expen / 1e6
  )

data_table

km_data = data_table %>%
  mutate(
    event = ifelse(first_event == "recovered", 1, 0)
  )

# ============================================================
# OVERALL KAPLAN-MEIER CURVE
# ============================================================

km_fit = survfit(
  Surv(recovery_days, event) ~ 1,
  data = km_data
)

p1 = ggsurvplot(
  km_fit,
  data = km_data,
  conf.int = TRUE,
  risk.table = TRUE,
  xlab = "Days after earthquake",
  ylab = "Probability of not yet recovered",
  title = "Kaplan–Meier Recovery Curve of Nighttime Light After Earthquakes",
  surv.scale = "percent",
  ggtheme = theme_minimal(base_size = 10)
)

# ============================================================
# KAPLAN-MEIER CURVES BY PGA GROUP
# ============================================================

km_data = data_table %>%
  filter(
    !is.na(mean_pga),
    !is.na(recovery_days)
  ) %>%
  mutate(
    event = ifelse(first_event == "recovered", 1, 0),
    pga_group = ntile(mean_pga, 2),
    pga_group = factor(
      pga_group,
      levels = c(1, 2),
      labels = c("Low PGA", "High PGA")
    )
  )

km_fit = survfit(
  Surv(recovery_days, event) ~ pga_group,
  data = km_data
)

p2 = ggsurvplot(
  km_fit,
  data = km_data,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Days after earthquake",
  ylab = "Probability of not yet recovered",
  title = "Kaplan–Meier Recovery Curves Stratified by Prefecture-Level Average PGA",
  legend.title = "mean PGA group",
  legend.labs = c("Low (1.22-3.50)", "High (3.69-30)"),
  surv.scale = "percent",
  ggtheme = theme_minimal(base_size = 10)
)

km_data %>%
  group_by(pga_group) %>%
  summarise(
    n = n(),
    min_pga = min(mean_pga),
    max_pga = max(mean_pga),
    median_pga = median(mean_pga)
  )

arrange_ggsurvplots(
  list(p1, p2),
  ncol = 2,
  nrow = 1
)

# ============================================================
# UNIVARIATE ANALYSIS
# ============================================================

surv_obj <- Surv(
  km_data$recovery_days,
  km_data$event
)

economic_vars <- c(
  "population",
  "gdp2015",
  "dis_expen",
  "maintain_expen",
  "restore_expen"
)

km_data %>%
  summarise(
    across(
      all_of(economic_vars),
      list(
        minimum = ~ min(.x, na.rm = TRUE),
        zeros = ~ sum(.x == 0, na.rm = TRUE)
      )
    )
  )

vars <- c(
  "ntl_decline_ratio",
  "mean_pga",
  "max_pga",
  "magnitude",
  "depth",
  "population",
  "gdp2015",
  "dis_expen",
  "maintain_expen",
  "restore_expen"
)

# Create log-transformed variables
km_data <- km_data %>%
  mutate(
    across(
      all_of(vars),
      log,
      .names = "log_{.col}"
    )
  )

# ============================================================
# UNIVARIATE COX MODELS
# ============================================================

cox_uni <- map_dfr(
  vars,
  function(x) {
    model <- coxph(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~",
          x
        )
      ),
      data = km_data
    )
    
    tidy(model) %>%
      mutate(variable = x)
  }
)

cox_uni

# ============================================================
# UNIVARIATE AFT MODELS
# ============================================================

aft_weibull <- map_dfr(
  vars,
  function(x) {
    model <- survreg(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~",
          x
        )
      ),
      data = km_data,
      dist = "weibull"
    )
    
    tidy(model) %>%
      mutate(
        variable = x,
        distribution = "Weibull"
      )
  }
)

aft_weibull

aft_lognormal <- map_dfr(
  vars,
  function(x) {
    model <- survreg(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~",
          x
        )
      ),
      data = km_data,
      dist = "lognormal"
    )
    
    tidy(model) %>%
      mutate(
        variable = x,
        distribution = "Log-normal"
      )
  }
)

aft_lognormal

aft_loglogistic <- map_dfr(
  vars,
  function(x) {
    model <- survreg(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~",
          x
        )
      ),
      data = km_data,
      dist = "loglogistic"
    )
    
    tidy(model) %>%
      mutate(
        variable = x,
        distribution = "Log-logistic"
      )
  }
)

aft_loglogistic

# Combine all AFT distributions
aft_uni <- map_dfr(
  vars,
  function(x) {
    map_dfr(
      c("weibull", "lognormal", "loglogistic"),
      function(dist) {
        model <- survreg(
          as.formula(
            paste(
              "Surv(recovery_days,event) ~",
              x
            )
          ),
          data = km_data,
          dist = dist
        )
        
        tidy(model) %>%
          filter(term == x) %>%
          mutate(
            variable = x,
            distribution = dist
          )
      }
    )
  }
)

print(aft_uni, n = 30)

# ============================================================
# NONLINEARITY TEST: COX MODELS
# ============================================================

nonlinear_test <- map_dfr(
  vars,
  function(x) {
    m1 <- coxph(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~",
          x
        )
      ),
      data = km_data
    )
    
    m2 <- coxph(
      as.formula(
        paste(
          "Surv(recovery_days,event) ~ ns(",
          x,
          ",df=3)"
        )
      ),
      data = km_data
    )
    
    test <- anova(
      m1,
      m2,
      test = "LRT"
    )
    
    data.frame(
      variable = x,
      Chi_square = test$Chisq[2],
      df = test$Df[2],
      p_value = test$`Pr(>|Chi|)`[2]
    )
  }
)

nonlinear_test

# ============================================================
# NONLINEARITY TEST: AFT MODELS
# ============================================================

dists <- c(
  "weibull",
  "lognormal",
  "loglogistic"
)

aft_nonlinear_test <- map_dfr(
  dists,
  function(dist) {
    map_dfr(
      vars,
      function(x) {
        linear <- survreg(
          as.formula(
            paste(
              "Surv(recovery_days,event) ~",
              x
            )
          ),
          data = km_data,
          dist = dist
        )
        
        nonlinear <- survreg(
          as.formula(
            paste(
              "Surv(recovery_days,event) ~ ns(",
              x,
              ",df=3)"
            )
          ),
          data = km_data,
          dist = dist
        )
        
        # Likelihood ratio test
        LR <- 2 * (
          as.numeric(logLik(nonlinear)) -
            as.numeric(logLik(linear))
        )
        
        df_diff <- attr(logLik(nonlinear), "df") -
          attr(logLik(linear), "df")
        
        p_value <- pchisq(
          LR,
          df = df_diff,
          lower.tail = FALSE
        )
        
        data.frame(
          distribution = dist,
          variable = x,
          Chi_square = LR,
          df = df_diff,
          p_value = p_value
        )
      }
    )
  }
)

aft_nonlinear_test

# ============================================================
# COX RESULTS TABLE
# ============================================================

options(scipen = 0)

cox_uni %>%
  mutate(
    estimate = formatC(
      estimate,
      format = "e",
      digits = 2
    ),
    std.error = formatC(
      std.error,
      format = "e",
      digits = 2
    )
  )

cox_report <- cox_uni %>%
  mutate(
    HR = exp(estimate),
    HR_p = paste0(
      round(HR, 3),
      " (p=",
      round(p.value, 3),
      ")"
    )
  ) %>%
  dplyr::select(
    variable,
    HR_p
  ) %>%
  left_join(
    nonlinear_test %>%
      dplyr::select(
        variable,
        p_value
      ),
    by = "variable"
  )

cox_report

# ============================================================
# AFT RESULTS TABLE
# ============================================================

sci_format <- function(x) {
  ifelse(
    abs(x) < 0.001,
    formatC(
      x,
      format = "e",
      digits = 2
    ),
    round(x, 3)
  )
}

aft_report <- aft_uni %>%
  filter(
    !term %in% c(
      "(Intercept)",
      "Log(scale)"
    )
  ) %>%
  mutate(
    beta_p = paste0(
      sci_format(estimate),
      " (p=",
      round(p.value, 3),
      ")"
    )
  ) %>%
  dplyr::select(
    variable,
    distribution,
    beta_p
  ) %>%
  pivot_wider(
    names_from = distribution,
    values_from = beta_p
  ) %>%
  left_join(
    aft_nonlinear_test %>%
      filter(
        distribution == "weibull"
      ) %>%
      dplyr::select(
        variable,
        p_value
      ),
    by = "variable"
  )

aft_report

aft_report %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Univariate AFT model results and nonlinear tests"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

cox_report %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Univariate Cox model results and nonlinear tests"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

# ============================================================
# COMPARE LINEAR AND LOG FORMS: COX
# ============================================================

compare_cox_forms <- function(x, data) {
  log_x <- paste0("log_", x)
  
  model_linear <- coxph(
    as.formula(
      paste(
        "Surv(recovery_days, event) ~",
        x
      )
    ),
    data = data,
    na.action = na.omit
  )
  
  model_log <- coxph(
    as.formula(
      paste(
        "Surv(recovery_days, event) ~",
        log_x
      )
    ),
    data = data,
    na.action = na.omit
  )
  
  tibble(
    variable = x,
    linear_AIC = AIC(model_linear),
    log_AIC = AIC(model_log),
    delta_AIC = AIC(model_linear) - AIC(model_log),
    preferred_form = case_when(
      AIC(model_log) + 2 < AIC(model_linear) ~ "Log",
      AIC(model_linear) + 2 < AIC(model_log) ~ "Linear",
      TRUE ~ "Similar"
    )
  )
}

cox_form_comparison <- map_dfr(
  vars,
  compare_cox_forms,
  data = km_data
)

cox_form_comparison

# ============================================================
# COMPARE LINEAR AND LOG FORMS: AFT
# ============================================================

compare_aft_forms <- function(
    x,
    data,
    dist) {
  
  log_x <- paste0(
    "log_",
    x
  )
  
  model_linear <- survreg(
    as.formula(
      paste(
        "Surv(recovery_days, event) ~",
        x
      )
    ),
    data = data,
    dist = dist,
    na.action = na.omit
  )
  
  model_log <- survreg(
    as.formula(
      paste(
        "Surv(recovery_days, event) ~",
        log_x
      )
    ),
    data = data,
    dist = dist,
    na.action = na.omit
  )
  
  tibble(
    variable = x,
    distribution = dist,
    linear_AIC = AIC(model_linear),
    log_AIC = AIC(model_log),
    delta_AIC = AIC(model_linear) - AIC(model_log),
    preferred_form = case_when(
      AIC(model_log) + 2 < AIC(model_linear) ~ "Log",
      AIC(model_linear) + 2 < AIC(model_log) ~ "Linear",
      TRUE ~ "Similar"
    )
  )
}

aft_form_comparison <- map_dfr(
  c(
    "weibull",
    "lognormal",
    "loglogistic"
  ),
  function(dist) {
    map_dfr(
      vars,
      compare_aft_forms,
      data = km_data,
      dist = dist
    )
  }
)

print(
  aft_form_comparison,
  n = 30
)

# ============================================================
# PROPORTIONAL HAZARDS ASSUMPTION
# ============================================================

ph_results <- lapply(
  vars,
  function(v) {
    fit <- coxph(
      as.formula(
        paste0(
          "Surv(recovery_days, event) ~ ",
          v
        )
      ),
      data = km_data
    )
    
    zph <- cox.zph(fit)
    
    data.frame(
      Variable = v,
      PH_p = zph$table[1, "p"]
    )
  }
)

do.call(
  rbind,
  ph_results
)