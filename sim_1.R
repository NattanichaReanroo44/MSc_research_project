# ============================================================
# PACKAGES
# ============================================================

library(dplyr)
library(purrr)
library(tibble)


# ============================================================
# FIXED SIMULATION PARAMETERS
# ============================================================

simulation_parameters <- list(
  n = 61,
  n_prefectures = 24,
  n_earthquakes = 17,
  beta0 = 4.00,
  beta1 = 0.50,
  beta2 = 0.35,
  sigma = 0.8,
  sigma_prefecture = 0.40,
  sigma_earthquake = 0.50
)


# ============================================================
# GENERATE BASE RECOVERY DATA
# ============================================================

simulate_recovery_data <- function(
    n = 61,
    n_prefectures = 24,
    n_earthquakes = 17,
    beta0 = 4.00,
    beta1 = 0.50,
    beta2 = 0.35,
    sigma = 0.8,
    sigma_prefecture = 0.40,
    sigma_earthquake = 0.50) {
  
  # Assign prefecture and earthquake clusters
  prefecture_id <- sample(seq_len(n_prefectures), size = n, replace = TRUE)
  earthquake_id <- sample(seq_len(n_earthquakes), size = n, replace = TRUE)
  
  # Generate covariates
  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  
  # Generate cluster-level random effects
  prefecture_effects <- rnorm(n_prefectures, mean = 0, sd = sigma_prefecture)
  earthquake_effects <- rnorm(n_earthquakes, mean = 0, sd = sigma_earthquake)
  
  # Linear predictor for log recovery time
  mu <- beta0 + beta1 * X1 + beta2 * X2 +
    prefecture_effects[prefecture_id] +
    earthquake_effects[earthquake_id]
  
  # Generate lognormal recovery times
  true_recovery_time <- rlnorm(n, meanlog = mu, sdlog = sigma)
  
  tibble(
    observation_id = seq_len(n),
    prefecture_id = factor(prefecture_id),
    earthquake_id = factor(earthquake_id),
    X1 = X1,
    X2 = X2,
    prefecture_effect = prefecture_effects[prefecture_id],
    earthquake_effect = earthquake_effects[earthquake_id],
    mu = mu,
    true_recovery_time = true_recovery_time
  )
}

set.seed(2026)

# Generate 1,000 base recovery datasets
base_simulations <- map(
  seq_len(1000),
  ~ do.call(
    simulate_recovery_data,
    simulation_parameters
  ) %>%
    mutate(simulation_id = .x)
)

length(base_simulations) # 1000


# ============================================================
# SCENARIO 1: NON-INFORMATIVE CENSORING
# ============================================================

apply_noninformative_censoring <- function(recovery_data, gamma) {
  n <- nrow(recovery_data)
  
  censoring_time <- rexp(n, rate = gamma)
  
  recovery_data %>%
    mutate(
      censoring_time = censoring_time,
      observed_time = pmin(true_recovery_time, censoring_time),
      event = as.integer(true_recovery_time <= censoring_time),
      censored = 1L - event
    )
}


# ============================================================
# CALIBRATE NON-INFORMATIVE CENSORING
# ============================================================

estimate_censoring_proportion <- function(
    gamma,
    parameters,
    pilot_repetitions = 300) {
  
  censoring_proportions <- replicate(
    pilot_repetitions,
    {
      recovery_data <- do.call(
        simulate_recovery_data,
        parameters
      )
      
      censored_data <- apply_noninformative_censoring(
        recovery_data = recovery_data,
        gamma = gamma
      )
      
      mean(censored_data$censored)
    }
  )
  
  mean(censoring_proportions)
}

calibrate_gamma <- function(
    target_censoring,
    parameters,
    pilot_repetitions = 300,
    lower_gamma = 1e-6,
    upper_gamma = 1) {
  
  objective_function <- function(gamma) {
    estimated_censoring <- estimate_censoring_proportion(
      gamma = gamma,
      parameters = parameters,
      pilot_repetitions = pilot_repetitions
    )
    
    estimated_censoring - target_censoring
  }
  
  calibrated_result <- uniroot(
    objective_function,
    interval = c(lower_gamma, upper_gamma),
    tol = 1e-4
  )
  
  calibrated_result$root
}

set.seed(2026)

gamma_10 <- calibrate_gamma(
  target_censoring = 0.10,
  parameters = simulation_parameters
)

gamma_25 <- calibrate_gamma(
  target_censoring = 0.25,
  parameters = simulation_parameters
)

gamma_40 <- calibrate_gamma(
  target_censoring = 0.40,
  parameters = simulation_parameters
)

gamma_table <- tibble(
  target_censoring = c(0.10, 0.25, 0.40),
  gamma = c(gamma_10, gamma_25, gamma_40)
)

gamma_table


# ============================================================
# CHECK CENSORING CALIBRATION
# ============================================================

set.seed(2027)

calibration_check <- gamma_table %>%
  mutate(
    estimated_censoring = map_dbl(
      gamma,
      ~ estimate_censoring_proportion(
        gamma = .x,
        parameters = simulation_parameters,
        pilot_repetitions = 1000
      )
    )
  )

calibration_check


# ============================================================
# EXAMPLE NON-INFORMATIVE DATASET
# ============================================================

set.seed(2028)

example_recovery_data <- do.call(
  simulate_recovery_data,
  simulation_parameters
)

example_scenario1_data <- apply_noninformative_censoring(
  recovery_data = example_recovery_data,
  gamma = gamma_25
)

example_scenario1_data

example_scenario1_data %>%
  summarise(
    n = n(),
    n_events = sum(event),
    n_censored = sum(censored),
    censoring_proportion = mean(censored)
  ) # 0.197


# ============================================================
# GENERATE NON-INFORMATIVE CENSORING SCENARIOS
# ============================================================

scenario1_10 <- map(
  base_simulations,
  ~ apply_noninformative_censoring(
    recovery_data = .x,
    gamma = gamma_10
  )
)

scenario1_25 <- map(
  base_simulations,
  ~ apply_noninformative_censoring(
    recovery_data = .x,
    gamma = gamma_25
  )
)

scenario1_40 <- map(
  base_simulations,
  ~ apply_noninformative_censoring(
    recovery_data = .x,
    gamma = gamma_40
  )
)

length(scenario1_10) # 1000
length(scenario1_25)
length(scenario1_40)


# ============================================================
# NON-INFORMATIVE CENSORING SUMMARY
# ============================================================

censoring_summary <- tibble(
  Scenario = c("10%", "25%", "40%"),
  Mean_Censoring = c(
    mean(map_dbl(scenario1_10, ~ mean(.x$censored))), # 0.103
    mean(map_dbl(scenario1_25, ~ mean(.x$censored))), # 0.253
    mean(map_dbl(scenario1_40, ~ mean(.x$censored)))  # 0.401
  ),
  SD_Censoring = c(
    sd(map_dbl(scenario1_10, ~ mean(.x$censored))), # 0.0431
    sd(map_dbl(scenario1_25, ~ mean(.x$censored))), # 0.0597
    sd(map_dbl(scenario1_40, ~ mean(.x$censored)))  # 0.0719
  )
)

print(censoring_summary)


# ============================================================
# SCENARIO 2: INFORMATIVE CENSORING
# ============================================================

# Generate informative-censoring components once
generate_informative_components <- function(
    recovery_data,
    sigma_b = 0.30,
    sigma_c = 0.30,
    sigma_error = 0.60) {
  
  prefecture_index <- as.integer(recovery_data$prefecture_id)
  earthquake_index <- as.integer(recovery_data$earthquake_id)
  
  n_prefectures <- max(prefecture_index)
  n_earthquakes <- max(earthquake_index)
  n <- nrow(recovery_data)
  
  b_prefecture <- rnorm(n_prefectures, mean = 0, sd = sigma_b)
  c_earthquake <- rnorm(n_earthquakes, mean = 0, sd = sigma_c)
  epsilon <- rnorm(n, mean = 0, sd = sigma_error)
  
  tibble(
    b_prefecture = b_prefecture[prefecture_index],
    c_earthquake = c_earthquake[earthquake_index],
    censoring_error = epsilon
  )
}

set.seed(2031)

informative_components <- map(
  base_simulations,
  generate_informative_components,
  sigma_b = 0.30,
  sigma_c = 0.30,
  sigma_error = 0.60
)

length(informative_components) # 1000


# ============================================================
# APPLY INFORMATIVE CENSORING
# ============================================================

apply_informative_censoring <- function(
    recovery_data,
    censoring_components,
    alpha0,
    alpha1,
    alpha2) {
  
  log_censoring_time <-
    alpha0 +
    alpha1 * recovery_data$X1 +
    alpha2 * recovery_data$X2 +
    censoring_components$b_prefecture +
    censoring_components$c_earthquake +
    censoring_components$censoring_error
  
  censoring_time <- exp(log_censoring_time)
  
  recovery_data %>%
    mutate(
      censoring_time = censoring_time,
      observed_time = pmin(true_recovery_time, censoring_time),
      event = as.integer(true_recovery_time <= censoring_time),
      censored = 1L - event,
      censoring_scenario = "Informative"
    )
}


# ============================================================
# CALIBRATE INFORMATIVE CENSORING
# ============================================================

estimate_informative_censoring <- function(
    alpha0,
    alpha1,
    alpha2,
    base_data,
    components) {
  
  censoring_rates <- map2_dbl(
    base_data,
    components,
    function(recovery_data, censoring_components) {
      
      simulated_data <- apply_informative_censoring(
        recovery_data = recovery_data,
        censoring_components = censoring_components,
        alpha0 = alpha0,
        alpha1 = alpha1,
        alpha2 = alpha2
      )
      
      mean(simulated_data$censored)
    }
  )
  
  mean(censoring_rates)
}

calibrate_alpha0 <- function(
    alpha1,
    alpha2,
    target_censoring = 0.25,
    base_data,
    components,
    interval = c(-5, 10)) {
  
  objective_function <- function(alpha0) {
    achieved_censoring <- estimate_informative_censoring(
      alpha0 = alpha0,
      alpha1 = alpha1,
      alpha2 = alpha2,
      base_data = base_data,
      components = components
    )
    
    achieved_censoring - target_censoring
  }
  
  uniroot(
    objective_function,
    interval = interval,
    tol = 1e-5
  )$root
}


# ============================================================
# INFORMATIVE CENSORING SETTINGS
# ============================================================

informative_settings <- tribble(
  ~strength, ~alpha1, ~alpha2,
  "Weak", -0.20, -0.15,
  "Moderate", -0.40, -0.30,
  "Strong", -0.60, -0.45
)

informative_settings <- informative_settings %>%
  mutate(
    alpha0 = map2_dbl(
      alpha1,
      alpha2,
      ~ calibrate_alpha0(
        alpha1 = .x,
        alpha2 = .y,
        target_censoring = 0.25,
        base_data = base_simulations,
        components = informative_components
      )
    )
  )

informative_settings


# ============================================================
# GENERATE INFORMATIVE CENSORING SCENARIOS
# ============================================================

scenario2_informative <- pmap(
  informative_settings,
  function(strength, alpha1, alpha2, alpha0) {
    map2(
      base_simulations,
      informative_components,
      function(recovery_data, censoring_components) {
        apply_informative_censoring(
          recovery_data = recovery_data,
          censoring_components = censoring_components,
          alpha0 = alpha0,
          alpha1 = alpha1,
          alpha2 = alpha2
        ) %>%
          mutate(
            informative_strength = strength,
            alpha0 = alpha0,
            alpha1 = alpha1,
            alpha2 = alpha2
          )
      }
    )
  }
)

names(scenario2_informative) <- informative_settings$strength

length(scenario2_informative) # 3
length(scenario2_informative$Weak) # 1000
length(scenario2_informative$Moderate) # 1000
length(scenario2_informative$Strong) # 1000


# ============================================================
# INFORMATIVE CENSORING SUMMARY
# ============================================================

scenario2_censoring_summary <- imap_dfr(
  scenario2_informative,
  function(dataset_list, strength_name) {
    
    realised_rates <- map_dbl(
      dataset_list,
      ~ mean(.x$censored)
    )
    
    tibble(
      strength = strength_name,
      mean_censoring = mean(realised_rates),
      sd_censoring = sd(realised_rates),
      minimum_censoring = min(realised_rates),
      maximum_censoring = max(realised_rates)
    )
  }
)

scenario2_censoring_summary


# ============================================================
# SCENARIO 3: COMPETING RISKS
# ============================================================

generate_competing_components <- function(
    recovery_data,
    sigma_d = 0.30) {
  
  earthquake_index <- as.integer(recovery_data$earthquake_id)
  n_earthquakes <- max(earthquake_index)
  n <- nrow(recovery_data)
  
  d_earthquake <- rnorm(
    n_earthquakes,
    mean = 0,
    sd = sigma_d
  )
  
  standard_exponential <- rexp(
    n,
    rate = 1
  )
  
  tibble(
    d_earthquake = d_earthquake[earthquake_index],
    standard_exponential = standard_exponential
  )
}

set.seed(2032)

competing_components <- map(
  base_simulations,
  generate_competing_components,
  sigma_d = 0.30
)

length(competing_components) # 1000


# ============================================================
# APPLY COMPETING RISK
# ============================================================

apply_competing_risk <- function(
    recovery_data,
    competing_components,
    gamma0,
    gamma1 = 0.30) {
  
  log_competing_rate <-
    gamma0 +
    gamma1 * recovery_data$X1 +
    competing_components$d_earthquake
  
  competing_rate <- exp(log_competing_rate)
  
  competing_time <-
    competing_components$standard_exponential /
    competing_rate
  
  # 1 = recovery, 2 = subsequent earthquake
  event_type <- ifelse(
    recovery_data$true_recovery_time <= competing_time,
    1L,
    2L
  )
  
  recovery_data %>%
    mutate(
      competing_rate = competing_rate,
      competing_time = competing_time,
      observed_time = pmin(true_recovery_time, competing_time),
      event_type = event_type,
      recovery_event = as.integer(event_type == 1L),
      competing_event = as.integer(event_type == 2L),
      scenario = "Competing risks"
    )
}


# ============================================================
# CALIBRATE COMPETING-RISK PROPORTIONS
# ============================================================

estimate_competing_proportion <- function(
    gamma0,
    gamma1 = 0.30,
    base_data,
    components) {
  
  competing_proportions <- map2_dbl(
    base_data,
    components,
    function(recovery_data, competing_component) {
      
      simulated_data <- apply_competing_risk(
        recovery_data = recovery_data,
        competing_components = competing_component,
        gamma0 = gamma0,
        gamma1 = gamma1
      )
      
      mean(simulated_data$competing_event)
    }
  )
  
  mean(competing_proportions)
}

calibrate_gamma0 <- function(
    target_competing,
    gamma1 = 0.30,
    base_data,
    components,
    interval = c(-15, 5)) {
  
  objective_function <- function(gamma0) {
    achieved_competing <- estimate_competing_proportion(
      gamma0 = gamma0,
      gamma1 = gamma1,
      base_data = base_data,
      components = components
    )
    
    achieved_competing - target_competing
  }
  
  uniroot(
    objective_function,
    interval = interval,
    tol = 1e-5
  )$root
}


# ============================================================
# COMPETING-RISK SETTINGS
# ============================================================

competing_settings <- tibble(
  target_competing = c(0.10, 0.25, 0.40)
)

competing_settings <- competing_settings %>%
  mutate(
    gamma1 = 0.30,
    gamma0 = map_dbl(
      target_competing,
      ~ calibrate_gamma0(
        target_competing = .x,
        gamma1 = 0.30,
        base_data = base_simulations,
        components = competing_components
      )
    )
  )

competing_settings


# ============================================================
# GENERATE COMPETING-RISK SCENARIOS
# ============================================================

scenario3_competing <- pmap(
  competing_settings,
  function(target_competing, gamma1, gamma0) {
    map2(
      base_simulations,
      competing_components,
      function(recovery_data, competing_component) {
        apply_competing_risk(
          recovery_data = recovery_data,
          competing_components = competing_component,
          gamma0 = gamma0,
          gamma1 = gamma1
        ) %>%
          mutate(
            target_competing = target_competing,
            gamma0 = gamma0,
            gamma1 = gamma1
          )
      }
    )
  }
)

names(scenario3_competing) <- c(
  "competing_10",
  "competing_25",
  "competing_40"
)

length(scenario3_competing) # 3
length(scenario3_competing$competing_10) # 1000
length(scenario3_competing$competing_25) # 1000
length(scenario3_competing$competing_40) # 1000


# ============================================================
# COMPETING-RISK SUMMARY
# ============================================================

scenario3_competing_summary <- imap_dfr(
  scenario3_competing,
  function(dataset_list, setting_name) {
    
    realised_proportions <- map_dbl(
      dataset_list,
      ~ mean(.x$competing_event)
    )
    
    tibble(
      setting = setting_name,
      mean_competing_proportion = mean(realised_proportions),
      sd_competing_proportion = sd(realised_proportions),
      minimum_competing_proportion = min(realised_proportions),
      maximum_competing_proportion = max(realised_proportions)
    )
  }
)

scenario3_competing_summary


# ============================================================
# COMBINE ALL SIMULATION SCENARIOS
# ============================================================

all_simulation_scenarios <- list(
  base_recovery = base_simulations,
  noninformative = list(
    censoring_10 = scenario1_10,
    censoring_25 = scenario1_25,
    censoring_40 = scenario1_40
  ),
  informative = scenario2_informative,
  competing_risks = scenario3_competing
)


# ============================================================
# SAVE SIMULATION DATA
# ============================================================

saveRDS(
  all_simulation_scenarios,
  file = "all_simulation_scenarios.rds"
)