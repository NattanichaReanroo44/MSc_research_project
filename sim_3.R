# ============================================================
# GENERATE ALL SCENARIOS FOR A GIVEN SAMPLE SIZE
# ============================================================

generate_all_scenarios <- function(
    sample_size,
    n_simulations = 1000,
    n_prefectures = 24,
    n_earthquakes = 17,
    seed = 2026) {
  
  # Simulation parameters
  parameters <- list(
    n = sample_size,
    n_prefectures = n_prefectures,
    n_earthquakes = n_earthquakes,
    beta0 = 4.00,
    beta1 = 0.50,
    beta2 = 0.35,
    sigma = 0.80,
    sigma_prefecture = 0.40,
    sigma_earthquake = 0.50
  )
  
  # Generate common underlying recovery datasets
  set.seed(seed)
  
  base_data <- map(
    seq_len(n_simulations),
    function(simulation_id) {
      do.call(
        simulate_recovery_data,
        parameters
      ) %>%
        mutate(
          simulation_id = simulation_id,
          sample_size = sample_size
        )
    }
  )
  
  
  # ==========================================================
  # SCENARIO 1: NON-INFORMATIVE CENSORING
  # ==========================================================
  
  set.seed(seed + 1)
  
  gamma_values <- tibble(
    target_censoring = c(0.10, 0.25, 0.40)
  ) %>%
    mutate(
      gamma = map_dbl(
        target_censoring,
        function(target) {
          calibrate_gamma(
            target_censoring = target,
            parameters = parameters,
            pilot_repetitions = 300
          )
        }
      )
    )
  
  set.seed(seed + 2)
  
  noninformative <- map(
    gamma_values$gamma,
    function(current_gamma) {
      map(
        base_data,
        function(recovery_data) {
          apply_noninformative_censoring(
            recovery_data = recovery_data,
            gamma = current_gamma
          )
        }
      )
    }
  )
  
  names(noninformative) <- c(
    "censoring_10",
    "censoring_25",
    "censoring_40"
  )
  
  
  # ==========================================================
  # SCENARIO 2: INFORMATIVE CENSORING
  # ==========================================================
  
  set.seed(seed + 3)
  
  informative_components_n <- map(
    base_data,
    generate_informative_components,
    sigma_b = 0.30,
    sigma_c = 0.30,
    sigma_error = 0.60
  )
  
  informative_settings_n <- tribble(
    ~strength, ~alpha1, ~alpha2,
    "Weak", -0.20, -0.15,
    "Moderate", -0.40, -0.30,
    "Strong", -0.60, -0.45
  ) %>%
    mutate(
      alpha0 = map2_dbl(
        alpha1,
        alpha2,
        function(a1, a2) {
          calibrate_alpha0(
            alpha1 = a1,
            alpha2 = a2,
            target_censoring = 0.25,
            base_data = base_data,
            components = informative_components_n
          )
        }
      )
    )
  
  informative <- pmap(
    informative_settings_n,
    function(strength, alpha1, alpha2, alpha0) {
      map2(
        base_data,
        informative_components_n,
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
  
  names(informative) <- informative_settings_n$strength
  
  
  # ==========================================================
  # SCENARIO 3: COMPETING RISKS
  # ==========================================================
  
  set.seed(seed + 4)
  
  competing_components_n <- map(
    base_data,
    generate_competing_components,
    sigma_d = 0.30
  )
  
  competing_settings_n <- tibble(
    target_competing = c(0.10, 0.25, 0.40)
  ) %>%
    mutate(
      gamma1 = 0.30,
      gamma0 = map_dbl(
        target_competing,
        function(target) {
          calibrate_gamma0(
            target_competing = target,
            gamma1 = 0.30,
            base_data = base_data,
            components = competing_components_n
          )
        }
      )
    )
  
  competing_risks <- pmap(
    competing_settings_n,
    function(target_competing, gamma1, gamma0) {
      map2(
        base_data,
        competing_components_n,
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
  
  names(competing_risks) <- c(
    "competing_10",
    "competing_25",
    "competing_40"
  )
  
  
  # ==========================================================
  # RETURN SIMULATION OBJECT
  # ==========================================================
  
  list(
    sample_size = sample_size,
    parameters = parameters,
    
    calibration = list(
      gamma = gamma_values,
      informative = informative_settings_n,
      competing = competing_settings_n
    ),
    
    base_recovery = base_data,
    noninformative = noninformative,
    informative = informative,
    competing_risks = competing_risks
  )
}


# ============================================================
# GENERATE SAMPLE SIZE N = 120
# ============================================================

simulation_scenarios_n120 <- generate_all_scenarios(
  sample_size = 120,
  n_simulations = 1000,
  n_prefectures = 24,
  n_earthquakes = 17,
  seed = 1202026
)

saveRDS(
  simulation_scenarios_n120,
  file = "all_simulation_scenarios_n120.rds"
)


# ============================================================
# GENERATE SAMPLE SIZE N = 250
# ============================================================

simulation_scenarios_n250 <- generate_all_scenarios(
  sample_size = 250,
  n_simulations = 1000,
  n_prefectures = 24,
  n_earthquakes = 17,
  seed = 2502026
)

saveRDS(
  simulation_scenarios_n250,
  file = "all_simulation_scenarios_n250.rds"
)


# ============================================================
# RUN EXISTING ANALYSIS
# ============================================================

run_existing_analysis <- function(simulation_object, sample_size) {
  
  
  # ==========================================================
  # SCENARIO 1: NON-INFORMATIVE CENSORING
  # ==========================================================
  
  scenario1_results_10 <- run_standard_simulation_list(
    dataset_list = simulation_object$noninformative$censoring_10,
    scenario_name = "Non-informative 10%"
  )
  
  scenario1_results_25 <- run_standard_simulation_list(
    dataset_list = simulation_object$noninformative$censoring_25,
    scenario_name = "Non-informative 25%"
  )
  
  scenario1_results_40 <- run_standard_simulation_list(
    dataset_list = simulation_object$noninformative$censoring_40,
    scenario_name = "Non-informative 40%"
  )
  
  
  # ==========================================================
  # SCENARIO 2: INFORMATIVE CENSORING
  # ==========================================================
  
  scenario2_results_weak <- run_standard_simulation_list(
    dataset_list = simulation_object$informative$Weak,
    scenario_name = "Informative weak"
  )
  
  scenario2_results_moderate <- run_standard_simulation_list(
    dataset_list = simulation_object$informative$Moderate,
    scenario_name = "Informative moderate"
  )
  
  scenario2_results_strong <- run_standard_simulation_list(
    dataset_list = simulation_object$informative$Strong,
    scenario_name = "Informative strong"
  )
  
  
  # ==========================================================
  # SCENARIO 3: COMPETING RISKS
  # ==========================================================
  
  # Cox and AFT treat competing events as censoring
  scenario3_standard_10 <- run_standard_simulation_list(
    dataset_list = map(
      simulation_object$competing_risks$competing_10,
      prepare_competing_as_censoring
    ),
    scenario_name = "Competing risk 10%"
  )
  
  scenario3_standard_25 <- run_standard_simulation_list(
    dataset_list = map(
      simulation_object$competing_risks$competing_25,
      prepare_competing_as_censoring
    ),
    scenario_name = "Competing risk 25%"
  )
  
  scenario3_standard_40 <- run_standard_simulation_list(
    dataset_list = map(
      simulation_object$competing_risks$competing_40,
      prepare_competing_as_censoring
    ),
    scenario_name = "Competing risk 40%"
  )
  
  
  # ==========================================================
  # FINE-GRAY MODELS
  # ==========================================================
  
  fine_gray_10 <- run_fine_gray_list(
    dataset_list = simulation_object$competing_risks$competing_10,
    scenario_name = "Competing risk 10%"
  )
  
  fine_gray_25 <- run_fine_gray_list(
    dataset_list = simulation_object$competing_risks$competing_25,
    scenario_name = "Competing risk 25%"
  )
  
  fine_gray_40 <- run_fine_gray_list(
    dataset_list = simulation_object$competing_risks$competing_40,
    scenario_name = "Competing risk 40%"
  )
  
  
  # ==========================================================
  # COMBINE COEFFICIENT RESULTS
  # ==========================================================
  
  coefficient_results <- bind_rows(
    scenario1_results_10$coefficient_results,
    scenario1_results_25$coefficient_results,
    scenario1_results_40$coefficient_results,
    
    scenario2_results_weak$coefficient_results,
    scenario2_results_moderate$coefficient_results,
    scenario2_results_strong$coefficient_results,
    
    scenario3_standard_10$coefficient_results,
    scenario3_standard_25$coefficient_results,
    scenario3_standard_40$coefficient_results
  ) %>%
    mutate(
      sample_size = sample_size,
      .before = 1
    )
  
  
  # ==========================================================
  # COMBINE PREDICTION RESULTS
  # ==========================================================
  
  prediction_results <- bind_rows(
    scenario1_results_10$prediction_results,
    scenario1_results_25$prediction_results,
    scenario1_results_40$prediction_results,
    
    scenario2_results_weak$prediction_results,
    scenario2_results_moderate$prediction_results,
    scenario2_results_strong$prediction_results,
    
    scenario3_standard_10$prediction_results,
    scenario3_standard_25$prediction_results,
    scenario3_standard_40$prediction_results
  ) %>%
    mutate(
      sample_size = sample_size,
      .before = 1
    )
  
  
  # ==========================================================
  # COMBINE FINE-GRAY RESULTS
  # ==========================================================
  
  fine_gray_results <- bind_rows(
    fine_gray_10,
    fine_gray_25,
    fine_gray_40
  ) %>%
    mutate(
      sample_size = sample_size,
      .before = 1
    )
  
  
  # ==========================================================
  # COMBINE FAILURE RESULTS
  # ==========================================================
  
  failure_results <- bind_rows(
    scenario1_results_10$failures,
    scenario1_results_25$failures,
    scenario1_results_40$failures,
    
    scenario2_results_weak$failures,
    scenario2_results_moderate$failures,
    scenario2_results_strong$failures,
    
    scenario3_standard_10$failures,
    scenario3_standard_25$failures,
    scenario3_standard_40$failures
  ) %>%
    mutate(
      sample_size = sample_size,
      .before = 1
    )
  
  
  # ==========================================================
  # COEFFICIENT PERFORMANCE SUMMARY
  # ==========================================================
  
  coefficient_summary <- coefficient_results %>%
    group_by(
      sample_size,
      scenario,
      model,
      parameter
    ) %>%
    summarise(
      n_successful = sum(is.finite(estimate)),
      mean_estimate = mean(estimate, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      estimation_rmse = sqrt(
        mean(squared_error, na.rm = TRUE)
      ),
      coverage_95 = mean(covered, na.rm = TRUE),
      mean_standard_error = mean(
        standard_error,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # ==========================================================
  # PREDICTIVE PERFORMANCE SUMMARY
  # ==========================================================
  
  prediction_summary <- prediction_results %>%
    group_by(
      sample_size,
      scenario,
      model
    ) %>%
    summarise(
      n_successful = sum(is.finite(median_mae)),
      mean_mae = mean(median_mae, na.rm = TRUE),
      sd_mae = sd(median_mae, na.rm = TRUE),
      mean_prediction_rmse = mean(
        median_rmse,
        na.rm = TRUE
      ),
      sd_prediction_rmse = sd(
        median_rmse,
        na.rm = TRUE
      ),
      mean_observed_c_index = mean(
        observed_c_index,
        na.rm = TRUE
      ),
      sd_observed_c_index = sd(
        observed_c_index,
        na.rm = TRUE
      ),
      mean_true_c_index = mean(
        true_c_index,
        na.rm = TRUE
      ),
      sd_true_c_index = sd(
        true_c_index,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # ==========================================================
  # FINE-GRAY SUMMARY
  # ==========================================================
  
  fine_gray_summary <- fine_gray_results %>%
    filter(
      !is.na(parameter)
    ) %>%
    group_by(
      sample_size,
      scenario,
      parameter
    ) %>%
    summarise(
      n_successful = sum(is.finite(estimate)),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      mean_subdistribution_hazard_ratio = mean(
        subdistribution_hazard_ratio,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # ==========================================================
  # FAILURE SUMMARY
  # ==========================================================
  
  failure_summary <- failure_results %>%
    group_by(
      sample_size,
      scenario
    ) %>%
    summarise(
      successful_replicates = sum(successful),
      failed_replicates = sum(!successful),
      success_rate = mean(successful),
      .groups = "drop"
    )
  
  
  # ==========================================================
  # RETURN ANALYSIS RESULTS
  # ==========================================================
  
  list(
    coefficient_results = coefficient_results,
    coefficient_summary = coefficient_summary,
    prediction_results = prediction_results,
    prediction_summary = prediction_summary,
    fine_gray_results = fine_gray_results,
    fine_gray_summary = fine_gray_summary,
    failure_results = failure_results,
    failure_summary = failure_summary
  )
}


# ============================================================
# ANALYSIS FOR SAMPLE SIZE N = 120
# ============================================================

results_n120 <- run_existing_analysis(
  simulation_object = simulation_scenarios_n120,
  sample_size = 120
)

saveRDS(results_n120, file = "simulation_model_results_n120.rds")

print(results_n120$coefficient_summary, n = 36)

print(results_n120$prediction_summary, n = 27)

results_n120$fine_gray_summary


# ============================================================
# ANALYSIS FOR SAMPLE SIZE N = 250
# ============================================================

results_n250 <- run_existing_analysis(
  simulation_object = simulation_scenarios_n250,
  sample_size = 250
)

saveRDS(results_n250, file = "simulation_model_results_n250.rds")

print(results_n250$coefficient_summary, n = 36)

print(results_n250$prediction_summary, n = 27)

results_n250$fine_gray_summary