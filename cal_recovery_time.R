# ============================================================
# PACKAGES
# ============================================================

library(dplyr)
library(ggplot2)
library(grid)
library(lubridate)
library(reshape2)
library(patchwork)
library(knitr)
library(kableExtra)

# ============================================================
# LOAD CLEANED DATASETS
# ============================================================

ntl = read.csv('/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL_clean.csv')
pga_clean = read.csv('/Users/nattanicha/Desktop/MSc_project/pga_clean_selected.csv')
econ_df = read.csv('/Users/nattanicha/Desktop/MSc_project/econ_df.csv')
eq = read.csv('/Users/nattanicha/Desktop/MSc_project/eq_df.csv')

colnames(ntl)
colnames(pga_clean)
colnames(econ_df)
colnames(eq)

hist(pga_clean$mean_pga)
hist(pga_clean$max_pga)

# ============================================================
# NTL VALIDATION AGAINST GDP
# ============================================================

econ_df$year

ntl = ntl %>%
  mutate(year = year(month))

ntl_yearly = ntl %>%
  group_by(prefecture, year) %>%
  summarise(
    mean_ntl = mean(mean_ntl, na.rm = TRUE),
    .groups = "drop"
  )

ntl_gdp_df = inner_join(
  econ_df,
  ntl_yearly,
  by = c("prefecture", "year")
)

summary(ntl_gdp_df)

cor.test(
  ntl_gdp_df$mean_ntl,
  ntl_gdp_df$gdp2015,
  method = "spearman"
) # rho 0.6598156, p-value < 2.2e-16

ggplot(
  ntl_gdp_df,
  aes(x = gdp2015, y = mean_ntl)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Relationship Between Nighttime Light Intensity and Gross Prefectural Product",
    x = "Gross Prefectural Product (2015 Base year, million yen, log scale)",
    y = "Mean Nighttime Light Intensity (log scale)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear on original scale
m1 = lm(gdp2015 ~ mean_ntl, data = ntl_gdp_df)

# Log-log
m2 = lm(log(gdp2015) ~ log(mean_ntl), data = ntl_gdp_df)

# Exponential
m3 = lm(log(gdp2015) ~ mean_ntl, data = ntl_gdp_df)

AIC(m1, m2, m3)

# AIC:
# m1  13063.7600
# m2    559.7216
# m3    621.8569

summary(m1) # Adjusted R-squared: 0.7658
summary(m2) # Adjusted R-squared: 0.6488
summary(m3) # Adjusted R-squared: 0.5857

# ============================================================
# PREPARE EARTHQUAKE AND NTL DATES
# ============================================================

ntl = ntl %>%
  mutate(date = as.Date(month))

pga_clean = pga_clean %>%
  mutate(date = as.Date(date))

# Order earthquakes within prefecture
eq_events = pga_clean %>%
  mutate(date = as.Date(date)) %>%
  arrange(prefecture, date) %>%
  group_by(prefecture) %>%
  mutate(eq_order = row_number()) %>%
  ungroup()

unique(eq_events$prefecture) # 35
max(eq_events$eq_order) # 10

ntl$month = as.Date(ntl$month)
eq_events$date = as.Date(eq_events$date)

# ============================================================
# RECOVERY STATE FUNCTION
# ============================================================

# Calculate recovery for one earthquake while tracking the NTL baseline.
#
# Baseline:
#   Mean NTL during the 3 months before the earthquake.
#
# Recovery:
#   First month when NTL reaches at least 95% of the baseline.
#
# Sequential earthquake rule:
#   If another earthquake occurs before recovery, the original baseline is
#   carried forward to the next earthquake. The baseline is reset only after
#   recovery is observed.
#
# Returns:
#   result        = recovery information for the current earthquake
#   next_baseline = baseline to carry forward, or NULL after recovery

calc_recovery_state <- function(
    eq_row,
    current_baseline = NULL,
    threshold = 0.95) {
  
  eq_id <- eq_row$eq_id
  pref <- eq_row$prefecture
  eq_date <- as.Date(eq_row$date)
  
  # NTL observations for the affected prefecture
  ntl_pref <- ntl %>%
    filter(prefecture == pref) %>%
    mutate(month = as.Date(month)) %>%
    arrange(month)
  
  baseline_type <- "new_baseline"
  
  # Use a new 3-month baseline unless one is carried forward
  if (is.null(current_baseline)) {
    baseline_months <- ntl_pref %>%
      filter(
        month < floor_date(eq_date, "month"),
        month >= floor_date(eq_date, "month") %m-% months(3)
      )
    
    baseline <- baseline_months %>%
      summarise(
        baseline = mean(mean_ntl, na.rm = TRUE)
      ) %>%
      pull(baseline)
  } else {
    baseline <- current_baseline
    baseline_type <- "carry_forward"
  }
  
  # Stop if no baseline can be calculated
  if (is.na(baseline)) {
    return(
      list(
        result = data.frame(
          eq_id = eq_id,
          prefecture = pref,
          date = eq_date,
          baseline = NA,
          baseline_type = baseline_type,
          recovery_days = NA,
          first_event = "no_baseline",
          ntl_decline_ratio = NA
        ),
        next_baseline = NULL
      )
    )
  }
  
  # First available NTL observation from the earthquake month onward
  first_after <- ntl_pref %>%
    filter(
      month >= floor_date(eq_date, "month")
    ) %>%
    arrange(month) %>%
    slice(1)
  
  # Stop if no post-earthquake NTL observation is available
  if (nrow(first_after) == 0) {
    return(
      list(
        result = data.frame(
          eq_id = eq_id,
          prefecture = pref,
          date = eq_date,
          baseline = baseline,
          baseline_type = baseline_type,
          recovery_days = NA,
          first_event = "no_ntl_after",
          ntl_decline_ratio = NA
        ),
        next_baseline = baseline
      )
    )
  }
  
  # Initial NTL decline relative to baseline
  ntl_decline_ratio = (baseline - first_after$mean_ntl) / baseline
  
  # Recovery threshold
  recovery_level <- threshold * baseline
  
  # Find the next earthquake in the same prefecture
  next_eq <- eq_events %>%
    filter(
      prefecture == pref,
      date > eq_date
    ) %>%
    arrange(date) %>%
    slice(1)
  
  # NTL observations after the current earthquake
  ntl_after <- ntl_pref %>%
    filter(
      month > floor_date(eq_date, "month")
    )
  
  # Stop searching for recovery when the next earthquake occurs
  if (nrow(next_eq) > 0) {
    ntl_after <- ntl_after %>%
      filter(
        month < floor_date(next_eq$date, "month")
      )
  }
  
  # First month reaching the recovery threshold
  recovery <- ntl_after %>%
    filter(
      mean_ntl >= recovery_level
    ) %>%
    arrange(month) %>%
    slice(1)
  
  # New earthquake occurs before recovery:
  # keep the same baseline for the next earthquake
  if (
    nrow(next_eq) > 0 &&
    (
      nrow(recovery) == 0 ||
      as.Date(next_eq$date) < as.Date(recovery$month)
    )
  ) {
    return(
      list(
        result = data.frame(
          eq_id = eq_id,
          prefecture = pref,
          date = eq_date,
          baseline = baseline,
          baseline_type = baseline_type,
          recovery_days = as.numeric(as.Date(next_eq$date) - eq_date),
          first_event = "new_earthquake",
          ntl_decline_ratio = ntl_decline_ratio
        ),
        next_baseline = baseline
      )
    )
  }
  
  # Recovery observed:
  # reset the baseline for the next independent earthquake
  if (nrow(recovery) > 0) {
    return(
      list(
        result = data.frame(
          eq_id = eq_id,
          prefecture = pref,
          date = eq_date,
          baseline = baseline,
          baseline_type = baseline_type,
          recovery_days = as.numeric(as.Date(recovery$month) - eq_date),
          first_event = "recovered",
          ntl_decline_ratio = ntl_decline_ratio
        ),
        next_baseline = NULL
      )
    )
  }
  
  # Neither recovery nor another earthquake is observed:
  # keep the baseline because recovery remains incomplete
  return(
    list(
      result = data.frame(
        eq_id = eq_id,
        prefecture = pref,
        date = eq_date,
        baseline = baseline,
        baseline_type = baseline_type,
        recovery_days = NA,
        first_event = "no_recovery_no_new_eq",
        ntl_decline_ratio = ntl_decline_ratio
      ),
      next_baseline = baseline
    )
  )
}

# ============================================================
# RUN RECOVERY ANALYSIS BY PREFECTURE
# ============================================================

results = list()

for (pref in unique(eq_events$prefecture)) {
  pref_events <- eq_events %>%
    filter(prefecture == pref) %>%
    arrange(date)
  
  current_baseline <- NULL
  
  for (i in 1:nrow(pref_events)) {
    out <- calc_recovery_state(
      pref_events[i, ],
      current_baseline = current_baseline,
      threshold = 0.95
    )
    
    results[[length(results) + 1]] <- out$result
    current_baseline <- out$next_baseline
  }
}

pga_recovery_state <- bind_rows(results)

pga_recovery_state
pga_recovery_state %>% arrange(prefecture, date)

table(pga_recovery_state$baseline_type)
table(pga_recovery_state$first_event)

write_csv(
  pga_recovery_state,
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_state.csv'
)

colnames(pga_recovery_state)

# ============================================================
# MERGE PGA VARIABLES
# ============================================================

pga_recovery_state = pga_recovery_state %>%
  left_join(
    pga_clean %>%
      dplyr::select(
        eq_id,
        prefecture,
        mean_pga,
        max_pga,
        magnitude,
        depth
      ),
    by = c("eq_id", "prefecture")
  )

pga_recovery_decline = pga_recovery_state %>%
  filter(ntl_decline_ratio > 0)

pga_recovery_decline

ggplot(
  pga_recovery_decline,
  aes(
    x = recovery_days,
    y = ntl_decline_ratio,
    colour = mean_pga
  )
) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    colour = "red",
    se = TRUE
  ) +
  scale_colour_viridis_c() +
  theme_minimal() +
  labs(
    title = "Relationship between NTL Decline and Recovery Time",
    subtitle = "Earthquakes with NTL decline ratio > 0 only",
    x = "Recovery time (days)",
    y = "NTL decline ratio",
    colour = "Mean PGA"
  )

# ============================================================
# MERGE ECONOMIC DATA
# ============================================================

pga_recovery_econ = pga_recovery_state %>%
  mutate(
    econ_year = year(date) - 1
  ) %>%
  left_join(
    econ_df,
    by = c(
      "prefecture" = "prefecture",
      "econ_year" = "year"
    )
  )

colnames(pga_recovery_econ)

summary(pga_recovery_econ$gdp2015)
summary(pga_recovery_econ$population)

colSums(
  is.na(
    pga_recovery_econ[, c(
      "population",
      "gdp2015",
      "dis_expen",
      "maintain_expen",
      "restore_expen"
    )]
  )
)

# ============================================================
# FINAL ANALYSIS DATASET
# ============================================================

pga_recovery_econ_clean = pga_recovery_econ %>%
  filter(
    !is.na(recovery_days),
    !is.na(ntl_decline_ratio),
    !is.na(mean_pga),
    !is.na(gdp2015),
    !is.na(population),
    !is.na(dis_expen),
    !is.na(maintain_expen),
    !is.na(restore_expen)
  )

pga_recovery_econ_clean

write_csv(
  pga_recovery_econ_clean,
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_clean_df.csv'
)

pga_recovery_econ_clean = read_csv(
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_clean_df.csv'
)

pga_recovery_econ_clean
nrow(pga_recovery_econ_clean)

pga_recovery_econ_clean1 = pga_recovery_econ_clean %>%
  filter(ntl_decline_ratio > 0)

nrow(pga_recovery_econ_clean1)

write_csv(
  pga_recovery_econ_clean1,
  '/Users/nattanicha/Desktop/MSc_project/pga_recovery_econ_onlydecline_df.csv'
)

colnames(pga_recovery_econ_clean1)

# ============================================================
# RECOVERY TIMELINE PLOT
# ============================================================

plot_data = pga_recovery_econ_clean1 %>%
  mutate(
    end_date = date + recovery_days,
    line_type = ifelse(
      first_event == "new_earthquake",
      "New earthquake",
      "Recovered"
    ),
    prefecture = factor(
      prefecture,
      levels = sort(unique(prefecture))
    ),
    line_type = factor(
      line_type,
      levels = c("Recovered", "New earthquake")
    )
  )

ggplot(
  plot_data,
  aes(
    x = date,
    xend = end_date,
    y = prefecture,
    yend = prefecture,
    colour = line_type,
    size = mean_pga
  )
) +
  geom_segment(lineend = "round") +
  scale_colour_manual(
    values = c(
      "New earthquake" = "#F8766D",
      "Recovered" = "#00BFC4"
    ),
    name = "Outcome"
  ) +
  scale_size_continuous(
    range = c(0.5, 3),
    name = "Mean PGA"
  ) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Earthquake Recovery Timeline by Prefecture",
    x = "Time period (Earthquake occurrence to recovery or subsequent earthquake)",
    y = "Prefecture"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

# ============================================================
# RECOVERY TIME DISTRIBUTION
# ============================================================

table(pga_recovery_econ_clean1$first_event)

pga_recovery_econ_clean2 = pga_recovery_econ_clean1 %>%
  filter(first_event == "recovered")

pga_recovery_econ_clean3 = pga_recovery_econ_clean1 %>%
  filter(first_event != "recovered")

ggplot(
  pga_recovery_econ_clean2,
  aes(x = recovery_days)
) +
  geom_histogram(
    bins = 20,
    fill = "#00BFC4",
    colour = "white"
  ) +
  labs(
    title = "Distribution of recovery time for recovered events",
    x = "Recovery time (days)",
    y = "Number of earthquakes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

pga_recovery_econ_clean1

table(pga_recovery_econ_clean1$first_event)
summary(pga_recovery_econ_clean2$recovery_days)
summary(pga_recovery_econ_clean3$recovery_days)

# ============================================================
# RECOVERY PLOTS
# ============================================================

p1 <- ggplot(
  pga_recovery_econ_clean2,
  aes(x = recovery_days)
) +
  geom_histogram(
    bins = 20,
    fill = "#00BFC4",
    colour = "white"
  ) +
  labs(
    title = "Distribution of Earthquake Recovery Time",
    x = "Recovery time (days)",
    y = "Number of earthquakes"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p2 <- ggplot(
  pga_recovery_econ_clean2,
  aes(
    x = recovery_days,
    y = ntl_decline_ratio,
    colour = mean_pga
  )
) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    colour = "red",
    se = TRUE
  ) +
  scale_colour_viridis_c() +
  theme_minimal(base_size = 10) +
  labs(
    title = "Relationship between NTL Decline and Recovery Time",
    subtitle = "Earthquakes with NTL decline ratio > 0 only",
    x = "Recovery time (days)",
    y = "NTL decline ratio",
    colour = "Mean PGA"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p1 + p2 + plot_layout(ncol = 2)

# ============================================================
# CORRELATION MATRIX
# ============================================================

cor_data = pga_recovery_econ_clean1 %>%
  dplyr::select(
    recovery_days,
    ntl_decline_ratio,
    mean_pga,
    max_pga,
    magnitude,
    depth,
    population,
    gdp2015,
    dis_expen,
    maintain_expen,
    restore_expen
  )

cor_matrix = cor(
  cor_data,
  use = "complete.obs"
)

cor_matrix

cor_long = melt(cor_matrix)

ggplot(
  cor_long,
  aes(
    x = Var1,
    y = Var2,
    fill = value
  )
) +
  geom_tile() +
  geom_text(
    aes(label = round(value, 2)),
    size = 3
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Correlation Matrix",
    x = "",
    y = "",
    fill = "Correlation"
  )

# ============================================================
# DESCRIPTIVE STATISTICS TABLE
# ============================================================

vars <- c(
  "recovery_days",
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

table1 <- pga_recovery_econ_clean1 %>%
  summarise(
    across(
      all_of(vars),
      list(
        Min = min,
        Mean = mean,
        Median = median,
        SD = sd,
        Max = max
      ),
      na.rm = TRUE
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "^(.*)_(Min|Mean|Median|SD|Max)$"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(., 2)
    )
  )

table1