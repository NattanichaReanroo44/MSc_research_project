# ============================================================
# PACKAGES
# ============================================================

library(tidyverse)
library(readxl)
library(sf)
library(terra)
library(raster)
library(exactextractr)
library(lubridate)
library(dplyr)

# ============================================================
# EARTHQUAKE DATA PRE-PROCESSING
# ============================================================

raw = read_excel('/Users/nattanicha/Desktop/MSc_project/6up.xlsx')

eq_clean = raw %>%
  dplyr::select(id, time, latitude, longitude, depth, mag, place)

eq_clean = eq_clean %>%
  rename(
    eq_id = id,
    date = time,
    lat = latitude,
    lon = longitude,
    magnitude = mag
  ) %>%
  arrange(date)

write_csv(eq_clean, '/Users/nattanicha/Desktop/MSc_project/eq_catalog_6up_clean.csv')

hist(eq_clean$magnitude)

# Manually cleaned according to available PGA raster files
eq = read_excel('/Users/nattanicha/Desktop/MSc_project/eq_catalog6.xlsx')

colnames(eq)

eq$date = as.Date(eq$date)
eq$year = format(eq$date, "%Y")

write_csv(eq, '/Users/nattanicha/Desktop/MSc_project/eq_df.csv')

# ============================================================
# PREFECTURE SHAPEFILE
# ============================================================

jp_pref = st_read('/Users/nattanicha/Desktop/MSc_project/GEE/jp_prefecture.shp')

eq_id = eq$eq_id # earthquake IDs

# ============================================================
# ECONOMIC DATA PRE-PROCESSING
# ============================================================

econ_data = read_excel('/Users/nattanicha/Desktop/MSc_project/econ_proxy.xlsx')

econ_data_clean = econ_data %>%
  mutate(
    prefecture = AREA %>%
      str_replace("-ken", "") %>%
      str_replace("-fu", "") %>%
      str_replace("-to", "") %>%
      str_replace("-do", "")
  )

unique(econ_data_clean$prefecture)

econ_data_clean$prefecture = gsub(
  "Gumma",
  "Gunma",
  econ_data_clean$prefecture
)

colnames(econ_data_clean)

econ_df = econ_data_clean %>%
  dplyr::select(
    year = YEAR,
    prefecture = prefecture,
    population = `A1101_Total population (Both sexes)[person]`,
    gdp2015 = `C1121_Gross prefectural product (2015 base)[million yen]`,
    dis_expen = `D310312_Disaster relief expenditure (Prefecture)[thousand yen]`,
    maintain_expen = `D310403_Expenditure for maintenance and repairs (Prefecture)[thousand yen]`,
    restore_expen = `D310407_Expenditure for disaster restoration (Prefecture)[thousand yen]`
  )

colnames(econ_df)

econ_df = econ_df %>%
  mutate(across(where(is.character), ~ na_if(., "***"))) %>%
  mutate(
    across(
      c(population, gdp2015, dis_expen, maintain_expen, restore_expen),
      as.numeric
    )
  )

write_csv(econ_df, '/Users/nattanicha/Desktop/MSc_project/econ_df.csv')

# ============================================================
# PGA EXTRACTION BY PREFECTURE
# ============================================================

all_results = list()


# Calculate mean and maximum PGA for each prefecture
for (i in eq_id) {
  cat('processing eq id:', i)
  
  path = paste0(
    '/Users/nattanicha/Desktop/MSc_project/shakemap_PGA/',
    i,
    '/pga_mean.flt'
  )
  
  if (!file.exists(path)) {
    cat('missing raster file for eq id:', i, '\n')
    next
  }
  
  pga = rast(path)
  jp_pref_ = st_transform(jp_pref, st_crs(pga))
  extract_val = exact_extract(pga, jp_pref_)
  
  all_val = values(pga)
  all_val = all_val[!is.na(all_val)]
  
  is_log = min(all_val) < 0
  cat(i, 'log scale:', is_log, '\n')
  
  if (is_log) {
    pga = exp(pga) * 100
  }
  
  pga_r = raster(pga)
  
  mean_val = exact_extract(pga, jp_pref, "mean")
  max_val = exact_extract(pga, jp_pref, "max")
  
  result = data.frame(
    eq_id = i,
    prefecture = jp_pref_$nam,
    mean_pga = mean_val,
    max_pga = max_val
  )
  
  all_results[[i]] = result
}

pga_all = bind_rows(all_results)

head(pga_all)

length(unique(pga_all$eq_id)) # 75
nrow(pga_all) # 47 * 75 = 3525

summary(pga_all$mean_pga)
summary(pga_all$max_pga)

# ============================================================
# CLEAN PGA DATA
# ============================================================

pga_all_clean = pga_all %>%
  mutate(
    prefecture = prefecture %>%
      str_replace(" Ken", "") %>%
      str_replace(" Fu", "") %>%
      str_replace(" To", "") %>%
      str_replace(" Do", "")
  )

unique(pga_all_clean$prefecture)

pga_all_clean$prefecture = gsub(
  "Hokkai",
  "Hokkaido",
  pga_all_clean$prefecture
)

pga_all_clean

write.csv(
  pga_all_clean,
  '/Users/nattanicha/Desktop/MSc_project/pga_all.csv',
  row.names = F
)

pga_all_clean = read.csv('/Users/nattanicha/Desktop/MSc_project/pga_all.csv')

pga_count = pga_all_clean %>%
  group_by(eq_id) %>%
  summarize(
    n_prefecture = sum(mean_pga > 0, na.rm = TRUE)
  )

pga_count

write.csv(
  pga_count,
  '/Users/nattanicha/Desktop/MSc_project/pga_count.csv',
  row.names = F
)

never_affected = pga_all_clean %>%
  group_by(prefecture) %>%
  summarize(
    ever_affected = any(mean_pga > 0, na.rm = TRUE)
  ) %>%
  filter(!ever_affected)

never_affected # 0 rows

pga_all_clean = pga_all_clean %>%
  left_join(
    eq %>%
      dplyr::select(eq_id, date, depth, magnitude),
    by = "eq_id"
  )

# pga_clean = pga_all_clean %>% filter(mean_pga > 1 | max_pga > 5)

pga_clean = pga_all_clean %>%
  filter(max_pga > 5)

nrow(pga_clean)

write_csv(
  pga_clean,
  '/Users/nattanicha/Desktop/MSc_project/pga_clean_selected.csv'
)

# ============================================================
# NIGHTTIME LIGHT DATA
# ============================================================

ntl = read_csv('/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL.csv')

head(ntl)

ntl_clean = ntl %>%
  dplyr::select(
    prefecture = nam,
    month,
    mean_ntl = mean
  ) %>%
  mutate(
    prefecture = prefecture %>%
      str_replace(" Ken", "") %>%
      str_replace(" Fu", "") %>%
      str_replace(" To", "") %>%
      str_replace(" Do", ""),
    month = ym(month)
  ) %>%
  arrange(prefecture, month)

head(ntl_clean)
nrow(ntl_clean) 
unique(ntl_clean$prefecture)

ntl_clean$prefecture = gsub("Hokkai", "Hokkaido", ntl_clean$prefecture)

write_csv(ntl_clean, '/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL_clean.csv')

summary(ntl_clean)

sum(is.na(ntl_clean$mean_ntl)) # 0
length(unique(ntl_clean$prefecture)) # 47
length(unique(ntl_clean$month)) # 147

ntl = read.csv('/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL_clean.csv')

ntl

# ============================================================
# PGA DISTRIBUTION CHECK
# ============================================================

quantile(pga_all_clean$mean_pga, probs = c(0.75, 0.9, 0.95), na.rm = T)

quantile( pga_all_clean$max_pga, probs = c(0.75, 0.9, 0.95), na.rm = T)

# ============================================================
# LOAD CLEANED DATASETS
# ============================================================

ntl = read.csv('/Users/nattanicha/Desktop/MSc_project/Japan_prefecture_NTL_clean.csv')
pga_clean = read.csv('/Users/nattanicha/Desktop/MSc_project/pga_clean_selected.csv')
econ_df = read.csv('/Users/nattanicha/Desktop/MSc_project/econ_df.csv')
eq = read.csv('/Users/nattanicha/Desktop/MSc_project/eq_df.csv')