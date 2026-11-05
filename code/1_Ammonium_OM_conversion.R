#### Step 1: NH4 conversion, OM per vial, and QC flags ####
# Adopted from Meredith Brehob 2021; Krause and Loria 2023 edits
# NH4:      NH4 = (1 - f_NH3) x measured ammonia (YSI temp and pH at outlet)
# sediment: OM_g = wet g in vial x f_dry x f_OM
# biofilm:  OM_g = slurry OM (g L-1) x vial volume (L)
# QC:       flags on lab notes and measured concentrations only

library(tidyverse)
library(lubridate)
library(readxl)

setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/")

# site names -> incubation names (NS1-3 pooled to 3 m site)
fix_site <- function(x) x %>%
  str_replace("^(BW|GB|SS|SH)NS\\d?$", "\\1_3m") %>%
  str_replace("^(BW|GB)10m$", "\\1_10m")


#### 1. Read in data ####
inc  <- read_csv("Step_1/Step1_N_incubation_raw_data.csv", show_col_types = FALSE)  # all months
sed  <- read_csv("Step_1/Step1_sediment_OM.csv", show_col_types = FALSE)
bio  <- read_csv("Step_1/Step1_biofilm_OM.csv", show_col_types = FALSE)
lake <- read_excel("Step_1/Step1_Lake_data_ammonium_conversion.xlsx")  # Shore, Inc_month, Temp_degC, pH

inc <- inc %>%
  filter(!is.na(Type)) %>% # drop empty rows
  mutate(Sample_Date = as.Date(Sample_Date),
         Location    = str_sub(Site, 1, 2),
         Shore       = case_when(Location %in% c("GB", "SH") ~ "E",
                                 Location %in% c("BW", "SS") ~ "W"))

#### 2. NH4 conversion ####
inc <- inc %>%
  left_join(lake %>% select(Shore, Inc_month, Temp_degC, pH) %>%
              mutate(across(c(Temp_degC, pH), as.numeric)),
            by = c("Shore", "Inc_month")) %>%
  mutate(pKa          = 0.09018 + 2727.92 / (Temp_degC + 273.15),
         f_NH3        = 1 / (10^(pKa - pH) + 1),
         NH4_mgNL     = if_else(Analyte == "NH3", (1 - f_NH3) * Conc_mgNL, NA_real_),
         NO2_NO3_mgNL = if_else(Analyte == "NO3", Conc_mgNL, NA_real_)) %>%
  select(-pKa, -f_NH3)

# NH4 should fall just below the 1:1 line
inc %>% filter(Analyte == "NH3") %>%
  ggplot(aes(Conc_mgNL, NH4_mgNL, color = Type)) +
  geom_abline(linetype = "dashed") + geom_point() +
  facet_wrap(~ Inc_month) + theme_bw()


#### 3. Sediment OM ####
sed_core <- sed %>%
  mutate(Date = mdy(Date), Site = fix_site(Site),
         Inc_month = format(Date, "%B"),
         f_dry = Dry_mass_g / Wet_mass_g,  # g dry / g wet
         f_OM  = (Dry_weight_plus_tin_g - Combusted_weight_plus_tin_g) /
           (Dry_weight_plus_tin_g - Mean_tin_weight_g), # g OM / g dry
         f_OM  = if_else(f_OM < 0, NA_real_, f_OM),   # GB0.5m May: combusted > dry
         f_OM  = if_else(Sample_ID == "NS23_F48", NA_real_, f_OM)) %>% # SHNS2 July outlier
  filter(Date != as.Date("2023-06-16"))  # extra BW sediment, not incubated

# mean per site x month (NS1-3 pooled); GB0.5m May f_OM = mean of June, July
sed_mean <- sed_core %>%
  group_by(Site, Inc_month) %>%
  summarise(f_dry = mean(f_dry), f_OM = mean(f_OM, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(f_OM = if_else(is.nan(f_OM), mean(f_OM[!is.nan(f_OM)]), f_OM)) %>%
  ungroup()


#### 4. Biofilm OM ####
bio_mean <- bio %>%
  mutate(Date = mdy(Date), Inc_month = format(Date, "%B")) %>%
  filter(Site %in% c("BW0.5m", "GB0.5m")) %>% # SSNS not incubated
  select(Site, Inc_month, OM_g_L)


#### 5. OM per vial ####
inc <- inc %>%
  left_join(sed_mean, by = c("Site", "Inc_month")) %>%
  left_join(bio_mean, by = c("Site", "Inc_month")) %>%
  mutate(OM_g = case_when(Type == "sediment" ~ Weight_g * f_dry * f_OM, # grams of sediment
                          Type == "biofilm"  ~ OM_g_L * Vial_vol_L)) # grams of sediment

# n_NA should be 0
inc %>% filter(Type %in% c("sediment", "biofilm")) %>%
  group_by(Type, Site, Inc_month) %>%
  summarise(n = n(), n_NA = sum(is.na(OM_g)), OM_g = mean(OM_g, na.rm = TRUE), .groups = "drop") %>%
  print(n = Inf)


#### 6. QC ####
# a) lab notes physically compromised vials -> removed
bad_vials <- tribble(
  ~Sample_ID,   ~Inc_month, ~reason,
  "Sed_NA_5",   "May",  "cracked, not in experiment",
  "Sed_NA_146", "May",  "cracked, spilled",
  "Sed_NA_243", "June", "half volume",
  "Sed_NA_213", "June", "filter broke"
)

inc <- inc %>% anti_join(bad_vials, by = c("Sample_ID", "Inc_month"))

# b1) concentration (µg N L-1); below detection -> half the detection limit
inc <- inc %>%
  mutate(Conc_mgNL = case_when(Analyte == "NH3" ~ NH4_mgNL,
                               Analyte == "NO3" ~ NO2_NO3_mgNL),
         Conc_mgNL = case_when(Analyte == "NH3" & Conc_mgNL < 0.002 ~ 0.001,
                               Analyte == "NO3" & Conc_mgNL < 0.003 ~ 0.0015,
                               TRUE ~ Conc_mgNL),
         Conc_ugNL = 1000 * Conc_mgNL,
         grp       = if_else(Type == "spike", "spike", "sample"))

# b2) triplicate median; only the vial farthest from it can be an outlier (n = 3)
inc <- inc %>%
  group_by(Site, Inc_month, Type, Analyte, Spike_µg_L) %>%
  mutate(n_rep   = sum(!is.na(Conc_ugNL)),
         rep_med = median(Conc_ugNL, na.rm = TRUE),
         is_max  = abs(Conc_ugNL - rep_med) == max(abs(Conc_ugNL - rep_med), na.rm = TRUE)) %>%
  ungroup()

# b3) deviation from triplicate median
inc <- inc %>%
  mutate(r        = (Conc_ugNL - rep_med) / rep_med,          # relative deviation
         abs_dev  = abs(Conc_ugNL - rep_med),                  # µg N L-1
         dl_floor = if_else(Analyte == "NH3", 5 * 2, 5 * 3))   # 5 x detection limit

# c) flag outliers (MAD x 1.4826 = robust SD)
# samples: absolute, by analyte x spike level
# spikes > 0: relative, pooled; 0 spike: absolute, by analyte

inc <- inc %>%
  mutate(rel = grp == "spike" & Spike_µg_L > 0,
         key = case_when(rel            ~ "spike_rel",
                         grp == "spike" ~ paste("spike0", Analyte),
                         TRUE           ~ paste("sample", Analyte, Spike_µg_L)))

mad_tab <- inc %>%
  filter(n_rep >= 3, abs_dev > 0) %>%
  group_by(key) %>%
  summarise(mad = 1.4826 * median(if_else(rel, abs(r), abs_dev)), n = n(), .groups = "drop") %>%
  mutate(cutoff = 3.5 * mad)  # µg N L-1, or fraction of median for spike_rel
print(mad_tab, n = Inf)

# d) flag outlier vials, > 3.5 MAD and > 5 x DL
inc <- inc %>%
  left_join(mad_tab %>% select(key, mad), by = "key") %>%
  mutate(cut = 3.5 * mad * if_else(rel, rep_med, 1),
         out = is_max & n_rep >= 3 & abs_dev > cut & abs_dev > dl_floor,
         qc_flag = case_when(out & grp == "spike"  ~ "spike_outlier",
                             out & grp == "sample" ~ "sample_outlier")) %>%
  select(-grp, -n_rep, -rep_med, -r, -abs_dev, -dl_floor, -is_max, -rel, -key, -mad, -cut, -out)

# review
flags <- inc %>% filter(!is.na(qc_flag)) %>%
  select(Site, Sample_ID, Inc_month, Type, Analyte, Spike_µg_L, Conc_ugNL, qc_flag, Notes) %>%
  arrange(qc_flag, Site, Inc_month)
print(flags, n = Inf)
count(inc, qc_flag)

# flagged vials by site
plot_flags <- function(type) inc %>%
  filter(Type == type, !is.na(Analyte)) %>%
  ggplot(aes(Spike_µg_L, Conc_ugNL, color = replace_na(qc_flag, "ok"))) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_grid(Site ~ Inc_month + Analyte, scales = "free_y") +
  scale_color_manual(values = flag_colors) +
  labs(x = "Spike (µg N L⁻¹)", y = "Measured (µg N L⁻¹)", color = NULL, title = type)

outliers_sed <- plot_flags("sediment")
outliers_bio <- plot_flags("biofilm")

#ggsave(plot = outliers_sed, "Step_1/outliers_sed.png",width=12,height=8)
#ggsave(plot = outliers_bio, "Step_1/outliers_bio.png",width=12,height=8)


#### 7. Write data ####
inc_export <- inc %>%
  select(Inc_month, Sample_Date, Site, Location, Shore, Type, Composite_num,
         Analyte, Spike_µg_L, Volume_lake_water_L, Weight_g, Vial_vol_L, 
         Incubation_start_time, Incubation_end_time,
         Actual_incubation_hrs, Temp_degC, pH, NH4_mgNL, NO2_NO3_mgNL,
         f_dry, f_OM, OM_g, qc_flag, Notes) %>%
  mutate(Inc_month = factor(Inc_month, levels = c("May", "June", "July"))) %>%
  arrange(Inc_month, Site, Type, Analyte, Spike_µg_L, Composite_num)

#write_csv(inc,"Step_2/Step2_NH4_OM_Ninc_data.csv")
