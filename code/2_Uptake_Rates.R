#### Step 2: Uptake rate calculation ####
# Script created by: Heili Lowman, KAL, JAK
#
# 1.   Read and check data
# 2.   Incubation time
# 3.   Spike standards: check outliers, then drop; C_pre = mean spike standard
#      by shore, month, analyte, spike level (0 spike = ambient lake water)
# 4.   Vial rate = (C_post - C_pre) / incubation h
# 5.   Sample outliers: look, check influence, then drop
# 6.   Pelagic rate = mean lake-water vials (SH uses GB water, SS uses BW water);
#      check spread and influence on net uptake
# 7.   Net = vial rate - pelagic rate; OM uptake = -net x V / g OM (ug N g OM-1 h-1)

library(tidyverse)
library(patchwork)

setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/")

month_levels <- c("May", "June", "July")
loc_colors   <- c(BW = "#3283a8", GB = "#a67d17", SS = "#86bedb", SH = "#dbbb6e")
flag_colors  <- c(ok = "grey60", spike_outlier = "red", sample_outlier = "blue")
theme_set(theme_bw(base_size = 13))

# "6:08" or hms 06:08:00 -> 6.13 h
to_hr <- function(x) {
  p <- str_split_fixed(as.character(x), ":", 3)
  as.numeric(p[, 1]) + as.numeric(p[, 2]) / 60
}


#### 1. Read data ####
dat <- read_csv("Step_1/Step1_NH4_OM_Ninc_data.csv", show_col_types = FALSE) %>%
  filter(!is.na(Analyte)) %>%
  mutate(Inc_month = factor(Inc_month, levels = month_levels),
         flag      = replace_na(qc_flag, "ok"))

count(dat, Inc_month, Type) %>% pivot_wider(names_from = Type, values_from = n)   # spikes each month?
count(dat, flag)

# check vials all have a concentration
dat %>% filter(is.na(Conc_ugNL)) %>% select(Sample_ID, Inc_month, Type, Analyte, Spike_µg_L)


#### 2. Incubation time ####
# some missing values; use end - start
dat <- dat %>%
  mutate(end_fix = is.na(Incubation_end_time) & !is.na(Actual_incubation_hrs) &
           to_hr(Actual_incubation_hrs) > 12,
         inc_hr  = if_else(end_fix,
                           to_hr(Actual_incubation_hrs) - to_hr(Incubation_start_time),
                           to_hr(Actual_incubation_hrs)))

# check: corrected vials (~6 h)
dat %>% filter(end_fix) %>%
  select(Sample_ID, Inc_month, Incubation_start_time, Actual_incubation_hrs, inc_hr)

# check: distribution
summary(dat$inc_hr[dat$Type != "spike"])
ggplot(dat %>% filter(Type != "spike"), aes(inc_hr, fill = Inc_month)) +
  geom_histogram(bins = 30) +
  labs(x = "Incubation time (h)", y = "Vials", title = "Incubation times")

# fill missing with mean of incubated vials
mean_inc_hr <- mean(dat$inc_hr[dat$Type != "spike"], na.rm = TRUE)
mean_inc_hr
dat <- dat %>% mutate(inc_hr_c = if_else(is.na(inc_hr), mean_inc_hr, inc_hr))


#### 3. Spike standards ####
spikes <- dat %>% filter(Type == "spike")

# check: measured vs nominal, outliers in red
ggplot(spikes, aes(Spike_µg_L, Conc_ugNL, color = flag)) +
  geom_abline(linetype = "dashed", color = "grey40") +
  geom_point(size = 2, alpha = 0.8) +
  facet_grid(Analyte ~ Inc_month + Shore) +
  scale_color_manual(values = flag_colors) +
  labs(x = "Nominal spike (µg N L-1)", y = "Measured (µg N L-1)")

# spike outliers
spikes %>% filter(flag == "spike_outlier") %>%
  select(Shore, Sample_ID, Inc_month, Analyte, Spike_µg_L, Conc_ugNL, any_of(c("Flag", "Notes")))

spike_means <- function(d) d %>%
  filter(Type == "spike") %>%
  group_by(Shore, Inc_month, Analyte, Spike_µg_L) %>%
  summarise(Mean_Spike_Conc = mean(Conc_ugNL), SD_Spike = sd(Conc_ugNL), n_spike = n(),
            .groups = "drop")

spike_cmp <- spike_means(spikes) %>%
  left_join(spike_means(filter(spikes, flag != "spike_outlier")),
            by = c("Shore", "Inc_month", "Analyte", "Spike_µg_L"), suffix = c("_all", "_clean")) %>%
  mutate(diff_ug  = Mean_Spike_Conc_clean - Mean_Spike_Conc_all,
         diff_pct = 100 * diff_ug / Mean_Spike_Conc_all)

# check: spike means that change when outliers are dropped
spike_cmp %>% filter(diff_ug != 0) %>%
  select(Shore, Inc_month, Analyte, Spike_µg_L, n_spike_all, n_spike_clean,
         Mean_Spike_Conc_all, Mean_Spike_Conc_clean, diff_ug, diff_pct)

# check: 0-spike background (ambient lake water)
ggplot(spike_cmp %>% filter(Spike_µg_L == 0),
       aes(Inc_month, Mean_Spike_Conc_clean, color = Shore, group = Shore)) +
  geom_line() + geom_point(size = 2.5) +
  facet_wrap(~ Analyte, scales = "free_y") +
  labs(x = NULL, y = "Background (µg N L-1)", title = "0-spike standards")

# DROP spike outliers; C_pre = mean spike standard (0 spike = ambient lake water)
dat   <- dat %>% filter(flag != "spike_outlier")
c_pre <- spike_means(dat)

# check: every incubated vial has a C_pre
dat %>% filter(Type != "spike") %>% distinct(Shore, Inc_month, Analyte, Spike_µg_L) %>%
  anti_join(c_pre, by = c("Shore", "Inc_month", "Analyte", "Spike_µg_L"))


#### 4. Vial rates of change ####
vials <- dat %>%
  filter(Type %in% c("sediment", "biofilm", "water")) %>%
  left_join(c_pre %>% select(Shore, Inc_month, Analyte, Spike_µg_L, Mean_Spike_Conc),
            by = c("Shore", "Inc_month", "Analyte", "Spike_µg_L")) %>%
  mutate(delta_Conc_ugNL   = Conc_ugNL - Mean_Spike_Conc,
         delta_Conc_ugNLhr = delta_Conc_ugNL / inc_hr_c,
         Assigned_Location = dplyr::recode(Location, SH = "GB", SS = "BW"))

# net uptake from a set of vials (water and substrate)
calc_net <- function(v) {
  water <- v %>%
    filter(Type == "water") %>%
    group_by(Inc_month, Location, Analyte, Spike_µg_L) %>%
    summarise(water_delta_Conc_ugNLhr = mean(delta_Conc_ugNLhr),
              water_sd = sd(delta_Conc_ugNLhr), n_water = n(), .groups = "drop")
  v %>%
    filter(Type %in% c("sediment", "biofilm")) %>%
    left_join(water, by = c("Assigned_Location" = "Location", "Inc_month", "Analyte", "Spike_µg_L")) %>%
    mutate(net_delta_Conc_ugNLhr    = delta_Conc_ugNLhr - water_delta_Conc_ugNLhr,
           net_delta_Conc_ugNLhr_OM = -net_delta_Conc_ugNLhr * Volume_lake_water_L / OM_g)
}

# means per site x month x substrate x analyte x spike level
group_means <- function(u) u %>%
  group_by(Site, Location, Inc_month, Type, Analyte, Spike_µg_L) %>%
  summarise(raw = mean(-net_delta_Conc_ugNLhr, na.rm = TRUE),
            OM  = mean(net_delta_Conc_ugNLhr_OM, na.rm = TRUE), .groups = "drop")

# visualize with and without organic matter correction
group_means(upt) %>%
  pivot_longer(c(raw, OM), names_to = "units") %>%
  ggplot(aes(Site, value, fill = Location)) +
  geom_boxplot() +
  facet_grid(units ~ Type + Analyte, scales = "free") +
  scale_fill_manual(values = loc_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### 5. Sample outliers: look, check influence, then drop ####
# % of incubated vials flagged (sediment, biofilm, water)
vials %>%
  summarise(n_vials   = n(),
            n_flagged = sum(flag == "sample_outlier"),
            pct       = round(100 * n_flagged / n_vials, 1))

# by substrate
vials %>%
  group_by(Type) %>%
  summarise(n_vials   = n(),
            n_flagged = sum(flag == "sample_outlier"),
            pct       = round(100 * n_flagged / n_vials, 1), .groups = "drop")

# where are the sample outliers?
vials %>% filter(flag == "sample_outlier") %>%
  count(Type, Site, Inc_month, Analyte) %>% arrange(desc(n)) %>% print(n = Inf)

# all vials in groups with at least one sample outlier
out_groups <- vials %>%
  semi_join(vials %>% filter(flag == "sample_outlier") %>%
              distinct(Type, Site, Inc_month, Analyte),
            by = c("Type", "Site", "Inc_month", "Analyte")) %>%
  mutate(panel = paste(Type, Site, Inc_month, Analyte))

ggplot(out_groups, aes(Spike_µg_L, Conc_ugNL, color = flag)) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_color_manual(values = flag_colors) +
  labs(x = "Spike (µg N L⁻¹)", y = "Measured (µg N L⁻¹)", color = NULL)

# check: GB0.5m May biofilm (many outliers, very low OM)
vials %>% filter(Site == "GB0.5m", Inc_month == "May", Type == "biofilm") %>%
  ggplot(aes(Spike_µg_L, Conc_ugNL, color = flag)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 2.5) +
  facet_wrap(~ Analyte) +
  scale_color_manual(values = flag_colors) +
  labs(x = "Spike (µg N L⁻¹)", y = "Measured (µg N L⁻¹)", title = "GB0.5m May biofilm")

# influence: uptake with vs without sample outliers
upt_keep <- calc_net(vials)
upt_drop <- calc_net(filter(vials, flag != "sample_outlier"))

samp_cmp <- group_means(upt_keep) %>%
  left_join(group_means(upt_drop),
            by = c("Site", "Location", "Inc_month", "Type", "Analyte", "Spike_µg_L"),
            suffix = c("_keep", "_drop")) %>%
  mutate(diff_OM  = OM_drop - OM_keep,
         diff_pct = 100 * diff_OM / abs(OM_keep),
         big      = abs(diff_pct) > 20 & abs(diff_OM) > 5)   # > 20% and > 5 ug N g OM-1 h-1

# check: groups that change a lot
samp_cmp %>% filter(big) %>%
  select(Site, Inc_month, Type, Analyte, Spike_µg_L, OM_keep, OM_drop, diff_OM, diff_pct) %>%
  arrange(desc(abs(diff_OM))) %>% print(n = Inf)
samp_cmp %>% count(Type, Analyte, big)

# check: do key results hold either way?
bind_rows(group_means(upt_keep) %>% mutate(version = "keep"),
          group_means(upt_drop) %>% mutate(version = "drop")) %>%
  group_by(version, Type, Analyte) %>%
  summarise(median_OM = median(OM, na.rm = TRUE), median_raw = median(raw, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(Type, Analyte, version)

##### drop sample outliers for uptake calculations
vials_use <- vials %>% filter(flag != "sample_outlier")



#### 6. Pelagic (lake water) uptake ####
water_vials <- vials_use %>% filter(Type == "water")

# visualize spread
pel_plot <- ggplot(water_vials, aes(factor(Spike_µg_L), delta_Conc_ugNLhr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, color = "grey40") +
  facet_grid(Analyte ~ Inc_month + Location, scales = "free_y") +
  labs(x = "Spike (µg N L-1)", y = expression(Pelagic~rate~(µg~N~L^-1~h^-1)))

ggsave(plot=pel_plot, "Step_1/figures/pel_plot.png",width=16, height=8)

# influence: pelagic correction vs substrate rate at the same level
upt <- calc_net(vials_use)

pel_infl <- upt %>%
  group_by(Inc_month, Assigned_Location, Site, Type, Analyte, Spike_µg_L) %>%
  summarise(sub_rate = mean(delta_Conc_ugNLhr),
            pel_rate = first(water_delta_Conc_ugNLhr),
            pel_sd   = first(water_sd),
            n_water  = first(n_water),
            net_rate = mean(net_delta_Conc_ugNLhr), .groups = "drop") %>%
  mutate(pel_share = abs(pel_rate) / (abs(sub_rate) + abs(pel_rate)),   # 0-1, share of correction
         sd_ratio  = pel_sd / abs(net_rate))                            # pelagic noise vs net signal

# check: groups where pelagic correction is > half the signal, or noise > net
pel_infl %>% filter(pel_share > 0.5 | sd_ratio > 1) %>%
  arrange(desc(pel_share)) %>% print(n = Inf)
pel_infl %>% summarise(median_share   = median(pel_share, na.rm = TRUE),
                       prop_over_half = mean(pel_share > 0.5, na.rm = TRUE))

# check: gross (no correction) vs net uptake per site, month, and spike level
gross_net <- upt %>%
  group_by(Site, Location, Inc_month, Type, Analyte, Spike_µg_L) %>%
  summarise(gross = mean(-delta_Conc_ugNLhr), net = mean(-net_delta_Conc_ugNLhr), .groups = "drop")

pelagic_net_plt <- ggplot(gross_net, aes(gross, net, color = Location, shape = Type)) +
  geom_abline(linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_point(size = 2.5) +
  facet_grid(Analyte ~ Inc_month) +
  scale_color_manual(values = loc_colors) +
  labs(x = expression(Gross~uptake~(µg~N~L^-1~h^-1)), y = "Net uptake (pelagic subtracted)",
       title = "Influence of pelagic correction")

#ggsave(plot=pelagic_net_plt, "Step_1/figures/pelagic_net_plt.png",width=10, height=8)

# does the pelagic correction change the key results?
gross_net %>%
  group_by(Type, Analyte) %>%
  summarise(median_gross = median(gross), median_net = median(net),
            rho = cor(gross, net, method = "spearman"), .groups = "drop")


#### 7. Final rates and export ####
dat_tidy <- upt %>%
  select(Inc_month, Site, Location, Analyte, Spike_µg_L, Mean_Spike_Conc,
         Volume_lake_water_L, Type, NO2_NO3_mgNL, NH4_mgNL,
         Conc_ugNLK = Conc_ugNL, pH, Weight_g, OM_g,
         real_icu_hrc = inc_hr_c,
         water_delta_Conc_ugNLhr, net_delta_Conc_ugNLhr, net_delta_Conc_ugNLhr_OM)

# check: missing rates (no C_pre, no water, or no OM)
dat_tidy %>% filter(is.na(net_delta_Conc_ugNLhr_OM)) %>%
  count(Inc_month, Site, Type, Analyte, Spike_µg_L)

# check: final rates
dat_tidy %>%
  pivot_longer(c(net_delta_Conc_ugNLhr, net_delta_Conc_ugNLhr_OM), names_to = "metric") %>%
  ggplot(aes(Location, value, fill = Location)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(metric ~ Inc_month + Analyte, scales = "free_y") +
  scale_fill_manual(values = loc_colors) +
  labs(x = NULL, y = NULL, title = "Final rates") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

#write_csv(dat_tidy,  "Step_2/Step2_Final_Uptake_Rates.csv")