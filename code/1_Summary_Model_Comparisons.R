## Summarize model fit comparisons
# JRB

## Import packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse"), require, character.only=T)

###########################
## Extract best model fit
###########################

## Import data
loo_df <- read.csv("../data/loo_comparisons_df_2026_10_29.csv", header = T)

# Extract the best-fitting model (elpd_diff == 0) for each .id
best_models <- loo_df |>
  group_by(.id) |>
  filter(elpd_diff == 0) |>
  ungroup()

# Sanity check: each .id should have exactly one best model
stopifnot(n_distinct(best_models$.id) == n_distinct(loo_df$.id),
          nrow(best_models) == n_distinct(loo_df$.id))
nrow(print(best_models, n = Inf))

# Quick tally of which model wins most often
best_models |> count(model, sort = TRUE)

#################################################################################
## Link best model to model diagnostics to confirm acceptable model
#################################################################################












