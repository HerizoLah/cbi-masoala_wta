library(interval)
library(survival)
library(broom)
library(dplyr)
library(forcats)

# load cleaned data
df <- readRDS("data/processed/df_clean.rds")

# fit Turnbull model
turnbull_data <- df %>%
  filter(!is.na(lower), !is.na(upper))

L <- turnbull_data$lower
R <- turnbull_data$upper

stopifnot(all(L <= R))

res_turnbull <- icfit(L, R)
turnbull_summary <- summary(res_turnbull)

# parse interval labels from summary output
parse_interval <- function(x) {
  x <- gsub("\\s", "", x)
  x <- sub("^\\(|^\\[", "", x)
  x <- sub("\\]$|\\)$", "", x)
  parts <- strsplit(x, ",")[[1]]

  c(
    L = as.numeric(parts[1]),
    U = ifelse(grepl("Inf", parts[2], ignore.case = TRUE), Inf, as.numeric(parts[2]))
  )
}

cdf <- cumsum(turnbull_summary$Probability)

interval_bounds <- t(vapply(as.character(turnbull_summary$Interval), parse_interval, numeric(2)))
L_int <- interval_bounds[, "L"]
R_int <- interval_bounds[, "U"]

# median interval
k_med <- which(cdf >= 0.5)[1]
med_L <- L_int[k_med]
med_U <- R_int[k_med]
median_interval <- c(lower = med_L, upper = med_U)

# survival objects for restricted mean
S_prev <- 1 - c(0, cdf[-length(cdf)])
S_next <- 1 - cdf
R_fin <- R_int[is.finite(R_int)]
S_fin <- S_next[is.finite(R_int)]

# restricted mean WTA up to the maximum offered bid
tau <- 103.9501

if (length(R_fin) == 0) {
  mean_wta_rmst <- tau
} else {
  t_pts <- c(0, R_fin[R_fin < tau], tau)
  s_match <- S_fin[R_fin < tau]
  last_s <- if (length(s_match) == 0) 1 else tail(s_match, 1)
  S_pts <- c(1, s_match, last_s)
  mean_wta_rmst <- sum(S_pts[-length(S_pts)] * diff(t_pts))
}

# parametric interval regression
surv_obj <- Surv(time = df$lower, time2 = df$upper, type = "interval2")

model_interval1 <- survreg(
  surv_obj ~ Head_Age + Head_Education + Food_Security +
    HH_Size + Rice_Consumption + FKT + Marital_Status + Head_Gender + total_annual_income +
    Tanety_Surf + Forest_Surf + IrrRice_Surf + Pasture_surface + Pasture_Surf + Fallow_surf,
  data = df,
  dist = "lognormal"
)

label_map <- c(
  "FKTFizono" = "Village Fizono",
  "FKTMahalevona Nord" = "Village Mahalevona Nord",
  "FKTMahalevona Sud" = "Village Mahalevona Sud",
  "FKTMasindrano" = "Village Masindrano",
  "FKTTanambao" = "Village Tanambao",
  "Head_Age" = "Age of household head, per year",
  "Head_Education" = "Education of head, per level",
  "Head_GenderWoman" = "Head is woman (ref man)",
  "Marital_StatusMarried" = "Married (ref cohabiting)",
  "Marital_StatusWidowed" = "Widowed (ref cohabiting)",
  "Marital_StatusSingle" = "Single (ref cohabiting)",
  "Marital_StatusDivorced" = "Divorced (ref cohabiting)",
  "Food_Security" = "Food security, per month",
  "total_annual_income" = "Annual income, 10% increase",
  "HH_Size" = "Household size, per person",
  "Rice_Consumption" = "Daily rice consumption, per unit",
  "Tanety_Surf" = "Upland field surface",
  "Forest_Surf" = "Forest surface",
  "IrrRice_Surf" = "Irrigated rice surface",
  "Pasture_Surf" = "Pasture without cloves surface",
  "Pasture_surface" = "Pasture with cloves surface",
  "Fallow_surf" = "Fallow surface"
)

tidy_df <- tidy(model_interval1, conf.int = TRUE) %>%
  filter(term != "(Intercept)", !grepl("Log\\(scale\\)", term)) %>%
  mutate(
    pct = (exp(estimate) - 1) * 100,
    pct_lo = (exp(conf.low) - 1) * 100,
    pct_hi = (exp(conf.high) - 1) * 100,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    category = case_when(
      grepl("^FKT", term) ~ "Location (ref Ankovana)",
      grepl("Marital_Status|Head_Gender|Head_", term) ~ "Demographic",
      term %in% c("Food_Security", "total_annual_income", "HH_Size", "Rice_Consumption") ~ "Socio-economic",
      grepl("surf|surface", term, ignore.case = TRUE) ~ "Land use, percent change per hectare",
      TRUE ~ "Other"
    ),
    label = recode(term, !!!label_map, .default = term)
  ) %>%
  arrange(category, pct) %>%
  mutate(label = fct_inorder(label))

# save analysis outputs
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

saveRDS(res_turnbull, "results/res_turnbull.rds")
saveRDS(turnbull_summary, "results/turnbull_summary.rds")
saveRDS(model_interval1, "results/model_interval1.rds")
write.csv(
  data.frame(
    statistic = c("median_lower", "median_upper", "restricted_mean_tau", "restricted_mean_wta"),
    value = c(med_L, med_U, tau, mean_wta_rmst)
  ),
  "results/turnbull_key_results.csv",
  row.names = FALSE
)
write.csv(tidy_df, "results/model_interval1_tidy.csv", row.names = FALSE)
