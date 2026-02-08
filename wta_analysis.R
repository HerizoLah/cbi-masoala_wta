library(readxl)
library(interval)
library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(forcats)

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
# Use a relative path or an environment variable so the script is portable.
# Example: DATA_PATH="data/df_17-12-2025.xlsx" Rscript wta_analysis.R

data_path <- Sys.getenv("DATA_PATH", "data/df_17-12-2025.xlsx")
if (!file.exists(data_path)) {
  stop(
    "Data file not found. Set DATA_PATH or place the file at ",
    data_path,
    call. = FALSE
  )
}

df <- read_excel(data_path, col_names = TRUE, na = c("", "NA"))

# Convert all character variables into factors.
df <- df %>%
  mutate(across(where(is.character), as.factor))

required_cols <- c(
  "lower", "upper", "Head_Age", "Head_Education", "Food_Security",
  "HH_Size", "Rice_Consumption", "FKT", "Marital_Status", "Head_Gender",
  "total_annual_income", "Tanety_Surf", "Forest_Surf", "IrrRice_Surf",
  "Pasture_surface", "Pasture_Surf", "Fallow_surf"
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1) Turnbull NPLME
# -----------------------------------------------------------------------------
L <- as.numeric(df$lower)
R <- as.numeric(df$upper)

# Keep complete cases only.
keep <- !is.na(L) & !is.na(R)
L <- L[keep]
R <- R[keep]

# Basic validity check.
stopifnot(all(L <= R))

res_turnbull <- icfit(L, R)
summary(res_turnbull)

# Plot survival curve.
plot(res_turnbull)

# -----------------------------------------------------------------------------
# 2) Parametric model (interval regression)
# -----------------------------------------------------------------------------
surv_obj <- Surv(time = df$lower, time2 = df$upper, type = "interval2")

model_interval1 <- survreg(
  surv_obj ~ Head_Age + Head_Education + Food_Security +
    HH_Size + Rice_Consumption + FKT + Marital_Status + Head_Gender +
    total_annual_income + Tanety_Surf + Forest_Surf + IrrRice_Surf +
    Pasture_surface + Pasture_Surf + Fallow_surf,
  data = df,
  dist = "lognormal"
)

# Tidy and transform to percent change in predicted median WTA.
# Scale income to reflect a 10% increase.
tidy_df <- tidy(model_interval1, conf.int = TRUE) %>%
  filter(term != "(Intercept)", !grepl("Log\\(scale\\)", term)) %>%
  mutate(
    effect = if_else(term == "total_annual_income", estimate * log(1.10), estimate),
    effect_lo = if_else(term == "total_annual_income", conf.low * log(1.10), conf.low),
    effect_hi = if_else(term == "total_annual_income", conf.high * log(1.10), conf.high),
    pct = (exp(effect) - 1) * 100,
    pct_lo = (exp(effect_lo) - 1) * 100,
    pct_hi = (exp(effect_hi) - 1) * 100,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Categories.
tidy_df <- tidy_df %>%
  mutate(category = case_when(
    grepl("^FKT", term) ~ "Location (ref Ankovana)",
    grepl("Marital_Status|Head_Gender|Head_", term) ~ "Demographic",
    term %in% c("Food_Security", "total_annual_income", "HH_Size", "Rice_Consumption") ~
      "Socio-economic",
    grepl("surf|surface", term, ignore.case = TRUE) ~
      "Land use, percent change per hectare",
    TRUE ~ "Other"
  ))

# Relabel with units.
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

tidy_df <- tidy_df %>%
  mutate(label = recode(term, !!!label_map, .default = term)) %>%
  arrange(category, pct) %>%
  mutate(label = fct_inorder(label))

# Color by direction on the percent scale.
tidy_df <- tidy_df %>%
  mutate(color = ifelse(pct < 0, "red", "blue"))

# Plot.
ggplot(tidy_df, aes(x = pct, y = label, color = color)) +
  geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.25) +
  geom_point() +
  geom_text(
    aes(label = paste0(round(pct, 1), "%", sig)),
    vjust = -0.8,
    size = 3,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_identity() +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Percent change in predicted median WTA, relative to reference",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 10)
  )
