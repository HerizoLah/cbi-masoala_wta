library(readxl)
library(dplyr)

# read data
df <- read_excel("data/df_17-12-2025.xlsx", col_names = TRUE, na = c("", "NA"))

# convert bid variables to numeric
df <- df %>%
  mutate(
    across(
      c(bid1_amount, bid2_amount, bid1_consent, bid2_consent),
      ~ as.numeric(as.character(.x))
    )
  )

# build interval bounds for Turnbull estimation
get_bounds <- function(b1, b2, y1, y2) {
  if (any(is.na(c(b1, b2, y1, y2)))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  lo <- min(b1, b2)
  hi <- max(b1, b2)

  if (y1 == 1 && y2 == 1) {
    return(c(lower = 0, upper = lo))
  }

  if ((y1 == 1 && y2 == 0) || (y1 == 0 && y2 == 1)) {
    return(c(lower = lo, upper = hi))
  }

  if (y1 == 0 && y2 == 0) {
    return(c(lower = hi, upper = Inf))
  }

  c(lower = NA_real_, upper = NA_real_)
}

bounds <- t(
  mapply(
    get_bounds,
    df$bid1_amount,
    df$bid2_amount,
    df$bid1_consent,
    df$bid2_consent
  )
)

df_clean <- df %>%
  mutate(
    lower = bounds[, "lower"],
    upper = bounds[, "upper"]
  ) %>%
  mutate(across(where(is.character), as.factor))

# save cleaned data for the next scripts
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)

saveRDS(df_clean, "data/processed/df_clean.rds")
write.csv(df_clean, "data/processed/df_clean.csv", row.names = FALSE)
