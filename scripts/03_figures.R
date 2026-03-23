library(ggplot2)
library(dplyr)

# load cleaned data and analysis objects
df <- readRDS("data/processed/df_clean.rds")
res_turnbull <- readRDS("results/res_turnbull.rds")
turnbull_summary <- readRDS("results/turnbull_summary.rds")
tidy_df <- read.csv("results/model_interval1_tidy.csv")

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

k_med <- which(cdf >= 0.5)[1]
med_L <- L_int[k_med]
med_U <- R_int[k_med]

S_prev <- 1 - c(0, cdf[-length(cdf)])
S_next <- 1 - cdf
R_fin <- R_int[is.finite(R_int)]
S_fin <- S_next[is.finite(R_int)]

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

# Turnbull survival plot
x_max_auto <- max(R_int[is.finite(R_int)], na.rm = TRUE)
x_max <- max(10, x_max_auto)
x_breaks <- seq(0, x_max, by = 20)

rect_df <- tibble(
  xmin = L_int,
  xmax = ifelse(is.finite(R_int), R_int, x_max),
  ymin = S_next,
  ymax = S_prev
) %>%
  mutate(xmax = ifelse(xmax == xmin, xmin + 1e-6, xmax))

surv_df <- tibble(
  bid = c(0, R_fin),
  S = c(1, S_fin)
) %>%
  arrange(bid)

if (any(is.infinite(R_int))) {
  surv_df <- bind_rows(
    surv_df,
    tibble(bid = x_max, S = tail(surv_df$S, 1))
  ) %>%
    arrange(bid)
}

shift_med <- 8
shift_mean <- 12
y_med_text <- 0.92
y_mean_text <- 0.83
x_med_mid <- (med_L + med_U) / 2
x_mean <- mean_wta_rmst
x_lab <- "Monthly bid offer (USD per household)"

p_turnbull <- ggplot() +
  geom_rect(
    data = rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey60",
    alpha = 0.35,
    color = NA
  ) +
  geom_step(
    data = surv_df,
    aes(x = bid, y = S),
    linewidth = 0.8
  ) +
  geom_vline(xintercept = med_L, linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = med_U, linetype = "dashed", linewidth = 0.6) +
  geom_segment(aes(x = med_L, xend = med_U, y = 0.5, yend = 0.5), linewidth = 0.6) +
  geom_vline(xintercept = mean_wta_rmst, linetype = "dotdash", linewidth = 0.7) +
  annotate(
    "segment",
    x = x_med_mid + shift_med,
    xend = x_med_mid,
    y = y_med_text,
    yend = 0.52,
    arrow = arrow(length = grid::unit(0.18, "cm"))
  ) +
  annotate(
    "text",
    x = x_med_mid + shift_med,
    y = y_med_text,
    label = paste0(
      "WTA median interval: [",
      round(med_L, 2), ", ", round(med_U, 2), "]"
    ),
    hjust = 0
  ) +
  annotate(
    "segment",
    x = x_mean + shift_mean,
    xend = x_mean,
    y = y_mean_text,
    yend = 0.58,
    arrow = arrow(length = grid::unit(0.18, "cm"))
  ) +
  annotate(
    "text",
    x = x_mean + shift_mean,
    y = y_mean_text,
    label = paste0(
      "Mean WTA truncated to max bid offered (~",
      round(tau, 0), "): ", round(mean_wta_rmst, 1)
    ),
    hjust = 0
  ) +
  scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "",
    x = x_lab,
    y = "Proportion rejecting offer (S(B))"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "oblique"),
    panel.grid.minor = element_blank()
  )

# coefficient plot
tidy_df$color <- ifelse(tidy_df$pct < 0, "red", "blue")

p_coeff <- ggplot(tidy_df, aes(x = pct, y = label, color = color)) +
  geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.25) +
  geom_point() +
  geom_text(
    aes(label = paste0(round(pct, 1), "%", sig)),
    vjust = -0.8,
    size = 4,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_identity() +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 14)
  )

# save figures
if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)

ggsave("results/figures/turnbull_survival_plot.png", p_turnbull, width = 10, height = 7, dpi = 300)
ggsave("results/figures/interval_regression_coefficients.png", p_coeff, width = 10, height = 9, dpi = 300)
