# question_3_tlg/02_create_visualizations.R
# Q3 - TLG Adverse Events Visualizations
# Roche ADS Programmer Assessment

# Purpose:
#   Create two AE visualizations using pharmaverseadam::adae and pharmaverseadam::adsl

# Outputs:
#   question_3_tlg/output/ae_severity_by_treatment.png
#   question_3_tlg/output/top10_ae_incidence_ci.png
#   question_3_tlg/logs/q3_plots_log.txt

# -----------------------------
# 0) Setup
# -----------------------------

options(stringsAsFactors = FALSE)

dir.create("question_3_tlg/output", recursive = TRUE, showWarnings = FALSE)
dir.create("question_3_tlg/logs", recursive = TRUE, showWarnings = FALSE)

log_file <- "question_3_tlg/logs/q3_plots_log.txt"
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Q3 AE Visualizations Log\n")
cat("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------
# 1) Load packages
# -----------------------------
required_pkgs <- c("dplyr", "tidyr", "stringr", "tibble", "ggplot2", "pharmaverseadam")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them before running this script.")
}

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(pharmaverseadam)

cat("Loaded packages successfully.\n\n")

# -----------------------------
# 2) Load datasets
# -----------------------------
data("adae", package = "pharmaverseadam")
data("adsl", package = "pharmaverseadam")

cat("Input dimensions:\n")
cat("ADAE:", paste(dim(adae), collapse = " x "), "\n")
cat("ADSL:", paste(dim(adsl), collapse = " x "), "\n\n")

# -----------------------------
# 3) Helper: Wilson CI
# -----------------------------
wilson_ci <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  denom <- 1 + (z^2 / n)
  center <- (p + (z^2 / (2 * n))) / denom
  half <- (z / denom) * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))
  c(max(0, center - half), min(1, center + half))
}

# -----------------------------
# 4) TEAE data prep
# -----------------------------
adae_teae <- adae %>%
  mutate(
    TRTEMFL = as.character(TRTEMFL),
    ACTARM  = as.character(ACTARM),
    AETERM  = as.character(AETERM),
    AESEV   = as.character(AESEV)
  ) %>%
  filter(TRTEMFL == "Y") %>%
  mutate(
    ACTARM = ifelse(is.na(ACTARM) | str_trim(ACTARM) == "", "Missing", ACTARM),
    AETERM = ifelse(is.na(AETERM) | str_trim(AETERM) == "", "Missing", AETERM),
    AESEV  = ifelse(is.na(AESEV)  | str_trim(AESEV)  == "", "Missing", AESEV)
  )

cat("TEAE records:", nrow(adae_teae), "\n\n")

# -----------------------------
# 5) Plot 1 - AE severity distribution by treatment (EVENT COUNTS)
# -----------------------------
# PDF sample matches event counts, not unique-subject percentages
sev_counts <- adae_teae %>%
  count(ACTARM, AESEV, name = "n_ae")

sev_order <- c("MILD", "MODERATE", "SEVERE", "Missing")
sev_counts <- sev_counts %>%
  mutate(
    AESEV = factor(AESEV, levels = unique(c(sev_order, sort(setdiff(unique(AESEV), sev_order)))))
  )

cat("Plot 1 data preview (event counts):\n")
print(sev_counts)
cat("\n")

p1 <- ggplot(sev_counts, aes(x = ACTARM, y = n_ae, fill = AESEV)) +
  geom_col(position = "stack") +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count of AEs",
    fill = "Severity/Intensity"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = "question_3_tlg/output/ae_severity_by_treatment.png",
  plot = p1,
  width = 8, height = 6, dpi = 300
)

cat("Saved Plot 1: question_3_tlg/output/ae_severity_by_treatment.png\n\n")

# -----------------------------
# 6) Plot 2 - Top 10 most frequent AEs (overall incidence with 95% CI)
# -----------------------------
# Use unique-subject incidence for top 10 AEs
ae_subj_term_overall <- adae_teae %>%
  distinct(USUBJID, AETERM)

n_teae_subj <- adae_teae %>%
  distinct(USUBJID) %>%
  nrow()

cat("Overall TEAE subject denominator (n):", n_teae_subj, "\n")

top10_overall <- ae_subj_term_overall %>%
  count(AETERM, name = "n_subj") %>%
  arrange(desc(n_subj), AETERM) %>%
  slice_head(n = 10)

cat("Top 10 AETERM overall:\n")
print(top10_overall)
cat("\n")

top10_plot_df <- top10_overall %>%
  mutate(
    rate = n_subj / n_teae_subj,
    ci_low = vapply(n_subj, function(x) wilson_ci(x, n_teae_subj)[1], numeric(1)),
    ci_high = vapply(n_subj, function(x) wilson_ci(x, n_teae_subj)[2], numeric(1)),
    rate_pct = 100 * rate,
    ci_low_pct = 100 * ci_low,
    ci_high_pct = 100 * ci_high
  ) %>%
  arrange(rate_pct) %>%
  mutate(
    AETERM = factor(AETERM, levels = AETERM),
    y_pos = seq_len(n())
  )

cat("Plot 2 data preview:\n")
print(top10_plot_df)
cat("\n")

p2 <- ggplot(top10_plot_df, aes(x = rate_pct, y = y_pos)) +
  geom_segment(
    aes(x = ci_low_pct, xend = ci_high_pct, y = y_pos, yend = y_pos),
    linewidth = 0.8
  ) +
  geom_segment(
    aes(x = ci_low_pct, xend = ci_low_pct, y = y_pos - 0.12, yend = y_pos + 0.12),
    linewidth = 0.8
  ) +
  geom_segment(
    aes(x = ci_high_pct, xend = ci_high_pct, y = y_pos - 0.12, yend = y_pos + 0.12),
    linewidth = 0.8
  ) +
  geom_point(size = 2.8) +
  scale_y_continuous(
    breaks = top10_plot_df$y_pos,
    labels = top10_plot_df$AETERM
  ) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    subtitle = paste0("n = ", n_teae_subj, " subjects; 95% Wilson CIs"),
    x = "Percentage of Patients (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = "question_3_tlg/output/top10_ae_incidence_ci.png",
  plot = p2,
  width = 10, height = 6.5, dpi = 300
)

cat("Saved Plot 2: question_3_tlg/output/top10_ae_incidence_ci.png\n\n")

cat("Q3 visualization script completed successfully.\n")