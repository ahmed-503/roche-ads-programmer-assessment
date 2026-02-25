# question_3_tlg/01_create_ae_summary_table.R
# Q3 - TLG Adverse Events Summary Table
# Roche ADS Programmer Assessment
#
# Purpose:
#   Create a summary table of treatment-emergent adverse events (TEAEs)
#   using pharmaverseadam::adae and pharmaverseadam::adsl
#
# Table specs:
#   - TEAEs only: TRTEMFL == "Y"
#   - Rows: AETERM (or AESOC; configurable)
#   - Columns: Treatment groups (ACTARM)
#   - Cell values: Count (n) and percentage (%)
#   - Include total column
#   - Sort by descending frequency
#
# Output:
#   question_3_tlg/output/ae_summary_table.html
#   question_3_tlg/logs/q3_table_log.txt

# -----------------------------
# 0) Setup
# -----------------------------
options(stringsAsFactors = FALSE)

dir.create("question_3_tlg/output", recursive = TRUE, showWarnings = FALSE)
dir.create("question_3_tlg/logs", recursive = TRUE, showWarnings = FALSE)

log_file <- "question_3_tlg/logs/q3_table_log.txt"
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Q3 AE Summary Table Log\n")
cat("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------
# 1) Load packages
# -----------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "tibble",
  "pharmaverseadam", "gtsummary", "gt"
)

missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(pharmaverseadam)
library(gtsummary)
library(gt)

cat("Loaded packages successfully.\n\n")

# -----------------------------
# 2) Load input datasets
# -----------------------------
data("adae", package = "pharmaverseadam")
data("adsl", package = "pharmaverseadam")

cat("Input dimensions:\n")
cat("ADAE:", paste(dim(adae), collapse = " x "), "\n")
cat("ADSL:", paste(dim(adsl), collapse = " x "), "\n\n")

# -----------------------------
# 3) Parameters (easy to change)
# -----------------------------
ROW_VAR <- "AETERM"  # change to "AESOC" if you want SOC rows instead
OUT_HTML <- "question_3_tlg/output/ae_summary_table.html"

if (!ROW_VAR %in% c("AETERM", "AESOC")) {
  stop("ROW_VAR must be 'AETERM' or 'AESOC'")
}

cat("Row variable selected:", ROW_VAR, "\n\n")

# -----------------------------
# 4) Prepare denominators from ADSL (one row per subject)
# -----------------------------
# Use ACTARM from ADSL as treatment groups / denominator
adsl_denom <- adsl %>%
  mutate(
    ACTARM = as.character(ACTARM),
    ACTARM = ifelse(is.na(ACTARM) | str_trim(ACTARM) == "", "Missing", ACTARM)
  ) %>%
  distinct(USUBJID, ACTARM)

denom_by_arm <- adsl_denom %>%
  count(ACTARM, name = "N_DENOM") %>%
  arrange(ACTARM)

n_total <- nrow(adsl_denom)

cat("Treatment denominators:\n")
print(denom_by_arm)
cat("Total N:", n_total, "\n\n")

# Optional: gtsummary denominator snapshot (useful for reviewer/log)
denom_tbl <- adsl_denom %>%
  tbl_summary(
    by = ACTARM,
    include = USUBJID,
    statistic = everything() ~ "{n}"
  )

cat("Created denominator summary with gtsummary (for QA/logging).\n\n")

# -----------------------------
# 5) Filter TEAEs and build subject-level event incidence
# -----------------------------
# TEAE records only
# Use distinct subject-event rows so counts are unique subjects (not event records)
adae_teae <- adae %>%
  mutate(
    TRTEMFL = as.character(TRTEMFL),
    ACTARM = as.character(ACTARM),
    ACTARM = ifelse(is.na(ACTARM) | str_trim(ACTARM) == "", "Missing", ACTARM),
    AETERM = as.character(AETERM),
    AESOC = as.character(AESOC)
  ) %>%
  filter(TRTEMFL == "Y")

cat("TEAE records (TRTEMFL == 'Y'):", nrow(adae_teae), "\n")

# Keep only needed columns and non-missing row variable
ae_subj_level <- adae_teae %>%
  transmute(
    USUBJID,
    ACTARM,
    ROW_LABEL = .data[[ROW_VAR]]
  ) %>%
  mutate(ROW_LABEL = ifelse(is.na(ROW_LABEL) | str_trim(ROW_LABEL) == "", "Missing", ROW_LABEL)) %>%
  distinct(USUBJID, ACTARM, ROW_LABEL)

cat("Distinct subject-", ROW_VAR, " rows: ", nrow(ae_subj_level), "\n\n", sep = "")

# -----------------------------
# 6) Count unique subjects by row and treatment
# -----------------------------
counts_by_arm <- ae_subj_level %>%
  count(ROW_LABEL, ACTARM, name = "n")

# Ensure all treatment columns exist for every row (fill 0)
counts_wide <- counts_by_arm %>%
  tidyr::complete(
    ROW_LABEL,
    ACTARM = denom_by_arm$ACTARM,
    fill = list(n = 0)
  ) %>%
  left_join(denom_by_arm, by = "ACTARM") %>%
  mutate(
    pct = ifelse(N_DENOM > 0, 100 * n / N_DENOM, NA_real_),
    cell = sprintf("%d (%.1f%%)", n, pct)
  ) %>%
  select(ROW_LABEL, ACTARM, n, cell)

# Total column across all subjects
counts_total <- ae_subj_level %>%
  count(ROW_LABEL, name = "n_total") %>%
  mutate(
    pct_total = ifelse(n_total > 0, 100 * n_total / n_total, NA_real_),
    Total = sprintf("%d (%.1f%%)", n_total, pct_total)
  ) %>%
  select(ROW_LABEL, n_total, Total)

# Row sorting by descending total frequency
row_order <- counts_total %>%
  arrange(desc(n_total), ROW_LABEL) %>%
  pull(ROW_LABEL)

# Pivot treatment columns to wide display format
tbl_display <- counts_wide %>%
  select(ROW_LABEL, ACTARM, cell) %>%
  tidyr::pivot_wider(
    names_from = ACTARM,
    values_from = cell
  ) %>%
  left_join(counts_total, by = "ROW_LABEL") %>%
  mutate(ROW_LABEL = factor(ROW_LABEL, levels = row_order)) %>%
  arrange(ROW_LABEL) %>%
  mutate(ROW_LABEL = as.character(ROW_LABEL))

# Put row label first, Total second, then treatment columns
treatment_cols <- sort(setdiff(names(tbl_display), c("ROW_LABEL", "Total", "n_total")))
tbl_display <- tbl_display %>%
  select(ROW_LABEL, Total, all_of(treatment_cols))

# -----------------------------
# 7) Create final table (gt) and save HTML
# -----------------------------
# (gtsummary is loaded and used above for denominator QA table;
#  final regulatory-style display here is built as a precise incidence table)
title_txt <- if (ROW_VAR == "AETERM") {
  "Treatment-Emergent Adverse Events by Preferred Term"
} else {
  "Treatment-Emergent Adverse Events by System Organ Class"
}

subtitle_txt <- paste0(
  "Cell values are unique subjects with TEAE: n (%) | Denominator = ADSL subjects by ACTARM | Total N = ",
  n_total
)

gt_tbl <- tbl_display %>%
  gt(rowname_col = "ROW_LABEL") %>%
  tab_header(
    title = md(paste0("**", title_txt, "**")),
    subtitle = subtitle_txt
  ) %>%
  cols_label(
    Total = md(paste0("**Total**<br>N=", n_total))
  )

# Add N to treatment column headers
for (arm in treatment_cols) {
  arm_n <- denom_by_arm %>% filter(ACTARM == arm) %>% pull(N_DENOM)
  gt_tbl <- gt_tbl %>%
    cols_label(!!arm := md(paste0("**", arm, "**<br>N=", arm_n)))
}

gt_tbl <- gt_tbl %>%
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(3)
  ) %>%
  cols_align(align = "left", columns = everything())

gtsave(gt_tbl, OUT_HTML)

cat("AE summary table written:\n")
cat("-", OUT_HTML, "\n\n")

# -----------------------------
# 8) Print preview to log
# -----------------------------
cat("Preview (first 10 rows):\n")
print(utils::head(tbl_display, 10))

cat("\nQ3 summary table script completed successfully.\n")