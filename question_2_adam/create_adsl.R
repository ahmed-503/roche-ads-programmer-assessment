# question_2_adam/create_adsl.R
# Q2 - ADaM ADSL Dataset Creation
# Roche ADS Programmer Assessment
#
# Purpose:
#   Create an ADSL (subject-level) dataset using SDTM source datasets
#   from pharmaversesdtm and derive the requested variables:
#   AGEGR9, AGEGR9N, TRTSDTM, TRTSTMF, ITTFL, LSTAVLDT
#
# Inputs:
#   pharmaversesdtm::dm, pharmaversesdtm::vs, pharmaversesdtm::ex,
#   pharmaversesdtm::ds, pharmaversesdtm::ae
#
# Outputs:
#   question_2_adam/output/adsl.csv
#   question_2_adam/output/adsl.rds
#   question_2_adam/logs/q2_run_log.txt

# -----------------------------
# 0) Setup
# -----------------------------
options(stringsAsFactors = FALSE)

dir.create("question_2_adam/output", recursive = TRUE, showWarnings = FALSE)
dir.create("question_2_adam/logs", recursive = TRUE, showWarnings = FALSE)

log_file <- "question_2_adam/logs/q2_run_log.txt"
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Q2 ADSL Creation Log\n")
cat("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------
# 1) Load packages
# -----------------------------
required_pkgs <- c("dplyr", "stringr", "tibble", "tidyr", "lubridate", "pharmaversesdtm", "readr")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them before running this script.")
}

library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(lubridate)
library(pharmaversesdtm)
library(readr)

cat("Loaded packages successfully.\n\n")

# -----------------------------
# 2) Load input SDTM datasets
# -----------------------------
data("dm", package = "pharmaversesdtm")
data("vs", package = "pharmaversesdtm")
data("ex", package = "pharmaversesdtm")
data("ds", package = "pharmaversesdtm")
data("ae", package = "pharmaversesdtm")

cat("Input dimensions:\n")
cat("DM:", paste(dim(dm), collapse = " x "), "\n")
cat("VS:", paste(dim(vs), collapse = " x "), "\n")
cat("EX:", paste(dim(ex), collapse = " x "), "\n")
cat("DS:", paste(dim(ds), collapse = " x "), "\n")
cat("AE:", paste(dim(ae), collapse = " x "), "\n\n")

# -----------------------------
# 3) Helper functions
# -----------------------------

# Parse ISO-like SDTM datetime strings and impute missing time components for EXSTDTC
# Rules (per assignment):
# - Date must be complete (YYYY-MM-DD) to derive TRTSDTM
# - If no time: impute 00:00:00, flag "H"
# - If HH only: impute :00:00, flag "M"
# - If HH:MM only: impute :00, no TRTSTMF if only seconds are missing
# - If HH:MM:SS complete: no flag
# Returns a list with:
#   dtm_chr = ISO datetime string
#   tmf     = imputation flag (custom: H or M, else NA)
impute_exstdtc <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NULL")] <- NA_character_
  
  out_dtm <- rep(NA_character_, length(x))
  out_tmf <- rep(NA_character_, length(x))
  
  # Require complete datepart YYYY-MM-DD
  has_complete_date <- !is.na(x) & grepl("^\\d{4}-\\d{2}-\\d{2}", x)
  
  # Full datetime HH:MM:SS
  idx_full <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}$", x)
  out_dtm[idx_full] <- x[idx_full]
  
  # Missing seconds only HH:MM
  idx_hm <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}$", x)
  out_dtm[idx_hm] <- paste0(x[idx_hm], ":00")
  # Per spec: if only seconds missing, do NOT populate TRTSTMF
  
  # Hour only HH
  idx_h <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}$", x)
  out_dtm[idx_h] <- paste0(x[idx_h], ":00:00")
  out_tmf[idx_h] <- "M"
  
  # Date only (no time)
  idx_date <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  out_dtm[idx_date] <- paste0(x[idx_date], "T00:00:00")
  out_tmf[idx_date] <- "H"
  
  # Also handle space separator instead of T (rare fallback)
  idx_space_hms <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}$", x)
  out_dtm[idx_space_hms] <- sub(" ", "T", x[idx_space_hms])
  
  idx_space_hm <- has_complete_date & grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}$", x)
  out_dtm[idx_space_hm] <- paste0(sub(" ", "T", x[idx_space_hm]), ":00")
  # no flag for seconds-only imputation
  
  tibble(dtm_chr = out_dtm, tmf = out_tmf)
}

# Convert ISO-like datetime/date strings to Date using complete datepart only
extract_datepart <- function(x) {
  x <- as.character(x)
  date_part <- ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}", x), substr(x, 1, 10), NA_character_)
  suppressWarnings(as.Date(date_part))
}

# -----------------------------
# 4) Start ADSL from DM (base dataset)
# -----------------------------
# Per assignment: use DM as the basis of ADSL
adsl <- dm

cat("Initialized ADSL from DM.\n\n")

# -----------------------------
# 5) Derive AGEGR9 and AGEGR9N
# -----------------------------
adsl <- adsl %>%
  mutate(
    AGE_NUM = suppressWarnings(as.numeric(AGE)),
    AGEGR9 = case_when(
      !is.na(AGE_NUM) & AGE_NUM < 18 ~ "<18",
      !is.na(AGE_NUM) & AGE_NUM <= 50 ~ "18 - 50",
      !is.na(AGE_NUM) & AGE_NUM > 50 ~ ">50",
      TRUE ~ NA_character_
    ),
    AGEGR9N = case_when(
      AGEGR9 == "<18" ~ 1,
      AGEGR9 == "18 - 50" ~ 2,
      AGEGR9 == ">50" ~ 3,
      TRUE ~ NA_real_
    )
  )

cat("Derived AGEGR9 and AGEGR9N.\n")

# -----------------------------
# 6) Derive ITTFL (randomized flag)
# -----------------------------
# Assignment rule: Y if DM.ARM populated, else N
adsl <- adsl %>%
  mutate(
    ITTFL = ifelse(!is.na(ARM) & str_trim(as.character(ARM)) != "", "Y", "N")
  )

cat("Derived ITTFL.\n")

# -----------------------------
# 7) Derive TRTSDTM and TRTSTMF from EX (first valid dose)
# -----------------------------
# Valid dose:
#   EXDOSE > 0 OR (EXDOSE == 0 AND EXTRT contains 'PLACEBO')
ex_prep <- ex %>%
  mutate(
    EXDOSE_NUM = suppressWarnings(as.numeric(EXDOSE)),
    valid_dose = (!is.na(EXDOSE_NUM) & EXDOSE_NUM > 0) |
      (!is.na(EXDOSE_NUM) & EXDOSE_NUM == 0 &
         str_detect(str_to_upper(coalesce(as.character(EXTRT), "")), "PLACEBO"))
  ) %>%
  filter(valid_dose)

# Impute time for EXSTDTC only when datepart is complete
ex_imputed <- ex_prep %>%
  bind_cols(impute_exstdtc(ex_prep$EXSTDTC)) %>%
  mutate(
    ex_date = extract_datepart(dtm_chr),
    ex_dtm_posix = suppressWarnings(ymd_hms(dtm_chr, tz = "UTC"))
  ) %>%
  filter(!is.na(ex_date), !is.na(ex_dtm_posix))

# First valid exposure per subject -> TRTSDTM/TRTSTMF
trts <- ex_imputed %>%
  arrange(USUBJID, ex_dtm_posix) %>%
  group_by(USUBJID) %>%
  summarise(
    TRTSDTM = first(dtm_chr),
    TRTSTMF = first(tmf),
    .groups = "drop"
  )

# Also derive last valid exposure datetime/date for LSTAVLDT support
trte <- ex_imputed %>%
  arrange(USUBJID, ex_dtm_posix) %>%
  group_by(USUBJID) %>%
  summarise(
    TRTEDTM = last(dtm_chr),
    TRTEDT = last(ex_date),
    .groups = "drop"
  )

adsl <- adsl %>%
  left_join(trts, by = "USUBJID") %>%
  left_join(trte, by = "USUBJID")

cat("Derived TRTSDTM, TRTSTMF, and TRTEDTM/TRTEDT from EX.\n")

# -----------------------------
# 8) Derive LSTAVLDT components from VS, AE, DS, EX
# -----------------------------
# (1) VS: last complete date with a valid result
# valid result = VSSTRESN and VSSTRESC not both missing
vs_last <- vs %>%
  mutate(
    vs_date = extract_datepart(VSDTC),
    has_vs_result = !(is.na(VSSTRESN) & (is.na(VSSTRESC) | str_trim(as.character(VSSTRESC)) == ""))
  ) %>%
  filter(has_vs_result, !is.na(vs_date)) %>%
  group_by(USUBJID) %>%
  summarise(VS_LASTDT = max(vs_date, na.rm = TRUE), .groups = "drop")

# (2) AE: last complete AE start date
ae_last <- ae %>%
  mutate(ae_date = extract_datepart(AESTDTC)) %>%
  filter(!is.na(ae_date)) %>%
  group_by(USUBJID) %>%
  summarise(AE_LASTDT = max(ae_date, na.rm = TRUE), .groups = "drop")

# (3) DS: last complete disposition date
ds_last <- ds %>%
  mutate(ds_date = extract_datepart(DSSTDTC)) %>%
  filter(!is.na(ds_date)) %>%
  group_by(USUBJID) %>%
  summarise(DS_LASTDT = max(ds_date, na.rm = TRUE), .groups = "drop")

# (4) EX: use TRTEDT (already from valid-dose treatment records)
ex_last <- trte %>%
  select(USUBJID, EX_LASTDT = TRTEDT)

# Merge and set LSTAVLDT = max of available dates
lstavldt_df <- adsl %>%
  select(USUBJID) %>%
  distinct() %>%
  left_join(vs_last, by = "USUBJID") %>%
  left_join(ae_last, by = "USUBJID") %>%
  left_join(ds_last, by = "USUBJID") %>%
  left_join(ex_last, by = "USUBJID") %>%
  rowwise() %>%
  mutate(
    LSTAVLDT = {
      candidates <- c(VS_LASTDT, AE_LASTDT, DS_LASTDT, EX_LASTDT)
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates) == 0) as.Date(NA) else max(candidates)
    }
  ) %>%
  ungroup() %>%
  select(USUBJID, LSTAVLDT)

adsl <- adsl %>%
  left_join(lstavldt_df, by = "USUBJID")

cat("Derived LSTAVLDT from VS, AE, DS, and EX.\n")

# -----------------------------
# 9) Final checks and output selection
# -----------------------------
# Keep all DM variables plus requested derived variables
# (reviewers usually like seeing full ADSL base + key derivations)
required_deriv <- c("AGEGR9", "AGEGR9N", "TRTSDTM", "TRTSTMF", "ITTFL", "LSTAVLDT")
missing_deriv <- setdiff(required_deriv, names(adsl))
if (length(missing_deriv) > 0) {
  stop("Missing required derived variables: ", paste(missing_deriv, collapse = ", "))
}

cat("\nADSL dimensions: ", paste(dim(adsl), collapse = " x "), "\n")
cat("Preview of key derived variables:\n")
print(
  adsl %>%
    select(USUBJID, AGE, AGEGR9, AGEGR9N, ARM, ITTFL, TRTSDTM, TRTSTMF, TRTEDTM, LSTAVLDT) %>%
    head(10)
)

# -----------------------------
# 10) Save outputs
# -----------------------------
write.csv(adsl, "question_2_adam/output/adsl.csv", row.names = FALSE)
saveRDS(adsl, "question_2_adam/output/adsl.rds")

cat("\nOutput files written:\n")
cat("- question_2_adam/output/adsl.csv\n")
cat("- question_2_adam/output/adsl.rds\n")
cat("\nQ2 script completed successfully.\n")