# question_1_sdtm/01_create_ds_domain.R
# Q1 - SDTM DS Domain Creation
# Roche ADS Programmer Assessment
#
# Purpose:
#   Create SDTM DS domain from pharmaverseraw::ds_raw
#   and study-specific controlled terminology loaded from CSV (sdtm_ct.csv).
#
# Outputs:
#   question_1_sdtm/output/ds.csv
#   question_1_sdtm/output/ds.rds
#   question_1_sdtm/logs/q1_run_log.txt

# -----------------------------
# 0) Setup
# -----------------------------
options(stringsAsFactors = FALSE)

# Create directories
dir.create("question_1_sdtm/output", showWarnings = FALSE)
dir.create("question_1_sdtm/logs", showWarnings = FALSE)

# Print messages in both terminal and save to log file
log_file <- "question_1_sdtm/logs/q1_run_log.txt"
sink(log_file, split = TRUE)

# Add log file title and time stamp
cat("Q1 DS Domain Creation Log\n")
cat("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------
# 1) Load packages
# -----------------------------

# Extra check that required packages are installed before running the script
required_pkgs <- c("dplyr", "stringr", "tibble", "lubridate", "readr", "pharmaverseraw")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them before running this script.")
}

# Load libraries
library(dplyr)
library(stringr)
library(tibble)
library(lubridate)
library(readr)
library(pharmaverseraw)

# Successful loading added to log file
cat("Loaded packages successfully.\n\n")

# -----------------------------
# 2) Load raw input data
# -----------------------------
data("ds_raw", package = "pharmaverseraw")

cat("ds_raw dimensions: ", paste(dim(ds_raw), collapse = " x "), "\n")
cat("ds_raw columns:\n")
print(names(ds_raw))
cat("\nHead of ds_raw:\n")
print(utils::head(ds_raw))

# -----------------------------
# 3) Load study controlled terminology from CSV
# -----------------------------

# file path
ct_file <- "question_1_sdtm/input/sdtm_ct.csv"

# check file availability
if (!file.exists(ct_file)) {
  stop("Controlled terminology file not found: ", ct_file,
       "\nPlace sdtm_ct.csv at: question_1_sdtm/input/sdtm_ct.csv")
}

# read csv file
study_ct <- read.csv(ct_file, stringsAsFactors = FALSE)

# Add logs to log file
cat("\nstudy_ct loaded from:", ct_file, "\n")
cat("study_ct dimensions: ", paste(dim(study_ct), collapse = " x "), "\n")
cat("study_ct columns:\n")
print(names(study_ct))

required_ct_cols <- c(
  "codelist_code",
  "term_code",
  "term_value",
  "collected_value",
  "term_preferred_term",
  "term_synonyms"
)

# validation check of required ct columns
missing_ct_cols <- setdiff(required_ct_cols, names(study_ct))
if (length(missing_ct_cols) > 0) {
  stop("sdtm_ct.csv is missing required columns: ",
       paste(missing_ct_cols, collapse = ", "))
}

# -----------------------------
# 4) Helper functions
# -----------------------------

# transform datetime to SDTM standard
to_iso8601_datetime <- function(date_chr, time_chr = NULL) {
  date_chr <- trimws(as.character(date_chr))
  if (!is.null(time_chr)) time_chr <- trimws(as.character(time_chr))
  
  date_chr[date_chr %in% c("", "NA", "NULL")] <- NA_character_
  if (!is.null(time_chr)) {
    time_chr[time_chr %in% c("", "NA", "NULL")] <- NA_character_
    dt_input <- ifelse(is.na(time_chr), date_chr, paste(date_chr, time_chr))
  } else {
    dt_input <- date_chr
  }
  
  parsed <- suppressWarnings(parse_date_time(
    dt_input,
    orders = c(
      "Ymd HMS", "Ymd HM", "Ymd",
      "dmY HMS", "dmY HM", "dmY",
      "dmy HMS", "dmy HM", "dmy",
      "d b Y HMS", "d b Y HM", "d b Y",
      "d B Y HMS", "d B Y HM", "d B Y",
      "mdY HMS", "mdY HM", "mdY"
    ),
    tz = "UTC"
  ))
  
  out <- rep(NA_character_, length(parsed))
  has_time <- !is.na(parsed) & !is.na(time_chr)
  
  out[!is.na(parsed) & !has_time] <- format(as.Date(parsed[!is.na(parsed) & !has_time]), "%Y-%m-%d")
  out[has_time] <- format(parsed[has_time], "%Y-%m-%dT%H:%M:%S")
  
  out
}

# transform date to SDTM standard
to_iso8601_date <- function(date_chr) {
  date_chr <- trimws(as.character(date_chr))
  date_chr[date_chr %in% c("", "NA", "NULL")] <- NA_character_
  
  parsed <- suppressWarnings(parse_date_time(
    date_chr,
    orders = c("Ymd", "dmY", "dmy", "d b Y", "d B Y", "mdY"),
    tz = "UTC"
  ))
  
  ifelse(is.na(parsed), NA_character_, format(as.Date(parsed), "%Y-%m-%d"))
}

derive_study_day <- function(date_chr, ref_chr) {
  d <- suppressWarnings(as.Date(date_chr))
  r <- suppressWarnings(as.Date(ref_chr))
  out <- rep(NA_integer_, length(d))
  ok <- !is.na(d) & !is.na(r)
  
  out[ok] <- ifelse(d[ok] >= r[ok],
                    as.integer(d[ok] - r[ok]) + 1L,
                    as.integer(d[ok] - r[ok]))
  out
}

# -----------------------------
# 5) Controlled terminology mapping prep for DSDECOD
# -----------------------------

# Lookup table to map to the same standard DSDECOD
ct_map <- bind_rows(
  study_ct %>% transmute(raw_key = str_to_upper(collected_value), DSDECOD_STD = term_value),
  study_ct %>% transmute(raw_key = str_to_upper(term_value), DSDECOD_STD = term_value),
  study_ct %>% transmute(raw_key = str_to_upper(term_synonyms), DSDECOD_STD = term_value),
  study_ct %>% transmute(raw_key = str_to_upper(term_preferred_term), DSDECOD_STD = term_value)
) %>%
  filter(!is.na(raw_key), raw_key != "") %>%
  distinct()

cat("\nControlled terminology mapping rows:", nrow(ct_map), "\n")

# -----------------------------
# 6) Create DS domain
# -----------------------------

# This is the main transformation step to create the SDTM variables
# Use lookup table to:
  # 1) create DS variables from raw data
  # 2) standardize DSDECOD
  # 3) assign DSSEQ number

ds <- ds_raw %>%
  mutate(
    STUDYID = as.character(STUDY),
    DOMAIN = "DS",
    USUBJID = paste0(STUDYID, "-", as.character(PATNUM)),
    
    DSTERM = dplyr::coalesce(as.character(`IT.DSTERM`), as.character(OTHERSP)),
    
    raw_dsdecod_key = str_to_upper(str_trim(
      dplyr::coalesce(
        as.character(`IT.DSDECOD`),
        as.character(`IT.DSTERM`),
        as.character(OTHERSP)
      )
    )),
    
    DSDTC = to_iso8601_datetime(DSDTCOL, DSTMCOL),
    DSSTDTC = to_iso8601_date(`IT.DSSTDAT`),
    
    VISIT = dplyr::coalesce(as.character(FORML), as.character(INSTANCE)),
    VISITNUM = suppressWarnings(readr::parse_number(
      dplyr::coalesce(as.character(INSTANCE), as.character(FORML))
    )),
    
    DSCAT = "DISPOSITION EVENT"
  ) %>%
  left_join(ct_map, by = c("raw_dsdecod_key" = "raw_key")) %>%
  mutate(
    DSDECOD = dplyr::coalesce(
      DSDECOD_STD,
      str_to_upper(as.character(`IT.DSDECOD`)),
      str_to_upper(DSTERM)
    )
  ) %>%
  group_by(USUBJID) %>%
  arrange(
    suppressWarnings(as.Date(substr(DSSTDTC, 1, 10))),
    suppressWarnings(as.Date(substr(DSDTC, 1, 10))),
    VISITNUM,
    .by_group = TRUE
  ) %>%
  mutate(DSSEQ = row_number()) %>%
  ungroup()

# -----------------------------
# 7) Derive DSSTDY
# -----------------------------
# Fallback reference date:
# Use first non-missing DSSTDTC within subject as subject-level ref date.
# (If DM RFSTDTC were available, we'd use that instead.)
# derives DSSTDY = the study day number for each DS record
# tells you how many days into the study the disposition event happened.

ds <- ds %>%
  group_by(USUBJID) %>%
  mutate(.ref_dt = min(as.Date(DSSTDTC), na.rm = TRUE)) %>%
  ungroup()

ds$.ref_dt[is.infinite(ds$.ref_dt)] <- NA

ds <- ds %>%
  mutate(
    DSSTDY = derive_study_day(DSSTDTC, as.character(.ref_dt))
  )

# -----------------------------
# 8) Final column order + checks
# -----------------------------
ds_final <- ds %>%
  select(
    STUDYID, DOMAIN, USUBJID, DSSEQ,
    DSTERM, DSDECOD, DSCAT,
    VISITNUM, VISIT,
    DSDTC, DSSTDTC, DSSTDY
  )

cat("\nFinal DS dimensions: ", paste(dim(ds_final), collapse = " x "), "\n")
cat("\nFinal DS preview:\n")
print(utils::head(ds_final, 10))

required_vars <- c("STUDYID","DOMAIN","USUBJID","DSSEQ","DSTERM","DSDECOD","DSCAT",
                   "VISITNUM","VISIT","DSDTC","DSSTDTC","DSSTDY")
missing_vars <- setdiff(required_vars, names(ds_final))
if (length(missing_vars) > 0) {
  stop("Missing required variables in final DS: ", paste(missing_vars, collapse = ", "))
}

cat("\nRequired variables present: YES\n")

# -----------------------------
# 9) Save outputs
# -----------------------------
write.csv(ds_final, "question_1_sdtm/output/ds.csv", row.names = FALSE)
saveRDS(ds_final, "question_1_sdtm/output/ds.rds")

cat("\nOutput files written:\n")
cat("- question_1_sdtm/output/ds.csv\n")
cat("- question_1_sdtm/output/ds.rds\n")

cat("\nQ1 script completed successfully.\n")
sink()