
library(tidyverse)
# ==================================================
# Prepare field metadata
# ==================================================

field_metadata <- read.csv(
  "raw/cohort3_epidata.csv",
  na.strings = c("", "NA")
) %>%
  mutate(
    Lab_id = trimws(Lab_id),
    Lab_id = gsub("_", "-", Lab_id),
    timepoint = as.numeric(`Days.at.room.temp`)
  ) %>%
  mutate(
    Lab_id = paste0(Lab_id, "-RDT")
  ) %>%
  select(
    Lab_id,
    timepoint
  ) 
    

# ==================================================
# Prepare sequencing data
# ==================================================

genome_length <- 11923

qc_raw <- read.csv(
  "raw/sequence_data/summary_stats_combined.csv"
) %>%
  # remove non-study samples and failed sample
  filter(!grepl("NTC|NEC|H26|RADDL4B-2026-013|S1-FTA-RT-D1", sample_id),
         # remove the contaminated run
         !(run_id == "run10")
         ) %>%
  mutate(
    # Correct sample naming
    sample_id= if_else(
      sample_id == "N14",
      "N14-D0",
      sample_id
    ))%>%
  mutate(
    # Some lab IDs use underscores instead of dashes
    # Normalize so they parse the same way as the rest
    sample_id = gsub("_", "-", sample_id),
    
    # Field samples: RADDL4B-2026-### (optionally _RDT suffixed)
    # and H##-###-# style IDs
    is_field = grepl("^RADDL|^H[0-9]", sample_id)
  )


# ==================================================
# Field samples
# ==================================================

qc_field <- qc_raw %>%
  filter(is_field) %>%
  mutate(
    # Explicitly identify RDT samples
    sample_type = case_when(
      grepl("-RDT$", sample_id) ~ "RDT",
      TRUE ~ "Brain"
    ),
    temp = "Field"
  ) %>%
  left_join(
    field_metadata,
    by = c("sample_id" = "Lab_id")
  )


# ==================================================
# Lab samples: existing dash-separated parsing
# ==================================================

qc_lab <- qc_raw %>%
  filter(!is_field) %>%
  separate(
    sample_id,
    into = c("sample_id", "sample_type", "temp", "timepoint"),
    sep = "-",
    fill = "right"
  ) %>%
  mutate(
    
    # Brain samples: N6-D0
    sample_type = case_when(
      sample_type == "D0" ~ "Brain",
      TRUE ~ sample_type
    ),
    
    # Experiment 2: N6-FTA-D3
    timepoint_from_temp = as.numeric(
      ifelse(
        grepl("^D", temp),
        sub("D", "", temp),
        NA_character_
      )
    ),
    
    timepoint_from_tp = as.numeric(
      ifelse(
        grepl("^D", timepoint),
        sub("D", "", timepoint),
        NA_character_
      )
    ),
    
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      is.na(timepoint) & grepl("^D", temp) ~ timepoint_from_temp,
      grepl("^D", timepoint) ~ timepoint_from_tp,
      TRUE ~ as.numeric(
        ifelse(
          grepl("^D", timepoint),
          NA_character_,
          timepoint
        )
      )
    ),
    
    # Move D3 out of temp column for Exp 2
    temp = case_when(
      sample_type == "Brain" ~ NA_character_,
      grepl("^D", temp) ~ "RT",
      TRUE ~ temp
    )
  ) %>%
  select(-timepoint_from_temp, -timepoint_from_tp) %>%
  
  # Drop bare IDs with no method/day info to parse
  filter(!is.na(sample_type))


# ==================================================
# Recombine and compute coverage
# ==================================================

qc_all <- bind_rows(qc_lab, qc_field) %>%
  mutate(
    coverage_20 = nonMaskedConsensusCov / genome_length * 100,
    coverage_100 = basesCovered_100 / genome_length * 100
  )


# Check cohort/timepoint structure
table(
  qc_all$cohort,
  qc_all$timepoint,
  useNA = "always"
)

# ---Save output -------------------------------
write.csv(qc_all, "processed/all-runs_sequencing_stats_prepared.csv", row.names=F)

