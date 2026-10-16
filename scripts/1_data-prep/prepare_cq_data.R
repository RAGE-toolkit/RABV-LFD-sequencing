
library(tidyverse)

# ==================================================
# Prepare Cq data- combine with sequencing stats data
# ==================================================
cq_raw <- read_csv(
  "raw/qPCR_lyssavirus_Cq.csv"
)
cq_raw <-cq_raw %>%
  mutate(
    # Correct sample naming
    sample_id= if_else(
      Sample_ID == "N14",
      "N14-D0",
      Sample_ID
    ))%>%
  mutate(
    # Some lab IDs use underscores instead of dashes
    # Normalize so they parse the same way as the rest
    sample_id = gsub("_", "-", sample_id),
    # Remove dash between letter and number at start of ID
    # e.g. N-10-D0 -> N10-D0
    sample_id = gsub("^([A-Za-z])-([0-9])", "\\1\\2", sample_id)
  )

cq_all  <- cq_raw  %>%
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

## import and correct the field sample ids
field <- read_csv("raw/cohort3_epidata.csv") 

cq_all <- cq_all %>%
  left_join(
    field %>%
      select(Case_id, Lab_id),
    by = c("sample_id" = "Case_id")
  ) %>%
  mutate(
    sample_id = coalesce(Lab_id, sample_id)
  ) %>%
  select(-Lab_id)%>%
  mutate(
    # Field samples: RADDL4B-2026-### (optionally _RDT suffixed)
    # and H##-###-# style IDs
    is_field = grepl("^RADDL|^H[0-9]", sample_id)
  )


# want to add to sequencing data
seq <- read.csv("processed/all-runs_sequencing_stats_prepared.csv")

seq_cq <- seq %>%
  left_join(
    cq_all %>%
      select(
        sample_id,
        sample_type,
        temp,
        timepoint,
        everything()
      ),
    by = c(
      "sample_id",
      "sample_type",
      "temp",
      "timepoint"
    )
  )

write.csv(seq_cq, "processed/all-runs_sequencing_stats_prepared_withCq.csv", row.names=F)