#--------------------------------------------------
# 1. Load packages
#--------------------------------------------------

library(tidyverse)
library(lubridate)


epi <- read.csv("raw/18aug26-mimaropa-outbreak-animal.csv")
gen <- read_csv("raw/18aug26-mimaropa-seq.csv", show_col_types = FALSE)
names(epi); names(gen)

# Remove columns that are completely empty
gen <- gen %>%
  select(where(~ !all(is.na(.))))

# Combine datsets based on case id
combined <- epi %>%
  left_join(gen, by = c("CASE_NO" = "case_number"))
combined

# Collect dates for key points in pipeline
# Sample collection - RDT test - FTA test - Sequencing
# DATE_SAMPLE_COLLECTION - DATE_LAB_SUBMISSION_FORM - DATE_TESTS_PVO - DATE_TESTS_SPEEDIER - TESTS_RADDL_RDT - DATE_TESTS_RADDL - DATE_TESTS_RITM- Date_sequenced

# Sequencing turnaround
# keep only records with sequencing data
combined_gen <- combined %>%
  filter(!is.na(Date_sequenced))
tat <- combined_gen %>%
  mutate(
    DATE_SAMPLE_COLLECTION = as.Date(DATE_SAMPLE_COLLECTION, format = "%m/%d/%y"),
    Date_sequenced = as.Date(Date_sequenced, format = "%d-%b-%y"),
    turnaround_days = as.numeric(
      Date_sequenced - DATE_SAMPLE_COLLECTION
    )
  )


## This is just cases that have been sequenced
# Stage-specific TATs, as well as total TAT

tat <- combined_gen %>%
  mutate(
    DATE_SAMPLE_COLLECTION = as.Date(DATE_SAMPLE_COLLECTION, format = "%m/%d/%y"),
    DATE_TESTS_PVO = as.Date(DATE_TESTS_PVO, format = "%m/%d/%y"),
    DATE_LAB_SUBMISSION_FORM = as.Date(DATE_LAB_SUBMISSION_FORM, format = "%m/%d/%y"),
    DATE_TESTS_SPEEDIER = as.Date(DATE_TESTS_SPEEDIER, format = "%m/%d/%y"),
    DATE_TESTS_RADDL = as.Date(DATE_TESTS_RADDL, format = "%m/%d/%y"),
    DATE_TESTS_RITM = as.Date(DATE_TESTS_RITM, format = "%m/%d/%y"),
    Date_sequenced = as.Date(Date_sequenced, format = "%d-%b-%y")
  )

# Between each stage TATs

tat <- tat %>%
  mutate(
    tat_collection_pvo =
      as.numeric(DATE_TESTS_PVO - DATE_SAMPLE_COLLECTION),
    
    tat_pvo_submission =
      as.numeric(DATE_LAB_SUBMISSION_FORM - DATE_TESTS_PVO),
    
    tat_submission_speedier =
      as.numeric(DATE_TESTS_SPEEDIER - DATE_LAB_SUBMISSION_FORM),
    
    tat_speedier_raddl =
      as.numeric(DATE_TESTS_RADDL - DATE_TESTS_SPEEDIER),
    
    tat_raddl_ritm =
      as.numeric(DATE_TESTS_RITM - DATE_TESTS_RADDL),
    
    tat_ritm_seq =
      as.numeric(Date_sequenced - DATE_TESTS_RITM),
    
    tat_total =
      as.numeric(Date_sequenced - DATE_SAMPLE_COLLECTION)
  )
summary(tat$tat_total)

# Individual sample timelines

timeline <- tat %>%
  select(
    CASE_NO,
    DATE_SAMPLE_COLLECTION,
    DATE_TESTS_PVO,
    DATE_LAB_SUBMISSION_FORM,
    DATE_TESTS_SPEEDIER,
    DATE_TESTS_RADDL,
    DATE_TESTS_RITM,
    Date_sequenced
  ) %>%
  pivot_longer(
    cols = -CASE_NO,
    names_to = "stage",
    values_to = "date"
  ) %>%
  filter(!is.na(date)) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "DATE_SAMPLE_COLLECTION",
        "DATE_TESTS_PVO",
        "DATE_LAB_SUBMISSION_FORM",
        "DATE_TESTS_SPEEDIER",
        "DATE_TESTS_RADDL",
        "DATE_TESTS_RITM",
        "Date_sequenced"
      ),
      labels = c(
        "Sample collected",
        "PVO",
        "Lab submission",
        "Speedier",
        "Raddl",
        "RITM",
        "Sequencing"
      )
    )
  )


# Plot individual sample timelines

# first order by collection dates
timeline <- timeline %>%
  group_by(CASE_NO) %>%
  mutate(
    collection_date = min(date, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    CASE_NO = reorder(CASE_NO, collection_date)
  )

# 50 samples with sequence data:
nrow(timeline)
# 227 stages recorded

n_distinct(timeline$CASE_NO)
# over 50 unique samples

# then plot
ggplot(
  timeline,
  aes(x = date, y = CASE_NO, group = CASE_NO)
) +
  geom_line(alpha = 0.4) +
  geom_point(aes(shape = stage), size = 3) +
  labs(
    x = "Date",
    y = "Sample",
    shape = "Pipeline stage"
  ) +
  theme_classic()


# TAT data in long format

tat_long <- tat %>%
  select(
    CASE_NO,
    tat_collection_pvo,
    tat_pvo_submission,
    tat_submission_speedier,
    tat_speedier_raddl,
    tat_raddl_ritm,
    tat_ritm_seq
  ) %>%
  pivot_longer(
    cols = -CASE_NO,
    names_to = "stage",
    values_to = "days"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "tat_collection_pvo",
        "tat_pvo_submission",
        "tat_submission_speedier",
        "tat_speedier_raddl",
        "tat_raddl_ritm",
        "tat_ritm_seq"
      ),
      labels = c(
        "Collection → PVO",
        "PVO → Lab submission",
        "Lab submission → Speedier",
        "Speedier → Raddl",
        "Raddl → RITM",
        "RITM → Sequencing"
      )
    )
  ) %>%
  # Remove stages with no observed TATs
  group_by(stage) %>%
  filter(any(!is.na(days))) %>%
  ungroup()

tat_summary <- tat_long %>%
  group_by(stage) %>%
  summarise(
    n = sum(!is.na(days)),
    median = median(days, na.rm = TRUE),
    mean = mean(days, na.rm = TRUE),
    IQR = IQR(days, na.rm = TRUE),
    SD = sd(days, na.rm = TRUE),
    min = min(days, na.rm = TRUE),
    max = max(days, na.rm = TRUE),
    .groups = "drop"
  )

tat_summary

# Plot TAT between stages

ggplot(tat_long, aes(x = stage, y = days)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    na.rm = TRUE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Turnaround time (days)"
  ) +
  theme_classic()
