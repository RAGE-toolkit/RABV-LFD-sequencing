#--------------------------------------------------
# 1. Load packages
#--------------------------------------------------

library(tidyverse)
library(lubridate)


epi <- read_csv("processed/RDT_validation_cases.csv")
gen <- read_csv("raw/26aug26_Genomic_surveillance_R4B_sequencing_data - Clean_sequence_data.csv", show_col_types = FALSE)
names(epi); names(gen)

#--------------------------------------------------
# 2. PREP DATA
#--------------------------------------------------

# Remove columns that are completely empty
gen <- gen %>%
  select(where(~ !all(is.na(.))))

# Correct for matching issue, gen has 3 digit code, gen has 1-3 digit 
gen <- gen %>%
  mutate(
    prefix = str_extract(case_number, "^[A-Za-z]+-"),
    number = str_extract(case_number, "\\d+$"),
    case_number = paste0(prefix, str_pad(number, width = 3, pad = "0"))
  ) %>%
  select(-prefix, -number)

# Combine datsets based on case id
combined <- epi %>%
  left_join(gen, by = c("case_no" = "case_number"))
dim(combined)


#--------------------------------------------------
# 3. Collect dates for key points in pipeline
#--------------------------------------------------
# DATE_SAMPLE_COLLECTION - DATE_LAB_SUBMISSION_FORM - DATE_TESTS_PVO - DATE_TESTS_SPEEDIER - TESTS_RADDL_RDT - DATE_TESTS_RADDL - DATE_TESTS_RITM- Date_sequenced

# Split into:
# a. Date sample collected
# b. Date of RDT (1st performed across PVO/MAO/RADDL)
# c. Date of DFAT (either RADDL or RITM)

# All diag turnaround stats (not just sequenced cases)

tat <- combined %>%
  mutate(
    date_sample_collection = as.Date(date_sample_collection, format = "%m/%d/%y"),
    date_tests_pvo= as.Date(date_tests_pvo, format = "%m/%d/%y"),
    date_tests_mao = as.Date(date_tests_mao, format = "%m/%d/%y"),
    date_lab_submission_form = as.Date(date_lab_submission_form, format = "%m/%d/%y"),
    date_tests_speedier = as.Date(date_tests_speedier, format = "%m/%d/%y"),
    Date_test_RADDL = as.Date(Date_test_RADDL, format = "%m/%d/%y"),
    date_tests_ritm = as.Date(date_tests_ritm, format = "%m/%d/%y"),
    Date_sequenced = as.Date(Date_sequenced, format = "%d-%b-%y")
  ) %>%
  mutate(
    date_rdt = pmin(
      if_else(
        tests_pvo %in% c("RDT+", "RDT-"),
        date_tests_pvo,
        as.Date(NA)
      ),
      if_else(
        tests_mao %in% c("RDT+", "RDT-"),
        date_tests_mao,
        as.Date(NA)
      ),
      if_else(
        tests_raddl_rdt %in% c("RDT+", "RDT-"),
        date_tests_raddl,
        as.Date(NA)
      ),
      na.rm = TRUE
    ),
    
    date_dfat = pmin(
      if_else(
        tests_raddl_d_fat %in% c("dFAT+", "dFAT-"),
        date_tests_raddl,
        as.Date(NA)
      ),
      if_else(
        tests_ritm_d_fat %in% c("dFAT+", "dFAT-"),
        date_tests_ritm,
        as.Date(NA)
      ),
      na.rm = TRUE
    )
  )

#--------------------------------------------------
# Between each stage TATs
# Collection - RDT - DFAST
#--------------------------------------------------

tat <- tat %>%
  mutate(
    tat_collection_rdt =
      as.numeric(date_rdt - date_sample_collection),
    tat_rdt_dfat =
      as.numeric(  date_dfat  - date_rdt),
  )
#--------------------------------------------------
# TAT data in long format
#--------------------------------------------------
tat_long <- tat %>%
  select(
    case_no,
    tat_collection_rdt,
    tat_rdt_dfat
  ) %>%
  pivot_longer(
    cols = -case_no,
    names_to = "stage",
    values_to = "days"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "tat_collection_rdt",
        "tat_rdt_dfat"
      ),
      labels = c(
        "Collection → RDT",
        "RDT → DFAT confirmation"
      )
    )
  ) %>%
  group_by(stage) %>%
  filter(any(!is.na(days))) %>%
  ungroup()

# For now assume that cases where RDT results before collection date are incorrect and replace with 0 (just two cases)
tat_long$days[which(tat_long$days<0)]=0


# Individual sample timelines

timeline <- tat %>%
  select(
    case_no,
    date_sample_collection,
    date_rdt,
    date_dfat
  ) %>%
  pivot_longer(
    cols = -case_no,
    names_to = "stage",
    values_to = "date"
  ) %>%
  filter(!is.na(date)) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "date_sample_collection",
        "date_rdt",
        "date_dfat"
      ),
      labels = c(
        "Sample collected",
        "RDT",
        "Confirmatory DFAT"
      )
    )
  )


# Plot individual sample timelines

# first order by collection dates
timeline <- timeline %>%
  group_by(case_no) %>%
  mutate(
    collection_date = min(date, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    case_no = reorder(case_no, collection_date)
  )


nrow(timeline)
# 933 stages recorded

n_distinct(timeline$case_no)
# 319 unique samples

# then plot
ggplot(
  timeline,
  aes(x = date, y = case_no, group = case_no)
) +
  geom_line(alpha = 0.3) +
  geom_point(aes(shape = stage, col = stage), size = 1) +
  scale_colour_manual(
    values = c(
      "Sample collected" = "#009E73",
      "RDT" = "#CC79A7",
      "Confirmatory DFAT" = "#E69F00"
    )
  ) +
  labs(
    x = "Date",
    y = "Sample",
    colour = "Pipeline stage",
    shape = "Pipeline stage"
  ) +
  theme_classic()+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

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

