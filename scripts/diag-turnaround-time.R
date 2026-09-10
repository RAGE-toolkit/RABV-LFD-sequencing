#--------------------------------------------------
# 1. Load packages
#--------------------------------------------------

library(tidyverse)
library(lubridate)


epi <- read_csv("raw/26aug26_MIMAROPA outbreak 2022-2026 - ANIMAL.csv",  skip = 1)
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
  left_join(gen, by = c("CASE_NO" = "case_number"))
dim(combined)

# Limit to just cases with RDT and DFAT tests
combined <- combined %>%
  filter(
    # RDT result
    (
      TESTS_PVO %in% c("RDT+", "RDT-") |
        TESTS_RADDL_RDT %in% c("RDT+", "RDT-") |
        TESTS_MAO %in% c("RDT+", "RDT-") 
    ) &
      
      # dFAT result
      # dFAT/FAT result
      (
        TESTS_RADDL_dFAT %in% c("dFAT+", "dFAT-", "FAT+", "FAT-") |
          TESTS_RITM_dFAT %in% c("dFAT+", "dFAT-", "FAT+", "FAT-")
      )
  )

dim(combined)

# Then remove the dfat unfit cases
combined <- combined %>%
  filter(!str_detect(LAB_RESULT_COPY, "unfit"))

dim(combined)

write.csv(combined,"processed/diagnostic-validation-cases_v3.csv", row.names=F)

## SHOULD BE 320 - only getting 318

# need to remove OM-134 as lab dates are incorrect
combined <- combined %>%
  filter(!CASE_NO == "OM-134")
# 317


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
    DATE_SAMPLE_COLLECTION = as.Date(DATE_SAMPLE_COLLECTION, format = "%m/%d/%y"),
    DATE_TESTS_PVO = as.Date(DATE_TESTS_PVO, format = "%m/%d/%y"),
    DATE_TESTS_MAO = as.Date(DATE_TESTS_MAO, format = "%m/%d/%y"),
    DATE_LAB_SUBMISSION_FORM = as.Date(DATE_LAB_SUBMISSION_FORM, format = "%m/%d/%y"),
    DATE_TESTS_SPEEDIER = as.Date(DATE_TESTS_SPEEDIER, format = "%m/%d/%y"),
    DATE_TESTS_RADDL = as.Date(DATE_TESTS_RADDL, format = "%m/%d/%y"),
    DATE_TESTS_RITM = as.Date(DATE_TESTS_RITM, format = "%m/%d/%y"),
    Date_sequenced = as.Date(Date_sequenced, format = "%d-%b-%y")
  ) %>%
  mutate(
    DATE_RDT = pmin(
      if_else(
        TESTS_PVO %in% c("RDT+", "RDT-"),
        DATE_TESTS_PVO,
        as.Date(NA)
      ),
      if_else(
        TESTS_MAO %in% c("RDT+", "RDT-"),
        DATE_TESTS_MAO,
        as.Date(NA)
      ),
      if_else(
        TESTS_RADDL_RDT %in% c("RDT+", "RDT-"),
        DATE_TESTS_RADDL,
        as.Date(NA)
      ),
      na.rm = TRUE
    ),
    
    DATE_DFAT = pmin(
      if_else(
        TESTS_RADDL_dFAT %in% c("dFAT+", "dFAT-"),
        DATE_TESTS_RADDL,
        as.Date(NA)
      ),
      if_else(
        TESTS_RITM_dFAT %in% c("dFAT+", "dFAT-"),
        DATE_TESTS_RITM,
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
      as.numeric(DATE_RDT - DATE_SAMPLE_COLLECTION),
    tat_rdt_dfat =
      as.numeric(  DATE_DFAT  - DATE_RDT),
  )
#--------------------------------------------------
# TAT data in long format
#--------------------------------------------------
tat_long <- tat %>%
  select(
    CASE_NO,
    tat_collection_rdt,
    tat_rdt_dfat
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
    CASE_NO,
    DATE_SAMPLE_COLLECTION,
    DATE_RDT,
    DATE_DFAT
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
        "DATE_RDT",
        "DATE_DFAT"
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
  group_by(CASE_NO) %>%
  mutate(
    collection_date = min(date, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    CASE_NO = reorder(CASE_NO, collection_date)
  )

# 512 samples with sequence data:
nrow(timeline)
# 1753 stages recorded

n_distinct(timeline$CASE_NO)
# over 512 unique samples

# then plot
ggplot(
  timeline,
  aes(x = date, y = CASE_NO, group = CASE_NO)
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
