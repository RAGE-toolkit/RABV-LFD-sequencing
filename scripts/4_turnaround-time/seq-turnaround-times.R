#--------------------------------------------------
# 1. Load packages
#--------------------------------------------------

library(tidyverse)
library(lubridate)
library(Polychrome)
library(patchwork)

epi <- read.csv("raw/21sep26-mimaropa-outbreak-animal.csv")
gen <- read_csv("raw/18aug26-mimaropa-seq.csv", show_col_types = FALSE)
names(epi); names(gen)

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
combined


# Collect dates for key points in pipeline
# Sample collection - RDT test - FTA test - Sequencing
# DATE_SAMPLE_COLLECTION - DATE_LAB_SUBMISSION_FORM - DATE_TESTS_PVO - DATE_TESTS_SPEEDIER - TESTS_RADDL_RDT - DATE_TESTS_RADDL - DATE_TESTS_RITM- Date_sequenced

# Sequencing turnaround
# keep only records with sequencing data
combined_gen <- combined %>%
  filter(!is.na(Date_sequenced))


## This is just cases that have been sequenced
# Stage-specific TATs, as well as total TAT

tat <- combined_gen %>%
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


# Between each stage TATs

tat <- tat %>%
  mutate(
    tat_collection_rdt =
      as.numeric(DATE_RDT - DATE_SAMPLE_COLLECTION),
    tat_rdt_dfat =
      as.numeric(  DATE_DFAT  - DATE_RDT),
    tat_dfat_seq =
      as.numeric( Date_sequenced  - DATE_DFAT),
    # Overall
    tat_total =
      as.numeric(Date_sequenced - DATE_SAMPLE_COLLECTION)
  )

#
summary(tat$tat_dfat_seq)
summary(tat$tat_total)

# Individual sample timelines

timeline <- tat %>%
  select(
    CASE_NO,
    run_name,
    DATE_SAMPLE_COLLECTION,
    DATE_RDT,
    DATE_DFAT,
    Date_sequenced
  ) %>%
  pivot_longer(
    cols = c(
      DATE_SAMPLE_COLLECTION,
      DATE_RDT,
      DATE_DFAT,
      Date_sequenced
    ),
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
        "DATE_DFAT",
        "Date_sequenced"
      ),
      labels = c(
        "Sample collected",
        "RDT",
        "DFAT",
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


nrow(timeline)
# 697 stages recorded

n_distinct(timeline$CASE_NO)
# over 180 unique samples

# timeline <- timeline %>%
#   mutate(year = year(date))

seq_runs <- tat %>%
  filter(!is.na(Date_sequenced), !is.na(run_name)) %>%
  group_by(run_name) %>%
  summarise(
    run_date = min(Date_sequenced),
    n_samples = n_distinct(CASE_NO),
    .groups = "drop"
  )

seq_runs <- tat %>%
  filter(
    !is.na(Date_sequenced),
    !is.na(run_name),
    !is.na(DATE_DFAT)
  ) %>%
  group_by(run_name) %>%
  summarise(
    run_date = min(Date_sequenced),
    n_samples = n_distinct(CASE_NO),
    median_dfat_seq = median(
      as.numeric(Date_sequenced - DATE_DFAT),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(run_date)

seq_runs <- seq_runs %>%
  arrange(run_date) %>%
  mutate(
    days_since_previous_run =
      as.numeric(run_date - lag(run_date))
  )

# batcg size over time
ggplot(seq_runs, aes(x = run_date, y = n_samples)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Sequencing run date",
    y = "Samples sequenced per run"
  ) +
  theme_classic()

# days between runs
ggplot(seq_runs, aes(x = run_date, y = days_since_previous_run)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Sequencing run date",
    y = "Days since previous sequencing run"
  ) +
  theme_classic()

timeline <- timeline %>%
  group_by(CASE_NO) %>%
  mutate(
    collection_date = min(date, na.rm = TRUE),
    collection_year = year(collection_date)
  ) %>%
  ungroup()

# Set a contrasting colour palette for seq runs
run_names <- sort(unique(timeline$run_name[timeline$stage == "Sequencing"]))
run_colours <- setNames(
  palette36.colors(length(run_names)),
  run_names
)

# then plot
timeline_plot <- ggplot(
  timeline,
  aes(
    x = date,
    y = CASE_NO,
    group = CASE_NO
  )
) +
  geom_line(
    colour = "grey70",
    alpha = 0.5,
    linewidth = 0.3
  ) +
  
  # Collection, PVO and Raddl
  geom_point(
    data = timeline %>%
    # filter(stage != "Sequencing"), # to show all diagnostic stages in plot 
    filter(stage == "DFAT"), # show just from point of confirmatory DFAT
    aes(shape = stage),
    size = 1.2
  ) +
  
  # Sequencing points coloured by sequencing run
  geom_point(
    data = timeline %>%
      filter(stage == "Sequencing"),
    aes(
      shape = stage,
      colour = run_name
    ),
    size = 1.2,
    show.legend = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "DFAT" = 17,        
      "Sequencing" = 16     
    )
    
  )+
  scale_colour_manual(
    values = run_colours
  ) +
  labs(
    x = "Date",
    y = "Per sample",
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.02, 0.96),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.5
    ),
    legend.key = element_blank(),
    legend.key.width = unit(0.5, "lines"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2)
  )+
  guides(colour = "none");timeline_plot

#--------------------------------------------------
# TAT data in long format
#--------------------------------------------------

tat_long <- tat %>%
  select(
    CASE_NO,
    tat_collection_rdt,
    tat_rdt_dfat,
    tat_dfat_seq,
    tat_total
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
        "tat_rdt_dfat", 
        "tat_dfat_seq",
        "tat_total"
      ),
      labels = c(
        "Collection → RDT",
        "RDT → DFAT confirmation",
        "DFAT confirmation → Sequence",
        "Sample → Sequence"
      )
    )
  )

tat_summary <- tat_long %>%
  group_by(stage) %>%
  summarise(
    n = sum(!is.na(days)),
    median = median(days, na.rm = TRUE),
    IQR = IQR(days, na.rm = TRUE),
    min = min(days, na.rm = TRUE),
    max = max(days, na.rm = TRUE),
    .groups = "drop"
  )

tat_summary

tat_yearly <- tat %>%
  mutate(
    year = year(DATE_SAMPLE_COLLECTION)
  ) %>%
  group_by(year) %>%
  summarise(
    n = sum(!is.na(tat_total)),
    median_TAT = median(tat_total, na.rm = TRUE),
    Q1 = quantile(tat_total, 0.25, na.rm = TRUE),
    Q3 = quantile(tat_total, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

tat_yearly

tat_quarterly <- tat %>%
  mutate(
    quarter = floor_date(DATE_SAMPLE_COLLECTION, "quarter")
  ) %>%
  group_by(quarter) %>%
  summarise(
    n = sum(!is.na(tat_total)),
    median_TAT = median(tat_total, na.rm = TRUE),
    Q1 = quantile(tat_total, 0.25, na.rm = TRUE),
    Q3 = quantile(tat_total, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


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

# improvement over time:
ggplot(tat, aes(x = DATE_DFAT, y = tat_dfat_seq)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "DFAT confirmation date",
    y = "Overall turnaround time (days)"
  ) +
  theme_classic()

tat_time <- tat %>%
  mutate(
    collection_quarter = floor_date(
      DATE_DFAT,
      unit = "quarter"
    )
  )

tat_quarterly <- tat_time %>%
  group_by(collection_quarter) %>%
  summarise(
    n = n(),
    median_TAT = median(tat_dfat_seq, na.rm = TRUE),
    IQR_low = quantile(tat_dfat_seq, 0.25, na.rm = TRUE),
    IQR_high = quantile(tat_dfat_seq, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  tat_quarterly,
  aes(x = collection_quarter, y = median_TAT)
) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = IQR_low,
      ymax = IQR_high
    ),
    width = 20
  ) +
  labs(
    x = "Sample collection quarter",
    y = "Median overall TAT (days)"
  ) +
  theme_classic()

model <- lm(
  tat_total ~ DATE_DFAT,
  data = tat
)

summary(model)

tat %>%
  mutate(year = year(DATE_DFAT)) %>%
  group_by(year) %>%
  summarise(
    n = n(),
    median_TAT = median(tat_dfat_seq, na.rm = TRUE),
    IQR = IQR(tat_dfat_seq, na.rm = TRUE)
  )

ggplot(
  tat %>% mutate(year = year(DATE_DFAT)),
  aes(x = factor(year), y = tat_dfat_seq)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4
  ) +
  labs(
    x = "Year of sample collection",
    y = "Overall turnaround time (days)"
  ) +
  theme_classic()

tat_panel <- ggplot(
  tat,
  aes(
    x = DATE_DFAT,
    y = tat_dfat_seq
  )
) +
  geom_point(aes(col=run_name),
    alpha = 0.6,
    size = 1.2,
  ) +
  scale_colour_manual(
    values = run_colours
  ) +
  # geom_smooth(
  #   method = "loess",
  #   se = TRUE
  # ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.7,
    col="black"
  )+
  labs(
    x = "Date",
    y = "Turnaround time (days)"
  ) +
  guides(colour = "none")+
  theme_classic(base_size = 14); tat_panel

tat_quarterly <- tat %>%
  mutate(
    quarter = floor_date(DATE_DFAT, "quarter")
  ) %>%
  group_by(quarter) %>%
  summarise(
    n = n(),
    median_TAT = median(tat_dfat_seq, na.rm = TRUE),
    q25 = quantile(tat_dfat_seq, 0.25, na.rm = TRUE),
    q75 = quantile(tat_dfat_seq, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  tat,
  aes(x = DATE_DFAT, y = tat_dfat_seq)
) +
  geom_point(alpha = 0.3) +
  geom_line(
    data = tat_quarterly,
    aes(x = quarter, y = median_TAT),
    linewidth = 1
  ) +
  geom_point(
    data = tat_quarterly,
    aes(x = quarter, y = median_TAT),
    size = 3
  ) +
  labs(
    x = "Sample collection date",
    y = "Sequencing turnaround time (days)"
  ) +
  theme_classic()



figure <- timeline_plot + tat_panel +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      size = 14,
      face = "bold"
    )
  );figure

# tiff
ggsave(
  "figures/mimaropa_timeline_TAT.tiff",
  figure,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

# png
ggsave(
  "figures/mimaropa_timeline_TAT.png",
  figure,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600
)
