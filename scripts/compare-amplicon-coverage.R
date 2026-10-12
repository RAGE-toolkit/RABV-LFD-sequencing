library(tidyverse)

# Read the output made by Python script
depths <- read_tsv("processed/all_seq_combined_amplicon_depths.tsv")


###############################################################
##### Individual depth plots per cohort 
#################################################################

# Turn sample columns into rows and extract metadata from their names
depths_long <- depths %>%
  pivot_longer(
    cols = -c(chrom, amplicon),
    names_to = "sample_name",
    values_to = "mean_depth"
  ) %>%
  filter(
    !grepl("NTC|NEC|H26", sample_name)
  ) %>%
  mutate(
    # Correct sample naming
    sample_name = if_else(
      sample_name == "N14",
      "N14-D0",
      sample_name
    ),
    
    # Normalise underscores to dashes
    sample_name_norm = str_replace_all(sample_name, "_", "-"),
    
    is_field = str_detect(
      sample_name_norm,
      "^RADDL|^H[0-9]"
    )
  ) %>%
  extract(
    sample_name_norm,
    into = c("sample", "condition", "temp", "day"),
    regex = "^(N\\d+|S\\d+)(?:-(FTA|RDT))?(?:-(RT|37))?(?:-(D\\d+))?$",
    remove = FALSE
  ) %>%
  mutate(
    
    # ---------------------------------------------------------
    # FIELD SAMPLES — Cohort 3
    # Brain and RDT from the same underlying field sample
    # e.g. RADDL4B-2026-024
    #      RADDL4B-2026-024_RDT
    # ---------------------------------------------------------
    
    sample = if_else(
      is_field,
      str_remove(sample_name, "_RDT$"),
      sample
    ),
    
    condition = if_else(
      is_field,
      if_else(
        str_detect(sample_name, "_RDT$"),
        "RDT",
        "Brain"
      ),
      condition
    ),
    
    temp = if_else(
      is_field,
      "Field",
      temp
    ),
    
    # ---------------------------------------------------------
    # General metadata
    # ---------------------------------------------------------
    
    condition = if_else(
      is.na(condition) | condition == "",
      "Baseline",
      condition
    ),
    
    day = factor(
      day,
      levels = c("D0", "D1", "D3", "D7", "D21")
    ),
    
    # ---------------------------------------------------------
    # Define cohorts
    # ---------------------------------------------------------
    
    cohort = case_when(
      str_detect(sample, "^S") ~ "Cohort 1",
      str_detect(sample, "^N") ~ "Cohort 2",
      is_field ~ "Cohort 3",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-sample_name_norm)

# define the cohorts 
cohort1 <- depths_long %>% filter(cohort == "Cohort 1")
cohort2 <- depths_long %>% filter(cohort == "Cohort 2")
cohort3 <- depths_long %>% filter(cohort == "Cohort 3")


# cohort 1 plots
ggplot(
  cohort1,
  aes(
    amplicon,
    mean_depth,
    colour = condition,
    linetype = day,
    group = sample_name
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  scale_linetype_manual(
    values = c(
      "D0" = "solid",
      "D1" = "dashed",
      "D3" = "dotted",
      "D7" = "dotdash",
      "D21" = "longdash"
    ),
    na.value = "solid"
  ) +
  geom_line() +
  facet_wrap(~ sample, scales = "free_y") +
  theme_classic()


#cohort 3 plots
ggplot(
  cohort3,
  aes(
    amplicon,
    mean_depth,
    colour = condition,
    group = sample_name
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_line() +
  facet_wrap(~ sample, scales = "free_y") +
  theme_classic()

# cohort 2 plots
ggplot(
  cohort2,
  aes(
    amplicon,
    mean_depth,
    colour = condition,
    linetype = day,
    group = sample_name
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  scale_linetype_manual(
    values = c(
      "D0" = "solid",
      "D1" = "dashed",
      "D3" = "dotted",
      "D7" = "dotdash",
      "D21" = "longdash"
    ),
    na.value = "solid"
  ) +
  geom_line() +
  facet_wrap(~ sample, scales = "free_y") +
  theme_classic()



###################################
#Consistency of amplicon failure 
####################################
amplicon_failure <- depths_long %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  group_by(cohort, amplicon) %>%
  summarise(
    n_samples = n_distinct(sample_name),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  )

# amplicon failure by cohort
ggplot(
  amplicon_failure,
  aes(
    x = factor(amplicon),
    y = prop_below_20x
  )
) +
  geom_col() +
  facet_wrap(~ cohort) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Amplicon",
    y = "Proportion of samples <20×"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


## Facet by time
amplicon_failure <- depths_long %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  filter(!is.na(day)) %>%
  group_by(day, amplicon) %>%
  summarise(
    n_samples = n_distinct(sample_name),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  )

# amplicon failure by time
ggplot(
  amplicon_failure,
  aes(
    x = factor(amplicon),
    y = prop_below_20x
  )
) +
  geom_col() +
  facet_wrap(~ day) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Amplicon",
    y = "Proportion of samples <20×"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

# amplicons with >20% of samples failing through time
amplicon_failure_time <- depths_long %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  filter(
    !is.na(day)
   # ,
    # day != "D0"
  ) %>%
  group_by(day, amplicon) %>%
  summarise(
    n_samples = sum(!is.na(mean_depth)),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  ) %>%
  group_by(day) %>%
  summarise(
    n_amplicons_gt20 = sum(prop_below_20x > 0.20),
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(day, levels = c("D0","D1", "D3", "D7", "D21"))
  )

ggplot(
  amplicon_failure_time,
  aes(x = day, y = n_amplicons_gt20, group = 1)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(
    breaks = 0:41
  ) +
  labs(
    x = "Storage time (days)",
    y = "Amplicons with >20% of samples failing"
  ) +
  theme_classic()


# overall amplicon failure proportions
amplicon_failure_all <- depths_long %>%
  # filter(condition != "Baseline", !is.na(condition)) %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  group_by(amplicon) %>%
  summarise(
    n_samples = n_distinct(sample_name),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  )

amplicon_failure_all %>%
  summarise(
    n_amplicons_gt20 = sum(prop_below_20x > 0.20),
    n_amplicons = n(),
    proportion = mean(prop_below_20x > 0.20)
  )

p_amplicon_failure <- ggplot(
  amplicon_failure_all,
  aes(
    x = factor(amplicon),
    y = prop_below_20x
  )
) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Amplicon",
    y = "Samples with <20× coverage"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  );p_amplicon_failure 

# Publication-ready plot
p_amplicon_failure <- p_amplicon_failure  +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ); p_amplicon_failure 


ggsave(
  "figures/amplicon_failure.pdf",
  plot = p_amplicon_failure,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  "figures/amplicon_failure.png",
  plot = p_amplicon_failure,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# proportion of >20 versus <20 coverage per amplicon

amplicon_coverage <- depths_long %>%
  mutate(
    coverage_category = case_when(
      mean_depth < 20 ~ "<20×",
      mean_depth > 20 ~ ">20x"
    )
  ) %>%
  count(amplicon, coverage_category) %>%
  group_by(amplicon) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()
ggplot(
  amplicon_coverage,
  aes(
    x = factor(amplicon),
    y = prop,
    fill = coverage_category
  )
) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Amplicon",
    y = "Proportion of samples",
    fill = "Mean depth"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


### Proportion of samples per preservation medium 

amplicon_failure <- depths_long %>%
  mutate(
    condition = case_when(
      condition %in% c("Brain", "Baseline") ~ "Brain",
      TRUE ~ condition
    ),
    day = case_when(
      condition == "Brain" ~ "D0",
      TRUE ~ day
    ),
    below_20x = mean_depth < 20
  ) %>%
  # filter(
  #   cohort == "Cohort 1",
  #   day != "D0",
  #   condition %in% c("FTA", "RDT")
  # ) %>%
  group_by(amplicon, condition, day) %>%
  summarise(
    n_samples = sum(!is.na(mean_depth)),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(
      day,
      levels = c("D1", "D3", "D7", "D21")
    )
  )


# dodged barplot
ggplot(
  amplicon_failure,
  aes(
    x = factor(amplicon),
    y = prop_below_20x,
    fill = condition
  )
) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Amplicon",
    y = "Proportion of samples <20×",
    fill = "Specimen type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


# look at amplicon failure over time, plot each amplicon as facet 
amplicon_failure <- depths_long %>%
  filter(
    cohort == "Cohort 1",
    day != "D0",
    condition %in% c("FTA", "RDT")
  ) %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  group_by(amplicon, condition, day) %>%
  summarise(
    n_samples = sum(!is.na(mean_depth)),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(
      day,
      levels = c("D1", "D3", "D7", "D21")
    )
  )

# amplicon failures over time, irrespective of medium
ggplot(
  amplicon_failure,
  aes(
    x = day,
    y = prop_below_20x,
    group = 1
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ amplicon) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Storage time",
    y = "Proportion of samples <20×"
  ) +
  theme_classic()


# RDT v FTA
ggplot(
  amplicon_failure,
  aes(
    x = day,
    y = prop_below_20x,
    colour = condition,
    group = condition
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ amplicon) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Storage time",
    y = "Proportion of samples <20×",
    colour = "Specimen type"
  ) +
  theme_classic()



##failure over time including baselines 
amplicon_failure <- depths_long %>%
  filter(
    cohort == "Cohort 1",
   # condition %in% c("FTA", "RDT")
  ) %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  group_by(amplicon, condition, day) %>%
  summarise(
    n_samples = sum(!is.na(mean_depth)),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_samples,
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(
      day,
      levels = c("D0", "D1", "D3", "D7", "D21")
    )
  )

# Manuscript-quality figure
p <- ggplot(
  amplicon_failure,
  aes(
    x = day,
    y = prop_below_20x,
    colour = condition,
    group = condition
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ amplicon) +
  scale_x_discrete(
    labels = c(
      "D0" = "0",
      "D1" = "1",
      "D3" = "3",
      "D7" = "7",
      "D21" = "21"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(-0.05, 1),
    breaks = seq(0, 1, 0.2)
  )+
  scale_colour_discrete(
    labels = c(
      "Baseline" = "Brain",
      "FTA" = "FTA",
      "RDT" = "RDT"
    )
  ) +
  labs(
    x = "Storage time (days)",
    y = "Proportion of samples <20×",
    colour = "Specimen type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  );p

# PNG — 600 dpi
ggsave(
  "figures/amplicon_failure_over_time.png",
  plot = p,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# TIFF — 600 dpi, LZW compressed
ggsave(
  "figures/amplicon_failure_over_time.tiff",
  plot = p,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# overall amplicon failure
amplicon_failure_summary <- depths_long %>%
  filter(
    cohort == "Cohort 1"
  ) %>%
  mutate(
    below_20x = mean_depth < 20
  ) %>%
  group_by(condition, day) %>%
  summarise(
    n_amplicons = sum(!is.na(mean_depth)),
    n_below_20x = sum(below_20x, na.rm = TRUE),
    prop_below_20x = n_below_20x / n_amplicons,
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(day, levels = c("D0", "D1", "D3", "D7", "D21"))
  )

p2 <- ggplot(
  amplicon_failure_summary,
  aes(
    x = day,
    y = prop_below_20x,
    colour = condition,
    group = condition
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_discrete(
    labels = c(
      "D0" = "0",
      "D1" = "1",
      "D3" = "3",
      "D7" = "7",
      "D21" = "21"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, 0.05)
  )+
  scale_colour_discrete(
    labels = c(
      "Baseline" = "Brain",
      "FTA" = "FTA",
      "RDT" = "LFD"
    )
  ) +
  labs(
    x = "Storage time (days)",
    y = "Amplicons with <20× mean depth",
    colour = "Specimen type"
  ) +
  theme_classic(base_size = 12)
p2

# PNG — 600 dpi
ggsave(
  "figures/summary_amplicon_failure_over_time.png",
  plot = p2,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# TIFF — 600 dpi, LZW compressed
ggsave(
  "figures/summary_amplicon_failure_over_time.tiff",
  plot = p2,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
