#==================================================
# Plot genome recovery across all experiments
#==================================================

library(tidyverse)


#==================================================
# Prepare sequencing data
#==================================================

genome_length <- 11923

qc_all <- read.csv(
  "raw/"
) %>%
  filter(!grepl("NTC|NEC", sample_id)) %>%
  select(-sequences) %>%
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
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      is.na(timepoint) & grepl("^D", temp) ~ as.numeric(sub("D", "", temp)),
      grepl("^D", timepoint) ~ as.numeric(sub("D", "", timepoint)),
      TRUE ~ as.numeric(timepoint)
    ),
    
    # Move D3 out of temp column for Exp 2
    temp = case_when(
      sample_type == "Brain" ~ NA_character_,
      grepl("^D", temp) ~ "RT",
      TRUE ~ temp
    ),
    
    coverage_20 = basesCovered_20 / genome_length * 100,
    coverage_100 = basesCovered_100 / genome_length * 100
  )

table(qc_all$Experiment_no, qc_all$timepoint, useNA="always")
#==================================================
# Remove known extraction failure
#==================================================

qc_all <- qc_all %>%
  filter(
    !(sample_id == "S1" &
        sample_type == "FTA" &
        temp == "RT" &
        timepoint == 1)
  )


#==================================================
# Summary across all experiments
#==================================================

recovery_summary_all <- qc_all %>%
  filter(sample_type != "Brain") %>%
  group_by(
    Experiment_no,
    sample_type,
    temp,
    timepoint
  ) %>%
  summarise(
    n = n(),
    mean_coverage20 = mean(coverage_20),
    sd_coverage20 = sd(coverage_20),
    median_coverage20 = median(coverage_20),
    success = sum(coverage_20 >= 90),
    success_percent = success/n*100,
    .groups = "drop"
  )

# write.csv(
#   recovery_summary_all,
#   "figures/all_experiments_genome_recovery_summary.csv",
#   row.names = FALSE
# )


#==================================================
# Plot ≥20x genome recovery
#==================================================

p_all_experiments <- ggplot(
  qc_all %>%
    filter(sample_type != "Brain"),
  aes(
    x = timepoint,
    y = coverage_20,
    colour = sample_type,
    group = interaction(sample_id, sample_type)
  )
) +
  geom_line(
    alpha = 0.1,
    linewidth = 0.4
  ) +
  geom_point(
    size = 2
  ) +
  stat_summary(
    aes(group = sample_type),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(group = sample_type),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  facet_wrap(~ temp) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_y_continuous(
    limits = c(70,100),
    breaks = seq(70,100,5)
  ) +
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method"
  ) +
  theme_classic(base_size = 11);p_all_experiments


ggsave(
  "figures/all_experiments_genome_recovery_coverage20.pdf",
  p_all_experiments,
  width = 16.5,
  height = 14,
  units = "cm",
  device = cairo_pdf
)

coverage_summary <- qc_all %>%
  filter(sample_type != "Brain") %>%
  group_by(temp, sample_type, timepoint) %>%
  summarise(
    mean_coverage = mean(coverage_20, na.rm = TRUE),
    sd_coverage = sd(coverage_20, na.rm = TRUE),
    .groups = "drop"
  )
p_all_experiments <- ggplot(
  coverage_summary,
  aes(
    x = timepoint,
    y = mean_coverage,
    colour = sample_type,
    fill = sample_type,
    group = sample_type
  )
) +
  
  geom_ribbon(
    aes(
      ymin = mean_coverage - sd_coverage,
      ymax = mean_coverage + sd_coverage
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(~ temp) +
  
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  scale_y_continuous(
    breaks = seq(70,100,5)
  ) +
  
  coord_cartesian(
    ylim = c(70,100)
  ) +
  
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method"
  ) +
  
  theme_classic(base_size = 11)

p_all_experiments
ggsave(
  "figures/all_experiments_mean&sd_genome_recovery_coverage20.pdf",
  p_all_experiments,
  width = 16.5,
  height = 14,
  units = "cm",
  device = cairo_pdf
)
#==================================================
# Plot ≥100x genome recovery
#==================================================
p_all_experiments <- ggplot(
  qc_all %>%
    filter(sample_type != "Brain"),
  aes(
    x = timepoint,
    y = coverage_100,
    colour = sample_type,
    group = interaction(sample_id, sample_type)
  )
) +
  geom_line(
    alpha = 0.1,
    linewidth = 0.4
  ) +
  geom_point(
    size = 2
  ) +
  stat_summary(
    aes(group = sample_type),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(group = sample_type),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  facet_wrap(~ temp) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_y_continuous(
    limits = c(70,100),
    breaks = seq(70,100,5)
  ) +
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥100× (%)",
    colour = "Preservation method"
  ) +
  theme_classic(base_size = 11);p_all_experiments


ggsave(
  "figures/all_experiments_genome_recovery_coverage100.pdf",
  p_all_experiments,
  width = 16.5,
  height = 14,
  units = "cm",
  device = cairo_pdf
)

coverage_summary100 <- qc_all %>%
  filter(sample_type != "Brain") %>%
  group_by(temp, sample_type, timepoint) %>%
  summarise(
    mean_coverage = mean(coverage_100, na.rm = TRUE),
    sd_coverage = sd(coverage_100, na.rm = TRUE),
    .groups = "drop"
  )
p_all_experiments <- ggplot(
  coverage_summary100,
  aes(
    x = timepoint,
    y = mean_coverage,
    colour = sample_type,
    fill = sample_type,
    group = sample_type
  )
) +
  
  geom_ribbon(
    aes(
      ymin = mean_coverage - sd_coverage,
      ymax = mean_coverage + sd_coverage
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(~ temp) +
  
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  scale_y_continuous(
    breaks = seq(50,100,5)
  ) +
  
  coord_cartesian(
    ylim = c(50,100)
  ) +
  
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥100× (%)",
    colour = "Preservation method",
    fill = "Preservation method"
  ) +
  
  theme_classic(base_size = 11)

p_all_experiments
ggsave(
  "figures/all_experiments_mean&sd_genome_recovery_coverage100.pdf",
  p_all_experiments,
  width = 16.5,
  height = 14,
  units = "cm",
  device = cairo_pdf
)



# Experiment 1 summary (reference)
exp1_summary <- qc_all %>%
  filter(
    sample_type != "Brain",
    temp == "RT",
    Experiment_no == 1
  ) %>%
  group_by(
    sample_type,
    timepoint
  ) %>%
  summarise(
    mean_coverage = mean(coverage_20, na.rm = TRUE),
    sd_coverage = sd(coverage_20, na.rm = TRUE),
    .groups = "drop"
  )


# Experiment 2 validation samples
exp2_validation <- qc_all %>%
  filter(
    sample_type != "Brain",
    temp == "RT",
    Experiment_no == 2
  )

p_RT_validation <- ggplot() +
  
  # Experiment 1 mean ± SD
  geom_ribbon(
    data = exp1_summary,
    aes(
      x = timepoint,
      ymin = mean_coverage - sd_coverage,
      ymax = mean_coverage + sd_coverage,
      fill = sample_type
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  geom_line(
    data = exp1_summary,
    aes(
      x = timepoint,
      y = mean_coverage,
      colour = sample_type,
      group = sample_type
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = exp1_summary,
    aes(
      x = timepoint,
      y = mean_coverage,
      colour = sample_type
    ),
    size = 3
  ) +
  
  # Experiment 2 validation points
  geom_point(
    data = exp2_validation,
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type
    ),
    shape = 21,
    fill = "white",
    size = 3,
    stroke = 1
  ) +
  
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(85,100,5)
  ) +
  
  coord_cartesian(
    ylim = c(85,100)
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method"
  ) +
  
  theme_classic(base_size = 11)


p_RT_validation
