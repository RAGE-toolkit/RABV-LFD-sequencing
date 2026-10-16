library(tidyverse)

#==================================================
# Plot genome recovery across all experiments
#==================================================

qc_all <- read.csv("processed/all-runs_sequencing_stats_prepared.csv")
pal <- wes_palette("Zissou1", 5, type = "discrete")

pal
fta_col <- pal[1]
rdt_col <- pal[4]

#==================================================
# Summary across all experiments
#==================================================

recovery_summary_all <- qc_all %>%
  filter(sample_type != "Brain") %>%
  group_by(
    cohort,
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
#   "tables/all_experiments_genome_recovery_summary.csv",
#   row.names = FALSE
# )


#==================================================
# Plot ≥20x genome recovery
#==================================================

p_all_experiments <- ggplot(
  qc_all %>%
    filter(sample_type != "Brain") %>%
    mutate(
      temp = factor(temp, levels = c("RT", "37", "Field"))
    ),
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
  scale_colour_manual(
    values = c(
      "FTA" = "#3B9AB2",
      "RDT" = "#E1AF00"
    )
  ) +
  facet_wrap(
    ~ temp,
    scales = "free_x"
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_y_continuous(
    limits = c(55, 100),
    breaks = seq(55, 100, 5)
  ) +
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method"
  ) +
  theme_classic(base_size = 11)

p_all_experiments


# ggsave(
#   "figures/all_experiments_genome_recovery_coverage20.pdf",
#   p_all_experiments,
#   width = 16.5,
#   height = 14,
#   units = "cm",
#   device = cairo_pdf
# )

coverage_summary <- qc_all %>%
  filter(sample_type != "Brain") %>%
  group_by(temp, sample_type, timepoint) %>%
  summarise(
    n = sum(!is.na(coverage_20)),
    mean_coverage = mean(coverage_20, na.rm = TRUE),
    sd_coverage = sd(coverage_20, na.rm = TRUE),
    .groups = "drop"
  )
p_all_experiments2 <- ggplot(
  coverage_summary %>%
    mutate(
      temp = factor(temp, levels = c("RT", "37", "Field"))
    ),
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
      ymin = pmax(mean_coverage - sd_coverage, 55),
      ymax = pmin(mean_coverage + sd_coverage, 100)
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
  # geom_point(
  #   aes(size = n)
  # ) +
  # scale_size_continuous(
  #   range = c(2.5, 5),
  #   name = "Number of samples"
  # ) +
  geom_text(
    aes(
      x = timepoint,
      y = mean_coverage,
      label = n,
      colour = sample_type
    ),
    vjust = -0.8,
    size = 3,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ temp,
    scales = "free_x"
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_y_continuous(
    breaks = seq(55, 100, 5)
  ) +
  coord_cartesian(
    ylim = c(55, 100)
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#3B9AB2",
      "RDT" = "#E1AF00"
    )
  ) +
  scale_fill_manual(
    values = c(
      "FTA" = "#3B9AB2",
      "RDT" = "#E1AF00"
    )
  ) +
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.08, 0.12),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      colour = "grey80"
    )
  );p_all_experiments2

# ggsave(
#   "figures/all_experiments_mean&sd_genome_recovery_coverage20-nLabels.pdf",
#   p_all_experiments2,
#   width = 16.5,
#   height = 14,
#   units = "cm",
#   device = cairo_pdf
# )


##have cohort 2 as validation points instead
# Cohort 1: summary
coverage_summary_c1 <- qc_all %>%
  filter(
    cohort == 1,
    sample_type != "Brain"
  ) %>%
  group_by(temp, sample_type, timepoint) %>%
  summarise(
    n = sum(!is.na(coverage_20)),
    mean_coverage = mean(coverage_20, na.rm = TRUE),
    sd_coverage = sd(coverage_20, na.rm = TRUE),
    .groups = "drop"
  )

# Cohort 2: individual validation observations
coverage_c2 <- qc_all %>%
  filter(
    cohort == 2,
    sample_type != "Brain"
  )

# Cohort 3: individual field observations
coverage_c3 <- qc_all %>%
  filter(
    cohort == 3,
    sample_type != "Brain"
  )

p_all_experiments3_pub <- ggplot(
  coverage_summary_c1 %>%
    mutate(
      temp = factor(temp, levels = c("RT", "37", "Field")),
      cohort = "Cohort 1"
    ),
  aes(
    x = timepoint,
    y = mean_coverage,
    colour = sample_type,
    fill = sample_type,
    group = sample_type
  )
) +
  
  # Cohort 1: mean ± SD
  geom_ribbon(
    aes(
      ymin = pmax(mean_coverage - sd_coverage, 55),
      ymax = pmin(mean_coverage + sd_coverage, 100)
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  # Cohort 2: validation points
  geom_point(
    data = coverage_c2 %>%
      mutate(
        temp = factor(temp, levels = c("RT", "37", "Field")),
        cohort = "Cohort 2"
      ),
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type,
      shape = cohort
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 1.5,
    fill = "white",
    alpha = 0.5,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  
  # Cohort 1: mean line
  geom_line(
    linewidth = 1
  ) +
  
  # Cohort 1: filled circles
  geom_point(
    aes(shape = cohort),
    size = 3,
    colour = "black"
  ) +
  
  # Cohort 3: field validation points
  geom_point(
    data = coverage_c3 %>%
      mutate(
        temp = factor(temp, levels = c("RT", "37", "Field")),
        cohort = "Cohort 3"
      ),
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type,
      shape = cohort
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2.5,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  
  # Panel labels
  facet_wrap(
    ~ temp,
    scales = "free_x",
    labeller = as_labeller(c(
      "RT" = "Ambient",
      "37" = "37°C",
      "Field" = "Field"
    ))
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  scale_y_continuous(
    breaks = seq(70, 100, 5)
  ) +
  
  coord_cartesian(
    ylim = c(70, 100)
  ) +
  
  # Preservation method
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    ),
    labels = c(
      "FTA" = "FTA card",
      "RDT" = "LFD"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    ),
    labels = c(
      "FTA" = "FTA card",
      "RDT" = "LFD"
    )
  ) +
  
  # Cohort shapes — preserved
  scale_shape_manual(
    values = c(
      "Cohort 1" = 21,  # filled circle
      "Cohort 2" = 1,   # open circle
      "Cohort 3" = 17   # triangle
    )
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method",
    shape = "Cohort"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    # Axis text
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    
    # Panel labels
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    # Legend
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0.12, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "transparent",
      colour = "grey80"
    ),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.8, "cm")
  ) +
  
  guides(
    colour = guide_legend(
      # title = "Preservation\nmedium",
      title =NULL,
      override.aes = list(
        shape = 16,
        fill = c("#0072B2", "#D55E00"),
        colour = c("#0072B2", "#D55E00"),
        alpha = 1,
        size = 3
      )
    ),
    fill = "none",
    shape = guide_legend(
      # title = "Cohort",
      title =NULL,
      override.aes = list(
        alpha = c(0.8, 0.4, 0.8)
      )
    )
  )

p_all_experiments3_pub
# ggsave(
#   "figures/p_all_experiments3_publication.pdf",
#   p_all_experiments3_pub,
#   width = 180,
#   height = 125,
#   units = "mm",
#   device = cairo_pdf
# )
# ggsave(
#   "figures/p_all_experiments3_publication.png",
#   p_all_experiments3_pub,
#   width = 180,
#   height = 125,
#   units = "mm",
#   dpi = 600,
#   bg = "white"
# )



## with brain 
# Brain: individual observations
coverage_brain <- qc_all %>%
  filter(
    sample_type == "Brain",
    timepoint == 0
  ) %>%
  mutate(
    temp = factor(temp, levels = c("RT", "37", "Field"))
  )

p_all_experiments4_pub <- ggplot(
  coverage_summary_c1 %>%
    mutate(
      temp = factor(temp, levels = c("RT", "37", "Field")),
      cohort = "Cohort 1"
    ),
  aes(
    x = timepoint,
    y = mean_coverage,
    colour = sample_type,
    fill = sample_type,
    group = sample_type
  )
) +
  
  # Cohort 1: mean ± SD
  geom_ribbon(
    aes(
      ymin = pmax(mean_coverage - sd_coverage, 55),
      ymax = pmin(mean_coverage + sd_coverage, 100)
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  # Cohort 2: validation points
  geom_point(
    data = coverage_c2 %>%
      mutate(
        temp = factor(temp, levels = c("RT", "37", "Field")),
        cohort = "Cohort 2"
      ),
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type,
      shape = cohort
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 1.5,
    fill = "white",
    alpha = 0.5,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  
  # Cohort 1: mean line
  geom_line(
    linewidth = 1
  ) +
  
  # Cohort 1: filled circles
  geom_point(
    aes(shape = cohort),
    size = 3,
    colour = "black"
  ) +
  
  # Cohort 3: field validation points
  geom_point(
    data = coverage_c3 %>%
      mutate(
        temp = factor(temp, levels = c("RT", "37", "Field")),
        cohort = "Cohort 3"
      ),
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type,
      shape = cohort
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2.5,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = coverage_brain %>%
      mutate(
        temp = factor(temp, levels = c("RT", "37", "Field"))
      ),
    aes(
      x = timepoint,
      y = coverage_20,
      shape = "Brain"
    ),
    shape = 23,
    size = 3,
    fill = "grey40",
    colour = "black",
    stroke = 0.8,
    inherit.aes = FALSE
  ) +
  
  # Panel labels
  facet_wrap(
    ~ temp,
    scales = "free_x",
    labeller = as_labeller(c(
      "RT" = "Room temperature",
      "37" = "37°C",
      "Field" = "Field"
    ))
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  scale_y_continuous(
    breaks = seq(55, 100, 5)
  ) +
  
  coord_cartesian(
    ylim = c(55, 100)
  ) +
  
  # Preservation method
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
  
  # Cohort shapes — preserved
  scale_shape_manual(
    values = c(
      "Cohort 1" = 21,  # filled circle
      "Cohort 2" = 1,   # open circle
      "Cohort 3" = 17   # triangle
    )
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method",
    shape = "Cohort"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    # Axis text
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    
    # Panel labels
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    # Legend
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0.12, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "transparent",
      colour = "grey80"
    ),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.8, "cm")
  ) +
  
  guides(
    colour = guide_legend(
      # title = "Preservation\nmedium",
      title =NULL,
      override.aes = list(
        shape = 16,
        fill = c("#0072B2", "#D55E00"),
        colour = c("#0072B2", "#D55E00"),
        alpha = 1,
        size = 3
      )
    ),
    fill = "none",
    shape = guide_legend(
      # title = "Cohort",
      title =NULL,
      override.aes = list(
        alpha = c(0.8, 0.4, 0.8)
      )
    )
  )

p_all_experiments4_pub 
# Experiment 1 summary (reference)
exp1_summary <- qc_all %>%
  filter(
    sample_type != "Brain",
    temp == "RT",
    cohort == 1
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
    cohort == 2
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
    breaks = seq(55,100,5)
  ) +
  
  coord_cartesian(
    ylim = c(55,100)
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation method",
    fill = "Preservation method"
  ) +
  
  theme_classic(base_size = 11)


p_RT_validation



coverage_by_type <- qc_all %>%
  filter(
    !is.na(coverage_20),
    sample_type %in% c("Brain", "FTA", "RDT")
  ) %>%
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    )
  )

p_brain_comparison <- ggplot(
  coverage_by_type,
  aes(
    x = sample_type,
    y = coverage_20,
    fill = sample_type
  )
) +
  
  geom_boxplot(
    width = 0.6,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      "Brain" = "grey50",
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_y_continuous(
    limits = c(55, 100),
    breaks = seq(0, 100, 20)
  ) +
  
  labs(
    x = "Specimen type",
    y = "Genome coverage ≥20× (%)"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

p_brain_comparison


coverage_by_type <- qc_all %>%
  filter(
    cohort == 1,
    !is.na(coverage_20),
    sample_type %in% c("Brain", "FTA", "RDT")
  ) %>%
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    )
  )
p_brain_comparison

coverage_by_type <- qc_all %>%
  filter(
    !is.na(coverage_20),
    sample_type %in% c("Brain", "FTA", "RDT")
  ) %>%
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    ),
    cohort = factor(
      cohort,
      levels = c(1, 2, 3),
      labels = c("Cohort 1", "Cohort 2", "Cohort 3")
    )
  )

p_brain_comparison <- ggplot(
  coverage_by_type,
  aes(
    x = sample_type,
    y = coverage_20,
    fill = sample_type
  )
) +
  
  geom_boxplot(
    width = 0.6,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  facet_wrap(~ cohort) +
  
  scale_fill_manual(
    values = c(
      "Brain" = "grey50",
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  
  scale_y_continuous(
    limits = c(55, 100),
    breaks = seq(0, 100, 20)
  ) +
  
  labs(
    x = "Specimen type",
    y = "Genome coverage ≥20× (%)"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "none"
  )

p_brain_comparison

coverage_thresholds <- qc_all %>%
  filter(
    sample_type %in% c("Brain", "FTA", "RDT"),
    !is.na(coverage_20)
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(1, 2, 3),
      labels = c("Cohort 1", "Cohort 2", "Cohort 3")
    ),
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    )
  ) %>%
  group_by(cohort, sample_type) %>%
  summarise(
    n = n(),
    n_ge80 = sum(coverage_20 >= 80),
    pct_ge80 = 100 * n_ge80 / n,
    n_ge90 = sum(coverage_20 >= 90),
    pct_ge90 = 100 * n_ge90 / n,
    .groups = "drop"
  ) %>%
  mutate(
    `≥80% coverage` = sprintf(
      "%d/%d (%.1f%%)",
      n_ge80, n, pct_ge80
    ),
    `≥90% coverage` = sprintf(
      "%d/%d (%.1f%%)",
      n_ge90, n, pct_ge90
    )
  ) %>%
  select(
    Cohort = cohort,
    `Specimen type` = sample_type,
    `≥80% coverage`,
    `≥90% coverage`
  )

coverage_thresholds
coverage_thresholds %>%
  select(
    Cohort,
    `Specimen type`,
    `≥90% coverage`,
    `≥80% coverage`
  )

coverage_success <- qc_all %>%
  filter(
    sample_type %in% c("Brain", "FTA", "RDT"),
    !is.na(coverage_20),
    !(cohort == 3 & sample_type == "FTA")
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(1, 2, 3),
      labels = c("Cohort 1", "Cohort 2", "Cohort 3")
    ),
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    )
  ) %>%
  group_by(cohort, sample_type) %>%
  summarise(
    n = n(),
    ge80_n = sum(coverage_20 >= 80),
    ge90_n = sum(coverage_20 >= 90),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(ge80_n, ge90_n),
    names_to = "threshold",
    values_to = "success"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("ge80_n", "ge90_n"),
      labels = c("≥80% coverage", "≥90% coverage")
    ),
    percentage = 100 * success / n
  )

coverage_success <- qc_all %>%
  filter(
    sample_type %in% c("Brain", "FTA", "RDT"),
    !is.na(coverage_20),
    !(cohort == 3 & sample_type == "FTA")
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(1, 2, 3),
      labels = c("Cohort 1", "Cohort 2", "Cohort 3")
    ),
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    )
  ) %>%
  group_by(cohort, sample_type) %>%
  summarise(
    n = n(),
    ge80_n = sum(coverage_20 >= 80),
    ge90_n = sum(coverage_20 >= 90),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(ge80_n, ge90_n),
    names_to = "threshold",
    values_to = "success"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("ge80_n", "ge90_n"),
      labels = c("≥80% coverage", "≥90% coverage")
    ),
    percentage = 100 * success / n
  )

p_coverage_success <- ggplot(
  coverage_success,
  aes(
    x = sample_type,
    y = percentage,
    fill = sample_type,
    alpha = threshold
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  # geom_text(
  #   aes(
  #     label = if_else(
  #       percentage == 0 & threshold == "≥90% coverage",
  #       "0%",
  #       ""
  #     )
  #   ),
  #   position = position_dodge(width = 0.75),
  #   vjust = -0.4,
  #   size = 3.5
  # ) +
  geom_text(
    aes(
      label = paste0(success, "/", n)
    ),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_grid(
    . ~ cohort,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "Brain" = "grey50",
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    ),
    labels = c(
      "Brain" = "Brain",
      "FTA" = "FTA card",
      "RDT" = "LFD"
    )
  ) +
  scale_alpha_manual(
    values = c(
      "≥80% coverage" = 0.45,
      "≥90% coverage" = 1
    )
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Specimen type",
    y = "Proportion of samples (%)",
    fill = "Specimen type",
    alpha = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "top"
  )
p_coverage_success

# # Publication-ready theme
# p_coverage_success_pub <- p_coverage_success +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 13),
#     strip.text = element_text(size = 15, face = "bold"),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 13)
#   )
# 
# # PNG
# ggsave(
#   "figures/p_coverage_success.png",
#   plot = p_coverage_success_pub,
#   width = 180,
#   height = 100,
#   units = "mm",
#   dpi = 600,
#   bg = "white"
# )
# 
# # PDF
# ggsave(
#   "figures/p_coverage_success.pdf",
#   plot = p_coverage_success_pub,
#   width = 180,
#   height = 100,
#   units = "mm",
#   device = cairo_pdf,
#   bg = "white"
# )
# 
# # TIFF
# ggsave(
#   "figures/p_coverage_success.tiff",
#   plot = p_coverage_success_pub,
#   width = 180,
#   height = 100,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw",
#   bg = "white"
# )

