library(tidyverse)
library(viridis)

#---------------------------
# 1. Read and parse
#---------------------------
df <- read.csv("raw/FTA_RDT_study_seq_depth_raw.csv")

df_parsed <- df %>%
  separate(Sample, into = c("sample_id", "sample_type", "temp", "timepoint"),
           sep = "-", fill = "right") %>%
  mutate(
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      !is.na(timepoint) ~ as.numeric(gsub("D", "", timepoint)),
      TRUE ~ NA_real_
    ),
    temp = if_else(sample_type == "Brain", NA_character_, temp),
    Depth = as.numeric(Depth)
  )

#---------------------------
# 2. Prepare baseline per sample
#---------------------------
baseline_df <- df_parsed %>%
  filter(sample_type == "Brain") %>%
  select(sample_id, Position, Depth) %>%
  crossing(temp = c("RT", "37")) %>%
  rename(Depth_baseline = Depth)

plot_df <- df_parsed %>%
  filter(sample_type != "Brain") %>%
  left_join(baseline_df, by = c("sample_id", "Position", "temp"))

#---------------------------
# 3. Summarise QC metrics
#---------------------------
qc_summary <- df_parsed %>%
  group_by(sample_id, sample_type, temp, timepoint) %>%
  summarise(
    genome_length = n(),
    mean_depth = mean(Depth, na.rm = TRUE),
    pct_10x = mean(Depth >= 10) * 100,
    pct_20x = mean(Depth >= 20) * 100,
    pct_50x = mean(Depth >= 50) * 100,
    pct_100x = mean(Depth >= 100) * 100,
    .groups = "drop"
  )

qc_summary_condition <- qc_summary %>%
  filter(sample_type != "Brain") %>%
  group_by(sample_type, temp, timepoint) %>%
  summarise(
    across(
      c(mean_depth, pct_10x:pct_100x),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      )
    ),
    n_samples = n(),
    .groups = "drop"
  )

#---------------------------
# 4. Reshape for plotting multiple coverage thresholds
#---------------------------
coverage_cols <- c("pct_10x_mean", "pct_20x_mean","pct_50x_mean", "pct_100x_mean")
sd_cols <- c("pct_10x_sd",  "pct_20x_sd", "pct_50x_sd", "pct_100x_sd")

qc_long <- qc_summary_condition %>%
  pivot_longer(
    cols = all_of(coverage_cols),
    names_to = "coverage_metric",
    values_to = "mean_value"
  ) %>%
  mutate(
    coverage = factor(str_extract(coverage_metric, "\\d+"),
                      levels = c("10", "20,","50", "100"),
                      labels = c("10x", "20x","50x", "100x"))
  ) %>%
  # add corresponding SD
  mutate(sd_value = case_when(
    coverage_metric == "pct_10x_mean" ~ pct_10x_sd,
    coverage_metric == "pct_20x_mean" ~ pct_20x_sd,
    coverage_metric == "pct_50x_mean" ~ pct_50x_sd,
    coverage_metric == "pct_100x_mean" ~ pct_100x_sd
  ))

#---------------------------
# 5. Define colour palette
#---------------------------
palette <- c("FTA" = "#1b9e77", "RDT" = "#d95f02")

#---------------------------
# 6. Plot function
#---------------------------
plot_mean_sd <- function(df, y_var, y_label, file_prefix, facet_var = "temp") {
  p <- ggplot(df, aes(x = timepoint, y = !!sym(y_var),
                      colour = sample_type, group = sample_type)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = !!sym(y_var) - sd_value,
                    ymax = !!sym(y_var) + sd_value,
                    fill = sample_type),
                alpha = 0.2, colour = NA) +
    facet_wrap(as.formula(paste("~", facet_var))) +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(
      x = "Days post sampling",
      y = y_label,
      colour = "Sample type",
      fill = "Sample type",
      title = paste(y_label, "over time")
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 12))
  
  ggsave(paste0("figures/", file_prefix, ".png"), p, width = 10, height = 6, dpi = 600)
  ggsave(paste0("figures/", file_prefix, ".pdf"), p, width = 10, height = 6)
  
  return(p)
}

#---------------------------
# 7. Generate plots
#---------------------------
# Coverage profiles per sample
p_coverage_profiles <- ggplot(plot_df,
                              aes(x = Position, y = Depth + 1,
                                  colour = sample_type, group = sample_type)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  geom_line(aes(y = Depth_baseline + 1), colour = "black", linewidth = 0.3) +
  facet_grid(temp ~ sample_id, scales = "free_y") +
  scale_y_log10() +
  labs(
    x = "Genome Position",
    y = "Sequencing Depth (log10)",
    colour = "Sample type",
    title = "Coverage profiles per sample with Brain-D0 baseline"
  ) +
  scale_colour_manual(values = palette) +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom")

# ggsave("figures/coverage_profiles.png", p_coverage_profiles, width = 14, height = 8, dpi = 600)
# ggsave("figures/coverage_profiles.pdf", p_coverage_profiles, width = 14, height = 8)

# Mean ± SD plots for each threshold
plot_mean_sd(qc_long %>% filter(coverage == "10x"), "mean_value", "% genome ≥10× coverage", "coverage_10x")
plot_mean_sd(qc_long %>% filter(coverage == "50x"), "mean_value", "% genome ≥50× coverage", "coverage_50x")
plot_mean_sd(qc_long %>% filter(coverage == "100x"), "mean_value", "% genome ≥100× coverage", "coverage_100x")

# All coverage thresholds in one faceted plot
p_all_cov <- ggplot(qc_long, aes(x = timepoint, y = mean_value,
                                 colour = sample_type, group = sample_type)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_value - sd_value,
                  ymax = mean_value + sd_value,
                  fill = sample_type),
              alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "black") +
  facet_grid(coverage ~ temp, scales = "free_y") +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Days post sampling",
    y = "Genome coverage (%) (mean ± SD)",
    colour = "Sample type",
    fill = "Sample type",
    title = "Sequencing coverage over time for multiple depth thresholds"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12))

# ggsave("figures/coverage_all_thresholds.png", p_all_cov, width = 12, height = 10, dpi = 600)
# ggsave("figures/coverage_all_thresholds.pdf", p_all_cov, width = 12, height = 10)


#==================================================
# 8. Statistical analysis
#==================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)

# Read sequencing summary data

seq_summary <- read.csv("raw/4aug26_FTA_RDT_study_seq_data - Summary.csv")

# Check columns

names(seq_summary)
qc_summary <- qc_summary %>%
  left_join(
    seq_summary %>%
      select(sample_id, mapped_reads),
    by = "sample_id"
  )
#--------------------------------------------------
# Prepare analysis dataset
#--------------------------------------------------
# Exclude fresh brain controls for preservation comparisons
# Brain is retained separately as the gold-standard baseline

analysis_conditions <- qc_summary %>%
  filter(sample_type != "Brain") %>%
  mutate(
    sample_type = factor(sample_type),
    temp = factor(temp),
    timepoint = factor(timepoint)
  )


# Set biologically meaningful reference groups

analysis_conditions$sample_type <- relevel(
  analysis_conditions$sample_type,
  ref = "FTA"
)

analysis_conditions$temp <- relevel(
  analysis_conditions$temp,
  ref = "RT"
)

analysis_conditions$timepoint <- relevel(
  analysis_conditions$timepoint,
  ref = "1"
)


#==================================================
# MODEL 1
# Sequencing yield: mapped reads
#==================================================

# Response should be log transformed because read counts are skewed
hist(analysis_conditions$mapped_read)
mapped_model <- lmer(
  log10(mapped_reads + 1) ~ sample_type * temp * timepoint +
    (1 | sample_id),
  data = analysis_conditions
)

summary(mapped_model)

anova(mapped_model)


# Pairwise comparisons

emmeans(
  mapped_model,
  pairwise ~ sample_type
)


# Plot mapped reads

p_reads <- ggplot(
  analysis_conditions,
  aes(
    x = timepoint,
    y = log10(mapped_reads + 1),
    colour = sample_type,
    group = sample_type
  )
) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1
  ) +
  facet_wrap(~temp) +
  labs(
    x = "Days post preservation",
    y = "log10(mapped reads + 1)",
    colour = "Sample type",
    title = "Sequencing yield across preservation conditions"
  ) +
  theme_minimal(base_size = 14)

p_reads



#==================================================
# MODEL 2
# Primary genome recovery endpoint
# % genome covered ≥20x
#==================================================

coverage20_model <- lmer(
  pct_20x ~ sample_type * temp * timepoint +
    (1 | sample_id),
  data = analysis_conditions
)

summary(coverage20_model)

anova(coverage20_model)


emmeans(
  coverage20_model,
  pairwise ~ sample_type
)


# Plot

p_20x <- ggplot(
  analysis_conditions,
  aes(
    x = timepoint,
    y = pct_20x,
    colour = sample_type,
    group = sample_type
  )
) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1
  ) +
  facet_wrap(~temp) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed"
  ) +
  labs(
    x = "Days post preservation",
    y = "% genome covered ≥20×",
    colour = "Sample type",
    title = "Genome recovery at 20× depth"
  ) +
  theme_minimal(base_size = 14)

p_20x



#==================================================
# MODEL 3
# High-depth genome recovery
#==================================================

coverage_threshold_models <- list(
  
  "50x" = lmer(
    pct_50x ~ sample_type * temp * timepoint +
      (1 | sample_id),
    data = analysis_conditions
  ),
  
  "100x" = lmer(
    pct_100x ~ sample_type * temp * timepoint +
      (1 | sample_id),
    data = analysis_conditions
  )
)


# Results

lapply(
  coverage_threshold_models,
  anova
)


# Pairwise FTA vs RDT

emmeans(
  coverage_threshold_models$`100x`,
  pairwise ~ sample_type
)


#==================================================
# Plot all depth thresholds
#==================================================

threshold_long <- analysis_conditions %>%
  pivot_longer(
    cols = c(pct_20x, pct_50x, pct_100x),
    names_to = "threshold",
    values_to = "coverage"
  )


p_thresholds <- ggplot(
  threshold_long,
  aes(
    x = timepoint,
    y = coverage,
    colour = sample_type,
    group = sample_type
  )
) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1
  ) +
  facet_grid(threshold ~ temp) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed"
  ) +
  labs(
    x = "Days post preservation",
    y = "% genome recovered",
    colour = "Sample type",
    title = "Genome recovery across sequencing depth thresholds"
  ) +
  theme_minimal(base_size = 14)

p_thresholds



#==================================================
# Brain-normalised recovery
# Secondary analysis
#==================================================

brain_comparison <- qc_summary %>%
  filter(sample_type != "Brain") %>%
  left_join(
    qc_summary %>%
      filter(sample_type == "Brain") %>%
      select(
        sample_id,
        brain_pct20x = pct_20x
      ),
    by="sample_id"
  ) %>%
  mutate(
    relative_recovery = pct_20x / brain_pct20x * 100,
    sample_type = factor(sample_type),
    temp = factor(temp),
    timepoint = factor(timepoint)
  )


brain_model <- lmer(
  relative_recovery ~ sample_type * temp * timepoint +
    (1 | sample_id),
  data = brain_comparison
)

summary(brain_model)

anova(brain_model)