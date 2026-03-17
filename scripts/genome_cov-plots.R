library(tidyverse)

#---------------------------
# 1. Read and parse
#---------------------------
df <- read.csv("raw/FTA_RDT_study_seq_depth_raw.csv")

# Remove specific sample if needed
df_parsed <- df %>%
  filter(Sample != "S1-FTA-RT-D1") %>%
  separate(Sample, into = c("sample_id", "sample_type", "temp", "timepoint"),
           sep = "-", fill = "right")

# Clean timepoints
df_parsed <- df_parsed %>%
  mutate(
    timepoint = as.numeric(gsub("D", "", timepoint)),
    Depth = as.numeric(Depth)
  )

# Fix brain baseline: assume D0 or Baseline samples
df_parsed <- df_parsed %>%
  mutate(
    timepoint = ifelse(sample_type == "Brain", 0, timepoint),
    temp = ifelse(sample_type == "Brain", NA, temp)
  )

#---------------------------
# 2. Plot coverage per sample
#---------------------------
ggplot(df_parsed,
       aes(x = Position,
           y = Depth,
           colour = factor(timepoint),
           group = interaction(sample_id, timepoint))) +
  geom_line(alpha = 0.7) +
  facet_grid(sample_type ~ temp, scales = "free_y") +
  labs(
    x = "Genome Position",
    y = "Read Depth",
    colour = "Days post sampling",
    title = "Coverage profiles per sample with baseline included"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  scale_colour_viridis_d(option = "C")   # nicer palette

#--- 
# summarise the xcoverage
qc_summary <- df_parsed %>%
  group_by(sample_id, sample_type, temp, timepoint) %>%
  summarise(
    genome_length = n(),
    mean_depth = mean(Depth, na.rm = TRUE),
    pct_10x = mean(Depth >= 10) * 100,
    pct_50x = mean(Depth >= 50) * 100,
    pct_100x = mean(Depth >= 100) * 100,
    .groups = "drop"
  )

# extract baseline brain values
baseline <- qc_summary %>%
  filter(sample_type == "Brain") %>%
  select(sample_id,
         baseline_depth = mean_depth,
         baseline_cov = pct_50x)

# calculate change relative to baseline
qc_summary <- qc_summary %>%
  left_join(baseline, by = "sample_id") %>%
  mutate(
    depth_ratio = mean_depth / baseline_depth,
    cov_loss = baseline_cov - pct_50x
  )


# plot degradation 
ggplot(qc_summary %>% filter(sample_type != "Brain"),
       aes(timepoint, pct_10x,
           colour = sample_type,
           group = interaction(sample_id, sample_type))) +
  geom_line(alpha = 0.6) +
  geom_point() +
  facet_wrap(~temp) +
  labs(
    x = "Days post sampling",
    y = "% genome ≥50× coverage",
    colour = "Sample type"
  )

condition_summary <- qc_summary %>%
  filter(sample_type != "Brain") %>%
  group_by(sample_type, temp, timepoint) %>%
  summarise(
    median_depth = median(mean_depth),
    median_cov50 = median(pct_50x),
    n_samples = n(),
    .groups = "drop"
  )

# heatmap 
ggplot(qc_summary %>% filter(sample_type != "Brain"),
       aes(timepoint, sample_id, fill = pct_50x)) +
  geom_tile() +
  facet_grid(sample_type ~ temp) +
  scale_fill_viridis_c(name = "% genome ≥50×") +
  labs(x = "Days", y = "Sample")
