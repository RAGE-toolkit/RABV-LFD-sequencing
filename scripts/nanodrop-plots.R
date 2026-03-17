library(tidyverse)
library(zoo)

#---------------------------
# 1. READ CSV WITH MULTI-HEADER
#---------------------------
raw <- read_csv("raw/nanodrop_result_summary.csv", col_names = FALSE)

# Combine multi-row headers
header1 <- zoo::na.locf(as.character(unlist(raw[1, ])), na.rm = FALSE)
header2 <- as.character(unlist(raw[2, ]))
colnames(raw) <- ifelse(is.na(header2) | header2 == "", header1, paste(header1, header2, sep = "_"))

# Remove header rows
df <- raw[-c(1,2), ]
colnames(df)[1] <- "Sample"

#---------------------------
# 2. CLEAN COLUMN NAMES
#---------------------------
df <- df %>%
  rename_with(~str_replace_all(., " ", "")) %>%
  rename_with(~str_replace_all(., "\\(ng/uL\\)", "conc"))

#---------------------------
# 3. PIVOT LONG AND CLEAN
#---------------------------
df_long <- df %>%
  pivot_longer(-Sample, names_to = "metric_time", values_to = "value") %>%
  drop_na(value) %>%
  mutate(
    value = as.numeric(str_replace_all(str_trim(str_replace_all(value, ",", "")), "")),
    Sample = str_trim(Sample)
  ) %>%
  separate(metric_time, into = c("metric", "timepoint"), sep = "_") %>%
  mutate(
    timepoint = as.numeric(str_remove(timepoint, "D")),
    sample_type = case_when(
      str_detect(Sample, "^FTA") ~ "FTA",
      str_detect(Sample, "^RDT") ~ "RDT",
      TRUE ~ "Baseline"
    ),
    temp = case_when(
      str_detect(Sample, "RT") ~ "RT",
      str_detect(Sample, "37") ~ "37C",
      TRUE ~ NA_character_
    ),
    sample_id = str_extract(Sample, "S\\d+")
  ) %>%
  filter(Sample != "NEC")  # remove negative controls

#---------------------------
# 4. PREPARE DATA FOR PLOTTING
#---------------------------
# Baseline for all temps
baseline_df <- df_long %>%
  filter(sample_type == "Baseline", str_detect(metric, "Concentration")) %>%
  mutate(timepoint = 0)

baseline_expanded <- bind_rows(
  baseline_df %>% mutate(temp = "RT"),
  baseline_df %>% mutate(temp = "37C")
)

# Non-baseline
non_baseline_df <- df_long %>%
  filter(sample_type != "Baseline", str_detect(metric, "Concentration"))

# Combine
plot_df <- bind_rows(non_baseline_df, baseline_expanded)

# Summary for mean ± SD
summary_df <- plot_df %>%
  group_by(sample_type, temp, timepoint) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            .groups = "drop")

#---------------------------
# 5. PER-SAMPLE PLOT
#---------------------------
per_sample_plot <- ggplot(plot_df,
                          aes(x = timepoint,
                              y = value,
                              colour = sample_type,
                              linetype = temp,
                              group = interaction(sample_type, temp))) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~sample_id, ncol = 2, scales = "free_y") +
  labs(
    x = "Days post sampling",
    y = "RNA concentration (ng/µL)",
    colour = "Sample type",
    linetype = "Temperature",
    title = "RNA degradation over time per sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_colour_manual(values = c(
    "Baseline" = "black",
    "FTA" = "#1b9e77",
    "RDT" = "#d95f02"
  ))

# Save per-sample plot
ggsave("figures/RNA_per_sample.pdf", per_sample_plot, width = 10, height = 8, dpi = 300)
ggsave("figures/RNA_per_sample.png", per_sample_plot, width = 10, height = 8, dpi = 600)

#---------------------------
# 6. MEAN ± SD SUMMARY PLOT
#---------------------------
mean_sd_plot <- ggplot(summary_df,
                       aes(x = timepoint, y = mean,
                           colour = sample_type,
                           group = sample_type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.3, alpha = 0.8, size = 0.8) +
  facet_wrap(~temp) +
  labs(
    x = "Time (days)",
    y = "RNA concentration (ng/µL)",
    colour = "Sample type",
    title = "RNA degradation over time (mean ± SD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_colour_manual(values = c(
    "Baseline" = "black",
    "FTA" = "#1b9e77",
    "RDT" = "#d95f02"
  ))

# Save summary plot
ggsave("figures/RNA_mean_sd.pdf", mean_sd_plot, width = 8, height = 6, dpi = 300)
ggsave("figures/RNA_mean_sd.png", mean_sd_plot, width = 8, height = 6, dpi = 600)