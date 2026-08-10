#==================================================
# Statistical analysis
#==================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)


#==================================================
# Prepare sequencing data
#==================================================

genome_length <- 11923

qc_summary <- read.csv(
  "raw/5aug26_FTA_RDT_study_seq_data - Summary (1).csv"
) %>%
  filter(Experiment_no == 1) %>%
  filter(
    Experiment_no == 1,
    !grepl("NTC|NEC", sample_id)
  ) %>%
  select(-sequences) %>%
  separate(
    sample_id,
    into = c("sample_id", "sample_type", "temp", "timepoint"),
    sep = "-",
    fill = "right"
  ) %>%
  mutate(
    
    # Correct brain controls (e.g. S4-D0)
    sample_type = if_else(
      sample_type == "D0",
      "Brain",
      sample_type
    ),
    
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      TRUE ~ as.numeric(gsub("D", "", timepoint))
    ),
    
    temp = if_else(
      sample_type == "Brain",
      NA_character_,
      temp
    ),
    
    coverage_20 =
      basesCovered_20 / genome_length * 100,
    
    coverage_100 =
      basesCovered_100 / genome_length * 100
  )


#==================================================
# Preservation conditions only
#==================================================

analysis_conditions <- qc_summary %>%
  filter(sample_type != "Brain") %>%
  mutate(
    sample_type = factor(sample_type),
    temp = factor(temp),
    timepoint = factor(timepoint)
  )

# set the reference for intercept (day 1, fta, rt)
analysis_conditions$sample_type <-
  relevel(
    analysis_conditions$sample_type,
    ref = "FTA"
  )

analysis_conditions$temp <-
  relevel(
    analysis_conditions$temp,
    ref = "RT"
  )

analysis_conditions$timepoint <-
  relevel(
    analysis_conditions$timepoint,
    ref = "1"
  )


#==================================================
# DESCRIPTIVE OUTCOME
# Successful genome recovery
# ≥90% genome coverage at ≥20x
#==================================================

analysis_conditions %>%
  summarise(
    successful = sum(coverage_20 >= 90),
    total = n(),
    percent_success = successful / total * 100
  )

# Remove the sample extraction failure S1-FTA-RT-D1
analysis_conditions <- analysis_conditions %>%
  filter(!(sample_id == "S1" &
             sample_type == "FTA" &
             temp == "RT" &
             timepoint == 1))
#==================================================
# Summary table of genome success across conditions
#==================================================

recovery_summary <- analysis_conditions %>%
group_by(sample_type, temp, timepoint) %>%
  summarise(
    n = n(),
    mean_coverage20 = mean(coverage_20),
    sd_coverage20 = sd(coverage_20),
    median_coverage20 = median(coverage_20),
    successful = sum(coverage_20 >= 90),
    success_percent = successful/n*100,
    .groups = "drop"
  )

recovery_summary


ggplot(
  analysis_conditions,
  aes(
    x = timepoint,
    y = coverage_20,
    colour = sample_id,
    shape = sample_type,
    group = interaction(sample_id, sample_type)
  )
) +
  geom_line(
    alpha = 0.5,
    linewidth = 0.8
  ) +
  geom_point(
    aes(size = sample_type),
    fill = NA,
    stroke = 1
  ) +
  facet_wrap(~temp) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +
  scale_shape_manual(
    values = c(
      "FTA" = 21,
      "RDT" = 24
    )
  ) +
  scale_size_manual(
    values = c(
      "FTA" = 3.5,
      "RDT" = 2.8
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "S1" = "#0072B2",
      "S2" = "#D55E00",
      "S3" = "#009E73",
      "S4" = "#CC79A7",
      "S5" = "#E69F00"
    )
  ) +
  scale_y_continuous(
    limits = c(70, 100),
    breaks = seq(70, 100, 5)
  ) +
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Sample",
    shape = "Preservation"
  ) +
  theme_classic(base_size = 10)

ggsave(
  filename = "figures/fig_genome_recovery_coverage20.pdf",
  width = 16.5,
  height = 9,
  units = "cm",
  device = cairo_pdf
)
ggsave(
  filename = "figures/fig_genome_recovery_coverage20.tiff",
  width = 16.5,
  height = 9,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)
ggsave(
  filename = "figures/fig_genome_recovery_coverage20.jpeg",
  width = 16.5,
  height = 9,
  units = "cm",
  dpi = 600
)

#==================================================
# Genome recovery summary
#==================================================

analysis_conditions %>%
  summarise(
    median_coverage20 = median(coverage_20),
    IQR_coverage20 = IQR(coverage_20),
    median_coverage100 = median(coverage_100),
    IQR_coverage100 = IQR(coverage_100),
    successful_genomes = sum(coverage_20 >= 90),
    total = n()
  )

#==================================================
# Model 
# look at loss of genome coverage to assess what factors increase loss of high-depth genome recovery?
#==================================================

# model reduction in genome coverage at ≥100× depth
analysis_conditions <- analysis_conditions %>%
  mutate(
    coverage100_loss = 100 - coverage_100,
    log_coverage100_loss = log10(coverage100_loss)
  )

hist(analysis_conditions$coverage100_loss, breaks=20)
hist(analysis_conditions$log_coverage100_loss, breaks=20)

# save df
write.csv(analysis_conditions, "processed/genome_data.csv", row.names=F)

coverage100_loss_model <- lmer(
  log_coverage100_loss ~
    sample_type * temp * timepoint +
    (1 | sample_id),
  data = analysis_conditions
)

summary(coverage100_loss_model)
anova(coverage100_loss_model)

anova.res <- anova(coverage100_loss_model)
anova_table <- as.data.frame(anova.res)
write.csv(
  anova_table,
  "figures/coverage100_loss_model_anova.csv",
  row.names = TRUE
)


# model diagnostics 
plot(coverage100_loss_model)
qqnorm(residuals(coverage100_loss_model))
qqline(residuals(coverage100_loss_model))

simulationOutput <- simulateResiduals(
  coverage100_loss_model
)
plot(simulationOutput)


library(emmeans)

# Preservation effect
emmeans(
  coverage100_loss_model,
  pairwise ~ sample_type | timepoint,
  adjust = "tukey"
)

# Temperature effect
emmeans(
  coverage100_loss_model,
  pairwise ~ temp | timepoint,
  adjust = "tukey"
)

# Time effect
emmeans(
  coverage100_loss_model,
  pairwise ~ timepoint | sample_type,
  adjust = "tukey"
)

