#==================================================
# Statistical analysis
#==================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(wesanderson)
library(cowplot)
library(dplyr)
library(ggplot2)

#==================================================
# Prepare sequencing data
#==================================================

qc_all <- read.csv("processed/all-runs_sequencing_stats_prepared.csv")

#==================================================
# Brain only
#==================================================
brain_summary <- qc_all %>%
  filter(sample_type == "Brain") %>%
  mutate(
    brain_source = case_when(
      temp == "Field" ~ "Field brain",
      is.na(temp) & cohort == 1 ~ "Cohort 1 archived brain",
      is.na(temp) & cohort == 2 ~ "Cohort 2 archived brain",
      is.na(temp) ~ "Archived brain",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(brain_source) %>%
  summarise(
    n = n(),
    
    successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
    percent_success_90 = successful_90 / n * 100,
    
    successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
    percent_success_80 = successful_80 / n * 100,
    
    .groups = "drop"
  ) %>%
  bind_rows(
    qc_all %>%
      filter(sample_type == "Brain") %>%
      summarise(
        brain_source = "Overall",
        n = n(),
        
        successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
        percent_success_90 = successful_90 / n * 100,
        
        successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
        percent_success_80 = successful_80 / n * 100
      )
  ) %>%
  mutate(
    brain_source = factor(
      brain_source,
      levels = c(
        "Overall",
        "Cohort 1 archived brain",
        "Cohort 2 archived brain",
        "Field brain"
      )
    )
  ) %>%
  arrange(brain_source)

brain_summary
#write.csv(brain_summary,"processed/ngs_summary_stats_brain.csv", row.names=F)
#==================================================
# Preservation conditions only
#==================================================

analysis_conditions <- qc_all %>%
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
# ≥80% genome coverage at ≥20x
#==================================================
analysis_conditions %>%
  summarise(
    successful = sum(coverage_20 >= 80),
    total = n(),
    percent_success = successful / total * 100
  )

success_summary <- analysis_conditions %>%
  mutate(cohort = as.character(cohort)) %>%
  group_by(cohort) %>%
  summarise(
    n = sum(!is.na(coverage_20)),
    
    successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
    percent_success_90 = successful_90 / n * 100,
    
    successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
    percent_success_80 = successful_80 / n * 100,
    
    .groups = "drop"
  ) %>%
  bind_rows(
    analysis_conditions %>%
      summarise(
        cohort = "Overall",
        n = sum(!is.na(coverage_20)),
        
        successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
        percent_success_90 = successful_90 / n * 100,
        
        successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
        percent_success_80 = successful_80 / n * 100
      )
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c("Overall", "1", "2", "3")
    )
  ) %>%
  arrange(cohort)

success_summary


# Success by cohort and sample type (RDT / FTA)
success_summary <- analysis_conditions %>%
  mutate(
    cohort = as.character(cohort)
  ) %>%
  group_by(cohort, sample_type) %>%
  summarise(
    n = sum(!is.na(coverage_20)),
    
    successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
    percent_success_90 = successful_90 / n * 100,
    
    successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
    percent_success_80 = successful_80 / n * 100,
    
    .groups = "drop"
  ) %>%
  bind_rows(
    # Cohort-specific overall (RDT + FTA)
    analysis_conditions %>%
      mutate(cohort = as.character(cohort)) %>%
      group_by(cohort) %>%
      summarise(
        sample_type = "Overall",
        n = sum(!is.na(coverage_20)),
        
        successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
        percent_success_90 = successful_90 / n * 100,
        
        successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
        percent_success_80 = successful_80 / n * 100,
        
        .groups = "drop"
      )
  ) %>%
  bind_rows(
    # Overall (all cohorts + RDT + FTA)
    analysis_conditions %>%
      summarise(
        cohort = "Overall",
        sample_type = "Overall",
        n = sum(!is.na(coverage_20)),
        
        successful_90 = sum(coverage_20 >= 90, na.rm = TRUE),
        percent_success_90 = successful_90 / n * 100,
        
        successful_80 = sum(coverage_20 >= 80, na.rm = TRUE),
        percent_success_80 = successful_80 / n * 100
      )
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c("Overall", "1", "2", "3")
    ),
    sample_type = factor(
      sample_type,
      levels = c("Overall", "FTA", "RDT")
    )
  ) %>%
  arrange(cohort, sample_type)

success_summary



#write.csv(success_summary,"processed/ngs_summary_stats_sample-mediums.csv", row.names=F)
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

# ggsave(
#   filename = "figures/fig_genome_recovery_coverage20.pdf",
#   width = 16.5,
#   height = 9,
#   units = "cm",
#   device = cairo_pdf
# )
# ggsave(
#   filename = "figures/fig_genome_recovery_coverage20.tiff",
#   width = 16.5,
#   height = 9,
#   units = "cm",
#   dpi = 600,
#   compression = "lzw"
# )
# ggsave(
#   filename = "figures/fig_genome_recovery_coverage20.jpeg",
#   width = 16.5,
#   height = 9,
#   units = "cm",
#   dpi = 600
# )

# Plot medium as the panel 

ggplot(
  analysis_conditions,
  aes(
    x = timepoint,
    y = coverage_20,
    colour = sample_id,
    shape = temp,
    group = interaction(sample_id, temp)
  )
) +
  geom_line(
    alpha = 0.7,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3,
    fill = NA,
    stroke = 1
  ) +
  facet_wrap(~sample_type) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +
  scale_shape_manual(
    values = c(
      "RT" = 21,
      "37" = 24
    )
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
    shape = "Temperature"
  ) +
  theme_classic(base_size = 10)


## Plot 
analysis_conditions <- analysis_conditions %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

# Summary statistics
summary_conditions <- analysis_conditions %>%
  group_by(sample_type, temp, timepoint) %>%
  summarise(
    median_coverage = median(coverage_20, na.rm = TRUE),
    lower = quantile(coverage_20, 0.25, na.rm = TRUE),
    upper = quantile(coverage_20, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  
  # Individual samples
  geom_point(
    data = analysis_conditions,
    aes(
      x = timepoint,
      y = coverage_20,
      colour = sample_type
    ),
    alpha = 0.25,
    size = 1.8,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  # IQR
  geom_ribbon(
    data = summary_conditions,
    aes(
      x = timepoint,
      ymin = lower,
      ymax = upper,
      fill = sample_type,
      group = sample_type
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Median
  geom_line(
    data = summary_conditions,
    aes(
      x = timepoint,
      y = median_coverage,
      colour = sample_type,
      group = sample_type
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = summary_conditions,
    aes(
      x = timepoint,
      y = median_coverage,
      colour = sample_type
    ),
    size = 3
  ) +
  
  # 90% threshold
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +
  
  # Temperature panels
  facet_wrap(~temp) +
  
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
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 3, 7, 21)
  ) +
  
  scale_y_continuous(
    limits = c(70, 100),
    breaks = seq(70, 100, 5)
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Preservation"
  ) +
  
  theme_classic(base_size = 10)

# Summary statistics

summary_conditions <- analysis_conditions %>%
  group_by(sample_type, temp, timepoint) %>%
  summarise(
    median_coverage = median(coverage_20, na.rm = TRUE),
    lower = quantile(coverage_20, 0.25, na.rm = TRUE),
    upper = quantile(coverage_20, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  
  # Individual samples
  geom_point(
    data = analysis_conditions,
    aes(
      x = timepoint,
      y = coverage_20,
      colour = temp
    ),
    alpha = 0.25,
    size = 1.8,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  # IQR
  geom_ribbon(
    data = summary_conditions,
    aes(
      x = timepoint,
      ymin = lower,
      ymax = upper,
      fill = temp,
      group = temp
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Median
  geom_line(
    data = summary_conditions,
    aes(
      x = timepoint,
      y = median_coverage,
      colour = temp,
      group = temp
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = summary_conditions,
    aes(
      x = timepoint,
      y = median_coverage,
      colour = temp
    ),
    size = 3
  ) +
  
  # 90% threshold
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +
  
  # 80% threshold
  geom_hline(
    yintercept = 80,
    linetype = "dotted",
    colour = "grey60",
    linewidth = 0.5
  ) +
  
  # Sample type panels
  facet_wrap(~sample_type) +
  
  scale_colour_manual(
    values = c(
      "RT" = "#0072B2",
      "37" = "#D55E00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "RT" = "#0072B2",
      "37" = "#D55E00"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 3, 7, 21)
  ) +
  
  scale_y_continuous(
    limits = c(60, 100),
    breaks = seq(60, 100, 5)
  ) +
  
  labs(
    x = "Storage time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = "Temperature"
  ) +
  
  theme_classic(base_size = 10)


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
# What factors increase loss of high-depth (>=20x) genome coverage?
#==================================================
#==================================================
# 1. Prepare data
#==================================================

# Safety filter: analysis is FTA vs RDT only (baseline/brain-tissue excluded)
analysis_conditions <- analysis_conditions %>%
  filter(sample_type %in% c("FTA", "RDT")) %>%
  filter(cohort %in% "1")

# Convert coverage_20 (0-100 scale) to a proportion
n_obs <- nrow(analysis_conditions)
analysis_conditions <- analysis_conditions %>%
  mutate(
    coverage_prop     = coverage_20 / 100,
  )

hist(analysis_conditions$coverage_prop, breaks = 20,
     main = "Adjusted coverage proportion (>=20x)", xlab = "Proportion")

#==================================================
# 2. Primary model: beta regression mixed model
#==================================================

analysis_conditions$timepoint <- factor(analysis_conditions$timepoint)

model_beta <- glmmTMB(
  coverage_prop ~ sample_type * temp * timepoint + (1 | sample_id),
  family = beta_family(),
  data   = analysis_conditions
)
summary(model_beta)

# Test whether the 3-way interaction is actually justified before
# interpreting it (likelihood ratio test, reduced vs full model)
model_beta_2way <- update(model_beta, . ~ . - sample_type:temp:timepoint)
anova(model_beta_2way, model_beta)


summary(model_beta_2way)
confint(model_beta_2way)
#Type III Wald chi-square
car::Anova(model_beta_2way, type = "III")
#==================================================
# 3. Model diagnostics (DHARMa)
#==================================================
sim_res_beta <- simulateResiduals(model_beta_2way)
plot(sim_res_beta)
testUniformity(sim_res_beta)
testDispersion(sim_res_beta)
testOutliers(sim_res_beta)
testZeroInflation(sim_res_beta)
 # all ok

#==================================================
# 4. Post-hoc pairwise comparisons
#==================================================

# Preservation medium (FTA vs RDT) effect, at each timepoint
emmeans(model_beta_2way, pairwise ~ sample_type | timepoint, adjust = "tukey")

# Temperature effect, at each timepoint
emmeans(model_beta_2way, pairwise ~ temp | timepoint, adjust = "tukey")

# Time effect, within each medium
emmeans(model_beta_2way, pairwise ~ timepoint | sample_type, adjust = "tukey")

#fta cov
plogis(2.03)
#rdt cov 
plogis(1.80)

#
exp(0.766)
exp(0.649)
exp(0.762)


# Specimen type at each timepoint: OR FTA:RDT, plus back-transformed means
emm_type <- emmeans(model_beta_2way, ~ sample_type | timepoint)
summary(emm_type, type = "response")                # mean coverage (proportion) with 95% CI
summary(pairs(emm_type), type = "response", infer = TRUE)   # FTA - RDT , OR FTA:RDT, 95% CI

# Temperature at each timepoint
emm_temp <- emmeans(model_beta_2way, ~ temp | timepoint)
summary(pairs(emm_temp), type = "response", infer = TRUE)                    # RT - 37 , OR RT:37
summary(pairs(emm_temp, reverse = TRUE), type = "response", infer = TRUE)    # 37 - RT , OR 37:RT
#==================================================
# 5. Save processed data
#==================================================
# write.csv(analysis_conditions, "processed/genome_data.csv", row.names = FALSE)

#==================================================
# 6. Save ANOVA-style output for the interaction test, for reporting
#==================================================
# anova_res <- anova(model_beta_2way, model_beta)
# write.csv(as.data.frame(anova_res), "tables/coverage20_beta_model_anova.csv",row.names = TRUE)


#### PLOTS to showcase coverage 


#heatmap 
# prepare data

qc_all %>%
  filter(sample_type == "Brain") %>%
  select(
    cohort,
    run_id,
    sample_id,
    sample_type,
    temp,
    timepoint,
    coverage_20,
    coverage_100,
    filepath
  ) %>%
  arrange(cohort, sample_id)

brain_d0 <- qc_all %>%
  filter(
    sample_type == "Brain",
    cohort == 1
  ) %>%
  mutate(
    sample_type = "Brain",
    timepoint = 0
  ) %>%
  select(
    sample_id,
    sample_type,
    temp,
    timepoint,
    coverage_20
  )

preservation_data <- analysis_conditions %>%
  filter(cohort == 1) %>%
  mutate(
    sample_id = as.character(sample_id),
    timepoint = as.numeric(as.character(timepoint)),
    temp = as.character(temp),
    sample_type = as.character(sample_type)
  )

preservation_plot <- preservation_data %>%
  mutate(
    x = case_when(
      temp == "RT" & timepoint == 1  ~ 2,
      temp == "RT" & timepoint == 3  ~ 3,
      temp == "RT" & timepoint == 7  ~ 4,
      temp == "RT" & timepoint == 21 ~ 5,
      
      temp == "37" & timepoint == 1  ~ 7,
      temp == "37" & timepoint == 3  ~ 8,
      temp == "37" & timepoint == 7  ~ 9,
      temp == "37" & timepoint == 21 ~ 10
    )
  )

brain_plot <- brain_d0 %>%
  mutate(
    x = 1
  )

heatmap_data <- bind_rows(
  brain_plot,
  preservation_plot
)

ggplot(
  heatmap_data,
  aes(
    x = x,
    y = sample_id,
    fill = coverage_20
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 1.5,
    colour = "black",
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = 6,
    colour = "black",
    linewidth = 0.6
  ) +
  
  facet_grid(
    sample_type ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 7, 8, 9, 10),
    labels = c(
      "D0",
      "D1", "D3", "D7", "D21",
      "D1", "D3", "D7", "D21"
    ),
    expand = expansion(add = 0.05)
  ) +
  
  scale_fill_viridis_c(
    limits = c(0, 100),
    name = "Genome coverage\n≥20× (%)"
  ) +
  
  labs(
    x = NULL,
    y = "Sample"
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

# Cohort 1 brain samples = common D0 baseline
brain_d0 <- qc_all %>%
  filter(
    sample_type == "Brain",
    cohort == 1
  ) %>%
  select(
    sample_id,
    coverage_20
  ) %>%
  mutate(
    timepoint = 0,
    temp = "D0"
  )

# Duplicate the same D0 baseline into RDT and FTA facets
brain_d0_plot <- bind_rows(
  brain_d0 %>% mutate(sample_type = "RDT"),
  brain_d0 %>% mutate(sample_type = "FTA")
)

# Preservation samples
preservation_plot <- qc_all %>%
  filter(
    cohort == 1,
    sample_type %in% c("FTA", "RDT")
  ) %>%
  mutate(
    timepoint = as.numeric(as.character(timepoint)),
    temp = as.character(temp)
  ) %>%
  select(
    sample_id,
    sample_type,
    timepoint,
    temp,
    coverage_20
  )

# Combine
heatmap_data <- bind_rows(
  brain_d0_plot,
  preservation_plot
)

heatmap_data <- heatmap_data %>%
  mutate(
    x = case_when(
      temp == "D0" ~ 1,
      
      temp == "RT" & timepoint == 1  ~ 2,
      temp == "RT" & timepoint == 3  ~ 3,
      temp == "RT" & timepoint == 7  ~ 4,
      temp == "RT" & timepoint == 21 ~ 5,
      
      temp == "37" & timepoint == 1  ~ 7,
      temp == "37" & timepoint == 3  ~ 8,
      temp == "37" & timepoint == 7  ~ 9,
      temp == "37" & timepoint == 21 ~ 10
    ),
    
    sample_id = factor(
      sample_id,
      levels = rev(paste0("S", 1:5))
    ),
    
    sample_type = factor(
      sample_type,
      levels = c("RDT", "FTA")
    )
  )



pal <- wes_palette(
  "Zissou1",
  100,
  type = "continuous"
)
heatmap_data <- heatmap_data %>%
  mutate(
    temperature_group = case_when(
      x == 1 ~ "",
      x %in% 2:5 ~ "RT",
      x %in% 7:10 ~ "37°C"
    )
  )

temp_labels <- tibble(
  sample_type = "RDT",
  x = c(3.5, 8.5),
  label = c("RT", "37°C")
)

heatmap <- ggplot(
  heatmap_data,
  aes(
    x = x,
    y = sample_id,
    fill = coverage_20
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  # geom_text(
  #   data = temp_labels,
  #   aes(
  #     x = x,
  #     y = Inf,
  #     label = label
  #   ),
  #   inherit.aes = FALSE,
  #   vjust = -0.8,
  #   fontface = "bold",
  #   size = 3.5
  # ) +
  
  # D0 | RT boundary
  geom_vline(
    xintercept = 1.5,
    colour = "black",
    linewidth = 0.6
  ) +
  
  # RT | 37°C boundary
  geom_vline(
    xintercept = 6,
    colour = "black",
    linewidth = 0.8
  ) +
  facet_grid(sample_type ~ .)+
  
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 7, 8, 9, 10),
    labels = c(
      "Brain",
      "D1", "D3", "D7", "D21",
      "D1", "D3", "D7", "D21"
    ),
    expand = expansion(add = 0.05)
  )+
  scale_fill_gradientn(
    colours = rev(pal),
    limits = c(0, 100),
    breaks = c(0, 20, 40, 60, 80, 90, 100),
    name = "Genome coverage\n≥20× (%)"
  )+
  labs(
    x = NULL,
    y = "Sample"
  ) +
  
  theme_classic(base_size = 13) +
  coord_cartesian(
    clip = "off"
  ) +
  theme(
    plot.margin = margin(
      t = 15,
      r = 5,
      b = 5,
      l = 5
    ),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(
      angle = 0,
      face = "bold",
      size = 13
    )
  )+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
temp_plot <- ggplot() +
  annotate(
    "text",
    x = c(3, 7),
    y = 1,
    label = c("RT", "37°C"),
    fontface = "bold",
    size = 4.5
  ) +
  scale_x_continuous(
    limits = c(0.5, 10.5),
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

p <- plot_grid(
  temp_plot,
  heatmap,
  ncol = 1,
  rel_heights = c(0.05, 1)
);p

# ggsave(
#   "figures/genome_coverage_heatmap.tiff",
#   plot = p,
#   width = 180,
#   height = 120,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )
# 
# ggsave(
#   "figures/genome_coverage_heatmap.png",
#   plot = p,
#   width = 180,
#   height = 120,
#   units = "mm",
#   dpi = 600
# )


# trajectory plots- lines over time
# Cohort 1 brain D0 values
brain_d0 <- qc_all %>%
  filter(
    sample_type == "Brain",
    cohort == 1,
    sample_id %in% paste0("S", 1:5)
  ) %>%
  select(
    sample_id,
    coverage_20
  )

# Preservation samples
preservation <- qc_all %>%
  filter(
    cohort == 1,
    sample_type %in% c("RDT", "FTA")
  ) %>%
  mutate(
    timepoint = as.numeric(as.character(timepoint)),
    temp = case_when(
      temp == "37" ~ "37°C",
      TRUE ~ as.character(temp)
    )
  ) %>%
  select(
    sample_id,
    sample_type,
    timepoint,
    temp,
    coverage_20
  )

# Duplicate the common brain D0 baseline into BOTH temperature panels
brain_d0 <- bind_rows(
  brain_d0 %>%
    mutate(
      sample_type = "RDT",
      timepoint = 0,
      temp = "RT"
    ),
  
  brain_d0 %>%
    mutate(
      sample_type = "RDT",
      timepoint = 0,
      temp = "37°C"
    ),
  
  brain_d0 %>%
    mutate(
      sample_type = "FTA",
      timepoint = 0,
      temp = "RT"
    ),
  
  brain_d0 %>%
    mutate(
      sample_type = "FTA",
      timepoint = 0,
      temp = "37°C"
    )
)

# Combine
trajectory_data <- bind_rows(
  brain_d0,
  preservation
) %>%
  mutate(
    temp = factor(
      temp,
      levels = c("RT", "37°C")
    ),
    
    point_type = case_when(
      timepoint == 0 ~ "Brain",
      temp == "RT" ~ "RT",
      temp == "37°C" ~ "37°C"
    )
  )


panel_labels <- data.frame(
  sample_type = c("RDT", "RDT", "FTA", "FTA"),
  temp = c("RT", "37°C", "RT", "37°C"),
  label = c("d)", "c)", "b)", "a)"),
  x = -0.5,
  y = 105
)

p <- ggplot(
  trajectory_data,
  aes(
    x = timepoint,
    y = coverage_20,
    group = sample_id
  )
) +
  
  # Individual trajectories
  geom_line(
    colour = "grey60",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.5
  ) +
  
  # D0 Brain — OPEN CIRCLE
  geom_point(
    data = trajectory_data %>%
      filter(point_type == "Brain"),
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8,
    size = 2
  ) +
  geom_point(
    data = trajectory_data %>% filter(point_type != "Brain"),
    aes(colour = sample_type),
    shape = 16,
    size = 2.8
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  )+
  geom_text(
    data = panel_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 12 / ggplot2:::.pt,
    fontface = "bold"
  )+
  # # RT — coloured circle
  # geom_point(
  #   data = trajectory_data %>%
  #     filter(point_type == "RT"),
  #   shape = 16,
  #   colour = rt_col,
  #   size = 2.8
  # ) +
  # 
  # # 37°C — coloured circle
  # geom_point(
  #   data = trajectory_data %>%
  #     filter(point_type == "37°C"),
  #   shape = 16,
  #   colour = temp37_col,
  #   size = 2.8
  # ) +
  
  facet_grid(
    sample_type ~ temp
  ) +
  
  scale_x_continuous(
    breaks = c(0, 1, 3, 7, 21),
    labels = c("0", "1", "3", "7", "21"),
    limits = c(-0.5, 21)
  ) +
  
  scale_y_continuous(
    limits = c(60, 105),
    breaks = seq(60, 100, 20),
    expand = expansion(mult = c(0, 0.01))
  ) +
  
  labs(
    x = NULL,
    y = "Genome coverage ≥20× (%)"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = 12
    ),
    
    strip.text.y = element_text(
      angle = 0,
      face = "bold",
      size = 12
    ),
    
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    panel.spacing = unit(1, "lines")
  )+
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 12
    ),
    strip.text.y = element_text(
      angle = 0,
      face = "bold",
      size = 12
    )
  )

p

# ggsave(
#   "figures/genome_coverage_trajectories_cohort1.tiff",
#   plot = p,
#   width = 180,
#   height = 120,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )
# 
# ggsave(
#   "figures/genome_coverage_trajectories_cohort1.png",
#   plot = p,
#   width = 180,
#   height = 120,
#   units = "mm",
#   dpi = 600
# )


# cohrot 2 and 3 
library(ggplot2)
library(dplyr)

# Cohorts 2 and 3
trajectory_validation <- qc_all %>%
  filter(
    cohort %in% c(2, 3),
    sample_type %in% c("Brain", "FTA", "RDT")
  ) %>%
  mutate(
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      TRUE ~ as.numeric(as.character(timepoint))
    ),
    sample_group = str_remove(sample_id, "-RDT$"),
    sample_type = factor(
      sample_type,
      levels = c("Brain", "FTA", "RDT")
    ),
    cohort = factor(cohort)
  )

p_validation <- ggplot(
  trajectory_validation,
  aes(
    x = timepoint,
    y = coverage_20,
    group = sample_group,
    colour = sample_type
  )
) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.5
  ) +
  geom_line(
    colour = "grey65",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  geom_point(
    data = trajectory_validation %>% filter(sample_type == "Brain"),
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8,
    size = 2
  ) +
  geom_point(
    data = trajectory_validation %>% filter(sample_type != "Brain"),
    shape = 16,
    size = 2.8
  ) +
  facet_wrap(
    ~ cohort,
    scales = "free_x"
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 3, 7, 14, 21, 28, 35, 42, 50, 56, 70, 77),
    labels = c(
      "0", "1", "3", "7", "14", "21",
      "28", "35", "42", "50", "56", "70", "77"
    ),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    limits = c(50, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = c(
      "Brain" = "black",
      "FTA" = fta_col,
      "RDT" = rdt_col
    )
  ) +
  labs(
    x = "Time (days)",
    y = "Genome coverage ≥20× (%)",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1
      )
    )

p_validation
# ggsave(
#   "figures/genome_coverage_validation.png",
#   plot = p_validation,
#   width = 180,
#   height = 100,
#   units = "mm",
#   dpi = 600
# )
# 
# ggsave(
#   "figures/genome_coverage_validation.tiff",
#   plot = p_validation,
#   width = 180,
#   height = 100,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )
