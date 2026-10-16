library(tidyverse)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(dplyr)
library(tibble)


seq_cq <- read.csv("processed/all-runs_sequencing_stats_prepared_withCq.csv")

## Prepare the cq dat
# Restrict analysis to cohort 1
cq_analysis <- seq_cq  %>%
  filter(cohort == 1) %>%
  mutate(
    sample_type = factor(sample_type),
    temp = factor(temp),
    timepoint = factor(timepoint)
  )

# brain as baseline
baseline <- cq_analysis %>%
  filter(
    sample_type == "Brain",
    !is.na(Cq)
  ) %>%
  select(
    sample_id,
    Baseline_Cq = Cq
  )

# calculate change in Cq
cq_analysis <- cq_analysis %>%
  left_join(
    baseline,
    by = "sample_id"
  ) %>%
  mutate(
    Delta_Cq = Cq - Baseline_Cq
  )

#then remove brain samples
analysis_df <- cq_analysis %>%
filter(sample_type != "Brain")

# model 
model <- lmer(
  Delta_Cq ~ sample_type * temp * timepoint +
    (1 | sample_id),
  data = analysis_df
)

summary(model)
anova(model, type = 3)

# random effect
rand_var=8.048/(8.048+1.512)*100;rand_var
#84.18%

# check residuals

sim_res <- simulateResiduals(
  fittedModel = model,
  n = 1000
)

plot(sim_res)
testUniformity(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

par(mfrow = c(2, 2))
plot(model)

qqnorm(residuals(model))
qqline(residuals(model))
plot(
  fitted(model),
  residuals(model),
  xlab = "Fitted values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)


# emeans
# effect size for sample mediyum
emm_medium <- emmeans(
  model,
  ~ sample_type
)

emm_medium
pairs(
  emm_medium,
  infer = c(TRUE, TRUE)
)
contrast_medium <- contrast(emm_medium, method = "revpairwise")
summary(
  contrast_medium,
  infer = c(TRUE, TRUE)
)

# effect size for sample mediyum
emm_temp <- emmeans(
  model,
  ~ temp
)
emm_temp
pairs(
  emm_temp,
  infer = c(TRUE, TRUE)
)

# effect size for storage duration
emm_timepoint <- emmeans(
  model,
  ~ timepoint
)
emm_temp
pairs(
  emm_timepoint,
  infer = c(TRUE, TRUE)
)


# plot cq versus genome cov 

cohort1 <- seq_cq %>%
  filter(cohort == 1) 
ggplot(
  seq_cq,
  aes(
    x = coverage_20,
    y = Cq ,
    colour = sample_type
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  scale_x_continuous(
    limits = c(65, 100),
    labels = \(x) paste0(x, "%")
  ) +
  #scale_y_reverse() +
  labs(
    x = "Genome coverage ≥20×",
    y = "Cq value"
   # colour = "Preservation medium"
  ) +
  theme_classic()

# ggsave(
#   filename = "figures/fig_Cq_vs_genome_coverage20_sampleType.tiff",
#   width = 16.5,
#   height = 9,
#   units = "cm",
#   dpi = 600,
#   compression = "lzw"
# )
# 
# ggsave(
#   filename = "figures/fig_Cq_vs_genome_coverage20_sampleType.png",
#   width = 16.5,
#   height = 9,
#   units = "cm",
#   dpi = 600
# )


#### Extract results and export all - chatgpt code
library(broom.mixed)
library(writexl)


model_coefficients <- broom.mixed::tidy(
  model,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p.value = signif(p.value, 3)
  ) %>%
  select(
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    p.value
  )

model_coefficients

#anova results 
# Extract Type III ANOVA results
anova_raw <- as.data.frame(anova(model, type = 3))

# Inspect the actual column names
names(anova_raw)

# Convert to a manuscript-ready table
anova_results <- anova_raw %>%
  rownames_to_column("Term") %>%
  rename_with(~ gsub("\\.", " ", .x)) %>%
  rename_with(~ gsub("Pr\\(>F\\)", "p-value", .x)) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

anova_results

# Estimated marginal means
emm_sample_type <- as.data.frame(
  emmeans(model, ~ sample_type)
)

emm_temp <- as.data.frame(
  emmeans(model, ~ temp)
)

emm_timepoint <- as.data.frame(
  emmeans(model, ~ timepoint)
)

# Pairwise comparisons
pairs_sample_type <- as.data.frame(
  pairs(emmeans(model, ~ sample_type), adjust = "tukey")
)

pairs_temp <- as.data.frame(
  pairs(emmeans(model, ~ temp), adjust = "tukey")
)

pairs_timepoint <- as.data.frame(
  pairs(emmeans(model, ~ timepoint), adjust = "tukey")
)


# this collects everyhting in an excel 
write_xlsx(
  list(
    Fixed_effects = model_coefficients,
    Type_III_ANOVA = anova_results,
    EMM_Sample_type = emm_sample_type,
    EMM_Temperature = emm_temp,
    EMM_Timepoint = emm_timepoint,
    Pairs_Sample_type = pairs_sample_type,
    Pairs_Temperature = pairs_temp,
    Pairs_Timepoint = pairs_timepoint
  ),
  path = "tables/Cq_model_results.xlsx"
)
