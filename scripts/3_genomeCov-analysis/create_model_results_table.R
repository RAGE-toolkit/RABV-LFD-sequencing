library(glmmTMB)
library(emmeans)
library(car)
library(dplyr)
library(tibble)
library(writexl)


###Export the model results for gen seq coverage analysis
source("scripts/3_genomeCov-analysis/stats-seq.R")

# Type III Wald chi-square tests
anova_res <- car::Anova(model_beta_2way, type = "III")

anova_table <- as.data.frame(anova_res) %>%
  rownames_to_column("Term") %>%
  dplyr::rename(
    `Wald chi-square` = Chisq,
    `Degrees of freedom` = Df,
    `p-value` = `Pr(>Chisq)`
  ) %>%
  mutate(
    `Wald chi-square` = round(`Wald chi-square`, 3),
    `p-value` = signif(`p-value`, 3)
  )

anova_table

# Fixed effects and confidence intervals
coef_table <- as.data.frame(
  summary(model_beta_2way)$coefficients$cond
) %>%
  rownames_to_column("Term") %>%
  dplyr::rename(
    Estimate = Estimate,
    SE = `Std. Error`,
    `z-statistic` = `z value`,
    `p-value` = `Pr(>|z|)`
  )

# Add Wald 95% confidence intervals
coef_ci <- as.data.frame(
  confint(model_beta_2way, parm = "beta_", method = "Wald")
) %>%
  rownames_to_column("Term") %>%
  dplyr::rename(
    `95% CI lower` = `2.5 %`,
    `95% CI upper` = `97.5 %`
  )

coef_table <- coef_table %>%
  left_join(coef_ci, by = "Term") %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    ),
    `p-value` = signif(`p-value`, 3)
  )

coef_table

# testing dropping 3 way itneraction
model_comparison <- as.data.frame(
  anova(model_beta_2way, model_beta)
) %>%
  rownames_to_column("Model")

model_comparison

# EMM
emm_type <- emmeans(
  model_beta_2way,
  ~ sample_type | timepoint
)

emm_type_table <- as.data.frame(
  summary(emm_type, type = "response", infer = c(TRUE, TRUE))
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

emm_type_table

# coverage by temp and timpoiint
emm_temp <- emmeans(
  model_beta_2way,
  ~ temp | timepoint
)

emm_temp_table <- as.data.frame(
  summary(emm_temp, type = "response", infer = c(TRUE, TRUE))
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

emm_temp_table

# pairwise comparisons 
pair_type_table <- as.data.frame(
  summary(
    pairs(emm_type),
    type = "response",
    infer = c(TRUE, TRUE),
    adjust = "tukey"
  )
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

pair_type_table

# tm comparisons
pair_temp_table <- as.data.frame(
  summary(
    pairs(emm_temp),
    type = "response",
    infer = c(TRUE, TRUE),
    adjust = "tukey"
  )
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

pair_temp_table

# random effect variance 
# Extract random-effect variance
variance_table <- as.data.frame(
  VarCorr(model_beta_2way)$cond$sample_id
)

variance_table


# save all to one excel 
write_xlsx(
  list(
    Type_III_tests = anova_table,
    Fixed_effects = coef_table,
    Model_comparison = model_comparison,
    EMM_sample_type = emm_type_table,
    EMM_temperature = emm_temp_table,
    Pairwise_sample_type = pair_type_table,
    Pairwise_temperature = pair_temp_table,
    Random_effects = variance_table
  ),
  path = "tables/coverage20_beta_model_results.xlsx"
)
