library(glmmTMB)
library(emmeans)
library(car)
library(dplyr)
library(tibble)
library(writexl)


###Export the model results for gen seq coverage analysis
source("scripts/3_genomeCov-analysis/stats-seq.R")

#### Export results for genome coverage analysis

library(glmmTMB)
library(emmeans)
library(car)
library(dplyr)
library(tibble)
library(writexl)

# Collate and export model results
source("scripts/3_genomeCov-analysis/stats-seq.R")


# ============================================================
# Type III Wald chi-square tests
# ============================================================

anova_res <- car::Anova(model_beta_2way, type = "III")

anova_table <- as.data.frame(anova_res) %>%
  rownames_to_column("Term") %>%
  dplyr::dplyr::rename(
    `Wald chi-square` = Chisq,
    `Degrees of freedom` = Df,
    `p-value` = `Pr(>Chisq)`
  ) %>%
  mutate(
    `Wald chi-square` = round(`Wald chi-square`, 3),
    `p-value` = signif(`p-value`, 3)
  )


# ============================================================
# Fixed effects and 95% confidence intervals
# ============================================================

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
    across(where(is.numeric), ~ round(.x, 3)),
    `p-value` = signif(`p-value`, 3)
  )


# ============================================================
# Likelihood ratio test: full vs reduced model
# ============================================================

model_comparison <- as.data.frame(
  anova(model_beta_2way, model_beta)
) %>%
  rownames_to_column("Model")


# ============================================================
# Estimated marginal means: specimen type by timepoint
# ============================================================

emm_type <- emmeans(
  model_beta_2way,
  ~ sample_type | timepoint
)

emm_type_table <- as.data.frame(
  summary(emm_type, type = "response", infer = c(TRUE, TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ============================================================
# Pairwise specimen-type comparisons at each timepoint
# ============================================================

pair_type_table <- as.data.frame(
  summary(
    pairs(emm_type),
    type = "response",
    infer = c(TRUE, TRUE),
    adjust = "tukey"
  )
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ============================================================
# Estimated marginal means: temperature by timepoint
# ============================================================

emm_temp <- emmeans(
  model_beta_2way,
  ~ temp | timepoint
)

emm_temp_table <- as.data.frame(
  summary(emm_temp, type = "response", infer = c(TRUE, TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ============================================================
# Pairwise temperature comparisons at each timepoint
# ============================================================

pair_temp_table <- as.data.frame(
  summary(
    pairs(emm_temp),
    type = "response",
    infer = c(TRUE, TRUE),
    adjust = "tukey"
  )
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ============================================================
# Timepoint comparisons within each specimen type
# ============================================================

emm_time <- emmeans(
  model_beta_2way,
  ~ timepoint | sample_type
)

emm_time_table <- as.data.frame(
  summary(emm_time, type = "response", infer = c(TRUE, TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

pair_time_table <- as.data.frame(
  summary(
    pairs(emm_time, adjust = "tukey"),
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ============================================================
# Random-effect variance and SD
# ============================================================

vc <- VarCorr(model_beta_2way)$cond$sample_id

variance_table <- tibble(
  Random_effect = "sample_id (random intercept)",
  Variance = as.numeric(vc[1, 1]),
  SD = sqrt(as.numeric(vc[1, 1]))
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )


# ============================================================
# Model diagnostics
# ============================================================

# Use the diagnostics from stats-seq.R if already available.
# Otherwise, calculate them here.

sim_res_beta <- simulateResiduals(model_beta_2way)

test_uniformity <- testUniformity(sim_res_beta)
test_dispersion <- testDispersion(sim_res_beta)
test_outliers <- testOutliers(sim_res_beta)
test_zero_inflation <- testZeroInflation(sim_res_beta)

diagnostic_row <- function(test, test_name) {
  tibble(
    Test = test_name,
    Statistic = if (length(test$statistic) > 0) {
      as.numeric(test$statistic[1])
    } else {
      NA_real_
    },
    `p-value` = if (length(test$p.value) > 0) {
      as.numeric(test$p.value[1])
    } else {
      NA_real_
    }
  )
}

diagnostics_table <- bind_rows(
  diagnostic_row(test_uniformity, "DHARMa uniformity"),
  diagnostic_row(test_dispersion, "DHARMa dispersion"),
  diagnostic_row(test_outliers, "DHARMa outliers"),
  diagnostic_row(test_zero_inflation, "DHARMa zero inflation")
) %>%
  mutate(
    Statistic = round(Statistic, 3),
    `p-value` = signif(`p-value`, 3)
  )


# ============================================================
# Export all results to Excel
# ============================================================

write_xlsx(
  list(
    Type_III_tests = anova_table,
    Fixed_effects = coef_table,
    Model_comparison = model_comparison,
    EMM_sample_type = emm_type_table,
    Pairwise_sample_type = pair_type_table,
    EMM_temperature = emm_temp_table,
    Pairwise_temperature = pair_temp_table,
    EMM_timepoint = emm_time_table,
    Pairwise_timepoint = pair_time_table,
    Random_effects = variance_table,
    Model_diagnostics = diagnostics_table
  ),
  path = "tables/coverage20_beta_model_results.xlsx"
)