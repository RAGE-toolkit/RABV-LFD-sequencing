#### Extract results and export all - chatgpt code
library(broom.mixed)
library(writexl)

# collate and export all the model results into an excel- chatgpt code
source("scripts/2_cq-analysis/stats-qc.R")

# ============================================================
# Fixed effects
# ============================================================

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


# ============================================================
# Type III ANOVA
# ============================================================

anova_raw <- as.data.frame(anova(model, type = 3))

anova_results <- anova_raw %>%
  rownames_to_column("Term") %>%
  rename_with(~ gsub("\\.", " ", .x)) %>%
  rename_with(~ gsub("Pr\\(>F\\)", "p-value", .x)) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )


# ============================================================
# Estimated marginal means
# Include 95% confidence intervals
# ============================================================

emm_sample_type <- as.data.frame(
  summary(
    emmeans(model, ~ sample_type),
    infer = c(TRUE, TRUE)
  )
)

emm_temp <- as.data.frame(
  summary(
    emmeans(model, ~ temp),
    infer = c(TRUE, TRUE)
  )
)

emm_timepoint <- as.data.frame(
  summary(
    emmeans(model, ~ timepoint),
    infer = c(TRUE, TRUE)
  )
)


# ============================================================
# Pairwise comparisons
# Include 95% confidence intervals and adjusted P values
# ============================================================

pairs_sample_type <- as.data.frame(
  summary(
    pairs(emmeans(model, ~ sample_type), adjust = "tukey"),
    infer = c(TRUE, TRUE)
  )
)

pairs_temp <- as.data.frame(
  summary(
    pairs(emmeans(model, ~ temp), adjust = "tukey"),
    infer = c(TRUE, TRUE)
  )
)

pairs_timepoint <- as.data.frame(
  summary(
    pairs(emmeans(model, ~ timepoint), adjust = "tukey"),
    infer = c(TRUE, TRUE)
  )
)


# ============================================================
# Random effects and ICC
# ============================================================

random_variance <- as.numeric(
  VarCorr(model)$sample_id[1, 1]
)

residual_variance <- sigma(model)^2

icc <- random_variance / (random_variance + residual_variance)

variance_components <- tibble(
  Component = c(
    "Sample ID random-intercept variance",
    "Residual variance",
    "ICC",
    "ICC (%)"
  ),
  Value = c(
    random_variance,
    residual_variance,
    icc,
    icc * 100
  )
) %>%
  mutate(Value = round(Value, 3))


# ============================================================
# Model diagnostics
# ============================================================

sim_res <- simulateResiduals(
  fittedModel = model,
  n = 1000
)

test_uniformity <- testUniformity(sim_res)
test_dispersion <- testDispersion(sim_res)
test_outliers <- testOutliers(sim_res)

diagnostic_row <- function(test, test_name) {
  tibble(
    Test = test_name,
    Statistic = if (length(test$statistic) > 0) {
      as.numeric(test$statistic[1])
    } else {
      NA_real_
    },
    P_value = if (length(test$p.value) > 0) {
      as.numeric(test$p.value[1])
    } else {
      NA_real_
    },
    Method = if (length(test$method) > 0) {
      paste(test$method, collapse = " ")
    } else {
      NA_character_
    }
  )
}

model_diagnostics <- bind_rows(
  diagnostic_row(test_uniformity, "DHARMa uniformity"),
  diagnostic_row(test_dispersion, "DHARMa dispersion"),
  diagnostic_row(test_outliers, "DHARMa outliers")
) %>%
  mutate(
    Statistic = round(Statistic, 3),
    P_value = signif(P_value, 3)
  )


# ============================================================
# Model summary
# ============================================================

model_summary <- tibble(
  Metric = c(
    "Number of observations",
    "Number of sample IDs",
    "AIC",
    "BIC",
    "Log-likelihood",
    "Residual standard deviation"
  ),
  Value = c(
    nobs(model),
    n_distinct(model.frame(model)$sample_id),
    AIC(model),
    BIC(model),
    as.numeric(logLik(model)),
    sigma(model)
  )
) %>%
  mutate(Value = round(Value, 3))


# ============================================================
# Cq detection summary
# ============================================================

cq_detection <- seq_cq %>%
  summarise(
    Total_specimens = n(),
    RABV_RNA_detected = sum(!is.na(Cq)),
    Not_detected = sum(is.na(Cq)),
    Detection_percentage = mean(!is.na(Cq)) * 100,
    Minimum_Cq = min(Cq, na.rm = TRUE),
    Maximum_Cq = max(Cq, na.rm = TRUE)
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )


# ============================================================
# Export all results to Excel
# ============================================================

write_xlsx(
  list(
    Cq_detection = cq_detection,
    Fixed_effects = model_coefficients,
    Type_III_ANOVA = anova_results,
    EMM_Sample_type = emm_sample_type,
    EMM_Temperature = emm_temp,
    EMM_Timepoint = emm_timepoint,
    Pairs_Sample_type = pairs_sample_type,
    Pairs_Temperature = pairs_temp,
    Pairs_Timepoint = pairs_timepoint,
    Variance_components = variance_components,
    Model_diagnostics = model_diagnostics,
    Model_summary = model_summary
  ),
  path = "tables/Cq_model_results.xlsx"
)
