#--------------------------------------------------
# 1. Load packages
#--------------------------------------------------

library(tidyverse)
library(lmerTest)
library(emmeans)
library(afex)
library(stargazer)

# referemces
#https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html
#https://ourcodingclub.github.io/tutorials/mixed-models/
#https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html
#--------------------------------------------------
# 2. Import qPCR data with two-row header
#--------------------------------------------------

ct_raw <- read_csv(
  "raw/4aug26_FTA_RDT_study_qPCR_data - Cq_lyssavirus (1).csv",
  col_names = FALSE
)

# Create column names from second header row
new_names <- c(
  "Sample_Name",
  as.character(ct_raw[2, -1])
)

# Remove header rows
ct <- ct_raw[-c(1, 2), ]
colnames(ct) <- new_names


#--------------------------------------------------
# 3. Convert to long format and assign metadata
# (single row per observation)
#--------------------------------------------------

ct_long <- ct %>%
  pivot_longer(
    cols = starts_with("D"),
    names_to = "Time",
    values_to = "Cq"
  ) %>%
  mutate(
    Cq = as.numeric(Cq),
    
    # Link storage samples back to original brain sample
    Sample_ID = case_when(
      str_detect(Sample_Name, "^S[0-9]+$") ~ Sample_Name,
      TRUE ~ str_extract(Sample_Name, "S[0-9]+$")
    ),
    
    # Preservation medium
    Medium = case_when(
      str_detect(Sample_Name, "FTA") ~ "FTA",
      str_detect(Sample_Name, "RDT") ~ "RDT",
      TRUE ~ "Brain tissue"
    ),
    
    # Storage temperature
    Temperature = case_when(
      str_detect(Sample_Name, "37") ~ "37°C",
      str_detect(Sample_Name, "RT") ~ "Ambient",
      TRUE ~ "Baseline"
    )
  )


#--------------------------------------------------
# 4. Calculate ΔCq relative to baseline brain tissue
#--------------------------------------------------

baseline <- ct_long %>%
  filter(
    Medium == "Brain tissue",
    !is.na(Cq)
  ) %>%
  select(
    Sample_ID,
    Baseline_Cq = Cq
  )


ct_long <- ct_long %>%
  left_join(
    baseline,
    by = "Sample_ID"
  ) %>%
  mutate(
    Delta_Cq = Cq - Baseline_Cq
  )


#--------------------------------------------------
# 5. Remove baseline samples for statistical analysis
#--------------------------------------------------

analysis_df <- ct_long %>%
  filter(
    Medium != "Brain tissue"
  )


#--------------------------------------------------
# 6. Explore raw ΔCq patterns
#--------------------------------------------------
hist(analysis_df$Delta_Cq, breaks=20) 

analysis_df %>%
  group_by(Temperature) %>%
  summarise(
    mean_DeltaCq = mean(Delta_Cq, na.rm = TRUE),
    sd = sd(Delta_Cq, na.rm = TRUE),
    n = n()
  )


analysis_df %>%
  group_by(Medium, Temperature, Time) %>%
  summarise(
    mean_DeltaCq = mean(Delta_Cq, na.rm = TRUE),
    sd = sd(Delta_Cq, na.rm = TRUE),
    n = n()
  )

# save df
# write.csv(analysis_df, "processed/Cq_data.csv", row.names=F)
#--------------------------------------------------
# 7. Fit linear mixed-effects model
#--------------------------------------------------
# note this uses lmerTest lmer, which produces p-values
model <- lmer(
  Delta_Cq ~ Medium * Temperature * Time +
    (1 | Sample_ID),
  data = analysis_df
)

# Full model summary
summary(model)
anova(model)
# sample explains ~84% of variance

#check residuals
par(mfrow=c(2,2))
plot(model)
qqnorm(residuals(model))
qqline(residuals(model))

# dharma diagnostics

#bit bimodal for medium (more well preserved for FTA so less of a change)
ggplot(analysis_df, aes(Delta_Cq, fill = Medium)) +
  geom_histogram(alpha=0.5)

#table of results with sjjplot
tab_model(
  model,
  show.re.var = TRUE,
  show.icc = TRUE,
  digits = 3
)
#--------------------------------------------------
# 8. Test effects
#--------------------------------------------------

anova(
  model,
  type = 3
)

#--------------------------------------------------
# 9. Extract estimated effects
# estimated marginal means (least squares)
#--------------------------------------------------

# Overall medium effect
medium.emm.s <- emmeans(
  model,
  pairwise ~ Medium
); medium.emm.s

# Medium effect within each storage time
medium.time.emm.s <- emmeans(
  model,
  pairwise ~ Medium | Time
);medium.time.emm.s 


# Temperature effect
temp.emm.s <- emmeans(
  model,
  pairwise ~ Temperature
); temp.emm.s


# Time effect
time.emm.s <- emmeans(
  model,
  pairwise ~ Time
);time.emm.s 

