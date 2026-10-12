library(tidyverse)
library(lmerTest)
library(emmeans)
library(DHARMa)
# ==================================================
# Prepare Cq data- combine with sequencing stats data
# ==================================================
cq_raw <- read_csv(
  "raw/qPCR_lyssavirus_Cq.csv"
)
cq_raw <-cq_raw %>%
  mutate(
    # Correct sample naming
    sample_id= if_else(
      Sample_ID == "N14",
      "N14-D0",
      Sample_ID
    ))%>%
  mutate(
    # Some lab IDs use underscores instead of dashes
    # Normalize so they parse the same way as the rest
    sample_id = gsub("_", "-", sample_id),
    # Remove dash between letter and number at start of ID
    # e.g. N-10-D0 -> N10-D0
    sample_id = gsub("^([A-Za-z])-([0-9])", "\\1\\2", sample_id)
  )

cq_all  <- cq_raw  %>%
  separate(
    sample_id,
    into = c("sample_id", "sample_type", "temp", "timepoint"),
    sep = "-",
    fill = "right"
  ) %>%
  mutate(
    
    # Brain samples: N6-D0
    sample_type = case_when(
      sample_type == "D0" ~ "Brain",
      TRUE ~ sample_type
    ),
    
    # Experiment 2: N6-FTA-D3
    timepoint_from_temp = as.numeric(
      ifelse(
        grepl("^D", temp),
        sub("D", "", temp),
        NA_character_
      )
    ),
    
    timepoint_from_tp = as.numeric(
      ifelse(
        grepl("^D", timepoint),
        sub("D", "", timepoint),
        NA_character_
      )
    ),
    
    timepoint = case_when(
      sample_type == "Brain" ~ 0,
      is.na(timepoint) & grepl("^D", temp) ~ timepoint_from_temp,
      grepl("^D", timepoint) ~ timepoint_from_tp,
      TRUE ~ as.numeric(
        ifelse(
          grepl("^D", timepoint),
          NA_character_,
          timepoint
        )
      )
    ),
    
    # Move D3 out of temp column for Exp 2
    temp = case_when(
      sample_type == "Brain" ~ NA_character_,
      grepl("^D", temp) ~ "RT",
      TRUE ~ temp
    )
  ) %>%
  select(-timepoint_from_temp, -timepoint_from_tp) %>%
  
  # Drop bare IDs with no method/day info to parse
  filter(!is.na(sample_type))

## import and correct the field sample ids
field <- read_csv("raw/cohort3_epidata.csv") 

cq_all <- cq_all %>%
  left_join(
    field %>%
      select(Case_id, Lab_id),
    by = c("sample_id" = "Case_id")
  ) %>%
  mutate(
    sample_id = coalesce(Lab_id, sample_id)
  ) %>%
  select(-Lab_id)%>%
  mutate(
    # Field samples: RADDL4B-2026-### (optionally _RDT suffixed)
    # and H##-###-# style IDs
    is_field = grepl("^RADDL|^H[0-9]", sample_id)
  )


# want to add to sequencing data
seq <- read.csv("processed/all-runs_sequencing_stats_prepared.csv")

seq_cq <- seq %>%
  left_join(
    cq_all %>%
      select(
        sample_id,
        sample_type,
        temp,
        timepoint,
        everything()
      ),
    by = c(
      "sample_id",
      "sample_type",
      "temp",
      "timepoint"
    )
  )

write.csv(seq_cq, "processed/all-runs_sequencing_stats_prepared_withCq.csv", row.names=F)


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

ggsave(
  filename = "figures/fig_Cq_vs_genome_coverage20_sampleType.tiff",
  width = 16.5,
  height = 9,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "figures/fig_Cq_vs_genome_coverage20_sampleType.png",
  width = 16.5,
  height = 9,
  units = "cm",
  dpi = 600
)
