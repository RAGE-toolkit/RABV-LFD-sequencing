
library(tidyverse)
## Cq versus genome cov 

cq<- read_csv("processed/Cq_data.csv")

# fix to match
qpcr_analysis <- cq %>%
  rename(
    sample_id = Sample_ID,
    sample_type = Medium,
    temp = Temperature,
    timepoint = Time
  ) %>%
  mutate(
    temp = recode(
      temp,
      "Ambient" = "RT",
      "37°C" = "37"
    ),
    timepoint = as.numeric(gsub("D", "", timepoint))  )

genome<- read_csv("processed/genome_data.csv")


seq_qpcr <- genome %>%
  left_join(
    qpcr_analysis,
    by = c(
      "sample_id",
      "sample_type",
      "temp",
      "timepoint"
    )
  )


# plot of cq by genome cov
ggplot(
  seq_qpcr,
  aes(
    x = Cq,
    y = coverage_20,
    colour = sample_type
  )
) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00"
    )
  ) +
  geom_point(size=3) +
  
  
  geom_hline(
    yintercept = 90,
    
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +
  theme_classic() +
  labs(
    x="Cq",
    y="Genome coverage ≥100× (%)",
    colour="Preservation method"
  )




