
library(tidyverse)
## Cq versus genome cov 

genome<- read_csv("processed/all-runs_sequencing_stats_prepared_withCq.csv")


# plot of cq by genome cov
ggplot(
  genome,
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




