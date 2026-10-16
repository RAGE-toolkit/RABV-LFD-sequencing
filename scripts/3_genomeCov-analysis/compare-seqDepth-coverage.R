library(tidyverse)
library(patchwork)
# =========================================================
# Read the combined sequencing depth output from Python
# =========================================================

depths <- read_tsv(
  "processed/all_seq_combined_post_normalisation_coverage.tsv"
)


# =========================================================
# Turn sample columns into rows and extract metadata
# =========================================================

depths_long <- depths %>%
  
  pivot_longer(
    cols = -c(chrom, pos),
    names_to = "sample_name",
    values_to = "depth"
  ) %>%
  
  # Remove controls / unwanted samples
  filter(
    !grepl("NTC|NEC|H26", sample_name)
  ) %>%
  
  mutate(
    
    # -------------------------------------------------------
    # Correct sample naming
    # -------------------------------------------------------
    
    sample_name = if_else(
      sample_name == "N14",
      "N14-D0",
      sample_name
    ),
    
    # -------------------------------------------------------
    # Normalise underscores to dashes
    # -------------------------------------------------------
    
    sample_name_norm = str_replace_all(
      sample_name,
      "_",
      "-"
    ),
    
    # -------------------------------------------------------
    # Identify field samples
    # -------------------------------------------------------
    
    is_field = str_detect(
      sample_name_norm,
      "^RADDL|^H[0-9]"
    )
  ) %>%
  
  # =======================================================
# Extract metadata from sample names
# =======================================================

extract(
  sample_name_norm,
  into = c(
    "sample",
    "condition",
    "temp",
    "day"
  ),
  regex = "^(N\\d+|S\\d+)(?:-(FTA|RDT))?(?:-(RT|37))?(?:-(D\\d+))?$",
  remove = FALSE
) %>%
  
  mutate(
    
    # =====================================================
    # FIELD SAMPLES — Cohort 3
    # =====================================================
    
    sample = if_else(
      is_field,
      str_remove(sample_name, "_RDT$"),
      sample
    ),
    
    condition = if_else(
      is_field,
      if_else(
        str_detect(sample_name, "_RDT$"),
        "RDT",
        "Brain"
      ),
      condition
    ),
    
    temp = if_else(
      is_field,
      "Field",
      temp
    ),
    
    # =====================================================
    # General metadata
    # =====================================================
    
    condition = if_else(
      is.na(condition) | condition == "",
      "Baseline",
      condition
    ),
    
    day = factor(
      day,
      levels = c(
        "D0",
        "D1",
        "D3",
        "D7",
        "D21"
      )
    ),
    
    # =====================================================
    # Define cohorts
    # =====================================================
    
    cohort = case_when(
      str_detect(sample, "^S") ~ "Cohort 1",
      str_detect(sample, "^N") ~ "Cohort 2",
      is_field ~ "Cohort 3",
      TRUE ~ NA_character_
    )
  ) %>%
  
  select(-sample_name_norm)


# =========================================================
# Define cohorts
# =========================================================

cohort1 <- depths_long %>%
  filter(cohort == "Cohort 1")

cohort2 <- depths_long %>%
  filter(cohort == "Cohort 2")

cohort3 <- depths_long %>%
  filter(cohort == "Cohort 3")


# =========================================================
# COHORT 1 — genome-wide sequencing depth
# =========================================================

ggplot(
  cohort1,
  aes(
    x = pos,
    y = depth,
    colour = condition,
    linetype = day,
    group = sample_name
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  scale_linetype_manual(
    values = c(
      "D0" = "solid",
      "D1" = "dashed",
      "D3" = "dotted",
      "D7" = "dotdash",
      "D21" = "longdash"
    ),
    na.value = "solid"
  ) +
  
  geom_line() +
  
  facet_wrap(
    ~ sample,
    scales = "free_y"
  ) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    linetype = "Day"
  ) +
  
  theme_classic()

# COHORT 1 — show collapsed version for RDT/brain/FTA comparison
# =========================================================

cohort1_summary <- cohort1 %>%
  mutate(
    condition = recode(condition, "Baseline" = "Brain")
  ) %>%
  group_by(condition, pos) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    sd_depth = sd(depth, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  cohort1_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    fill = condition
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth
    ),
    alpha = 0.2,
    colour = NA
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  
  geom_line(linewidth = 0.9) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    fill = "Condition"
  ) +
  
  theme_classic()

# =========================================================
# COHORT 2 — genome-wide sequencing depth
# =========================================================

ggplot(
  cohort2,
  aes(
    x = pos,
    y = depth,
    colour = condition,
    linetype = day,
    group = sample_name
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  scale_linetype_manual(
    values = c(
      "D0" = "solid",
      "D3" = "dotted",
      "D7" = "dotdash"
    ),
    na.value = "solid"
  ) +
  
  geom_line() +
  
  facet_wrap(
    ~ sample,
    scales = "free_y"
  ) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    linetype = "Day"
  ) +
  
  theme_classic()

# COHORT 2 — show collapsed version for RDT/brain/FTA comparison
# =========================================================

cohort2_summary <- cohort2 %>%
  mutate(
    condition = recode(condition, "Baseline" = "Brain")
  ) %>%
  group_by(condition, pos) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    sd_depth = sd(depth, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  cohort2_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    fill = condition
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth
    ),
    alpha = 0.2,
    colour = NA
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  
  geom_line(linewidth = 0.9) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    fill = "Condition"
  ) +
  
  theme_classic()

# =========================================================
# COHORT 3 — genome-wide sequencing depth
# =========================================================

ggplot(
  cohort3,
  aes(
    x = pos,
    y = depth,
    colour = condition,
    group = sample_name
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  geom_line() +
  
  facet_wrap(
    ~ sample,
    scales = "free_y"
  ) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition"
  ) +
  
  theme_classic()

# COHORT 3 — show collapsed version for RDT/brain/FTA comparison
# =========================================================

cohort3_summary <- cohort3 %>%
  mutate(
    condition = recode(condition, "Baseline" = "Brain")
  ) %>%
  group_by(condition, pos) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    sd_depth = sd(depth, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  cohort3_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    fill = condition
  )
) +
  
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth
    ),
    alpha = 0.2,
    colour = NA
  ) +
  scale_colour_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "FTA" = "#0072B2",
      "RDT" = "#D55E00",
      "Brain" = "black"
    )
  ) +
  
  geom_line(linewidth = 0.9) +
  
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    fill = "Condition"
  ) +
  
  theme_classic()


condition_cols <- c(
  "Brain" = "black",
  "FTA" = "#0072B2",
  "RDT" = "#D55E00"
)

condition_levels <- c("Brain", "FTA", "RDT")
cohort1_summary <- cohort1_summary %>%
  mutate(condition = factor(condition, levels = condition_levels))

cohort2_summary <- cohort2_summary %>%
  mutate(condition = factor(condition, levels = condition_levels))

cohort3_summary <- cohort3_summary %>%
  mutate(condition = factor(condition, levels = condition_levels))


# Define shared genome position axis
x_limits <- range(cohort3_summary$pos, na.rm = TRUE)
x_breaks <- scales::pretty_breaks(n = 8)(x_limits)

p_cohort1 <- ggplot(
  cohort1_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    fill = condition
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth,
      fill = condition
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  )+
  geom_line(linewidth = 0.3) +
  scale_colour_manual(
    values = condition_cols,
    limits = c("Brain", "FTA", "RDT"),
    breaks = c("Brain", "FTA", "RDT"),
    labels = c("Brain", "FTA", "LFD"),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = condition_cols,
    guide = "none"
  )+
  scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = c(0.01, 0.01))
  )+
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Specimen type",
    fill = "Specimen type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )


p_cohort2 <- ggplot(
  cohort2_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    fill = condition
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth,
      fill = condition
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  )+
  geom_line(linewidth = 0.3) +
  scale_colour_manual(
    values = condition_cols,
    guide = "none"
  ) +
  scale_fill_manual(
    values = condition_cols,
    guide = "none"
  )+
  scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = c(0.01, 0.01))
  )+
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition",
    fill = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )


p_cohort3 <- ggplot(
  cohort3_summary,
  aes(
    x = pos,
    y = mean_depth,
    colour = condition,
    group = condition
  )
) +
  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_depth - sd_depth, 0),
      ymax = mean_depth + sd_depth,
      fill = condition
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.3
  ) +
  scale_colour_manual(
    values = condition_cols,
    guide = "none"
  ) +
  scale_fill_manual(
    values = condition_cols,
    guide = "none"
  )+
  scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = c(0.01, 0.01))
  )+
  labs(
    x = "Genome position",
    y = "Sequencing depth",
    colour = "Condition"
  ) +
  theme_classic(base_size = 12) 

source("scripts/3_genomeCov-analysis/amplicons_versus_snps.R")

p_snp_add <- p_snp_add+scale_x_continuous(
  limits = x_limits,
  breaks=x_breaks,
  expand = c(0, 0)
)

p_depth_combined <- (
  p_cohort1 + p_cohort2 + p_cohort3 +p_snp_add
) +
  plot_layout(
    ncol = 2
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "top"
  )

p_depth_combined

# Publication-quality outputs

# ggsave(
#   "figures/genome_depth_by_cohort.pdf",
#   plot = p_depth_combined,
#   width = 180,
#   height = 150,
#   units = "mm",
#   device = cairo_pdf,
#   bg = "white"
# )
# 
# ggsave(
#   "figures/genome_depth_by_cohort.tiff",
#   plot = p_depth_combined,
#   width = 180,
#   height = 150,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw",
#   bg = "white"
# )
# 
# ggsave(
#   "figures/genome_depth_by_cohort.png",
#   plot = p_depth_combined,
#   width = 180,
#   height = 150,
#   units = "mm",
#   dpi = 600,
#   bg = "white"
# )

