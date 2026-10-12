library(ggplot2)
library(Biostrings)
library(tidyverse)
library(dplyr)

# ============================================================
# Settings
# ============================================================

alignment_file <- "processed/all-phl-specimens-n160.aln.fasta"
metadata_file  <- "processed/all-runs_sequencing_stats_prepared.csv"

# ============================================================
# 1. Import sequencing QC metadata
# ============================================================

qc_all <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

# Create sequence IDs matching the FASTA
sequence_metadata <- qc_all %>%
  mutate(
    sequence = basename(filepath),
    sequence = str_remove(sequence, "\\.sorted\\.bam$")
  ) %>%
  select(
    sequence,
    coverage_20,
    sample_type,
    timepoint,
    temp
  )

print(sequence_metadata)


# Check for duplicate sequence IDs
duplicate_sequences <- sequence_metadata %>%
  count(sequence) %>%
  filter(n > 1)

if (nrow(duplicate_sequences) > 0) {
  print(duplicate_sequences)
  stop("Duplicate sequence IDs found in sequence_metadata.")
}


# ============================================================
# 2. Import alignment
# ============================================================

aln <- readDNAStringSet(alignment_file)

seq_names <- names(aln)

# Remove everything after first space in FASTA header
seq_names <- gsub(" .*", "", seq_names)

cat("Number of sequences:", length(aln), "\n")
cat("Alignment length:", unique(width(aln)), "\n")


# ============================================================
# 3. Extract underlying sample ID
# ============================================================

metadata <- tibble(
  sequence = seq_names
) %>%
  mutate(
    sample_id = case_when(
      
      # S1-D0, S1-FTA-37-D1 etc. -> S1
      str_detect(sequence, "^[SN]\\d+") ~
        str_extract(sequence, "^[SN]\\d+"),
      
      # RADDL4B-2026-016_RDT -> RADDL4B-2026-016
      str_detect(sequence, "_RDT$") ~
        str_remove(sequence, "_RDT$"),
      
      # Anything else remains unchanged
      TRUE ~ sequence
    )
  )

print(metadata)


# ============================================================
# 4. Convert alignment to matrix
# ============================================================

seq_matrix <- do.call(
  rbind,
  strsplit(
    as.character(aln),
    split = ""
  )
)

rownames(seq_matrix) <- seq_names


# ============================================================
# 5. Function to calculate pairwise differences
# ============================================================

calculate_pairwise_differences <- function(
    sequence_names,
    seq_matrix
) {
  
  n <- length(sequence_names)
  
  distance_matrix <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(
      sequence_names,
      sequence_names
    )
  )
  
  for (i in seq_len(n)) {
    
    for (j in i:n) {
      
      x <- seq_matrix[sequence_names[i], ]
      y <- seq_matrix[sequence_names[j], ]
      
      # Only compare positions where BOTH sequences
      # have an unambiguous nucleotide.
      valid <- x %in% c("A", "C", "G", "T") &
        y %in% c("A", "C", "G", "T")
      
      if (sum(valid) == 0) {
        
        d <- NA_real_
        
      } else {
        
        # Number of nucleotide differences
        d <- sum(
          x[valid] != y[valid]
        )
      }
      
      distance_matrix[i, j] <- d
      distance_matrix[j, i] <- d
    }
  }
  
  distance_matrix
}


# ============================================================
# 6. Calculate distances within each sample
# ============================================================

all_distances <- list()

for (sample in unique(metadata$sample_id)) {
  
  sample_sequences <- metadata %>%
    filter(sample_id == sample) %>%
    pull(sequence)
  
  # Need at least two sequences
  if (length(sample_sequences) < 2) {
    next
  }
  
  dist_matrix <- calculate_pairwise_differences(
    sample_sequences,
    seq_matrix
  )
  
  # Convert to long format
  dist_long <- as.data.frame(dist_matrix) %>%
    rownames_to_column("sequence_1") %>%
    pivot_longer(
      cols = -sequence_1,
      names_to = "sequence_2",
      values_to = "differences"
    ) %>%
    mutate(
      sample_id = sample
    )
  
  all_distances[[sample]] <- dist_long
}

all_distances <- bind_rows(all_distances)


# ============================================================
# 7. Keep each pair only once
# ============================================================

distances_unique <- all_distances %>%
  filter(sequence_1 != sequence_2) %>%
  mutate(
    seq_a = pmin(sequence_1, sequence_2),
    seq_b = pmax(sequence_1, sequence_2)
  ) %>%
  distinct(
    sample_id,
    seq_a,
    seq_b,
    .keep_all = TRUE
  )


# ============================================================
# 8. Add metadata for BOTH sequences
# ============================================================

distances_with_metadata <- distances_unique %>%
  
  # Metadata for sequence 1
  left_join(
    sequence_metadata %>%
      dplyr::rename(
        sequence_1 = sequence,
        coverage_20_1 = coverage_20,
        sample_type_1 = sample_type,
        timepoint_1 = timepoint,
        temp_1 = temp
      ),
    by = "sequence_1"
  ) %>%
  
  # Metadata for sequence 2
  left_join(
    sequence_metadata %>%
      dplyr::rename(
        sequence_2 = sequence,
        coverage_20_2 = coverage_20,
        sample_type_2 = sample_type,
        timepoint_2 = timepoint,
        temp_2 = temp
      ),
    by = "sequence_2"
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------
    
    min_coverage =
      pmin(
        coverage_20_1,
        coverage_20_2,
        na.rm = TRUE
      ),
    
    low_coverage_pair =
      coverage_20_1 < 90 |
      coverage_20_2 < 90,
    
    # --------------------------------------------------------
    # Pair-level temperature
    # --------------------------------------------------------
    
    temp_pair = case_when(
      
      temp_1 == temp_2 ~ temp_1,
      
      TRUE ~ paste(
        sort(
          c(temp_1, temp_2)
        ),
        collapse = " + "
      )
    ),
    
    # --------------------------------------------------------
    # Pair-level sample type
    # --------------------------------------------------------
    
    sample_type_pair = case_when(
      
      sample_type_1 == sample_type_2 ~ sample_type_1,
      
      TRUE ~ paste(
        sort(
          c(sample_type_1, sample_type_2)
        ),
        collapse = " + "
      )
    ),
    
    # --------------------------------------------------------
    # Storage time
    # --------------------------------------------------------
    
    day_1 =
      as.numeric(
        str_remove(timepoint_1, "^D")
      ),
    
    day_2 =
      as.numeric(
        str_remove(timepoint_2, "^D")
      ),
    
    time_difference =
      abs(day_1 - day_2)
  )


# ============================================================
# 9. Check metadata matching
# ============================================================

cat(
  "Pairs with missing metadata:",
  sum(
    is.na(distances_with_metadata$coverage_20_1) |
      is.na(distances_with_metadata$coverage_20_2)
  ),
  "\n"
)

if (
  any(
    is.na(distances_with_metadata$coverage_20_1) |
    is.na(distances_with_metadata$coverage_20_2)
  )
) {
  
  print(
    distances_with_metadata %>%
      filter(
        is.na(coverage_20_1) |
          is.na(coverage_20_2)
      ) %>%
      select(
        sequence_1,
        sequence_2,
        coverage_20_1,
        coverage_20_2
      ) %>%
      distinct()
  )
}


# ============================================================
# 10. Save pairwise dataset
# ============================================================

write_csv(
  distances_with_metadata,
  "processed/within_sample_pairwise_differences_with_metadata.csv"
)


# ============================================================
# 11. Overall summary
# ============================================================

distances_with_metadata %>%
  summarise(
    n_comparisons = n(),
    mean = mean(differences, na.rm = TRUE),
    sd = sd(differences, na.rm = TRUE),
    median = median(differences, na.rm = TRUE),
    IQR = IQR(differences, na.rm = TRUE),
    min = min(differences, na.rm = TRUE),
    max = max(differences, na.rm = TRUE),
    n_identical =
      sum(
        differences == 0,
        na.rm = TRUE
      ),
    prop_identical =
      mean(
        differences == 0,
        na.rm = TRUE
      )
  )


# ============================================================
# 12. Summary by coverage
# ============================================================

distances_with_metadata %>%
  group_by(low_coverage_pair) %>%
  summarise(
    n_comparisons = n(),
    median = median(differences, na.rm = TRUE),
    IQR = IQR(differences, na.rm = TRUE),
    mean = mean(differences, na.rm = TRUE),
    max = max(differences, na.rm = TRUE),
    prop_identical =
      mean(
        differences == 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


# ============================================================
# 13. Overall histogram
# ============================================================

ggplot(
  distances_with_metadata,
  aes(
    x = differences
  )
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    colour = "white"
  ) +
  labs(
    x = "Pairwise nucleotide differences",
    y = "Number of comparisons"
  ) +
  theme_classic(base_size = 11)

# ============================================================
# 14. Pairwise differences by sample
# ============================================================

plot_distances <- distances_with_metadata %>%
  mutate(
    cohort_order = case_when(
      str_detect(sample_id, "^S\\d+") ~ 1,
      str_detect(sample_id, "^N\\d+") ~ 2,
      str_detect(sample_id, "^RADDL") ~ 3,
      TRUE ~ 4
    ),
    
    sample_number = case_when(
      str_detect(sample_id, "^[SN]\\d+") ~
        as.numeric(
          str_extract(sample_id, "\\d+")
        ),
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(
    cohort_order,
    sample_number,
    sample_id
  ) %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = unique(sample_id)
    )
  )


p_pairwise <- ggplot(
  plot_distances,
  aes(
    x = sample_id,
    y = differences
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.6,
    size = 1.8
  ) +
  labs(
    x = "Sample",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
     axis.text.y = element_text(size = 12),
     axis.title = element_text(size = 12),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  )

p_pairwise
ggsave(
  "figures/pairwise_nucleotide_differences.pdf",
  plot = p_pairwise,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  "figures/pairwise_nucleotide_differences.tiff",
  plot = p_pairwise,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)
ggsave(
  "pairwise_nucleotide_differences.png",
  plot = p_pairwise,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# ============================================================
# 15. Distance vs minimum genome coverage
# ============================================================

ggplot(
  distances_with_metadata,
  aes(
    x = min_coverage,
    y = differences
  )
) +
  geom_jitter(
    width = 0,
    height = 0.1,
    alpha = 0.5
  ) +
  # geom_smooth(
  #   method = "loess",
  #   se = TRUE
  # ) +
  labs(
    x = "Minimum genome coverage of pair (%)",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 11)


# ============================================================
# 16. Distance vs minimum coverage, coloured by temperature
# ============================================================
distances_with_metadata <- distances_with_metadata %>%
  mutate(
    temp_comparison = case_when(
      temp_1 == "RT" & temp_2 == "RT" ~ "RT vs RT",
      temp_1 == "37" & temp_2 == "37" ~ "37°C vs 37°C",
      temp_1 %in% c("RT", "37") &
        temp_2 %in% c("RT", "37") &
        temp_1 != temp_2 ~ "RT vs 37°C",
      temp_1 == "Field" & temp_2 == "Field" ~ "Field vs Field",
      TRUE ~ NA_character_
    )
  )
ggplot(
  distances_with_metadata %>%
    filter(!is.na(temp_comparison)),
  aes(
    x = min_coverage,
    y = differences,
    colour = temp_comparison
  )
) +
  geom_jitter(
    width = 0,
    height = 0.1,
    alpha = 0.6
  ) +
  labs(
    x = "Minimum genome coverage of pair (%)",
    y = "Pairwise nucleotide differences",
    colour = "Temperature comparison"
  ) +
  theme_classic(base_size = 11)


# ============================================================
# 17. Distance vs minimum coverage, coloured by sample type
# ============================================================
distances_with_metadata <- distances_with_metadata %>%
  mutate(
    specimen_comparison = case_when(
      sample_type_1 == "Brain" &
        sample_type_2 == "Brain" ~ "Brain vs Brain",
      
      sample_type_1 == "FTA" &
        sample_type_2 == "FTA" ~ "FTA vs FTA",
      
      sample_type_1 == "RDT" &
        sample_type_2 == "RDT" ~ "RDT vs RDT",
      
      sample_type_1 %in% c("FTA", "RDT") &
        sample_type_2 %in% c("FTA", "RDT") &
        sample_type_1 != sample_type_2 ~ "FTA vs RDT",
      
      sample_type_1 %in% c("Brain", "FTA") &
        sample_type_2 %in% c("Brain", "FTA") &
        sample_type_1 != sample_type_2 ~ "Brain vs FTA",
      
      sample_type_1 %in% c("Brain", "RDT") &
        sample_type_2 %in% c("Brain", "RDT") &
        sample_type_1 != sample_type_2 ~ "Brain vs RDT",
      
      TRUE ~ NA_character_
    )
  )
ggplot(
  distances_with_metadata %>%
    filter(!is.na(specimen_comparison)),
  aes(
    x = min_coverage,
    y = differences
  )
) +
  geom_jitter(
    width = 0,
    height = 0.1,
    alpha = 0.6
  ) +
  facet_wrap(~ specimen_comparison) +
  labs(
    x = "Minimum genome coverage of pair (%)",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 11)

# ============================================================
# 18. Distance vs difference in storage time
# ============================================================

ggplot(
  distances_with_metadata,
  aes(
    x = time_difference,
    y = differences,
    colour = low_coverage_pair
  )
) +
  geom_jitter(
    width = 0.15,
    height = 0.1,
    alpha = 0.6
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "red"
    ),
    labels = c(
      "FALSE" = "Both sequences ≥90% coverage",
      "TRUE" = "At least one sequence <90% coverage"
    ),
    name = "Genome coverage"
  ) +
  labs(
    x = "Difference in storage time (days)",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 11)


# ============================================================
# 19. Pairwise heatmaps
# ============================================================

dir.create(
  "processed/pairwise_heatmaps",
  showWarnings = FALSE
)

for (sample in unique(all_distances$sample_id)) {
  
  plot_data <- all_distances %>%
    filter(sample_id == sample)
  
  # Preserve sequence ordering
  sequence_order <- unique(
    c(
      plot_data$sequence_1,
      plot_data$sequence_2
    )
  )
  
  plot_data <- plot_data %>%
    mutate(
      sequence_1 = factor(
        sequence_1,
        levels = sequence_order
      ),
      
      sequence_2 = factor(
        sequence_2,
        levels = rev(sequence_order)
      )
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = sequence_1,
      y = sequence_2,
      fill = differences
    )
  ) +
    geom_tile(
      colour = "white"
    ) +
    scale_fill_viridis_c(
      name = "Nucleotide\ndifferences",
      na.value = "grey90"
    ) +
    geom_text(
      aes(
        label = differences
      ),
      size = 3
    ) +
    coord_fixed() +
    labs(
      title = sample,
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.text.y = element_text(),
      plot.title = element_text(
        face = "bold"
      )
    )
  
  print(p)
  
  ggsave(
    filename = file.path(
      "processed/pairwise_heatmaps",
      paste0(
        sample,
        "_pairwise_differences.pdf"
      )
    ),
    plot = p,
    width = 7,
    height = 6
  )
}

sequence_metadata <- qc_all %>%
  mutate(
    sequence = basename(filepath),
    sequence = str_remove(sequence, "\\.sorted\\.bam$")
  ) %>%
  select(
    sequence,
    coverage_20,
    meanReads,
    sample_type,
    timepoint,
    temp,
    cohort
  )

distances_with_metadata <- distances_unique %>%
  
  # Metadata for sequence 1
  left_join(
    sequence_metadata %>%
      dplyr::rename(
        sequence_1 = sequence,
        coverage_20_1 = coverage_20,
        mean_depth_1 = meanReads,
        sample_type_1 = sample_type,
        timepoint_1 = timepoint,
        temp_1 = temp,
        cohort_1 = cohort
      ),
    by = "sequence_1"
  ) %>%
  
  # Metadata for sequence 2
  left_join(
    sequence_metadata %>%
      dplyr::rename(
        sequence_2 = sequence,
        coverage_20_2 = coverage_20,
        mean_depth_2 = meanReads,
        sample_type_2 = sample_type,
        timepoint_2 = timepoint,
        temp_2 = temp,
        cohort_2 = cohort
      ),
    by = "sequence_2"
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------
    
    min_coverage = pmin(
      coverage_20_1,
      coverage_20_2,
      na.rm = TRUE
    ),
    
    low_coverage_pair =
      coverage_20_1 < 90 |
      coverage_20_2 < 90,
    
    # --------------------------------------------------------
    # Mean depth of the pair
    # --------------------------------------------------------
    
    mean_depth_pair = rowMeans(
      cbind(mean_depth_1, mean_depth_2),
      na.rm = TRUE
    ),
    
    # --------------------------------------------------------
    # Cohort
    # --------------------------------------------------------
    
    cohort = case_when(
      cohort_1 == cohort_2 ~ as.character(cohort_1),
      TRUE ~ paste(
        sort(c(cohort_1, cohort_2)),
        collapse = " + "
      )
    ),
    
    # --------------------------------------------------------
    # Pair-level temperature
    # --------------------------------------------------------
    
    temp_pair = case_when(
      temp_1 == temp_2 ~ temp_1,
      TRUE ~ paste(
        sort(c(temp_1, temp_2)),
        collapse = " + "
      )
    ),
    
    # --------------------------------------------------------
    # Pair-level sample type
    # --------------------------------------------------------
    
    sample_type_pair = case_when(
      sample_type_1 == sample_type_2 ~ sample_type_1,
      TRUE ~ paste(
        sort(c(sample_type_1, sample_type_2)),
        collapse = " + "
      )
    ),
    
    # --------------------------------------------------------
    # Storage time
    # --------------------------------------------------------
    
    day_1 = as.numeric(
      str_remove(timepoint_1, "^D")
    ),
    
    day_2 = as.numeric(
      str_remove(timepoint_2, "^D")
    ),
    
    time_difference = abs(day_1 - day_2)
  )

ggplot(
  distances_with_metadata,
  aes(
    x = mean_depth_pair,
    y = differences
  )
) +
  geom_jitter(
    width = 0,
    height = 0.1,
    alpha = 0.6,
    size = 2
  ) +
  labs(
    x = "Mean sequencing depth of pair",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 13)

distances_with_metadata %>%
  group_by(cohort) %>%
  summarise(
    n = sum(complete.cases(mean_depth_pair, differences)),
    rho = cor(
      mean_depth_pair,
      differences,
      method = "spearman",
      use = "complete.obs"
    ),
    p = cor.test(
      mean_depth_pair,
      differences,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  )

ggplot(
  distances_with_metadata,
  aes(
    x = mean_depth_pair,
    y = differences
  )
) +
  geom_jitter(
    width = 0,
    height = 0.1,
    alpha = 0.6,
    size = 2
  ) +
  facet_wrap(~ cohort) +
  labs(
    x = "Mean sequencing depth of pair",
    y = "Pairwise nucleotide differences"
  ) +
  theme_classic(base_size = 13)

