## Plot of SNP positions and amplicon posisitons

snps <- read.csv("processed/snp_positions_to_baseline.csv")

snp_counts <- snps %>%
  count(
    position_1_based,
    baseline_base,
    specimen_base,
    name = "n"
  ) %>%
  arrange(desc(n))


# Read primer BED

primers <- read_tsv(
  "raw/rabvSEasia.scheme.bed",
  col_names = FALSE,
  comment = "#"
) %>%
  select(
    chrom = X1,
    start = X2,
    end = X3,
    primer = X4,
    strand = X6
  ) 
# %>%
#   mutate(
#     amplicon = as.integer(str_extract(primer, "(?<=SEasia_rabv_)\\d+")),
#     pool = if_else(amplicon %% 2 == 1, "Pool 1", "Pool 2")
#   )



# Derive amplicon boundaries and pool from primer coordinates
amplicons <- primers %>%
  mutate(
    amplicon = as.integer(
      str_extract(primer, "(?<=SEasia_rabv_)\\d+")
    ),
    pool = if_else(amplicon %% 2 == 1, "Pool 1", "Pool 2")
  ) %>%
  group_by(amplicon, pool) %>%
  summarise(
    start = start[strand == "+"],
    end = end[strand == "-"],
    amplicon_length = end - start + 1,
    .groups = "drop"
  )

# Position amplicons above the SNP data
y_max <- max(snp_counts$n, na.rm = TRUE)

amplicons <- amplicons %>%
  mutate(
    y_min = if_else(pool == "Pool 1", y_max + 1, y_max + 3),
    y_max_block = if_else(pool == "Pool 1", y_max + 1.7, y_max + 3.7),
    y_label = if_else(pool == "Pool 1", y_max + 1.35, y_max + 3.35)
  )

p_snp_amplicons <- ggplot() +
  
  # Lollipop stems
  geom_segment(
    data = snp_counts,
    aes(
      x = position_1_based,
      xend = position_1_based,
      y = 0,
      yend = n
    ),
    linewidth = 0.5,
    colour = "grey40"
  ) +
  
  # Lollipop heads
  geom_point(
    data = snp_counts,
    aes(
      x = position_1_based,
      y = n
    ),
    size = 2.5
  ) +
  
  # Amplicon regions
  geom_rect(
    data = amplicons,
    aes(
      xmin = start,
      xmax = end,
      ymin = y_min,
      ymax = y_max_block
    ),
    fill = "grey75",
    colour = "grey30"
  ) +
  
  # Amplicon labels
  geom_text(
    data = amplicons,
    aes(
      x = (start + end) / 2,
      y = y_label,
      label = amplicon
    ),
    size = 3,
    fontface = "bold"
  ) +
  
  labs(
    x = "Genome position",
    y = "Number of nucleotide differences"
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 13)


# Publication-ready plot
p_snp_amplicons <- p_snp_amplicons +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ); p_snp_amplicons

# PNG
ggsave(
  "figures/snp_amplicon_plot.png",
  plot = p_snp_amplicons,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# PDF
ggsave(
  "figures/snp_amplicon_plot.pdf",
  plot = p_snp_amplicons,
  width = 180,
  height = 110,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# TIFF
ggsave(
  "figures/snp_amplicon_plot.tiff",
  plot = p_snp_amplicons,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Publication-ready plot
p_snp_amplicons <- p_snp_amplicons +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ); p_snp_amplicons
 



p_snp <- ggplot() +
  
  # Lollipop stems
  geom_segment(
    data = snp_counts,
    aes(
      x = position_1_based,
      xend = position_1_based,
      y = 0,
      yend = n
    ),
    linewidth = 0.5,
    colour = "grey40"
  ) +
  
  # Lollipop heads
  geom_point(
    data = snp_counts,
    aes(
      x = position_1_based,
      y = n
    ),
    size = 2.5
  ) +
  
  # # Amplicon regions
  # geom_rect(
  #   data = amplicons,
  #   aes(
  #     xmin = start,
  #     xmax = end,
  #     ymin = y_min,
  #     ymax = y_max_block
  #   ),
  #   fill = "grey75",
  #   colour = "grey30"
  # ) +
  
  # # Amplicon labels
  # geom_text(
  #   data = amplicons,
  #   aes(
  #     x = (start + end) / 2,
  #     y = y_label,
  #     label = amplicon
  #   ),
  #   size = 3,
  #   fontface = "bold"
  # ) +
  
  labs(
    x = "Genome position",
    y = "No of SNPs"
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 13)


# Publication-ready plot
p_snp <- p_snp +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ); p_snp
p_snp_add <- p_snp +  theme_classic(base_size = 12)

# PNG
ggsave(
  "figures/snp_plot.png",
  plot = p_snp,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# PDF
ggsave(
  "figures/snp_plot.pdf",
  plot = p_snp,
  width = 180,
  height = 110,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# TIFF
ggsave(
  "figures/snp_plot.tiff",
  plot = p_snp,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
