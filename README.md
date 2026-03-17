# FTA and RDT Sample Sequencing for Rabies Virus Whole Genome Analysis

## Overview

This repository contains data and analysis code for experiments assessing the use of **FTA cards** and **Rapid Diagnostic Test (RDT) strips** as input materials for rabies virus whole genome sequencing.

The study evaluates how different **sample preservation methods and storage conditions** affect sequencing performance over time, using frozen brain tissue as the baseline comparator.

---

## Experimental Rationale

This experiment is designed as a **post-collection stress test** to assess RNA stability and sequencing success under different storage conditions.

- **Baseline comparator**: Defrosted frozen brain tissue (−80°C), extracted at Day 0
- **Test conditions**: Brain tissue applied to FTA cards and RDT strips
- **Objective**: Quantify degradation and sequencing performance relative to baseline

---

## Experimental Design

### Sample Processing

For each brain sample:

- One aliquot is extracted immediately (**Day 0 baseline**)
- Remaining material is applied to:
  - FTA cards
  - RDT strips

---

### Storage Conditions

Each sample type is stored under:

- Room Temperature (RT)
- 37°C

---

### Timepoints

| Timepoint | Description |
|----------|------------|
| T0       | Baseline brain extraction |
| T1       | Day 1 |
| T3       | Day 3 |
| T7       | Day 7 |
| T21      | Day 21 |

Note: Day 0 extraction is only feasible for brain tissue. FTA/RDT begin at Day 1.

---

### Experimental Matrix (per sample)

| Material | Temperature | T0 | T1 | T3 | T7 | T21 |
|----------|------------|----|----|----|----|-----|
| Brain tissue | −80°C | ✓ | – | – | – | – |
| FTA | RT | – | ✓ | ✓ | ✓ | ✓ |
| FTA | 37°C | – | ✓ | ✓ | ✓ | ✓ |
| RDT | RT | – | ✓ | ✓ | ✓ | ✓ |
| RDT | 37°C | – | ✓ | ✓ | ✓ | ✓ |

---

### Extraction Strategy

- One extract per **brain × temperature × timepoint**
- Total extracts per sample: **17**
- RDT strips:
  - ~2 extracts per strip
  - ~4 strips required per sample

---

### Sequencing Runs

- 5 sequencing runs (one per timepoint)
- T0: ~10 samples + controls
- Other timepoints: ~40 samples + controls

---

## QC Metrics

The following metrics are recorded for each sample:

### RNA Quality
- Nanodrop concentration
- RNA Integrity Number (RIN)

### Viral Load
- RT-qPCR Ct value

### Sequencing Performance
- % genome covered at ≥50×
- Mean sequencing depth
- Number of Ns
- Ambiguous base fraction

---

## Controls and Contamination Prevention

### Controls

- Blank FTA cards (processed alongside samples)
- Extraction negatives (per batch)
- Positive extraction control (e.g. CVS strain or high viral load sample)
- Negative control (dFAT-negative brain sample)

### Laboratory Practices

- Dedicated clean area for FTA processing
- Single-use sterile blades for card punching
- Gloves changed between samples
- Surfaces decontaminated between samples

---

## Data Structure

### Raw Data

Per-position sequencing depth:

| Sample | Position | Depth |
|--------|----------|-------|

Sample naming convention:

```
SampleID-Material-Temperature-Timepoint
```

Examples:

- `S1-Brain-D0`
- `S1-FTA-RT-D1`
- `S1-RDT-37-D21`

---

### Parsed Metadata

| Column | Description |
|--------|------------|
| sample_id | Biological sample ID |
| sample_type | Brain / FTA / RDT |
| temp | Storage temperature |
| timepoint | Days post sampling |

---

### Derived QC Metrics

Generated per sample:

| Metric | Description |
|--------|------------|
| mean_depth | Average sequencing depth |
| pct_10x | % genome ≥10× |
| pct_50x | % genome ≥50× |
| pct_100x | % genome ≥100× |
| depth_ratio | Relative to baseline |
| cov_loss | Loss vs baseline |

---

## Analysis Overview

The analysis pipeline:

1. Parse sample metadata from sample names
2. Filter failed sequencing libraries (optional)
3. Calculate per-sample QC metrics
4. Compare each condition to baseline brain sample
5. Visualise degradation trends across:
   - Timepoints
   - Storage conditions
   - Sample types

---

## Key Outputs

### Plots

- Genome coverage profiles
- % genome ≥50× vs time
- Depth decay vs time
- Heatmaps of sequencing success

### Tables

- Condition-level summaries (median depth, coverage)
- Success rates across conditions

---

## Key Questions Addressed

- Can FTA cards and RDT strips preserve rabies RNA for whole genome sequencing?
- How does sequencing success degrade over time?
- What is the impact of elevated temperature (37°C)?
- How do FTA and RDT compare as preservation methods?

---

## Reproducibility

All analysis is performed using R (tidyverse-based workflow).  
Scripts are provided to:

- Parse and clean data
- Compute QC metrics
- Generate figures and summary tables

---

## Notes

- Baseline brain samples serve as the **reference consensus**
- All comparisons are made **relative to baseline per sample**
- Failed sequencing libraries may be excluded from analysis where appropriate

---

## Contact

For questions or collaboration, please contact the repository maintainer.