# LFD and FTA Card Rabies Virus Sequencing

This repository contains analysis code and data-processing workflows for evaluating lateral flow devices (LFDs) and Flinders Technology Associates (FTA) cards as specimen types for rabies virus (RABV) whole-genome sequencing.

The study investigates the effects of specimen type, storage duration and temperature on RNA detection, genome coverage and consensus sequence concordance, using controlled laboratory experiments and field-collected specimens from the Philippines.

## Study design

Three cohorts were evaluated:

- **Cohort 1:** Five biological samples, with matched FTA-card and LFD specimens stored at ambient temperature (25–30°C) or 37°C for 1, 3, 7 or 21 days.
- **Cohort 2:** Fifteen additional biological samples, with matched FTA-card and LFD specimens stored at ambient temperature for 3 or 7 days.
- **Cohort 3:** Eight field-collected LFDs evaluated alongside matched brain specimens.

Baseline brain tissue was used as the reference specimen for the laboratory experiments.

## Analysis

Analyses include:

- RABV RNA detection by RT-qPCR.
- Genome coverage and sequencing success at ≥20× read depth.
- Effects of specimen type, storage duration and temperature on genome coverage.
- Amplicon-level sequencing performance.
- Consensus sequence concordance between specimens from the same animal.
- Diagnostic performance of LFDs against direct fluorescent antibody testing (DFAT).
- Diagnostic and sequencing turnaround times.

Genome coverage was analysed using beta mixed-effects regression, with biological sample included as a random effect.

## Key findings

- LFDs showed 99.1% sensitivity and 100% specificity against DFAT in 319 paired diagnostic results.
- Of 160 specimens retained for downstream sequencing analysis, 81.9% achieved ≥90% genome coverage at ≥20× depth.
- Genome coverage declined with storage duration, with lower coverage at 37°C than ambient temperature at days 7 and 21.
- Consensus sequences from matched field-collected LFDs and brain specimens were identical in seven of eight pairs.

## Tools and protocols

- **Sequencing pipeline:** [Artic-nf](https://github.com/RAGE-toolkit/Artic-nf)
- **Rabies lineage assignment:** [MADDOG V2](https://github.com/RAGE-toolkit/MADDOG-V2)
- **Laboratory protocol:** [RNA extraction from rapid diagnostic tests and FTA cards](https://doi.org/10.17504/protocols.io.14egn5oz6g5d/v1)

## Citation

*Lateral flow devices can fast-track rabies virus genomic surveillance.*

Please cite the associated manuscript when using this repository.