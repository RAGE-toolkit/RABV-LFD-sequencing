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

## Scripts   
### 1_data-prep
Collate and standardise seq and cq data into dataframes suitable for analyses
R scripts: 
- prepare_sequencing_stats.R : combines sequencing summary statistics (read depth etc) for all cohorts, removes 1 failed sample (S1-FTA-RT-D1), non-study samples and negatives. Standardises temp and storage info for all specimen types.
- prepare_cq_data.R : run after above script to standardise data and add cq values to the sequencing info
Python scripts:
- combine_amplicon_depths.py : used to collate all the amplicon depth files ("*.amplicon_depths.tsv") from bioinformatics pipeline output into one file for analysis
- combine_seq_depths.py : used to collate all the post normalisation (i.e. capped to 200x) per base depth files   

### 2_cq-analysis  
Analysis of change in real-time qPCR cq values from baseline brain measurement
- stats-qc.R : calculate delta_cq (change in cq from baseline) and apply lmer model to measure effect of specimen type, temp and storage duration
- cqVcoverage.R : plot cq values versus genome coverage

### 3_genomeCov-analysis   
Analysis of genome coverage (at read depth of 20x) in response to temp, storage duration and specimen type
Python script:
- compare_sample_specimens.py : compares same-sample specimens to identify SNPs relative to baseline brain sample, using n160 alignment as input. Produces comparison_to_baseline.csv and snp_positions_to_baseline.csv
R scripts:   
- stats-seq.R : summaries of samples achieving >20x, range of exploratory and final plots, beta regression mixed model to assess proportion of genome covered in relation to temp, storage duration and specimen type
- create_model_results_table.R : export all model results to an excel
- all_experiments_plots.R : more plots, including main manuscript figures
- compare-seqDepth-coverage.R : summarise and plot read depth across genome for each cohort and produce plot, including cumulative snp frequency across genome
- amplicons_versus_snps.R : lollipop plots of SNP frequencies across the genome, with amplicons mapped
- compare-amplicon-coverage.R : explore amplicons failures across cohorts and by specimen type
- pairwise_distances.R : calculate pairwise SNP differences for the n160 sequences used in study, outputs within sample distance heatmaps and exploratory plots to look at associations with exp variables

### 4_genomeCov-analysis   
Calculating turnaround times for diagnostic workflows (n=319 RDT samples) and sequencing workflows in MIMAROPA, Philippines
- diag-turnaround-time.R : calculates TAT for LFD and FAT diagnostic results fro samples processed in RDT validation in MIMAROPA
- seq-turnaround-times.R : calculate TAT from point of confirmatory FAT test for sequenced samples from SPEEDIER project


## Tools and protocols

- **Sequencing pipeline:** [Artic-nf](https://github.com/RAGE-toolkit/Artic-nf)
- **Rabies lineage assignment:** [MADDOG V2](https://github.com/RAGE-toolkit/MADDOG-V2)
- **Laboratory protocol:** [RNA extraction from rapid diagnostic tests and FTA cards](https://doi.org/10.17504/protocols.io.14egn5oz6g5d/v1)

## Citation

*Lateral flow devices can fast-track rabies virus genomic surveillance.*

Please cite the associated manuscript when using this repository.