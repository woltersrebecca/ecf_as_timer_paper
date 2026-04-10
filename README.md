# ECF Sigma and Anti-Sigma Factor Timer Paper — Data Analysis

**Authors:** Rebecca Wolters, Stefano Vecchione, Angelika Diel, Georg Fritz  
**Institution:** The University of Western Australia, University of Marburg
**Status:** In preparation  

---

## Overview

This repository contains all data analysis scripts and raw plate reader data
for the paper:

> Wolters R. et al. (202X). *Title of paper*. Journal Name.

The paper characterises transcriptional timer circuits built from ECF sigma
and anti-sigma factors. Analysis covers:

- **Toxicity analysis** — growth effects of sigma and anti-sigma
  factor expression in E. coli
- **Zero-step timer** — characterisation of inducible promoters (Ptet, Pbad, PpcaU)
  driving luciferase directly
- **Single-step timer** — one-step ECF cascade with anti-sigma threshold control
- **Two-step timer** — two-step ECF cascade with combined aTc and DHBA thresholds


---

## Raw data

Raw data are plate reader measurements from a BioTek plate reader exported
as .csv files using Gen5 software. Each experiment folder contains:

- `*_exp#.csv` — raw OD600 and luminescence measurements
- `metainfo.csv` — plate layout in plater format describing strain,
  inducer concentrations, and well assignments per experiment

---

## Reproducing the analysis

### Requirements

R >= 4.2.0

Install required packages:

```r
install.packages(c(
  "tidyverse", "janitor", "MetBrewer", "broom",
  "scales", "wesanderson", "patchwork", "here",
  "plater", "rstatix", "reluxr"
))

# wellr must be installed from GitHub
install.packages("wellr", repos = "bradyajohnston.r-universe.dev")
```

### Steps

1. Clone this repository
```bash
git clone https://github.com/woltersrebecca/ecf_as_timer_paper.git
```

2. Open `ecf_as_timer_paper.Rproj` in RStudio

3. Run the analysis notebooks in `analysis/` — each is self-contained:
   - `zerostep/Zerostep_timer_analysis_clean.Rmd`
   - `singlestep/Singlestep_timer_analysis_clean.Rmd`
   - `twostep/Twostep_timer_analysis_clean.Rmd`
   - `toxicity/as_and_ecf_tox_analysis_clean.Rmd`

All file paths use the `here` package and resolve relative to the
`.Rproj` file, so no manual path editing is needed.

---

## Key analysis decisions

A few methodological choices worth noting for reproducibility:

**Time delay calculation** — time delays are calculated as the time point at
which the luminescence output of the ECF promoter is fourfold over the pre-induction 
level.

**Fold-change flooring** — luminescence fold-change values below 1 are set
to 1 prior to analysis, as values at or below baseline represent measurement
noise rather than a biological signal.

**OD600 normalisation** — all luminescence values are normalised by OD600
to account for differences in cell density between wells and time points.

**Background correction** — OD600 is corrected for pathlength using
empirically determined correction factors (96-well: 4.898, 384-well: 3.09).

---

## Statistical tests

- **Fold-change significance** — paired t-test on log10(fold-change) across
  biological replicates, comparing threshold ON vs OFF per strain
- **Growth rate significance** — one-way ANOVA with Tukey HSD post-hoc test
  across inducer combinations per strain

---
