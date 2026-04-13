# Do damaging earthquakes worldwide respond to lake storage changes?

Code for the manuscript:

**Li, Linxuan & Dong, Junjie** (submitted)
*Do damaging earthquakes worldwide respond to lake storage changes?*

This repository contains code to reproduce the statistical analyses and figures for the paper on the relationship between global lake storage changes and nearby damaging earthquakes.

## Overview

This study tests whether large and shallow earthquakes (M 5.5+ and depth <= 30 km) are statistically associated with lake and reservoir storage changes at the global scale. This repository includes:

- global lake storage time series for 1,972 large inland lakes
- a declustered global earthquake catalog
- statistical tests comparing earthquake occurrence during:
  - generally high vs. generally low storage periods
  - high vs. low storage
  - increasing vs. decreasing storage

The main finding is that earthquakes within **20 km of reservoirs** preferentially occur at **high-storage** and **decreasing-storage**, while no comparable signal is observed for natural lakes.

## Repository structure

```text
.
├── Processed_Data/
│   └── processed global lake data and declustered earthqauke catalogue used by the analysis
├── figs_submitted/
│   └── manuscript figures
├── A0_mc.m
├── A_Decluster.m
├── B_read_Lake_general_info.m
├── C_read_water_time_series.m
├── D_Look_water_time_series.m
├── E_Number_EQ.m
├── F_EF_longterm.m
├── G_EF_storgae.m
├── G_EF_storgae_rate.m
├── LICENSE
└── README.md
