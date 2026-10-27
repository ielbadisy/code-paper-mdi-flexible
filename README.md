# Simulation code for missing data imputation in flexible survival models

Code and scripts for the simulation study and the application in the paper
**"Classical and machine learning-based imputation of missing covariates in
flexible survival models with time-varying and non-linear effects"**
(El Badisy and Giorgi, under consideration in _Biometrical Journal_).

## Contents

| Folder | Content | Paper output |
|---|---|---|
| `sim_PH_scenario/` | Scenario 1 (proportional hazards) | Section 4.1, Appendix A |
| `sim_TVE_bathub_scenario/` | Scenario 2 (bathtub TVE of x3) | Section 4.2, Appendix B |
| `sim_TVE_decreasing_scenario/` | Scenario 3 (decreasing TVE of x3) | Section 4.2, Appendix B |
| `sim_TVE_increasing_scenario/` | Scenario 4 (increasing TVE of x3) | Figures 1-3, Appendix B |
| `sim_NLE_scenario/` | Scenario 5 (NLE of x2) | Figures 4-6, Appendix C |
| `case_study/` | Colorectal cancer application | Figures 7-8, Appendix D |

Each simulation folder contains:

- `utils_*.R`: data generation, amputation (MCAR, MAR, MNAR), imputation
  methods, and estimation of the flexible parametric survival model
  (`rstpm2::stpm2`).
- `sim_*.R`: runs the scenario (R = 500 replicates, N = 1000 observations per
  replicate) and writes `res_*.csv`.
- `res_*.R` / `results_PH.R`: read `res_*.csv` and produce the figures of
  hazard ratios, empirical coverage, and MSE.

`sim_NLE_scenario/original_NLE.csv` holds the true hazard ratio curve used by
the NLE result scripts.

The application uses `case_study/CRC.csv` (343 patients), a previously
published dataset openly available on Figshare:
<https://figshare.com/articles/dataset/dataset_CRC_Fes_csv/21821439?file=38720988>.
`case_study_NLE_Age.R` prints the AIC of the age models (linear, natural
spline df 1 to 6) and produces Figure 7, `case_study_TVE_Surgery.R` produces
Figure 8, and `case_study_Models.R` fits the models on the full dataset.

## Installation

All dependencies are on CRAN except `imputer`, the wrapper used to call the
imputation methods, which is installed from GitHub:

```r
install.packages(c("pacman", "tidyverse", "vroom", "data.table", "rstpm2",
  "survival", "survAUC", "biostat3", "mice", "miceadds", "mixgb",
  "missCforest", "gower", "fastDummies", "purrr", "pbapply", "tictoc"))
remotes::install_github("ielbadisy/imputer")
```

## Running

Open a scenario folder as the working directory (each has its own `.Rproj`)
and run `sim_*.R`, then the `res_*.R` scripts. The simulation scripts run the
replicates in parallel through `pbapply` with `nCores` cores (set at the top of
each script; adjust to your machine). Random seeds are fixed in each script.

## License

MIT License, see `LICENSE.md`.
