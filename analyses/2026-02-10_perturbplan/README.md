# perturbplan workflow

This directory contains the steps used to:

1. convert a SCEPTRE object into MuData / H5MU for perturbplan
2. run perturbplan power analysis on the resulting MuData object
3. run perturbplan power analysis directly from a final SCEPTRE object

## Prerequisites

Install required R packages:

```r
devtools::install_github("Katsevich-Lab/sceptreIGVF")
devtools::install_github("Katsevich-Lab/perturbplan")
```

If you render the R Markdown notebook, `pandoc` must be available because the notebook produces HTML output in addition to the H5MU file.

## Create MuData from a SCEPTRE object

Render the notebook once per day. The H5MU output is written under `results/2026-04-24_perturbplan/dayX/h5mu/`.

### Day 0

```sh
Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_perturbplan/00_prepare_matrix_and_run_perturbplan.Rmd",
  params = list(selected_day = 0),
  output_dir = "results/2026-04-24_perturbplan/day0/"
)'
```

### Day 2

```sh
Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_perturbplan/00_prepare_matrix_and_run_perturbplan.Rmd",
  params = list(selected_day = 2),
  output_dir = "results/2026-04-24_perturbplan/day2/"
)'
```

### Day 4

```sh
Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_perturbplan/00_prepare_matrix_and_run_perturbplan.Rmd",
  params = list(selected_day = 4),
  output_dir = "results/2026-04-24_perturbplan/day4/"
)'
```

Expected H5MU outputs:

- `results/2026-04-24_perturbplan/day0/h5mu/day0_grna20_mudata.h5mu`
- `results/2026-04-24_perturbplan/day2/h5mu/day2_grna20_mudata.h5mu`
- `results/2026-04-24_perturbplan/day4/h5mu/day4_grna20_mudata.h5mu`

## Run the power analysis

If `--pvalue_cutoff` is omitted, `run_power_analysis.R` now derives it from `metadata(mudata)$test_results` by taking the largest discovery `p_value` among rows marked `significant == TRUE`.

### Day 0

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis.R \
  -i results/2026-04-24_perturbplan/day0/h5mu/day0_grna20_mudata.h5mu \
  -o results/2026-04-24_perturbplan/day0/h5mu/perturbplan_es_0.85_significant_cutoff_both_mean_theta.h5mu \
  --expression_stats_output results/2026-04-24_perturbplan/day0/expression_stats_mean_theta.tsv \
  --power_results_output results/2026-04-24_perturbplan/day0/individual_power.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  -t 8 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```

### Day 2

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis.R \
  -i results/2026-04-24_perturbplan/day2/h5mu/day2_grna20_mudata.h5mu \
  -o results/2026-04-24_perturbplan/day2/h5mu/perturbplan_es_0.85_significant_cutoff_both_mean_theta.h5mu \
  --expression_stats_output results/2026-04-24_perturbplan/day2/expression_stats_mean_theta.tsv \
  --power_results_output results/2026-04-24_perturbplan/day2/individual_power.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  -t 8 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```

### Day 4

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis.R \
  -i results/2026-04-24_perturbplan/day4/h5mu/day4_grna20_mudata.h5mu \
  -o results/2026-04-24_perturbplan/day4/h5mu/perturbplan_es_0.85_significant_cutoff_both_mean_theta.h5mu \
  --expression_stats_output results/2026-04-24_perturbplan/day4/expression_stats_mean_theta.tsv \
  --power_results_output results/2026-04-24_perturbplan/day4/individual_power.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  -t 8 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```

## Run the power analysis directly from a final SCEPTRE object

`run_power_analysis_from_sceptre_object.R` uses the final `sceptre_object.rds` directly, applies `cells_in_use`, derives the default cutoff from `sceptre_object@discovery_result`, and uses the stored theta values in `sceptre_object@response_precomputations` with a fitted mean reconstructed from `fitted_coefs` and the filtered covariate matrix. Discovery pairs whose guide targets are absent from the guide assignments or whose response genes lack stored theta are dropped with a message. The script writes a standalone RDS, a gene-level expression-stats TSV, and a per-pair power TSV.

### Day 0

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis_from_sceptre_object.R \
  -i results/sceptre/day0_grna20/sceptre_object.rds \
  -o results/2026-04-24_perturbplan/day0/perturbplan_es_0.85_significant_cutoff_from_sceptre_fitted_only.rds \
  --expression_stats_output results/2026-04-24_perturbplan/day0/expression_stats_from_sceptre_fitted_only.tsv \
  --power_results_output results/2026-04-24_perturbplan/day0/individual_power_from_sceptre_fitted_only.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```

### Day 2

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis_from_sceptre_object.R \
  -i results/sceptre/day2_grna20/sceptre_object.rds \
  -o results/2026-04-24_perturbplan/day2/perturbplan_es_0.85_significant_cutoff_from_sceptre_fitted_only.rds \
  --expression_stats_output results/2026-04-24_perturbplan/day2/expression_stats_from_sceptre_fitted_only.tsv \
  --power_results_output results/2026-04-24_perturbplan/day2/individual_power_from_sceptre_fitted_only.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```

### Day 4

```sh
Rscript analyses/2026-02-10_perturbplan/run_power_analysis_from_sceptre_object.R \
  -i results/sceptre/day4_grna20/sceptre_object.rds \
  -o results/2026-04-24_perturbplan/day4/perturbplan_es_0.85_significant_cutoff_from_sceptre_fitted_only.rds \
  --expression_stats_output results/2026-04-24_perturbplan/day4/expression_stats_from_sceptre_fitted_only.tsv \
  --power_results_output results/2026-04-24_perturbplan/day4/individual_power_from_sceptre_fitted_only.tsv \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  --n_nonzero_trt_thresh 15 \
  --n_nonzero_cntrl_thresh 15 \
  --side both
```