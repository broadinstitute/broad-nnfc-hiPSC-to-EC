# How to run perturbplan

## Create mudata from sceptre object

Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_perturbplan/prepare_matrix_and_run_perturbplan.Rmd",
  params = list(
    selected_day = 0
  ),
  output_dir = "results/2025-09-24_power_analysis/"
)'


## Run the power analysis
Requires:
devtools::install_github("Katsevich-Lab/sceptreIGVF")
devtools::install_github("Katsevich-Lab/perturbplan")

Rscript analyses/2026-02-10_perturbplan/run_power_analysis.R \
  -i results/2026-02-10_perturbplan/day0/mudata.h5mu \
  -o results/2026-02-10_perturbplan/day0/perturbplan_es_15_pval_0.1_both.h5mu \
  --effect_size 0.85 \
  --effect_size_sd 0.13 \
  -t 8 \
  --pvalue_cutoff 0.1 \
  --side both