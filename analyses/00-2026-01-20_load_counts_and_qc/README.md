
# Render day 0 analysis

Rscript -e 'rmarkdown::render(
  "analyses/00-2026-01-20_load_counts_and_qc/cdna-load-counts-and-qc.Rmd",
  params = list(
    selected_day = 0,
    ribo_threshold = 20,
    target_umi_threshold = 55,
    umi_threshold = 1500
  ),
  output_dir = "results/2026-01-20_load_counts_and_qc/day0_ribo20_tumi55_umi1500"
)'


# Render day 2 analysis
Rscript -e 'rmarkdown::render(
  "analyses/00-2026-01-20_load_counts_and_qc/cdna-load-counts-and-qc.Rmd",
  params = list(
    selected_day = 2,
    ribo_threshold = 16,
    target_umi_threshold = 45,
    umi_threshold = 1000
  ),
  output_dir = "results/2026-01-20_load_counts_and_qc/day2_ribo16_tumi45_umi1000"
)'

# Render day 4 analysis
Rscript -e 'rmarkdown::render(
  "analyses/00-2026-01-20_load_counts_and_qc/cdna-load-counts-and-qc.Rmd",
  params = list(
    selected_day = 4,
    ribo_threshold = 6,
    target_umi_threshold = 25,
    umi_threshold = 1000
  ),
  output_dir = "results/2026-01-20_load_counts_and_qc/day4_ribo6_tumi25_umi1000"
)'