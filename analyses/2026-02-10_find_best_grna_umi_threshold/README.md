Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_find_best_grna_umi_threshold/kd_effect_size_at_different_threshold.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 10,
    threads = 8
  ),
  output_dir = "results/2026-02-10_find_best_grna_umi_threshold/day0_grna10"
)'


Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_find_best_grna_umi_threshold/kd_effect_size_at_different_threshold.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 20,
    threads = 8
  ),
  output_dir = "results/2026-02-10_find_best_grna_umi_threshold/day0_grna20"
)'


Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_find_best_grna_umi_threshold/kd_effect_size_at_different_threshold.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 40,
    threads = 8
  ),
  output_dir = "results/2026-02-10_find_best_grna_umi_threshold/day0_grna40"
)'


Rscript -e 'rmarkdown::render(
  "analyses/2026-02-10_find_best_grna_umi_threshold/kd_effect_size_at_different_threshold.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 100,
    threads = 8
  ),
  output_dir = "results/2026-02-10_find_best_grna_umi_threshold/day0_grna100"
)'