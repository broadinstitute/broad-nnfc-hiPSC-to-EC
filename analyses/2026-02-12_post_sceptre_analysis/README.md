# Render day 0 analysis

Rscript -e 'rmarkdown::render(
  "analyses/2026-02-12_post_sceptre_analysis/create_distance_by_effect_size_plot.Rmd",
  params = list(
    selected_day = 0
  ),
  output_dir = "results/2026-02-12_post_sceptre_analysis/day0"
)'

# Render day 2 analysis

Rscript -e 'rmarkdown::render(
  "analyses/2026-02-12_post_sceptre_analysis/create_distance_by_effect_size_plot.Rmd",
  params = list(
    selected_day = 2
  ),
  output_dir = "results/2026-02-12_post_sceptre_analysis/day2"
)'

# Render day 4 analysis

Rscript -e 'rmarkdown::render(
  "analyses/2026-02-12_post_sceptre_analysis/create_distance_by_effect_size_plot.Rmd",
  params = list(
    selected_day = 4
  ),
  output_dir = "results/2026-02-12_post_sceptre_analysis/day4"
)'