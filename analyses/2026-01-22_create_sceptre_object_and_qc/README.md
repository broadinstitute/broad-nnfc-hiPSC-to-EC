
# Render day 0 analysis

Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 20,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day0_grna20"
)'

Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 0,
    grna_threshold = 40,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day0_grna40"
)'


# Render day 2 analysis
Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 2,
    grna_threshold = 20,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day2_grna20"
)'

Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 2,
    grna_threshold = 40,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day2_grna40"
)'

# Render day 4 analysis
Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 4,
    grna_threshold = 20,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day4_grna20"
)'

Rscript -e 'rmarkdown::render(
  "analyses/2026-01-22_create_sceptre_object_and_qc/create-sceptre-object-and-qc.Rmd",
  params = list(
    selected_day = 4,
    grna_threshold = 40,
    threads = 8
  ),
  output_dir = "results/2026-01-22_create_sceptre_object_and_qc/day4_grna40"
)'

