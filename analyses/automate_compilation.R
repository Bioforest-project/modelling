library(tidyverse)

# automate qmd compilation
raw_data_sites <- list.files("data/raw_data", "harmonized") |>
  gsub(
    pattern = "harmonized_data_|_v1.csv|_temoin|_t[0-9]|[0-9]",
    replacement = ""
  ) |>
  unique()

in_plot_info <- read.csv("data/raw_data/bioforest-plot-information.csv") |>
  subset(!is.na(longitude) & !is.na(Treatment) & Treatment != "") |>
  select(site) |>
  unlist() |>
  tolower() |>
  gsub(pattern = " ", replacement = "_") |>
  gsub(pattern = "_km_|[0-9]", replacement = "") |>
  gsub(pattern = "sg_", replacement = "sungai_") |>
  iconv(to = "ASCII//TRANSLIT") |>
  unique()

compile_sites <- intersect(raw_data_sites, in_plot_info)

# cache: if we don't want to redo the compilation for files that already exist

cache <- TRUE

if (cache) {
  done <- list.files("reports/") |>
    gsub(pattern = "data_aggregation_|.pdf", replacement = "")

  compile_sites <- setdiff(compile_sites, done)
}

if (length(compile_sites) > 0) {
  for (s in compile_sites) {
    # render document
    quarto::quarto_render(
      input = "analyses/data_aggregation.qmd",
      output_format = "all",
      output_file = paste0("data_aggregation_", s, ".pdf"),
      execute_params = list(site = s)
    )
    # move to "outputs/data_aggregation_reports" folder
    file.rename(
      from = paste0("data_aggregation_", s, ".pdf"),
      to = paste0("reports/data_aggregation_", s, ".pdf")
    )
  }
}
