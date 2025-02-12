# sites <- "Mbaiki" #nolint
sites <- unique(readr::read_tsv("data/derived_data/data.tsv")$site) # nolint
for (site in sites) {
  print(site)
  file_name <- paste0(site, "_trajectories.pdf")
  quarto::quarto_render(
    input = "analyses/trajectories.qmd",
    output_file = file_name,
    execute_params = list(site = site)
  )
  file.rename(
    from = file_name,
    to = file.path(
      "outputs",
      "trajectories",
      file_name
    )
  )
}

sites <- c(
  "Kabo", "Lesong", "Malinau", "Mbaiki", "Misiones",
  "Paracou", "Sg Lalang", "Ulu Muda"
)
for (site in sites) {
  print(site)
  file_name <- paste0(site, "_wg1.pdf")
  quarto::quarto_render(
    input = "analyses/wg1.qmd",
    output_file = file_name,
    execute_params = list(site = site)
  )
  file.rename(
    from = file_name,
    to = file.path(
      "outputs",
      "wg1",
      file_name
    )
  )
}
