sites <- "Mbaiki"
# sites <- unique(readr::read_tsv("data/derived_data/data.tsv")$site) #nolint
for (site in sites) {
  file_name <- paste0(site, "_trajectories.pdf")
  quarto::quarto_render(
    input = "analyses/trajectories.qmd",
    output_file = file_name,
    execute_params = list(site = site)
  )
  file.rename(
    from = file_name,
    to = file.path("outputs",
                   "trajectories",
                   file_name)
  )
}
